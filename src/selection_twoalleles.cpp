//
//  selection.cpp
//  
//
//  Created by Thijs Janzen on 28/02/2018.
//
//

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <numeric>
#include <cmath>

#include <vector>
#include <algorithm>

#include "Fish.h"
#include "main.h"
#include "random_functions.h"
#include "selection.h"


//#include "randomc.h"

#include <Rcpp.h>
using namespace Rcpp;

double calculate_fitness_twoAllele(const Fish& focal,
                                   const NumericMatrix& select) {
    int number_of_markers = select.nrow();

    std::vector< int > num_alleles(number_of_markers, 0);
    
    int focal_marker = 0;
    double pos = select(focal_marker, 0);
    // loc aa  Aa  AA ancestor
    //  0  1   2  3  4


    double anc = select(focal_marker, 4);
   
    for(auto it = (focal.chromosome1.begin()+1); it != focal.chromosome1.end(); ++it) {
        if((*it).pos > pos) {
            if((*(it-1)).right == anc) num_alleles[focal_marker]++;
            focal_marker++;
            if(focal_marker >= number_of_markers) {
                break;
            }
            pos = select(focal_marker, 0);
            anc = select(focal_marker, 4);
        }

    }

    focal_marker = 0;
    pos = select(focal_marker, 0);
    anc = select(focal_marker, 4);

    for(auto it = (focal.chromosome2.begin()+1); it != focal.chromosome2.end(); ++it) {
        if((*it).pos > pos) {
            if((*(it-1)).right == anc) num_alleles[focal_marker]++;
            focal_marker++;
            if(focal_marker >= number_of_markers) {
                break;
            }
            pos = select(focal_marker, 0);
            anc = select(focal_marker, 4);
        }
    }

    double fitness = 0.0;
    for(int i = 0; i < num_alleles.size(); ++i) {
        int fitness_index = 1 + num_alleles[i];
        fitness += select(i, fitness_index);
    }

    return(fitness);
}

std::vector< Fish > selectPopulation_twoAlleles(const std::vector< Fish>& sourcePop,
                                                const NumericMatrix& select,
                                                int pop_size,
                                                int total_runtime,
                                                double morgan,
                                                bool progress_bar,
                                                NumericMatrix& frequencies,
                                                bool track_frequency) {

    double expected_max_fitness = select(1,1) + 1e-6;
    for(int j = 0; j < select.nrow(); ++j) {
        for(int i = 1; i < 4; ++i) {
            if(select(j, i) > expected_max_fitness) {
                expected_max_fitness = select(j, i);
            }
        }
    }

    std::vector<Fish> Pop = sourcePop;
    std::vector<double> fitness;
    double maxFitness = -1;
    for(auto it = Pop.begin(); it != Pop.end(); ++it){
        double fit = calculate_fitness_markers((*it), select);
        if(fit > maxFitness) maxFitness = fit;

        if(fit > (expected_max_fitness)) { // little fix to avoid numerical problems
            Rcout << "Expected maximum " << expected_max_fitness << " found " << fit << "\n";
            Rcpp::stop("ERROR in calculating fitness, fitness too large\n");
        }

        fitness.push_back(fit);
    }

    int updateFreq = total_runtime / 20;
    if(updateFreq < 1) updateFreq = 1;

    if(progress_bar) {
        Rcout << "0--------25--------50--------75--------100\n";
        Rcout << "*";
    }

    for(int t = 0; t < total_runtime; ++t) {

        if(track_frequency) {
            frequencies(t,_) = update_frequency(Pop, select(0, 0), frequencies.ncol());
        }

        std::vector<Fish> newGeneration;
        std::vector<double> newFitness;
        double newMaxFitness = -1.0;

        for(int i = 0; i < pop_size; ++i)  {
            int index1 = draw_prop_fitness(fitness, maxFitness);
            int index2 = draw_prop_fitness(fitness, maxFitness);
            while(index2 == index1) index2 = draw_prop_fitness(fitness, maxFitness);

            Fish kid = mate(Pop[index1], Pop[index2], morgan);

            newGeneration.push_back(kid);

            double fit = calculate_fitness_markers(kid, select);
            if(fit > newMaxFitness) newMaxFitness = fit;

            if(fit > expected_max_fitness) {
                Rcout << "Expected maximum " << expected_max_fitness << " found " << fit << "\n";
                Rcpp::stop("ERROR in calculating fitness, fitness too large\n");
            }

            newFitness.push_back(fit);
        }

        if(t % updateFreq == 0 && progress_bar) {
            Rcout << "**";
        }
        Rcpp::checkUserInterrupt();

        Pop = newGeneration;
        newGeneration.clear();
        fitness = newFitness;
        maxFitness = newMaxFitness;
    }
    Rcout << "\n";
    return(Pop);
}

// [[Rcpp::export]]
List create_population_selection_twoalleles_cpp(NumericMatrix select,
                                                 int pop_size,
                                                 int number_of_founders,
                                                 int total_runtime,
                                                 double morgan,
                                                 bool progress_bar,
                                                 bool track_frequency)
{
    std::vector< Fish > Pop;
    for(int i = 0; i < pop_size; ++i) {
        Fish p1 = Fish( random_number( number_of_founders ) );
        Fish p2 = Fish( random_number( number_of_founders ) );

        Pop.push_back(mate(p1,p2, morgan));
    }

    NumericMatrix frequencies_table;
    if(track_frequency) {
        frequencies_table = NumericMatrix(total_runtime, number_of_founders);
    }

    std::vector<Fish> outputPop = selectPopulation_twoAlleles(Pop,
                                                          select,
                                                          pop_size,
                                                          total_runtime,
                                                          morgan,
                                                          progress_bar,
                                                          frequencies_table,
                                                          track_frequency);
    
    return List::create( Named("population") = convert_to_list(outputPop),
                         Named("frequencies") = frequencies_table);
}

// [[Rcpp::export]]
List select_population_twoalleles_cpp(Rcpp::NumericVector v1,
                           Rcpp::NumericMatrix selectM,
                           int population_size,
                           int run_time,
                           double morgan,
                           bool progress_bar,
                           bool track_frequency) {

    std::vector< Fish > Pop = convert_NumericVector_to_fishVector(v1);

    NumericMatrix frequencies_table;
    if(track_frequency) {
        int number_of_founders = 0;
        for(auto it = Pop.begin(); it != Pop.end(); ++it) {
            for(auto i = (*it).chromosome1.begin(); i != (*it).chromosome1.end(); ++i) {
                if((*i).right > number_of_founders) {
                    number_of_founders = (*i).right;
                }
            }
            for(auto i = (*it).chromosome2.begin(); i != (*it).chromosome2.end(); ++i) {
                if((*i).right > number_of_founders) {
                    number_of_founders = (*i).right;
                }
            }
        }

        frequencies_table = NumericMatrix(run_time, 1 + number_of_founders);
    }


    std::vector<Fish> outputPop = selectPopulation_twoAlleles(Pop,
                                                   selectM,
                                                   population_size,
                                                   run_time,
                                                   morgan,
                                                   progress_bar,
                                                   frequencies_table,
                                                   track_frequency);

    return List::create( Named("population") = convert_to_list(outputPop),
                        Named("frequencies") = frequencies_table);
}
