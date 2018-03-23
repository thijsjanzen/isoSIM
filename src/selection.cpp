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

//#include "randomc.h"

#include <Rcpp.h>
using namespace Rcpp;

int draw_prop_fitness(const std::vector<double> fitness,
                      double maxFitness) {

    if(maxFitness <= 0.0) {
        Rcpp::stop("Cannot draw fitness if maxFitness <= 0");
        return(-1);
    }

    if(maxFitness > 100.0) {
        Rcout << "maxFitness\n";
        Rcpp::stop("It appears maxfitness has encountered a memory access violation\n");
        return(-1);
    }

    for(int i = 0; i < 1e6; ++i) {
        int index = random_number(fitness.size());
        double prob = 1.0 * fitness[index] / maxFitness;
        if(uniform() < prob) {
            return index;
        }
    }
    Rcout << maxFitness << "\n";
    Rcpp::stop("ERROR!Couldn't pick proportional to fitness");
    return -1;
}


double calculate_fitness_markers(const Fish& focal,
                                 const NumericMatrix& select) {

    double fitness = 2.0;
    int number_of_markers = select.nrow();

    
    int focal_marker = 0;
    double pos = select(focal_marker, 0);
    double anc = select(focal_marker, 1);
    double s = select(focal_marker, 2);

    for(auto it = (focal.chromosome1.begin()+1); it != focal.chromosome1.end(); ++it) {
        if((*it).pos > pos) {
            if((*(it-1)).right == anc) fitness += s;
            focal_marker++;
            if(focal_marker >= number_of_markers) {
                break;
            }
            pos = select(focal_marker, 0);
            anc = select(focal_marker, 1);
            s = select(focal_marker, 2);
        }

    }

    focal_marker = 0;
    pos = select(focal_marker, 0);
    anc = select(focal_marker, 1);
    s = select(focal_marker, 2);

    for(auto it = (focal.chromosome2.begin()+1); it != focal.chromosome2.end(); ++it) {
        if((*it).pos > pos) {
            if((*(it-1)).right == anc) fitness += s;
            focal_marker++;
            if(focal_marker >= number_of_markers) {
                break;
            }
            pos = select(focal_marker, 0);
            anc = select(focal_marker, 1);
            s = select(focal_marker, 2);

        }
    }

    return(fitness / 2.0);
}

NumericVector update_frequency(const std::vector< Fish >& v, double m, int num_alleles) {
    NumericVector freq(num_alleles);

    for(auto it = v.begin(); it != v.end(); ++it) {
        for(auto i = (*it).chromosome1.begin(); i != (*it).chromosome1.end(); ++i) {
            if((*i).pos > m) {
                int index = (*(i-1)).right;
                freq(index)++;
                break;
            }
        }

        for(auto i = (*it).chromosome2.begin(); i != (*it).chromosome2.end(); ++i) {
            if((*i).pos > m) {
                int index = (*(i-1)).right;
                freq(index)++;
                break;
            }
        }
    }

    for(int i = 0; i < freq.size(); ++i) {
        freq(i) = freq(i) * 1.0 / (2*v.size());
    }

    return(freq);
}

std::vector< Fish > selectPopulation_vector(const std::vector< Fish>& sourcePop,
                                            const NumericMatrix& select,
                                            int pop_size,
                                            int total_runtime,
                                            double morgan,
                                            bool progress_bar,
                                            NumericMatrix& frequencies,
                                            bool track_frequency) {

    double expected_max_fitness = 1.0;

    Rcout << "Applying selection for:\n";
    Rcout << "Location" << "\t" << "Ancestor" << "\t" << "Selection coefficient" << "\n";
    for(int i = 0; i < select.nrow(); ++i) {
        Rcout << select(i, 0) << "\t\t" << select(i, 1) << "\t\t\t" << select(i, 2) << "\n";
        expected_max_fitness += select(i, 2);
    }

    std::vector<Fish> Pop = sourcePop;
    std::vector<double> fitness;
    double maxFitness = -1;
    for(auto it = Pop.begin(); it != Pop.end(); ++it){
        double fit = calculate_fitness_markers((*it), select);
        if(fit > maxFitness) maxFitness = fit;

        if(fit > expected_max_fitness) {
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
List create_population_selection_markers_cpp(NumericMatrix select,
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

    std::vector<Fish> outputPop = selectPopulation_vector(Pop,
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
List select_population_markers_cpp(Rcpp::NumericVector v1,
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


    std::vector<Fish> outputPop = selectPopulation_vector(Pop,
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
