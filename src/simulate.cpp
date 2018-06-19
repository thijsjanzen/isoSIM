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

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

std::vector< Fish > simulate_Population(const std::vector< Fish>& sourcePop,
                                        const NumericMatrix& select,
                                        int pop_size,
                                        int total_runtime,
                                        double morgan,
                                        bool progress_bar,
                                        arma::cube& frequencies,
                                        bool track_frequency,
                                        bool track_junctions,
                                        std::vector<double>& junctions) {

    bool use_selection = FALSE;
    if(select(1, 1) > -1e4) use_selection = TRUE;


    double expected_max_fitness = 1e-6;
    std::vector<Fish> Pop = sourcePop;
    std::vector<double> fitness;
    double maxFitness = -1;

    if(use_selection) {
        for(int j = 0; j < select.nrow(); ++j) {
            if(select(j, 4) < 0) break; // these entries are only for tracking, not for selection calculations
            double local_max_fitness = 0.0;
            for(int i = 1; i < 4; ++i) {
                if(select(j, i) > local_max_fitness) {
                    local_max_fitness = select(j, i);
                }
            }
            expected_max_fitness += local_max_fitness;
        }

        for(auto it = Pop.begin(); it != Pop.end(); ++it){
            double fit = calculate_fitness_twoAllele((*it), select);
            if(fit > maxFitness) maxFitness = fit;

            if(fit > (expected_max_fitness)) { // little fix to avoid numerical problems
                Rcout << "Expected maximum " << expected_max_fitness << " found " << fit << "\n";
                Rcpp::stop("ERROR in calculating fitness, fitness too large\n");
            }

            fitness.push_back(fit);
        }
    }

    int updateFreq = total_runtime / 20;
    if(updateFreq < 1) updateFreq = 1;

    if(progress_bar) {
        Rcout << "0--------25--------50--------75--------100\n";
        Rcout << "*";
    }

    for(int t = 0; t < total_runtime; ++t) {

        if(track_junctions) junctions.push_back(calc_mean_junctions(Pop));

        if(track_frequency) {
            for(int i = 0; i < select.nrow(); ++i) {
                Rcout << "updating frequencies\n";
                arma::mat x = frequencies.slice(i);
                NumericVector v = update_frequency(Pop, select(i, 0), x.n_cols);

                for(int j = 0; j < v.size(); ++j) {
                    x(t, j) = v(j);
                }

                frequencies.slice(i) = x;
                Rcout << "frequencies updated\n";
            }
        }

        std::vector<Fish> newGeneration;
        std::vector<double> newFitness;
        double newMaxFitness = -1.0;
        Rcout << "updating fish\n";
        for(int i = 0; i < pop_size; ++i)  {
            int index1 = 0;
            int index2 = 0;
            if(use_selection) {
                index1 =  draw_prop_fitness(fitness, maxFitness);
                index2 = draw_prop_fitness(fitness, maxFitness);
                while(index2 == index1) index2 = draw_prop_fitness(fitness, maxFitness);
            } else {
                index1 = random_number( (int)Pop.size() );
                index2 = random_number( (int)Pop.size() );
                while(index2 == index1) index2 = random_number( (int)Pop.size() );
            }

            Fish kid = mate(Pop[index1], Pop[index2], morgan);

            newGeneration.push_back(kid);

            double fit = -2.0;
            if(use_selection) fit = calculate_fitness_twoAllele(kid, select);
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
        Rcout << "done updating, again!\n";
    }
    if(progress_bar) Rcout << "\n";
    return(Pop);
}


// [[Rcpp::export]]
List simulate_cpp(Rcpp::NumericVector input_population,
              NumericMatrix select,
              int pop_size,
              int number_of_founders,
              int total_runtime,
              double morgan,
              bool progress_bar,
              bool track_frequency,
              bool track_junctions)
{
    std::vector< Fish > Pop;
    int number_of_alleles = number_of_founders;

    if(input_population[0] > -1e4) {
        Rcout << "Found input population! converting!\n";
        Pop = convert_NumericVector_to_fishVector(input_population);

        number_of_founders = 0;

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
        number_of_alleles = number_of_founders + 1;
        Rcout << "Number of alleles calculated\n";
    } else {
         for(int i = 0; i < pop_size; ++i) {
            Fish p1 = Fish( random_number( number_of_founders ) );
            Fish p2 = Fish( random_number( number_of_founders ) );

            Pop.push_back(mate(p1,p2, morgan));
        }
    }

    arma::cube frequencies_table;

    if(track_frequency) {
        int number_entries = select.nrow();
        arma::cube x(total_runtime, number_of_alleles, number_entries); // n_row, n_col, n_slices, type
        frequencies_table = x;
    }

   arma::mat initial_frequencies = update_all_frequencies(Pop, select, number_of_alleles);

    std::vector<double> junctions;

    std::vector<Fish> outputPop = simulate_Population(Pop,
                                                      select,
                                                      pop_size,
                                                      total_runtime,
                                                      morgan,
                                                      progress_bar,
                                                      frequencies_table,
                                                      track_frequency,
                                                      track_junctions,
                                                      junctions);

    arma::mat final_frequencies = update_all_frequencies(outputPop, select, number_of_alleles);

    return List::create( Named("population") = convert_to_list(outputPop),
                        Named("frequencies") = frequencies_table,
                        Named("initial_frequencies") = initial_frequencies,
                        Named("final_frequencies") = final_frequencies,
                        Named("junctions") = junctions);
}
