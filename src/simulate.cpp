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
                                        const NumericVector& track_markers,
                                        bool track_junctions,
                                        std::vector<double>& junctions,
                                        bool multiplicative_selection,
                                        int num_alleles) {

    bool use_selection = FALSE;
    if(select(1, 1) >= 0) use_selection = TRUE;


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
            double fit = calculate_fitness_twoAllele((*it), select, multiplicative_selection);
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
            for(int i = 0; i < track_markers.size(); ++i) {
                arma::mat x = frequencies.slice(i);
                if(track_markers[i] < 0) break;
                NumericVector v = update_frequency(Pop, track_markers[i], num_alleles);
                for(int j = 0; j < v.size(); ++j) {
                    x(t, j) = v(j);
                }
                frequencies.slice(i) = x;
            }
        }

        std::vector<Fish> newGeneration;
        std::vector<double> newFitness;
        double newMaxFitness = -1.0;
        //Rcout << "updating fish\n";
        for(int i = 0; i < pop_size; ++i)  {
            int index1 = 0;
            int index2 = 0;
            if(use_selection) {
                index1 = draw_prop_fitness(fitness, maxFitness);
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
            if(use_selection) fit = calculate_fitness_twoAllele(kid, select, multiplicative_selection);
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
   //     Rcout << "done updating, again!\n";
    }
    if(progress_bar) Rcout << "\n";
    return(Pop);
}

int draw_random_founder(const std::vector<double>& v) {
    double r = uniform();
    for(int i = 0; i < v.size(); ++i) {
        r -= v[i];
        if(r <= 0) {
            return(i);
        }
    }
    return(v.back());
}

// [[Rcpp::export]]
List simulate_cpp(Rcpp::NumericVector input_population,
              NumericMatrix select,
              int pop_size,
              int number_of_founders,
              Rcpp::NumericVector starting_proportions,
              int total_runtime,
              double morgan,
              bool progress_bar,
              bool track_frequency,
              NumericVector track_markers,
              bool track_junctions,
              bool multiplicative_selection)
{
    std::vector< Fish > Pop;
    int number_of_alleles = number_of_founders;

    if(input_population[0] > -1e4) {
     //   Rcout << "Found input population! converting!\n";
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
       // Rcout << "Number of alleles calculated\n";
    } else {
        std::vector<double> starting_freqs = as< std::vector<double> >(starting_proportions);
        for(int i = 0; i < pop_size; ++i) {
            int founder_1 = draw_random_founder(starting_freqs);
            int founder_2 = draw_random_founder(starting_freqs);

            Fish p1 = Fish( founder_1 );
            Fish p2 = Fish( founder_2 );

            Pop.push_back(mate(p1,p2, morgan));
        }
    }

    arma::cube frequencies_table;

    if(track_frequency) {
        //Rcout << "Preparing frequencies_table\n";
        int number_entries = track_markers.size();
        arma::cube x(total_runtime, number_of_alleles, number_entries); // n_row, n_col, n_slices, type
        frequencies_table = x;
    }

    arma::mat initial_frequencies = update_all_frequencies(Pop, track_markers, number_of_alleles);

    std::vector<double> junctions;
  //  Rcout << "starting simulation\n";
    std::vector<Fish> outputPop = simulate_Population(Pop,
                                                      select,
                                                      pop_size,
                                                      total_runtime,
                                                      morgan,
                                                      progress_bar,
                                                      frequencies_table,
                                                      track_frequency,
                                                      track_markers,
                                                      track_junctions,
                                                      junctions,
                                                      multiplicative_selection,
                                                      number_of_alleles);

    arma::mat final_frequencies = update_all_frequencies(outputPop, track_markers, number_of_alleles);

    return List::create( Named("population") = convert_to_list(outputPop),
                         Named("frequencies") = frequencies_table,
                         Named("initial_frequencies") = initial_frequencies,
                         Named("final_frequencies") = final_frequencies,
                         Named("junctions") = junctions);
}

// [[Rcpp::export]]
List create_pop_admixed_cpp(int num_individuals,
                  int num_ancestors,
                  int population_size,
                  double size_in_morgan) {

    double p = 1.0 / num_ancestors;
    double init_heterozygosity = 2*p*(1-p);

    int max_num_j = 2*init_heterozygosity * population_size * size_in_morgan;

    std::vector< Fish > output;
    for(int i = 0; i < num_individuals; ++i) {
        Fish focal(random_number(num_ancestors));
        focal.chromosome1.pop_back();
        focal.chromosome2.pop_back();
        double pos = 0.0;
        while(pos < 1) {
            double u = uniform();
            double lambda = max_num_j;
            double exp_u = (-1.0 / lambda) * log(u);
            pos += exp_u;
            if(pos < 1) {
                junction to_add(pos, random_number(num_ancestors));
                focal.chromosome1.push_back(to_add);
            } else {
                junction to_add(1.0, -1);
                focal.chromosome1.push_back(to_add);
            }
        }

        pos = 0.0;
        while(pos < 1) {
            double u = uniform();
            double lambda = max_num_j;
            double exp_u = (-1.0 / lambda) * log(u);
            pos += exp_u;
            if(pos < 1) {
                junction to_add(pos, random_number(num_ancestors));
                focal.chromosome2.push_back(to_add);
            } else {
                junction to_add(1.0, -1);
                focal.chromosome2.push_back(to_add);
            }
        }

        output.push_back(focal);
    }

    return List::create( Named("population") = convert_to_list(output));
}




