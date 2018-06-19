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




double calculate_fitness_twoAllele(const Fish& focal,
                                   const NumericMatrix& select) {


    int number_of_markers = select.nrow();
    std::vector< int > num_alleles(number_of_markers, 0);

    int focal_marker = 0;
    double pos = select(focal_marker, 0);
    double anc = select(focal_marker, 4);
    // loc aa  Aa  AA ancestor
    //  0  1   2  3  4

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
        if(anc < 0) break; // these entries are only for tracking alleles over time, not for selection calculation
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
        if(anc < 0) break; // these entries are only for tracking alleles over time, not for selection calculation
    }

    double fitness = 0.0;
    for(int i = 0; i < num_alleles.size(); ++i) {
        if(select(i, 4) < 0) break; // these entries are only for tracking alleles over time, not for selection calculation

        int fitness_index = 1 + num_alleles[i];
        fitness += select(i, fitness_index);
    }

    return(fitness);
}

std::vector< Fish > selectPopulation(const std::vector< Fish>& sourcePop,
                                                const NumericMatrix& select,
                                                int pop_size,
                                                int total_runtime,
                                                double morgan,
                                                bool progress_bar,
                                                arma::cube& frequencies,
                                                bool track_frequency) {

    double expected_max_fitness = 1e-6;
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

    std::vector<Fish> Pop = sourcePop;
    std::vector<double> fitness;
    double maxFitness = -1;
    for(auto it = Pop.begin(); it != Pop.end(); ++it){
        double fit = calculate_fitness_twoAllele((*it), select);
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
           for(int i = 0; i < select.nrow(); ++i) {
               if(select(i, 4) < 0) break;
                arma::mat x = frequencies.slice(i);
                NumericVector v = update_frequency(Pop, select(i, 0), x.n_cols);

                for(int j = 0; j < v.size(); ++j) {
                    x(t, j) = v(j);
                }
           
                frequencies.slice(i) = x;
            }
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

            double fit = calculate_fitness_twoAllele(kid, select);
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

arma::mat update_all_frequencies(const std::vector< Fish >& pop,
                                 const NumericMatrix& select_matrix,
                                 int number_of_founders) {

    arma::mat output(select_matrix.nrow(), number_of_founders);

    for(int i = 0; i < select_matrix.nrow(); ++i) {
        NumericVector v = update_frequency(pop,
                                           select_matrix(i, 0),
                                           number_of_founders);
        for(int j = 0; j < v.size(); ++j) {
            output(i, j) = v(j);
        }
    }
    return(output);
}



// [[Rcpp::export]]
List create_population_selection_cpp(NumericMatrix select,
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

    arma::cube frequencies_table;

    if(track_frequency) {
        int number_entries = select.nrow();
        arma::cube x(total_runtime, number_of_founders, number_entries); // n_row, n_col, n_slices, type
        frequencies_table = x;
    }

    arma::mat initial_frequencies = update_all_frequencies(Pop, select, number_of_founders);

    std::vector<Fish> outputPop = selectPopulation(Pop,
                                                          select,
                                                          pop_size,
                                                          total_runtime,
                                                          morgan,
                                                          progress_bar,
                                                          frequencies_table,
                                                          track_frequency);

    arma::mat final_frequencies = update_all_frequencies(outputPop, select, number_of_founders);



    return List::create( Named("population") = convert_to_list(outputPop),
                         Named("frequencies") = frequencies_table,
                         Named("initial_frequencies") = initial_frequencies,
                         Named("final_frequencies") = final_frequencies);
}

// [[Rcpp::export]]
List select_population_cpp(Rcpp::NumericVector v1,
                           Rcpp::NumericMatrix selectM,
                           int population_size,
                           int run_time,
                           double morgan,
                           bool progress_bar,
                           bool track_frequency) {

    std::vector< Fish > Pop = convert_NumericVector_to_fishVector(v1);

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

    int num_alleles = 1 + number_of_founders; // because of counting from 0

    arma::cube frequencies_table;
    if(track_frequency) {

        int number_entries = selectM.nrow();
        arma::cube x(run_time, num_alleles, number_entries); // n_row, n_col, n_slices, type
        frequencies_table = x;
    }

    arma::mat initial_frequencies = update_all_frequencies(Pop, selectM, num_alleles);

    std::vector<Fish> outputPop = selectPopulation(Pop,
                                                   selectM,
                                                   population_size,
                                                   run_time,
                                                   morgan,
                                                   progress_bar,
                                                   frequencies_table,
                                                   track_frequency);

  //  Rcout << "simulation done, starting calculation of final frequencies\n";
    arma::mat final_frequencies = update_all_frequencies(outputPop, selectM, num_alleles);

    return List::create( Named("population") = convert_to_list(outputPop),
                        Named("frequencies") = frequencies_table,
                        Named("initial_frequencies") = initial_frequencies,
                        Named("final_frequencies") = final_frequencies);
}


int draw_prop_fitness(const std::vector<double> fitness,
                      double maxFitness) {

    if(maxFitness <= 0.0) {
        Rcout << "maxFitness = " << maxFitness << "\n";
        Rcpp::stop("Cannot draw fitness if maxFitness <= 0");
        return(-1);
    }

    if(maxFitness > 100.0) {
        Rcout << "maxFitness = " << maxFitness << "\n";
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

NumericVector update_frequency(const std::vector< Fish >& v,
                               double m,
                               int num_alleles) {

    NumericVector freq(num_alleles, 0.0);

    for(auto it = v.begin(); it != v.end(); ++it) {
        for(auto i = ((*it).chromosome1.begin()+1); i != (*it).chromosome1.end(); ++i) {
            if((*i).pos > m) {
                int index = (*(i-1)).right;
                if(index >= num_alleles || index < 0) {
                    Rcout << index << "\t" << num_alleles << "\t" << freq.size() << "\n";
                }
                freq(index)++;
                break;
            }
        }

        for(auto i = ((*it).chromosome2.begin()+1); i != (*it).chromosome2.end(); ++i) {
            if((*i).pos > m) {
                int index = (*(i-1)).right;
                if(index >= num_alleles || index < 0) {
                     Rcout << index << "\t" << num_alleles << "\t" << freq.size() << "\n";
                }
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
