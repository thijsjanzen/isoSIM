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
#include "random_functions.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;

NumericVector update_frequency(const std::vector< Fish >& v,
                               double m,
                               int num_alleles) {

    NumericVector freq(num_alleles, 0.0);

    for(auto it = v.begin(); it != v.end(); ++it) {
        for(auto i = ((*it).chromosome1.begin()+1); i != (*it).chromosome1.end(); ++i) {
            if((*i).pos > m) {
                int index = (*(i-1)).right;
                if(index >= num_alleles || index < 0) {
                    Rcout << "ERROR!!\n";
                    Rcout << "trying to access NumericVector freq outside bounds\n";
                    Rcout << "in update_frequency\n";
                    Rcout << index << "\t" << num_alleles << "\t" << freq.size() << "\n";
                    Rcout << (*i).pos << "\t" << m << "\t" << (*it).chromosome1.size() << "\n";
                    Rcout << (*(i-1)).pos << "\t" << (*(i-1)).right << "\t" << (*it).chromosome1.empty() << "\n";
                }
                freq(index)++;
                break;
            }
        }

        for(auto i = ((*it).chromosome2.begin()+1); i != (*it).chromosome2.end(); ++i) {
            if((*i).pos > m) {
                int index = (*(i-1)).right;
                if(index >= num_alleles || index < 0) {
                    Rcout << "ERROR!!\n";
                    Rcout << "trying to access NumericVector freq outside bounds\n";
                    Rcout << "in update_frequency\n";
                    Rcout << index << "\t" << num_alleles << "\t" << freq.size() << "\n";
                    Rcout << (*i).pos << "\t" << m << "\t" << (*it).chromosome2.size() << "\n";
                    Rcout << (*(i-1)).pos << "\t" << (*(i-1)).right << "\t" << (*it).chromosome2.empty() << "\n";
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

arma::mat update_all_frequencies(const std::vector< Fish >& pop,
                                 const NumericVector& markers,
                                 int number_of_founders) {

    arma::mat output(markers.size(), number_of_founders);

    for(int i = 0; i < markers.size(); ++i) {
        NumericVector v = update_frequency(pop,
                                           markers[i],
                                           number_of_founders);
        for(int j = 0; j < v.size(); ++j) {
            output(i, j) = v(j);
        }
    }
    return(output);
}

double calc_mean_junctions(const std::vector< Fish> & pop) {

    double mean_junctions = 0.0;
    for(auto it = pop.begin(); it != pop.end(); ++it) {
        mean_junctions += (*it).chromosome1.size() - 2; // start and end don't count
        mean_junctions += (*it).chromosome2.size() - 2;
    }
    mean_junctions *= 1.0 / (pop.size() * 2); // diploid

    return(mean_junctions);
}

int draw_prop_fitness(const std::vector<double> fitness,
                      double maxFitness) {

    if(maxFitness <= 0.0) {
        Rcout << "maxFitness = " << maxFitness << "\n";
        Rcpp::stop("Cannot draw fitness if maxFitness <= 0");
        return(-1);
    }

    if(maxFitness > 10000.0) {
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




std::vector< Fish > convert_NumericVector_to_fishVector(const NumericVector v) {
    std::vector< Fish > output;

    Fish temp;
    int indic_chrom = 1;
    bool add_indiv = false;

    for(int i = 0; i < v.size(); i += 2) {
        junction temp_j;
        temp_j.pos = v[i];
        temp_j.right = v[i+1];

        if(indic_chrom == 1) {
            temp.chromosome1.push_back(temp_j);
        } else {
            temp.chromosome2.push_back(temp_j);
        }

        if(temp_j.right == -1) {
            if(indic_chrom == 1) {
                indic_chrom = 2;
            } else {
                add_indiv = true;
            }
        }

        if(add_indiv) {
            output.push_back(temp);
            add_indiv = false;
            indic_chrom = 1;
            temp.chromosome1.clear();
            temp.chromosome2.clear();
        }
    }

    return(output);
}

List convert_to_list(const std::vector<Fish>& v) {
    int list_size = (int)v.size();
    List output(list_size);

    for(int i = 0; i < v.size(); ++i) {

        Fish focal = v[i];

        NumericMatrix chrom1(focal.chromosome1.size(), 2); // nrow = number of junctions, ncol = 2
        for(int j = 0; j < focal.chromosome1.size(); ++j) {
            chrom1(j, 0) = focal.chromosome1[j].pos;
            chrom1(j, 1) = focal.chromosome1[j].right;
        }

        NumericMatrix chrom2(focal.chromosome2.size(), 2); // nrow = number of junctions, ncol = 2
        for(int j = 0; j < focal.chromosome2.size(); ++j) {
            chrom2(j, 0) = focal.chromosome2[j].pos;
            chrom2(j, 1) = focal.chromosome2[j].right;
        }

        List toAdd = List::create( Named("chromosome1") = chrom1,
                                   Named("chromosome2") = chrom2
                                  );

        output(i) = toAdd;
    }

    return output;
}




double calculate_fitness(const Fish& focal,
                         const NumericMatrix& select,
                         bool multiplicative_selection) {

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

    if(!multiplicative_selection) {
        for(int i = 0; i < num_alleles.size(); ++i) {
            if(select(i, 4) < 0) break; // these entries are only for tracking alleles over time, not for selection calculation

            int fitness_index = 1 + num_alleles[i];
            fitness += select(i, fitness_index);
        }
    }

    if(multiplicative_selection) {
        fitness = 1.0;
        for(int i = 0; i < num_alleles.size(); ++i) {
            if(select(i, 4) < 0) break; // these entries are only for tracking alleles over time, not for selection calculation

            int fitness_index = 1 + num_alleles[i];
            fitness *= select(i, fitness_index);
        }
    }

    return(fitness);
}

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
            double fit = calculate_fitness((*it), select, multiplicative_selection);
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
            if(use_selection) fit = calculate_fitness(kid, select, multiplicative_selection);
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




