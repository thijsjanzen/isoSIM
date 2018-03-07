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
// #include "random_functions.h"

#include "randomc.h"

#include <Rcpp.h>
using namespace Rcpp;

int draw_prop_fitness(const std::vector<double> fitness,
                      double maxFitness) {

    if(maxFitness < 0.0) {
        //   Rcout << "Cannot draw fitness if maxFitness < 0, terminating\n";
        return(-1);
    }

    for(int i = 0; i < 1e6; ++i) {
        int index = random_number(fitness.size());
        double prob = 1.0 * fitness[index] / maxFitness;
        if(uniform() < prob) {
            return index;
        }
    }
    Rcout << "\nERROR! ERROR! Couldn't pick proportional to fitness\t" << maxFitness << "\n";
    return -1;
}

double assess_match(const std::vector<junction>& chrom,
                    double start,
                    double end,
                    int ancestor) {

    std::vector< junction > block;

    for(int i = 0; i < chrom.size(); ++i) {
        if(chrom[i].pos > end) { // stop after block
            break;
        } else {
            if(chrom[i].pos > start) {
                block.push_back(chrom[i]); //junctions within the block should be stored
            } else {
                if(i+1 < chrom.size()) {
                    if(chrom[i+1].pos > start) { //one junction before the end as well.
                        block.push_back(chrom[i]);
                    }
                }
            }
        }
    }

    if(block.size() == 1) {
        if(block[0].right == ancestor) {
            return 1.0;
        } else {
            return 0.0;
        }
    }


    double match = 0.0;

    for(int i = 0; i < block.size(); ++i) {
        if(block[i].right == ancestor) {

            double local_right = end;
            if(i+1 < block.size()) {
                local_right = block[i+1].pos;
            }

            double local_left = block[i].pos;
            if(local_left < start) local_left = start;

            match += local_right - local_left;
        }
    }

    match *= 1.0 / (end - start);

    return(match);
}

double calculate_fitness(const Fish& focal,
                         const std::vector< std::vector< double > >& select,
                         double s) {

    // double fitness = 2.0 * select[0][0];
    double fitness = 2.0;

    for(int i = 0; i < select.size(); ++i) {
        double start = select[i][0];
        double end = select[i][1];
        int ancestor = select[i][2];

        double a1 = assess_match(focal.chromosome1, start, end, ancestor);
        double a2 = assess_match(focal.chromosome2, start, end, ancestor);

        fitness += (end - start) * (s * (a1 + a2));

        // double nextStart = 1.0;
        // if((i+1) < select.size()) {
        //     nextStart = select[i+1][0];
        // }
        // fitness += 2.0 * (nextStart - end);
    }
    
    fitness = fitness / 2.0;
    return fitness;
}



std::vector< Fish > selectPopulation(const std::vector< Fish>& sourcePop,
                                     const std::vector< std::vector< double > >& select,
                                     double s,
                                     int popSize,
                                     int maxTime,
                                     double Morgan)
{

    std::vector<Fish> Pop = sourcePop;
    std::vector<double> fitness;
    double maxFitness = -1e6;

    for(auto it = Pop.begin(); it != Pop.end(); ++it){
        double fit = calculate_fitness((*it), select, s);
        if(fit < 0.0) {
            Rcout << "ERROR in calculating fitness\n"; 
        }

        if(fit > maxFitness) maxFitness = fit;
        fitness.push_back(fit);
    }

    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";

    int updateFreq = maxTime / 20;
    if(updateFreq < 1) updateFreq = 1;



    for(int t = 0; t < maxTime; ++t) {

        std::vector<Fish> newGeneration;
        std::vector<double> newFitness;
        double newMaxFitness = - 1.0;

        for(int i = 0; i < popSize; ++i)  {
            int index1 = draw_prop_fitness(fitness, maxFitness);
            int index2 = draw_prop_fitness(fitness, maxFitness);
            while(index2 == index1) index2 = draw_prop_fitness(fitness, maxFitness);

            Fish kid = mate(Pop[index1], Pop[index2], Morgan);

            newGeneration.push_back(kid);

            double fit = calculate_fitness(kid, select, s);
            if(fit > newMaxFitness) newMaxFitness = fit;
            newFitness.push_back(fit);
        }

        Pop = newGeneration;
        newGeneration.clear();
        fitness = newFitness;
        maxFitness = newMaxFitness;

        if(t % updateFreq == 0) {
            Rcout << "**";
        }
        Rcpp::checkUserInterrupt();
    }
    
    return(Pop);
}



// [[Rcpp::export]]
List select_population_cpp(Rcpp::NumericVector v1,
                           Rcpp::NumericVector selectM,
                           double s,
                           int population_size,
                           int run_time,
                           double morgan,
                           int seed) {

    set_seed(seed);
    std::vector< Fish > Pop = convert_NumericVector_to_fishVector(v1);

    std::vector< std::vector< double > > select;
    std::vector<double> temp_select;
    for(int i = 0; i < selectM.size(); ++i) {
        temp_select.push_back(selectM[i]);
        if(temp_select.size() == 3) {
            select.push_back(temp_select);
            temp_select.clear();
        }
    }

    std::vector<Fish> outputPop = selectPopulation(Pop,
                                                   select,
                                                   s,
                                                   population_size,
                                                   run_time,
                                                   morgan);

    return List::create( Named("population") = convert_to_list(outputPop) );
}

// [[Rcpp::export]]
List create_population_selection_cpp(int pop_size,
                                 int number_of_founders,
                                 int total_runtime,
                                 double morgan,
                                 Rcpp::NumericVector select_matrix,
                                 double selection,
                                 int seed) {

    set_seed(seed);
    std::vector< Fish > Pop;
    for(int i = 0; i < pop_size; ++i) {
        Fish p1 = Fish( random_number( number_of_founders ) );
        Fish p2 = Fish( random_number( number_of_founders ) );

        Pop.push_back(mate(p1,p2, morgan));
    }

    std::vector< std::vector< double > > select;
    std::vector<double> temp_select;
    for(int i = 0; i < select_matrix.size(); ++i) {
        temp_select.push_back(select_matrix[i]);
        if(temp_select.size() == 3) {
            select.push_back(temp_select);
            temp_select.clear();
        }
    }

    std::vector<Fish> outputPop = selectPopulation(Pop,
                                                   select,
                                                   selection,
                                                   pop_size,
                                                   total_runtime,
                                                   morgan);

    return List::create( Named("population") = convert_to_list(outputPop) );
}

