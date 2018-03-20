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

double assess_match(const std::vector<junction> chrom,
                    double start,
                    double end,
                    int ancestor) {

   /* std::vector< junction > block;
    for(int i = 1; i < chrom.size(); ++i) {
        if(chrom[i].pos > start) {
            block.push_back(chrom[i-1]);

            if(chrom[i].pos > end) {
                break;
            }


        }
    }

    // now we have:
    //  | |  |   |   |     |  |
    // j0 | j1  j2  j3 .. jn  |
    //    s                   e

    // we add a junction indicating the start, and remove the one
    // just before the start

    junction start_corner;
    start_corner.pos = start;
    start_corner.right = block[0].right;
    block[0] = start_corner;

    // we add a junction indicating the end

    junction end_corner;
    end_corner.pos = end;
    end_corner.right = -1;

    block.push_back(end_corner);

    // now we have:
    //  |   |   |  |     | |
    //  s  j1  j2 j3 .. jn e


    double match = 0.0;
    for(int i = 1; i < block.size(); ++i) {
        double stretch = block[i].pos - block[i-1].pos;
        if(block[i-1].right == ancestor) match += stretch;
    }
    
    match *= 1.0 / (end - start);
    
    return(match);
    */
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
                       //  const std::vector< std::vector< double > >& select,
                         const NumericMatrix& select,
                         double s) {

    double fitness = 2.0;

    for(int i = 0; i < select.nrow(); ++i) {
        //double start = select[i][0];
        //double end = select[i][1];
        //int ancestor = select[i][2];
        double start = select(i, 0);
        double end = select(i, 1);
        int ancestor = select(i, 2);


        double a1 = assess_match(focal.chromosome1, start, end, ancestor);

        if(a1 < 0) {
            Rcout << "ERROR! assess_match with chromosome 1 failure\n";
            Rcout << start <<"\t" << end << "\t" << ancestor << "\n";
            for(int i = 0; i < focal.chromosome1.size(); ++i) {
                Rcout << focal.chromosome1[i].pos << "\t" << focal.chromosome1[i].right << "\n";
            }
            Rcpp::stop("ERROR! assess_match with chromosome 1 failure\n");
        }


        double a2 = assess_match(focal.chromosome2, start, end, ancestor);

        if(a2 < 0) {
            Rcout << "ERROR! assess_match with chromosome 2 failure\n";
            Rcout << start <<"\t" << end << "\t" << ancestor << "\n";
            for(int i = 0; i < focal.chromosome2.size(); ++i) {
                Rcout << focal.chromosome2[i].pos << "\t" << focal.chromosome2[i].right << "\n";
            }
            Rcpp::stop("ERROR! assess_match with chromosome 2 failure");
        }

     //   if((end - start) < 0.0) {
     //       Rcout << end << "\t" << start << "\n";
     //       Rcpp::stop("ERROR! end - start < 0");
     //   }

        double to_add = (end - start) * (s * (a1 + a2));
        if(to_add < 0.0) {
       //     Rcout << select[i][0] << "\t" << select[i][1] << "\t" << select[i][2] << "\n";
            Rcout << select(i, 0) << "\t" << select(i, 1) << "\t" << select(i, 2) << "\t" << i << "\n";
            Rcout << start << "\t" << end << "\t" << s << "\t" << a1 << "\t" << a2 << "\n";
            Rcpp::stop("ERROR! Fitness increase negative!");
        }

        fitness += (end - start) * (s * (a1 + a2));
    }
    
    fitness = fitness / 2.0;
    return fitness;
}



std::vector< Fish > selectPopulation(const std::vector< Fish>& sourcePop,
               //                      const std::vector< std::vector< double > >& select,
                                     const NumericMatrix& select,
                                     double s,
                                     int popSize,
                                     int maxTime,
                                     double Morgan)
{
    Rcout << "Applying matrix with:\n";
    Rcout << "Start" << "\t" << "End" << "\t" << "Ancestor" << "\n";
    for(int i = 0; i < select.nrow(); ++i) {
        Rcout << select(i, 0) << "\t" << select(i, 1) << "\t" << select(i, 2) << "\n";
    }



    std::vector<Fish> Pop = sourcePop;
    std::vector<double> fitness;
    double maxFitness = -1.0;

    for(auto it = Pop.begin(); it != Pop.end(); ++it){
        double fit = calculate_fitness((*it), select, s);
        if(fit < 0.0) {
            Rcpp::stop("ERROR in calculating fitness");
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

            double fit = calculate_fitness(kid, select, s);
            if(fit > newMaxFitness) newMaxFitness = fit;

            newGeneration.push_back(kid);
            newFitness.push_back(fit);

            if(fit < 0.0) {
                Rcpp::stop("ERROR in calculating fitness");
            }
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
                           Rcpp::NumericMatrix selectM,
                           double s,
                           int population_size,
                           int run_time,
                           double morgan) {

    std::vector< Fish > Pop = convert_NumericVector_to_fishVector(v1);

   // std::vector< std::vector< double > > select;
   // std::vector<double> temp_select;
  //  for(int i = 0; i < selectM.size(); ++i) {
   //     temp_select.push_back(selectM[i]);
   //     if(temp_select.size() == 3) {
   //         select.push_back(temp_select);
   //         temp_select.clear();
   //     }
   // }

    std::vector<Fish> outputPop = selectPopulation(Pop,
                                                   selectM,
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
                                 Rcpp::NumericMatrix select_matrix,
                                 double selection) {

    std::vector< Fish > Pop;
    for(int i = 0; i < pop_size; ++i) {
        Fish p1 = Fish( random_number( number_of_founders ) );
        Fish p2 = Fish( random_number( number_of_founders ) );

        Pop.push_back(mate(p1,p2, morgan));
    }

  //  std::vector< std::vector< double > > select;
  //  std::vector<double> temp_select;
  //  for(int i = 0; i < select_matrix.size(); ++i) {
  //      temp_select.push_back(select_matrix[i]);
  //      if(temp_select.size() == 3) {
  //          select.push_back(temp_select);
  //          temp_select.clear();
  //      }
  //  }

    std::vector<Fish> outputPop = selectPopulation(Pop,
                                                   select_matrix,
                                                   selection,
                                                   pop_size,
                                                   total_runtime,
                                                   morgan);

    return List::create( Named("population") = convert_to_list(outputPop) );
}





double calculate_fitness_markers(const Fish& focal,
                                 const NumericMatrix& select) {

    double fitness = 2.0;

    
    int focal_marker = 0;
    double pos = select(focal_marker, 0);
    double anc = select(focal_marker, 1);
    double s = select(focal_marker, 2);

    for(auto it = (focal.chromosome1.begin()+1); it != focal.chromosome1.end(); ++it) {
        if((*it).pos > pos) {
            if((*(it-1)).right == anc) fitness += s;
            focal_marker++;
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
            pos = select(focal_marker, 0);
            anc = select(focal_marker, 1);
            s = select(focal_marker, 2);

        }
    }

    return(fitness / 2.0);
}


std::vector< Fish > selectPopulation_vector(const std::vector< Fish>& sourcePop,
                                            const NumericMatrix& select,
                                            int pop_size,
                                            int total_runtime,
                                            double morgan) {

    double expected_max_fitness = 1.0;

    Rcout << "Applying matrix with:\n";
    Rcout << "Location" << "\t" << "Ancestor" << "\t" << "Selection coefficient" << "\n";
    for(int i = 0; i < select.nrow(); ++i) {
        Rcout << select(i, 0) << "\t" << select(i, 1) << "\t" << select(i, 2) << "\n";
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

    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";

    for(int t = 0; t < total_runtime; ++t) {

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

        if(t % updateFreq == 0) {
            Rcout << "**";
        }

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
                                                 double morgan)
{
    std::vector< Fish > Pop;
    for(int i = 0; i < pop_size; ++i) {
        Fish p1 = Fish( random_number( number_of_founders ) );
        Fish p2 = Fish( random_number( number_of_founders ) );

        Pop.push_back(mate(p1,p2, morgan));
    }

    std::vector<Fish> outputPop = selectPopulation_vector(Pop,
                                                          select,
                                                          pop_size,
                                                          total_runtime,
                                                          morgan);
    
    return List::create( Named("population") = convert_to_list(outputPop) );
}

// [[Rcpp::export]]
List select_population_markers_cpp(Rcpp::NumericVector v1,
                           Rcpp::NumericMatrix selectM,
                           int population_size,
                           int run_time,
                           double morgan) {

    std::vector< Fish > Pop = convert_NumericVector_to_fishVector(v1);

    std::vector<Fish> outputPop = selectPopulation_vector(Pop,
                                                   selectM,
                                                   population_size,
                                                   run_time,
                                                   morgan);

    return List::create( Named("population") = convert_to_list(outputPop) );
}
