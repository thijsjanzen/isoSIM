//
//  main.cpp
//  SMC_freqs
//
//  Created by Thijs Janzen on 29/01/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <numeric>
#include <cmath>

#include <vector>
#include <algorithm>

#include "randomc.h"
#include "Fish.h"
#include "Output.h"

#include <Rcpp.h>
using namespace Rcpp;





/***********************************************************************
 Start of simulation code
 ***********************************************************************/


Output doSimulation(int popSize,
                    double initRatio,
                    int maxTime,
                    double numRecombinations,
                    int numberOfMarkers
                    )    {

    Output O;
    std::vector<Fish> Pop;
    std::vector<double> markers;
    if(numberOfMarkers > 0) {
        for(int i = 0; i < numberOfMarkers; ) {
            double pos = uniform();
            if(pos > 0 && pos < 1.0) {
                ++i;
                markers.push_back(pos);
            }
        }
        std::sort(markers.begin(), markers.end());
    }


    Fish parent1 = Fish(0);
    Fish parent2 = Fish(1);

    for(int i = 0; i < popSize; ++i) {
        Fish p1 = parent1;
        Fish p2 = parent1;

        if(uniform() < initRatio) {
            p1 = parent2;
        }
        if(uniform() < initRatio) {
            p2 = parent2;
        }


        Pop.push_back(mate(p1,p2, numRecombinations));
    }

    for(int t = 0; t < maxTime; ++t) {
        O.update(Pop);
        if(numberOfMarkers > 0) O.detectNumJunctions(Pop, markers);

        std::vector<Fish> newGeneration;

        for(int i = 0; i < popSize; ++i)  {
            int index1 = random_number(popSize);
            int index2 = random_number(popSize);

            Fish kid = mate(Pop[index1], Pop[index2], numRecombinations);

            newGeneration.push_back(kid);
        }

        Pop = newGeneration;
        newGeneration.clear();
    }

    return O;
}

Output continue_simulation(std::vector< Fish > Pop,
                           int max_time,
                           double morgan,
                           int numberOfMarkers)
{
    Output O;
    int popSize = (int)Pop.size();

    std::vector<double> markers;
    if(numberOfMarkers > 0) {
        for(int i = 0; i < numberOfMarkers; ) {
            double pos = uniform();
            if(pos > 0 && pos < 1.0) {
                ++i;
                markers.push_back(pos);
            }
        }
        std::sort(markers.begin(), markers.end());
    }



    for(int t = 0; t < max_time; ++t) {
        O.update(Pop);
        if(numberOfMarkers > 0) O.detectNumJunctions(Pop, markers);

        std::vector<Fish> newGeneration;

        for(int i = 0; i < popSize; ++i)  {
            int index1 = random_number(popSize);
            int index2 = random_number(popSize);

            Fish kid = mate(Pop[index1], Pop[index2], morgan);

            newGeneration.push_back(kid);
        }

        Pop = newGeneration;
        newGeneration.clear();
    }
    
    return O;
}

std::vector<Fish> createPopulation(int popSize,
                      int numFounders,
                      int maxTime,
                      double Morgan,
                      double overlap,
                      double populationIndicator)
{
    std::vector<Fish> Pop;
    std::vector<Fish> parents;
    int i = 0;

    if(populationIndicator == 0) {
        for(int i = 0; i < 1000; ++i) {
            parents.push_back(Fish(random_number(numFounders)));
        }
    }

    if(populationIndicator == 1) {
        //second population, we will have to take
        //into account overlap
        //i = numFounders - overlap * numFounders;
        //numFounders += i;
        for(int i = 0; i < 1000; ++i) {
            int index = numFounders + random_number(numFounders);
            if(uniform() < overlap) {
                index = random_number(numFounders);
            }
            parents.push_back(Fish(index));
        }
    }



    for(int i = 0; i < popSize; ++i) {
        Fish p1 = parents[ random_number( parents.size() ) ];
        Fish p2 = parents[ random_number( parents.size() ) ];

        Pop.push_back(mate(p1,p2, Morgan));
    }

    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";

    int updateFreq = maxTime / 20;
    if(updateFreq < 1) updateFreq = 1;

    for(int t = 0; t < maxTime; ++t) {

        std::vector<Fish> newGeneration;

        for(int i = 0; i < popSize; ++i)  {
            int index1 = random_number(popSize);
            int index2 = random_number(popSize);

            Fish kid = mate(Pop[index1], Pop[index2], Morgan);

            newGeneration.push_back(kid);
        }

        Pop = newGeneration;
        newGeneration.clear();

        if(t % updateFreq == 0) {
            Rcout << "**";
        }
    }
    Rcout << "\n";
    return(Pop);
}


bool is_fixed(const std::vector< Fish >& v) {

    if(v[0].chromosome1 != v[0].chromosome2) return false;

    for(auto it = v.begin(); it != v.end(); ++it) {
        if((*it).chromosome1 != v[0].chromosome1) return false;
        if((*it).chromosome2 != v[0].chromosome2) return false;
    }

    return true;
}


std::vector<Fish> create_line(const std::vector< Fish >& founders,
                             int popSize,
                             int maxTime,
                             double Morgan,
                             double overlap,
                             double populationIndicator)
{
    std::vector<Fish> Pop;


    for(int i = 0; i < popSize; ++i) {
        Pop.push_back(mate(founders[0], founders[1], Morgan));
    }

    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";

    int updateFreq = maxTime / 20;
    if(updateFreq < 1) updateFreq = 1;

    for(int t = 0; t < maxTime; ++t) {

        std::vector<Fish> newGeneration;

        for(int i = 0; i < popSize; ++i)  {
            int index1 = random_number(popSize);
            int index2 = random_number(popSize);

            Fish kid = mate(Pop[index1], Pop[index2], Morgan);

            newGeneration.push_back(kid);
        }

        Pop = newGeneration;
        newGeneration.clear();
        
        if(t % updateFreq == 0) {
            Rcout << "**";
        }

        if(is_fixed(Pop)) {
            Rcout << "\nPreliminary exit because the population is already completely homozygous\n";
            break;
        }
    }
    Rcout << "\n";
    return(Pop);
}



std::vector<double> createPopVector(const std::vector< Fish >& v) {
    std::vector<double> output;
    for(auto it = v.begin(); it != v.end(); ++it) {

        for(auto i = (*it).chromosome1.begin(); i != (*it).chromosome1.end(); ++i) {
            output.push_back((*i).pos);
            output.push_back((*i).right);
        }

        for(auto j = (*it).chromosome2.begin(); j != (*it).chromosome2.end(); ++j) {
            output.push_back((*j).pos);
            output.push_back((*j).right);
        }
    }
    return(output);
}





/*
 //c++ code for independent compiling
 int main(int argc, const char * argv[]) {

 macstart(argv);

	set_seed(seed);

 Output O = doSimulation(popSize, initRatio, maxTime, numRecombinations);

	return 0;
 }
 */


// [[Rcpp::export]]
Numeric calc_heterozygosity(NumericVector v) {
    std::vector< Fish > pop;
    Fish temp;
    int indic_chrom = 1
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
            pop.push_back(temp);
            add_indiv = false;
            indic_chrom = 1;
            temp.chromosome1.clear();
            temp.chromosome2.clear();
        }
    }

    double heterozygosity = calculate_heterozygosity(pop);
    return(heterozygosity);
}




// [[Rcpp::export]]
List simulate_from_population(std::string file_name,
                              int total_runtime,
                              double morgan,
                              int number_of_markers,
                              int seed)
{
    set_seed(seed);
    std::vector< Fish > Pop;
    readPopfromFile(Pop, file_name);

    Output O = continue_simulation(Pop, total_runtime, number_of_markers, morgan);
    return List::create(Named("avgJunctions") = O.avgJunct,
                        Named("detectedJunctions") = O.avg_detected_Junctions);
}

// [[Rcpp::export]]
List create_population(int pop_size,
                       int number_of_founders,
                       int total_runtime,
                       double morgan,
                       int seed,
                       bool writeToFile)
{
    set_seed(seed);
    std::vector<Fish> Pop = createPopulation(pop_size, number_of_founders,
                                             total_runtime, morgan, 0.0, 0);
    if(writeToFile) {
        writePoptoFile(Pop, "population_1.pop");
    }

    return List::create( Named("population") = createPopVector(Pop) );
}

// [[Rcpp::export]]
List create_femaleLine(NumericVector indiv,
                       int pop_size,
                       int total_runtime,
                       double morgan,
                       int seed)
{
    set_seed(seed);

    std::vector< Fish > founders;

    for(int f = 0; f < 2; ++f) {
        Fish founder;
        int chromIndicator = 0;
        for(int i = 0; i < indiv.size();  i+=2) {
            junction temp;
            temp.pos = indiv[i];
            temp.right = indiv[i+1];

            if(chromIndicator == 0) {
                founder.chromosome1.push_back(temp);
            }
            if(chromIndicator == 1) {
                founder.chromosome2.push_back(temp);
            }

            if(temp.right == -1) {
                chromIndicator++;
            }
            if(chromIndicator == 2) {
                founders.push_back(founder);
                founder.chromosome1.clear();
                founder.chromosome2.clear();
                chromIndicator = 0;
            }
        }
    }

    std::vector<Fish> Pop = create_line(founders, pop_size,
                                        total_runtime, morgan, 0.0, 0);

    return List::create( Named("population") = createPopVector(Pop) );
}




// [[Rcpp::export]]
List create_two_populations(int pop_size,
                            int number_of_founders,
                            int total_runtime,
                            double morgan,
                            int seed,
                            double overlap,
                            bool writeToFile) {
    set_seed(seed);
    std::vector<Fish> Pop1 = createPopulation(pop_size, number_of_founders,
                                              total_runtime, morgan, overlap, 0);
    if(writeToFile) {
        writePoptoFile(Pop1, "population_1.pop");
    }

    std::vector<Fish> Pop2 = createPopulation(pop_size, number_of_founders,
                                              total_runtime, morgan, overlap, 1);
    if(writeToFile) {
        writePoptoFile(Pop2, "population_2.pop");
    }

    return List::create( Named("population_1") = createPopVector(Pop1),
                         Named("population_2") = createPopVector(Pop2)
                       );
}

// [[Rcpp::export]]
List sim_inf_chrom(int pop_size,
                   double initial_heterozygosity,
                   int total_runtime,
                   double morgan,
                   int markers,
                   int seed) {
    set_seed(seed);
    double p = 0.5 * (1 - sqrt(1 - 2 * initial_heterozygosity));

    Output O = doSimulation(pop_size, p, total_runtime, morgan, markers);
    return List::create(Named("avgJunctions") = O.avgJunct,
                        Named("detectedJunctions") = O.avg_detected_Junctions);
}
