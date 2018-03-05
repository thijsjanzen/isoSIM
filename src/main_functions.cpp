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

#include <unistd.h> //for sleep

#include "Fish.h"
#include "Output.h"
#include "random_functions.h"

#include <Rcpp.h>
using namespace Rcpp;


/***********************************************************************
 Start of simulation code
 ***********************************************************************/

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
        Rcpp::checkUserInterrupt();
    }
    
    return O;
}

std::vector< Fish > simulate(const std::vector< Fish >& input_pop,
                             int popSize,
                             int maxTime,
                             double Morgan)
{
    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";

    std::vector< Fish > Pop = input_pop;

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
            Rcout << "\n After " << t << " generations, the population has become completely homozygous and fixed\n iso-females are ready!\n";
            return(Pop);
        }

        Rcpp::checkUserInterrupt();
    }
    Rcout << "\n";
    return(Pop);
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

    Pop = simulate(Pop, popSize, maxTime, Morgan);

    return(Pop);
}

std::vector<Fish> create_line(const std::vector< Fish >& founders,
                             int popSize,
                             int maxTime,
                             double Morgan)
{
    std::vector<Fish> Pop;

    if(founders.size() == 2) {
        for(int i = 0; i < popSize; ++i) {
            Fish temp = mate( founders[0], founders[1], Morgan);
            std::vector< Fish > to_print;
            to_print.push_back(temp);
           Pop.push_back( temp);
        }
    } else {
        for(int i = 0; i < popSize; ++i) {
            int index1 = random_number(founders.size());
            int index2 = random_number(founders.size());
            while(index1 == index2) index2 = random_number(founders.size());

            Fish temp = mate( founders[index1], founders[index2], Morgan);
            std::vector< Fish > to_print;
            to_print.push_back(temp);
            Pop.push_back( temp);
        }

    }
    Pop = simulate(Pop, popSize, maxTime, Morgan);
    return(Pop);
}

std::vector<Fish> update_pop(const std::vector< Fish>& main_pop,
                             const std::vector< Fish>& immigrant_pop,
                             double Morgan,
                             double m) {

    std::vector< Fish > offspring;
    int pop_size = main_pop.size();

    for(int i = 0; i < pop_size; ++i)  {

        int index1 = random_number(pop_size);
        int index2 = random_number(pop_size);
        while(index2 == index1) index2 = random_number(pop_size);

        Fish kid;

        if(uniform() < m) {
            kid = mate(immigrant_pop[index1], immigrant_pop[index2], Morgan);
        } else {
            kid = mate(main_pop[index1], main_pop[index2], Morgan);
        }

        offspring.push_back(kid);
    }
    return offspring;
}

void create_two_pop_migration( std::vector< Fish >& p1,
                              std::vector< Fish >& p2,
                              int num_ancestors_per_pop,
                              int pop_size,
                              int max_time,
                              double Morgan,
                              double m)
{
    std::vector<Fish> Pop1;
    std::vector<Fish> Pop2;

    std::vector<Fish> parents1;
    std::vector<Fish> parents2;

    for(int i = 0; i < num_ancestors_per_pop; ++i) {
        parents1.push_back(Fish(i));
        parents2.push_back(Fish(i + num_ancestors_per_pop));

    }

    for(int i = 0; i < pop_size; ++i) {
        Fish parent1 = parents1[ random_number( num_ancestors_per_pop ) ];
        Fish parent2 = parents1[ random_number( num_ancestors_per_pop ) ];

        Pop1.push_back(mate(parent1, parent2, Morgan));
    }

    for(int i = 0; i < pop_size; ++i) {
        Fish parent1 = parents2[ random_number( num_ancestors_per_pop ) ];
        Fish parent2 = parents2[ random_number( num_ancestors_per_pop ) ];

        Pop2.push_back(mate(parent1, parent2, Morgan));
    }

    int updateFreq = max_time / 20;
    if(updateFreq < 1) updateFreq = 1;

    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";

    for(int t = 0; t < max_time; ++t) {

        std::vector<Fish> newGeneration1 = update_pop(Pop1, Pop2, Morgan, m);
        std::vector<Fish> newGeneration2 = update_pop(Pop2, Pop1, Morgan, m);

        Pop1 = newGeneration1;
        newGeneration1.clear();

        Pop2 = newGeneration2;
        newGeneration2.clear();

        if(t % updateFreq == 0) {
            Rcout << "**";
        }
        Rcpp::checkUserInterrupt();
    }
    
    p1 = Pop1;
    p2 = Pop2;
    
    return;
}


// [[Rcpp::export]]
List simulate_from_population(std::string file_name,
                              int total_runtime,
                              double morgan,
                              int number_of_markers)
{
    std::vector< Fish > Pop;
    readPopfromFile(Pop, file_name);

    Output O = continue_simulation(Pop, total_runtime, number_of_markers, morgan);
    return List::create(Named("avgJunctions") = O.avgJunct,
                        Named("detectedJunctions") = O.avg_detected_Junctions);
}

// [[Rcpp::export]]
List create_population_cpp(int pop_size,
                       int number_of_founders,
                       int total_runtime,
                       double morgan,
                       bool writeToFile)
{
    std::vector<Fish> Pop = createPopulation(pop_size, number_of_founders,
                                             total_runtime, morgan, 0.0, 0);
    if(writeToFile) {
        writePoptoFile(Pop, "population_1.pop");
    }

    return List::create( Named("population") = createPopVector(Pop) );
}

// [[Rcpp::export]]
List create_isofemale_line_cpp(NumericVector v,
                       int pop_size,
                       int total_runtime,
                       double morgan)
{

    std::vector< Fish > founders = convert_NumericVector_to_fishVector(v);

    std::vector<Fish> Pop = create_line(founders, pop_size,
                                        total_runtime, morgan);

    return List::create( Named("population") = createPopVector(Pop) );
}


// [[Rcpp::export]]
List create_two_populations_cpp(int pop_size,
                            int number_of_founders,
                            int total_runtime,
                            double morgan,
                            double overlap,
                            bool writeToFile) {

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
List create_two_populations_migration_cpp(int pop_size,
                                          int number_of_founders,
                                          int total_runtime,
                                          double morgan,
                                          double migration,
                                          bool writeToFile) {
    std::vector< Fish > Pop1;
    std::vector< Fish > Pop2;

    create_two_pop_migration(Pop1, Pop2, number_of_founders, pop_size, total_runtime, morgan, migration);

    if(writeToFile) {
        writePoptoFile(Pop1, "population_1.pop");
        writePoptoFile(Pop2, "population_2.pop");
    }

    return List::create( Named("population_1") = createPopVector(Pop1),
                        Named("population_2") = createPopVector(Pop2)
                        );
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

void flush_console() {
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();
}

bool is_fixed(const std::vector< Fish >& v) {

    if(v[0].chromosome1 != v[0].chromosome2) return false;

    for(auto it = v.begin(); it != v.end(); ++it) {
        if((*it).chromosome1 != v[0].chromosome1) return false;
        if((*it).chromosome2 != v[0].chromosome2) return false;
    }

    return true;
}


// [[Rcpp::export]]
void test_fish_functions() {
    Fish test_fish;

    junction temp;
    junction temp2(0.5, -1, 0);
    junction temp3(0.5, -1, 0);
    if(temp2 == temp3) {
        temp = temp2;
    }
    if(temp2 != temp3) {
        temp = temp3;
    }

    junction temp4(temp);

    test_fish.chromosome1.push_back(temp);
    test_fish.chromosome1.push_back(temp2);
    test_fish.chromosome1.push_back(temp3);
    test_fish.chromosome1.push_back(temp4);

    Fish test_fish2 = test_fish;

    if(test_fish == test_fish2) {
        Rcout << "fishes are equal!\n";
    }

    Fish test_fish3(5);
    if(test_fish == test_fish3) {
        Rcout << "these fishes were supposed to be different!\n";
    }

    std::vector< junction > chrom;
    chrom.push_back(temp);

    Fish test_fish4(chrom, chrom);
    
    return;
}


/*
 void output_alleles(const std::vector<Fish>& v) {
 std::vector< int > alleles;
 for(int i = 0; i < v.size(); ++i) {
 for(int j = 0; j < v[i].chromosome1.size(); ++j) {
 int a = v[i].chromosome1[j].right;
 alleles.push_back(a);
 }
 for(int j = 0; j < v[i].chromosome2.size(); ++j) {
 int a = v[i].chromosome2[j].right;
 alleles.push_back(a);
 }
 }

 std::sort(alleles.begin(), alleles.end() );
 alleles.erase(std::unique(alleles.begin(), alleles.end()), alleles.end());

 for(int i = 0; i < alleles.size(); ++i) {
 Rcout << alleles[i] << "\t";
 }
 Rcout << "\n";

 return;
 }
 */







