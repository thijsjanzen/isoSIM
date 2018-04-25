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
#include "main.h"

#include "random_functions.h"

//#include "randomc.h"


#include <Rcpp.h>
using namespace Rcpp;

bool verify_individual_cpp(const Fish& Nemo) {
    for(int i = 0; i < Nemo.chromosome1.size(); ++i) {
        if(Nemo.chromosome1[i].right >  1000 |
           Nemo.chromosome1[i].right < -1000) {
            return false;
        }
    }

    for(int i = 0; i < Nemo.chromosome2.size(); ++i) {
        if(Nemo.chromosome2[i].right >  1000 |
           Nemo.chromosome2[i].right < -1000) {
            return false;
        }
    }

    return true;
}


bool verify_pop_cpp(const std::vector< Fish >& pop) {
    for(auto it = pop.begin(); it != pop.end(); ++it) {
        if(!verify_individual_cpp((*it))) {
            return false;
        }
    }

    return true;
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



/***********************************************************************
 Start of simulation code
 ***********************************************************************/

std::vector< Fish > simulate(const std::vector< Fish >& input_pop,
                             int popSize,
                             int maxTime,
                             double Morgan,
                             bool progress_bar,
                             bool track_junctions,
                             std::vector<double>& junctions)
{
    if(progress_bar) {
        Rcout << "0--------25--------50--------75--------100\n";
        Rcout << "*";
    }

    std::vector< Fish > Pop = input_pop;

    int updateFreq = maxTime / 20;
    if(updateFreq < 1) updateFreq = 1;

    for(int t = 0; t < maxTime; ++t) {

        if(track_junctions) junctions.push_back(calc_mean_junctions(Pop));

        std::vector<Fish> newGeneration;

        for(int i = 0; i < popSize; ++i)  {
            int index1 = random_number(Pop.size());
            while(index1 >= Pop.size()) index1 = random_number(Pop.size());

            int index2 = random_number(Pop.size());
            while(index2 >= Pop.size() || index1 == index2) index2 = random_number(Pop.size());

            Fish kid = mate(Pop[index1], Pop[index2], Morgan);
            
           // if(verify_individual_cpp(kid)) {
                newGeneration.push_back(kid);
           // } else {
            //    Rcout << "\n After " << t << " generations, verify individual failed\n"; R_FlushConsole();
            //}
        }

        Pop = newGeneration;
        newGeneration.clear();

        if(t % updateFreq == 0 && progress_bar) {
            Rcout << "**";
        }

        if(is_fixed(Pop)) {
            Rcout << "\n After " << t << " generations, the population has become completely homozygous and fixed\n"; R_FlushConsole();
            return(Pop);
        }

      //  if(!verify_pop_cpp(Pop)) {
      //      Rcout << "\n After " << t << " generations, verify population failed\n"; R_FlushConsole();
      //  }

        Rcpp::checkUserInterrupt();
    }
   if(progress_bar) Rcout << "\n";
   // Rcout << "No fixation before maximum run time was reached\n"; R_FlushConsole();
    return(Pop);
}

std::vector<Fish> create_line(const std::vector< Fish >& founders,
                             int popSize,
                             int maxTime,
                             double Morgan,
                             bool progress_bar)
{
    std::vector<Fish> Pop;

    // if(!verify_pop_cpp(founders)) {
    //    Rcout << "Verify founders in create_line failed\n"; R_FlushConsole();
    // }

    if(founders.size() == 2) {
        for(int i = 0; i < popSize; ++i) {
            Fish temp = mate( founders[0], founders[1], Morgan);
            Pop.push_back( temp);
        }
    } else {
        for(int i = 0; i < popSize; ++i) {
            int index1 = random_number(founders.size());
            while(index1 >= (int)founders.size()) index1 = random_number(founders.size());

            int index2 = random_number(founders.size());
            while(index2 >= (int)founders.size()) index2 = random_number(founders.size());

            while(index1 == index2) index2 = random_number(founders.size());

            Fish temp = mate( founders[index1], founders[index2], Morgan);
            Pop.push_back( temp);
        }

    }

  //  if(!verify_pop_cpp(Pop)) {
  //      Rcout << "Creation in create_line failed\n"; R_FlushConsole();
  //  }
    std::vector<double> junctions;
    Pop = simulate(Pop, popSize, maxTime, Morgan, progress_bar, FALSE, junctions);
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
            kid = mate(main_pop[index1], immigrant_pop[index2], Morgan);
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
                              double m,
                              bool progress_bar)
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
        Fish parent1 = parents1[ random_number( parents1.size() ) ];
        Fish parent2 = parents1[ random_number( parents1.size() ) ];

        Pop1.push_back(mate(parent1, parent2, Morgan));
    }

    for(int i = 0; i < pop_size; ++i) {
        Fish parent1 = parents2[ random_number( parents2.size() ) ];
        Fish parent2 = parents2[ random_number( parents2.size() ) ];

        Pop2.push_back(mate(parent1, parent2, Morgan));
    }

    int updateFreq = max_time / 20;
    if(updateFreq < 1) updateFreq = 1;

    if(progress_bar) {
        Rcout << "0--------25--------50--------75--------100\n";
        Rcout << "*";
    }

    for(int t = 0; t < max_time; ++t) {

        std::vector<Fish> newGeneration1 = update_pop(Pop1, Pop2, Morgan, m);
        std::vector<Fish> newGeneration2 = update_pop(Pop2, Pop1, Morgan, m);

        Pop1 = newGeneration1;
        newGeneration1.clear();

        Pop2 = newGeneration2;
        newGeneration2.clear();

        if(t % updateFreq == 0 && progress_bar) {
            Rcout << "**";
        }
        Rcpp::checkUserInterrupt();
    }
    
    p1 = Pop1;
    p2 = Pop2;
    
    return;
}


// [[Rcpp::export]]
List create_population_cpp(int pop_size,
                       int number_of_founders,
                       int total_runtime,
                       double morgan,
                       bool progress_bar,
                       bool track_junctions)
{
    std::vector<Fish> Pop;

    for(int i = 0; i < pop_size; ++i) {
        Fish p1 = Fish(random_number(number_of_founders));
        Fish p2 = Fish(random_number(number_of_founders));

        Pop.push_back(mate(p1,p2, morgan));
    }
    std::vector<double> junctions;

    Pop = simulate(Pop, pop_size, total_runtime, morgan, progress_bar,
                   track_junctions, junctions);

    if(track_junctions) {
        return List::create( Named("population") = convert_to_list(Pop) ,
                             Named("junctions") = junctions);
    }

    return List::create( Named("population") = convert_to_list(Pop) );

}

// [[Rcpp::export]]
List create_isofemale_line_cpp(NumericVector v,
                       int pop_size,
                       int total_runtime,
                       double morgan,
                       bool progress_bar)
{
    std::vector< Fish > founders = convert_NumericVector_to_fishVector(v);

    std::vector<Fish> Pop = create_line(founders, pop_size,
                                        total_runtime, morgan, progress_bar);

    return List::create( Named("population") = convert_to_list(Pop) );
}

// [[Rcpp::export]]
List create_two_populations_migration_cpp(int pop_size,
                                          int number_of_founders,
                                          int total_runtime,
                                          double morgan,
                                          double migration,
                                          bool progress_bar) {
    std::vector< Fish > Pop1;
    std::vector< Fish > Pop2;

    create_two_pop_migration(Pop1, Pop2, number_of_founders,
                             pop_size, total_runtime, morgan, migration, progress_bar);

    return List::create( Named("population_1")  = convert_to_list(Pop1),
                         Named("population_2")  = convert_to_list(Pop2)
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

bool matching_chromosomes(const std::vector< junction >& v1,
                          const std::vector< junction >& v2)
{
    if(v1.size() != v2.size()) {
        return false;
    }
    for(int i = 0; i < v1.size(); ++i) {
        if(v1[i] != v2[i]) {
            return false;
        }
    }
    return true;
}



bool is_fixed(const std::vector< Fish >& v) {

    if(!matching_chromosomes(v[0].chromosome1, v[0].chromosome2)) {
        return false;
    }

    for(auto it = v.begin(); it != v.end(); ++it) {
        if(!matching_chromosomes((*it).chromosome1, v[0].chromosome1)) {
            return false;
        }
        if(!matching_chromosomes((*it).chromosome1, (*it).chromosome2)) {
            return false;
        }
    }
    return true;
}

// [[Rcpp::export]]
void test_fish_functions() {
    Fish test_fish;

    junction temp;
    junction temp2(0.5, 0);
    junction temp3(0.5, 0);
    if(temp2 == temp3) {
        temp = temp2;
    }

    bool temp400 = (temp2 != temp3);
    Rcout << temp400 << "\t" << "this is only for testing\n";

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
    bool b = (test_fish == test_fish3);
    Rcout << b << "\t" << "this is only for testing\n";

    std::vector< junction > chrom;
    chrom.push_back(temp);

    Fish test_fish4(chrom, chrom);

    std::vector< Fish > pop;
    pop.push_back(test_fish);
    pop.push_back(test_fish2);
    pop.push_back(test_fish3);
    pop.push_back(test_fish4);

    verify_pop_cpp(pop);


    return;
}
