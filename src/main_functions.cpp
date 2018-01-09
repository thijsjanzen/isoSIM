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



std::vector<Fish> create_line(const std::vector< Fish >& founders,
                             int popSize,
                             int maxTime,
                             double Morgan)
{
    std::vector<Fish> Pop;

    for(int i = 0; i < popSize; ++i) {
        Fish temp = mate( founders[0], founders[1], Morgan);
        std::vector< Fish > to_print;
        to_print.push_back(temp);
       // Rcout << i << "\t";
     //   output_alleles(to_print);
        Pop.push_back( temp);
    }
   // output_alleles(Pop);

    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";

    int updateFreq = maxTime / 20;
    if(updateFreq < 1) updateFreq = 1;

    for(int t = 0; t < maxTime; ++t) {

        std::vector<Fish> newGeneration;

        for(int i = 0; i < popSize; ++i)  {
            int index1 = random_number(popSize);
            int index2 = random_number(popSize);

            while(index2 == index1) index2 = random_number(popSize);

            Fish kid = mate(Pop[index1], Pop[index2], Morgan);

            newGeneration.push_back(kid);
        }

        Pop = newGeneration;
        newGeneration.clear();
        
        if(t % updateFreq == 0) {
            Rcout << "**";
        }

        if(is_fixed(Pop)) {
            //Rcout << "\nPreliminary exit because the population is already completely homozygous\n";
            Rcout << "\n After " << t << " generations, the population has become completely homozygous and fixed\n iso-females are ready!\n";
            return(Pop);
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

void update_allele_freq(const Fish& indiv, std::vector<double>& allele_freq) {


    for(auto it = indiv.chromosome1.begin(); it != indiv.chromosome1.end(); ++it) {
        double left = (*it).pos;
        double right = 1.0;
        if( ((it+1)) != indiv.chromosome1.end()) {
            right = (*(it+1)).pos;
        }

        allele_freq[ (*it).right] += (right - left);
    }

    for(auto it = indiv.chromosome2.begin(); it != indiv.chromosome2.end(); ++it) {
        double left = (*it).pos;
        double right = 1.0;
        if( ((it+1)) != indiv.chromosome2.end()) {
            right = (*(it+1)).pos;
        }

        allele_freq[ (*it).right] += (right - left);
    }
    return;
}



std::vector< double > calculate_mean_allelefreq(const std::vector< Fish >& pop,
                                                       int number_of_founders) {

    std::vector< double > allele_freq(number_of_founders * 2, 0.0);

    for(auto it = pop.begin(); it != pop.end(); ++it) {
        update_allele_freq((*it), allele_freq);
    }

    for(int i = 0; i < allele_freq.size(); ++i) {
        allele_freq[i] = 1.0 * allele_freq[i] / (2 * pop.size());
    }
    return(allele_freq);
}

int draw_prop_fitness(const std::vector<double> fitness,
                      double maxFitness) {
    for(int i = 0; i < 1e6; ++i) {
        int index = random_number(fitness.size());
        if(uniform() < fitness[index] / maxFitness) {
            return index;
        }
    }
    Rcout << "\nERROR! ERROR! Couldn't pick proportional to fitness\n";
    return -1;
}

double assess_match(const std::vector<junction> chrom,
                    double start,
                    double end,
                    int ancestor) {

    std::vector< junction > block;
    bool first_one = true;
    for(int i = 0; i < chrom.size(); ++i) {
        if(chrom[i].pos > start) {

            if(chrom[i].pos > end) {
                if(first_one) {
                    block.push_back(chrom[i]);
                    first_one = false;
                } else {
                    break;
                }
            }
            block.push_back(chrom[i]);
        }
    }

    if(block.size() == 1) {
        if(block[0].left == ancestor) return 1.0;

        return 0.0;
    }

    double match = 0.0;
    if(block[0].left == ancestor) {
        match += (block[0].pos - start);
    }
    for(int i = 1; i < block.size(); ++i) {
        double local_left = block[i-1].pos;
        double local_right = block[i].pos;
        if(local_right > end) local_right = end;

        if(block[i].left == ancestor) match += (local_right - local_left);
    }
    match *= 1.0 / (end - start);
    return(match);
}



double calculate_fitness(const Fish& focal,
                         const std::vector< std::vector< double > >& select,
                         double s) {

    double fitness = 2.0 * select[0][0];

    for(int i = 0; i < select.size(); ++i) {
        double start = select[i][0];
        double end = select[i][1];
        int ancestor = select[i][2];

        double a1 = assess_match(focal.chromosome1, start, end, ancestor);
        double a2 = assess_match(focal.chromosome2, start, end, ancestor);

        //fitness += (end - start) * a1 * (1+s) + (end - start) * (1 - a1) * (1);
        //fitness += (end - start) * a2 * (1+s) + (end - start) * (1 - a2) * (1);
        //fitness += (end - start) * (1 + s * a1);
        //fitness += (end - start) * (1 + s * a2);
        fitness += (end - start) * (2 + s * (a1 + a2));

        double nextStart = 1.0;
        if((i+1) < select.size()) {
            nextStart = select[i+1][0];
        }
        fitness += 2.0 * (nextStart - end);
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
    double maxFitness = -1;
    for(auto it = Pop.begin(); it != Pop.end(); ++it){
        double fit = calculate_fitness((*it), select, s);
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
    }
    return(Pop);
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
double calc_heterozygosity_cpp(NumericVector v) {
    std::vector< Fish > pop;
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

void flush_console() {
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();

}


// [[Rcpp::export]]
List select_population_cpp(Rcpp::NumericVector v1,
                       Rcpp::NumericVector selectM,
                       double s,
                       int population_size,
                       int run_time,
                       double morgan,
                       int seed,
                       bool writeToFile) {

    Rcout << "CPP: start\n"; flush_console();
    set_seed(seed);
    std::vector< Fish > Pop;

    usleep(100);

    Rcout << "CPP: converting population\n"; flush_console();

    std::vector<double> v = Rcpp::as<std::vector<double> >(v1);
    Rcout << "CPP: RCPP vector conversion done\n"; flush_console();

    Fish temp;
    int indic_chrom = 1;
    bool add_indiv = false;

    usleep(100);

    if(1 == 2) {


    for(int i = 0; i < (v.size() - 1); i += 2) {
        junction temp_j;
        temp_j.pos = v[i];
        if(i+1 > v.size()) break;
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
            Rcout << v.size() << "\t" << v1.size() << "\t" << pop.size() << "\n"; flush_console();
        }
    }

    Rcout << "CPP: loaded individuals\n"; flush_console();

    usleep(100);

    std::vector<double> selectMatrix = Rcpp::as<std::vector<double>>(selectM);
    Rcout << "CPP: converting select\n";
    std::vector<std::vector<double>> select;
    for(int i = 0; i < selectMatrix.size(); ++i) {
        std::vector<double> temp;
        temp.push_back(selectMatrix[i]);
        if(temp.size() == 3) {
            select.push_back(temp);
            temp.clear();
        }
    }


    Rcout << "CPP: starting simulation\n"; flush_console();
    std::vector<Fish> Pop = selectPopulation( pop,
                                              select,
                                              s,
                                              population_size,
                                              run_time,
                                              morgan);
    if(writeToFile) {
        writePoptoFile(Pop, "population_1.pop");
    }
    }

    Rcout << "CPP: Done\n"; flush_console();
    return List::create( Named("population") = createPopVector(Pop) );
}


// [[Rcpp::export]]
List calculate_summaryStats(NumericVector v,
                            int number_of_founders) {
    std::vector< Fish > pop;
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
            pop.push_back(temp);
            add_indiv = false;
            indic_chrom = 1;
            temp.chromosome1.clear();
            temp.chromosome2.clear();
        }
    }

    double heterozygosity = calculate_heterozygosity(pop);

    std::vector < double > allele_freq = calculate_mean_allelefreq(pop,
                                                                   number_of_founders);

    return List::create( Named("Hst") = heterozygosity,
                         Named("freq_pop") = allele_freq );
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
List create_femaleLine(NumericVector v,
                       int pop_size,
                       int total_runtime,
                       double morgan,
                       int seed)
{
    set_seed(seed);

    std::vector< Fish > founders;

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
            founders.push_back(temp);
            add_indiv = false;
            indic_chrom = 1;
            temp.chromosome1.clear();
            temp.chromosome2.clear();
        }
    }

    //Rcout << "Founders: ";
    //output_alleles(founders);


    std::vector<Fish> Pop = create_line(founders, pop_size,
                                        total_runtime, morgan);


    //Rcout << "Final: ";
    //output_alleles(Pop);

    return List::create( Named("population") = createPopVector(Pop) );
}




// [[Rcpp::export]]
List create_two_populations_cpp(int pop_size,
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
