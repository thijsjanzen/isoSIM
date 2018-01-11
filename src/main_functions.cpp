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


void flush_console() {
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();

}



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
            Rcout << "ERROR in calculating fitness\n"; flush_console();
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



// [[Rcpp::export]]
List select_population_cpp(Rcpp::NumericVector v1,
                       Rcpp::NumericVector selectM,
                       double s,
                       int population_size,
                       int run_time,
                       double morgan,
                       int seed,
                       bool writeToFile) {

    set_seed(seed);
    std::vector< Fish > Pop;

    std::vector<double> v = Rcpp::as<std::vector<double> >(v1);

    Fish temp;
    int indic_chrom = 1;
    bool add_indiv = false;

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
            Pop.push_back(temp);
            add_indiv = false;
            indic_chrom = 1;
            temp.chromosome1.clear();
            temp.chromosome2.clear();
           // Rcout << v.size() << "\t" << v1.size() << "\t" << Pop.size() << "\n"; flush_console();
        }
    }

    std::vector< std::vector< double > > select;
    std::vector<double> temp_select;
    for(int i = 0; i < selectM.size(); ++i) {
        temp_select.push_back(selectM[i]);
        if(temp_select.size() == 3) {
            select.push_back(temp_select);
            temp_select.clear();
        }
    }

    //for(int i = 0; i < select.size(); ++i){
    //    for(int j = 0; j < 3; ++j) {
    //        Rcout << select[i][j] << " "; flush_console();
    //    }
    //    Rcout << "\n"; flush_console();
   // }



    std::vector<Fish> outputPop = selectPopulation(Pop,
                                              select,
                                              s,
                                              population_size,
                                              run_time,
                                              morgan);


    if(writeToFile) {
       writePoptoFile(outputPop, "population_1.pop");
    }

    return List::create( Named("population") = createPopVector(outputPop) );
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

    std::vector<Fish> Pop = create_line(founders, pop_size,
                                        total_runtime, morgan);

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


NumericMatrix allele_spectrum(const std::vector<Fish>& v,
                              double step_size,
                              int numAncestors) {

    int numSteps = 1.0 / step_size;

    NumericMatrix spectrum(numSteps * numAncestors, 3);

    double left = 0.0;
    double right = step_size;

    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";

    int updateFreq = numSteps / 20;
    if(updateFreq < 1) updateFreq = 1;


    for(int i = 0; i < numSteps; ++i) {

        if(i % updateFreq == 0) {
            Rcout << "**";
        }


        for(int ancestor = 0; ancestor < numAncestors; ++ancestor) {
            double local_freq = 0.0;
            for(auto it = v.begin(); it != v.end(); ++it) {

                double a = assess_match( (*it).chromosome1, left, right, ancestor);
                double b = assess_match( (*it).chromosome2, left, right, ancestor);

                double freq =  (a+b);

                local_freq += freq;
            }
            int index = ancestor * numSteps + i;
            spectrum(index, 0) = right;
            spectrum(index, 1) = ancestor + 1;
            spectrum(index, 2) = local_freq / (2.0 * v.size());
        }
        left = right;
        right += step_size;
    }
    
    return spectrum;
}

// [[Rcpp::export]]
NumericMatrix calculate_allele_spectrum_cpp(NumericVector v1,
                                            int numFounders,
                                            double step_size)
{
    std::vector< Fish > Pop;

    std::vector<double> v = Rcpp::as<std::vector<double> >(v1);

    Fish temp;
    int indic_chrom = 1;
    bool add_indiv = false;

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
            Pop.push_back(temp);
            add_indiv = false;
            indic_chrom = 1;
            temp.chromosome1.clear();
            temp.chromosome2.clear();
        }
    }

    NumericMatrix output = allele_spectrum(Pop, step_size, numFounders);

    return output;
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







