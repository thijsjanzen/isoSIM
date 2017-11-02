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


bool isFixed(const std::vector<Fish>& V);
bool same(const Fish& A, const Fish& B);
template <typename T>
double calculateMean(const std::vector<T>& v);
template <typename T>
double calculateSD(const std::vector<T>& v);
double calcMeanFreq(const std::vector< Fish >& P, double L);


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

    //let's check!!!!
    writePoptoFile(Pop, "pop1.txt");

    std::vector<Fish> Pop2;
    readPopfromFile(Pop2, "pop1.txt");

    for(int i = 0; i < Pop.size(); ++i) {
        bool are_the_same = (Pop[i] == Pop2[i]);
        std::cout << are_the_same << "\n";
    }


    return O;
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
List sim_inf_chrom(int popSize,
                   double Hzero,
                   int maxTime,
                   double size_in_Morgan,
                   int markers,
                   int seed) {
    set_seed(seed);
    double p = 0.5 * (1 - sqrt(1 - 2*Hzero));

    Output O = doSimulation(popSize, p, maxTime, size_in_Morgan, markers);
    return List::create(Named("avgJunctions") = O.avgJunct,
                        Named("detectedJunctions") = O.avg_detected_Junctions);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// HELPER FUNCTIONS  /////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
double calculateMean(const std::vector<T>& v)
{
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = 1.0 * sum / v.size();
    return(mean);
}

template <typename T>
double calculateSD(const std::vector<T>& v)
{
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();


    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / v.size() - mean * mean);
    return stdev;
}

bool same(const Fish& A, const Fish& B) {

    for(int i = 0; i < A.chromosome1.size(); ++i) {
        if(A.chromosome1[i] != B.chromosome1[i]) {
            return false;
        }
        if(A.chromosome2[i] != B.chromosome2[i]) {
            return false;
        }
    }
    return true;

}
