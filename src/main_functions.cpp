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

#include <Rcpp.h>
using namespace Rcpp;

/***********************************************************************
 Start of simulation code
 ***********************************************************************/



struct Output {
    std::vector<double> avgJunct;
    std::vector<double> avg_detected_Junctions;
    void update(const std::vector< Fish >& Pop);
    void detectNumJunctions(const std::vector<Fish> &Pop,
                            const std::vector<double> &markers);
};

double getRecomPos();

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
////////////////////////   MATING FUNCTIONS ////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void Output::update(const std::vector<Fish>& Pop) {
    double averageNumJunctions = 0;
    for(std::vector<Fish>::const_iterator i = Pop.begin(); i != Pop.end(); ++i)
    {
        int numJ  = (int)(*i).chromosome1.size() - 2; //exclude the ends
        numJ += (int)(*i).chromosome2.size() - 2; //exclude the ends
        averageNumJunctions += numJ;
    }

    averageNumJunctions = 1.0 * averageNumJunctions / (2 * (int)Pop.size());

    avgJunct.push_back(averageNumJunctions);

    return;
}

int countJunctions(const std::vector<bool>& B) {
    int numJunctions = 0;
    for(std::size_t i = 1; i < B.size(); ++i) {
        if(B[i] != B[i-1]) {
            numJunctions++;
        }
    }
    return numJunctions;
}

std::vector<bool> detectJunctions(const std::vector<junction>& G,
                                  const std::vector<double>& markers) {
    std::vector<bool> output(markers.size());

    int j = 0;
    for(int i = 0; i < markers.size(); ++i) {
        double focalPos = markers[i];
        for(; j <= (G.size()-1); ++j) {
            double left = G[j].pos;
            double right = G[j+1].pos;
            if(left <= focalPos && right >= focalPos) {

                output[i] = ((bool)G[j].right);
                break;
            }
        }
    }
    return output;
}

void Output::detectNumJunctions(const std::vector<Fish> &Pop,
                                const std::vector<double> &markers) {
    double averageNumJunctions = 0;
    for(std::vector<Fish>::const_iterator i = Pop.begin(); i != Pop.end(); ++i) {
        std::vector<bool> genome1 = detectJunctions((*i).chromosome1, markers);
        averageNumJunctions += countJunctions(genome1);

        std::vector<bool> genome2 = detectJunctions((*i).chromosome2, markers);
        averageNumJunctions += countJunctions(genome2);
    }

    averageNumJunctions = 1.0 * averageNumJunctions / (2 * Pop.size()); //diploid
    avg_detected_Junctions.push_back(averageNumJunctions);
    return;
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
