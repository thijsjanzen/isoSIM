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

#include <Rcpp.h>
using namespace Rcpp;

/***********************************************************************
 Start of simulation code
 ***********************************************************************/

struct junction {
    double pos;
    int left;
    int right;

    junction()  {}

    junction(double loc, int A, int B)  {
        pos = loc;
        left = A;
        right = B;
    }

    junction(const junction& other) {
        pos = other.pos;
        left = other.left;
        right = other.right;
    }

    bool operator ==(const junction& other) const {
        if(pos != other.pos) return false;
        if(left != other.left) return false;
        if(right != other.right) return false;

        return true;
    }

    bool operator <(const junction& other) const {
        return(pos < other.pos);
    }

    bool operator !=(const junction& other) const {
        return( !( (*this) == other) );
    }
};

struct Fish {
    std::vector< junction > chromosome1;
    std::vector< junction > chromosome2;

    Fish()
    {}

    Fish(int initLoc)    {
        junction left(0.0, -1, initLoc);
        junction right(1, initLoc, -1);
        chromosome1.push_back( left  );
        chromosome1.push_back( right );
        chromosome2.push_back( left  );
        chromosome2.push_back( right );
    }

    Fish(const std::vector<junction>& A,
         const std::vector<junction>& B)    {
        chromosome1 = A;
        chromosome2 = B;
    }
};


struct Output {
    std::vector<double> avgJunct;
    std::vector<double> avg_detected_Junctions;
    void update(const std::vector< Fish >& Pop);
    void detectNumJunctions(const std::vector<Fish> &Pop,
                            const std::vector<double> &markers);
};

Fish mate(const Fish& A, const Fish& B,
          double numRecombinations);

double getRecomPos();

void Recombine(std::vector<junction>& offspring,
               std::vector<junction> chromosome1,
               std::vector<junction> chromosome2,
               double numRecombinations);

bool isFixed(const std::vector<Fish>& V);
void macstart(const char * argv[]);  // forward declaration
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

double getRecomPos() {
    double pos = uniform();
    while(pos == 0 || pos == 1.0) {
        pos = uniform(); //code to avoid drawing exactly the borders of the chromosome
    }

    return pos;
}

void addJunction(std::vector<junction>& offspring,
                 const std::vector<junction>& focalChrom,
                 double recomPos,
                 double nextRecomPos) {

    bool added = false;
    for(int i = 0; i < focalChrom.size(); ++i)
    {
        if(focalChrom[i].pos > recomPos & focalChrom[i].pos < nextRecomPos) {
            if(!added) {
                added = true;
                junction toAdd(recomPos,
                               offspring.back().right ,
                               focalChrom[i].left);
                if(toAdd.left != toAdd.right) offspring.push_back(toAdd); //only add true recombinations
            }
            offspring.push_back(focalChrom[i]);
        }
    }

    return;
}

void Recombine(std::vector<junction>& offspring,
               std::vector<junction> chromosome1,
               std::vector<junction> chromosome2,
               double MORGAN)  {

    std::vector<double> recomPos;

    int numRecombinations = poisson(MORGAN);

    while (recomPos.size() < numRecombinations) {
        double pos = getRecomPos();
        recomPos.push_back(pos);
        // sort them, in case they are not sorted yet
        // we need this to remove duplicates, and later
        // to apply crossover
        std::sort(recomPos.begin(), recomPos.end() );
        // remove duplicate recombination sites
        recomPos.erase(std::unique(recomPos.begin(), recomPos.end()), recomPos.end());
    }

    std::vector< junction > toAdd; //first create junctions on exactly the recombination positions
    for(int i = 0; i < recomPos.size(); ++i) {
        junction temp;
        temp.pos = recomPos[i];
        toAdd.push_back(temp);
    }

    for(int i = 1; i < chromosome1.size(); ++i) {
        double leftpos = chromosome1[i-1].pos;
        double rightpos = chromosome1[i].pos;

        for(int j = 0; j < recomPos.size(); ++j) {
            if(recomPos[j] > leftpos && recomPos[j] < rightpos) {
                if(j % 2 == 0) { //even, so chrom1 = L, chrom2 = R
                    toAdd[j].left = chromosome1[i].left;
                }
                if(j % 2 == 1) { //uneven so chrom1 = R, chrom2 = L
                    toAdd[j].right = chromosome1[i].left;
                }
            }
        }
    }

    for(int i = 1; i < chromosome2.size(); ++i) {
        double leftpos = chromosome2[i-1].pos;
        double rightpos = chromosome2[i].pos;

        for(int j = 0; j < recomPos.size(); ++j) {
            if(recomPos[j] > leftpos && recomPos[j] < rightpos) {
                if(j % 2 == 0) { //even, so chrom1 = L, chrom2 = R
                    toAdd[j].right = chromosome2[i].left;
                }
                if(j % 2 == 1) { //uneven so chrom1 = R, chrom2 = L
                    toAdd[j].left = chromosome2[i].left;
                }
            }
        }
    }

    for(int i = 0; i < toAdd.size(); ++i) {
        if(toAdd[i].left != toAdd[i].right) {
            offspring.push_back(toAdd[i]);
        }
    }

    //now we have to add the other junctions from chrom1 and chrom2.
    double leftpos = 0;
    double rightpos = 0;


    for(int i = 0; i < (recomPos.size() + 1); ++i) {
        rightpos = 1.0;
        if(i < recomPos.size()) rightpos = recomPos[i];
        if(i % 2 == 0) { //even, so take from chromosome 1
            for(std::vector<junction>::iterator it = chromosome1.begin(); it != chromosome1.end(); ++it) {
                if((*it).pos >= leftpos && (*it).pos <= rightpos) {
                    offspring.push_back((*it));
                }
                if((*it).pos > rightpos) { //we are past the recombination section
                    break;
                }
            }
        }

        if(i % 2 == 1) { //odd, so take from chromosome 2
            for(std::vector<junction>::iterator it = chromosome2.begin(); it != chromosome2.end(); ++it) {
                if((*it).pos >= leftpos && (*it).pos <= rightpos) {
                    offspring.push_back((*it));
                }
                if((*it).pos > rightpos) { //we are past the recombination section
                    break;
                }
            }
        }

        //move forward
        leftpos = rightpos;
    }

    std::sort(offspring.begin(), offspring.end());
    offspring.erase(std::unique(offspring.begin(), offspring.end()), offspring.end());

    return;
}



Fish mate(const Fish& A, const Fish& B, double numRecombinations)
{
    Fish offspring;
    offspring.chromosome1.clear();
    offspring.chromosome2.clear(); //just to be sure.

    //first the father chromosome
    int event = random_number(2);
    switch(event) {
        case 0:  {
            Recombine(offspring.chromosome1, A.chromosome1, A.chromosome2, numRecombinations);
            break;
        }
        case 1: {
            Recombine(offspring.chromosome1, A.chromosome2, A.chromosome1, numRecombinations);
            break;
        }
    }


    //then the mother chromosome
    event = random_number(2);
    switch(event) {
        case 0:  {
            Recombine(offspring.chromosome2, B.chromosome1, B.chromosome2, numRecombinations);
            break;
        }
        case 1: {
            Recombine(offspring.chromosome2, B.chromosome2, B.chromosome1, numRecombinations);
            break;
        }
    }

    return offspring;
}


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
