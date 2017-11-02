//
//  Fish.cpp
//  
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#include "Fish.h"
#include "randomc.h"
#include <algorithm>



double getRecomPos() {
    double pos = uniform();
    while(pos == 0 || pos == 1.0) {
        pos = uniform(); //code to avoid drawing exactly the borders of the chromosome
    }

    return pos;
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
