//
//  calculate_allele_spectrum.cpp
//  
//
//  Created by Thijs Janzen on 11/01/2018.
//
//

#include <stdio.h>
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

#include <Rcpp.h>
using namespace Rcpp;

void flush_console2() {
    R_FlushConsole();
    R_ProcessEvents();
    R_CheckUserInterrupt();

}

void assess_matches(const std::vector<junction>& chrom,
                    double start,
                    double end,
                    NumericMatrix& spectrum,
                    int numSteps,
                    int progress,
                    double pop_correction) {

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

    double stretch_size_normalization = pop_correction * 1.0 / (end - start);

    for(int i = 0; i < block.size(); ++i) {
        double local_right = end;
        if(i+1 < block.size()) {
            local_right = block[i+1].pos;
        }

        double local_left = block[i].pos;
        if(local_left < start) local_left = start;

        double stretch = (local_right - local_left) * stretch_size_normalization;
        int ancestor = block[i].right;
        int index = ancestor * numSteps + progress;
        if(ancestor >= 0) {
            spectrum(index, 2) += stretch;
        }
    }
    
    return;
}


NumericMatrix allele_spectrum(const std::vector<Fish>& v,
                              double step_size,
                              int numAncestors) {

    int numSteps = 1.0 / step_size;

    NumericMatrix spectrum(1 + numSteps * numAncestors, 3);

    for(int a = 0; a < numAncestors; ++a) {
        for(int i = 0; i < numSteps; ++i) {

            int index = a * numSteps + i;
            spectrum(index, 0) = i * step_size;
            spectrum(index, 1) = a;
            spectrum(index, 2) = 0;
        }
    }


    double left = 0.0;
    double right = step_size;
    double correction = 1.0 / (2.0 * v.size());

    Rcout << "0--------25--------50--------75--------100\n";
    Rcout << "*";
    int updateFreq = numSteps / 20;
    if(updateFreq < 1) updateFreq = 1;

    for(int i = 0; i < numSteps; ++i) {

        for(auto it = v.begin(); it != v.end(); ++it) {
            assess_matches((*it).chromosome1, left , right, spectrum, numSteps, i, correction);
            assess_matches((*it).chromosome2, left , right, spectrum, numSteps, i, correction);
        }
        left = right;
        right += step_size;

        if(i % updateFreq == 0) {
            Rcout << "**"; flush_console2();
        }
        Rcpp::checkUserInterrupt();
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

/*

double assess_match_old(const std::vector<junction>& chrom,
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


NumericMatrix allele_spectrum_old(const std::vector<Fish>& v,
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
        Rcpp::checkUserInterrupt();


        for(int ancestor = 0; ancestor < numAncestors; ++ancestor) {
            double local_freq = 0.0;
            for(auto it = v.begin(); it != v.end(); ++it) {

                double a = assess_match_old( (*it).chromosome1, left, right, ancestor);
                double b = assess_match_old( (*it).chromosome2, left, right, ancestor);

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

NumericMatrix calculate_allele_spectrum_cpp_old(NumericVector v1,
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

    NumericMatrix output = allele_spectrum_old(Pop, step_size, numFounders);
    
    return output;
}
*/
