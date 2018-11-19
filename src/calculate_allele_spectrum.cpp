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

#include <Rcpp.h>
using namespace Rcpp;

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
                              int max_num_ancestors,
                              bool progress_bar) {

    int numSteps = 1.0 / step_size;

    NumericMatrix spectrum(1 + numSteps * (1 + max_num_ancestors), 3);

    for(int a = 0; a < (1 + max_num_ancestors); ++a) {
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

    if(progress_bar) {
        Rcout << "0--------25--------50--------75--------100\n";
        Rcout << "*";
    }

    int updateFreq = numSteps / 20;
    if(updateFreq < 1) updateFreq = 1;

    for(int i = 0; i < numSteps; ++i) {

        for(auto it = v.begin(); it != v.end(); ++it) {
            assess_matches((*it).chromosome1, left , right, spectrum, numSteps, i, correction);
            assess_matches((*it).chromosome2, left , right, spectrum, numSteps, i, correction);
        }
        left = right;
        right += step_size;

        if(i % updateFreq == 0 && progress_bar) {
            Rcout << "**";
        }
        Rcpp::checkUserInterrupt();
    }
    return spectrum;
}


// [[Rcpp::export]]
NumericMatrix calculate_allele_spectrum_cpp(NumericVector v1,
                                            double step_size,
                                            bool progress_bar)
{
    std::vector< Fish > Pop;

    std::vector<double> v = Rcpp::as<std::vector<double> >(v1);

    Fish temp;
    int indic_chrom = 1;
    bool add_indiv = false;

    int max_num_ancestor = -1;

    for(int i = 0; i < (v.size() - 1); i += 2) {
        junction temp_j;
        temp_j.pos = v[i];
        if(i+1 > v.size()) break;
        temp_j.right = v[i+1];

        if(temp_j.right > max_num_ancestor) max_num_ancestor = temp_j.right;

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
    NumericMatrix output = allele_spectrum(Pop, step_size, max_num_ancestor, progress_bar);
    
    return output;
}
