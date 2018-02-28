//
//  calculate_stats.cpp
//  
//
//  Created by Thijs Janzen on 28/02/2018.
//
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
