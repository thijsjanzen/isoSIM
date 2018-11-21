//
//  Fish.hpp
//  
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#ifndef selection_hpp
#define selection_hpp

#include <vector>
#include <algorithm>

#include "Fish.h"
#include "main.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

int draw_prop_fitness(const std::vector<double> fitness,
                      double maxFitness);

double calculate_fitness_markers(const Fish& focal,
                                 const NumericMatrix& select);

double calculate_fitness_twoAllele(const Fish& focal,
                                   const NumericMatrix& select,
                                   bool multiplicative_selection);

NumericVector update_frequency(const std::vector< Fish >& v,
                               double m,
                               int num_alleles);

arma::mat update_all_frequencies(const std::vector< Fish >& pop,
                                 const NumericVector& select_matrix,
                                 int number_of_founders);


#endif /* Fish_hpp */
