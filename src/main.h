//
//  main.h
//  
//
//  Created by Thijs Janzen on 06/03/2018.
//
//

#ifndef main_h
#define main_h

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "Fish.h"

bool is_fixed(const std::vector< Fish >& v);
std::vector< Fish > convert_NumericVector_to_fishVector(const NumericVector v);
std::vector<double> createPopVector(const std::vector< Fish >& v);
List convert_to_list(const std::vector<Fish>& v);
double calc_mean_junctions(const std::vector< Fish> & pop);

#endif /* main_h */
