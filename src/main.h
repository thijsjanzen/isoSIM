//
//  main.h
//  
//
//  Created by Thijs Janzen on 06/03/2018.
//
//

#ifndef main_h
#define main_h

#include <Rcpp.h>
using namespace Rcpp;

#include "Fish.h"

bool is_fixed(const std::vector< Fish >& v);
std::vector< Fish > convert_NumericVector_to_fishVector(const NumericVector v);
std::vector<double> createPopVector(const std::vector< Fish >& v);
void flush_console();


#endif /* main_h */
