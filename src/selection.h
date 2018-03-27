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


#include <Rcpp.h>
using namespace Rcpp;

int draw_prop_fitness(const std::vector<double> fitness,
                      double maxFitness);

double calculate_fitness_markers(const Fish& focal,
                                 const NumericMatrix& select);

NumericVector update_frequency(const std::vector< Fish >& v,
                               double m,
                               int num_alleles);



#endif /* Fish_hpp */
