//
//  random_functions.cpp
//  
//
//  Created by Thijs Janzen on 05/03/2018.
//
//

#include "random_functions.h"
#include <Rcpp.h>
using namespace Rcpp;

double uniform()
{
    return R::runif(1, 0, 1);
}

int random_number(int n)
{
    return (int)(R::runif(1, 0, n));
}

double poisson(double lambda)
{
    return R::rpois(lambda);
}

