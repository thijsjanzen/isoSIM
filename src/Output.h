//
//  Output.hpp
//  
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#ifndef Output_hpp
#define Output_hpp

#include <stdio.h>
#include <vector>
#include <sstream>
#include <string>

#include "Fish.h"
#include <Rcpp.h>
using namespace Rcpp;

struct Output {
    std::vector<double> avgJunct;
    std::vector<double> avg_detected_Junctions;
    
    void update(const std::vector< Fish >& Pop);
    void detectNumJunctions(const std::vector< Fish>    &Pop,
                            const std::vector< double > &markers);
};

void writePoptoFile( const std::vector< Fish >& Pop,
                           std::string filename);

void readPopfromFile(      std::vector< Fish >& Pop,
                           std::string filename);

double calculate_heterozygosity(const std::vector < Fish >& v);
void flush_console();
bool is_fixed(const std::vector< Fish >& v);
std::vector< Fish > convert_NumericVector_to_fishVector(const NumericVector v);


#endif /* Output_hpp */
