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

#endif /* Output_hpp */
