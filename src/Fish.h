//
//  Fish.hpp
//  
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#ifndef Fish_hpp
#define Fish_hpp

#include <stdio.h>
#include <vector>

struct junction {
    double pos;
    int right;

    junction()  {}

    junction(double loc, int A)  {
        pos = loc;
        right = A;
    }

    junction(const junction& other) {
        pos = other.pos;
        right = other.right;
    }

    bool operator ==(const junction& other) const {
        if(pos != other.pos) return false;
        if(right != other.right) return false;

        return true;
    }

    bool operator <(const junction& other) const {
        return(pos < other.pos);
    }

    bool operator !=(const junction& other) const {
        return( !( (*this) == other) );
    }
};

struct Fish {
    std::vector< junction > chromosome1;
    std::vector< junction > chromosome2;

    Fish()
    {}

    Fish(int initLoc)    {
        junction left(0.0, initLoc);
        junction right(1.0, -1);
        chromosome1.push_back( left  );
        chromosome1.push_back( right );
        chromosome2.push_back( left  );
        chromosome2.push_back( right );
    }

    Fish(const std::vector<junction>& A,
         const std::vector<junction>& B)    {
        chromosome1 = A;
        chromosome2 = B;
    }

    bool operator ==(const Fish& other) const {
        if(chromosome1.size() != other.chromosome1.size()) return false;
        if(chromosome2.size() != other.chromosome2.size()) return false;

        double eps = 1e-4;

        for(int i = 0; i < (int)chromosome1.size(); ++i) {
            double diffPos = chromosome1[i].pos - other.chromosome1[i].pos;
            if(diffPos < -eps || diffPos > eps) return false;
            if(chromosome1[i].right != other.chromosome1[i].right) return false;
        }

        for(int i = 0; i < (int)chromosome2.size(); ++i) {
            double diffPos = chromosome2[i].pos - other.chromosome2[i].pos;
            if(diffPos < -eps || diffPos > eps) return false;
            if(chromosome2[i].right != other.chromosome2[i].right) return false;
        }
        
        return true;
    }
    
    
};


Fish mate(const Fish& A, const Fish& B, double numRecombinations);

#endif /* Fish_hpp */
