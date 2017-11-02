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
    int left;
    int right;

    junction()  {}

    junction(double loc, int A, int B)  {
        pos = loc;
        left = A;
        right = B;
    }

    junction(const junction& other) {
        pos = other.pos;
        left = other.left;
        right = other.right;
    }

    bool operator ==(const junction& other) const {
        if(pos != other.pos) return false;
        if(left != other.left) return false;
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
        junction left(0.0, -1, initLoc);
        junction right(1, initLoc, -1);
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
};

Fish mate(const Fish& A, const Fish& B, double numRecombinations);

#endif /* Fish_hpp */
