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

    junction& operator =(const junction& other) {
        if(*this != other) {
            pos = other.pos;
            left = other.left;
            right = other.right;
        }
        return *this;
    }

    junction(const junction& other) {
        (*this) = other;
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

    Fish& operator=(const Fish& other) {
        chromosome1 = other.chromosome1;
        chromosome2 = other.chromosome2;
        return *this;
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
std::vector<double> createPopVector(const std::vector< Fish >& v);

#endif /* Fish_hpp */
