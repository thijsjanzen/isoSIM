//
//  Output.cpp
//  
//
//  Created by Thijs Janzen on 02/11/2017.
//
//

#include "Output.h"
#include "Fish.h"
#include <vector>



void Output::update(const std::vector<Fish>& Pop) {
    double averageNumJunctions = 0;
    for(std::vector<Fish>::const_iterator i = Pop.begin(); i != Pop.end(); ++i)
    {
        int numJ  = (int)(*i).chromosome1.size() - 2; //exclude the ends
        numJ += (int)(*i).chromosome2.size() - 2; //exclude the ends
        averageNumJunctions += numJ;
    }

    averageNumJunctions = 1.0 * averageNumJunctions / (2 * (int)Pop.size());

    avgJunct.push_back(averageNumJunctions);

    return;
}

int countJunctions(const std::vector<bool>& B) {
    int numJunctions = 0;
    for(std::size_t i = 1; i < B.size(); ++i) {
        if(B[i] != B[i-1]) {
            numJunctions++;
        }
    }
    return numJunctions;
}

std::vector<bool> detectJunctions(const std::vector<junction>& G,
                                  const std::vector<double>& markers) {
    std::vector<bool> output(markers.size());

    int j = 0;
    for(int i = 0; i < markers.size(); ++i) {
        double focalPos = markers[i];
        for(; j <= (G.size()-1); ++j) {
            double left = G[j].pos;
            double right = G[j+1].pos;
            if(left <= focalPos && right >= focalPos) {

                output[i] = ((bool)G[j].right);
                break;
            }
        }
    }
    return output;
}

void Output::detectNumJunctions(const std::vector<Fish> &Pop,
                                const std::vector<double> &markers) {
    double averageNumJunctions = 0;
    for(std::vector<Fish>::const_iterator i = Pop.begin(); i != Pop.end(); ++i) {
        std::vector<bool> genome1 = detectJunctions((*i).chromosome1, markers);
        averageNumJunctions += countJunctions(genome1);

        std::vector<bool> genome2 = detectJunctions((*i).chromosome2, markers);
        averageNumJunctions += countJunctions(genome2);
    }

    averageNumJunctions = 1.0 * averageNumJunctions / (2 * Pop.size()); //diploid
    avg_detected_Junctions.push_back(averageNumJunctions);
    return;
}

std::istream& operator >> (std::istream& is, junction& j)
{
    is >> j.left;
    is >> j.pos;
    is >> j.right;

    return is;
}

std::ostream& operator << (std::ostream& os, const junction& j)
{
    os << j.left << "\t";
    os << j.pos << "\t";
    os << j.right << "\t";

    return os;
}

std::ostream& operator << (std::ostream& os, const std::vector< junction >& chrom)
{
    for(auto it = chrom.begin(); it != chrom.end(); ++it) {
        os << (*it);
    }

    return os;
}

void writePoptoFile(const std::vector<Fish>& Pop, std::string filename)
{
    std::ofstream outFile(filename.c_str());
    for(auto indiv = Pop.begin(); indiv != Pop.end(); ++indiv) {
        outFile << (*indiv).chromosome1 << "\n";
        outFile << (*indiv).chromosome2 << "\n";
    }
    outFile.close();
    return;
}

void readPopfromFile(std::vector<Fish>& Pop, std::string filename)
{
    Pop.clear();
    std::ifstream inFile(filename.c_str());
    if(!inFile.is_open()) {
        std::cout << "Could not open file\n";
        return;
    }
    std::string line1;
    std::string line2;
    while(std::getline(inFile, line1) &&
          std::getline(inFile, line2)) {

        Fish temp;
        std::istringstream iss1(line1);
        while(!iss1.eof()) {
            junction temp_junction;
            iss1 >> temp_junction;
            temp.chromosome1.push_back(temp_junction);
        }
        temp.chromosome1.pop_back(); //!.eof always reads last line twice for some fucked up reason

        std::istringstream iss2(line2);
        while(!iss2.eof()) {
            junction temp_junction;
            iss2 >> temp_junction;
            temp.chromosome2.push_back(temp_junction);
        }
        temp.chromosome2.pop_back(); //!.eof always reads last line twice for some fucked up reason
        Pop.push_back(temp);
    }
    return;
}
