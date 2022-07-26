#include "LongImpSingalBunch.h"


#include <vector>
#include <complex>
#include <iostream>
#include <fstream>

using namespace std;
using std::vector;
using std::complex;



LongImpSingalBunch::LongImpSingalBunch()
{

}


LongImpSingalBunch::~LongImpSingalBunch()
{
    
}

void LongImpSingalBunch::Initial()
{
    ifstream fin;
    fin.open("ZL_WithoutIDs_20180713_All.txt");
    
    double tempR;       //Ohm
    double tempI;       //Ohm
    double tempF;       //GHz
    
    
    int i=0;
    while(fin>>tempF && fin>>tempR && fin>>tempI)
    {
        
        tempF = tempF *1.0E9;           //Hz
        
        longImpedFre.push_back(tempF);
        longImpedR.push_back(tempR);
        longImpedI.push_back(tempI);
        i++;
    }

    fin.close();
    
    
}
