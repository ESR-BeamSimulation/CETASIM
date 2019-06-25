#ifndef LongImpSingalBunch_H
#define LongImpSingalBunch_H

#include <vector>
#include <complex>
#include "Global.h"



using namespace std;
using std::vector;
using std::complex;


class LongImpSingalBunch
{
public:
    LongImpSingalBunch();
    ~LongImpSingalBunch();
    
    vector<double> longImpedR;   //ohm 
    vector<double> longImpedI;   //ohm
    vector<double> longImpedFre; //Hz
    void Initial();
    double binImpedFre;

private:

};




#endif
