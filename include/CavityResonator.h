#ifndef CAVITYRESONATOR_H
#define CAVITYRESONATOR_H


#include "Global.h"
#include <vector>
#include <complex>
#include "ReadInputSettings.h"
#include "Resonator.h"


using namespace std;
using std::vector;
using std::complex;


class CavityResonator
{


public:
    CavityResonator();
    ~CavityResonator();
    
    void Initial(ReadInputSettings &inputParameter);
       
    //void Initial(Beam &beam, ReadInputSettings &inputParameter);
    //void GetGenPar(Beam &beam, ReadInputSettings &inputParameter);

    vector<Resonator> resonatorVec;
    
    double GetWakefunction(double tauij);      // used to get the longitudinal wakefunction. (short range) BBR model as usuall
    

private:

};




#endif
