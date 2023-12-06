//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
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
    
    void GetInitalCavityPowerInfo(ReadInputSettings &inputParameter); 
    void GetIntialGenIg(ReadInputSettings &inputParameter);
private:

};




#endif
