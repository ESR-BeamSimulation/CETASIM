//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#ifndef RAMPING_H
#define RAMPING_H


#include "Global.h"
#include <vector>
#include <complex>
#include "ReadInputSettings.h"
#include "LatticeInterActionPoint.h"


using namespace std;
using std::vector;
using std::complex;


class Ramping
{


public:
    Ramping();
    ~Ramping();
    
    void RampingPara(ReadInputSettings &inputParameter,  LatticeInterActionPoint &latticeInterActionPoint,int n);
    void RampSKQ(LatticeInterActionPoint &latticeInterActionPoint,int n);

private:

};




#endif
