//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************

#ifndef BEAMION2DPIC_H
#define BEAMION2DPIC_H

#include <vector>
#include <complex>
#include <fftw3.h>
#include "ReadInputSettings.h"
#include "PIC2D.h"
using namespace std;
using std::complex;


using std::vector;
using v1i= vector<int> ;
using v1d= vector<double> ;
using v2d= vector<vector<double> > ;
using v3d= vector<vector<vector<double> > > ;


class BeamIon2DPIC
{

public:
    BeamIon2DPIC();
    ~BeamIon2DPIC();

    vector<int> meshNumBeam{128,128};
    vector<int> meshNumIon{128,128};

    PIC2D pic2DBeam;
    PIC2D pic2DIon; 

    void InitialPIC2D();

};




#endif
