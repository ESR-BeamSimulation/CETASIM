//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************

#ifndef BoardBandImp_H
#define BoardBandImp_H

#include <vector>
#include <complex>
#include "Global.h"
#include "LatticeInterActionPoint.h"
#include "ReadInputSettings.h"
#include "Spline.h"
using namespace std;
using std::vector;
using std::complex;

class BoardBandImp
{


public:
    BoardBandImp();
    ~BoardBandImp();

    vector<complex<double> > zZImp;
    vector<complex<double> > zDxImp;
    vector<complex<double> > zDyImp;
    vector<complex<double> > zQxImp;
    vector<complex<double> > zQyImp;
    vector<double>           freq   ;
    int                      freSampN;

    int nBins;
    double freqMax;  //   impedance is from 0 to Fmax
    double dt;       //   dt  = 1/ (2* fMax) -- naquist certeria
    double tMax;
    double zMax;
    double dz;     
    vector<double> binPosZ;

    void ReadInImp(const ReadInputSettings &inputParameter);

private:

};




#endif
