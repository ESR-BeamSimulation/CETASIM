//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#ifndef FIRFeedBack_H
#define FIRFeedBack_H

#include <vector>
#include <complex>
#include "Global.h"
#include "ReadInputSettings.h"


using namespace std;
using std::vector;
using std::complex;


class FIRFeedBack
{


public:
    FIRFeedBack();
    ~FIRFeedBack();

    int delay;
    int taps;
    double kickStrengthKx;		// Eq. 116
    double kickStrengthKy;		// Eq. 116
    double kickStrengthF;		// Eq. 117
    double kickerDisp;			// [m] disperson function at kicker
    double kickerDispP;			// D(disperson)/Ds function at kicker
    double gain;

    vector<double >  firCoeffx;
    vector<double >  firCoeffy;
    vector<double >  firCoeffz;
    vector<double >  firCoeffxy;   //Nakamura's TDLSF approaches, one group pickup and kicker for 2D  
    vector<vector<double> > posxData;
    vector<vector<double> > posyData;
    vector<vector<double> > poszData;
    

    double fIRBunchByBunchFeedbackPowerLimit;// =1000;                 // power wat limit on feedback
    double fIRBunchByBunchFeedbackKickerImped;// =123E+3;              // Ohm
    double fIRBunchByBunchFeedbackKickLimit;// =0.E0;
    
    void Initial(ReadInputSettings &inputParameter);
    
      
private:

};




#endif
