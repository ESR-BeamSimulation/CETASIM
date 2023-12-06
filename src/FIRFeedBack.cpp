//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             

// Ref. Nakamura paper "Single loop multi-dimensional digitla feedback by fir filter"
// y[0] = K \sum_0^{N} a_k x[n-k] 

#include "FIRFeedBack.h"
#include <vector>
#include <iostream>
#include<numeric>
#include <fstream>
 #include<iomanip>


using namespace std;
using std::vector;
using std::complex;


FIRFeedBack::FIRFeedBack()
{
      
}

FIRFeedBack::~FIRFeedBack()
{

}

void FIRFeedBack::Initial(ReadInputSettings &inputParameter)
{
    double beamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    int totBunchNum   = inputParameter.ringFillPatt->totBunchNumber;
    delay = inputParameter.ringBBFB->delay;
    taps  = inputParameter.ringBBFB->taps;
    gain  = inputParameter.ringBBFB->gain;
    kickerDisp = inputParameter.ringBBFB->kickerDisp;
    kickerDispP= inputParameter.ringBBFB->kickerDispP;
    fIRBunchByBunchFeedbackPowerLimit = inputParameter.ringBBFB->fIRBunchByBunchFeedbackPowerLimit;
    fIRBunchByBunchFeedbackKickerImped = inputParameter.ringBBFB->fIRBunchByBunchFeedbackKickerImped;
    kickStrengthKx = inputParameter.ringBBFB->kickStrengthK[0];
    kickStrengthKy = inputParameter.ringBBFB->kickStrengthK[1];   
    kickStrengthF  = inputParameter.ringBBFB->kickStrengthK[2];      

       
    fIRBunchByBunchFeedbackKickLimit = sqrt(2*fIRBunchByBunchFeedbackPowerLimit * fIRBunchByBunchFeedbackKickerImped)/beamEnergy; // rad unit // p=u^2/z ->u=sqrt(p*z)
     
    firCoeffx  = inputParameter.ringBBFB->firCoeffx;
    firCoeffy  = inputParameter.ringBBFB->firCoeffy;
    firCoeffz  = inputParameter.ringBBFB->firCoeffz;
    firCoeffxy = inputParameter.ringBBFB->firCoeffxy;


    int firOrder = delay + taps;
    posxData.resize(firOrder);
    posyData.resize(firOrder);
    poszData.resize(firOrder);
    

    for(int i=0;i<firOrder;i++)
    {
        posxData[i].resize(totBunchNum);
        posyData[i].resize(totBunchNum);
        poszData[i].resize(totBunchNum);
    }

    for(int i=0;i<firOrder;i++)
    {
        for(int j=0;j<totBunchNum;j++)
        {
            posxData[i][j] = 0.E0;
            posyData[i][j] = 0.E0;
            poszData[i][j] = 0.E0;
        }
    }
    
//    double tempCoefx=1;
//    double tempCoefy=1;

//    kickStrengthKx = tempCoefx * 1.0E-1; 
//    kickStrengthKy = tempCoefy * 1.0E-1;
//    kickStrengthF = 0.E0;

    cout<<"--------Sum of bunch-by-bunch feedbakc coeff--DC rejection------------"<<endl;
    double sumx  =  accumulate(begin(firCoeffx), end(firCoeffx), 0.0);
    double sumy  =  accumulate(begin(firCoeffy), end(firCoeffy), 0.0);
    double sumz  =  accumulate(begin(firCoeffz), end(firCoeffz), 0.0);
    double sumxy =  accumulate(begin(firCoeffxy),end(firCoeffxy), 0.0);
    cout<<setw(40)<<left<<"bunch-by-bunch sum coefx "<<setw(10)<<left<< sumx<<endl;
    cout<<setw(40)<<left<<"bunch-by-bunch sum coefy "<<setw(10)<<left<< sumy<<endl;
    cout<<setw(40)<<left<<"bunch-by-bunch sum coefz "<<setw(10)<<left<< sumz<<endl;
    cout<<setw(40)<<left<<"bunch-by-bunch sum coefxy "<<setw(10)<<left<< sumxy<<endl;


}
