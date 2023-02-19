//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             


#include <vector>
#include <complex>
#include "Resonator.h"
#include <iostream>
#include<fstream>
#include<iomanip>
	
#include <numeric>
#include <cmath>

using namespace std;
using std::vector;
using std::complex;


// resonator plays roles as cavities 

Resonator::Resonator()
{

}


Resonator::~Resonator()
{
    delete dirCavFB;
    delete filterCavFB; 
}

void Resonator::GetInitialGenIg()
{
    double   argVgr = arg(resGenVol) - resDeTunePsi;
    double   absVgr = abs(resGenVol) / abs(cos(resDeTunePsi));
    
    // Ref. P.B. Wilson, Eq. (3.2.2), SLAC-PUB-6062 
    resGenVgr =  absVgr * exp(li * argVgr);
    resGenIg  =  resGenVgr * (1.0 + resCouplingBeta) / resShuntImpRs;

}

// make sure this function is only called at t=m*Trf.  
void Resonator::GetInitialCavityPower(const ReadInputSettings &inputParameter)
{
    double beamCurr = inputParameter.ringParBasic->ringCurrent;
    
    // Ref. N.G. Bill Eq. (7.40) and Eq. (7.41)
    complex<double> cavVoltage = vbAccum + resGenVol;

    resCavPower   = pow(abs(cavVoltage),2) / 2.0 / resShuntImpRs;
    resBeamPower  = beamCurr * cavVoltage.real();  
     
    // Ref. Eq. (4.1.2 Pb Wilsond)
    double coef  = pow(abs(cavVoltage),2) / 2.0 / resShuntImpRs * pow(1.0 + resCouplingBeta,2) / 4.0 / resCouplingBeta / pow(cos(resDeTunePsi),2);
    double coef0 = 2 * beamCurr * resShuntImpRs / abs(cavVoltage)  / (1.0 + resCouplingBeta) ;
    double term1 = cos(arg(cavVoltage)) + coef0 * cos(resDeTunePsi) * cos(resDeTunePsi);
    double term2 = sin(arg(cavVoltage)) + coef0 * cos(resDeTunePsi) * sin(resDeTunePsi);
    
    resGenPower   = (pow(term1,2) + pow(term2,2)) * coef; 
    resGenPowerReflect = resGenPower - resBeamPower - resCavPower;
 
}


void Resonator::CavityResonatorDynamics(double time)
{
    //Ref. T. Berenc and Borland IPAC 2015 pape, MOPMA006, Eq.(2)
    // Notice: factor of k in Eq.4  is the accelerator definition R_a/Q. 
    //cavity voltage is only solved at t=m*tRF, m=0,1,2,3.., that the golable TrackingTime=m*tRF
    
    double tB = time;
    double sigma =  2.0 * PI * resFre / (2.0 * resQualityQL);
    double deltaOmega = 2.0 * PI * resDetuneFre;
 
    complex<double> genIg =  resGenIg ;   

    // the same as matirx multiplying by matrix A of Eq.(3) in PAC 2015-MOPMA006
    resGenVol *= exp( - sigma * tB ) * exp (li * deltaOmega * tB);

    double alpha = deltaOmega * exp(- sigma * tB ) * sin(deltaOmega * tB) -  sigma      * exp(- sigma * tB ) * cos(deltaOmega * tB) + sigma;
    double beta  = sigma      * exp(- sigma * tB ) * sin(deltaOmega * tB) +  deltaOmega * exp(- sigma * tB ) * cos(deltaOmega * tB) - deltaOmega;

    double coefB = 2 * PI * resFre / 2.0 * resShuntImpRs / resQualityQ0 / (pow(sigma,2) + pow(deltaOmega,2) ); 

    double resGenVolReal = coefB * ( alpha * genIg.real() + beta  * genIg.imag());
    double resGenVolImag = coefB * (- beta * genIg.real() + alpha * genIg.imag());

    resGenVol += complex<double>(resGenVolReal, resGenVolImag); 

}







