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

void Resonator::GetInitialResonatorGenIg(const ReadInputSettings &inputParameter)
{
    if(resType==0)     // passive cavity
    {
        resGenIg = complex<double>(0.E0, 0.E0);
    }
    else        // active cavity
    {
        double   argVgr = arg(resGenVol) - resDeTunePsi;
        double   absVgr = abs(resGenVol) / abs(cos(resDeTunePsi));
        
        // Ref. P.B. Wilson, Eq. (3.2.2), SLAC-PUB-6062  -- simple but not accurate. 
        resGenVgr =  absVgr * exp(li * argVgr);
        resGenIg  =  resGenVgr * (1.0 + resCouplingBeta) / resShuntImpRs;
    }

   

    // // try to get resGenIg by iteration process 
    // int nturns = 20;
    // int ringHarm        = inputParameter.ringParRf->ringHarm;
    // double t0           = inputParameter.ringParBasic->t0;
    // double tRF          = t0 / ringHarm;



    // // get initial resGenVol as reference 
    // resGenVol = complex<double>(0.E0,0.E0);
    // for (int n=0;n<nturns;n++)
    // {
    //     for(int h=0;h<ringHarm;h++)
    //     {
    //         ResonatorDynamics(tRF);
    //     }
    // }
    // // resGenVol is updated by cavity-dynamics

    // complex<double> resGenIg_0, resGenIg_1;
    // double temp;

    // for(int k=0;k<10;k++)
    // {
    //     temp = abs(resGenVol + vbAccum) - abs(resCavVolReq);
    //     if(temp < 0)
    //     {   
    //         resGenIg_0 =       resGenIg;
    //         resGenIg_1 = 1.1 * resGenIg ;
    //     }
    //     else
    //     {
    //         resGenIg_0 = 0.9 * resGenIg;
    //         resGenIg_1 =       resGenIg;
    //     }

    //     while( abs(temp) > 1.0E1 )     // target to ensure amplitude 
    //     {     
    //         resGenVol = complex<double>(0.E0,0.E0);
    //         resGenIg = (resGenIg_0 + resGenIg_1)/2.0;

    //         for (int n=0;n<nturns;n++)
    //         {
    //             for(int h=0;h<ringHarm;h++)
    //             {
    //                 ResonatorDynamics(tRF);
    //             }
    //         }

    //         temp = abs(resGenVol + vbAccum) - abs(resCavVolReq);
            
    //         if (temp >0)
    //         {
    //             resGenIg_0 = resGenIg_0;
    //             resGenIg_1 = resGenIg;
    //         }
    //         else
    //         {
    //             resGenIg_0 = resGenIg;
    //             resGenIg_1 = resGenIg_1;
    //         }
    //         cout<<k<<"  "<<temp<<endl;
    //     }

    //     double resGenIgAmp = abs(resGenIg);
    //     double resGenIgArg = arg(resGenIg);

    //     double resGenIgPhase_0,resGenIgPhase_1,resGenIgPhase;
    //     temp = arg(resGenVol + vbAccum) - arg(resCavVolReq);
        
    //     if(temp<0)
    //     {
    //         resGenIgPhase_0 = arg(resGenIg);
    //         resGenIgPhase_1 = arg(resGenIg) * 1.1;
    //     }
    //     else
    //     {
    //         resGenIgPhase_0 = 0.9 * arg(resGenIg);
    //         resGenIgPhase_1 = arg(resGenIg);
    //     }

    //     temp = 1;
    //     while (abs(temp)>1.E-5)   // target to ensure phase 
    //     {
    //         resGenVol = complex<double>(0.E0,0.E0);

    //         resGenIg_0 = resGenIgAmp * exp(li * resGenIgPhase_0);
    //         resGenIg_1 = resGenIgAmp * exp(li * resGenIgPhase_1);
    //         resGenIg = (resGenIg_0 + resGenIg_1)/2.0;
            
    //         for (int n=0;n<nturns;n++)
    //         {
    //             for(int h=0;h<ringHarm;h++)
    //             {
    //                 ResonatorDynamics(tRF);
    //             }
    //         }

    //         temp = arg(resGenVol + vbAccum) - arg(resCavVolReq);

    //         if (temp >0)
    //         {
    //             resGenIgPhase_0 = resGenIgPhase_0;
    //             resGenIgPhase_1 = arg(resGenIg);
    //         }
    //         else
    //         {
    //             resGenIgPhase_0 = arg(resGenIg);
    //             resGenIgPhase_1 = resGenIgPhase_1;
    //         }
            
    //     }
    // }

    // cout<<left<<setw(21)<<abs(resGenIg_0)
    //     <<left<<setw(21)<<arg(resGenIg_0)
    //     <<left<<setw(21)<<abs(resGenIg_1)
    //     <<left<<setw(21)<<arg(resGenIg_1)
    //     <<left<<setw(21)<<abs(resGenIg)
    //     <<left<<setw(21)<<arg(resGenIg)
    //     <<left<<setw(21)<<abs(resGenVol + vbAccum)
    //     <<left<<setw(21)<<arg(resGenVol + vbAccum)<<endl;

    // getchar();

}

// make sure this function is only called at t=m*Trf.  
void Resonator::GetInitialResonatorPower(const ReadInputSettings &inputParameter)
{
    double beamCurr = inputParameter.ringParBasic->ringCurrent;
    
    // Ref. N.G. Bill Eq. (7.40) and Eq. (7.41)
    complex<double> cavVoltage = vbAccum + resGenVol;

    resCavPower   = pow(abs(cavVoltage),2) / 2.0 / resShuntImpRs;
    resBeamPower  = beamCurr * cavVoltage.real();  
    
    // Ref. Eq. (4.1.2 P. B. Wilson)
    double coef  = pow(abs(cavVoltage),2) / 2.0 / resShuntImpRs * pow(1.0 + resCouplingBeta,2) / 4.0 / resCouplingBeta / pow(cos(resDeTunePsi),2);
    double coef0 = 2 * beamCurr * resShuntImpRs / abs(cavVoltage)  / (1.0 + resCouplingBeta);
    double term1 = cos(arg(cavVoltage)) + coef0 * cos(resDeTunePsi) * cos(resDeTunePsi);
    double term2 = sin(arg(cavVoltage)) + coef0 * cos(resDeTunePsi) * sin(resDeTunePsi);
    
    resGenPower   = (pow(term1,2) + pow(term2,2)) * coef; 
    resGenPowerReflect = resGenPower - resBeamPower - resCavPower; 
}


void Resonator::GetBeamInducedVol(double timeToNextBunch)
{    
  
    double deltaL = timeToNextBunch / tF;        
    double cPsi   = 2.0 * PI * resFre * timeToNextBunch;       
    vbAccum *= exp( - deltaL ) * exp (li * cPsi);
    
    // deltaL = tB / tF ;
    // cPsi   = deltaL * tan(cavityResonator.resonatorVec[i].resDeTunePsi )

    // cout<<left<<setw(15)<<timeToNextBunch
    //     <<left<<setw(15)<<tF
    //     <<left<<setw(15)<<deltaL
    //     <<left<<setw(15)<<cPsi
    //     <<left<<setw(15)<<abs(exp( - deltaL ) * exp (li * cPsi))
    //     <<left<<setw(15)<<arg(exp( - deltaL ) * exp (li * cPsi))
    //     <<left<<setw(15)<<abs(vbAccum)
    //     <<left<<setw(15)<<arg(vbAccum)
    //     <<endl;

    
    // getchar();

    // cout<<abs(vbAccum)<<endl;
    // cout<<arg(vbAccum)<<endl;
    
    
}


void Resonator::ResonatorDynamics(double time)
{
    //Ref. T. Berenc and Borland IPAC 2015 paper, MOPMA006, Eq.(2)
    // Notice: factor of k in Eq.4  is the accelerator definition R_a/Q. 
    //cavity voltage is only solved at t=m*tRF, m=0,1,2,3.., that the golable TrackingTime=m*tRF
    
    double tB = time;
    double sigma =  2.0 * PI * resFre / (2.0 * resQualityQL);
    double deltaOmega = 2.0 * PI * resDetuneFre;
 
    complex<double> genIg =  resGenIg;   

    // the same as matirx multiplying by matrix A of Eq.(3) in PAC 2015-MOPMA006
    resGenVol *= exp( - sigma * tB ) * exp (li * deltaOmega * tB);

    double alpha = deltaOmega * exp(- sigma * tB ) * sin(deltaOmega * tB) -  sigma      * exp(- sigma * tB ) * cos(deltaOmega * tB) + sigma;
    double beta  = sigma      * exp(- sigma * tB ) * sin(deltaOmega * tB) +  deltaOmega * exp(- sigma * tB ) * cos(deltaOmega * tB) - deltaOmega;

    double coefB = 2 * PI * resFre / 2.0 * resShuntImpRs / resQualityQ0 / (pow(sigma,2) + pow(deltaOmega,2) ); 

    double resGenVolReal = coefB * ( alpha * genIg.real() + beta  * genIg.imag());
    double resGenVolImag = coefB * (- beta * genIg.real() + alpha * genIg.imag());

    resGenVol += complex<double>(resGenVolReal, resGenVolImag); 

}


void Resonator::GetResonatorInfoAtNTrf(int harmonicNum,double dt)
{
    double tB = dt;
    double deltaL = tB / tF;        
    double cPsi   = 2.0 * PI * resFre * tB;       

    vBSample[harmonicNum]   = vbAccum *  exp( - deltaL ) * exp (li * cPsi);
    
    vGenSample[harmonicNum] = resGenVol;
    vCavSample[harmonicNum] = vBSample[harmonicNum] + vGenSample[harmonicNum];
    deltaVCavSample[harmonicNum] = vCavSample[harmonicNum] - resCavVolReq;
    vCavDueToDirFB[harmonicNum] = - deltaVCavSample[harmonicNum];  
}


    
    

    





