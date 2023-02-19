//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once

#include "SPBunch.h"
#include "Global.h"
#include "Faddeeva.h"
#include "WakeFunction.h"
#include "Spline.h"
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <cmath>
#include <random>
#include <gsl/gsl_fft_complex.h>
#include <algorithm>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_histogram.h>


SPBunch::SPBunch()
{
}

SPBunch::~SPBunch()
{ 
}

void SPBunch::InitialSPBunch(const  ReadInputSettings &inputParameter)
{    
    haissinski->cavAmp.resize(inputParameter.ringParRf->resNum,0);
    haissinski->cavPhase.resize(inputParameter.ringParRf->resNum,0);   
    //call the Initial at base class (Bunch.h)
    Bunch::Initial(inputParameter);           
}


void SPBunch::DistriGenerator(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter, int bunchIndex)
{
    // single particle always located at (0,0,0,0,0,0), Dis errors are defined by input
 
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> doffset{0,1};

    double initialStaticOffSet[6];
    double initialDynamicOffSet[6];

    for(int i=0;i<6;i++)
    {
       initialStaticOffSet[i]   = inputParameter.ringBunchPara->initialStaticOffSet[i];
       initialDynamicOffSet[i]  = inputParameter.ringBunchPara->initialDynamicOffSet[i];       
    }

    srand(time(0)+bunchIndex);
    
    double temp;
    double disDx,disDy,disDz,disMx,disMy,disMz;
    
    do{
        temp  = doffset(gen);            
        disDx = temp * initialDynamicOffSet[0] + initialStaticOffSet[0];
    }while(temp>3); 
    
    do{
        temp  = doffset(gen);            
        disDy = temp * initialDynamicOffSet[1] + initialStaticOffSet[1];
    }while(temp>3); 

    do{
        temp  = doffset(gen);            
        disDz = temp * initialDynamicOffSet[2] + initialStaticOffSet[2];
    }while(temp>3);         
    
    do{
        temp  = doffset(gen);            
        disMx = temp * initialDynamicOffSet[3] + initialStaticOffSet[3];
    }while(temp>3);
    
    do{
        temp  = doffset(gen);            
        disMy = temp * initialDynamicOffSet[4] + initialStaticOffSet[4];
    }while(temp>3);
    
    do{
        temp  = doffset(gen);            
        disMz = temp * initialDynamicOffSet[5] + initialStaticOffSet[5];
    }while(temp>3);
    

    ePositionZ[0] =  disDz;
    eMomentumZ[0] =  disMz;

  	double dispersionX  =  latticeInterActionPoint.twissDispX[0];
  	double dispersionY  =  latticeInterActionPoint.twissDispY[0];
  	double dispersionPX =  latticeInterActionPoint.twissDispPX[0];
  	double dispersionPY =  latticeInterActionPoint.twissDispPY[0];

	ePositionX[0] =  disDx + dispersionX  * eMomentumZ[0];
	ePositionY[0] =  disDy + dispersionY  * eMomentumZ[0];
	eMomentumX[0] =  disMx + dispersionPX * eMomentumZ[0];
	eMomentumY[0] =  disMy + dispersionPY * eMomentumZ[0];

    zAverLastTurn = ePositionZ[0];

    GetSPBunchRMS(latticeInterActionPoint, 0);

}

void SPBunch::GetSPBunchRMS(const LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    if (eSurive[0] != 0)
    {
        cerr<<"electron beam is lost at certain bunch in SP model "<<eSurive[0] <<endl;
        exit(0);
    } 
    
    rmsRx = sqrt(emittanceX * latticeInterActionPoint.twissBetaX[k] +  pow(rmsEnergySpread *latticeInterActionPoint.twissDispX[k],2) );
    rmsRy = sqrt(emittanceY * latticeInterActionPoint.twissBetaY[k] +  pow(rmsEnergySpread *latticeInterActionPoint.twissDispY[k],2) );

    xAver  = ePositionX[0];
    yAver  = ePositionY[0];
    zAver  = ePositionZ[0];
    pxAver = eMomentumX[0];
    pyAver = eMomentumY[0];
    pzAver = eMomentumZ[0];

    actionJx = 0.E0;
    actionJy = 0.E0;

    double twissAlphaXTemp = latticeInterActionPoint.twissAlphaX[k];
    double twissAlphaYTemp = latticeInterActionPoint.twissAlphaY[k];
    double twissBetaXTemp  = latticeInterActionPoint.twissBetaX[k];
    double twissBetaYTemp  = latticeInterActionPoint.twissBetaY[k];

    actionJx = (1+pow(twissAlphaXTemp,2))/twissBetaXTemp * pow(xAver,2)
             +  2*twissAlphaXTemp* xAver * pxAver
             +   twissBetaXTemp * pow(pxAver,2);

    actionJy = (1+pow(twissAlphaYTemp,2))/twissBetaYTemp * pow(yAver,2)
            + 2*twissAlphaYTemp* yAver * pyAver
            +   twissBetaYTemp * pow(pyAver,2);

    actionJx = actionJx /2.;
    actionJy = actionJy /2.;
    
}

void SPBunch::WSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    double nE0 = electronNumPerBunch;
    double nI0;
    double ionMassNumber;
    double omegaE, omegaI;
    double coeffI, coeffE;
    double tempFx=0.E0;
    double tempFy=0.E0;
    double posx=0.E0;
    double posy=0.E0;
    double rmsRxTemp = rmsRx;
    double rmsRyTemp = rmsRy;
    double eFxTemp;
    double eFyTemp;

    eFxDueToIon[0] =0.E0;
    eFyDueToIon[0] =0.E0;
    eFxTemp = 0.E0;
    eFyTemp = 0.E0;


    for(int p=0;p<latticeInterActionPoint.gasSpec;p++)
    {
        nI0 = latticeInterActionPoint.macroIonCharge[k][p];
        ionMassNumber = latticeInterActionPoint.ionMassNumber[p];

        coeffI = 2.0*nE0*ElecClassicRadius*ElectronMassEV/IonMassEV/ionMassNumber * CLight;   // [m * m/s]
        coeffE = 2.0*ElecClassicRadius/rGamma;                                                // [m]

        for(int j=0;j<latticeInterActionPoint.ionAccumuNumber[k][p];j++)
        {
            latticeInterActionPoint.ionAccumuFx[k][p][j]=0.E0;
            latticeInterActionPoint.ionAccumuFy[k][p][j]=0.E0;

            posx    =  latticeInterActionPoint.ionAccumuPositionX[k][p][j] - xAver;
            posy    =  latticeInterActionPoint.ionAccumuPositionY[k][p][j] - yAver;

            if(abs(posx/rmsRxTemp) + abs(posy/rmsRyTemp)<1.0E-5)
            {
                tempFx=0;
                tempFy=0;
            }
            else
            {
                if( (rmsRxTemp - rmsRyTemp)/rmsRyTemp > 1.e-4 )
                {
                    BassettiErskine1(posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy);         // tempFx ->[1/m]; // assume electron rms as gaussion-- get force at ion position
                }
                else if ( (rmsRxTemp - rmsRyTemp)/rmsRyTemp < -1.e-4  )
                {
                    BassettiErskine1(posy,posx,rmsRyTemp,rmsRxTemp,tempFy,tempFx);
                }
                else
                {
                    GaussianField(posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy);
                }
            }

            latticeInterActionPoint.ionAccumuFx[k][p][j]= coeffI*tempFx;                  //[1/m] * [m * m/s] - > [m/s]; -> (13) integrate along dt gives ion velovity change
            latticeInterActionPoint.ionAccumuFy[k][p][j]= coeffI*tempFy;

            eFxTemp +=  tempFx;
            eFyTemp +=  tempFy;
        }

        eFxDueToIon[0] += -1*eFxTemp * coeffE * nI0;                                      // since the e and ion with opposite charge state  [1/m * m ]-> [rad]  integrage (12) along ds -- beam dpx change.
        eFyDueToIon[0] += -1*eFyTemp * coeffE * nI0;

    }

    latticeInterActionPoint.GetTotIonCharge();
    totIonCharge = latticeInterActionPoint.totIonCharge;

}


void SPBunch::BunchTransferDueToLatticeL(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    // Ref. bunch.h that ePositionZ = - ePositionT * c. head pariticles: deltaT<0, ePositionZ[i]>0.
    // During the tracking, from head to tail means ePositionZMin from [+,-]; 
    
    // Notification: thanks to info exchange with Naoto, Yamamoto. In below, the equation to get cPsi,
    // In original PB Wilson and N.G. Bill's paper,
    // cPsi = 2 * PI * (cavityResonator.resonatorVec[j].resFre - ringHarmH * f0 * resHarm) * tB;     
    // It is the picture in the phasor frame and can not apply to real frame. 
    // In instead, the equaiton have to be modified as, which goes into the real frame.
    // cPsi = 2 * PI * cavityResonator.resonatorVec[j].resFre * tB;
    // Then the instability can be excited... 


    // in the longitudinal tracking, the (z,pz) rotate in phase with anti-colckwise
    double *synchRadDampTime = inputParameter.ringParBasic->synchRadDampTime;

    int resNum        = inputParameter.ringParRf->resNum;
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double t0         = inputParameter.ringParBasic->t0;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double eta        = inputParameter.ringParBasic->eta;
    double u0         = inputParameter.ringParBasic->u0;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    double fRF         = f0 * ringHarmH;
  
    double tB = 0.E0;
    double deltaL = 0.E0;
    double cPsi   = 0.E0;
    double resFre       = 0.E0;
    int resHarm = 0; 
    
    complex<double> cavVoltage =(0.E0,0.E0);
    complex<double> genVoltage =(0.E0,0.E0);
    complex<double> vb0=(0,0);

 

    // loop of cavities. 
    for(int j=0;j<resNum;j++)
    {                              
        resHarm = cavityResonator.resonatorVec[j].resHarm;
        resFre  = cavityResonator.resonatorVec[j].resFre; 

        //-----------------------------------------------------------------------------------------------------------------------------------------
        vb0  = complex<double>( -1 * cavityResonator.resonatorVec[j].resFre * 2 * PI * cavityResonator.resonatorVec[j].resShuntImpRs /
                                     cavityResonator.resonatorVec[j].resQualityQ0, 0.E0) * electronNumPerBunch * ElectronCharge;       //[Volt]

        //Elegant method    
        // genVoltage = cavityResonator.resonatorVec[j].resGenVol * exp( - li * ePositionZ[0] / CLight * 2. * PI * double(resHarm) * fRF );        
        // cavVoltage = cavityResonator.resonatorVec[j].vbAccum + genVoltage;
        // eMomentumZ[0] += cavVoltage.real()  /electronBeamEnergy / pow(rBeta,2) ;
        // eMomentumZ[0] += vb0.real()/2.0     /electronBeamEnergy / pow(rBeta,2) ;
        // cavityResonator.resonatorVec[j].vbAccum += vb0;  
        //---------------------------------------------------------------------------------------------------

        //Yamamoto;s approaches, can be used to get stable solution
        genVoltage = cavityResonator.resonatorVec[j].resGenVol; 
        cavVoltage = cavityResonator.resonatorVec[j].vbAccum + cavityResonator.resonatorVec[j].resGenVol; 
        eMomentumZ[0] += (cavVoltage * exp(  - li * ePositionZ[0] / CLight / rBeta * 2. * PI * double(resHarm) * fRF ) ).real() / electronBeamEnergy / pow(rBeta,2);
        eMomentumZ[0] += vb0.real()/2.0     / electronBeamEnergy / pow(rBeta,2);
        cavityResonator.resonatorVec[j].vbAccum += vb0; 
        //---------------------------------------------------------------------------------------------------


        
        // get the beam induced voltage sampled at the time t = m*Trf according to the ePositionZ[0] = - c dt[0],
        // if dt[0]>0, or ePositionZ[0]<0,  phasor grow and rotate inversely to t = m*Trf 
        // if dt[0]<0, or ePositionZ[0]>0,  show grow and rotate as it is
        // for main cavity, -- opposite direction 
        // head particle, dt[0]<0, or ePositionZ[0]>0, dtTemp < 0, deltaL<0: growth; cPsi<0,  rotate clockwise 
        // tail particle, dt[0]>0, or ePositionZ[0]<0, dtTemp > 0, deltaL>0: decay;  cPsi>0,  rotate anticlockwise, just as cavity evolove extra dt[0]                
        cavityResonator.resonatorVec[j].vBSampleTemp = cavityResonator.resonatorVec[j].vbAccum;  
        double dtTemp = - ePositionZ[0] / CLight / rBeta; 
        deltaL = dtTemp / cavityResonator.resonatorVec[j].tF;        
        cPsi   = 2.0 * PI *  cavityResonator.resonatorVec[j].resFre * dtTemp;      
        cavityResonator.resonatorVec[j].vBSampleTemp *=  exp( - deltaL ) * exp (li * cPsi);        
        //---------------------------------------------------------------------------------------------------



        // set the to pinrt, what bunch feels.
        bunchRFModeInfo->induceVolBunchCen[j]   = cavityResonator.resonatorVec[j].vbAccum;
        bunchRFModeInfo->genVolBunchAver[j]     = genVoltage;               
        bunchRFModeInfo->cavVolBunchCen[j]      = cavVoltage;
        bunchRFModeInfo->selfLossVolBunchCen[j] = vb0/2.0;
        ///////////////////////////////////////////////////////////////////////
        
        if(cavityResonator.resonatorVec[j].resExciteIntability==0)
        {
            tB = bunchGap * 1.0 / f0 / ringHarmH;    // instablity is not exited
        }
        else
        {
            tB    = timeFromCurrnetBunchToNextBunch;  // instablity is exited during beam loading simulation
        }
     
        deltaL = tB / cavityResonator.resonatorVec[j].tF;        
        cPsi   = 2.0 * PI *  cavityResonator.resonatorVec[j].resFre * tB; // - 2 * PI * resHarm * bunchGap;      
        cavityResonator.resonatorVec[j].vbAccum *= exp( - deltaL ) * exp (li * cPsi); 
        //-----------------------------------------------------------------------------------------------------------------------------------------
    }
       

    eMomentumZ[0] -= u0 / electronBeamEnergy / pow(rBeta,2);
    if(inputParameter.ringRun->synRadDampingFlag[1]==1)
    {
        eMomentumZ[0] *= (1.0-2.0/synchRadDampTime[2]);                
    }
    ePositionZ[0]  -= eta * t0  * CLight  * eMomentumZ[0];    // deltaZ = -deltaT * CLight = -eta * T0 * deltaPOverP * CLight   
}




// for this one , it is used for coupled bunch and single bunch growth rate calculaiton, the stastic potential well term is removed 
void SPBunch:: BunchTransferDueToLatticeLMatarix(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    // in the longitudinal tracking, the (z,pz) rotate in phase with anti-colckwise
    double *synchRadDampTime = inputParameter.ringParBasic->synchRadDampTime;

    int resNum        = inputParameter.ringParRf->resNum;
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double t0         = inputParameter.ringParBasic->t0;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double eta        = inputParameter.ringParBasic->eta;
    double u0         = inputParameter.ringParBasic->u0;
    double alphac     = inputParameter.ringParBasic->alphac;
    double workQz     = inputParameter.ringParBasic->workQz;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    double circRing = inputParameter.ringParBasic->circRing;
    double rfbucketLen = CLight * t0 / ringHarmH;
    double fRF         = f0 * ringHarmH;

    // matarix transfer in longitudinal
    // double phaseAdvanceZ = -2 * PI * workQz;    //in longitudinal particle rotate anti-clockwise  in phase space
    // double betaz  = inputParameter.ringBunchPara->rmsBunchLength / inputParameter.ringBunchPara->rmsEnergySpread; 

    // double a00 =              cos(phaseAdvanceZ);
    // double a01 =      betaz * sin(phaseAdvanceZ);
    // double a10 = -1 / betaz * sin(phaseAdvanceZ);
    // double a11 =              cos(phaseAdvanceZ);
    // double tempz  = ePositionZ[0];
    // double temppz = eMomentumZ[0];
    // ePositionZ[0]  =  a00 * tempz + a01 * temppz;
    // eMomentumZ[0]  =  a10 * tempz + a11 * temppz; 
    //--------------------------------------------------------------------------------------------------------


    // here for coupled bunc benchmark with phasor notation.
    // in below to benchmark with analytical growth rate prediction -- and stastic potential well has to be substracted. 
    double tF = 0.E0;
    double tB = 0.E0;
    double tB0 = 0.E0;
    double deltaL = 0.E0;
    double cPsi   = 0.E0;
    double resGenVolAmp = 0.E0;
    double resGenVolArg = 0.E0;
    double resFre       = 0.E0;
    double resPhaseReq;
    double resVoltageReq;
    int resHarm = 0;
    
    double cavVolAmp;
    double cavVolArg;
    complex<double> cavVoltage =(0.E0,0.E0);
    complex<double> vb0=(0,0);
    double temp = 0;

    for(int j=0;j<resNum;j++)
    {                       
        resHarm = cavityResonator.resonatorVec[j].resHarm;
        resFre  = cavityResonator.resonatorVec[j].resFre; 

        vb0  = complex<double>( -1 * cavityResonator.resonatorVec[j].resFre * 2 * PI * cavityResonator.resonatorVec[j].resShuntImpRs /
                                     cavityResonator.resonatorVec[j].resQualityQ0, 0.E0) * electronNumPerBunch * ElectronCharge;       //[Volt]
        
        cavVoltage = cavityResonator.resonatorVec[j].vbAccum + cavityResonator.resonatorVec[j].resGenVol;                  
        eMomentumZ[0] += (cavVoltage * exp( - li * ePositionZ[0] / CLight * 2. * PI * double(resHarm) * fRF )/ electronBeamEnergy / pow(rBeta,2) ).real(); 

        tB    = timeFromCurrnetBunchToNextBunch;  // instablity is exited during beam loading simulation                
        deltaL = tB / cavityResonator.resonatorVec[j].tF;        
        cPsi   = 2.0 * PI *  cavityResonator.resonatorVec[j].resFre * tB; // - 2 * PI * resHarm * bunchGap;
        cavityResonator.resonatorVec[j].vbAccum += vb0;  
        cavityResonator.resonatorVec[j].vbAccum *= exp(- deltaL ) * exp( li * cPsi);
      
        cavVoltage = cavityResonator.resonatorVec[j].vbAccum0 + cavityResonator.resonatorVec[j].resGenVol;  
        eMomentumZ[0] -= (cavVoltage * exp( - li * ePositionZ[0] / CLight * 2. * PI * double(resHarm) * fRF )/ electronBeamEnergy / pow(rBeta,2)).real(); 
        tB  = bunchGap * 1.0 / f0 / ringHarmH; 
        deltaL = tB / cavityResonator.resonatorVec[j].tF; 
        cPsi   = 2.0 * PI *  cavityResonator.resonatorVec[j].resFre * tB;
        cavityResonator.resonatorVec[j].vbAccum0 += vb0;  
        cavityResonator.resonatorVec[j].vbAccum0 *= exp(- deltaL ) * exp( li * cPsi);

    }
    //--------------------------------------------------------------------------------------------------------

    // leapfrog integration process.  Alex Chao 4.2 and 4.1 -- same as SYL 3.42 
    // (1) get the f(q) at the n step
    double fqn = pow(2 * PI *  workQz ,2) /  eta / circRing *  ePositionZ[0] ;
    // (2) get  p at n + 1/2 step     
    eMomentumZ[0] +=  fqn / 2; 
    // (3) get  q at n + 1 step
    ePositionZ[0] -= eta * circRing * eMomentumZ[0];
    // (4) get the f(q) at n + 1 step 
    fqn = pow(2 * PI *  workQz ,2) /  eta / circRing *  ePositionZ[0];
    // (5) get the p at n + 1 step
    eMomentumZ[0] +=  fqn / 2; 
}


// retate reference to the n*trf frame, in which the generator voltate always rotate n*2*pi.  
// In simulation it is not a right frame  for simulation 
void SPBunch::BunchTransferDueToLatticeLYamamoto(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    // in the longitudinal tracking, the (z,pz) rotate in phase with anti-colckwise
    double *synchRadDampTime = inputParameter.ringParBasic->synchRadDampTime;

    int resNum        = inputParameter.ringParRf->resNum;
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double t0         = inputParameter.ringParBasic->t0;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double eta        = inputParameter.ringParBasic->eta;
    double u0         = inputParameter.ringParBasic->u0;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    double fRF         = f0 * ringHarmH;

    double tF = 0.E0;
    double tB = 0.E0;
    double deltaL = 0.E0;
    double cPsi   = 0.E0;
    double resFre = 0.E0;
    int resHarm = 0;
    
    complex<double> cavVoltage =(0.E0,0.E0);
    complex<double> genVoltage =(0.E0,0.E0);
    complex<double> vb0=(0,0);

    // loop of cavities. 
    // Get Yamomoto's approaches to get his fig.10 -------------------------------------------------
    for(int j=0;j<resNum;j++)
    {                       
        resHarm = cavityResonator.resonatorVec[j].resHarm;
        resFre  = cavityResonator.resonatorVec[j].resFre;

        vb0  = complex<double>( -1 * cavityResonator.resonatorVec[j].resFre * 2 * PI * cavityResonator.resonatorVec[j].resShuntImpRs /
                                     cavityResonator.resonatorVec[j].resQualityQ0, 0.E0) * electronNumPerBunch * ElectronCharge;       //[Volt]
    
        cavVoltage = cavityResonator.resonatorVec[j].vbAccum + cavityResonator.resonatorVec[j].resGenVol; 
        eMomentumZ[0] += (cavVoltage * exp(  - li * ePositionZ[0] / CLight / rBeta * 2. * PI * double(resHarm) * fRF ) ).real() / electronBeamEnergy / pow(rBeta,2);
        eMomentumZ[0] += vb0.real()/2.0     / electronBeamEnergy / pow(rBeta,2);
        cavityResonator.resonatorVec[j].vbAccum += vb0;  
        
        // set the data used for feedback 
        bunchRFModeInfo->induceVolBunchCen[j]   = cavityResonator.resonatorVec[j].vbAccum;                     
        bunchRFModeInfo->cavVolBunchCen[j]      = cavVoltage;
        bunchRFModeInfo->selfLossVolBunchCen[j] = vb0/2.0;   
        //--------------------------------------------------------------------------------------------------


        // time to next bunch
        if(cavityResonator.resonatorVec[j].resExciteIntability==0)
        {
            tB = bunchGap * 1.0 / f0 / ringHarmH;    // instablity is not exited
        }
        else
        {
            tB    = timeFromCurrnetBunchToNextBunch;  // instablity is exited during beam loading simulation
        }
                
        deltaL = tB / cavityResonator.resonatorVec[j].tF;        
        cPsi   = 2.0 * PI *  cavityResonator.resonatorVec[j].resFre * tB;     
        cavityResonator.resonatorVec[j].vbAccum *= exp(- deltaL ) * exp(li * cPsi); 
    }
    
    eMomentumZ[0] -= u0 / electronBeamEnergy / pow(rBeta,2);
    if(inputParameter.ringRun->synRadDampingFlag[1]==1)
    {
        eMomentumZ[0] *= (1.0-2.0/synchRadDampTime[2]);                
    }
    ePositionZ[0]  -= eta * t0  * CLight  * eMomentumZ[0];             // deltaZ = - deltaT * CLight = - eta * T0 * deltaPOverP * CLight 
   
}

