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
    delete haissinski;
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
    

    ePositionZ[0] +=  disDz;
    eMomentumZ[0] +=  disMz;

  	double dispersionX  =  latticeInterActionPoint.twissDispX[0];
  	double dispersionY  =  latticeInterActionPoint.twissDispY[0];
  	double dispersionPX =  latticeInterActionPoint.twissDispPX[0];
  	double dispersionPY =  latticeInterActionPoint.twissDispPY[0];

	ePositionX[0] +=  disDx + dispersionX  * eMomentumZ[0];
	ePositionY[0] +=  disDy + dispersionY  * eMomentumZ[0];
	eMomentumX[0] +=  disMx + dispersionPX * eMomentumZ[0];
	eMomentumY[0] +=  disMy + dispersionPY * eMomentumZ[0];


    // set the initial coupled bunch mode according mode index /mu 
    
    //int N = inputParameter.ringFillPatt->totBunchNumber; 
    
    //double alphaX=latticeInterActionPoint.twissAlphaX[0];
    //double alphaY=latticeInterActionPoint.twissAlphaY[0];
    //double betaX =latticeInterActionPoint.twissBetaX[0];
    //double betaY =latticeInterActionPoint.twissBetaY[0];
    //int mu = 3;
    //double phase = 2 * PI * mu * bunchIndex / N;

    //
    //ePositionX[0] =  sqrt(emittanceX * betaX) * cos(phase) ;
    //ePositionY[0] =  sqrt(emittanceY * betaY) * cos(phase) ;
    //eMomentumX[0] = -sqrt(emittanceX / betaX) * (alphaX * cos(phase) + sin(phase) );
    //eMomentumY[0] = -sqrt(emittanceY / betaY) * (alphaY * cos(phase) + sin(phase) );



}

void SPBunch::GetSPBunchRMS(const LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    if (eSurive[0] = 0)
    {
        cerr<<"electron beam is lost at certain bunch in SP model"<<endl;
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
    
    double *synchRadDampTime = inputParameter.ringParBasic->synchRadDampTime;

    int resNum        = inputParameter.ringParRf->resNum;
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double t0         = inputParameter.ringParBasic->t0;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double eta        = inputParameter.ringParBasic->eta;
    double u0         = inputParameter.ringParBasic->u0;
    double alphac    = inputParameter.ringParBasic->alphac;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;

    complex<double> vb0=(0,0);

    double tF = 0.E0;
    double tB = 0.E0;
    double deltaL = 0.E0;
    double cPsi   = 0.E0;
    double resGenVolAmp = 0.E0;
    double resGenVolArg = 0.E0;
    double resFre       = 0.E0;
    double resPhaseReq;
    int resHarm = 0;
    
    double cavVolAmp;
    double cavVolArg;
    complex<double> cavVoltage=(0.E0,0.E0);

    for(int j=0;j<resNum;j++)
    {
        resHarm = cavityResonator.resonatorVec[j].resHarm;
        resFre  = cavityResonator.resonatorVec[j].resFre;
        resPhaseReq = cavityResonator.resonatorVec[j].resPhaseReq;

        vb0  = complex<double>( -1 * cavityResonator.resonatorVec[j].resFre * 2 * PI * cavityResonator.resonatorVec[j].resShuntImpRs /
                                     cavityResonator.resonatorVec[j].resQualityQ0, 0.E0) * electronNumPerBunch * ElectronCharge;            //[Volt]

        // induced voltage vbAccum is at time just before the next bunch comes,  genVolBunchAver[j] includes average one turn feedback.
        cavVoltage = cavityResonator.resonatorVec[j].vbAccum + cavFBCenInfo->genVolBunchAver[j];
  
        //cavVoltage = cavityResonator.resonatorVec[j].resCavVolReq;
        cavVolAmp  = abs(cavVoltage);
        cavVolArg  = arg(cavVoltage);

        haissinski->cavAmp[j]   = cavVolAmp;
        haissinski->cavPhase[j] = cavVolArg;

        //particle momentum due to cavity
        eMomentumZ[0] += cavVolAmp  * cos( -2. * PI * ringHarmH * f0 * resHarm * ePositionZ[0] / (rBeta * CLight)  +  cavVolArg )/electronBeamEnergy;

        //particle momentum due to self-field
        eMomentumZ[0] += vb0.real()/2.0/electronBeamEnergy;

        // beam induced voltage accumulated  
        cavityResonator.resonatorVec[j].vbAccum += vb0;

        //tB = timeToFirstBinOfNextBunch;
        tB = bunchGap * t0 / ringHarmH;
        deltaL = tB / cavityResonator.resonatorVec[j].tF;
        cPsi   = deltaL * tan(cavityResonator.resonatorVec[j].resDeTunePsi);


        // beam induced voltage rotate and decay
        cavityResonator.resonatorVec[j].vbAccum *= exp(- deltaL ) * exp (li * cPsi);   

        // get the information used in cavity feedback, get the average cavity voltage, beam induced vol and selflose voltage of this bunch.
        cavFBCenInfo->induceVolBunchCen[j]   = cavityResonator.resonatorVec[j].vbAccum ;            
        cavFBCenInfo->cavAmpBunchCen[j]      = cavVolAmp;
        cavFBCenInfo->cavPhaseBunchCen[j]    = cavVolArg;
        cavFBCenInfo->cavVolBunchCen[j]      = cavVoltage;
        cavFBCenInfo->selfLossVolBunchCen[j] = vb0/2.0;

    
    }
    eMomentumZ[0]  = eMomentumZ[0] - u0 / electronBeamEnergy;
    eMomentumZ[0]  = eMomentumZ[0] * (1.0-2.0/synchRadDampTime[2]);                   // to maintain the symetric,
    ePositionZ[0] -= eta * t0  * CLight * rBeta * eMomentumZ[0];                      // deltaZ = -deltaT * CLight = -eta * T0 * deltaPOverP * CLight         
}

void SPBunch::GetBunchHaissinski(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,WakeFunction &sRWakeFunction)
{         
    // Head to tail means bunchPosZ from [+,-]; 
    //cout<<"get the Haissinski solution for the specified bunch"<<endl;
    int ringHarm          = inputParameter.ringParRf->ringHarm;
    double f0             = inputParameter.ringParBasic->f0;
    double sigmaT0        = inputParameter.ringParBasic->sigmaT0 ;
    double sigmaZ0        = sigmaT0 * CLight;   

    int nz                = sigmaT0 * 50 / haissinski->dt;      // +- 50 rms bunch size as the range for haissinksi binsize=1 ps here used as default setting
                                                                // usually is enougth, double RF and impedance increase bunch at most 10 times: sigmaT0 -> 10 * sigmaT0    
   
    haissinski->zMax      =   nz * haissinski->dt * CLight;
    haissinski->zMin      =  -haissinski->zMax;    
    haissinski->nz        =   2 * nz + 1;    
    nz                    = haissinski->nz; 
                      
    haissinski->bunchPosZ.resize(haissinski->nz);                            // (-,+)-> tail to head   
    haissinski->bunchProfile.resize(haissinski->nz);                         // [1/m] 
    haissinski->totWakePoten.resize(haissinski->nz);                         // [1/m]
    haissinski->rwWakePoten.resize(haissinski->nz);
    haissinski->bbrWakePoten.resize(haissinski->nz);
    haissinski->wakeHamiltonian.resize(haissinski->nz);  
    haissinski->rfHamiltonian.resize(haissinski->nz);
    haissinski->totHamiltonian.resize(haissinski->nz);                      // dimensionless

	double coef=electronNumPerBunch * 1.0/sqrt(2*PI)/sigmaZ0;
	double temp=0;		


    for(int i=0;i<haissinski->nz;i++)
    {        
        haissinski->bunchPosZ[i]    = i * haissinski->dz + haissinski->zMin;
        haissinski->bunchProfile[i] = coef * exp(-pow(haissinski->bunchPosZ[i],2)/2/pow(sigmaZ0,2));   //[electroNumber/m]~[1/m]                
    }

    // double loop to find out the cavity phase, short range wakes from the haissinski solution 
    vector<double> bunchProfile1(haissinski->nz,0);      
    for(int k=0;k<1001;k++)     
    {
        //cout<<"haissinski loop: " <<k<<endl;
        GetRFHamiltonian(inputParameter,cavityResonator);     

        if(inputParameter.ringRun->sRWakeFlag)
        {
            GetWakeHamiltonian(inputParameter,sRWakeFunction);  
        }
          
        GetTotHamiltonian();
        bunchProfile1=GetProfile(inputParameter);
        
        double temp=0;
        for(int i=0;i<haissinski->nz;i++)
        {
            haissinski->bunchProfile[i] = haissinski->bunchProfile[i] *0.9 + bunchProfile1[i] * 0.1;            
            temp += abs(haissinski->bunchProfile[i] - bunchProfile1[i]) * haissinski->dz;
        }        
                    
        if (temp/electronNumPerBunch < 1.e-9)        
        {
            break;
        }              
    }
    
    GetBunchAverAndBunchLengthFromHaissinskiSolution();
    
    cout<<"haissinki, bunch.bucketIndex:="<<setw(15)<<left<<bunchHarmNum 
        <<"Bunch Center (m): "            <<setw(15)<<left<<haissinski->averZ
        <<"Bunch Length (m): "            <<setw(15)<<left<<haissinski->rmsZ
        <<endl; 
        
}

void SPBunch::GetBunchAverAndBunchLengthFromHaissinskiSolution()
{
    double temp=0;
    double norm=0;
    for(int i=0;i<haissinski->nz;i++)
    {
        temp += haissinski->bunchProfile[i] * haissinski->dz * haissinski->bunchPosZ[i]  ;
        norm += haissinski->bunchProfile[i] * haissinski->dz ;
    }
    
    haissinski->averZ = temp / norm;
    
    temp = 0;
    
    for(int i=0;i<haissinski->nz;i++)
    {
        temp += haissinski->bunchProfile[i] * haissinski->dz * pow(haissinski->bunchPosZ[i] - haissinski->averZ, 2);
    }
    
    haissinski->rmsZ = sqrt(temp / norm);
    
}


vector<double> SPBunch::GetProfile(const ReadInputSettings &inputParameter)
{
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double f0         = inputParameter.ringParBasic->f0;
    double eta        = inputParameter.ringParBasic->eta;
    double sdelta0    = inputParameter.ringParBasic->sdelta0;

    double coef = 2 * PI * ringHarmH * f0 * eta * pow(sdelta0,2);  //  [1/s]


    vector<double> profile(haissinski->nz,0);
    double norm=0;
    for(int i=0;i<haissinski->nz;i++)
    {
        profile[i] = exp(-haissinski->totHamiltonian[i]/coef);
        norm += profile[i] * haissinski->dz;
    }
    
    for(int i=0;i<haissinski->nz;i++)
    {
        profile[i] =  profile[i] / norm * electronNumPerBunch; 
    }
    /*
    ofstream fout("test.dat");
    for(int i=0;i<haissinski->nz;i++)
    {    
        fout<<setw(15)<<left<< haissinski->bunchPosZ[i] 
            <<setw(15)<<left<< haissinski->totWakePoten[i] 
            <<setw(15)<<left<< haissinski->wakeHamiltonian[i]
            <<setw(15)<<left<< haissinski->rfHamiltonian[i] 
            <<setw(15)<<left<< haissinski->totHamiltonian[i]
            <<setw(15)<<left<<profile[i]
            <<endl;
    }
    cout<<__LINE__<<endl;
    getchar();
    */
    return profile;    
        
} 

void SPBunch::GetTotHamiltonian()
{    
    
    for(int i=0;i<haissinski->nz;i++)
    {
        haissinski->totHamiltonian[i] = haissinski->rfHamiltonian[i] + haissinski->wakeHamiltonian[i];
        haissinski->totWakePoten[i]   = haissinski->bbrWakePoten[i]  + haissinski->rwWakePoten[i];
    }
                   
    auto iter0 = min_element(haissinski->totHamiltonian.begin(),haissinski->totHamiltonian.end() );
    double minTemp = *iter0;
    for(int i=0;i<haissinski->nz;i++)
    {
        haissinski->totHamiltonian[i]  -=  minTemp;
        haissinski->rfHamiltonian[i]   -=  minTemp;        
    }

}

void SPBunch::GetRFHamiltonian(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator)
{
    vector<double> vTotRF(haissinski->nz,0);  

    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double t0         = inputParameter.ringParBasic->t0;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double eta        = inputParameter.ringParBasic->eta;
    double u0         = inputParameter.ringParBasic->u0;
    double alphac     = inputParameter.ringParBasic->alphac;
    double sdelta0    = inputParameter.ringParBasic->sdelta0;
    double circRing   = inputParameter.ringParBasic->circRing;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy; 
    int resHarm;  

    double vb0;
     
    // the rf voltage and phase are the same as initialized for haissinski.  
    for(int i=0;i<haissinski->nz;i++)
    {
        for(int j=0;j<cavityResonator.resonatorVec.size();j++)
        {
            vb0  =  -1 * cavityResonator.resonatorVec[j].resFre * 2 * PI * cavityResonator.resonatorVec[j].resShuntImpRs /
                         cavityResonator.resonatorVec[j].resQualityQ0 * haissinski->dz * haissinski->bunchProfile[i] * ElectronCharge;     // [Volt]
            
            resHarm = cavityResonator.resonatorVec[j].resHarm;   
            vTotRF[i] += haissinski->cavAmp[j] * cos( - 2. * PI * ringHarmH * f0 * resHarm * haissinski->bunchPosZ[i] /  ( rBeta * CLight)  + haissinski->cavPhase[j] ) 
                       + vb0 / 2 ; // [Volt]
        }  
    }    
    
    double coeffDelta = f0 / pow(rBeta,2) / electronBeamEnergy;   // 1/[V]/[s];

    //Eq.(2~8), PRAB 17 064401, and S. Y. Lee Eq. (3.35 and 3.36)  ddelta/dt = - dH/dphi  and (dz=- c * dt = - c * dphi / (2 * PI * h * f0) )
    haissinski->rfHamiltonian[0] = 0.e0;
    for(int i=1;i<haissinski->nz;i++)
    {              
        haissinski->rfHamiltonian[i] = haissinski->rfHamiltonian[i-1] - ((vTotRF[i] + vTotRF[i-1]) / 2.0 -  u0 ) * coeffDelta 
                                     * (-1) *  2. * PI * ringHarmH * f0 * haissinski->dz / ( rBeta * CLight);                        // [1/s]          
    }    
    
}

void SPBunch::GetWakeHamiltonian(const ReadInputSettings &inputParameter, WakeFunction &sRWakeFunction)
{
    for(int i=0;i<haissinski->nz;i++)
    {
        haissinski->wakeHamiltonian[i] = 0.E0; 
        haissinski->rwWakePoten[i]     = 0.E0;
        haissinski->bbrWakePoten[i]    = 0.E0;         
    }
     
    if(!inputParameter.ringSRWake->bbrInput.empty())
    {   
        GetWakeHamiltonianFromBBR(inputParameter,sRWakeFunction); 
    }    
    
    if(!inputParameter.ringSRWake->pipeGeoInput.empty())
    {   
        GetWakeHamiltonianFromRW(inputParameter,sRWakeFunction);    
    }
      
}


void SPBunch::GetWakeHamiltonianFromBBR(const ReadInputSettings &inputParameter, WakeFunction &sRWakeFunction)
{
    // Head to tail means bunchPosZ from [+,-]; 

    vector<double> wakeHamiltonTemp(haissinski->nz,0.E0);
    vector<double> wakePotenTemp(haissinski->nz,0.E0);
    
    vector<double> wakeFunji;
    double tauji;   
    
    for(int i=0;i<haissinski->nz;i++)
    {
        for(int j=0;j<haissinski->nz;j++)
        {
            tauji = (i-j) * haissinski->dz / CLight;                      // [s]  must .LE. 0             
            wakeFunji = sRWakeFunction.GetBBRWakeFun1(tauji);              // [V/C]                                                    
            wakePotenTemp[i] += wakeFunji[2] *  haissinski->bunchProfile[j] * haissinski->dz ;  // [V/C] * [1/m] * [m]  = [V/C]                                       
        }
        wakePotenTemp[i] = -wakePotenTemp[i] * ElectronCharge;                           // [V]  Eq. (3.7) -- multiplty -1;  
    }
    
    
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;  

    double coeffDelta = f0 / pow(rBeta,2) / electronBeamEnergy;   // 1/[V]/[s];

    // get the Hamiltonian due to wakePoten.   dp/dt = - dH/dq  
    wakeHamiltonTemp[0] = 0;
    for(int i=1;i<haissinski->nz;i++)
    {
        wakeHamiltonTemp[i] = wakeHamiltonTemp[i-1] - ( wakePotenTemp[i] + wakePotenTemp[i-1] ) / 2.0   * coeffDelta 
                             *  (-1) * 2. * PI * ringHarmH * f0 * haissinski->dz / ( rBeta * CLight);                       //  [V]  *  1/[V]/[s] = 1/[s]         
    }

    for(int i=0;i<haissinski->nz;i++)
    {
        haissinski->wakeHamiltonian[i] += wakeHamiltonTemp[i];      //[1/s]
        haissinski->bbrWakePoten[i]    += wakePotenTemp[i];         //[V]
    }                 
}



void SPBunch::GetWakeHamiltonianFromRW(const ReadInputSettings &inputParameter, WakeFunction &sRWakeFunction)
{
    vector<double> wakeHamiltonTemp(haissinski->nz,0.E0);            
    vector<double> wakePotenTemp(haissinski->nz,0.E0);
    vector<double> wakeFunji(3,0);
    double tauji;   
     
    for(int i=0;i<haissinski->nz;i++)                 // with the persudo green wake function down and up limit have to change!
    {
        for(int j=0;j<haissinski->nz;j++)
        {
            tauji = (i-j) * haissinski->dz / CLight;                                     
            wakeFunji = sRWakeFunction.GetRWSRWakeFun(tauji);                                                                             
            wakePotenTemp[i] += wakeFunji[2] * haissinski->bunchProfile[j] * haissinski->dz;  // [V/C] * [1/m] * [m]  =  [V/C]                        
        }
        wakePotenTemp[i] = -wakePotenTemp[i] * ElectronCharge;        // [V]  Eq. (3.7) -- multiplty -1;     
    }


    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;  
    
    double coeffDelta = f0 / pow(rBeta,2) / electronBeamEnergy;   // 1/[V]/[s];

   // get the Hamiltonian due to wakePoten     dp/dt = - dH/dq  
    wakeHamiltonTemp[0] = 0;
    for(int i=1;i<haissinski->nz;i++)
    {
        //wakeHamiltonTemp[i] = wakeHamiltonTemp[i-1] - coef * (wakePotenTemp[i] + wakePotenTemp[i-1])/2.0 * haissinski->dz;  //[dimenionless]                          

        wakeHamiltonTemp[i] = wakeHamiltonTemp[i-1] - ( wakePotenTemp[i] + wakePotenTemp[i-1] ) / 2.0  * coeffDelta 
                             *  (-1) * 2. * PI * ringHarmH * f0 * haissinski->dz / ( rBeta * CLight);                       //  [V]  *  1/[V]/[s] = 1/[s]   
   
    }

    for(int i=0;i<haissinski->nz;i++)
    {
        haissinski->wakeHamiltonian[i] += wakeHamiltonTemp[i];
        haissinski->rwWakePoten[i]     += wakePotenTemp[i];
    } 
}



void SPBunch::GetParticleLongitudinalPhaseSpace(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,int bunchIndex)
{
    tk::spline wakePotenFit, totHamiltonianFit;
    wakePotenFit.set_points(haissinski->bunchPosZ,haissinski->totWakePoten,tk::spline::cspline);           // fitting the wakePoten 1/[m]
    totHamiltonianFit.set_points(haissinski->bunchPosZ,haissinski->totHamiltonian,tk::spline::cspline);    // fitting the totHamilton 1/[m]


    double zMax = haissinski->averZ + haissinski->rmsZ * 5;
    double zMin = haissinski->averZ - haissinski->rmsZ * 5;
    
    int nz      = 41;                               // total have 41 partilces at differet initial conditions 
    double dz   = (zMax - zMin) / nz;     
    
    double circRing   = inputParameter.ringParBasic->circRing;
    double eta        = inputParameter.ringParBasic->eta;
    double workQz     = inputParameter.ringParBasic->workQz;
    int turnsLongiOscilation = int(1/workQz);
    
    
    vector<vector<vector<double> > > longiTrajZeta;
    longiTrajZeta.resize(nz);
    vector<double> actionJ(nz,0);
    vector<double> nus(nz,0); 
    vector<double> totHamilton(nz,0); 
    // print data and calculate the hamilotnion and action J. 
    

    string filename=inputParameter.ringRun->TBTBunchLongTraj + to_string(bunchIndex) + ".sdds";
    ofstream fout(filename); 
    fout<<"SDDS1"<<endl;
    fout<<"&parameter name=z,           units=m             type=float,  &end"<<endl;
    fout<<"&parameter name=hamilton,    units=1/s           type=float,  &end"<<endl;  // S Y L 3.36
    fout<<"&parameter name=action,      units=m,            type=float,  &end"<<endl;
    fout<<"&parameter name=nus,                             type=float,  &end"<<endl;
    
    fout<<"&column name=z,              units=m,            type=float,  &end"<<endl;
    fout<<"&column name=delta,          units=rad,          type=float,  &end"<<endl;
    fout<<"&data mode=ascii, &end"<<endl;
    

    for(int i=0;i<nz;i++)
    {
        vector<double> zeta0={i*dz+zMin, 0};            
        vector<double> zeta1=zeta0;
        vector<double> zeta2(2,0);
        totHamilton[i] = totHamiltonianFit(i*dz+zMin);          //[1/s]
        double deltaS=0;
        double p0,p1;
        int counter=0;
    
        zeta2 = LeapFrog(inputParameter,cavityResonator,zeta0,wakePotenFit);
        p0 = zeta2[1];     

        for(int k=0;k<10000*turnsLongiOscilation;k++)              
        {            
            longiTrajZeta[i].push_back(zeta1);
            zeta2 = LeapFrog(inputParameter,cavityResonator,zeta1,wakePotenFit);                      
            actionJ[i] +=  abs(zeta2[0] - zeta1[0]) * abs(zeta2[1] + zeta1[1]) / 2  / (2 * PI);            
            deltaS     +=  abs(zeta2[0] - zeta1[0]) / abs(zeta2[1] + zeta1[1]) * 2  / eta ; 
            
            p1 = zeta2[1]; 
            zeta1 =zeta2;
            if(p0 * p1<0)    //ensure the longitudinal phase space only rotate one-trun.
            {
                counter +=1;
                p0=p1;
            }
            if( counter==2 ) break;                         
        }


        nus[i]     = circRing / deltaS;
        

        fout<<"! page number "<<i + 1<<endl;
        fout<<i*dz+zMin<<endl;
        fout<<totHamilton[i]<<endl;
        fout<<actionJ[i]<<endl;
        fout<<nus[i]<<endl;
        fout<<longiTrajZeta[i].size()<<endl;;

        for(int k=0;k<longiTrajZeta[i].size();k++)
        {
            fout<<setw(15)<<left<<longiTrajZeta[i][k][0]
                <<setw(15)<<left<<longiTrajZeta[i][k][1]
                <<endl;
        }
    
        
    }
    

    fout.close();
    
    
}

vector<double> SPBunch::LeapFrog(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,vector<double> zeta,const tk::spline &wakePotenFit)
{
    
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double t0         = inputParameter.ringParBasic->t0;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double eta        = inputParameter.ringParBasic->eta;
    double u0         = inputParameter.ringParBasic->u0;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;  
    int resHarm;     
    double temp;
    double wakePoten=0;
    
    
    double dT = t0;                        // 1/[s]
    double coeffDelta = f0 / pow(rBeta,2) / electronBeamEnergy;   // 1/[V]/[s];
    double coeffZ     = eta * rBeta * CLight;                     // [m]/[s]

    double q=zeta[0];
    double p=zeta[1];
    double fq;    

    // Ref. S.Y. Lee Eq. (3.35 and 3.36) -- symplectic leap-frog integration process.
    
    // leapfrog integration process.
    // (1) get the f(q) at the n step
    temp = 0;
    for(int j=0;j<cavityResonator.resonatorVec.size();j++)
    {
        resHarm = cavityResonator.resonatorVec[j].resHarm;                 
        temp   += haissinski->cavAmp[j] * cos( haissinski->cavPhase[j] - 2. * PI * ringHarmH * f0 * resHarm * q / ( rBeta * CLight));       //[V]          
    }
    temp  = temp - u0 + wakePotenFit(q) ;                                                                                                   //[v]                            
    fq = temp * coeffDelta;                                                                                                                 //1/[s]

    // (2) get  p at n + 1/2 step     
    double pHalf = p + fq * dT / 2;                                                                                                         // [rad]
    
    // (3) get  q at n + 1 step
    q = q -  coeffZ * pHalf * dT;                                                                                                           //[m]/[s] * [s] ->[m]
    
    // (4) get the f(q) at n + 1 step 
    temp = 0;
    for(int j=0;j<cavityResonator.resonatorVec.size();j++)
    {
        resHarm = cavityResonator.resonatorVec[j].resHarm;                 
        temp   +=  haissinski->cavAmp[j] * cos( haissinski->cavPhase[j] - 2. * PI * ringHarmH * f0 * resHarm * q / ( rBeta * CLight));                 
    }
    temp  = temp - u0 + wakePotenFit(q) ;                                                                                                   //[v]                           
    fq = temp * coeffDelta;                                                                                                                 //1/[s]

    // (5) get the p at n + 1 step
    p = pHalf +  fq * dT / 2;

    zeta[0]=q;
    zeta[1]=p;

    return zeta;
}
