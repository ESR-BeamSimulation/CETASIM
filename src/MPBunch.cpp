//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************                                                           
#pragma once

#include "MPBunch.h"
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
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <memory.h>


MPBunch::MPBunch()
{
}

MPBunch::~MPBunch()
{
           
}


void MPBunch::InitialMPBunch(const  ReadInputSettings &inputParameter)
{    
    //call the Initial at base class (Bunch.h)
    haissinski->cavAmp.resize(inputParameter.ringParRf->resNum,0);
    haissinski->cavPhase.resize(inputParameter.ringParRf->resNum,0); 
    Bunch::Initial(inputParameter); 

    srWakePoten.resize(3);
    int bunchBinNumberZ = inputParameter.ringParRf->rfBunchBinNum;
    
    beamCurDenZProf.resize(bunchBinNumberZ+1);
    posZBins.resize(bunchBinNumberZ);
    densProfVsBin.resize(bunchBinNumberZ+1);
    densProfVsBinAnalytical.resize(bunchBinNumberZ+1);
    hamiltonPotenWell.resize(bunchBinNumberZ+1);

    cavVolInfoVsLongBins.resize(inputParameter.ringParRf->resNum);
    cavForceInfoVsLongBins.resize(inputParameter.ringParRf->resNum);

    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        cavVolInfoVsLongBins[i].resize(bunchBinNumberZ+1);
        cavForceInfoVsLongBins[i].resize(bunchBinNumberZ+1);
    }
   
}


void MPBunch::DistriGenerator(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter, int bunchIndex)
{

// longitudial bunch phase space generatetion --simple rms in both z and z' phase space.
    double rBeta      = inputParameter.ringParBasic->rBeta;

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

    srand(time(0) + bunchIndex);
    double disDx = doffset(gen) * initialDynamicOffSet[0] + initialStaticOffSet[0];
    double disDy = doffset(gen) * initialDynamicOffSet[1] + initialStaticOffSet[1];
    double disDz = doffset(gen) * initialDynamicOffSet[2] + initialStaticOffSet[2];
    double disMx = doffset(gen) * initialDynamicOffSet[3] + initialStaticOffSet[3];
    double disMy = doffset(gen) * initialDynamicOffSet[4] + initialStaticOffSet[4];
    double disMz = doffset(gen) * initialDynamicOffSet[5] + initialStaticOffSet[5];

    std::normal_distribution<> dx{0,rmsBunchLength};
    std::normal_distribution<> dy{0,rmsEnergySpread};

    double tempx;
    double tempy;
    double temp;
    
    int i=0;
    while(i<macroEleNumPerBunch)
    {
        tempx = dx(gen);
        tempy = dy(gen);

        temp =  pow(tempx/rmsBunchLength,2) + pow(tempy/rmsEnergySpread,2);

        if( temp<pow(3,2))                                      // longitudinal is truncted at 3 sigma both in z and dp direction
        {
            ePositionZ[i] = tempx;			                    // m
            eMomentumZ[i] = tempy / (pow(rBeta,2));			    // rad dE/E ->dp/p
        }
        else
        {
            continue;
        }
        i++;
    }


    double gammaX;
    double gammaY;
    double sigmaX;
    double sigmaY;
    double dSigmaX;
    double dSigmaY;
    double alphaX;
    double alphaY;
    double betaX;
    double betaY;
  	double dispersionX;
  	double dispersionY;
  	double dispersionPX;
  	double dispersionPY;


    // bunch distribution is generated accroding the twiss parameter at the first interaction points.

    alphaX       = latticeInterActionPoint.twissAlphaX[0];
    alphaY       = latticeInterActionPoint.twissAlphaY[0];
    betaX        = latticeInterActionPoint.twissBetaX[0];
    betaY        = latticeInterActionPoint.twissBetaY[0];
  	dispersionX  = latticeInterActionPoint.twissDispX[0];
  	dispersionY  = latticeInterActionPoint.twissDispY[0];
  	dispersionPX = latticeInterActionPoint.twissDispPX[0];
  	dispersionPY = latticeInterActionPoint.twissDispPY[0];

    gammaX = (1+pow(alphaX,2))/betaX;
    gammaY = (1+pow(alphaY,2))/betaY;

    sigmaX  =  sqrt(betaX);
    sigmaY  =  sqrt(betaY);
    dSigmaX = -alphaX/sqrt(betaX);
    dSigmaY = -alphaY/sqrt(betaY);

    double rX;
    double rY;
    rX   = sqrt(emittanceX*betaX);
    rY   = sqrt(emittanceY*betaY);

    

    double f0;
    double if0;
    double fi;
    double axax;
    double ayay;
    double ax;
    double ay;
    double phaseX;
    double phaseY;
    double ix;

    int distributionType = inputParameter.ringBunchPara->distributionType;
    double kappa =  inputParameter.ringBunchPara->kappa;
    
    i=0;

    while(i<macroEleNumPerBunch)
    {
      switch(distributionType)
        {
            case 1:
                f0  = 4*emittanceX;
                if0 = double(std::rand())/RAND_MAX;
                fi = f0;
                break;
            case 2:
                f0  = 6*emittanceX;
                if0 = double(std::rand())/RAND_MAX;
                fi = f0 * sqrt(if0);
                break;
            case 3:
                f0  = 2*emittanceX;
                if0 = double(std::rand())/RAND_MAX;
                fi  = GSSlover(if0);
                fi  = fi * f0 ;
                break;
            default:
            cerr<<"Do Nothing, no distribution type in generation.  "<<endl;
        }


        ix   = double(rand())/RAND_MAX;

        axax =  fi * ix;
        ayay = (fi - axax) * kappa;

        ax  = sqrt(axax);
        ay  = sqrt(ayay);

        phaseX = 2 * PI* double(rand())/RAND_MAX;
        phaseY = 2 * PI* double(rand())/RAND_MAX;

        ePositionX[i] = ax *   sigmaX * cos( phaseX ) ;
        ePositionY[i] = ay *   sigmaY * cos( phaseY ) ;
        eMomentumX[i] = ax * (dSigmaX * cos( phaseX ) - sin(phaseX)/sigmaX );
        eMomentumY[i] = ay * (dSigmaY * cos( phaseY ) - sin(phaseY)/sigmaY );


        if(distributionType==3)
        {
            if(pow(ePositionX[i]/rX,2) + pow(ePositionY[i]/rY,2)>9)
            {
                continue;
            }
        }

        i++;
    }

	// take the dispersion and initial error into accout

	for(int i=0;i<macroEleNumPerBunch;i++)
	{
        ePositionZ[i] +=  disDz;
        eMomentumZ[i] +=  disMz;

		ePositionX[i] +=  disDx + dispersionX  * eMomentumZ[i];
		ePositionY[i] +=  disDy + dispersionY  * eMomentumZ[i];

		eMomentumX[i] +=  disMx + dispersionPX * eMomentumZ[i];
		eMomentumY[i] +=  disMy + dispersionPY * eMomentumZ[i];
	}

    GetMPBunchRMS(latticeInterActionPoint, 0);
    
    rmsBunchLengthLastTurn =  rmsBunchLength;
    zAverLastTurn          =  zAver;

    GetZMinMax();
  

}

double MPBunch::GSSlover(double if0)
{

    double fLow;
    double fUp;
    double fMid;
    double rangeLow;
    double rangeUP;
    double mid;

    rangeLow = 0.E0;
    rangeUP  = 8.E0;    // seems like 8 rms emittance truncated..

    fLow = 1 - if0 - (1+rangeLow) * exp(-rangeLow);    // compared with yuri' equation in NIMA BeamPath.
    fUp  = 1 - if0 - (1+rangeUP ) * exp(-rangeUP );


    while(rangeUP - rangeLow > 1.0E-5 )
    {
        mid  = (rangeLow + rangeUP)/2.E0;
        fMid = 1 - if0 - (1+mid) * exp(-mid);

        if(abs(fMid)<1.0E-6)
        {
            return mid;
        }

        if(fMid<0)
        {
            rangeLow = mid;
        }
        else if (fMid>0)
        {
            rangeUP = mid;
        }
    }

    mid  = (rangeLow + rangeUP)/2.E0;
    return mid;
}


void MPBunch::GetMPBunchRMS(const LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    pxAver=0.E0;
    pyAver=0.E0;
    pzAver=0.E0;
    xAver =0.E0;
    yAver =0.E0;
    zAver =0.E0;

    macroEleNumSurivePerBunch = 0;

    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        if(eSurive[i]!=0) continue;

        xAver   +=  ePositionX[i];
        yAver   +=  ePositionY[i];
        zAver   +=  ePositionZ[i];
        pxAver  +=  eMomentumX[i];
        pyAver  +=  eMomentumY[i];
        pzAver  +=  eMomentumZ[i];

        macroEleNumSurivePerBunch++;
    }

    xAver   /= macroEleNumSurivePerBunch;
    yAver   /= macroEleNumSurivePerBunch;
    zAver   /= macroEleNumSurivePerBunch;
    pxAver  /= macroEleNumSurivePerBunch;
    pyAver  /= macroEleNumSurivePerBunch;
    pzAver  /= macroEleNumSurivePerBunch;

    double x2Aver=0.E0;
    double y2Aver=0.E0;
    double z2Aver=0.E0;
    double px2Aver=0.E0;
    double py2Aver=0.E0;
    double pz2Aver=0.E0;
    double xpxAver=0.E0;
    double ypyAver=0.E0;
    double zpzAver=0.E0;

// The below section is used to calculate the effective emittance

    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        if(eSurive[i]!=0) continue;

        x2Aver  +=  pow(ePositionX[i],2);
        y2Aver  +=  pow(ePositionY[i],2);
        z2Aver  +=  pow(ePositionZ[i],2);
        px2Aver +=  pow(eMomentumX[i],2);
        py2Aver +=  pow(eMomentumY[i],2);
        pz2Aver +=  pow(eMomentumZ[i],2);

        xpxAver +=  (ePositionX[i]) * (eMomentumX[i]);
        ypyAver +=  (ePositionY[i]) * (eMomentumY[i]);
        zpzAver +=  (ePositionZ[i]) * (eMomentumZ[i]);
    }



    x2Aver  /= macroEleNumSurivePerBunch;
    y2Aver  /= macroEleNumSurivePerBunch;
    z2Aver  /= macroEleNumSurivePerBunch;
    px2Aver /= macroEleNumSurivePerBunch;
    py2Aver /= macroEleNumSurivePerBunch;
    pz2Aver /= macroEleNumSurivePerBunch;
    xpxAver /= macroEleNumSurivePerBunch;
    ypyAver /= macroEleNumSurivePerBunch;
    zpzAver /= macroEleNumSurivePerBunch;

    rmsEffectiveRingEmitX = sqrt(x2Aver * px2Aver - pow(xpxAver,2));
    rmsEffectiveRingEmitY = sqrt(y2Aver * py2Aver - pow(ypyAver,2));
    rmsEffectiveRingEmitZ = sqrt(z2Aver * pz2Aver - pow(zpzAver,2));


    rmsEffectiveRx = sqrt(rmsEffectiveRingEmitX * latticeInterActionPoint.twissBetaX[k]);
    rmsEffectiveRy = sqrt(rmsEffectiveRingEmitY * latticeInterActionPoint.twissBetaY[k]);



// The below section is used to calculate the  rms emittance
// the obtained rms size is used to calculate the interaction between beam and ion

    x2Aver=0.E0;
    y2Aver=0.E0;
    z2Aver=0.E0;
    px2Aver=0.E0;
    py2Aver=0.E0;
    pz2Aver=0.E0;
    xpxAver=0.E0;
    ypyAver=0.E0;
    zpzAver=0.E0;


    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        if(eSurive[i]!=0) continue;

        x2Aver  +=  pow(ePositionX[i]-xAver ,2);
        y2Aver  +=  pow(ePositionY[i]-yAver ,2);
        z2Aver  +=  pow(ePositionZ[i]-zAver ,2);

        px2Aver +=  pow(eMomentumX[i]-pxAver,2);
        py2Aver +=  pow(eMomentumY[i]-pyAver,2);
        pz2Aver +=  pow(eMomentumZ[i]-pzAver,2);

        xpxAver +=  (ePositionX[i]-xAver) * (eMomentumX[i]-pxAver);
        ypyAver +=  (ePositionY[i]-yAver) * (eMomentumY[i]-pyAver);
        zpzAver +=  (ePositionZ[i]-zAver) * (eMomentumZ[i]-pzAver);
    }



    x2Aver  /= macroEleNumSurivePerBunch;
    y2Aver  /= macroEleNumSurivePerBunch;
    z2Aver  /= macroEleNumSurivePerBunch;
    px2Aver /= macroEleNumSurivePerBunch;
    py2Aver /= macroEleNumSurivePerBunch;
    pz2Aver /= macroEleNumSurivePerBunch;
    xpxAver /= macroEleNumSurivePerBunch;
    ypyAver /= macroEleNumSurivePerBunch;
    zpzAver /= macroEleNumSurivePerBunch;

    emittanceX = sqrt(x2Aver * px2Aver - pow(xpxAver,2));
    emittanceY = sqrt(y2Aver * py2Aver - pow(ypyAver,2));
    emittanceZ = sqrt(z2Aver * pz2Aver - pow(zpzAver,2));

    rmsRx = sqrt(emittanceX * latticeInterActionPoint.twissBetaX[k]);
    rmsRy = sqrt(emittanceY * latticeInterActionPoint.twissBetaY[k]);

    rmsBunchLength  = sqrt(z2Aver) ;
    rmsEnergySpread = sqrt(pz2Aver);

    GetZMinMax();

}

void MPBunch::SSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    // (1) get the force of accumulated ions due to the bunch electron beam --- BassettiErskine model.

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

    for(int p=0;p<latticeInterActionPoint.gasSpec;p++)
    {
        nI0 = latticeInterActionPoint.macroIonCharge[k][p];
        ionMassNumber = latticeInterActionPoint.ionMassNumber[p];

        coeffI = 2.0 * nE0*ElecClassicRadius * ElectronMassEV/IonMassEV/ionMassNumber * CLight;   // [m * m/s]
        coeffE = 2.0 * ElecClassicRadius/rGamma;                                                  // [m]

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

            eFxTemp = eFxTemp + tempFx;
            eFyTemp = eFyTemp + tempFy;
        }
   
    } 



    // (2) get the force of certain bunched beam due to accumulated ions --- BassettiErskine model.
    // The process to get the force from accumulated ion beam is the similar to inverse "strong-weak" model.
    
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        if(eSurive[i] != 0) continue;

        eFxDueToIon[i] =0.E0;
        eFyDueToIon[i] =0.E0;

        for(int p=0;p<latticeInterActionPoint.gasSpec;p++)
        {
            nI0 = latticeInterActionPoint.macroIonCharge[k][p] * latticeInterActionPoint.ionAccumuNumber[k][p];
            coeffE = 2.0*nI0*ElecClassicRadius/rGamma;
            
            rmsRxTemp = latticeInterActionPoint.ionAccumuRMSX[k][p];
            rmsRyTemp = latticeInterActionPoint.ionAccumuRMSY[k][p];

            posx    =  ePositionX[i] - latticeInterActionPoint.ionAccumuAverX[k][p];
            posy    =  ePositionY[i] - latticeInterActionPoint.ionAccumuAverY[k][p];

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

            eFxDueToIon[i] +=    coeffE * tempFx;
            eFyDueToIon[i] +=    coeffE * tempFy;

        }          
    }    

    latticeInterActionPoint.GetTotIonCharge();
    totIonCharge = latticeInterActionPoint.totIonCharge;
}



void MPBunch::BunchTransferDueToLatticeLRigid(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    // Ref. bunch.h that ePositionZ = - ePositionT * c. head pariticles: deltaT<0, ePositionZ[i]>0.

    int resNum        = inputParameter.ringParRf->resNum;
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    double fRF         = f0 * ringHarmH;
    int bunchBinNumberZ = inputParameter.ringParRf->rfBunchBinNum;
     double u0         = inputParameter.ringParBasic->u0;
    GetZMinMax();
    // cut the bunch into  bins and store the particle index
    double dzBin = (zMaxCurrentTurn - zMinCurrentTurn) / bunchBinNumberZ;
    double dtBin = dzBin / CLight / rBeta;

    // bunch bins are alined from head to tail [zmax-> zMin] when i=0,1,2,...
    for(int i=0;i<posZBins.size();i++)
    {
        posZBins[i] = zMaxCurrentTurn - i * dzBin  - dzBin / 2.0;   
    }

    vector<vector<int>> histoParIndex;
    histoParIndex.resize(bunchBinNumberZ);
    for(int i=0;i<ePositionZ.size();i++)
    {
        if(eSurive[i]!=0) continue;
        int index = floor((  zMaxCurrentTurn - ePositionZ[i] ) / dzBin); 
        histoParIndex[index].push_back(i);    
    }

    for(int j=0;j<resNum;j++)
    {
        double resHarm = cavityResonator.resonatorVec[j].resHarm;
        double resFre  = cavityResonator.resonatorVec[j].resFre;
        double tB, deltaL, cPsi;

        // accumulate is along bin-by-bin--------
 
        complex<double> selfLossVolAccume=(0.0,0.0);
        // do the averaging process in  Cartesian system...
        // polar system average will bring trouble when all values around -pi        
        double cavVoltageReal = 0.E0;
        double cavVoltageImag = 0.E0;
        double induceVolReal = 0.E0;
        double induceVolImag = 0.E0;
        double genVolReal = 0.E0;
        double genVolImag = 0.E0;

        complex<double> vb0=(0,0);
        complex<double> cavVoltage=(0.E0,0.E0);
        complex<double> genVoltage=(0.E0,0.E0);
        int particleInBunch=0;

        // calculate the Vb0 once per bunch
        vb0  = complex<double>(-1 * 2 * PI * resFre * cavityResonator.resonatorVec[j].resShuntImpRs /
                            cavityResonator.resonatorVec[j].resQualityQ0,0.E0) * electronNumPerBunch * ElectronCharge;  // [Volt]s

        // loop for bin-by-bin in one bunch from head to tail...
        for(int k=0;k<histoParIndex.size();k++) 
        {
            genVoltage = cavityResonator.resonatorVec[j].resGenVol * exp( - li * posZBins[k]  / CLight * 2. * PI * double(resHarm) * fRF );   
            cavVoltage = cavityResonator.resonatorVec[j].vbAccum   +  genVoltage;
                        
            // loop for particle in each bin --- histoParIndex is ligned from head to tail [-t,t]          
            for(int i=0;i<histoParIndex[k].size();i++)
            {
                int index          = histoParIndex[k][i];
                eMomentumZ[index] += cavVoltage.real() / electronBeamEnergy / pow(rBeta,2);
                eMomentumZ[index] += vb0.real()/2.0    / electronBeamEnergy / pow(rBeta,2);
                eMomentumZ[index] -= u0                / electronBeamEnergy / pow(rBeta,2); 
            }

            // to get the weighing average of cavvity voltage, beam inuced voltage...
            particleInBunch   += histoParIndex[k].size();
            selfLossVolAccume += vb0/2.0    * double(histoParIndex[k].size());

            cavVoltageReal    += cavVoltage.real() * double(histoParIndex[k].size());
            cavVoltageImag    += cavVoltage.imag() * double(histoParIndex[k].size());
            genVolReal        += genVoltage.real() * double(histoParIndex[k].size());
            genVolImag        += genVoltage.imag() * double(histoParIndex[k].size());
            induceVolReal     += (cavityResonator.resonatorVec[j].vbAccum + vb0).real() * double(histoParIndex[k].size()); 
            induceVolImag     += (cavityResonator.resonatorVec[j].vbAccum + vb0).imag() * double(histoParIndex[k].size());         
        }
        
        cavVoltageReal /= double(particleInBunch);
        cavVoltageImag /= double(particleInBunch);
        genVolReal     /= double(particleInBunch);
        genVolImag     /= double(particleInBunch);
        induceVolReal  /= double(particleInBunch);
        induceVolImag  /= double(particleInBunch); 

        bunchRFModeInfo->cavVolBunchCen[j]    = complex<double>(cavVoltageReal,cavVoltageImag);
        bunchRFModeInfo->genVolBunchAver[j]   = complex<double>(genVolReal    ,genVolImag    ); 
        bunchRFModeInfo->induceVolBunchCen[j] = complex<double>(induceVolReal ,induceVolImag); 
        bunchRFModeInfo->selfLossVolBunchCen[j]  =  selfLossVolAccume / double(particleInBunch);
   

        // beam induced voltage rotate to next bunch   
        tB     = timeFromCurrnetBunchToNextBunch;
        deltaL = tB / cavityResonator.resonatorVec[j].tF;        
        cPsi   = 2.0 * PI * cavityResonator.resonatorVec[j].resFre * tB;
        cavityResonator.resonatorVec[j].vbAccum +=  vb0;   
        cavityResonator.resonatorVec[j].vbAccum *=  exp(- deltaL ) * exp (li * cPsi);
    }

}


void MPBunch:: BunchTransferDueToLatticeLBinByBin(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    // Ref. bunch.h that ePositionZ = - ePositionT * c. head pariticles: deltaT<0, ePositionZ[i]>0.
    double *synchRadDampTime = inputParameter.ringParBasic->synchRadDampTime;
    
    int resNum        = inputParameter.ringParRf->resNum;
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    int bunchBinNumberZ = inputParameter.ringParRf->rfBunchBinNum;
    double u0         = inputParameter.ringParBasic->u0;
    double fRF         = f0 * ringHarmH;

    // cut the bunch into  bins and store the particle index
    double dzBin = (zMaxCurrentTurn - zMinCurrentTurn) / bunchBinNumberZ;
    double dtBin = dzBin / CLight / rBeta;


    // bunch bins are alined from head to tail [zmax-> zMin] when i=0,1,2,...
    for(int i=0;i<posZBins.size();i++)
    {
        posZBins[i] = zMaxCurrentTurn - i * dzBin  - dzBin / 2.0;   
    }

    vector<vector<int>> histoParIndex;
    histoParIndex.resize(bunchBinNumberZ);
    for(int i=0;i<ePositionZ.size();i++)
    {
        if(eSurive[i]!=0) continue;
        int index = floor((  zMaxCurrentTurn - ePositionZ[i] ) / dzBin);         
        histoParIndex[index].push_back(i);    
    }


    for(int j=0;j<resNum;j++)
    {           
        double resHarm = cavityResonator.resonatorVec[j].resHarm;
        double resFre  = cavityResonator.resonatorVec[j].resFre;
        double tB, deltaL, cPsi;
        tB = deltaL = cPsi = 0.E0;    
        // accumulate is along bin-by-bin--------
        complex<double> selfLossVolAccume=(0.E0,0.E0);

        // do the averaging process in  Cartesian system...
        // polar system average will bring trouble when all values around -pi        
        double cavVoltageReal = 0.E0;
        double cavVoltageImag = 0.E0;
        double induceVolReal = 0.E0;
        double induceVolImag = 0.E0;
        double genVolReal = 0.E0;
        double genVolImag = 0.E0;

        complex<double> vb0=(0,0);
        complex<double> cavVoltage=(0.E0,0.E0);
        complex<double> genVoltage=(0.E0,0.E0);
        int particleInBunch=0;


        // loop for bin-by-bin in one bunch, bins is alined from head to tail
        // each bin excite beam induced voltage itself
        for(int k=0;k<histoParIndex.size();k++) 
        {
            // change in the real time frame for generator voltage calculation...
            
            genVoltage = cavityResonator.resonatorVec[j].resGenVol * exp(  - li * posZBins[k]  / CLight * 2. * PI * double(resHarm) * fRF); 
            cavVoltage = cavityResonator.resonatorVec[j].vbAccum   +  genVoltage;
            
            vb0  = complex<double>(-1 * 2 * PI * resFre * cavityResonator.resonatorVec[j].resShuntImpRs /
                    cavityResonator.resonatorVec[j].resQualityQ0,0.E0) * macroEleCharge * double(histoParIndex[k].size()) * ElectronCharge;  // [Volt]

            // loop for particle in each bin          
            for(int i=0;i<histoParIndex[k].size();i++)
            {
                int index          = histoParIndex[k][i];
                eMomentumZ[index] += cavVoltage.real() / electronBeamEnergy / pow(rBeta,2);
                eMomentumZ[index] += vb0.real()/2.0    / electronBeamEnergy / pow(rBeta,2);
                eMomentumZ[index] -= u0                / electronBeamEnergy / pow(rBeta,2); 
            }
  
            // to get the weighing average of cavvity voltage, beam inuced voltage...
            particleInBunch   += histoParIndex[k].size();
            selfLossVolAccume += vb0/2.0 * double(histoParIndex[k].size());
           
            cavVoltageReal    += cavVoltage.real() * double(histoParIndex[k].size());
            cavVoltageImag    += cavVoltage.imag() * double(histoParIndex[k].size());
            genVolReal        += genVoltage.real() * double(histoParIndex[k].size());
            genVolImag        += genVoltage.imag() * double(histoParIndex[k].size());
            induceVolReal     += (cavityResonator.resonatorVec[j].vbAccum + vb0).real() * double(histoParIndex[k].size()); 
            induceVolImag     += (cavityResonator.resonatorVec[j].vbAccum + vb0).imag() * double(histoParIndex[k].size());

            // beam induced voltage retotate and decay bin-by-bin            
            tB     = dtBin;
            deltaL = tB  / cavityResonator.resonatorVec[j].tF;
            cPsi   = 2.0 * PI *  cavityResonator.resonatorVec[j].resFre * tB;
            cavityResonator.resonatorVec[j].vbAccum += vb0;
            cavityResonator.resonatorVec[j].vbAccum = cavityResonator.resonatorVec[j].vbAccum * exp(- deltaL ) * exp (li * cPsi);
    
        }
        
        cavVoltageReal /= double(particleInBunch);
        cavVoltageImag /= double(particleInBunch);
        genVolReal     /= double(particleInBunch);
        genVolImag     /= double(particleInBunch);
        induceVolReal  /= double(particleInBunch);
        induceVolImag  /= double(particleInBunch); 

        bunchRFModeInfo->cavVolBunchCen[j]    = complex<double>(cavVoltageReal,cavVoltageImag);
        bunchRFModeInfo->genVolBunchAver[j]   = complex<double>(genVolReal    ,genVolImag    ); 
        bunchRFModeInfo->induceVolBunchCen[j] = complex<double>(induceVolReal ,induceVolImag); 
        bunchRFModeInfo->selfLossVolBunchCen[j]  =  selfLossVolAccume / double(particleInBunch) ;

        
        // beam induced voltage rotate to next bunch   
        tB     = timeFromCurrnetBunchToNextBunch;
        deltaL = tB / cavityResonator.resonatorVec[j].tF;        
        cPsi   = 2.0 * PI * cavityResonator.resonatorVec[j].resFre * tB;
        cavityResonator.resonatorVec[j].vbAccum +=  vb0;   
        cavityResonator.resonatorVec[j].vbAccum *=  exp(- deltaL ) * exp (li * cPsi);
    }

    
}

void MPBunch::BunchLongiMomentumUpdateDuetoRF(const ReadInputSettings &inputParameter)
{
    double *synchRadDampTime = inputParameter.ringParBasic->synchRadDampTime;
    double t0         = inputParameter.ringParBasic->t0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double eta        = inputParameter.ringParBasic->eta;
    double u0         = inputParameter.ringParBasic->u0;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;

    // momentum and position update
    double rmsEnergySpreadTemp = inputParameter.ringParBasic-> sdelta0; // quantum excitation  natural beam energy spread.
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> qEpsilon{0,2*rmsEnergySpreadTemp/sqrt(synchRadDampTime[2])};

    for (int i=0;i<eMomentumZ.size();i++)
    {
        eMomentumZ[i]  -= u0 / electronBeamEnergy / pow(rBeta,2); 
        if(inputParameter.ringRun->synRadDampingFlag[1]==1)
        {
            eMomentumZ[i] +=  qEpsilon(gen);
            eMomentumZ[i] *= (1.0-2.0/synchRadDampTime[2]);        
        }            
        // ePositionZ[i] -= eta * t0 * CLight  * eMomentumZ[i];
    }   
}



void MPBunch::BunchTransferDueToLatticeLNoInstability(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    // Ref. bunch.h that ePositionZ = - ePositionT * c. head pariticles: deltaT<0, ePositionZ[i]>0.
    // During the tracking, from head to tail means ePositionZMin from [+,-];   

    double *synchRadDampTime = inputParameter.ringParBasic->synchRadDampTime;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    double rmsEnergySpreadTemp = inputParameter.ringBunchPara-> rmsEnergySpread; // quantum excitation  natural beam energy spread.
    //rmsEnergySpreadTemp = rmsEnergySpread;                  


    std::normal_distribution<> qEpsilon{0,2*rmsEnergySpreadTemp/sqrt(synchRadDampTime[2])};

    int resNum        = inputParameter.ringParRf->resNum;
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double t0         = inputParameter.ringParBasic->t0;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double eta        = inputParameter.ringParBasic->eta;
    double u0         = inputParameter.ringParBasic->u0;
    double *alphac    = inputParameter.ringParBasic->alphac;
    double sdelta0    = inputParameter.ringParBasic->sdelta0;                   // natural beam energy spread
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
    
    int bunchBinNumberZ = inputParameter.ringParRf->rfBunchBinNum;
  
    vector< vector<double> > hamiltonPotenWellTemp(resNum);    
    for(int j=0;j<resNum;j++)
    {
        hamiltonPotenWellTemp[j].resize(bunchBinNumberZ);

       for(int i=0;i<bunchBinNumberZ;i++)
        {
            hamiltonPotenWellTemp[j][i]=0.E0;
        }
    }

    
    double zMin,zMax;   
    GetZMinMax();
    zMin = zMinCurrentTurn;
    zMax = zMaxCurrentTurn;

    double poszMin = zMin ;
    double poszMax = zMax ;
    //poszMin = -1.0/f0/ringHarmH * rBeta * CLight/10;  // set the \pm (-2 pi/10, 2 pi/10) simulation range.
    //poszMax =  1.0/f0/ringHarmH * rBeta * CLight/10;
    double dpzTemp;

    double dzBin = (poszMax - poszMin) / bunchBinNumberZ;
    double dtBin = dzBin / CLight;
    vector<vector<int>> histoParIndex;
    histoParIndex.resize(bunchBinNumberZ);

    int counter;

    for(int i=0;i<ePositionZ.size();i++)
    {
        int index = int( (ePositionZ[i] -poszMin  ) / dzBin );  // head particle index stores in histoParIndex[i]. Smaller i represents head particle.  
        histoParIndex[index].push_back(i);
    }


    for(int j=0;j<resNum;j++)
    {
        resHarm = cavityResonator.resonatorVec[j].resHarm;
        resFre  = cavityResonator.resonatorVec[j].resFre;
        resPhaseReq = cavityResonator.resonatorVec[j].resPhaseReq;

        complex<double> cavVoltageAccume=(0.0,0.0); // along one bunch, which is cutted into different bins
        complex<double> selfLossVolAccume=(0.0,0.0);
        complex<double> induceVolAccume=(0.0,0.0);
        int particleInBunch=0;

        counter=0;

        for(int k=0;k<bunchBinNumberZ;k++)  //bin by bin... From small bin index to large -> head to tail 
        {
            vb0  = complex<double>(-1 * cavityResonator.resonatorVec[j].resFre * 2 * PI * cavityResonator.resonatorVec[j].resShuntImpRs /
                                        cavityResonator.resonatorVec[j].resQualityQ0,0.E0) * macroEleCharge * double(histoParIndex[k].size()) * ElectronCharge;  // [Volt]

            cavVoltage = cavityResonator.resonatorVec[j].vbAccum + bunchRFModeInfo->genVolBunchAver[j]; 
            cavVolAmp  = abs(cavVoltage);
            cavVolArg  = arg(cavVoltage);

            cavVolInfoVsLongBins[j][k] = cavVoltage;
            posZBins[k]                = poszMin + k * dzBin;

            // get bunch intensity in current bin
            densProfVsBin[k] = histoParIndex[k].size()/dzBin;        //C/m  

            if(histoParIndex[k].size()==0)
            {
                // no particle in this bin;
                dpzTemp = cavVolAmp /electronBeamEnergy * (cos( - 2. * PI * ringHarmH * f0 * resHarm * posZBins[k] / ( rBeta * CLight)  +  cavVolArg ) - cos(resPhaseReq) );
                dpzTemp = dpzTemp * (1.0-2.0/synchRadDampTime[2]);
            }
            else
            {
                //particle exist in this bin;
                for(int i=0;i<histoParIndex[k].size();i++)
                {
                    int index          = histoParIndex[k][i];
                    dpzTemp            = cavVolAmp / electronBeamEnergy  * cos( - 2. * PI * ringHarmH * f0 * resHarm * ePositionZ[index] / ( rBeta * CLight)  +  cavVolArg );
                    eMomentumZ[index] += dpzTemp / pow(rBeta,2);
                    eMomentumZ[index] += vb0.real()/2.0/electronBeamEnergy / pow(rBeta,2);
                    eMomentumZ[index] -= u0            /electronBeamEnergy / pow(rBeta,2); 
                }
                counter++ ;

                particleInBunch   +=  histoParIndex[k].size();
                cavVoltageAccume  +=  cavVoltage * double(histoParIndex[k].size());
                selfLossVolAccume +=  vb0/2.0    * double(histoParIndex[k].size());
                induceVolAccume   += cavityResonator.resonatorVec[j].vbAccum *  double(histoParIndex[k].size());

                dpzTemp = cavVolAmp /electronBeamEnergy * (cos(  - 2. * PI * ringHarmH * f0 * resHarm * posZBins[k]  / ( rBeta * CLight)  +  cavVolArg ) - cos(resPhaseReq))
                                        + vb0.real()/2.0/electronBeamEnergy;

                dpzTemp = dpzTemp * (1.0-2.0/synchRadDampTime[2]);

                // beam induced voltage accumu  -- due to each bin
                cavityResonator.resonatorVec[j].vbAccum += vb0;

                //phasor rotate and decay due to each bin with time distance dtBin
                
                deltaL = dtBin  / cavityResonator.resonatorVec[j].tF;
                cPsi   = deltaL * tan(cavityResonator.resonatorVec[j].resDeTunePsi);

                cavityResonator.resonatorVec[j].vbAccum = cavityResonator.resonatorVec[j].vbAccum * exp(- deltaL ) * exp (li * cPsi);
            }

            cavForceInfoVsLongBins[j][k] = dpzTemp;   // rad
        }
        

        bunchRFModeInfo->selfLossVolBunchCen[j] =  selfLossVolAccume /double(particleInBunch) ;
        bunchRFModeInfo->induceVolBunchCen[j]   =  induceVolAccume   /double(particleInBunch) ;
        bunchRFModeInfo->cavVolBunchCen[j]      =  cavVoltageAccume  /double(particleInBunch) ;
                
        hamiltonPotenWellTemp[j][0]=0;
        for(int k=1;k<bunchBinNumberZ;k++)
        {
            hamiltonPotenWellTemp[j][k]  = hamiltonPotenWellTemp[j][k-1] + (cavForceInfoVsLongBins[j][k-1] + cavForceInfoVsLongBins[j][k])/2.0; // dimensionless (22) PRAB 2018,21 012001
        }                                                                                                                                       // dz is dealed with later

        tB     = timeFromCurrnetBunchToNextBunch;
        //tB     = bunchGap * t0 / ringHarmH;

        deltaL = tB / cavityResonator.resonatorVec[j].tF;
        cPsi   = deltaL * tan(cavityResonator.resonatorVec[j].resDeTunePsi);
        cavityResonator.resonatorVec[j].vbAccum = cavityResonator.resonatorVec[j].vbAccum * exp(- deltaL ) * exp (li * cPsi);       
    }

    for(int k=0;k<bunchBinNumberZ;k++)
    {
        hamiltonPotenWell[k]=0;
        for(int j=0;j<resNum;j++)
        {
            hamiltonPotenWell[k] +=  hamiltonPotenWellTemp[j][k];
        }
    }
   

    auto iter = max_element(hamiltonPotenWell.begin(),hamiltonPotenWell.end() );
    double temp = *iter;

    // hamiltonian refer.to Eq. 22, PRAB 2018-21-012001// get the normalized profile
    double densityNormAna=0.0;
    double densityNormTracking=0.0;
    for(int k=0;k<bunchBinNumberZ;k++)
    {
        hamiltonPotenWell[k] -= temp;
        hamiltonPotenWell[k] = hamiltonPotenWell[k] * alphac[0] / ( 2 * PI * ringHarmH ) / pow(alphac[0] * sdelta0,2) * 2 * PI * dzBin * ringHarmH * f0 / rBeta /CLight;
        densProfVsBinAnalytical[k] = exp(hamiltonPotenWell[k]);
        densityNormAna += densProfVsBinAnalytical[k];
        densityNormTracking += densProfVsBin[k];
    }
    for(int k=0;k<bunchBinNumberZ;k++)
    {
        densProfVsBinAnalytical[k] = densProfVsBinAnalytical[k] / densityNormAna;
        densProfVsBin[k]           = densProfVsBin[k]           / densityNormTracking;
    }

    //get the analytical bunchcenter and bunchLength;
    double tempNorm=0.E0;
    zAverAnalytical =0.E0;
    bunchLengthAnalytical =0.E0;
    for(int k=0;k<bunchBinNumberZ;k++)
    {
        zAverAnalytical       += densProfVsBinAnalytical[k] * posZBins[k];
    }
    for(int k=0;k<bunchBinNumberZ;k++)
    {
        bunchLengthAnalytical += densProfVsBinAnalytical[k] * pow(posZBins[k]-zAverAnalytical,2);
    }
    
    bunchLengthAnalytical = sqrt(bunchLengthAnalytical);


}

void MPBunch::GetZMinMax()
{

    double zMin = ePositionZ[0];
    double zMax = ePositionZ[0];

    for (int i=0;i<ePositionZ.size();i++)
    {        
        if(eSurive[i]!=0) continue;      //  eSurive[i]=0 when partilce is not lost 
    
        if(zMin>ePositionZ[i])
        {
           zMin = ePositionZ[i]; 
        }
        if(zMax<ePositionZ[i])
        {
           zMax = ePositionZ[i]; 
        } 
    }

    zMinCurrentTurn = zMin - 1.e-6;
    zMaxCurrentTurn = zMax + 1.e-6; 
}

void MPBunch::BunchTransferDueToSRWake(const  ReadInputSettings &inputParameter, WakeFunction &sRWakeFunction, const LatticeInterActionPoint &latticeInterActionPoint,int turns)
{
    // Ref. bunch.h that ePositionZ = - ePositionT * c. head pariticles: deltaT<0, ePositionZ[i]>0.
    // During the tracking, from head to tail means ePositionZMin from [+,-];   
     
    int bunchBinNumberZ = inputParameter.ringSRWake->SRWBunchBinNum;

    for(int i=0;i<3;i++)
    {
        srWakePoten[i].resize(bunchBinNumberZ,0.e0);
    }
    double zMin,zMax; 

    GetZMinMax();
    zMin = zMinCurrentTurn;
    zMax = zMaxCurrentTurn;

    double poszMin = zMin  - 2 * rmsBunchLength;
    double poszMax = zMax  + 2 * rmsBunchLength;

    double dzBin = (poszMax - poszMin) / bunchBinNumberZ;
    double dtBin = dzBin / CLight;
    vector<vector<int>> histoParIndex;
    histoParIndex.resize(bunchBinNumberZ);
    
    vector<double> averXAlongBunch(bunchBinNumberZ,0);
    vector<double> averYAlongBunch(bunchBinNumberZ,0);

    for(int i=0;i<ePositionZ.size();i++)
    {
        int index = int( (ePositionZ[i] -poszMin ) / dzBin );          // head particle index stores in histoParIndex[i]. Smaller i represents head particle.  
        histoParIndex[index].push_back(i);
    }
     
    int partID;
    
    for(int i=0;i<bunchBinNumberZ;i++)
    {
        for(int j=0;j<histoParIndex[i].size();j++)
        {
            partID = histoParIndex[i][j];
            if(eSurive[partID]!=0) continue;
           
            averXAlongBunch[i] += ePositionX[partID]; 
            averYAlongBunch[i] += ePositionY[partID]; 
        }

        if (histoParIndex[i].size()==0)
        {
            averXAlongBunch[i]=0;
            averYAlongBunch[i]=0;    
        }
        else
        {
            averXAlongBunch[i] /= histoParIndex[i].size();
            averYAlongBunch[i] /= histoParIndex[i].size();
        }
    }
    
    vector<double> wakeFunji;
    double tauji;   
    int partNumInBin;

    ofstream fout(inputParameter.ringSRWake->SRWWakePotenWriteTo+".sdds",ios_base::app);
    if(turns==0) 
    {
        fout<<"SDDS1"<<endl;
        fout<<"&column name=z,              units=m,              type=float,  &end" <<endl;
        fout<<"&column name=profile,                              type=float,  &end" <<endl;
        fout<<"&column name=wakePotenX,                        type=float,  &end" <<endl;
        fout<<"&column name=wakePotenY,                        type=float,  &end" <<endl;
        fout<<"&column name=wakePotenZ,                        type=float,  &end" <<endl;
        fout<<"&data mode=ascii, &end"                                               <<endl;
    }

    if(turns% (inputParameter.ringRun->bunchInfoPrintInterval)==0)
    {
        fout<<"! page number " << int(turns/100)+1 <<endl;
        fout<<bunchBinNumberZ<<endl;
    }
    for(int i=0;i<bunchBinNumberZ;i++)
    {
        // RW puesdo wake function -- integration range have to modified
        if(!inputParameter.ringSRWake->pipeGeoInput.empty())
        {
             for(int j=0;j<bunchBinNumberZ;j++)                         //integration range have to modified
             {
                tauji = (i-j) * dtBin;                                        // [s]  must .LE. 0             
                partNumInBin = histoParIndex[j].size();
                wakeFunji = sRWakeFunction.GetRWSRWakeFun(tauji);                                                             
                srWakePoten[0][i] += wakeFunji[0] * partNumInBin * averXAlongBunch[j];  // [V/C m] [m] ->[V/C] X
                srWakePoten[1][i] += wakeFunji[1] * partNumInBin * averXAlongBunch[j];  // [V/C m] [m] ->[V/C] Y
                srWakePoten[2][i] += wakeFunji[2] * partNumInBin;        
             }            
        }

        if(!inputParameter.ringSRWake->bbrInput.empty())
        {
            for(int j=0;j<bunchBinNumberZ;j++)
             {
                tauji = (i-j) * dtBin;                                        // [s]  must .LE. 0             
                partNumInBin = histoParIndex[j].size();
                wakeFunji = sRWakeFunction.GetBBRWakeFun1(tauji);                                                             
                srWakePoten[0][i] += wakeFunji[0] * partNumInBin * averXAlongBunch[j];  // [V/C m] [m] ->[V/C] X
                srWakePoten[1][i] += wakeFunji[1] * partNumInBin * averXAlongBunch[j];  // [V/C m] [m] ->[V/C] Y
                srWakePoten[2][i] += wakeFunji[2] * partNumInBin;                        
             }         
        }

        srWakePoten[0][i] *= (-1)  * ElectronCharge * macroEleCharge / electronEnergy;              // [V/V] [rad]  Eq. (3.7) -- multiplty -1; 
        srWakePoten[1][i] *= (-1)  * ElectronCharge * macroEleCharge / electronEnergy;              // [V/V] [rad]  Eq. (3.7) -- multiplty -1;  
        srWakePoten[2][i] *= (-1)  * ElectronCharge * macroEleCharge / electronEnergy;              // [V/V] [rad]  Eq. (3.7) -- multiplty -1;   

        if(turns% (inputParameter.ringRun->bunchInfoPrintInterval)==0) 
        {
            fout<<setw(15)<<left<< i*dtBin*CLight + poszMin
                <<setw(15)<<left<<histoParIndex[i].size()
                <<setw(15)<<left<<srWakePoten[0][i]
                <<setw(15)<<left<<srWakePoten[1][i]
                <<setw(15)<<left<<srWakePoten[2][i]
                <<endl;
        }
    }
    fout.close();
    

    // can be updated to include the quadrupole wakes. -- left for future. 
    
    for(int i=0;i<bunchBinNumberZ;i++)
    {
        for(int j=0;j<histoParIndex[i].size();j++)
        {
            partID = histoParIndex[i][j];
            if(eSurive[partID]!=0) continue;
            //eMomentumX[partID] += srWakePoten[0][i];            //rad
            //eMomentumY[partID] += srWakePoten[1][i];            //rad    
            eMomentumZ[partID] += srWakePoten[2][i];            //rad
        }
    }
}


void MPBunch::BBImpBunchInteraction(const ReadInputSettings &inputParameter, const BoardBandImp &boardBandImp )
{
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    double zMinBin = -boardBandImp.zMax;
    double dzBin   =  boardBandImp.dz;
    
    GetBunchProfileForBeamBBImpEffect(inputParameter,boardBandImp);
    int nBins = profileForBunchBBImp.size();

    gsl_fft_real_workspace          *workspace  = gsl_fft_real_workspace_alloc (nBins);
    gsl_fft_real_wavetable          *wavetable0 = gsl_fft_real_wavetable_alloc (nBins); 
    gsl_fft_halfcomplex_wavetable   *wavetable1 = gsl_fft_halfcomplex_wavetable_alloc (nBins);
   
	double temp[nBins];                     
    double bunchSpectrum[nBins];        // bunch spectrum only in the positive frequency side
	for(int i=0;i<nBins;i++)
	{
		temp[i] = profileForBunchBBImp[i];              //  [C/s]
    }
    gsl_fft_real_transform (temp,1, nBins, wavetable0, workspace);  //  rho(\tau) -> rho(f) [C/s] -> [C/s]  // refer to (2.69)   
    for(int i=0;i<nBins;i++) bunchSpectrum[i] = temp[i];

    // (0) z longitudial kick s	
    temp[0] = 0;  // in longitudial, DC term impedance is zZ(0) = (0 + 0I)    
    for(int i=1;i<boardBandImp.zZImp.size();i++)
    {
        complex<double> indVF = complex<double> (temp[2*i-1], temp[2*i] ) * boardBandImp.zZImp[i];
        temp[2*i-1] = indVF.real();
        temp[2*i  ] = indVF.imag();                                       // [C/s] * [Ohm] = [V]                          
    }
    gsl_fft_halfcomplex_inverse (temp, 1, nBins, wavetable1, workspace);  // [V] -> [V]
    for(int i=0;i<nBins;i++) beamIndVFromBBImpZ[i] = - temp[i];           // Eq. (3.7)  

    //-------------------------------------------------------------------------------    
    // (1) x kick 
    // for(int i=0;i<nBins;i++) temp[i] = bunchSpectrum[i]; 
    // temp[0] = - bunchSpectrum[0] * boardBandImp.zDxImp[0].imag();     // in x DC term impedance is zX(0) = (0 + I * Im);  
    // for(int i=1;i<boardBandImp.zDxImp.size();i++)
    // {
    //     complex<double> indVF = complex<double> (temp[2*i-1], temp[2*i] ) * boardBandImp.zDxImp[i] * li;  
    //     temp[2*i-1] = indVF.real();
    //     temp[2*i  ] = indVF.imag();                                       // [C/s] * [Ohm/m] = [V/m]                          
    // }
    // gsl_fft_halfcomplex_inverse (temp, 1, nBins, wavetable1, workspace);  // [V/m] -> [V/m]
    // for(int i=0;i<nBins;i++) beamIndVFromBBImpX[i] = temp[i];             // Eq. (3.51)   
    //---------------------------------------------------------------------------------------------------

    // (2) y kick 
    // for(int i=0;i<nBins;i++) temp[i] = bunchSpectrum[i];
    // temp[0] = - bunchSpectrum[0] * boardBandImp.zDyImp[0].imag();     
    // for(int i=1;i<boardBandImp.zDyImp.size();i++)
    // {
    //     complex<double> indVF = complex<double> (temp[2*i-1], temp[2*i] ) * boardBandImp.zDyImp[i] * li;  
    //     temp[2*i-1] = indVF.real();
    //     temp[2*i  ] = indVF.imag();                                       // [C/s] * [Ohm/m] = [V/m]                          
    // }
    // gsl_fft_halfcomplex_inverse (temp, 1, nBins, wavetable1, workspace);  // [V/m] -> [V/m]
    // for(int i=0;i<nBins;i++) beamIndVFromBBImpY[i] = temp[i];             // Eq. (3.51)  
    //---------------------------------------------------------------------------------------------------

    // kick particles in bunch
    int index;
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        index  =  (ePositionZ[i] - zMinBin + dzBin / 2 ) / dzBin; 
        eMomentumZ[i] += beamIndVFromBBImpZ[index] / electronBeamEnergy / pow(rBeta,2);
        // eMomentumX[i] += beamIndVFromBBImpX[index] / electronBeamEnergy / pow(rBeta,2) * ePositionX[i];
        // eMomentumY[i] += beamIndVFromBBImpY[index] / electronBeamEnergy / pow(rBeta,2) * ePositionY[i];
    }


    gsl_fft_real_workspace_free (workspace);
    gsl_fft_real_wavetable_free (wavetable0);
    gsl_fft_halfcomplex_wavetable_free (wavetable1);

    //test the longitudinal wakefield here. 
    // ofstream fout("test.dat");
    // fout<<"SDDS1"<<endl;
    // fout<<"&column name=z,              units=m,              type=float,  &end" <<endl;
    // fout<<"&column name=profile,                              type=float,  &end" <<endl;
    // fout<<"&column name=indVZ,                                type=float,  &end" <<endl;
    // fout<<"&column name=indVX,                                type=float,  &end" <<endl;
    // fout<<"&column name=indVY,                                type=float,  &end" <<endl;
    // fout<<"&data mode=ascii, &end"                                               <<endl;
    // fout<<"! page number " << 1 <<endl;
    // fout<<nBins<<endl;

    // for(int i=0; i<nBins;i++)
    // {
    //     fout<<setw(15)<<left<<boardBandImp.binPosZ[i]
    //         <<setw(15)<<left<<profileForBunchBBImp[i]
    //         <<setw(15)<<left<<beamIndVFromBBImpZ[i]
    //         <<setw(15)<<left<<beamIndVFromBBImpX[i]
    //         <<setw(15)<<left<<beamIndVFromBBImpY[i]
    //         <<endl;
    // }
    // cout<<"tets"<<endl;
    // getchar();
}


void MPBunch::GetBunchProfileForBeamBBImpEffect(const ReadInputSettings &inputParameter, const BoardBandImp &boardBandImp )
{
    double zMaxBin =  boardBandImp.zMax;
    double zMinBin = -boardBandImp.zMax;
    double dzBin   =  boardBandImp.dz;
    int    nBins   =  boardBandImp.nBins;   // from 0 to max
    double rBeta   = inputParameter.ringParBasic->rBeta;
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double rfLen   = inputParameter.ringParBasic->t0 / ringHarmH * CLight;
    int index; 

    fill(profileForBunchBBImp.begin(),profileForBunchBBImp.end(),0); // array[ 2* nBins -1 ] to store the profile
 
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        index  =  (ePositionZ[i] - zMinBin + dzBin / 2 ) / dzBin; 
        profileForBunchBBImp[index] +=  macroEleCharge * ElectronCharge / dzBin *  CLight;  // [C/s]  
    }

    // for(int i=0;i<profileForBunchBBImp.size();i++)
    // {
    //     profileForBunchBBImp[i] = 1.0/sqrt(2*PI)/rmsBunchLength * exp(-pow(boardBandImp.binPosZ[i] - 0,2)/2/pow(rmsBunchLength,2)) 
    //                             * macroEleCharge * ElectronCharge * CLight * macroEleNumPerBunch; // [C/s]
    // }

    int indexStart =   ( -rfLen / 2. - zMinBin + dzBin / 2 ) / dzBin;
    int indexEnd   =   (  rfLen / 2. - zMinBin + dzBin / 2 ) / dzBin;
    int N = indexEnd - indexStart;
  
    int K=51;
    const double alpha[3] = { 10, 3.0, 10.0 };   /* alpha values */
    gsl_vector *x = gsl_vector_alloc(N);         /* input vector */
    gsl_vector *y1 = gsl_vector_alloc(N);        /* filtered output vector for alpha1 */
    gsl_vector *k1 = gsl_vector_alloc(K);        /* Gaussian kernel for alpha1 */
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_filter_gaussian_workspace *gauss_p = gsl_filter_gaussian_alloc(K);
    double sum = 0.0;

    /* generate input signal */
    for(int i=0; i<N; i++) gsl_vector_set(x, i, profileForBunchBBImp[i+indexStart]);
    
    gsl_filter_gaussian_kernel(alpha[0], 0, 0, k1);    
    gsl_filter_gaussian(GSL_FILTER_END_PADVALUE, alpha[0], 0, x, y1, gauss_p);
    
    double sum0=0.E0; 
    double sum1=0.E0; 

    for (int i=0; i<N; i++)
    {
        double xi  = gsl_vector_get(x,  i);
        double y1i = gsl_vector_get(y1, i);

        sum0 += xi;
        sum1 += y1i;  
        profileForBunchBBImp[i+indexStart] =  y1i; 
    }
    sum0 = sum0 / sum1;

    for(int i=0; i<N; i++) profileForBunchBBImp[i+indexStart] *= sum0;
    
    gsl_vector_free(x);
    gsl_vector_free(y1);
    gsl_vector_free(k1);
    gsl_rng_free(r);
    gsl_filter_gaussian_free(gauss_p);

}

