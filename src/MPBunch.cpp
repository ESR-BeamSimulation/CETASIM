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


MPBunch::MPBunch()
{
}

MPBunch::~MPBunch()
{
  
}


void MPBunch::InitialMPBunch(const  ReadInputSettings &inputParameter)
{    
    //call the Initial at base class (Bunch.h)
    Bunch::Initial(inputParameter); 
    srWakePoten.resize(3);
    int bunchBinNumberZ = inputParameter.ringParRf->rfBunchBinNum;
    
    beamCurDenZProf.resize(bunchBinNumberZ);
    posZBins.resize(bunchBinNumberZ);
    densProfVsBin.resize(bunchBinNumberZ);
    densProfVsBinAnalytical.resize(bunchBinNumberZ);
    hamiltonPotenWell.resize(bunchBinNumberZ);

    cavVolInfoVsLongBins.resize(inputParameter.ringParRf->resNum);
    cavForceInfoVsLongBins.resize(inputParameter.ringParRf->resNum);

    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        cavVolInfoVsLongBins[i].resize(bunchBinNumberZ);
        cavForceInfoVsLongBins[i].resize(bunchBinNumberZ);
    }
   
}


void MPBunch::DistriGenerator(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter, int bunchIndex)
{

// longitudial bunch phase space generatetion --simple rms in both z and z' phase space.

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

        if( temp<pow(5,2))                                      // longitudinal is truncted at 5 sigma both in z and dp direction
        {
            ePositionZ[i]=tempx;			                    // m
            eMomentumZ[i]=tempy;			                    // rad
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

    // set the first particle without any error
    ePositionX[0]=0.E0;
    ePositionY[0]=0.E0;
    ePositionZ[0]=0.E0;
    eMomentumX[0]=0.E0;
    eMomentumY[0]=0.E0;
    eMomentumZ[0]=0.E0;

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
        if(eSurive[i]==0)
        {
            continue;
        }
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
        if(eSurive[i]==0)
        {
            continue;
        }

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
        if(eSurive[i]==0)
        {
            continue;
        }

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
        if(eSurive[i] == 0) continue;

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

void MPBunch::BunchTransferDueToLatticeL(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
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
    double alphac     = inputParameter.ringParBasic->alphac;
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

    //auto iter0 = min_element(ePositionZ.begin(),ePositionZ.end() );
    //auto iter1 = max_element(ePositionZ.begin(),ePositionZ.end() );
    //double poszMin = *iter0  - 2 * rmsBunchLength;
    //double poszMax = *iter1  + 2 * rmsBunchLength;
    
    double zMin,zMax;   
    vector<double> zMinMax = GetZMinMax();
    zMin = zMinMax[0];
    zMax = zMinMax[1];


    double poszMin = zMin  - 2 * rmsBunchLength;
    double poszMax = zMax  + 2 * rmsBunchLength;
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

            cavVoltage = cavityResonator.resonatorVec[j].vbAccum + cavFBCenInfo->genVolBunchAver[j]; 
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
                    eMomentumZ[index] += dpzTemp;
                    eMomentumZ[index] += vb0.real()/2.0/electronBeamEnergy;
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
        

        cavFBCenInfo->selfLossVolBunchCen[j] =  selfLossVolAccume /double(particleInBunch) ;
        cavFBCenInfo->induceVolBunchCen[j]   =  induceVolAccume   /double(particleInBunch) ;
        cavFBCenInfo->cavVolBunchCen[j]      =  cavVoltageAccume  /double(particleInBunch) ;
        cavFBCenInfo->cavAmpBunchCen[j]      =  abs(cavFBCenInfo->cavVolBunchCen[j]);
        cavFBCenInfo->cavPhaseBunchCen[j]    =  arg(cavFBCenInfo->cavVolBunchCen[j]);

        
        //cout<<abs(cavFBCenInfo->induceVolBunchCen[j])  <<"  "<<arg(cavFBCenInfo->induceVolBunchCen[j])<<endl;
        //cout<<abs(cavFBCenInfo->cavVolBunchCen[j])     <<"  "<<arg(cavFBCenInfo->cavVolBunchCen[j])<<endl;
        //cout<<abs(cavFBCenInfo->selfLossVolBunchCen[j])<<"  "<<arg(cavFBCenInfo->selfLossVolBunchCen[j])<<endl;
        //getchar();
        

        hamiltonPotenWellTemp[j][0]=0;
        for(int k=1;k<bunchBinNumberZ;k++)
        {
            hamiltonPotenWellTemp[j][k]  = hamiltonPotenWellTemp[j][k-1] + (cavForceInfoVsLongBins[j][k-1] + cavForceInfoVsLongBins[j][k])/2.0; // dimensionless (22) PRAB 2018,21 012001
        }                                                                                                                                       // dz is dealed with later

        //tB     = timeToFirstBinOfNextBunch;
        tB     = bunchGap * t0 / ringHarmH;

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
        hamiltonPotenWell[k] = hamiltonPotenWell[k] * alphac / ( 2 * PI * ringHarmH ) / pow(alphac*sdelta0,2) * 2 * PI * dzBin * ringHarmH * f0 / rBeta /CLight;
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


    // momentum and position rotate due to next trun
    for (int i=0;i<eMomentumZ.size();i++)
    {
        eMomentumZ[i]  = eMomentumZ[i] -  u0 / electronBeamEnergy;
        eMomentumZ[i]  = eMomentumZ[i] * (1.0-2.0/synchRadDampTime[2]) + qEpsilon(gen) ;
        ePositionZ[i] -= eta * t0 * CLight * rBeta * eMomentumZ[i];
    }

}

vector<double> MPBunch::GetZMinMax()
{
    vector<double> zMinMax(2,0);
    double zMin = ePositionZ[0];
    double zMax = ePositionZ[0];

    for (int i=0;i<ePositionZ.size();i++)
    {        
        if(eSurive[i]==0) continue;        
    
        if(zMin>ePositionZ[i])
        {
           zMin = ePositionZ[i]; 
        }
        if(zMax<ePositionZ[i])
        {
           zMax = ePositionZ[i]; 
        } 
    }
    zMinMax[0] = zMin;
    zMinMax[1] = zMax;
    return zMinMax;
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
    vector<double> zMinMax = GetZMinMax();
    zMin = zMinMax[0];
    zMax = zMinMax[1];

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
            if(eSurive[partID]==0) continue;
           
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
            if(eSurive[partID]==0) continue;
            //eMomentumX[partID] += srWakePoten[0][i];            //rad
            //eMomentumY[partID] += srWakePoten[1][i];            //rad    
            eMomentumZ[partID] += srWakePoten[2][i];            //rad
        }
    }
}
