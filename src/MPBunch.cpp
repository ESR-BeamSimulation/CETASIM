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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <memory.h>
#include <fftw3.h>
#include <complex.h>

using namespace std;
using std::vector;
using std::complex;
using v1d= vector<double> ;
using v2d= vector<vector<double> > ;
using v3d= vector<vector<vector<double> > >;
using v1i= vector<int> ;
using v2i= vector<vector<int> > ;

MPBunch::MPBunch()
{
}

MPBunch::~MPBunch()
{
    delete wakePotenFromBBI;          
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
    densProfVsBin.resize(bunchBinNumberZ);
    densProfVsBinAnalytical.resize(bunchBinNumberZ+1);
    hamiltonPotenWell.resize(bunchBinNumberZ+1);

    cavVolInfoVsLongBins.resize(inputParameter.ringParRf->resNum);
    cavForceInfoVsLongBins.resize(inputParameter.ringParRf->resNum);

    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        cavVolInfoVsLongBins[i].resize(bunchBinNumberZ+1);
        cavForceInfoVsLongBins[i].resize(bunchBinNumberZ+1);
    }

    // for space charge simulation
    // vector<vector<vector<double> > > slicedPIC2DBeam
    // slicedBunchDisInfo.resize(11);
    // slicedBunchChargeInfo.resize(11);
    // slicedBunchPartiIndex.resize(11);
    // for (int i=0;i<slicedBunchDisInfo.size();i++)
    // {
    //     slicedBunchDisInfo[i].resize(6);   
    // }
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
  
    // for(int i=0;i<macroEleNumPerBunch;i++)
    // {
        
    //     ePositionX[i] = (i+1) * 1.E-5;
    //     ePositionY[i] = (i+1) * 1.E-5; 
    //     eMomentumX[i] = 0.E0;
    //     eMomentumY[i] = 0.E0;
    //     ePositionZ[i] = 0.E0;
    //     eMomentumZ[i] = 0.E0;
    // }

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
    rangeUP  = 20.E0;    // seems like 10 rms emittance truncated..

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

    double etax   = latticeInterActionPoint.twissDispX[0];
    double etaxp  = latticeInterActionPoint.twissDispPX[0];  // \frac{disP}{ds} 
    double etay   = latticeInterActionPoint.twissDispY[0];
    double etayp  = latticeInterActionPoint.twissDispPY[0];  // \frac{disP}{ds} 

    // substracut the orbit shift due to the dispersion // Ref. Elegant ILMatrix Eq(56)
    // for(int i=0;i<macroEleNumPerBunch;i++)
    // {
    //     ePositionX[i]  -=  etax  * eMomentumZ[i];   // [m] 
    //     ePositionY[i]  -=  etay  * eMomentumZ[i];   // [m]
    //     eMomentumX[i]  -=  etaxp * eMomentumZ[i];   // [rad]
    //     eMomentumY[i]  -=  etayp * eMomentumZ[i];   // [rad]
    // }

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

    // // add the orbit shift due to the dispersion // Ref. Elegant ILMatrix Eq(56)
    // for(int i=0;i<macroEleNumPerBunch;i++)
    // {
    //     ePositionX[i]  +=  etax  * eMomentumZ[i];   // [m] 
    //     ePositionY[i]  +=  etay  * eMomentumZ[i];   // [m]
    //     eMomentumX[i]  +=  etaxp * eMomentumZ[i];   // [rad]
    //     eMomentumY[i]  +=  etayp * eMomentumZ[i];   // [rad]
    // }

    GetZMinMax();
    GetEigenEmit(latticeInterActionPoint);
}

void MPBunch::GetEigenEmit(const LatticeInterActionPoint &latticeInterActionPoint)
{
    double etax   = latticeInterActionPoint.twissDispX[0];
    double etaxp  = latticeInterActionPoint.twissDispPX[0];  // \frac{disP}{ds} 
    double etay   = latticeInterActionPoint.twissDispY[0];
    double etayp  = latticeInterActionPoint.twissDispPY[0];  // \frac{disP}{ds} 

    vector<double> x  = ePositionX;
    vector<double> y  = ePositionY;
    vector<double> xp = eMomentumX;
    vector<double> yp = eMomentumY;
    vector<double> zp = eMomentumZ;

    // substracut the orbit shift due to the dispersion // Ref. Elegant ILMatrix Eq(56)
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        x[i]   -=  etax  * zp[i];   // [m] 
        y[i]   -=  etay  * zp[i];   // [m]
        xp[i]  -=  etaxp * zp[i];   // [rad]
        yp[i]  -=  etayp * zp[i];   // [rad]
    }

    double averX  = std::accumulate(x.begin(),  x.end(),  0);
    double averY  = std::accumulate(y.begin(),  y.end(),  0);
    double averXP = std::accumulate(xp.begin(), xp.end(), 0);
    double averYP = std::accumulate(yp.begin(), yp.end(), 0);

    double aver2[10] = {0};  // x_x, x_xp, x_y,x_yp/ Ref. Lars, PRAB,16,044201 
    

    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        x[i]   -=  averX;   
        y[i]   -=  averY;   // [m]
        xp[i]  -=  averXP;   // [rad]
        yp[i]  -=  averYP;   // [rad]

        aver2[0] +=  x[i] *   x[i];
        aver2[1] +=  x[i] *  xp[i];
        aver2[2] +=  x[i] *   y[i];
        aver2[3] +=  x[i] *  yp[i];
        aver2[4] += xp[i] *  xp[i];
        aver2[5] += xp[i] *   y[i];
        aver2[6] += xp[i] *  yp[i];
        aver2[7] +=  y[i] *   y[i];
        aver2[8] +=  y[i] *  yp[i];
        aver2[9] += yp[i] *  yp[i]; 
    }

    for(int i=0;i<10;i++) aver2[i] /= macroEleNumPerBunch;
    
    // aver2[0] = 8.57;
    // aver2[1] = -4.34;
    // aver2[2] = -3.83;
    // aver2[3] = -1.15;
    // aver2[4] = 3.35;
    // aver2[5] = -0.54;
    // aver2[6] = 1.52;
    // aver2[7] = 11.2;
    // aver2[8] = -3.05;
    // aver2[9] = 1.87;

    gsl_matrix *sympleMarixJ4 = gsl_matrix_alloc (4, 4);
    gsl_matrix_set_zero(sympleMarixJ4);
    gsl_matrix_set(sympleMarixJ4,0,1,1);
    gsl_matrix_set(sympleMarixJ4,1,0,-1);
    gsl_matrix_set(sympleMarixJ4,2,3,1);
    gsl_matrix_set(sympleMarixJ4,3,2,-1);


    gsl_matrix *simgMatrix    = gsl_matrix_alloc (4, 4);
 
    gsl_matrix_set(simgMatrix,0,0,aver2[0]);
    gsl_matrix_set(simgMatrix,0,1,aver2[1]);
    gsl_matrix_set(simgMatrix,0,2,aver2[2]);
    gsl_matrix_set(simgMatrix,0,3,aver2[3]);

    gsl_matrix_set(simgMatrix,1,0,aver2[1]);
    gsl_matrix_set(simgMatrix,1,1,aver2[4]);
    gsl_matrix_set(simgMatrix,1,2,aver2[5]);
    gsl_matrix_set(simgMatrix,1,3,aver2[6]);    

    gsl_matrix_set(simgMatrix,2,0,aver2[2]);
    gsl_matrix_set(simgMatrix,2,1,aver2[5]);
    gsl_matrix_set(simgMatrix,2,2,aver2[7]);
    gsl_matrix_set(simgMatrix,2,3,aver2[8]);  

    gsl_matrix_set(simgMatrix,3,0,aver2[3]);
    gsl_matrix_set(simgMatrix,3,1,aver2[6]);
    gsl_matrix_set(simgMatrix,3,2,aver2[8]);
    gsl_matrix_set(simgMatrix,3,3,aver2[9]);  

  

    gsl_matrix *matrixCJ  = gsl_matrix_alloc (4, 4);
    gsl_matrix *matrixCJ2 = gsl_matrix_alloc (4, 4);
    gsl_matrix_set_zero(matrixCJ);
    gsl_matrix_set_zero(matrixCJ2);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,simgMatrix,sympleMarixJ4,0.0,matrixCJ); 
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,matrixCJ,  matrixCJ,     0.0,matrixCJ2); 

    
    double trCJ2=0;
    for(int i=0;i<4;i++) trCJ2 += gsl_matrix_get(matrixCJ2,i,i);
    
    double  det = get_det(simgMatrix);

    eigenEmitX = sqrt(- trCJ2 + sqrt( trCJ2 * trCJ2 -  16 * det ) )/2;
    eigenEmitY = sqrt(- trCJ2 - sqrt( trCJ2 * trCJ2 -  16 * det ) )/2; 
    emitXY4D = sqrt(sqrt(det));


    gsl_matrix_free(matrixCJ);
    gsl_matrix_free(matrixCJ2);
    gsl_matrix_free(simgMatrix);
    gsl_matrix_free(sympleMarixJ4);


    // get the coupling alpha at x-y plane // Andrea Franchi  phd theiss, 2006 Fig.9.5
    // double emitXY = sqrt(aver2[0] * aver2[7] -  aver2[2] * aver2[2] );   // m 
    // xyCouplingAlpha = - aver2[2] / emitXY;  // rad 

    // gsl_matrix *xyEllipse = gsl_matrix_alloc (2, 2); 
    // gsl_matrix_set(xyEllipse,0,0,aver2[0]);
    // gsl_matrix_set(xyEllipse,0,1,aver2[2]);
    // gsl_matrix_set(xyEllipse,1,0,aver2[2]);
    // gsl_matrix_set(xyEllipse,1,1,aver2[7]);

    // gsl_vector *eval = gsl_vector_alloc (2);
    // gsl_matrix *evec = gsl_matrix_alloc (2, 2);
     
    // gsl_eigen_symmv_workspace * w =  gsl_eigen_symmv_alloc (2);
    // gsl_eigen_symmv (xyEllipse, eval, evec, w);
    // gsl_eigen_symmv_sort (eval, evec,GSL_EIGEN_SORT_ABS_ASC);

    // double v1y = gsl_matrix_get(evec,0,1);
    // double v1x = gsl_matrix_get(evec,0,0);
    
    // xyCouplingAlpha = atan2(v1y,v1x) - PI /2;
    
    // gsl_eigen_symmv_free (w); 
    // gsl_vector_free(eval);
    // gsl_matrix_free(evec);
    // gsl_matrix_free(xyEllipse);

    double a = aver2[0];
    double b = aver2[2];
    double c = aver2[7];
    double lambda =(a + c) / 2 + sqrt( (a - c) * (a - c) / 4 + b * b  );
    if(b==0)
    {
        xyCouplingAlpha = a>=c? 0: PI/2;
    } 
    else
    {
        xyCouplingAlpha = atan2((lambda-a),b);
    }
}

void MPBunch::SetSlicedBunchInfo(int const nz)
{   
     //(1) set the slicedBuncDisInfo, 2.5D PIC model (-5 simgaz, 5*simgaz) 11 slices as default.   
    slicedBunchDisInfo    = v3d(nz,v2d(6,v1d()));  
    slicedBunchChargeInfo = v2d(nz,v1d());
    slicedBunchPartiIndex = v2i(nz,v1i());

    GetZMinMax();  
    double zMin = zMinCurrentTurn;
    double zMax = zMaxCurrentTurn;  
    double dz = (zMax - zMin) /  (nz - 1);

    //int range = 10;
    //double dz = 2 * range * rmsBunchLength / (nz - 1);
    
    int sliceIndex;
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        //sliceIndex =  floor( (ePositionZ[i] -  (zAver - range * rmsBunchLength)) / dz);
        sliceIndex =  floor( (ePositionZ[i] -  zMin) / dz);
        sliceIndex = sliceIndex < 0    ? 0    :  sliceIndex;
        sliceIndex = sliceIndex > nz-1 ? nz-1 :  sliceIndex; 
			
        slicedBunchDisInfo[sliceIndex][0].push_back(ePositionX[i]);
        slicedBunchDisInfo[sliceIndex][1].push_back(ePositionY[i]);
        slicedBunchDisInfo[sliceIndex][2].push_back(ePositionZ[i]);
        slicedBunchDisInfo[sliceIndex][3].push_back(eMomentumX[i]);
        slicedBunchDisInfo[sliceIndex][4].push_back(eMomentumY[i]);
        slicedBunchDisInfo[sliceIndex][5].push_back(eMomentumZ[i]);
        slicedBunchChargeInfo[sliceIndex].push_back(macroEleCharge * ElectronCharge / dz);  // [C/m]
        slicedBunchPartiIndex[sliceIndex].push_back(i);
    }

}

void MPBunch::UpdateTransverseMomentumBEModel(vector<vector<double>>  &particles, vector<double> eCharge, const double ds)
{
    int dims = 2; 
    vector<double> beamCen(dims,0.E0);
    vector<double> beamRms(dims,0.E0);
    for (int plane=0; plane<dims;plane++)
    { 
        vector<double> pos = particles[plane];
        beamCen[plane] = accumulate(pos.begin(), pos.end(), 0.E0) / pos.size();
        beamRms[plane] = 0.E0;
        for(int i = 0; i<pos.size();i++ )
        {
            beamRms[plane] += pow(pos[i] - beamCen[plane] ,2);
        }
        beamRms[plane] =  sqrt(beamRms[plane] / pos.size() );
    }
    double posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy;
    rmsRxTemp = beamRms[0];
    rmsRyTemp = beamRms[1];
    
    double lambda = accumulate(eCharge.begin(), eCharge.end(), 0.E0); // [C/m] 
    
    vector<vector<double>> partExEyField(2,v1d(particles[0].size(),0.E0));
    //(1) get the exEy field 
    for(int i=0;i<particles[0].size();++i)
    {
        posx = particles[0][i] - beamCen[0] ;  // reference to the center of the beam
        posy = particles[1][i] - beamCen[1] ;  // reference to the center of the beam

        if( (rmsRxTemp - rmsRyTemp) / rmsRyTemp > 1.e-4) 
        {
            BassettiErskine1(posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy);
        }
        else if ( (rmsRxTemp - rmsRyTemp)/rmsRyTemp < -1.e-4 )
        {
            BassettiErskine1(posy,posx,rmsRyTemp,rmsRxTemp,tempFy,tempFx); //tempFx, [1/m]
        }
        else
        {
            GaussianField(posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy);    //tempFx, [1/m]
        }
        partExEyField[0][i] = -lambda / PI / 2.0 / Epsilon * tempFx;          // [C/m] /  (C / (V m)) * [1/m] -> [V/m]
        partExEyField[1][i] = -lambda / PI / 2.0 / Epsilon * tempFy;          // [V/m]
    }

    //(2) kick beam transversely  
    double gamma0 = electronEnergy / ElectronMassEV; 
    double beta0  = sqrt(1 - 1 / pow(gamma0, 2));
    double p0     = beta0 * gamma0;
    double px,py,pz,p,gamma,beta,deltaPx,deltaPy;
    
    // ofstream fout("testBE.sdds");
    for(int i=0;i<particles[0].size();++i)
    {
        pz = (1 + particles[5][i]) * p0;
        px = pz * particles[1][i];
        py = pz * particles[3][i];
        p = sqrt(pow(pz,2) + pow(px,2) + pow(py,2) );
        gamma = sqrt(1 + pow(p,2) );
        beta  = p / gamma;

        deltaPx =  partExEyField[0][i] / pow(gamma,2) * ElectronCharge * ds / (beta * CLight) / (ElectronMass * p * CLight );
        deltaPy =  partExEyField[1][i] / pow(gamma,2) * ElectronCharge * ds / (beta * CLight) / (ElectronMass * p * CLight );

        particles[3][i] += deltaPx;
        particles[4][i] += deltaPy;

        // fout<<setw(15)<< particles[0][i] 
        //     <<setw(15)<< particles[1][i]
        //     <<setw(15)<< deltaPx
        //     <<setw(15)<< deltaPy
        //     <<endl; 
    }

    // cout<<__LINE__<<__FILE__<<endl;
    // getchar();
}


void MPBunch::UpdateTransverseMomentumLinear(vector<vector<double>>  &particles, vector<double> eCharge, const double ds)
{
    int dims = 2; 
    vector<double> beamCen(dims,0.E0);
    vector<double> beamRms(dims,0.E0);
    for (int plane=0; plane<dims;plane++)
    { 
        vector<double> pos = particles[plane];
        beamCen[plane] = accumulate(pos.begin(), pos.end(), 0.E0) / pos.size();
        beamRms[plane] = 0.E0;
        for(int i = 0; i<pos.size();i++ )
        {
            beamRms[plane] += pow(pos[i] - beamCen[plane] ,2);
        }
        beamRms[plane] =  sqrt(beamRms[plane] / pos.size() );
    }
    double posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy;
    rmsRxTemp = beamRms[0];   // assume unifrom in (-2,2) rms beam size // same idea as KV linear approximation
    rmsRyTemp = beamRms[1];   
    
    double lambda = accumulate(eCharge.begin(), eCharge.end(), 0.E0); // [C/m] 
    
    vector<vector<double>> partExEyField(2,v1d(particles[0].size(),0.E0));
    //(1) get the exEy field 
    for(int i=0;i<particles[0].size();++i)
    {
        posx = particles[0][i] - beamCen[0] ;  // reference to the center of the beam
        posy = particles[1][i] - beamCen[1] ;  // reference to the center of the beam

        tempFx = posx / (rmsRxTemp * (rmsRxTemp + rmsRyTemp));
        tempFy = posy / (rmsRyTemp * (rmsRxTemp + rmsRyTemp));  

        partExEyField[0][i] = lambda / PI / Epsilon * tempFx;          // [C/m] /  (C / (V m)) * [1/m] -> [V/m]
        partExEyField[1][i] = lambda / PI / Epsilon * tempFy;          // [V/m]
    }


    //(2) kick beam transversely  
    double gamma0 = electronEnergy / ElectronMassEV; 
    double beta0  = sqrt(1 - 1 / pow(gamma0, 2));
    double p0     = beta0 * gamma0;
    double px,py,pz,p,gamma,beta,deltaPx,deltaPy;
    
    // ofstream fout("testLinar.sdds");
    for(int i=0;i<particles[0].size();++i)
    {
        pz = (1 + particles[5][i]) * p0;
        px = pz * particles[1][i];
        py = pz * particles[3][i];
        p = sqrt(pow(pz,2) + pow(px,2) + pow(py,2) );
        gamma = sqrt(1 + pow(p,2) );
        beta  = p / gamma;

        deltaPx =  partExEyField[0][i] / pow(gamma,2) * ElectronCharge * ds / (beta * CLight) / (ElectronMass * p * CLight );
        deltaPy =  partExEyField[1][i] / pow(gamma,2) * ElectronCharge * ds / (beta * CLight) / (ElectronMass * p * CLight );

        particles[3][i] += deltaPx;
        particles[4][i] += deltaPy;

        // fout<<setw(15)<< particles[0][i] 
        //     <<setw(15)<< particles[1][i]
        //     <<setw(15)<< deltaPx
        //     <<setw(15)<< deltaPy
        //     <<endl; 
    }

    // cout<<__LINE__<<__FILE__<<endl;
    // getchar();
}




void MPBunch::BunchMomentumUpdateDueToSpaceChargeAnalytical(LatticeInterActionPoint &latticeInterActionPoint, int k, const ReadInputSettings &inputParameter)
{
    // 2.5D model for simulaiton with the ideal model
    int nz = inputParameter.ringRun->scMeshNum[2];
    int scFlag = inputParameter.ringRun->spaceChargeFlag;
    SetSlicedBunchInfo(nz);
    vector<vector<double> > particles(6,v1d());
    vector<double> eCharge;
    vector<int> partiIndexInSlice;
    double ds = latticeInterActionPoint.interactionLength[k];

    for (int slice = 0; slice<slicedBunchDisInfo.size(); ++slice )
    {
        particles =  slicedBunchDisInfo[slice];
        eCharge   =  slicedBunchChargeInfo[slice];
        partiIndexInSlice = slicedBunchPartiIndex[slice];             
        if(particles[0].size()<2) continue;
        
        if(scFlag==1) UpdateTransverseMomentumBEModel(particles,eCharge,ds);
        if(scFlag==3) UpdateTransverseMomentumLinear(particles,eCharge,ds);

        int partIndex;
        for(int i=0; i<partiIndexInSlice.size(); ++i)
        {
            partIndex = partiIndexInSlice[i];
            eMomentumX[partIndex] = particles[3][i];
            eMomentumY[partIndex] = particles[4][i];
        } 
    }
    // clear the slicebunchInfo
    for (int slice=0; slice< slicedBunchDisInfo.size(); ++slice )
    {
        slicedBunchChargeInfo[slice].clear();
        slicedBunchPartiIndex[slice].clear();
        for(int i=0; i<6; ++i)  slicedBunchDisInfo[slice][i].clear();
    }   


}

void MPBunch::BunchMomentumUpdateDueToSpaceChargePIC(PIC3D &picSCBeam3D,LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    // generate the mesh according to the rms beam size, and keep mesh size for in each slice. 
    double rmsXY[2]  = {rmsRx,rmsRy};
    // double rmsXY[2]  = {1.E-4,1.E-4};
    double averXY[2] = {xAver,yAver};
    if(rmsXY[0]==0 || rmsXY[1]==0)
    {
        cerr<<"rms beam size is zero, not able to generate mesh. RMSX:"<<rmsXY[0]<<"RMSY:"<<rmsXY[1]<<endl;
        exit(0);
    }

    // 2.5D model for simulaiton--since sigmaz~1.E-3, sigmax~1.E-5, sigmay~1.E-6, very un-symmetryic 
    vector<vector<double> > particles(6,v1d());
    vector<double> eCharge;
    vector<int> partiIndexInSlice;
    // mesh is generatote in XY accoding the rms beam size.
    picSCBeam3D.slicedBunchPIC2D->Set2DMesh(rmsXY,averXY);
    int nz = picSCBeam3D.numberOfGrid[2];
    SetSlicedBunchInfo(nz);  // sliced longitudianl bunch profile

    double ds = latticeInterActionPoint.interactionLength[k];

    for (int slice = 0; slice<slicedBunchDisInfo.size(); ++slice )
    {
        particles =  slicedBunchDisInfo[slice];
        eCharge   =  slicedBunchChargeInfo[slice];
        partiIndexInSlice = slicedBunchPartiIndex[slice];
             
        if(particles[0].size()<2) continue;
        picSCBeam3D.slicedBunchPIC2D->Set2DRho(particles,eCharge);
        picSCBeam3D.slicedBunchPIC2D->Set2DPhi(); 
        picSCBeam3D.slicedBunchPIC2D->Set2DEField();
        picSCBeam3D.slicedBunchPIC2D->SetPartSCField();
        picSCBeam3D.slicedBunchPIC2D->UpdateTransverseMomentum(particles, ds);

        // update the particle momentum
        int partIndex;
        for(int i=0; i<partiIndexInSlice.size(); ++i)
        {
            partIndex = partiIndexInSlice[i];
            eMomentumX[partIndex] = particles[3][i];
            eMomentumY[partIndex] = particles[4][i];
        } 
    }
    
    // // clear the slicebunchInfo
    for (int slice=0; slice< slicedBunchDisInfo.size(); ++slice )
    {
        slicedBunchChargeInfo[slice].clear();
        slicedBunchPartiIndex[slice].clear();
        for(int i=0; i<6; ++i)  slicedBunchDisInfo[slice][i].clear();
    }    

    // cout<<__LINE__<<__FILE__<<endl;
    // getchar();



    // 3D solver directly 
    // double ds = latticeInterActionPoint.interactionLength[k];
    // vector<vector<double> > particles(6,vector<double>() ) ;
    // for(int i=0;i<macroEleNumPerBunch;i++)
    // {
    //     if(eSurive[i]!=0) continue;  // only surive particles are set for PIC solver
    //     particles[0].push_back(ePositionX[i]);
    //     particles[1].push_back(ePositionY[i]);
    //     particles[2].push_back(ePositionZ[i]);
    //     particles[3].push_back(eMomentumX[i]); // px / p0
    //     particles[4].push_back(eMomentumY[i]); // py / p0
    //     particles[5].push_back(eMomentumZ[i]); // deltaP / p0;
    // }
    

    // cout<<macroEleCharge * ElectronCharge * macroEleNumPerBunch<<endl;
    // picSCBeam3D.Set3DMesh(particles);
    // picSCBeam3D.Set3DRho(particles, macroEleCharge * ElectronCharge);
    // picSCBeam3D.Set3DPhi1();
    // picSCBeam3D.Set3DEField();
    // picSCBeam3D.SetPartSCField(particles);
    // picSCBeam3D.UpdatMomentum(particles,ds);



    // // copy coodinates banck to member variables momentum
    // int j=0;
    // for(int i=0;i<particles[0].size();i++)
    // {
    //     if(eSurive[i]!=0)
    //     {
    //         j++;
    //         continue;  // reset the momentum of the particiles.
    //     } 
    //     eMomentumX[i+j] = particles[3][i];
    //     eMomentumY[i+j] = particles[4][i];
    //     eMomentumZ[i+j] = particles[5][i];
    // }   
}

void MPBunch::SSIonBunchInteractionPIC(BeamIon2DPIC &beamIon2DPIC, LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    // (1) get the E fiel from electron beam
    double ds = latticeInterActionPoint.interactionLength[k];
    double circRing = latticeInterActionPoint.circRing;
    vector<vector<double> > particles(4,vector<double>()) ;  // [x,y,px,py]
    vector<double> eCharge;

    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        if(eSurive[i]!=0) continue;  // only surive particles are set for PIC solver
        particles[0].push_back(ePositionX[i]);
        particles[1].push_back(ePositionY[i]);
        // particles[2].push_back(eMomentumX[i]); // px / p0
        // particles[3].push_back(eMomentumY[i]); // py / p0
        eCharge.push_back(macroEleCharge * ElectronCharge / circRing);   // [C/m] 2D electron beam line density
    }

    beamIon2DPIC.pic2DBeam.Set2DMesh(particles);
    beamIon2DPIC.pic2DBeam.Set2DRho(particles,eCharge); 
    beamIon2DPIC.pic2DBeam.Set2DPhi(); 
    beamIon2DPIC.pic2DBeam.Set2DEField();

    // here to compare with result from BasstiErskine formular 
    // int nx = beamIon2DPIC.pic2DBeam.numberOfGrid[0];
    // int ny = beamIon2DPIC.pic2DBeam.numberOfGrid[1]; 
    // cout<<nx<<" "<<ny<<endl;   
    // ofstream fout("test1.sdds");
    // double a = beamIon2DPIC.pic2DBeam.beamrmssize[0] ;
    // double b = beamIon2DPIC.pic2DBeam.beamrmssize[1] ;
    // double x,y,tempFx,tempFy;
    // for(int i=0;i<beamIon2DPIC.pic2DBeam.numberOfGrid[0];i++)
    // {
    //     for(int j=0;j<beamIon2DPIC.pic2DBeam.numberOfGrid[1];j++)
    //     {
    //         x = beamIon2DPIC.pic2DBeam.meshx[0] + i * beamIon2DPIC.pic2DBeam.meshWidth[0];
    //         y = beamIon2DPIC.pic2DBeam.meshy[0] + j * beamIon2DPIC.pic2DBeam.meshWidth[1];
            
    //         BassettiErskine1(x,y,a,b,tempFx,tempFy);

    //         fout<<setw(15)<< x / a
    //             <<setw(15)<< y / b
    //             <<setw(15)<< beamIon2DPIC.pic2DBeam.rho[i][j]
    //             <<setw(15)<< beamIon2DPIC.pic2DBeam.phi[i][j]
    //             <<setw(15)<< beamIon2DPIC.pic2DBeam.ex[i][j]
    //             <<setw(15)<< beamIon2DPIC.pic2DBeam.ey[i][j]
    //             <<setw(15)<< -tempFx * eCharge[0] * macroEleNumPerBunch / 2 / PI / Epsilon   
    //             <<setw(15)<< -tempFy * eCharge[0] * macroEleNumPerBunch / 2 / PI / Epsilon  
    //             <<endl;
    //     }
    // }
    // cout<<__LINE__<<__FILE__<<endl;
    // getchar();

    
    // (2) get the E fiel from ions on the mesh  
    vector<vector<double> > ions(4,vector<double>()) ;
    vector<double> ionCharge;
    for(int p=0;p<latticeInterActionPoint.gasSpec;p++)
    {
        for(int j=0;j<latticeInterActionPoint.ionAccumuNumber[k][p];j++)
        {
            ions[0].push_back(latticeInterActionPoint.ionAccumuPositionX[k][p][j]);
            ions[1].push_back(latticeInterActionPoint.ionAccumuPositionY[k][p][j]);
            // ions[2].push_back(latticeInterActionPoint.macroIonCharge[k][p] * ElectronCharge)
            // ions[3].push_back(latticeInterActionPoint.ionMassNumber[p]);
            ionCharge.push_back(latticeInterActionPoint.macroIonCharge[k][p] * ElectronCharge);   // [C] 
        }
    }

    beamIon2DPIC.pic2DIon.Set2DMesh(ions);
    beamIon2DPIC.pic2DIon.Set2DRho(ions,ionCharge);
    beamIon2DPIC.pic2DIon.Set2DPhi();
    beamIon2DPIC.pic2DIon.Set2DEField();

    // to be updated...
    //(3) Get the momentum kick of electron due to ions 
    vector<vector<double> > eFieldPart = beamIon2DPIC.pic2DIon.GetPartSCField(particles);
    int j;
    for(int i=0;i<particles[0].size();i++)
    {
        if(eSurive[i]!=0)
        {
            j++;
            continue;  // reset the momentum of the particiles.
        } 
        eFxDueToIon[i+j] = eFieldPart[0][i];
        eFxDueToIon[i+j] = eFieldPart[1][i];
    }   
    
    //(4) Get the momentum kick of ions due to elecrons
    vector<vector<double> > eFieldIons = beamIon2DPIC.pic2DBeam.GetPartSCField(ions);
    int count = 0;
    for(int p=0;p<latticeInterActionPoint.gasSpec;p++)
    {
        for(int j=0;j<latticeInterActionPoint.ionAccumuNumber[k][p];j++)
        {            
            latticeInterActionPoint.ionAccumuFx[k][p][j]= eFieldIons[0][count];
            latticeInterActionPoint.ionAccumuFy[k][p][j]= eFieldIons[1][count];
            count++;
        }
    }

    cout<<__LINE__<<__FILE__<<endl;
    getchar();

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



void MPBunch::BunchMomentumUpdateDuetoRFRigid(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
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
                            cavityResonator.resonatorVec[j].resQualityQ0,0.E0) * electronNumPerBunch * ElectronCharge;  // [Volt]

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
                if(j==0) eMomentumZ[index] -= u0       / electronBeamEnergy / pow(rBeta,2); 
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

void MPBunch::BunchMomentumUpdateDueToRFModeStable(const ReadInputSettings &inputParameter,Resonator &resonator, int resIndex)
{

    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    double fRF         = f0 * ringHarmH;
    int bunchBinNumberZ = inputParameter.ringParRf->rfBunchBinNum;
    double dzBin = (zMaxCurrentTurn - zMinCurrentTurn) / bunchBinNumberZ;
    double dtBin = dzBin / CLight / rBeta;
    double resHarm = resonator.resHarm;
    double resFre  = resonator.resFre;
    double tB, deltaL, cPsi;
    tB = deltaL = cPsi = 0.E0;

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
    for(int i=0;i<bunchBinNumberZ;i++)
    {
        densProfVsBin[i] = histoParIndex[i].size() / double(macroEleNumPerBunch) /  dzBin;
    }
    

    double cavVoltageReal = 0.E0;
    double cavVoltageImag = 0.E0;
    double induceVolReal = 0.E0;
    double induceVolImag = 0.E0;
    double genVolReal = 0.E0;
    double genVolImag = 0.E0;
    complex<double> vb0=(0,0);
    complex<double> cavVoltage=(0.E0,0.E0);
    complex<double> genVoltage=(0.E0,0.E0);
    complex<double> selfLossVolAccume=(0.E0,0.E0);
    int particleInBunch=0;

    // rigid beam Vb is calculated once per bunch
    vb0  = complex<double>(-1 * 2 * PI * resFre * resonator.resShuntImpRs / resonator.resQualityQ0, 0.E0) * electronNumPerBunch * ElectronCharge;  // [Volt]
 
    genVoltage = resonator.resGenVol;     
    cavVoltage = resonator.vbAccum + genVoltage;

    for(int k=0;k<histoParIndex.size();k++) 
    {                    
        // loop for particle in each bin          
        for(int i=0;i<histoParIndex[k].size();i++)
        {
            int index          = histoParIndex[k][i];
            eMomentumZ[index] += (cavVoltage * exp( - li * posZBins[k]  / CLight / rBeta * 2. * PI * double(resHarm) * fRF)).real() / electronBeamEnergy / pow(rBeta,2);
            eMomentumZ[index] += vb0.real()/2.0  / electronBeamEnergy / pow(rBeta,2);
        }

        particleInBunch   += histoParIndex[k].size();
        selfLossVolAccume += vb0/2.0           * double(histoParIndex[k].size());
        cavVoltageReal    += cavVoltage.real() * double(histoParIndex[k].size());
        cavVoltageImag    += cavVoltage.imag() * double(histoParIndex[k].size());
        genVolReal        += genVoltage.real() * double(histoParIndex[k].size());
        genVolImag        += genVoltage.imag() * double(histoParIndex[k].size());
        induceVolReal     += (resonator.vbAccum + vb0).real() * double(histoParIndex[k].size()); 
        induceVolImag     += (resonator.vbAccum + vb0).imag() * double(histoParIndex[k].size());

        // beam induced voltage retotate and decay bin-by-bin            
        // tB     = dtBin;
        // deltaL = tB  / resonator.tF;
        // cPsi   = 2.0 * PI *  resonator.resFre * tB;
        // resonator.vbAccum += vb0;
        // resonator.vbAccum *= exp(- deltaL ) * exp (li * cPsi);
    }
    resonator.vbAccum += vb0;

    cavVoltageReal /= double(particleInBunch);
    cavVoltageImag /= double(particleInBunch);
    genVolReal     /= double(particleInBunch);
    genVolImag     /= double(particleInBunch);
    induceVolReal  /= double(particleInBunch);
    induceVolImag  /= double(particleInBunch); 

    bunchRFModeInfo->cavVolBunchCen[resIndex]       = complex<double>(cavVoltageReal,cavVoltageImag);
    bunchRFModeInfo->genVolBunchAver[resIndex]      = complex<double>(genVolReal    ,genVolImag    ); 
    bunchRFModeInfo->induceVolBunchCen[resIndex]    = complex<double>(induceVolReal ,induceVolImag); 
    bunchRFModeInfo->selfLossVolBunchCen[resIndex]  =  selfLossVolAccume / double(particleInBunch) ;
}


void MPBunch::BunchMomentumUpdateDueToRFMode(const ReadInputSettings &inputParameter,Resonator &resonator, int resIndex)
{
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    double fRF         = f0 * ringHarmH;
    int bunchBinNumberZ = inputParameter.ringParRf->rfBunchBinNum;
    double dzBin = (zMaxCurrentTurn - zMinCurrentTurn) / bunchBinNumberZ;
    double dtBin = dzBin / CLight / rBeta;
    double resHarm = resonator.resHarm;
    double resFre  = resonator.resFre;
    double tB, deltaL, cPsi;
    tB = deltaL = cPsi = 0.E0;

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
    for(int i=0;i<bunchBinNumberZ;i++)
    {
        densProfVsBin[i] = histoParIndex[i].size() / double(macroEleNumPerBunch) / dzBin;
    }
   
    double cavVoltageReal = 0.E0;
    double cavVoltageImag = 0.E0;
    double induceVolReal = 0.E0;
    double induceVolImag = 0.E0;
    double genVolReal = 0.E0;
    double genVolImag = 0.E0;
    complex<double> vb0=(0,0);
    complex<double> cavVoltage=(0.E0,0.E0);
    complex<double> genVoltage=(0.E0,0.E0);
    complex<double> selfLossVolAccume=(0.E0,0.E0);
    int particleInBunch=0;

    // loop for bin-by-bin in one bunch, bins is alined from head to tail and each bin excite beam induced voltage itself

    for(int k=0;k<histoParIndex.size();k++) 
    {            
        genVoltage = resonator.resGenVol * exp( - li * posZBins[k]  / CLight / rBeta * 2. * PI * double(resHarm) * fRF); 
        cavVoltage = resonator.vbAccum   +  genVoltage;
        
        vb0  = complex<double>(-1 * 2 * PI * resFre * resonator.resShuntImpRs / resonator.resQualityQ0, 0.E0) * macroEleCharge*double(histoParIndex[k].size()) * ElectronCharge;  // [Volt]

        // loop for particle in each bin          
        for(int i=0;i<histoParIndex[k].size();i++)
        {
            int index          = histoParIndex[k][i];
            eMomentumZ[index] += cavVoltage.real() / electronBeamEnergy / pow(rBeta,2);
            eMomentumZ[index] += vb0.real()/2.0    / electronBeamEnergy / pow(rBeta,2);
        }

        particleInBunch   += histoParIndex[k].size();
        selfLossVolAccume += vb0/2.0           * double(histoParIndex[k].size());
        cavVoltageReal    += cavVoltage.real() * double(histoParIndex[k].size());
        cavVoltageImag    += cavVoltage.imag() * double(histoParIndex[k].size());
        genVolReal        += genVoltage.real() * double(histoParIndex[k].size());
        genVolImag        += genVoltage.imag() * double(histoParIndex[k].size());
        induceVolReal     += (resonator.vbAccum + vb0).real() * double(histoParIndex[k].size()); 
        induceVolImag     += (resonator.vbAccum + vb0).imag() * double(histoParIndex[k].size());

        // beam induced voltage retotate and decay bin-by-bin            
        tB     = dtBin;
        deltaL = tB  / resonator.tF;
        cPsi   = 2.0 * PI *  resonator.resFre * tB;
        resonator.vbAccum += vb0;
        resonator.vbAccum *= exp(- deltaL ) * exp (li * cPsi);
    }
    
    
    cavVoltageReal /= double(particleInBunch);
    cavVoltageImag /= double(particleInBunch);
    genVolReal     /= double(particleInBunch);
    genVolImag     /= double(particleInBunch);
    induceVolReal  /= double(particleInBunch);
    induceVolImag  /= double(particleInBunch); 
    selfLossVolAccume /= double(particleInBunch);

    bunchRFModeInfo->cavVolBunchCen[resIndex]       = complex<double>(cavVoltageReal,cavVoltageImag);
    bunchRFModeInfo->genVolBunchAver[resIndex]      = complex<double>(genVolReal    ,genVolImag    ); 
    bunchRFModeInfo->induceVolBunchCen[resIndex]    = complex<double>(induceVolReal ,induceVolImag) + selfLossVolAccume ; 
    bunchRFModeInfo->selfLossVolBunchCen[resIndex]  =  selfLossVolAccume ;

}


void MPBunch::BunchMomentumUpdateDuetoRFBinByBin(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    // Ref. bunch.h that ePositionZ = - ePositionT * c. head pariticles: deltaT<0, ePositionZ[i]>0.
    
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
                if(j==0) eMomentumZ[index] -= u0       / electronBeamEnergy / pow(rBeta,2); 
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
                    if(j==0)  eMomentumZ[index] -= u0            /electronBeamEnergy / pow(rBeta,2); 
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

    // zMin = *min_element(ePositionZ.begin(),ePositionZ.end());
    // zMax = *max_element(ePositionZ.begin(),ePositionZ.end());
 

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


void MPBunch::BBImpBunchInteractionTD(const ReadInputSettings &inputParameter, const BoardBandImp &boardBandImp,const LatticeInterActionPoint &latticeInterActionPoint,
                                      vector<vector<double>> wakePoten)
{
    double rBeta              = inputParameter.ringParBasic->rBeta;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    int ringHarmH             = inputParameter.ringParRf->ringHarm;
    double rfLen              = inputParameter.ringParBasic->t0 / ringHarmH * CLight; 
    double dzBin              = abs(wakePoten[0][1] - wakePoten[0][0]);

    double zMinBin            = - 1 * rfLen / 2;
    
    int nBinWake = wakePoten[0].size();
    int nBinBunchDen = int(1 * rfLen / dzBin);    // bunch density profile -- from [-rfLen/2,rfLen/2]
    
    double temp[5]; // wake force in wz wdx wdy wqx, wqy
    
    double densityProfile[nBinBunchDen];
    for(int i=0;i<nBinBunchDen;i++)densityProfile[i]=0;
    
    vector<vector<int>> histoParIndex;
    histoParIndex.resize(nBinBunchDen);
    for(int i=0;i<ePositionZ.size();i++)
    {
        int index = int( (ePositionZ[i] - zMinBin ) / dzBin );          // head particle index stores in histoParIndex[i]. Smaller i represents head particle.  
        histoParIndex[index].push_back(i);
        densityProfile[index] +=  macroEleCharge * ElectronCharge ;  // [C]
    } 
    GetSmoothedBunchProfileGassionFilter(densityProfile,nBinBunchDen);

    vector<double> averXAlongBunch(nBinBunchDen,0);
    vector<double> averYAlongBunch(nBinBunchDen,0);
    
    for(int i=0;i<nBinBunchDen;i++)
    {
        for(int j=0;j<histoParIndex[i].size();j++)
        {
            int partID = histoParIndex[i][j];
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

    // ofstream fout("wakePotenTimeDomain_fromGreenFun.dat");
    // fout<<"SDDS1"<<endl;
    // fout<<"&column name=z,              units=m,              type=float,  &end" <<endl;
    // fout<<"&column name=profile,        units=m,              type=float,  &end" <<endl;
    // fout<<"&column name=indVZ,                                type=float,  &end" <<endl;
    // fout<<"&column name=indVDX,                               type=float,  &end" <<endl;
    // fout<<"&column name=indVDY,                               type=float,  &end" <<endl;
    // fout<<"&column name=indVQX,                                type=float,  &end" <<endl;
    // fout<<"&column name=indVQY,                                type=float,  &end" <<endl;
    // fout<<"&data mode=ascii, &end"                                               <<endl;
    // fout<<"! page number " << 1 <<endl;
    // fout<<nBinBunchDen<<endl;

    int partID, tij;
    for(int i=0;i<nBinBunchDen;i++) // ith bins in front of jth bins. how ith bin affect jth bin
    {    
        for(int j=0;j<5;j++)temp[j] = 0.E0;
                
        for(int j=i;j<nBinBunchDen;j++)          // bins in front of the ith bin
        {
            tij = j - i;
            if(tij>nBinWake) break;
            if(inputParameter.ringBBImp->impedSimFlag[0]==1) temp[0] -= wakePoten[1][tij] * densityProfile[j];                      //[V/C] * [C] = V
            if(inputParameter.ringBBImp->impedSimFlag[1]==1) temp[1] -= wakePoten[2][tij] * densityProfile[j] * averXAlongBunch[j]; //[V/C] * [C] * [m] = [V m]
            if(inputParameter.ringBBImp->impedSimFlag[2]==1) temp[2] -= wakePoten[3][tij] * densityProfile[j] * averYAlongBunch[j];
            if(inputParameter.ringBBImp->impedSimFlag[3]==1) temp[3] -= wakePoten[4][tij] * densityProfile[j] ; //[V/C] * [C]  = [V]
            if(inputParameter.ringBBImp->impedSimFlag[4]==1) temp[4] -= wakePoten[5][tij] * densityProfile[j] ;
        }

        for(int k=0;k<histoParIndex[i].size();k++) // Eq3.7 and Eq3.49 
        {                    
            partID = histoParIndex[i][k];
            if(inputParameter.ringBBImp->impedSimFlag[0]==1) eMomentumZ[partID] += temp[0] / electronBeamEnergy / pow(rBeta,2);
            if(inputParameter.ringBBImp->impedSimFlag[1]==1) eMomentumX[partID] += temp[1] / electronBeamEnergy / pow(rBeta,2) / latticeInterActionPoint.twissBetaX[0];
            if(inputParameter.ringBBImp->impedSimFlag[2]==1) eMomentumY[partID] += temp[2] / electronBeamEnergy / pow(rBeta,2) / latticeInterActionPoint.twissBetaY[0];
            if(inputParameter.ringBBImp->impedSimFlag[3]==1) eMomentumX[partID] += temp[3] / electronBeamEnergy / pow(rBeta,2) / latticeInterActionPoint.twissBetaX[0] * ePositionX[partID];
            if(inputParameter.ringBBImp->impedSimFlag[4]==1) eMomentumY[partID] += temp[4] / electronBeamEnergy / pow(rBeta,2) / latticeInterActionPoint.twissBetaY[0] * ePositionY[partID];
        }
    
        
        // fout<<setw(15)<<left<<i * dzBin + zMinBin
        //     <<setw(15)<<left<<-densityProfile[i]
        //     <<setw(15)<<left<<temp[0]
        //     <<setw(15)<<left<<temp[1]
        //     <<setw(15)<<left<<temp[2]
        //     <<setw(15)<<left<<temp[3]
        //     <<setw(15)<<left<<temp[4]
        //     <<endl;
          
    }

   
    // fout.close();
    // cout<<__LINE__<<__FILE__<<endl;
    // getchar();

}


void MPBunch::BBImpBunchInteraction(const ReadInputSettings &inputParameter, const BoardBandImp &boardBandImp,const LatticeInterActionPoint &latticeInterActionPoint)
{
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    double zMinBin = -boardBandImp.zMax;
    double dzBin   =  boardBandImp.dz;
    int   nBins    =  profileForBunchBBImp.size(); // bins for bunch profile, that is 2*boardBandImp.zZImp.size() - 1
   
    int doTracking=0;
    for (int i=0;i<5;i++)
    {
        doTracking = doTracking || inputParameter.ringBBImp->impedSimFlag[i];
    }
    if(doTracking==0) return;
    

    fftw_complex *r2cout  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nBins /2 +1) );  // the same size as boardBandImp.zZImp.size() // bunch specturm
    fftw_complex *c2rin   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nBins /2 +1) );  // the same size as boardBandImp.zZImp.size()

    double temp[nBins];
    int    partBinIndex[macroEleNumPerBunch];

    // get the smoothed longi-BunchCharge-profile 
    fill(profileForBunchBBImp.begin(),profileForBunchBBImp.end(),0); // array[ 2* nBins -1 ] to store the profile
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        partBinIndex[i]  =  (ePositionZ[i] - zMinBin + dzBin / 2 ) / dzBin; 
        int index = partBinIndex[i];
        profileForBunchBBImp[index] +=  macroEleCharge * ElectronCharge / dzBin *  CLight;  // [C/s]  
    }
    GetSmoothedBunchProfileGassionFilter(inputParameter, boardBandImp);
    for(int i=0;i<nBins;i++) temp[i] = profileForBunchBBImp[i];

    fftw_plan p = fftw_plan_dft_r2c_1d(nBins, temp, r2cout, FFTW_ESTIMATE);    // ffw, according to g++/nvcc complile flag is FFTW_MEASURE/FFTW_ESTIMATE
    fftw_execute(p);

    
    // (0) in zz  directoin --------------------------------------------------- Eq3.10 
    if(inputParameter.ringBBImp->impedSimFlag[0]==1)  // beam longitudianl impedance interactoin
    {
        for(int i=0;i<boardBandImp.zZImp.size();i++)
        {
            complex<double> indVF = complex<double> (r2cout[i][0], r2cout[i][1] ) * boardBandImp.zZImp[i];
            c2rin[i][0] = indVF.real();
            c2rin[i][1] = indVF.imag();                              // [C/s] * [Ohm] = [V]   
        }
        p = fftw_plan_dft_c2r_1d(nBins, c2rin, temp, FFTW_ESTIMATE);
        fftw_execute(p);   
        for(int i=0;i<nBins;i++)  wakePotenFromBBI->wakePotenZ[i] = - temp[i] / nBins;
        
        for(int i=0;i<macroEleNumPerBunch;i++)
        {
            int index  = partBinIndex[i] ; 
            eMomentumZ[i] += wakePotenFromBBI->wakePotenZ[index] / electronBeamEnergy / pow(rBeta,2);
        }
    }

    // (3) in zQx  directoin --------------------------------------------------- 
    if(inputParameter.ringBBImp->impedSimFlag[3]==1)  // beam ZQx interactoin
    { 
        for(int i=0;i<boardBandImp.zZImp.size();i++)
        {
            complex<double> indVF = complex<double> (r2cout[i][0], r2cout[i][1] ) * boardBandImp.zQxImp[i] * li;
            c2rin[i][0] = indVF.real();
            c2rin[i][1] = indVF.imag();                              // [C/s] * [Ohm] = [V]   
        }
        p = fftw_plan_dft_c2r_1d(nBins, c2rin, temp, FFTW_ESTIMATE);
        fftw_execute(p);   
        for(int i=0;i<nBins;i++)  wakePotenFromBBI->wakePotenQx[i] = temp[i] / nBins;
        
        for(int i=0;i<macroEleNumPerBunch;i++)
        {
            int index  = partBinIndex[i] ; 
            eMomentumX[i] += wakePotenFromBBI->wakePotenQx[index] / electronBeamEnergy / pow(rBeta,2) *  ePositionX[i] / latticeInterActionPoint.twissBetaX[0];
        }
    }


    // (4) in zQy  directoin --------------------------------------------------- 
    if(inputParameter.ringBBImp->impedSimFlag[4]==1)  // beam ZQx interactoin
    {
        for(int i=0;i<boardBandImp.zZImp.size();i++)
        {
            complex<double> indVF = complex<double> (r2cout[i][0], r2cout[i][1] ) * boardBandImp.zQyImp[i] * li;
            c2rin[i][0] = indVF.real();
            c2rin[i][1] = indVF.imag();                              // [C/s] * [Ohm] = [V]   
        }
        p = fftw_plan_dft_c2r_1d(nBins, c2rin, temp, FFTW_ESTIMATE);
        fftw_execute(p);   
        for(int i=0;i<nBins;i++)  wakePotenFromBBI->wakePotenQy[i] = temp[i] / nBins;
        
        for(int i=0;i<macroEleNumPerBunch;i++)
        {
            int index  = partBinIndex[i] ; 
            eMomentumY[i] += wakePotenFromBBI->wakePotenQy[index] / electronBeamEnergy / pow(rBeta,2) *  ePositionY[i] / latticeInterActionPoint.twissBetaY[0];
        }
    }

    

    // (1) in zDx  directoin --------------------------------------------------- //Alex Chao Eq 3.51
    if(inputParameter.ringBBImp->impedSimFlag[1]==1)   // beam with ZDx
    {
        fill(profileForBunchBBImp.begin(),profileForBunchBBImp.end(),0); // array[ 2* nBins -1 ] to store the profile

        for(int i=0;i<macroEleNumPerBunch;i++)
        {
            partBinIndex[i]  =  (ePositionZ[i] - zMinBin + dzBin / 2 ) / dzBin; 
            int index        = partBinIndex[i];
            profileForBunchBBImp[index] +=  macroEleCharge * ElectronCharge / dzBin *  CLight * ePositionX[i];  // [C/s]  
        }
        GetSmoothedBunchProfileGassionFilter(inputParameter, boardBandImp);
        for(int i=0;i<nBins;i++) temp[i] = profileForBunchBBImp[i];

        p = fftw_plan_dft_r2c_1d(nBins, temp, r2cout, FFTW_ESTIMATE);    // ffw, according to g++/nvcc complile flag is FFTW_MEASURE/FFTW_ESTIMATE
        fftw_execute(p);
        for(int i=0;i<boardBandImp.zZImp.size();i++)
        {
            complex<double> indVF = complex<double> (r2cout[i][0], r2cout[i][1] ) * boardBandImp.zDxImp[i] * li;
            c2rin[i][0] = indVF.real();
            c2rin[i][1] = indVF.imag();                              // [C/s] * [Ohm] = [V]   
        }
        p = fftw_plan_dft_c2r_1d(nBins, c2rin, temp, FFTW_ESTIMATE);
        fftw_execute(p);   
        for(int i=0;i<nBins;i++)  wakePotenFromBBI->wakePotenDx[i] = temp[i] / nBins;           
            
        for(int i=0;i<macroEleNumPerBunch;i++)
        {
            int index      = partBinIndex[i] ; 
            eMomentumX[i] += wakePotenFromBBI->wakePotenDx[index] / electronBeamEnergy / pow(rBeta,2) / latticeInterActionPoint.twissBetaX[0];
        }
    }
     
     
     // (2) in zDy  directoin --------------------------------------------------- 
    if(inputParameter.ringBBImp->impedSimFlag[2]==1)   // beam beam with ZDy
    {
        fill(profileForBunchBBImp.begin(),profileForBunchBBImp.end(),0); // array[ 2* nBins -1 ] to store the profile

        for(int i=0;i<macroEleNumPerBunch;i++)
        {
            partBinIndex[i]  =  (ePositionZ[i] - zMinBin + dzBin / 2 ) / dzBin; 
            int index        = partBinIndex[i];
            profileForBunchBBImp[index] +=  macroEleCharge * ElectronCharge / dzBin *  CLight * ePositionY[i];  // [C/s]  
        }
        GetSmoothedBunchProfileGassionFilter(inputParameter, boardBandImp);
        for(int i=0;i<nBins;i++) temp[i] = profileForBunchBBImp[i];

        p = fftw_plan_dft_r2c_1d(nBins, temp, r2cout, FFTW_ESTIMATE);    // ffw, according to g++/nvcc complile flag is FFTW_MEASURE/FFTW_ESTIMATE
        fftw_execute(p);
        for(int i=0;i<boardBandImp.zZImp.size();i++)
        {
            complex<double> indVF = complex<double> (r2cout[i][0], r2cout[i][1] ) * boardBandImp.zDyImp[i] * li;
            c2rin[i][0] = indVF.real();
            c2rin[i][1] = indVF.imag();                              // [C/s] * [Ohm] = [V]   
        }
        p = fftw_plan_dft_c2r_1d(nBins, c2rin, temp, FFTW_ESTIMATE);
        fftw_execute(p);   
        for(int i=0;i<nBins;i++)  wakePotenFromBBI->wakePotenDy[i] = temp[i] / nBins;
            
        for(int i=0;i<macroEleNumPerBunch;i++)
        {
            int index      = partBinIndex[i] ; 
            eMomentumY[i] += wakePotenFromBBI->wakePotenDy[index] / electronBeamEnergy / pow(rBeta,2) / latticeInterActionPoint.twissBetaY[0];
        }
    }

    

    fftw_free(r2cout);
    fftw_free(c2rin);
    fftw_destroy_plan(p);


    // test the wakefield her. 
    // ofstream fout("wakePotenFreqDomain_23mm.dat");
    // fout<<"SDDS1"<<endl;
    // fout<<"&column name=z,              units=m,              type=float,  &end" <<endl;
    // fout<<"&column name=profile,                              type=float,  &end" <<endl;
    // fout<<"&column name=indVZ,                                type=float,  &end" <<endl;
    // fout<<"&column name=indVDX,                               type=float,  &end" <<endl;
    // fout<<"&column name=indVDY,                               type=float,  &end" <<endl;
    // fout<<"&column name=indVQX,                                type=float,  &end" <<endl;
    // fout<<"&column name=indVQY,                                type=float,  &end" <<endl;
    // fout<<"&data mode=ascii, &end"                                               <<endl;
    // fout<<"! page number " << 1 <<endl;
    // fout<<nBins<<endl;

    // for(int i=0; i<nBins;i++)
    // {
    //     fout<<setw(15)<<left<<boardBandImp.binPosZ[i]
    //         <<setw(15)<<left<<profileForBunchBBImp[i]
    //         <<setw(15)<<left<<wakePotenFromBBI->wakePotenZ[i]
    //         <<setw(15)<<left<<wakePotenFromBBI->wakePotenDx[i]
    //         <<setw(15)<<left<<wakePotenFromBBI->wakePotenDy[i]
    //         <<setw(15)<<left<<wakePotenFromBBI->wakePotenQx[i]
    //         <<setw(15)<<left<<wakePotenFromBBI->wakePotenQy[i]
    //         <<endl;
    // }
    // cout<<"wrie file"<<endl;
    // fout.close();
    // getchar();
}


void MPBunch::GetSmoothedBunchProfileGassionFilter(const ReadInputSettings &inputParameter, const BoardBandImp &boardBandImp )
{
    double zMaxBin =  boardBandImp.zMax;
    double zMinBin = -boardBandImp.zMax;
    double dzBin   =  boardBandImp.dz;
    int    nBins   =  boardBandImp.nBins;   // from 0 to max
    double rBeta   = inputParameter.ringParBasic->rBeta;
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double rfLen   = inputParameter.ringParBasic->t0 / ringHarmH * CLight;
    int index; 
   
    int indexStart =   ( -rfLen / 2. - zMinBin + dzBin / 2 ) / dzBin;
    int indexEnd   =   (  rfLen / 2. - zMinBin + dzBin / 2 ) / dzBin;
    int N = indexEnd - indexStart;
  

    double tempProfile[N];
    for(int i=0;i<N;i++) tempProfile[i] = profileForBunchBBImp[i+indexStart];  

    GetSmoothedBunchProfileGassionFilter(tempProfile,N);

    for(int i=0;i<N;i++) profileForBunchBBImp[i+indexStart] = tempProfile[i];

}

void MPBunch::GetSmoothedBunchProfileGassionFilter(double *profile,int N)
{
     
    int K=51;
    const double alpha[3] = { 10, 3.0, 10.0 };   /* alpha values */
    gsl_vector *x = gsl_vector_alloc(N);         /* input vector */
    gsl_vector *y1 = gsl_vector_alloc(N);        /* filtered output vector for alpha1 */
    gsl_vector *k1 = gsl_vector_alloc(K);        /* Gaussian kernel for alpha1 */
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_filter_gaussian_workspace *gauss_p = gsl_filter_gaussian_alloc(K);
    double sum = 0.0;

    // /* generate input signal */
    for(int i=0; i<N; i++) gsl_vector_set(x, i, profile[i]);
    
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
        profile[i] =  y1i; 
    }
    sum0 = sum0 / sum1;

    for(int i=0; i<N; i++) profile[i] *= sum0;
    
    gsl_vector_free(x);
    gsl_vector_free(y1);
    gsl_vector_free(k1);
    gsl_rng_free(r);
    gsl_filter_gaussian_free(gauss_p);
}
