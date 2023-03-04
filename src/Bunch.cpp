//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             


#include "Bunch.h"
#include "Global.h"
#include "Faddeeva.h"
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


Bunch::Bunch()
{
}

Bunch::~Bunch()
{
    delete bunchRFModeInfo; 
    delete haissinski;  
}


void Bunch::Initial(const  ReadInputSettings &inputParameter)
{    
    macroEleNumPerBunch = inputParameter.ringBunchPara->macroEleNumPerBunch;      // macroEleNumPerBunch=1 -- electron beam weak-strong model                                                                           
    electronEnergy      = inputParameter.ringParBasic->electronBeamEnergy;
    rmsEnergySpread     = inputParameter.ringBunchPara-> rmsEnergySpread;
    rmsBunchLength      = inputParameter.ringBunchPara-> rmsBunchLength;
    rGamma              = inputParameter.ringParBasic->rGamma;
    rBeta               = inputParameter.ringParBasic->rBeta;
    // current             = inputParameter.ringBunchPara->current; 
    // in sp model, it does not change during the tracking 
    emittanceZ      = inputParameter.ringBunchPara->emittanceZ;
    emittanceX      = inputParameter.ringBunchPara->emittanceX;
    emittanceY      = inputParameter.ringBunchPara->emittanceY;
     
    electronNumPerBunch = current  / inputParameter.ringParBasic->f0  / ElectronCharge;  // [dimmenless]
    macroEleCharge      = electronNumPerBunch / macroEleNumPerBunch;                     // [dimmenless]

    lRWakeForceAver.resize(3,0E0);
    
    ePositionX.resize(macroEleNumPerBunch,0E0);
    ePositionY.resize(macroEleNumPerBunch,0E0);
    ePositionZ.resize(macroEleNumPerBunch,0E0);
    eMomentumX.resize(macroEleNumPerBunch,0E0);
    eMomentumY.resize(macroEleNumPerBunch,0E0);
    eMomentumZ.resize(macroEleNumPerBunch,0E0);
    eFxDueToIon.resize(macroEleNumPerBunch,0E0);
    eFyDueToIon.resize(macroEleNumPerBunch,0E0);
    eFzDueToIon.resize(macroEleNumPerBunch,0E0);
    eSurive.resize(macroEleNumPerBunch,0);
    
    
    bunchRFModeInfo->cavVolBunchCen.resize(inputParameter.ringParRf->resNum);
    bunchRFModeInfo->genVolBunchAver.resize(inputParameter.ringParRf->resNum);
    bunchRFModeInfo->induceVolBunchCen.resize(inputParameter.ringParRf->resNum);
    bunchRFModeInfo->selfLossVolBunchCen.resize(inputParameter.ringParRf->resNum);
    bunchRFModeInfo->genIgBunchAver.resize(inputParameter.ringParRf->resNum); 
    bunchRFModeInfo->genPower.resize(inputParameter.ringParRf->resNum);
    bunchRFModeInfo->beamPower.resize(inputParameter.ringParRf->resNum);
    bunchRFModeInfo->cavPower.resize(inputParameter.ringParRf->resNum);
    bunchRFModeInfo->refPower.resize(inputParameter.ringParRf->resNum);

    
    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {        
        bunchRFModeInfo->genVolBunchAver[i]     =complex<double>(0.e0,0.e0);
        bunchRFModeInfo->induceVolBunchCen[i]   =complex<double>(0.e0,0.e0);
        bunchRFModeInfo->selfLossVolBunchCen[i] =complex<double>(0.e0,0.e0);
        bunchRFModeInfo->cavVolBunchCen[i]      =complex<double>(0.e0,0.e0);              
    }

    // refer to Maro -- 2015 PRAB 18 031001 --Eq(24) to get the coupled bunch mode growth rate. 
    xyzHistoryDataToFit.resize(3);
    for (int i=0;i<3;i++)
    {
        xyzHistoryDataToFit[i].resize(inputParameter.ringRun->bunchInfoPrintInterval);    
    }
}


void Bunch::BassettiErskine1(double posx,double posy,double rmsRxTemp, double rmsRyTemp,double &tempFx,double &tempFy)
{
    complex<double> z1(0.E0,0.E0);
    complex<double> z2(0.E0,0.E0);
    complex<double> w1(0.E0,0.E0);
    complex<double> w2(0.E0,0.E0);
    double z3;
    double w3;

    complex<double> w0(0.E0,0.E0);

    double sigma=0.E0;
    double ryOverRx=0.E0;
    double tempPosix=0.E0;
    double tempPosiy=0.E0;
    tempFx=0;
    tempFy=0;

    tempPosix    =  abs(posx);
    tempPosiy    =  abs(posy);

    sigma    = sqrt(2*pow(rmsRxTemp,2)-2*pow(rmsRyTemp,2));
    ryOverRx = rmsRyTemp/rmsRxTemp;

    double coeffBE;
    coeffBE    = -1 * sqrt(PI)/sigma;

    z1           =  {tempPosix/sigma,tempPosiy/sigma};
    w1           =  Faddeeva::w(z1);

    z2           =  {ryOverRx*tempPosix/sigma,tempPosiy/sigma/ryOverRx};
    w2           =  Faddeeva::w(z2);

    z3           =  - pow(tempPosix/rmsRxTemp,2)/2 - pow(tempPosiy/rmsRyTemp,2)/2;
    w3           =  - exp(z3);

    w0           =  coeffBE * (w1 + w3 * w2 );

    tempFx       = w0.imag();
    tempFy       = w0.real();


    if(posx<=0)
    {
        tempFx  = -tempFx;
    }

    if(posy<=0)
    {
        tempFy = -tempFy;
    }
}

void Bunch::MarkLostParticle(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint)
{
    int ringHarm        = inputParameter.ringParRf->ringHarm;
    double t0           = inputParameter.ringParBasic->t0;

    for(int i=0;i<macroEleNumPerBunch;i++)
    {    
        if( abs(ePositionZ[i]) > t0 * CLight / ringHarm ) eSurive[i] = 2;   // loss in longitudianl
    }

}

void Bunch::GaussianField(double posx,double posy,double rmsRxTemp, double rmsRyTemp,double &tempFx,double &tempFy)
{
    double r2 = pow(posx,2) + pow(posy,2);
    double sigma  = rmsRxTemp;
    double sigma2 = sigma * sigma;

    tempFx = - posx / r2 * (1 - exp(- r2 /2/sigma2 ) );
    tempFy = - posy / r2 * (1 - exp(- r2 /2/sigma2 ) );
}

void Bunch::BunchTransferDueToIon(const LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        eMomentumX[i] +=  eFxDueToIon[i]  ;    // rad
        eMomentumY[i] +=  eFyDueToIon[i]  ;    // rad
    }
}

void Bunch::BunchTransferDueToLatticeOneTurnT66(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint)
{
    // get the twiss parameters of from lattice
    double alphax,betax,alphay,betay,gammax,gammay,etax,etaxp,etay,etayp;
    alphax = latticeInterActionPoint.twissAlphaX[0];
    alphay = latticeInterActionPoint.twissAlphaY[0];
    betax  = latticeInterActionPoint.twissBetaX[0];
    betay  = latticeInterActionPoint.twissBetaY[0];
    
    etax   = latticeInterActionPoint.twissDispX[0];
    etaxp  = latticeInterActionPoint.twissDispPX[0];  // \frac{disP}{ds} 
    etay   = latticeInterActionPoint.twissDispY[0];
    etayp  = latticeInterActionPoint.twissDispPY[0];  // \frac{disP}{ds} 

    double *alphac = inputParameter.ringParBasic->alphac;
    double *aDTX  = inputParameter.ringParBasic->aDTX;
    double *aDTY  = inputParameter.ringParBasic->aDTY;
    double *aDTXY = inputParameter.ringParBasic->aDTXY;
    
    double circRing = inputParameter.ringParBasic->circRing;
    double nux = inputParameter.ringParBasic->workQx;
    double nuy = inputParameter.ringParBasic->workQy;
    double chromx = inputParameter.ringParBasic->chrom[0];
    double chromy = inputParameter.ringParBasic->chrom[1];
    double t0         = inputParameter.ringParBasic->t0;
    double eta        = inputParameter.ringParBasic->eta;

    gammax = (1 + pow(alphax,2))/betax;
    gammay = (1 + pow(alphay,2))/betay;


    // here follow what elegant approaches to build the R66 matix for each particle
    gsl_matrix *ILmatrix  = gsl_matrix_alloc (6, 6);
    gsl_matrix *cord0     = gsl_matrix_alloc (6, 1);
    gsl_matrix *cord1     = gsl_matrix_alloc (6, 1);
    gsl_matrix_set_zero(ILmatrix);
    gsl_matrix_set_zero(cord0);
    gsl_matrix_set_zero(cord1);

    //ref. Elegant ILMatrix element
    double xPosN,yPosN,xMomN,yMomN;
    double ampX,ampY;
    double nuxtmp,nuytmp,tmp;
    double phix,phiy;
    // generate the transfer matrix for each particle in bunch 
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        // only the first order is kept
        xPosN  = ePositionX[i] - etax  * eMomentumZ[i];   // [m] 
        yPosN  = ePositionY[i] - etay  * eMomentumZ[i];   // [m]
        xMomN  = eMomentumX[i] - etaxp * eMomentumZ[i];   // [rad]
        yMomN  = eMomentumY[i] - etayp * eMomentumZ[i];   // [rad]


        ampX   = ( pow(xPosN,2) + pow( alphax * xPosN + betax * xMomN,2) ) / betax;  //[m]
        ampY   = ( pow(yPosN,2) + pow( alphay * yPosN + betay * yMomN,2) ) / betay;  //[m]

        nuxtmp = nux + chromx * eMomentumZ[i] +  aDTX[0] * ampX + aDTX[1] * pow(ampX,2) / 2 + aDTXY[0] * ampX * ampY;  //[]
        nuytmp = nuy + chromy * eMomentumZ[i] +  aDTY[0] * ampY + aDTY[1] * pow(ampY,2) / 2 + aDTXY[1] * ampX * ampY;  //[]

        phix = 2 * PI * nuxtmp ;
        phiy = 2 * PI * nuytmp ;

        tmp = cos(phix) + alphax * sin(phix);   gsl_matrix_set(ILmatrix,0,0,tmp);   //R11
        tmp =             betax  * sin(phix);   gsl_matrix_set(ILmatrix,0,1,tmp);   //R12
        tmp =           - gammax * sin(phix);   gsl_matrix_set(ILmatrix,1,0,tmp);   //R21
        tmp = cos(phix) - alphax * sin(phix);   gsl_matrix_set(ILmatrix,1,1,tmp);   //R22

        tmp = cos(phiy) + alphay * sin(phiy);   gsl_matrix_set(ILmatrix,2,2,tmp);   //R33
        tmp =             betay  * sin(phiy);   gsl_matrix_set(ILmatrix,2,3,tmp);   //R34
        tmp =           - gammay * sin(phiy);   gsl_matrix_set(ILmatrix,3,2,tmp);   //R43
        tmp = cos(phiy) - alphay * sin(phiy);   gsl_matrix_set(ILmatrix,3,3,tmp);   //R44

        tmp = etax  - etax  * cos(phix) - (alphax * etax + betax * etaxp) *  sin(phix);                                    gsl_matrix_set(ILmatrix,0,5,tmp);  //R16
        tmp = etay  - etay  * cos(phiy) - (alphay * etay + betay * etayp) *  sin(phiy);                                    gsl_matrix_set(ILmatrix,2,5,tmp);  //R36  

        tmp = etaxp - etaxp * cos(phix)  + ( etax  + pow(alphax,2) * etax +  alphax * betax * etaxp) *  sin(phix) / betax; gsl_matrix_set(ILmatrix,1,5,tmp); //R26    
        tmp = etayp - etayp * cos(phiy)  + ( etay  + pow(alphay,2) * etay +  alphay * betay * etayp) *  sin(phiy) / betay; gsl_matrix_set(ILmatrix,3,5,tmp); //R46

        tmp = -etaxp + etaxp * cos(phix) + ( etax  + pow(alphax,2) * etax +  alphax * betax * etaxp) *  sin(phix) / betax; gsl_matrix_set(ILmatrix,4,0,tmp); //R51
        tmp = -etayp + etayp * cos(phiy) + ( etay  + pow(alphay,2) * etay +  alphay * betay * etayp) *  sin(phiy) / betay; gsl_matrix_set(ILmatrix,4,2,tmp); //R53 

        tmp =  etax  - etax  * cos(phix) + ( alphax * etax +  betax * etaxp) * sin(phix);                                  gsl_matrix_set(ILmatrix,4,1,tmp); //R52 
        tmp =  etay  - etay  * cos(phiy) + ( alphay * etay +  betay * etayp) * sin(phiy);                                  gsl_matrix_set(ILmatrix,4,3,tmp); //R54  

        // how to set R55, R56, R65, R66. checked it with yong-chul recently.
        gsl_matrix_set(ILmatrix,4,4,1); //R55
        gsl_matrix_set(ILmatrix,5,5,1); //R66

        // set the right particle coordinate (x_{beta} ) for one turn transfer. Kick the particle position here, but not in rf cavity subroutine.
        // also make the R56 kick extra due to high order momentum compact factor, ampllidude depended length...    
        ePositionX[i] = xPosN; 
        ePositionY[i] = yPosN;
        eMomentumX[i] = xMomN;
        eMomentumY[i] = yMomN;
        
        gsl_matrix_set(cord0,0,0,ePositionX[i]);
        gsl_matrix_set(cord0,1,0,eMomentumX[i]);
        gsl_matrix_set(cord0,2,0,ePositionY[i]);
        gsl_matrix_set(cord0,3,0,eMomentumY[i]);
        gsl_matrix_set(cord0,4,0,ePositionZ[i]);
        gsl_matrix_set(cord0,5,0,eMomentumZ[i]);


        gsl_matrix_mul(ILmatrix,cord0,cord1);
        ePositionX[i]  = gsl_matrix_get(cord1,0,0);
        eMomentumX[i]  = gsl_matrix_get(cord1,1,0);
        ePositionY[i]  = gsl_matrix_get(cord1,2,0);
        eMomentumY[i]  = gsl_matrix_get(cord1,3,0);
        ePositionZ[i]  = gsl_matrix_get(cord1,4,0);
        eMomentumZ[i]  = gsl_matrix_get(cord1,5,0);
        
        // longitudinal pozition updated in one turn -- momentum compactor of the whole ring. 
        ePositionZ[i] -= circRing * (alphac[0] * eMomentumZ[i]  - pow(alphac[1] * eMomentumZ[i] ,2) + pow(alphac[2] * eMomentumZ[i] ,3) ) ;    
    }

    gsl_matrix_free(ILmatrix);
    gsl_matrix_free (cord0);
    gsl_matrix_free (cord1);

}


void Bunch::BunchTransferDueToLatticeT(const LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    double xtemp=0.E0;
    double ytemp=0.E0;
    double xPtemp=0.E0;
    double yPtemp=0.E0;
    double ztemp=0.E0;
    double zPtemp=0.E0;
   
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        // refers to SY. Lee Eq. 2.67

        xtemp  = latticeInterActionPoint.xTransferMatrix[k][0] * ePositionX[i]
               + latticeInterActionPoint.xTransferMatrix[k][1] * eMomentumX[i];

        xPtemp = latticeInterActionPoint.xTransferMatrix[k][2] * ePositionX[i]
               + latticeInterActionPoint.xTransferMatrix[k][3] * eMomentumX[i];

        ytemp  = latticeInterActionPoint.yTransferMatrix[k][0] * ePositionY[i]
               + latticeInterActionPoint.yTransferMatrix[k][1] * eMomentumY[i];

        yPtemp = latticeInterActionPoint.yTransferMatrix[k][2] * ePositionY[i]
               + latticeInterActionPoint.yTransferMatrix[k][3] * eMomentumY[i];

        ePositionX[i] = xtemp;
        ePositionY[i] = ytemp;
        eMomentumX[i] = xPtemp;
        eMomentumY[i] = yPtemp;

        double lossTemp =pow(ePositionX[i]/latticeInterActionPoint.pipeAperatureX[k],2) + pow(ePositionY[i]/latticeInterActionPoint.pipeAperatureY[k],2) ; 
        if(lossTemp >1) eSurive[i] = 1;   // loss in transverse
    }
    
}



void Bunch::BunchSynRadDamping(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint)
{
//Note: the SynRadDamping and excitation is follow Yuan ZHang's approaches. The results is not consistent with approaches used in MBtrack.
//This section is only used for transverse direction. The longitudinla is dealt with in longitudinal tracking
//Bunch::BunchTransferDueToLatticeL() -- 2021-12-30

    // in the unit of number of truns for synchRadDampTime setting.
    vector<double> synchRadDampTime;
    synchRadDampTime.resize(3);
    for(int i=0;i<synchRadDampTime.size();i++)
    {
        synchRadDampTime[i] = inputParameter.ringParBasic->synchRadDampTime[i];
    }

    int k=0;
    //	 set the J2_sympeletic matrix
    gsl_matrix * sympleMarixJ = gsl_matrix_alloc (2, 2);
    gsl_matrix_set_zero(sympleMarixJ);
    gsl_matrix_set (sympleMarixJ, 0, 1, 1);
    gsl_matrix_set (sympleMarixJ, 1, 0,-1);

    //	// set the H matrix
    gsl_matrix * dispMatrixH  = gsl_matrix_alloc (6, 6);
    gsl_matrix * dispMatrixHX = gsl_matrix_alloc (2, 2);
    gsl_matrix * dispMatrixHY = gsl_matrix_alloc (2, 2);

    gsl_matrix * H31   = gsl_matrix_alloc (2, 2);
    gsl_matrix * H32   = gsl_matrix_alloc (2, 2);


    gsl_matrix_set_identity(dispMatrixH);
    gsl_matrix_set_zero(dispMatrixHX);
    gsl_matrix_set_zero(dispMatrixHY);
    gsl_matrix_set_zero(H31);
    gsl_matrix_set_zero(H32);



    gsl_matrix_set(dispMatrixHX,0,1,latticeInterActionPoint.twissDispX[k]);
    gsl_matrix_set(dispMatrixHX,1,1,latticeInterActionPoint.twissDispPX[k]);

    gsl_matrix_set(dispMatrixHY,0,1,latticeInterActionPoint.twissDispY[k]);
    gsl_matrix_set(dispMatrixHY,1,1,latticeInterActionPoint.twissDispPY[k]);


    gsl_matrix * tempMatrix1 = gsl_matrix_alloc (2, 2);
    gsl_matrix * tempMatrix2 = gsl_matrix_alloc (2, 2);
    gsl_matrix_set_zero(tempMatrix1);
    gsl_matrix_set_zero(tempMatrix2);


    //	// get H31 sub_matrix;
    gsl_matrix_transpose(dispMatrixHX);

    gsl_matrix_mul(dispMatrixHX,sympleMarixJ,tempMatrix1);
    gsl_matrix_mul(sympleMarixJ,tempMatrix1,tempMatrix2);
    gsl_matrix_scale (tempMatrix2, -1);
    gsl_matrix_memcpy (H31, tempMatrix2);

    gsl_matrix_transpose(dispMatrixHX);



    //	// get H32 sub_matrix;
    gsl_matrix_transpose(dispMatrixHY);

    gsl_matrix_mul(dispMatrixHY,sympleMarixJ,tempMatrix1);
    gsl_matrix_mul(sympleMarixJ,tempMatrix1,tempMatrix2);
    gsl_matrix_scale (tempMatrix2, -1);

    gsl_matrix_memcpy (H32, tempMatrix2);
    gsl_matrix_transpose(dispMatrixHY);


    //	// set H Marix
    gsl_matrix_set(dispMatrixH,0,4,-gsl_matrix_get(dispMatrixHX,0,0));
    gsl_matrix_set(dispMatrixH,0,5,-gsl_matrix_get(dispMatrixHX,0,1));
    gsl_matrix_set(dispMatrixH,1,4,-gsl_matrix_get(dispMatrixHX,1,0));
    gsl_matrix_set(dispMatrixH,1,5,-gsl_matrix_get(dispMatrixHX,1,1));

    gsl_matrix_set(dispMatrixH,2,4,-gsl_matrix_get(dispMatrixHY,0,0));
    gsl_matrix_set(dispMatrixH,2,5,-gsl_matrix_get(dispMatrixHY,0,1));
    gsl_matrix_set(dispMatrixH,3,4,-gsl_matrix_get(dispMatrixHY,1,0));
    gsl_matrix_set(dispMatrixH,3,5,-gsl_matrix_get(dispMatrixHY,1,1));

    gsl_matrix_set(dispMatrixH,4,0,gsl_matrix_get(H31,0,0));
    gsl_matrix_set(dispMatrixH,5,0,gsl_matrix_get(H31,1,0));
    gsl_matrix_set(dispMatrixH,4,1,gsl_matrix_get(H31,0,1));
    gsl_matrix_set(dispMatrixH,5,1,gsl_matrix_get(H31,1,1));

    gsl_matrix_set(dispMatrixH,4,2,gsl_matrix_get(H32,0,0));
    gsl_matrix_set(dispMatrixH,5,2,gsl_matrix_get(H32,1,0));
    gsl_matrix_set(dispMatrixH,4,3,gsl_matrix_get(H32,0,1));
    gsl_matrix_set(dispMatrixH,5,3,gsl_matrix_get(H32,1,1));


 
    //set Teng Marrix R

    gsl_matrix * tengMatrixR  = gsl_matrix_alloc (6, 6);
    gsl_matrix_set_identity(tengMatrixR);

    gsl_matrix * R2  = gsl_matrix_alloc (2, 2);     //coupling matrix Eq.(6)
    // gsl_matrix_set_identity(R2);                  // full coupling
    gsl_matrix_set_zero(R2);                        // no coupling at the point in SR calculation
    // gsl_matrix_set(R2,0,0,1);
    // gsl_matrix_set(R2,1,1,0.1);

    double b = get_det(R2);
    b = sqrt(1-b);

    gsl_matrix * R12  = gsl_matrix_alloc (2, 2);
    gsl_matrix * R21  = gsl_matrix_alloc (2, 2);


    //	//set R12
    gsl_matrix_transpose(R2);
    gsl_matrix_mul(R2,sympleMarixJ,tempMatrix1);
    gsl_matrix_mul(sympleMarixJ,tempMatrix1,tempMatrix2);
    gsl_matrix_memcpy (R12, tempMatrix2);
    gsl_matrix_transpose(R2);
    //set R21
    gsl_matrix_memcpy (R21, R2);


    gsl_matrix_set(tengMatrixR,0,0,b);
    gsl_matrix_set(tengMatrixR,1,1,b);
    gsl_matrix_set(tengMatrixR,2,2,b);
    gsl_matrix_set(tengMatrixR,3,3,b);

    gsl_matrix_set(tengMatrixR,0,2,gsl_matrix_get(R12,0,0));
    gsl_matrix_set(tengMatrixR,0,3,gsl_matrix_get(R12,0,1));
    gsl_matrix_set(tengMatrixR,1,2,gsl_matrix_get(R12,1,0));
    gsl_matrix_set(tengMatrixR,1,3,gsl_matrix_get(R12,1,1));

    gsl_matrix_set(tengMatrixR,2,0,gsl_matrix_get(R21,0,0));
    gsl_matrix_set(tengMatrixR,2,1,gsl_matrix_get(R21,0,1));
    gsl_matrix_set(tengMatrixR,3,0,gsl_matrix_get(R21,1,0));
    gsl_matrix_set(tengMatrixR,3,1,gsl_matrix_get(R21,1,1));

 
    //set Twiss B Matrix
    gsl_matrix * twissMatrixB = gsl_matrix_alloc (6, 6);
    gsl_matrix_set_zero(twissMatrixB);


    double a00,a01,a10,a11;

    a00 = 1.0/sqrt(latticeInterActionPoint.twissBetaX[k]);
    a10 = latticeInterActionPoint.twissAlphaX[k] / sqrt(latticeInterActionPoint.twissBetaX[k]);
    a11 = sqrt(latticeInterActionPoint.twissBetaX[k]);


    gsl_matrix_set(twissMatrixB,0,0,a00);
    gsl_matrix_set(twissMatrixB,1,0,a10);
    gsl_matrix_set(twissMatrixB,1,1,a11);


    a00 = 1.0/sqrt(latticeInterActionPoint.twissBetaY[k]);
    a10 = latticeInterActionPoint.twissAlphaY[k] / sqrt(latticeInterActionPoint.twissBetaY[k]);
    a11 = sqrt(latticeInterActionPoint.twissBetaY[k]);


    gsl_matrix_set(twissMatrixB,2,2,a00);
    gsl_matrix_set(twissMatrixB,3,2,a10);
    gsl_matrix_set(twissMatrixB,3,3,a11);


    a00 = 1.0/sqrt(latticeInterActionPoint.twissBetaZ[k]);
    a10 = latticeInterActionPoint.twissAlphaZ[k] / sqrt(latticeInterActionPoint.twissBetaZ[k]);
    a11 = sqrt(latticeInterActionPoint.twissBetaZ[k]);


    gsl_matrix_set(twissMatrixB,4,4,a00);
    gsl_matrix_set(twissMatrixB,5,4,a10);
    gsl_matrix_set(twissMatrixB,5,5,a11);



    gsl_matrix * cordTransfer  = gsl_matrix_alloc (6, 6);
    gsl_matrix * invCordTransfer  = gsl_matrix_alloc (6, 6);
    gsl_matrix * cordTransferTemp  = gsl_matrix_alloc (6, 6);
    gsl_matrix_mul(tengMatrixR,dispMatrixH,cordTransferTemp);
    gsl_matrix_mul(twissMatrixB,cordTransferTemp,cordTransfer);


    //get the   invCordTransfer =  (cordTransfer)^-1 ---- approach 1

    gsl_matrix *  sympleMarixJ6 = gsl_matrix_alloc (6, 6);
    gsl_matrix_set_zero(sympleMarixJ6);
    gsl_matrix_set(sympleMarixJ6,0,1, 1);
    gsl_matrix_set(sympleMarixJ6,1,0,-1);
    gsl_matrix_set(sympleMarixJ6,2,3, 1);
    gsl_matrix_set(sympleMarixJ6,3,2,-1);
    gsl_matrix_set(sympleMarixJ6,4,5, 1);
    gsl_matrix_set(sympleMarixJ6,5,4,-1);

    gsl_matrix_transpose(cordTransfer);

    gsl_matrix_mul(cordTransfer,sympleMarixJ6,cordTransferTemp);
    gsl_matrix_mul(sympleMarixJ6,cordTransferTemp,invCordTransfer);
    gsl_matrix_scale (invCordTransfer, -1);

    gsl_matrix_transpose(cordTransfer);
    //-------------------------------------------

    // the same as approache blow, which is much slower
    // get the   invCordTransfer =  (cordTransfer)^-1 
    // gsl_matrix_memcpy(invCordTransfer,cordTransfer);
    // gsl_matrix_inv(invCordTransfer);
    //---------------------------------------------

    double lambda[3];
    double coeff[3];

    for(int i=0;i<3;i++)
    {
        lambda[i] = exp(-1.0/synchRadDampTime[i]);
    }

    // replace here by natural emittance value from input, also set the synchrontron integration values.
    coeff[0] = sqrt(1 - pow(lambda[0],2)) * sqrt(inputParameter.ringParBasic->emitNat[0]);
    coeff[1] = sqrt(1 - pow(lambda[1],2)) * sqrt(inputParameter.ringParBasic->emitNat[1]);
    coeff[2] = sqrt(1 - pow(lambda[2],4)) * sqrt(inputParameter.ringParBasic->emitNat[2]);

    gsl_matrix * vecX   = gsl_matrix_alloc (6, 1);
    gsl_matrix * vecNX  = gsl_matrix_alloc (6, 1);

    double tempX,tempPX,tempY,tempPY,tempZ,tempPZ;
    double randR[6];

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> dx{0,1};

    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        //(0) transfer x to X;
        gsl_matrix_set(vecX,0,0,ePositionX[i]);
        gsl_matrix_set(vecX,1,0,eMomentumX[i]);
        gsl_matrix_set(vecX,2,0,ePositionY[i]);
        gsl_matrix_set(vecX,3,0,eMomentumY[i]);
        gsl_matrix_set(vecX,4,0,ePositionZ[i]);
        gsl_matrix_set(vecX,5,0,eMomentumZ[i]);
        
		gsl_matrix_mul(cordTransfer,vecX,vecNX);

   
        tempX  = gsl_matrix_get(vecNX,0,0);
        tempPX = gsl_matrix_get(vecNX,1,0);
        tempY  = gsl_matrix_get(vecNX,2,0);
        tempPY = gsl_matrix_get(vecNX,3,0);
        tempZ  = gsl_matrix_get(vecNX,4,0);
        tempPZ = gsl_matrix_get(vecNX,5,0);

		//(1)   Eq(8~10) damping and excitation

        //(1.1) synchron-radiation-damping --- only in transverse direction
        tempX  = tempX  * lambda[0] ;
        tempPX = tempPX * lambda[0] ;
        tempY  = tempY  * lambda[1] ;
        tempPY = tempPY * lambda[1] ;
        tempPZ = tempPZ * pow(lambda[2],2) ;


        //(1.2) synchron-radiation-excitiation (transverse only and multi-particle case only)

        if(macroEleNumPerBunch!=1)
        {
            for(int j=0;j<6;j++)
            {
                // randR[j]=Gaussrand(1,0,100.0*j);
                randR[j]=dx(gen);
            }

            tempX  +=  coeff[0] * randR[0];
            tempPX +=  coeff[0] * randR[1];
            tempY  +=  coeff[1] * randR[2];
            tempPY +=  coeff[1] * randR[3];
            tempPZ +=  coeff[2] * randR[5];           
        }


        gsl_matrix_set(vecNX,0,0,tempX);
        gsl_matrix_set(vecNX,1,0,tempPX);
        gsl_matrix_set(vecNX,2,0,tempY);
        gsl_matrix_set(vecNX,3,0,tempPY);
        gsl_matrix_set(vecNX,4,0,tempZ);
        gsl_matrix_set(vecNX,5,0,tempPZ);

		//(2) transfer X to x,  Eq(11)

        gsl_matrix_mul(invCordTransfer,vecNX,vecX);

  
        ePositionX[i] = gsl_matrix_get(vecX,0,0);
        eMomentumX[i] = gsl_matrix_get(vecX,1,0);
        ePositionY[i] = gsl_matrix_get(vecX,2,0);
        eMomentumY[i] = gsl_matrix_get(vecX,3,0);
        ePositionZ[i] = gsl_matrix_get(vecX,4,0);
        eMomentumZ[i] = gsl_matrix_get(vecX,5,0);
    }

    gsl_matrix_free (sympleMarixJ);
    gsl_matrix_free (dispMatrixH);
    gsl_matrix_free (dispMatrixHX);
    gsl_matrix_free (dispMatrixHY);
    gsl_matrix_free (H31);
    gsl_matrix_free (H32);
    gsl_matrix_free (tempMatrix1);
    gsl_matrix_free (tempMatrix2);

    gsl_matrix_free (tengMatrixR);
    gsl_matrix_free (R2);
    gsl_matrix_free (R21);
    gsl_matrix_free (R12);


    gsl_matrix_free (twissMatrixB);

    gsl_matrix_free (cordTransfer);
    gsl_matrix_free (invCordTransfer);
    gsl_matrix_free (cordTransferTemp);
    gsl_matrix_free (sympleMarixJ6);


    gsl_matrix_free (vecX);
    gsl_matrix_free (vecNX);

}

void Bunch::BunchTransferDueToDriveMode(const ReadInputSettings &inputParameter, const int n)
{
    // Ref to Alex Chao Eq.(2.90) in transverse and Eq.(2.86) in longitudinal
    // Can be translated to BBR model as well
    double time = 0.E0;
    int ringHarm        = inputParameter.ringParRf->ringHarm;
    double f0           = inputParameter.ringParBasic->f0;
    double t0           = inputParameter.ringParBasic->t0;
    double rBeta        = inputParameter.ringParBasic->rBeta;
    double tRF          = t0 / ringHarm;
    double driveAmp     = inputParameter.driveMode->driveAmp;
    double driveFre     = inputParameter.driveMode->driveFre;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;

    if(inputParameter.driveMode->drivePlane == 0) // mode is exicted in z direction
    {
        for(int i=0;i<ePositionX.size();i++)
        {
            time = n * t0 + bunchHarmNum * tRF - ePositionZ[i] /  CLight / rBeta;            
            eMomentumZ[i] += driveAmp * cos( 2 * PI * driveFre * time) / electronBeamEnergy / pow(rBeta,2);
        }
    }
    else if (inputParameter.driveMode->drivePlane == 1) // mode is exicted in x direction
    {
        for(int i=0;i<ePositionX.size();i++)
        {
            time = n * t0 + bunchHarmNum * tRF - ePositionZ[i] /  CLight / rBeta;
            eMomentumX[i] += driveAmp * sin( 2 * PI * driveFre * time);
        }
    }
    else if (inputParameter.driveMode->drivePlane == 2) // // mode is exicted in y direction
    {
        for(int i=0;i<ePositionX.size();i++)
        {
            time = n * t0 + bunchHarmNum * tRF - ePositionZ[i] /  CLight / rBeta;
            eMomentumY[i] += driveAmp * sin( 2 * PI * driveFre * time);
        }
    }
    else
    {
        cerr<<"wrong settings in the DRIVEMode->drivePlane, have to be 0,or 1 0r 2 "<<endl;
    }
}


void Bunch::BunchTransferDueToWake()
{
    for(int i=0;i<ePositionX.size();i++)
    {
        eMomentumX[i] += lRWakeForceAver[0];      //rad
        eMomentumY[i] += lRWakeForceAver[1];    
        eMomentumZ[i] += lRWakeForceAver[2];        
    }
}

void Bunch::GetLongiKickDueToCavFB(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double circRing   = inputParameter.ringParBasic->circRing;
    double rfLen      = circRing  / ringHarmH;
    double tRF        = inputParameter.ringParBasic->t0 / ringHarmH; 
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double t0         = inputParameter.ringParBasic->t0;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;

    for(int j=0;j<cavityResonator.resonatorVec.size();j++)
    {
        for(int k=0;k<bunchGap;k++)
        {

            int index = bunchHarmNum + k; 

            cavityResonator.resonatorVec[j].vBSample[index]   = cavityResonator.resonatorVec[j].vBSampleTemp;
            cavityResonator.resonatorVec[j].vGenSample[index] = cavityResonator.resonatorVec[j].resGenVol;
            cavityResonator.resonatorVec[j].vCavSample[index] = cavityResonator.resonatorVec[j].resGenVol + cavityResonator.resonatorVec[j].vBSampleTemp;
            cavityResonator.resonatorVec[j].deltaVCavSample[index] = cavityResonator.resonatorVec[j].vCavSample[index]
                                                                   - cavityResonator.resonatorVec[j].resCavVolReq ;

            //method 1: set this vCavVolDueToDirFB value directly as deltaVCavSample, not robust idea 
            if(cavityResonator.resonatorVec[j].resDirFB!=0)
            {
                cavityResonator.resonatorVec[j].vCavDueToDirFB[index] = cavityResonator.resonatorVec[j].deltaVCavSample[index];
            }

            // method 2: according to vCavDueToDirFB try to get Ig
            // to update the cavity voltage due to the generator currnt genIg, and with updated genIG track the cavity dynamics, 
            // Ref. https://www.aps.anl.gov/files/APS-Uploads/ASDSeminars/2015/2015-05-27_Berenc.pdf
            // slid 44. in his loop, deltaVCavSample is applied to build the control loop
            
            
            
            
            cavityResonator.resonatorVec[j].CavityResonatorDynamics(TrackingTime + k * tRF);
            GetVBSampledDueToIthBunch(j,k,inputParameter,cavityResonator);
        }
    }
   
    TrackingTime += bunchGap * tRF;


    complex<double> cavFB=(0.E0,0.E0);
    double absCavFB,argCavFB;

    // kick due to the cavity feedback, here in below only works with dirFB scheme. It seesm not a good idea for our case.
    for(int j=0;j<cavityResonator.resonatorVec.size();j++)
    {
        if(cavityResonator.resonatorVec[j].resDirFB!=0)
        {
            cavFB    = cavityResonator.resonatorVec[j].vCavDueToDirFB[bunchHarmNum];
            absCavFB = abs(cavFB) * cavityResonator.resonatorVec[j].dirCavFB->gain;
            argCavFB = arg(cavFB) + cavityResonator.resonatorVec[j].dirCavFB->phaseShift 
                                  + cavityResonator.resonatorVec[j].dirCavFB->delay * cavityResonator.resonatorVec[j].resFre * 2 * PI;
            
            cavFB    = absCavFB * exp(li * argCavFB);

            for(int i=0;i<ePositionX.size();i++)
            {
                eMomentumZ[i] -= cavFB.real() /electronBeamEnergy / pow(rBeta,2);
            }
        }
    }
    // end of the kick from the dirfb scheme----------------------------------------------------

}

void Bunch::GetVBSampledDueToIthBunch(const int j,const int k,const ReadInputSettings &inputParameter, CavityResonator &cavityResonator)
{
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double tRF        = inputParameter.ringParBasic->t0 / ringHarmH; 


    double deltaL = tRF / cavityResonator.resonatorVec[j].tF;        
    double cPsi   = 2.0 * PI *  cavityResonator.resonatorVec[j].resFre * tRF;
    cavityResonator.resonatorVec[j].vBSampleTemp *= exp( - deltaL ) * exp (li * cPsi);
    
}



void Bunch::SetBunchPosHistoryDataWithinWindow()
{
    for(int i=0;i<3;i++)
    {
        xyzHistoryDataToFit[i].erase(xyzHistoryDataToFit[i].begin());
    }
    
    xyzHistoryDataToFit[0].push_back(xAver);
    xyzHistoryDataToFit[1].push_back(yAver);
    xyzHistoryDataToFit[2].push_back(zAver);

}


void Bunch::GetBunchHaissinski(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,WakeFunction &sRWakeFunction)
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

void Bunch::GetBunchAverAndBunchLengthFromHaissinskiSolution()
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


vector<double> Bunch::GetProfile(const ReadInputSettings &inputParameter)
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

void Bunch::GetTotHamiltonian()
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

void Bunch::GetRFHamiltonian(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator)
{
    vector<double> vTotRF(haissinski->nz,0);  

    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double t0         = inputParameter.ringParBasic->t0;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double eta        = inputParameter.ringParBasic->eta;
    double u0         = inputParameter.ringParBasic->u0;
    double *alphac    = inputParameter.ringParBasic->alphac;
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

void Bunch::GetWakeHamiltonian(const ReadInputSettings &inputParameter, WakeFunction &sRWakeFunction)
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


void Bunch::GetWakeHamiltonianFromBBR(const ReadInputSettings &inputParameter, WakeFunction &sRWakeFunction)
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



void Bunch::GetWakeHamiltonianFromRW(const ReadInputSettings &inputParameter, WakeFunction &sRWakeFunction)
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



void Bunch::GetParticleLongitudinalPhaseSpace(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,int bunchIndex)
{
    tk::spline wakePotenFit, totHamiltonianFit;
    wakePotenFit.set_points(haissinski->bunchPosZ,haissinski->totWakePoten,tk::spline::cspline);           // fitting the wakePoten 1/[m]
    totHamiltonianFit.set_points(haissinski->bunchPosZ,haissinski->totHamiltonian,tk::spline::cspline);    // fitting the totHamilton 1/[m]


    double zMax = haissinski->averZ + haissinski->rmsZ * 5 + 0.1;
    double zMin = haissinski->averZ - haissinski->rmsZ * 5 - 0.1;
    zMax = + 0.075;
    zMin = - 0.075;
    
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
    

    fout<<"&column name=turns,                              type=float,  &end"<<endl;
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
            fout<<setw(15)<<left<<k
                <<setw(15)<<left<<longiTrajZeta[i][k][0]
                <<setw(15)<<left<<longiTrajZeta[i][k][1]
                <<endl;
        }
    
        
    }
    
    fout.close();        
}

vector<double> Bunch::LeapFrog(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,vector<double> zeta,const tk::spline &wakePotenFit)
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