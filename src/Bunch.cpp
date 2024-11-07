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
#include <gsl/gsl_blas.h>
// #include <cuda_runtime.h>
// #include <cufftXt.h>
// #include <cufft.h>
// #include "CUDAFunction.cuh"
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


using namespace std;
//using std::vector;
//using std::complex;


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

    
    eSurive.resize(macroEleNumPerBunch,0);  
    ePositionX.resize(macroEleNumPerBunch,0E0);
    ePositionY.resize(macroEleNumPerBunch,0E0);
    ePositionZ.resize(macroEleNumPerBunch,0E0);
    eMomentumX.resize(macroEleNumPerBunch,0E0);
    eMomentumY.resize(macroEleNumPerBunch,0E0);
    eMomentumZ.resize(macroEleNumPerBunch,0E0);
    eFxDueToIon.resize(macroEleNumPerBunch,0E0);
    eFyDueToIon.resize(macroEleNumPerBunch,0E0);
    eFzDueToIon.resize(macroEleNumPerBunch,0E0);
    
    lRWakeForceAver.resize(3,0E0);
    accPhaseAdvX = v2d(macroEleNumPerBunch,v1d(3,0.E0));
    accPhaseAdvY = v2d(macroEleNumPerBunch,v1d(3,0.E0));
    accPhaseAdvZ = v2d(macroEleNumPerBunch,v1d(3,0.E0));

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
    int count           = 0;
    int k = 0;
    for(int i=0;i<macroEleNumPerBunch;i++)
    {    
        if( abs(ePositionZ[i]) > t0 * CLight / ringHarm / 2 ) 
        {
            eSurive[i] = 2;   // loss in longitudianl
            count++;
        }
        
        double lossTemp =pow(ePositionX[i]/latticeInterActionPoint.pipeAperatureX[k],2) + pow(ePositionY[i]/latticeInterActionPoint.pipeAperatureY[k],2) ; 
        if(lossTemp >1) 
        {   
            eSurive[i] = 1;                                     // loss in transverse
            count++;
        }
    }
    transmission = (macroEleNumPerBunch - count) / double(macroEleNumPerBunch);
    if(transmission<0.5)
    {
        cout<<"bunch at harmoinc "<<bunchHarmNum<<", more than 50% particles are lost"<<endl;
        exit(0);
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



// void Bunch::BunchTransferDueToLatticeOneTurnT66GPU(const ReadInputSettings &inputParameter, LatticeInterActionPoint &latticeInterActionPoint)
// {
    
//     // prepare the data to munted to GPU 
//     int dimPart6 = macroEleNumPerBunch * 6;   
//     double partCord[dimPart6];
//     for(int i=0;i<macroEleNumPerBunch;i++)
//     {
//         partCord[6*i  ] =  ePositionX[i];
//         partCord[6*i+1] =  eMomentumX[i];
//         partCord[6*i+2] =  ePositionY[i];
//         partCord[6*i+3] =  eMomentumY[i];
//         partCord[6*i+4] =  ePositionZ[i];
//         partCord[6*i+5] =  eMomentumZ[i];        
//     }
  
//     int  paraNum = sizeof(latticeInterActionPoint.latticeParaForOneTurnMap)/sizeof(latticeInterActionPoint.latticeParaForOneTurnMap[0]);
//     GPU_PartiOneTurnTransfer(macroEleNumPerBunch,partCord,paraNum,latticeInterActionPoint.latticeParaForOneTurnMap);
    
//     for(int i=0;i<macroEleNumPerBunch;i++)
//     {
//         ePositionX[i]   = partCord[6*i  ];
//         eMomentumX[i]   = partCord[6*i+1]  ;
//         ePositionY[i]   = partCord[6*i+2]  ;
//         eMomentumY[i]   = partCord[6*i+3]  ;
//         ePositionZ[i]   = partCord[6*i+4]  ;
//         eMomentumZ[i]   = partCord[6*i+5]  ;        
//     }

// }





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
    
    //double *aDTXY = inputParameter.ringParBasic->aDTXY;
    
    double aDTXY[2]; 
    
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
    gsl_matrix *skewQuad  = gsl_matrix_alloc (6, 6);
    gsl_matrix *oneTurnMap  = gsl_matrix_alloc (6, 6);
    gsl_matrix_set_zero(ILmatrix);
    gsl_matrix_set_zero(cord0);
    gsl_matrix_set_zero(cord1);
    gsl_matrix_set_zero(oneTurnMap);
    gsl_matrix_set_identity(skewQuad);
    
    gsl_matrix_set(skewQuad,1,2,inputParameter.ringParBasic->skewQuadK);
    gsl_matrix_set(skewQuad,3,0,inputParameter.ringParBasic->skewQuadK);


    gsl_vector_complex *eval = gsl_vector_complex_alloc (6);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (6, 6);
    gsl_eigen_nonsymmv_workspace * w =  gsl_eigen_nonsymmv_alloc (6);


    //ref. Elegant ILMatrix element  compute_matrices.c line:2442
    //Sy. Lee, Eq.2.257-with dispersion to maintain the periodic condition. ILmatrix[0][5] ILmatrix[1][5] ILmatrix[2][5] ILmatrix[3][5].
    double xPosN,yPosN,xMomN,yMomN;
    double ampX,ampY;
    double nuxtmp,nuytmp,tmp;
    double phix,phiy;

    // generate the transfer matrix for each particle in bunch 
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        if(eSurive[i]!=0) continue;
        // elegant ILMATRIX Eq(56), only keeo the first order here. -- notice the unit of \frac{d eta}/{d delta}
        ePositionX[i]  -=  etax  * eMomentumZ[i];   // [m] 
        ePositionY[i]  -=  etay  * eMomentumZ[i];   // [m]
        eMomentumX[i]  -=  etaxp * eMomentumZ[i];   // [rad]
        eMomentumY[i]  -=  etayp * eMomentumZ[i];   // [rad]

        ampX   = ( pow(ePositionX[i],2) + pow( alphax * ePositionX[i] + betax * eMomentumX[i],2) ) / betax;  //[m]
        ampY   = ( pow(ePositionY[i],2) + pow( alphay * ePositionY[i] + betay * eMomentumY[i],2) ) / betay;  //[m]

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

        tmp = etax  - etax   * cos(phix) - (alphax * etax + betax * etaxp) *  sin(phix);                                   gsl_matrix_set(ILmatrix,0,5,tmp); //R16
        tmp = etay  - etay   * cos(phiy) - (alphay * etay + betay * etayp) *  sin(phiy);                                   gsl_matrix_set(ILmatrix,2,5,tmp); //R36  
        tmp = etaxp - etaxp  * cos(phix) + ( etax  + pow(alphax,2) * etax +  alphax * betax * etaxp) *  sin(phix) / betax; gsl_matrix_set(ILmatrix,1,5,tmp); //R26    
        tmp = etayp - etayp  * cos(phiy) + ( etay  + pow(alphay,2) * etay +  alphay * betay * etayp) *  sin(phiy) / betay; gsl_matrix_set(ILmatrix,3,5,tmp); //R46
        
        tmp = -etaxp + etaxp * cos(phix) + ( etax  + pow(alphax,2) * etax +  alphax * betax * etaxp) *  sin(phix) / betax; gsl_matrix_set(ILmatrix,4,0,tmp); //R51
        tmp = -etayp + etayp * cos(phiy) + ( etay  + pow(alphay,2) * etay +  alphay * betay * etayp) *  sin(phiy) / betay; gsl_matrix_set(ILmatrix,4,2,tmp); //R53 
        tmp =  etax  - etax  * cos(phix) + ( alphax * etax +  betax * etaxp) * sin(phix);                                  gsl_matrix_set(ILmatrix,4,1,tmp); //R52 
        tmp =  etay  - etay  * cos(phiy) + ( alphay * etay +  betay * etayp) * sin(phiy);                                  gsl_matrix_set(ILmatrix,4,3,tmp); //R54  
        // R55, R56, R65, R66 and R16 R26 R36 R46 is set to ensure the period one turn solution
        gsl_matrix_set(ILmatrix,4,4,1); //R55
        gsl_matrix_set(ILmatrix,5,5,1); //R66

        // set the right particle coordinate (x_{beta} ) for one turn transfer. 
        // also make the R56 kick extra due to high order momentum compact factor, ampllidude depended length...    
        // gsl matrix multi give almost the same perfromance in simulation speed
        gsl_matrix_set(cord0,0,0,ePositionX[i]);
        gsl_matrix_set(cord0,1,0,eMomentumX[i]);
        gsl_matrix_set(cord0,2,0,ePositionY[i]);
        gsl_matrix_set(cord0,3,0,eMomentumY[i]);
        gsl_matrix_set(cord0,4,0,ePositionZ[i]);
        gsl_matrix_set(cord0,5,0,eMomentumZ[i]);
       
        
        // gsl_matrix_mul(ILmatrix,cord0,cord1);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,ILmatrix,cord0,0.0,cord1);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,skewQuad,cord1,0.0,cord0);

        ePositionX[i]  = gsl_matrix_get(cord0,0,0);
        eMomentumX[i]  = gsl_matrix_get(cord0,1,0);
        ePositionY[i]  = gsl_matrix_get(cord0,2,0);
        eMomentumY[i]  = gsl_matrix_get(cord0,3,0);
        ePositionZ[i]  = gsl_matrix_get(cord0,4,0);
        eMomentumZ[i]  = gsl_matrix_get(cord0,5,0);


        // here try to get the eigen value of the one turn transfer map
        // ----start -------------------------------  
        // gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,ILmatrix,skewQuad,0.0,oneTurnMap);
        // gsl_eigen_nonsymmv (oneTurnMap, eval, evec, w);
        // gsl_eigen_nonsymmv_sort (eval, evec,GSL_EIGEN_SORT_ABS_DESC);
        
        // gsl_complex gsltemp;
        // GSL_REAL(gsltemp)=1;
        // GSL_IMAG(gsltemp)=0;
        // for(int i=0;i<6;i++)
        // {
        //     gsl_complex eval_i = gsl_vector_complex_get (eval, i);
        //     gsl_vector_complex_view evec_i  = gsl_matrix_complex_column (evec, i);
        //     gsltemp = gsl_complex_mul(eval_i,gsltemp);
        //     printf ("eigenvalue = %g + %gi   %g     %g \n", GSL_REAL(eval_i), GSL_IMAG(eval_i), gsl_complex_abs(eval_i),gsltemp );

        // }
        // ----end -------------------------------

        // for(int j=0;j<6;j++)
        // {
        //     for(int k=0;k<6;k++) printf("%15.7f\t",ILmatrix->data[j * ILmatrix->tda+k]);
        //     cout<<endl;       
        // }
        
        // for(int j=0;j<6;j++) printf("%15.10f\t", cord0->data[j * cord0->tda]);
        // cout<<endl;
        // for(int j=0;j<6;j++) printf("%15.10f\t", cord1->data[j * cord1->tda]);
        // cout<<endl;
                              
    }

    gsl_matrix_free(ILmatrix);
    gsl_matrix_free (cord0);
    gsl_matrix_free (cord1);
    gsl_matrix_free(skewQuad);
    gsl_matrix_free(oneTurnMap);

    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    gsl_eigen_nonsymmv_free (w);

}


void Bunch::GetLongiKickDueToCavFB(const ReadInputSettings &inputParameter,Resonator &resonator)
{
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
   
    complex<double> cavFB=(0.E0,0.E0);
    double absCavFB,argCavFB;

    // kick due to the cavity dirFB. -- not physical one        
    cavFB    = resonator.vCavDueToDirFB[bunchHarmNum];
    absCavFB = abs(cavFB) * resonator.dirCavFB->gain;
    argCavFB = arg(cavFB) + resonator.dirCavFB->phaseShift + resonator.dirCavFB->delay * resonator.resFre * 2 * PI;
    
    cavFB    = absCavFB * exp(li * argCavFB);

    for(int i=0;i<ePositionX.size();i++)
    {
        eMomentumZ[i] += cavFB.real() /electronBeamEnergy / pow(rBeta,2);
    }
}


void Bunch::BunchMomentumUpdateDueToRFCA(const ReadInputSettings &inputParameter,Resonator &resonator,int resIndex)
{
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double f0         = inputParameter.ringParBasic->f0;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    double fRF         = f0 * ringHarmH;
    
    int resHarm = resonator.resHarm;
    double resFre  = resonator.resFre;

    complex<double> cavVoltage =(0.E0,0.E0);
    complex<double> genVoltage =(0.E0,0.E0);

    for (int i=0;i<macroEleNumPerBunch;i++)
    {
        genVoltage = resonator.resCavVolReq * exp( - li * ePositionZ[i] / CLight / rBeta * 2. * PI * double(resHarm) * fRF );
        cavVoltage = genVoltage;
        eMomentumZ[i] += cavVoltage.real()  /electronBeamEnergy / pow(rBeta,2);
    }

    // the info cavity, generator and beam induced voltage bunch feels
    bunchRFModeInfo->induceVolBunchCen[resIndex]   = complex<double>(0.E0,0.E0);
    bunchRFModeInfo->genVolBunchAver[resIndex]     = genVoltage;               
    bunchRFModeInfo->cavVolBunchCen[resIndex]      = cavVoltage;
    bunchRFModeInfo->selfLossVolBunchCen[resIndex] = complex<double>(0.E0,0.E0);


}

void Bunch::BunchEnergyLossOneTurn(const ReadInputSettings &inputParameter)
{
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double u0         = inputParameter.ringParBasic->u0;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;

    for (int i=0;i<macroEleNumPerBunch;i++)
        eMomentumZ[i] -= u0 / electronBeamEnergy / pow(rBeta,2);   
}


void Bunch::BunchLongPosTransferOneTurn(const ReadInputSettings &inputParameter)
{
    double circRing = inputParameter.ringParBasic->circRing;
    double *alphac = inputParameter.ringParBasic->alphac;
    
    for(int i=0;i<macroEleNumPerBunch;i++)
        ePositionZ[i] -= circRing * (alphac[0] * eMomentumZ[i]  + alphac[1] * pow( eMomentumZ[i] ,2) + alphac[2] * pow( eMomentumZ[i], 3) ) ;
}


void Bunch::InitialAccumPhaseAdV(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter)
{
    int k = 0;
	gsl_matrix *vecX     = gsl_matrix_alloc (6, 1);
	gsl_matrix *vecX1    = gsl_matrix_alloc (6, 1);
    gsl_matrix *averVecX   = gsl_matrix_alloc (6, 1);
	gsl_matrix *avervecX1  = gsl_matrix_alloc (6, 1);   
	gsl_matrix *B1H1 = latticeInterActionPoint.symplecticMapB1H1[k].mat2D;

    averVecX->data[0 * averVecX->tda] = xAver;
    averVecX->data[1 * averVecX->tda] = pxAver;
    averVecX->data[2 * averVecX->tda] = yAver;
    averVecX->data[3 * averVecX->tda] = pyAver;
    averVecX->data[4 * averVecX->tda] = zAver;
    averVecX->data[5 * averVecX->tda] = pzAver;

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,B1H1,averVecX,0.0,avervecX1);

    double phaseAdvXMean = atan2(gsl_matrix_get(averVecX,1,0),gsl_matrix_get(averVecX,0,0));
    double phaseAdvYMean = atan2(gsl_matrix_get(averVecX,3,0),gsl_matrix_get(averVecX,2,0)); 

    double nux = inputParameter.ringParBasic->workQx;
    double nuy = inputParameter.ringParBasic->workQy;
    nux -= floor(nux);
    nuy -= floor(nuy);

    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        vecX->data[0 * vecX->tda] = ePositionX[i];
        vecX->data[1 * vecX->tda] = eMomentumX[i];
        vecX->data[2 * vecX->tda] = ePositionY[i];
        vecX->data[3 * vecX->tda] = eMomentumY[i];
        vecX->data[4 * vecX->tda] = ePositionZ[i]; 
        vecX->data[5 * vecX->tda] = eMomentumZ[i]; 
    	
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,B1H1,vecX,0.0,vecX1);

        accPhaseAdvX[i][0] = atan2(gsl_matrix_get(vecX1,1,0),gsl_matrix_get(vecX1,0,0));
        accPhaseAdvY[i][0] = atan2(gsl_matrix_get(vecX1,3,0),gsl_matrix_get(vecX1,2,0));
        accPhaseAdvZ[i][0] = atan2(gsl_matrix_get(vecX1,5,0),gsl_matrix_get(vecX1,4,0));
        // accPhaseAdvX[i][0] = atan2(gsl_matrix_get(vecX1,1,0),gsl_matrix_get(vecX1,0,0)) - phaseAdvXMean;
        // accPhaseAdvY[i][0] = atan2(gsl_matrix_get(vecX1,1,0),gsl_matrix_get(vecX1,0,0)) - phaseAdvYMean;
    }

	gsl_matrix_free (vecX);      
	gsl_matrix_free (vecX1); 
    gsl_matrix_free (averVecX);
	gsl_matrix_free (avervecX1);
    B1H1 = NULL;   
}

void Bunch::GetAccumuPhaseAdv(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter)
{
    // not used 
    int k = 0;
	gsl_matrix *vecX       = gsl_matrix_alloc (6, 1);
	gsl_matrix *vecX1      = gsl_matrix_alloc (6, 1); 
    gsl_matrix *averVecX   = gsl_matrix_alloc (6, 1);
	gsl_matrix *avervecX1  = gsl_matrix_alloc (6, 1); 
	gsl_matrix *B1H1 = latticeInterActionPoint.symplecticMapB1H1[k].mat2D;

    averVecX->data[0 * averVecX->tda] = xAver;
    averVecX->data[1 * averVecX->tda] = pxAver;
    averVecX->data[2 * averVecX->tda] = yAver;
    averVecX->data[3 * averVecX->tda] = pyAver;
    averVecX->data[4 * averVecX->tda] = zAver;
    averVecX->data[5 * averVecX->tda] = pzAver;

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,B1H1,averVecX,0.0,avervecX1);

    double phaseAdvXMean = atan2(gsl_matrix_get(averVecX,1,0),gsl_matrix_get(averVecX,0,0));
    double phaseAdvYMean = atan2(gsl_matrix_get(averVecX,3,0),gsl_matrix_get(averVecX,2,0)); 

    double nux = inputParameter.ringParBasic->workQx;
    double nuy = inputParameter.ringParBasic->workQy;
    nux -= floor(nux);
    nuy -= floor(nuy);
    double tmpx,tmpy;

    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        vecX->data[0 * vecX->tda] = ePositionX[i] - xAver;
        vecX->data[1 * vecX->tda] = eMomentumX[i] - pxAver;
        vecX->data[2 * vecX->tda] = ePositionY[i] - yAver;
        vecX->data[3 * vecX->tda] = eMomentumY[i] - pyAver;
        vecX->data[4 * vecX->tda] = ePositionZ[i] - zAver; 
        vecX->data[5 * vecX->tda] = eMomentumZ[i] - pzAver; 

        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,B1H1,vecX,0.0,vecX1);
        // here to get the phase advances reference to the bunch center 
        accPhaseAdvX[i][1] = atan2(gsl_matrix_get(vecX1,1,0),gsl_matrix_get(vecX1,0,0)) - phaseAdvXMean;
        accPhaseAdvY[i][1] = atan2(gsl_matrix_get(vecX1,3,0),gsl_matrix_get(vecX1,2,0)) - phaseAdvYMean;

        tmpx = accPhaseAdvX[i][0] - accPhaseAdvX[i][1];
        tmpy = accPhaseAdvY[i][0] - accPhaseAdvY[i][1]; 
        accPhaseAdvX[i][2] += tmpx >= 0? tmpx : tmpx + 2 * PI;
        accPhaseAdvY[i][2] += tmpy >= 0? tmpy : tmpy + 2 * PI;

        cout<<setw(15)<<i;
        for(int j=0;j<3;++j)
        {
            cout<<setw(15)<<accPhaseAdvX[i][j];
        }
        cout<<endl;

        accPhaseAdvX[i][0] = accPhaseAdvX[i][1];
        accPhaseAdvY[i][0] = accPhaseAdvY[i][1];

    }
    
	gsl_matrix_free (vecX);
	gsl_matrix_free (vecX1);
    gsl_matrix_free (averVecX);
	gsl_matrix_free (avervecX1);
    B1H1 = NULL;  
    cout<<__LINE__<<__FILE__<<endl;
    getchar();
}



void Bunch::BunchTransferDueToLatticeTSymplectic(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint, int k)
{
	latticeSetionPassedCount = currentTurnNum * inputParameter.ringParBasic->ringSectNum + k;
    double etax   = latticeInterActionPoint.twissDispX[k];
    double etaxp  = latticeInterActionPoint.twissDispPX[k];  // \frac{disP}{ds} 
    double etay   = latticeInterActionPoint.twissDispY[k];
    double etayp  = latticeInterActionPoint.twissDispPY[k];  // \frac{disP}{ds} 
	double alphax = latticeInterActionPoint.twissAlphaX[k];
    double alphay = latticeInterActionPoint.twissAlphaY[k];
    double betax  = latticeInterActionPoint.twissBetaX[k];
    double betay  = latticeInterActionPoint.twissBetaY[k];
    
    double gammax = (1 + pow(alphax,2))/betax;
    double gammay = (1 + pow(alphay,2))/betay;
	
    double *alphac = inputParameter.ringParBasic->alphac;
    // dnux/dAx dnux/dAy dnux^2/dAx^2 dnux^2/dAy^2 dnux^2/dAxdAy 
    // dnuy/dAx dnuy/dAy dnuy^2/dAx^2 dnuy^2/dAy^2 dnuy^2/dAxdAy
    
    double *aDTX  = inputParameter.ringParBasic->aDTX;
    double *aDTY  = inputParameter.ringParBasic->aDTY;
    
    double circRing = inputParameter.ringParBasic->circRing;
    double nux    = inputParameter.ringParBasic->workQx;
    double nuy 	  = inputParameter.ringParBasic->workQy;
    double chromx = inputParameter.ringParBasic->chrom[0];
    double chromy = inputParameter.ringParBasic->chrom[1];
    double t0     = inputParameter.ringParBasic->t0;
    double eta    = inputParameter.ringParBasic->eta;

	double phaseAdvX = latticeInterActionPoint.phaseAdvX12[k];
	double phaseAdvY = latticeInterActionPoint.phaseAdvY12[k];

	double weighX = phaseAdvX / (2 * PI * nux);
	double weighY = phaseAdvY / (2 * PI * nuy);
	
	double x,px,y,py,ampX,ampY,phiX,phiY;
	
	gsl_matrix *matRotat = gsl_matrix_alloc (6, 6);
	gsl_matrix *vecX     = gsl_matrix_alloc (6, 1);
	gsl_matrix *vecX1    = gsl_matrix_alloc (6, 1);  
	
	gsl_matrix *B1H1 = latticeInterActionPoint.symplecticMapB1H1[k].mat2D;
	gsl_matrix *H2B2 = latticeInterActionPoint.symplecticMapInvH2InvB2[k].mat2D;

    double tmpx,tmpy,tmpz;
    double vec0[2],vec1[2];

	for(int i=0;i<macroEleNumPerBunch;i++)
    {
        if(eSurive[i]!=0) continue;
        // elegant ILMATRIX Eq(56), only keeo the first order here. -- notice the unit of \frac{d eta}/{d delta}
   	
		x  = ePositionX[i] - etax  * eMomentumZ[i]; // [m] 
		px = eMomentumX[i] - etaxp * eMomentumZ[i]; // [rad]
		y  = ePositionY[i] - etay  * eMomentumZ[i]; // [m] 
		py = eMomentumY[i] - etayp * eMomentumZ[i]; // [rad]

        ampX   = ( pow(x,2) + pow( alphax * x + betax * px,2) ) / betax;  //[m]
        ampY   = ( pow(y,2) + pow( alphay * y + betay * py,2) ) / betay;  //[m]

        phiX = phaseAdvX + (chromx * eMomentumZ[i] + aDTX[0] * ampX + aDTX[1] * ampY + aDTX[2] * pow(ampX,2)/2 + aDTX[3] * pow(ampY,2)/2 + aDTX[4] * ampX * ampY ) * weighX;
        phiY = phaseAdvY + (chromy * eMomentumZ[i] + aDTY[0] * ampX + aDTY[1] * ampY + aDTY[2] * pow(ampX,2)/2 + aDTY[3] * pow(ampY,2)/2 + aDTY[4] * ampX * ampY ) * weighY;
    	
    	gsl_matrix_set_zero(matRotat);
    	
    	gsl_matrix_set(matRotat,0,0, cos(phiX)); gsl_matrix_set(matRotat,0,1, sin(phiX)); gsl_matrix_set(matRotat,1,0, -sin(phiX)); gsl_matrix_set(matRotat,1,1, cos(phiX));
    	gsl_matrix_set(matRotat,2,2, cos(phiY)); gsl_matrix_set(matRotat,2,3, sin(phiY)); gsl_matrix_set(matRotat,3,2, -sin(phiY)); gsl_matrix_set(matRotat,3,3, cos(phiY));
    	gsl_matrix_set(matRotat,4,4, 1		  );																					gsl_matrix_set(matRotat,5,5, 		1 );
    	
        vecX->data[0 * vecX->tda] = ePositionX[i];
        vecX->data[1 * vecX->tda] = eMomentumX[i];
        vecX->data[2 * vecX->tda] = ePositionY[i];
        vecX->data[3 * vecX->tda] = eMomentumY[i];
        vecX->data[4 * vecX->tda] = ePositionZ[i]; 
        vecX->data[5 * vecX->tda] = eMomentumZ[i]; 
    	
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,B1H1,vecX,0.0,vecX1);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,matRotat,vecX1,0.0,vecX);

        // here to get the phase advances based on turn-by-trun data 
        if(k==inputParameter.ringParBasic->ringSectNum-1)
        {
            accPhaseAdvX[i][1] = atan2(gsl_matrix_get(vecX,1,0),gsl_matrix_get(vecX,0,0));
            accPhaseAdvY[i][1] = atan2(gsl_matrix_get(vecX,3,0),gsl_matrix_get(vecX,2,0));
            accPhaseAdvZ[i][1] = atan2(gsl_matrix_get(vecX,5,0),gsl_matrix_get(vecX,4,0));
            
            tmpx = accPhaseAdvX[i][0] - accPhaseAdvX[i][1];
            tmpy = accPhaseAdvY[i][0] - accPhaseAdvY[i][1];
            tmpz = accPhaseAdvZ[i][0] - accPhaseAdvZ[i][1];
            
            accPhaseAdvX[i][2] += tmpx >= 0? tmpx : tmpx + 2 * PI;
            accPhaseAdvY[i][2] += tmpy >= 0? tmpy : tmpy + 2 * PI;
            accPhaseAdvZ[i][2] += tmpz >= 0? tmpz : tmpz + 2 * PI;

            accPhaseAdvX[i][0] = accPhaseAdvX[i][1];
            accPhaseAdvY[i][0] = accPhaseAdvY[i][1];
            accPhaseAdvZ[i][0] = accPhaseAdvZ[i][1];

        }

		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,H2B2,vecX,0.0,vecX1);
		
		ePositionX[i] =  gsl_matrix_get(vecX1,0,0);
		eMomentumX[i] =  gsl_matrix_get(vecX1,1,0);
		ePositionY[i] =  gsl_matrix_get(vecX1,2,0);
		eMomentumY[i] =  gsl_matrix_get(vecX1,3,0);
		ePositionZ[i] =  gsl_matrix_get(vecX1,4,0);
		eMomentumZ[i] =  gsl_matrix_get(vecX1,5,0);
	}
	
	gsl_matrix_free (vecX);
	gsl_matrix_free (vecX1);
	B1H1 = NULL;
	H2B2 = NULL;
    matRotat = NULL;	

    // getchar();
    
}

void Bunch::BunchTransferDuetoSkewQuad(const ReadInputSettings &inputParameter)
{
	gsl_matrix *skewQuad  = gsl_matrix_alloc (6, 6);
	gsl_matrix_set_identity(skewQuad);
    gsl_matrix_set(skewQuad,1,2,inputParameter.ringParBasic->skewQuadK);
    gsl_matrix_set(skewQuad,3,0,inputParameter.ringParBasic->skewQuadK);
   
    gsl_matrix *vecX     = gsl_matrix_alloc (6, 1);
    gsl_matrix *vecX1    = gsl_matrix_alloc (6, 1);
    gsl_matrix_set_zero(vecX);
    gsl_matrix_set_zero(vecX1);
    
	for(int i=0;i<macroEleNumPerBunch;i++)
    {
        if(eSurive[i]!=0) continue;
		
        vecX->data[0 * vecX->tda] = ePositionX[i];
        vecX->data[1 * vecX->tda] = eMomentumX[i];
        vecX->data[2 * vecX->tda] = ePositionY[i];
        vecX->data[3 * vecX->tda] = eMomentumY[i];
        vecX->data[4 * vecX->tda] = ePositionZ[i]; 
        vecX->data[5 * vecX->tda] = eMomentumZ[i]; 
		
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,skewQuad,vecX,0.0,vecX1);	
		
		ePositionX[i] =  gsl_matrix_get(vecX1,0,0);
		eMomentumX[i] =  gsl_matrix_get(vecX1,1,0);
		ePositionY[i] =  gsl_matrix_get(vecX1,2,0);
		eMomentumY[i] =  gsl_matrix_get(vecX1,3,0);
		ePositionZ[i] =  gsl_matrix_get(vecX1,4,0);
		eMomentumZ[i] =  gsl_matrix_get(vecX1,5,0);
	}
	
	
	gsl_matrix_free(skewQuad);
	gsl_matrix_free(vecX);
	gsl_matrix_free(vecX1);

}




void Bunch::BunchTransferDueToLatticeT(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint, int k)
{

    double circRing = inputParameter.ringParBasic->circRing;
    double *alphac = inputParameter.ringParBasic->alphac;
    double nux = inputParameter.ringParBasic->workQx;
    double nuy = inputParameter.ringParBasic->workQy;
    
    double betaX1, betaX2, alphaX1, alphaX2, phiX0, etaX, etaXP, chromX;
    double betaY1, betaY2, alphaY1, alphaY2, phiY0, etaY, etaYP, chromY;
    
    etaX  = latticeInterActionPoint.twissDispX[k];
    etaY  = latticeInterActionPoint.twissDispY[k];
    etaXP = latticeInterActionPoint.twissDispPX[k];
    etaYP = latticeInterActionPoint.twissDispPY[k];

    alphaX1 = latticeInterActionPoint.twissAlphaX[k];
    betaX1  = latticeInterActionPoint.twissBetaX[k];
    alphaY1 = latticeInterActionPoint.twissAlphaY[k];
    betaY1  = latticeInterActionPoint.twissBetaY[k];

    chromX = inputParameter.ringParBasic->chrom[0] * latticeInterActionPoint.interactionLength[k] / circRing;
    chromY = inputParameter.ringParBasic->chrom[1] * latticeInterActionPoint.interactionLength[k] / circRing; 


    if (k==latticeInterActionPoint.numberOfInteraction-1)
    {
        alphaX2 = latticeInterActionPoint.twissAlphaX[0];
        betaX2  = latticeInterActionPoint.twissBetaX[0];
        alphaY2 = latticeInterActionPoint.twissAlphaY[0];
        betaY2  = latticeInterActionPoint.twissBetaY[0]; 
        phiX0   = 2 * PI * nux  -   latticeInterActionPoint.xPhaseAdv[k]; 
        phiY0   = 2 * PI * nuy  -   latticeInterActionPoint.yPhaseAdv[k];      
    }
    else
    {
        alphaX2 = latticeInterActionPoint.twissAlphaX[k+1];
        betaX2  = latticeInterActionPoint.twissBetaX[k+1];
        alphaY2 = latticeInterActionPoint.twissAlphaY[k+1];
        betaY2  = latticeInterActionPoint.twissBetaY[k+1];
        phiX0   = latticeInterActionPoint.xPhaseAdv[k+1]  - latticeInterActionPoint.xPhaseAdv[k]; 
        phiY0   = latticeInterActionPoint.yPhaseAdv[k+1]  - latticeInterActionPoint.yPhaseAdv[k];   
    }

    gsl_matrix * vecX   = gsl_matrix_alloc (3, 1);
    gsl_matrix * vecNX  = gsl_matrix_alloc (3, 1);
    gsl_matrix * cordTransfer     = gsl_matrix_alloc (3, 3);
    
    double tmp, phiX,phiY;

    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        if(eSurive[i]!=0) continue;

        phiX = phiX0 + chromX * eMomentumZ[i];
        // transformation in the x direction:
        vecX->data[0 * vecX->tda] = ePositionX[i];
        vecX->data[1 * vecX->tda] = eMomentumX[i];
        vecX->data[2 * vecX->tda] = eMomentumZ[i]; 

        gsl_matrix_set_zero(cordTransfer);
        tmp = sqrt(betaX2 / betaX1) * (cos(phiX) + alphaX1 * sin(phiX));                                         gsl_matrix_set(cordTransfer,0,0,tmp);
        tmp = sqrt(betaX2 * betaX1) * (                      sin(phiX));                                         gsl_matrix_set(cordTransfer,0,1,tmp);
        tmp = etaX;                                                                                              gsl_matrix_set(cordTransfer,0,2,tmp);

        tmp = ((-1 - alphaX1 * alphaX2) * sin(phiX) + (alphaX1 - alphaX2) * cos(phiX)) / sqrt(betaX2 * betaX1);  gsl_matrix_set(cordTransfer,1,0,tmp);
        tmp = sqrt(betaX1 / betaX2) * (cos(phiX) - alphaX2 * sin(phiX));                                         gsl_matrix_set(cordTransfer,1,1,tmp);  
        tmp = etaXP ;                                                                                            gsl_matrix_set(cordTransfer,1,2,tmp);
        gsl_matrix_set(cordTransfer,2,2,1);

 
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,cordTransfer,vecX,0.0,vecNX);

        ePositionX[i] = gsl_matrix_get(vecNX,0,0);
        eMomentumX[i] = gsl_matrix_get(vecNX,1,0);
        eMomentumZ[i] = gsl_matrix_get(vecNX,2,0);

        // transformation in the y direction:
        phiY = phiY0 + chromY * eMomentumZ[i];
        vecX->data[0 * vecX->tda] = ePositionY[i];
        vecX->data[1 * vecX->tda] = eMomentumY[i];
        vecX->data[2 * vecX->tda] = eMomentumZ[i];

        gsl_matrix_set_zero(cordTransfer);
        tmp = sqrt(betaY2 / betaY1) * (cos(phiY) + alphaY1 * sin(phiY));                                         gsl_matrix_set(cordTransfer,0,0,tmp);
        tmp = sqrt(betaY2 * betaY1) * (                      sin(phiY));                                         gsl_matrix_set(cordTransfer,0,1,tmp);
        tmp = etaY;                                                                                              gsl_matrix_set(cordTransfer,0,2,tmp);

        tmp = ((-1 - alphaY1 * alphaY2) * sin(phiY) + (alphaY1 - alphaY2) * cos(phiY)) / sqrt(betaY2 * betaY1);  gsl_matrix_set(cordTransfer,1,0,tmp);
        tmp = sqrt(betaY1 / betaY2) * (cos(phiY) - alphaY2 * sin(phiY));                                         gsl_matrix_set(cordTransfer,1,1,tmp);  
        tmp = etaYP ;                                                                                            gsl_matrix_set(cordTransfer,1,2,tmp);

        gsl_matrix_set(cordTransfer,2,2,1);

        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,cordTransfer,vecX,0.0,vecNX);

        ePositionY[i] = gsl_matrix_get(vecNX,0,0);
        eMomentumY[i] = gsl_matrix_get(vecNX,1,0);
        eMomentumZ[i] = gsl_matrix_get(vecNX,2,0);

        double lossTemp =pow(ePositionX[i]/latticeInterActionPoint.pipeAperatureX[k],2) + pow(ePositionY[i]/latticeInterActionPoint.pipeAperatureY[k],2) ; 
        if(lossTemp >1) eSurive[i] = 1;   // loss in transverse
    }   

    gsl_matrix_free (cordTransfer);
    gsl_matrix_free (vecX);
    gsl_matrix_free (vecNX);


    // very initial version of particle tansfer in transverse
    // for(int i=0;i<macroEleNumPerBunch;i++)
    // {
    //     if(eSurive[i]!=0) continue;
    //     // refers to SY. Lee Eq. 2.67

    //     xtemp  = latticeInterActionPoint.xTransferMatrix[k][0] * ePositionX[i]
    //            + latticeInterActionPoint.xTransferMatrix[k][1] * eMomentumX[i];

    //     xPtemp = latticeInterActionPoint.xTransferMatrix[k][2] * ePositionX[i]
    //            + latticeInterActionPoint.xTransferMatrix[k][3] * eMomentumX[i];

    //     ytemp  = latticeInterActionPoint.yTransferMatrix[k][0] * ePositionY[i]
    //            + latticeInterActionPoint.yTransferMatrix[k][1] * eMomentumY[i];

    //     yPtemp = latticeInterActionPoint.yTransferMatrix[k][2] * ePositionY[i]
    //            + latticeInterActionPoint.yTransferMatrix[k][3] * eMomentumY[i];

    //     ePositionX[i] = xtemp;
    //     ePositionY[i] = ytemp;
    //     eMomentumX[i] = xPtemp;
    //     eMomentumY[i] = yPtemp;

    //     double lossTemp =pow(ePositionX[i]/latticeInterActionPoint.pipeAperatureX[k],2) + pow(ePositionY[i]/latticeInterActionPoint.pipeAperatureY[k],2) ; 
    //     if(lossTemp >1) eSurive[i] = 1;   // loss in transverse
    // }
   
}



void Bunch::BunchSynRadDamping(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint)
{
    //Note: the SynRadDamping and excitation is follow Yuan ZHang's PRAB paper. 

    // in the unit of number of truns for synchRadDampTime setting.
    vector<double> synchRadDampTime;
    synchRadDampTime.resize(3);
    for(int i=0;i<synchRadDampTime.size();i++)
    {
        synchRadDampTime[i] = inputParameter.ringParBasic->synchRadDampTime[i];   
    }
	
	gsl_matrix *B1H1 = latticeInterActionPoint.symplecticMapB1H1[0].mat2D;
	gsl_matrix *H1B1 = latticeInterActionPoint.symplecticMapInvH1InvB1[0].mat2D;

	/*
	gsl_matrix * B1H1  = gsl_matrix_alloc (6, 6);
    gsl_matrix * H2B2  = gsl_matrix_alloc (6, 6);
    for(int i=0;i<6;i++)
    {
        for(int j=0;j<6;j++)
        {
            gsl_matrix_set(cordTransfer,   i,j,latticeInterActionPoint.latticeSynRadBRH[6*i+j]);
            gsl_matrix_set(invCordTransfer,i,j,latticeInterActionPoint.latticeSynRadBRH[6*i+j+36]);
        }
    }
	*/
		
    double lambda[3];
    double coeff[3];

    for(int i=0;i<3;i++)
    {
        lambda[i] = exp(-1.0/synchRadDampTime[i]);
    }

    // replace here by natural emittance value from input, also set the synchrontron integration values.
    // which is actually are the eigen emittance from simulation
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
        if(eSurive[i]!=0) continue;
        
        vecX->data[0 * vecX->tda] = ePositionX[i];
        vecX->data[1 * vecX->tda] = eMomentumX[i];
        vecX->data[2 * vecX->tda] = ePositionY[i];
        vecX->data[3 * vecX->tda] = eMomentumY[i];
        vecX->data[4 * vecX->tda] = ePositionZ[i];
        vecX->data[5 * vecX->tda] = eMomentumZ[i];      
		
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,B1H1,vecX,0.0,vecNX);   
              
        vecNX->data[0 * vecNX->tda] *= lambda[0];
        vecNX->data[1 * vecNX->tda] *= lambda[0];
        vecNX->data[2 * vecNX->tda] *= lambda[1];
        vecNX->data[3 * vecNX->tda] *= lambda[1];
        vecNX->data[5 * vecNX->tda] *= lambda[2] * lambda[2];
    
        if(macroEleNumPerBunch!=1)
        {
            for(int j=0;j<6;j++)
            {
                randR[j]=dx(gen);
            }

            vecNX->data[0 * vecNX->tda] +=  coeff[0] * randR[0];
            vecNX->data[1 * vecNX->tda] +=  coeff[0] * randR[1];
            vecNX->data[2 * vecNX->tda] +=  coeff[1] * randR[2];
            vecNX->data[3 * vecNX->tda] +=  coeff[1] * randR[3];
            vecNX->data[5 * vecNX->tda] +=  coeff[2] * randR[5];           
        }

		//(2) transfer X to x,  Eq(11)
        
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,H1B1,vecNX,0.0,vecX); 
      
        ePositionX[i] = gsl_matrix_get(vecX,0,0);
        eMomentumX[i] = gsl_matrix_get(vecX,1,0);
        ePositionY[i] = gsl_matrix_get(vecX,2,0);
        eMomentumY[i] = gsl_matrix_get(vecX,3,0);
        ePositionZ[i] = gsl_matrix_get(vecX,4,0);
        eMomentumZ[i] = gsl_matrix_get(vecX,5,0);
    }

    H1B1 = NULL;
    B1H1 = NULL;
        
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

        // cout<<lRWakeForceAver[0]<<" "<<lRWakeForceAver[1]<<"    "<<lRWakeForceAver[2]<<endl;
    }

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

    int nz                = floor(sigmaT0 * 20 / haissinski->dt);      // +- 20 rms bunch size as the range for haissinksi binsize=1 ps here used as default setting
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
    haissinski->vRF.resize(haissinski->nz);
    haissinski->nus.resize(haissinski->nz);
    haissinski->actionJ.resize(haissinski->nz); 

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

    double coef = 2 * PI * ringHarmH * f0 * eta * pow(sdelta0,2);  //  [1/s]   // PRAB 17.064401 Eq.(8) 


    vector<double> profile(haissinski->nz,0);
    double norm=0;
    for(int i=0;i<haissinski->nz;i++)
    {
        profile[i] = exp(-haissinski->totHamiltonian[i]/coef);
        norm += profile[i] * haissinski->dz;
    }
    
    int electronNumPerBunchTemp;
    electronNumPerBunchTemp = electronNumPerBunch==0? 1 : electronNumPerBunch;
    

    for(int i=0;i<haissinski->nz;i++)
    {
        profile[i] =  profile[i] / norm * electronNumPerBunchTemp; 
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

     
    // the rf voltage and phase from tracking  
    for(int i=0;i<haissinski->nz;i++)
    {
        for(int j=0;j<cavityResonator.resonatorVec.size();j++)
        {          
            resHarm = cavityResonator.resonatorVec[j].resHarm;   
            vTotRF[i] += haissinski->cavAmp[j] * cos( - 2. * PI * ringHarmH * f0 * resHarm * haissinski->bunchPosZ[i] /  ( rBeta * CLight)  + haissinski->cavPhase[j] );            
        }
        haissinski->vRF[i]  = vTotRF[i] - u0;
    }    
    
    double coeffDelta = f0 / pow(rBeta,2) / electronBeamEnergy;   // 1/[V]/[s];

    //Eq.(2~8), PRAB 17 064401, and S. Y. Lee Eq. (3.35 and 3.36)  ddelta/dt = - dH/dphi  and (dz=- c * dt = - c * dphi / (2 * PI * h * f0) )
    haissinski->rfHamiltonian[0] = 0.e0;
    for(int i=1;i<haissinski->nz;i++)
    {              
        haissinski->rfHamiltonian[i] = haissinski->rfHamiltonian[i-1] - ((vTotRF[i] + vTotRF[i-1]) / 2.0 -  u0 ) * coeffDelta 
                                     * (-1) *  2. * PI * ringHarmH * f0 * haissinski->dz / ( rBeta * CLight);                        // [1/s]          
   

        // haissinski->rfHamiltonian[i] = haissinski->rfHamiltonian[i-1] - (vTotRF[i-1] -  u0 ) * coeffDelta 
        //                              * (-1) *  2. * PI * ringHarmH * f0 * haissinski->dz / ( rBeta * CLight); 
    
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


void Bunch::GetParticleLongitudinalPhaseSpace1(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,int bunchIndex)
{
    
    tk::spline wakePotenFit, totHamiltonianFit;
    wakePotenFit.set_points(haissinski->bunchPosZ,haissinski->totWakePoten,tk::spline::cspline);           // fitting the wakePoten 1/[m]
    totHamiltonianFit.set_points(haissinski->bunchPosZ,haissinski->totHamiltonian,tk::spline::cspline);    // fitting the totHamilton 1/[m]

    double zMax = haissinski->averZ + haissinski->rmsZ * 5 + 0.1;
    double zMin = haissinski->averZ - haissinski->rmsZ * 5 - 0.09;
    
    // zMin = haissinski->bunchPosZ[0];
    // zMax = haissinski->bunchPosZ.back(); 

    // zMax = + 0.075;
    // zMin = - 0.075;
    
    int nz      = 41;                               // total have 41 partilces at differet initial conditions 
    double dz   = (zMax - zMin) / nz;     

    // nz = haissinski->nz; 
    // dz = haissinski->dz;
 
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
    fout<<"&parameter name=deltaS,      units=m,            type=float,  &end"<<endl;

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

        for(int k=0;k<100*turnsLongiOscilation;k++)              
        {            
            longiTrajZeta[i].push_back(zeta1);
            zeta2 = LeapFrog(inputParameter,cavityResonator,zeta1,wakePotenFit);                      

            // actionJ[i] +=  abs(zeta2[0] - zeta1[0]) *  abs( zeta2[1] )      / (2 * PI);            
            // deltaS     +=  abs(zeta2[0] - zeta1[0]) /  abs( zeta2[1] )      / eta ;

            actionJ[i] +=  abs(zeta2[0] - zeta1[0]) * abs(zeta2[1] + zeta1[1]) / 2  / (2 * PI);         // m   
            deltaS     +=  abs(zeta2[0] - zeta1[0]) / abs(zeta2[1] + zeta1[1]) * 2  / eta ;             // m
            
            p1 = zeta2[1]; 
            zeta1 = zeta2;
            if(p0 * p1<0)    //ensure the longitudinal phase space only rotate one-trun (360 degree).
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
        fout<<deltaS<<endl;
        fout<<longiTrajZeta[i].size()<<endl;

        for(int k=0;k<longiTrajZeta[i].size();k++)
        {
            fout<<setw(15)<<left<<k
                <<setw(15)<<left<<longiTrajZeta[i][k][0]
                <<setw(15)<<left<<longiTrajZeta[i][k][1]
                <<endl;
        }
        haissinski->nus[i]     =  nus[i];
        haissinski->actionJ[i] =  actionJ[i];      
    }


    // get the average and rms longitudinal nus
    double averNus = 0;
    double norm=0 ;

    for(int i=0;i<nz;i++)
    {
        norm    +=  haissinski->bunchProfile[i] * dz;
        averNus +=  haissinski->bunchProfile[i] * nus[i] * dz; 
    }
    averNus /= norm;
    
    double rmsNus  = 0;     
    for(int i=0;i<nz;i++)
    {
        rmsNus  +=  haissinski->bunchProfile[i] * pow(nus[i] - averNus,2) * dz ; 
    }
    rmsNus = sqrt(rmsNus / norm);

    haissinski->averNus = averNus;
    haissinski->rmsNus  = rmsNus;

    fout.close();        
}




void Bunch::GetParticleLongitudinalPhaseSpace(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,int bunchIndex)
{
    tk::spline wakePotenFit, totHamiltonianFit;
    wakePotenFit.set_points(haissinski->bunchPosZ,haissinski->totWakePoten,tk::spline::cspline);           // fitting the wakePoten 1/[m]
    totHamiltonianFit.set_points(haissinski->bunchPosZ,haissinski->totHamiltonian,tk::spline::cspline);    // fitting the totHamilton 1/[m]

    double zMax = haissinski->averZ + haissinski->rmsZ * 5 + 0.1;
    double zMin = haissinski->averZ - haissinski->rmsZ * 5 - 0.09;
    
    // zMin = haissinski->bunchPosZ[0];
    // zMax = haissinski->bunchPosZ.back(); 

    // zMax = + 0.075;
    // zMin = - 0.075;
    
    int nz      = 41;                               // total have 41 partilces at differet initial conditions 
    double dz   = (zMax - zMin) / nz;     

    // nz = haissinski->nz; 
    // dz = haissinski->dz;
 
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
    fout<<"&parameter name=deltaS,      units=m,            type=float,  &end"<<endl;

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

        for(int k=0;k<100*turnsLongiOscilation;k++)              
        {            
            longiTrajZeta[i].push_back(zeta1);
            zeta2 = LeapFrog(inputParameter,cavityResonator,zeta1,wakePotenFit);                      

            // actionJ[i] +=  abs(zeta2[0] - zeta1[0]) *  abs( zeta2[1] )      / (2 * PI);            
            // deltaS     +=  abs(zeta2[0] - zeta1[0]) /  abs( zeta2[1] )      / eta ;

            actionJ[i] +=  abs(zeta2[0] - zeta1[0]) * abs(zeta2[1] + zeta1[1]) / 2  / (2 * PI);         // m   
            deltaS     +=  abs(zeta2[0] - zeta1[0]) / abs(zeta2[1] + zeta1[1]) * 2  / eta ;             // m
            
            p1 = zeta2[1]; 
            zeta1 = zeta2;
            if(p0 * p1<0)    //ensure the longitudinal phase space only rotate one-trun (360 degree).
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
        fout<<deltaS<<endl;
        fout<<longiTrajZeta[i].size()<<endl;

        for(int k=0;k<longiTrajZeta[i].size();k++)
        {
            fout<<setw(15)<<left<<k
                <<setw(15)<<left<<longiTrajZeta[i][k][0]
                <<setw(15)<<left<<longiTrajZeta[i][k][1]
                <<endl;
        }
        haissinski->nus[i]     =  nus[i];
        haissinski->actionJ[i] =  actionJ[i];      
    }


    // get the average and rms longitudinal nus
    double averNus = 0;
    double norm=0 ;

    for(int i=0;i<nz;i++)
    {
        norm    +=  haissinski->bunchProfile[i] * dz;
        averNus +=  haissinski->bunchProfile[i] * nus[i] * dz; 
    }
    averNus /= norm;
    
    double rmsNus  = 0;     
    for(int i=0;i<nz;i++)
    {
        rmsNus  +=  haissinski->bunchProfile[i] * pow(nus[i] - averNus,2) * dz ; 
    }
    rmsNus = sqrt(rmsNus / norm);

    haissinski->averNus = averNus;
    haissinski->rmsNus  = rmsNus;

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
