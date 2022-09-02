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
    delete cavFBCenInfo;   
}


void Bunch::Initial(const  ReadInputSettings &inputParameter)
{    
    macroEleNumPerBunch = inputParameter.ringBunchPara->macroEleNumPerBunch;      // macroEleNumPerBunch=1 -- electron beam weak-strong model
    current             = inputParameter.ringBunchPara->current;                                                                              
    electronEnergy      = inputParameter.ringParBasic->electronBeamEnergy;
    rmsEnergySpread     = inputParameter.ringBunchPara-> rmsEnergySpread;
    rmsBunchLength      = inputParameter.ringBunchPara-> rmsBunchLength;
    rGamma              = inputParameter.ringParBasic->rGamma;
    rBeta               = inputParameter.ringParBasic->rBeta;
    
    // in sp model, it does not change during the tracking 
    emittanceZ      = inputParameter.ringBunchPara->emittanceZ;
    emittanceX      = inputParameter.ringBunchPara->emittanceX;
    emittanceY      = inputParameter.ringBunchPara->emittanceY;
    
    
    electronNumPerBunch = current / inputParameter.ringParBasic->f0  / ElectronCharge;
    macroEleCharge      = electronNumPerBunch / macroEleNumPerBunch;

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
    eSurive.resize(macroEleNumPerBunch,1);
    
    
    cavFBCenInfo->cavVolBunchCen.resize(inputParameter.ringParRf->resNum);
    cavFBCenInfo->cavAmpBunchCen.resize(inputParameter.ringParRf->resNum);
    cavFBCenInfo->genVolBunchAver.resize(inputParameter.ringParRf->resNum);
    cavFBCenInfo->cavPhaseBunchCen.resize(inputParameter.ringParRf->resNum);
    cavFBCenInfo->induceVolBunchCen.resize(inputParameter.ringParRf->resNum);
    cavFBCenInfo->selfLossVolBunchCen.resize(inputParameter.ringParRf->resNum);


    
    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {        
        cavFBCenInfo->genVolBunchAver[i]     =complex<double>(0.e0,0.e0);
        cavFBCenInfo->induceVolBunchCen[i]   =complex<double>(0.e0,0.e0);
        cavFBCenInfo->selfLossVolBunchCen[i] =complex<double>(0.e0,0.e0);
        cavFBCenInfo->cavVolBunchCen[i]      =complex<double>(0.e0,0.e0);
        cavFBCenInfo->cavAmpBunchCen[i]      =0.E0;
        cavFBCenInfo->cavPhaseBunchCen[i]    =0.E0;                
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


        // ztemp  = latticeInterActionPoint.zTransferMatrix[k][0] * ePositionZ[i]
        //        + latticeInterActionPoint.zTransferMatrix[k][1] * eMomentumZ[i];

        // zPtemp = latticeInterActionPoint.zTransferMatrix[k][2] * ePositionZ[i]
        //        + latticeInterActionPoint.zTransferMatrix[k][3] * eMomentumZ[i];

        // ePositionZ[i] = ztemp;
        // eMomentumZ[i] = zPtemp;

        double lossTemp =pow(ePositionX[i]/latticeInterActionPoint.pipeAperatureX[k],2) + pow(ePositionY[i]/latticeInterActionPoint.pipeAperatureY[k],2);

        if( lossTemp > 1)
        {            
            eSurive[i] = 0;  // particle is lost in transverse when  eSurive[i]==0
        }
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
//    gsl_matrix_set_identity(R2);                  // full coupling
    gsl_matrix_set_zero(R2);                        // no coupling at the point in SR calculation

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



    coeff[0] = sqrt(1 - pow(lambda[0],2)) * sqrt(emittanceX);
    coeff[1] = sqrt(1 - pow(lambda[1],2)) * sqrt(emittanceY);
    coeff[2] = sqrt(1 - pow(lambda[2],4)) * sqrt(emittanceZ);

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
        tempZ  = tempZ  			;
        tempPZ = tempPZ             ;
        //tempPZ = tempPZ * pow(lambda[2],2) ;


        //(1.2) synchron-radiation-excitiation (transverse only and multi-particle case only)

        if(macroEleNumPerBunch!=1)
        {
            for(int j=0;j<6;j++)
            {
                randR[j]=Gaussrand(1,0,100.0*j);
            }

            tempX  +=  coeff[0] * randR[0];
            tempPX +=  coeff[0] * randR[1];
            tempY  +=  coeff[1] * randR[2];
            tempPY +=  coeff[1] * randR[3];
            tempZ   =  tempZ ;
            tempPZ  =  tempPZ;
            //tempPZ = tempPZ  + coeff[2] * randR[5];
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

void Bunch::BunchTransferDueToWake()
{
    for(int i=0;i<ePositionX.size();i++)
    {
        /*    
        if(wakeForceAver[0]>1.E-10) eMomentumX[i] -= wakeForceAver[0];
        if(wakeForceAver[1]>1.E-10) eMomentumY[i] -= wakeForceAver[1];
        if(wakeForceAver[2]>1.E-10) eMomentumZ[i] -= wakeForceAver[2];
        */
        //cout<<setw(15)<<left<<xAver
        //   <<setw(15)<<left<<lRWakeForceAver[0]
        //    <<setw(15)<<left<<lRWakeForceAver[1]
        //    <<setw(15)<<left<<lRWakeForceAver[2]
        //    <<setw(15)<<left<<__LINE__<<__FILE__<<endl;
        //getchar();
        eMomentumX[i] += lRWakeForceAver[0];      //rad
        eMomentumY[i] += lRWakeForceAver[1];    
        eMomentumZ[i] += lRWakeForceAver[2];
        
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


