//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             
#include "LatticeInterActionPoint.h"
#include "Global.h"
#include <vector>
#include <stdlib.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <time.h>
#include <random>
#include<string>
#include<iomanip>



LatticeInterActionPoint::LatticeInterActionPoint()
{
}


LatticeInterActionPoint::~LatticeInterActionPoint()
{
}

void LatticeInterActionPoint::Initial(const ReadInputSettings &inputParameter)
{

    //0) get the input from the inputParameter 
    numberOfInteraction       = inputParameter.ringIonEffPara->numberofIonBeamInterPoint;             
    ionMaxNumberOneInterPoint = inputParameter.ringIonEffPara->ionMaxNumber;
    ionLossBoundary           = inputParameter.ringIonEffPara->ionLossBoundary;

    gasSpec                   = inputParameter.ringIonEffPara->gasSpec;
    ionMassNumber             = inputParameter.ringIonEffPara->ionMassNumber; 
    corssSectionEI            = inputParameter.ringIonEffPara->corssSectionEI;
    gasPercent                = inputParameter.ringIonEffPara->gasPercent; 
    
    circRing        =  inputParameter.ringParBasic->circRing;
    harmonics       =  inputParameter.ringParBasic->harmonics;
 
    double temp=0;  
    for(int i=0;i<gasSpec;i++)
    {
        temp += gasPercent[i]; 
    }
    if(temp!=1.e0)
    {
        cerr<<"non 100 percent gas"<<endl;
        exit(0);
    }
    

   
    //1) twiss parameter is divided by number of interaction points    
    twissAlphaX.resize(numberOfInteraction);
    twissBetaX .resize(numberOfInteraction);
    twissAlphaY.resize(numberOfInteraction);
    twissBetaY .resize(numberOfInteraction);
    twissAlphaZ.resize(numberOfInteraction);
    twissBetaZ .resize(numberOfInteraction);
    pipeAperatureX.resize(numberOfInteraction);
    pipeAperatureY.resize(numberOfInteraction);
    
    
    twissDispX .resize(numberOfInteraction);				
    twissDispPX.resize(numberOfInteraction);;          
    twissDispY .resize(numberOfInteraction);;				
    twissDispPY.resize(numberOfInteraction);;             

    xPhaseAdv.resize(numberOfInteraction);
    yPhaseAdv.resize(numberOfInteraction);
    zPhaseAdv.resize(numberOfInteraction);
    
    xTransferMatrix.resize(numberOfInteraction);
    yTransferMatrix.resize(numberOfInteraction);
    zTransferMatrix.resize(numberOfInteraction);
    
    transferMatrix.resize(numberOfInteraction);
    

    // in current assumption the 6*6 transfer matrix does not includes couling between either two degrees of freedom   
    for(int i=0;i<numberOfInteraction;i++)
    {
        xTransferMatrix[i].resize(4);
        yTransferMatrix[i].resize(4);
        zTransferMatrix[i].resize(4);
    }
    
    // end

   
    //2.1) ion data at kth interaction point
    temperature      .resize(numberOfInteraction);
    vacuumPressure   .resize(numberOfInteraction); 
    interactionLength.resize(numberOfInteraction);
    totMacroIonsAtInterPoint.resize(numberOfInteraction);
    totIonChargeAtInterPoint.resize(numberOfInteraction);


    //2.2) ion data at kth interaction point and pth ion species
    vacuumPressureEachGas.resize(numberOfInteraction); 
    ionLineDensity   .resize(numberOfInteraction);        
    macroIonNumber.resize(numberOfInteraction);
    ionNumber     .resize(numberOfInteraction);
    macroIonCharge.resize(numberOfInteraction);
    
    ionPositionX.resize(numberOfInteraction);
    ionPositionY.resize(numberOfInteraction);
    ionVelocityX.resize(numberOfInteraction);
    ionVelocityY.resize(numberOfInteraction);
    
    ionAccumuFx.resize(numberOfInteraction);
    ionAccumuFy.resize(numberOfInteraction);
    
    ionAccumuNumber   .resize(numberOfInteraction);
    ionAccumuPositionX.resize(numberOfInteraction);
    ionAccumuPositionY.resize(numberOfInteraction);
    ionAccumuVelocityX.resize(numberOfInteraction);
    ionAccumuVelocityY.resize(numberOfInteraction);
    
    ionAccumuAverX.resize(numberOfInteraction);
    ionAccumuAverY.resize(numberOfInteraction);
    ionAccumuRMSX .resize(numberOfInteraction);
    ionAccumuRMSY .resize(numberOfInteraction); 
    ionAccumuAverVelX.resize(numberOfInteraction); 
    ionAccumuAverVelY.resize(numberOfInteraction); 

    for(int i=0;i<numberOfInteraction;i++)
    {
        vacuumPressureEachGas[i].resize(gasSpec);
        ionLineDensity[i]       .resize(gasSpec);
        ionNumber[i]            .resize(gasSpec);
        macroIonNumber[i]       .resize(gasSpec);
        macroIonCharge[i]       .resize(gasSpec);

        ionPositionX[i].resize(gasSpec);
        ionPositionY[i].resize(gasSpec);
        ionVelocityX[i].resize(gasSpec);
        ionVelocityY[i].resize(gasSpec);

        ionAccumuFx[i].resize(gasSpec);
        ionAccumuFy[i].resize(gasSpec);
        
        ionAccumuNumber[i].resize(gasSpec);
        ionAccumuPositionX[i].resize(gasSpec);
        ionAccumuPositionY[i].resize(gasSpec);
        ionAccumuVelocityX[i].resize(gasSpec);
        ionAccumuVelocityY[i].resize(gasSpec);   
        
        ionAccumuAverX[i].resize(gasSpec);
        ionAccumuAverY[i].resize(gasSpec);
        ionAccumuRMSX[i].resize(gasSpec);
        ionAccumuRMSY[i].resize(gasSpec);  
        ionAccumuAverVelX[i].resize(gasSpec);
        ionAccumuAverVelY[i].resize(gasSpec);                                                 
    }
    
    InitialLattice(inputParameter);   


}

void LatticeInterActionPoint::SetLatticeBRHForSynRad(const ReadInputSettings &inputParameter)
{
    for(int i=0;i<3;i++)
    {
        latticeSynRadCoff[i] = exp(-1.0/inputParameter.ringParBasic->synchRadDampTime[i]);
    }
    latticeSynRadCoff[2] = pow(latticeSynRadCoff[2],2);
  
    latticeSynRadCoff[3] = sqrt(1 - pow(latticeSynRadCoff[0],2)) * sqrt(inputParameter.ringParBasic->emitNat[0]);
    latticeSynRadCoff[4] = sqrt(1 - pow(latticeSynRadCoff[1],2)) * sqrt(inputParameter.ringParBasic->emitNat[1]);
    latticeSynRadCoff[5] = sqrt(1 - pow(latticeSynRadCoff[2],2)) * sqrt(inputParameter.ringParBasic->emitNat[2]);


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

    gsl_matrix_set(dispMatrixHX,0,1,twissDispX[k]);
    gsl_matrix_set(dispMatrixHX,1,1,twissDispPX[k]);

    gsl_matrix_set(dispMatrixHY,0,1,twissDispY[k]);
    gsl_matrix_set(dispMatrixHY,1,1,twissDispPY[k]);


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
    // gsl_matrix_set(R2,1,0,inputParameter.ringParBasic->tengR21);
    double det = gsl_matrix_get(R2,0,0) * gsl_matrix_get(R2,1,1) - gsl_matrix_get(R2,1,0) * gsl_matrix_get(R2,0,1);
    double b = sqrt(1 - det);
    
    gsl_matrix * R12  = gsl_matrix_alloc (2, 2);
    gsl_matrix * R21  = gsl_matrix_alloc (2, 2);

    // set R12
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

    a00 = 1.0/sqrt(twissBetaX[k]);
    a10 = twissAlphaX[k] / sqrt(twissBetaX[k]);
    a11 = sqrt(twissBetaX[k]);


    gsl_matrix_set(twissMatrixB,0,0,a00);
    gsl_matrix_set(twissMatrixB,1,0,a10);
    gsl_matrix_set(twissMatrixB,1,1,a11);


    a00 = 1.0/sqrt(twissBetaY[k]);
    a10 = twissAlphaY[k] / sqrt(twissBetaY[k]);
    a11 = sqrt(twissBetaY[k]);


    gsl_matrix_set(twissMatrixB,2,2,a00);
    gsl_matrix_set(twissMatrixB,3,2,a10);
    gsl_matrix_set(twissMatrixB,3,3,a11);


    a00 = 1.0/sqrt(twissBetaZ[k]);
    a10 = twissAlphaZ[k] / sqrt(twissBetaZ[k]);
    a11 = sqrt(twissBetaZ[k]);

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


    // test M^{T}.J6.M==J6
    // gsl_matrix *  mat66Temp = gsl_matrix_alloc (6, 6);
    // gsl_matrix_transpose(cordTransfer);
    // gsl_matrix_mul(cordTransfer,sympleMarixJ6,cordTransferTemp);
    // gsl_matrix_transpose(cordTransfer);
    // gsl_matrix_mul(cordTransferTemp,cordTransfer,mat66Temp);

    // for(int i=0;i<6;i++)
    // {
    //     for(int j=0;j<6;j++)
    //     {
            
    //         cout<<setw(15)<<left<<gsl_matrix_get(invCordTransfer,i,j) ;
    //     }
    //     cout<<endl;
    // }
    // cout<<__FILE__<<endl;
    // getchar();



    for(int i=0;i<6;i++)
    {
        for(int j=0;j<6;j++)
        {
            latticeSynRadBRH[6*i+j] = gsl_matrix_get(cordTransfer,i,j);
            // cout<<setw(15)<<left<<latticeSynRadBRH[6*i+j];
        }
        // cout<<endl;
    }
    // getchar();

    for(int i=0;i<6;i++)
    {
        for(int j=0;j<6;j++)
        {
            latticeSynRadBRH[6*i+j+36] = gsl_matrix_get(invCordTransfer,i,j);
        }
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
}

void LatticeInterActionPoint::SetLatticeParaForOneTurnMap(const ReadInputSettings &inputParameter)
{
    

    latticeParaForOneTurnMap[0] = twissAlphaX[0];                         latticeParaForOneTurnMap[1] = twissAlphaY[0];
    latticeParaForOneTurnMap[2] = twissBetaX[0] ;                         latticeParaForOneTurnMap[3] = twissBetaY[0] ;
    latticeParaForOneTurnMap[4] = inputParameter.ringParBasic->workQx;    latticeParaForOneTurnMap[5] = inputParameter.ringParBasic->workQy;
    latticeParaForOneTurnMap[6] = inputParameter.ringParBasic->chrom[0];  latticeParaForOneTurnMap[7] = inputParameter.ringParBasic->chrom[1];
    latticeParaForOneTurnMap[8] = twissDispX[0];                          latticeParaForOneTurnMap[9] = twissDispY[0];
    latticeParaForOneTurnMap[10]= twissDispPX[0];                         latticeParaForOneTurnMap[11]= twissDispPY[0];
    latticeParaForOneTurnMap[12]= inputParameter.ringParBasic->aDTX[0];   latticeParaForOneTurnMap[13]= inputParameter.ringParBasic->aDTX[1];
    latticeParaForOneTurnMap[14]= inputParameter.ringParBasic->aDTY[0];   latticeParaForOneTurnMap[15]= inputParameter.ringParBasic->aDTY[1];
    latticeParaForOneTurnMap[16]= inputParameter.ringParBasic->aDTXY[0];  latticeParaForOneTurnMap[17]= inputParameter.ringParBasic->aDTXY[1];
    latticeParaForOneTurnMap[18]= inputParameter.ringParBasic->alphac[0]; latticeParaForOneTurnMap[19]= inputParameter.ringParBasic->alphac[1]; 
    latticeParaForOneTurnMap[20]= inputParameter.ringParBasic->alphac[2];
    latticeParaForOneTurnMap[21]= inputParameter.ringParBasic->circRing;
  
    // for(int i=0;i<22;i++)
    //     cout<<i<<"  "<<latticeParaForOneTurnMap[i]<<endl;

    // getchar();
}


void LatticeInterActionPoint::InitialLattice(const ReadInputSettings &inputParameter)
{

    double workQx = inputParameter.ringParBasic->workQx;
    double workQy = inputParameter.ringParBasic->workQy;
    double workQz = inputParameter.ringParBasic->workQz;
    double xAperture = inputParameter.ringParBasic->pipeAperature[1];
    double yAperture = inputParameter.ringParBasic->pipeAperature[2];
    
    ifstream fin(inputParameter.ringIonEffPara->twissInput);
    if (! fin.is_open())
    {
        cerr<< "Error opening file "<<inputParameter.ringIonEffPara->twissInput<<endl;
        exit (1);
    }

    
    string str;
    vector<string> strVec;
    
    int index = 0;
    int i=0;
    while (!fin.eof())
    {
        getline(fin,str);

        if(index<=4)
        {
            index++;
            continue;
        }
        

        string stringTest;
        for(int i=0;i<str.size();i++)
        {
           stringTest.push_back(' ');		
        }		
        if(stringTest == str )  continue;
                    
        StringSplit2(str, strVec);
                
        twissBetaX[i]  = stod(strVec[1]);           
        twissAlphaX[i] = stod(strVec[2]);
        xPhaseAdv[i]   = stod(strVec[3]);
        twissDispX[i]  = stod(strVec[4]);  
        twissDispPX[i] = stod(strVec[5]);                    
        pipeAperatureX[i] = stod(strVec[6]);
                        
        twissBetaY[i]  = stod(strVec[7]);            
        twissAlphaY[i] = stod(strVec[8]);
        yPhaseAdv[i]   = stod(strVec[9]);
        twissDispY[i]  = stod(strVec[10]);
        twissDispPY[i] = stod(strVec[11]);
        pipeAperatureY[i] = stod(strVec[12]);
        interactionLength[i] = stod(strVec[13]);
        vacuumPressure[i] =  stod(strVec[14])* 1.0E-9 * 133.3224;               // Torr to Pascals    
        temperature[i] = stod(strVec[15]) ;
 

        twissBetaZ[i]  = inputParameter.ringParBasic->naturalBunchLength / inputParameter.ringParBasic->sdelta0;   //ref. Zhang Yuan's paper, have to equibrium value.  
        twissAlphaZ[i] = 0.0;
        zPhaseAdv[i]   = interactionLength[i] * 2 * PI * workQz;   
   
        xAperture>pipeAperatureX[i] ? (pipeAperatureX[i]=pipeAperatureX[i]) : (pipeAperatureX[i]=xAperture);
        yAperture>pipeAperatureY[i] ? (pipeAperatureY[i]=pipeAperatureY[i]) : (pipeAperatureY[i]=yAperture);

        i++;
        index++;
    }


    fin.close();


    vector<double> intLengthTemp=interactionLength;    
    for (int i=0; i<interactionLength.size();i++)
    {
        if(i != interactionLength.size()-1 )
        {
            interactionLength[i] = (intLengthTemp[i+1] - intLengthTemp[i]) * circRing ;        
        }
        else
        {
            interactionLength[i] = (1.0                - intLengthTemp[i]) * circRing ;  
        }                      
    }
 
    
    if((i)!=numberOfInteraction)
    {
        cerr<<"data of numberofIonBeamInterPoint in input.dat and InterPointParameter.dat files does not match";
        exit(0);
    }
    


    double betaX1;
    double betaX2;
    double alphaX1;
    double alphaX2;
    
    double betaY1;
    double betaY2;
    double alphaY1;
    double alphaY2;
  
	double betaZ1;
    double betaZ2;
    double alphaZ1;
    double alphaZ2;
  
    
    
    double phaseAdvanceX;
    double phaseAdvanceY;
    double phaseAdvanceZ;
    

    for(int i=0;i<numberOfInteraction;i++)
    {
	           
        alphaX1 = twissAlphaX[i];
        betaX1  = twissBetaX[i];
        alphaY1 = twissAlphaY[i];
        betaY1  = twissBetaY[i];
        alphaZ1 = twissAlphaZ[i];
        betaZ1  = twissBetaZ[i];


        if(i<numberOfInteraction-1)
        {

            alphaX2 = twissAlphaX[i+1];
            betaX2  = twissBetaX[i+1];
            alphaY2 = twissAlphaY[i+1];
            betaY2  = twissBetaY[i+1];
            alphaZ2 = twissAlphaZ[i+1];
            betaZ2  = twissBetaZ[i+1];
            
            phaseAdvanceX =   xPhaseAdv[i+1]  -   xPhaseAdv[i];
            phaseAdvanceY =   yPhaseAdv[i+1]  -   yPhaseAdv[i]; 
            phaseAdvanceZ =   zPhaseAdv[i+1]  -   zPhaseAdv[i]; 
               
        }
        else
        {
        
            alphaX2 = twissAlphaX[0];
            betaX2  = twissBetaX[0];
            alphaY2 = twissAlphaY[0];
            betaY2  = twissBetaY[0];
            alphaZ2 = twissAlphaZ[0];
            betaZ2  = twissBetaZ[0];
         
            phaseAdvanceX =   2*PI*workQx  -   xPhaseAdv[i];
            phaseAdvanceY =   2*PI*workQy  -   yPhaseAdv[i];
            phaseAdvanceZ =   2*PI*workQz  -   zPhaseAdv[i];
              
        }

                 
                                             
        xTransferMatrix[i][0]  = sqrt(betaX2 / betaX1) * (cos(phaseAdvanceX) + alphaX1 * sin(phaseAdvanceX));
        xTransferMatrix[i][1]  = sqrt(betaX2 * betaX1) *  sin(phaseAdvanceX);
        xTransferMatrix[i][2]  = -(1 + alphaX1 * alphaX2)/sqrt(betaX1 * betaX2) * sin(phaseAdvanceX) 
                                 +(    alphaX1 - alphaX2)/sqrt(betaX1 * betaX2) * cos(phaseAdvanceX);
        xTransferMatrix[i][3]  = sqrt(betaX1 / betaX2) * (cos(phaseAdvanceX) - alphaX2 * sin(phaseAdvanceX));

        yTransferMatrix[i][0]  = sqrt(betaY2 / betaY1) * (cos(phaseAdvanceY) + alphaY1 * sin(phaseAdvanceY));
        yTransferMatrix[i][1]  = sqrt(betaY2 * betaY1) * sin(phaseAdvanceY);
        yTransferMatrix[i][2]  = -(1 + alphaY1 * alphaY2)/sqrt(betaY1 * betaY2) * sin(phaseAdvanceY) 
                                 +(    alphaY1 - alphaY2)/sqrt(betaY1 * betaY2) * cos(phaseAdvanceY);
        yTransferMatrix[i][3]  = sqrt(betaY1 / betaY2) * (cos(phaseAdvanceY) - alphaY2 * sin(phaseAdvanceY));
			
        zTransferMatrix[i][0]  = sqrt(betaZ2 / betaZ1) * (cos(phaseAdvanceZ) + alphaZ1 * sin(phaseAdvanceZ));
        zTransferMatrix[i][1]  = sqrt(betaZ2 * betaZ1) * sin(phaseAdvanceZ);
        zTransferMatrix[i][2]  = -(1 + alphaZ1 * alphaZ2)/sqrt(betaZ1 * betaZ2) * sin(phaseAdvanceZ) 
                                 +(    alphaZ1 - alphaZ2)/sqrt(betaZ1 * betaZ2) * cos(phaseAdvanceZ);
        zTransferMatrix[i][3]  = sqrt(betaZ1 / betaZ2) * (cos(phaseAdvanceZ) - alphaZ2 * sin(phaseAdvanceZ));
        
        // cout<<setw(15)<<left<<zTransferMatrix[i][0]<<setw(15)<<left<<zTransferMatrix[i][1]<<endl;
        // cout<<setw(15)<<left<<zTransferMatrix[i][2]<<setw(15)<<left<<zTransferMatrix[i][3]<<endl;
        // cout<<setw(15)<<left<<yTransferMatrix[i][0]<<setw(15)<<left<<yTransferMatrix[i][1]<<endl;
        // cout<<setw(15)<<left<<yTransferMatrix[i][2]<<setw(15)<<left<<yTransferMatrix[i][3]<<endl;
        // cout<<setw(15)<<left<<xTransferMatrix[i][0]<<setw(15)<<left<<xTransferMatrix[i][1]<<endl;
        // cout<<setw(15)<<left<<xTransferMatrix[i][2]<<setw(15)<<left<<xTransferMatrix[i][3]<<endl;        
        // cout<< inputParameter.ringBunchPara->rmsBunchLength<<endl;
        // cout<< inputParameter.ringBunchPara->rmsEnergySpread<<endl;      
        // getchar();
        
    }

  
    for(int k=0; k<numberOfInteraction;k++)
    {
        for(int j=0;j<gasSpec;j++)
        {        
            macroIonNumber[k][j]=inputParameter.ringIonEffPara->macroIonNumberGeneratedPerIP;          // at each interaction point, macroIonNumber[k] macroIon are generated 
	    }
	}

    for(int k=0;k<numberOfInteraction;k++)
    {
        for(int j=0;j<gasSpec;j++)
        {
            ionPositionX[k][j].resize(macroIonNumber[k][j]);
            ionPositionY[k][j].resize(macroIonNumber[k][j]);
            ionVelocityX[k][j].resize(macroIonNumber[k][j]);
            ionVelocityY[k][j].resize(macroIonNumber[k][j]);
        }
    }       


    for(int k=0;k<numberOfInteraction;k++)     
    {
        for (int j=0;j<gasSpec;j++)
        {
            vacuumPressureEachGas[k][j] = gasPercent[j] * vacuumPressure[k];  // related to the gasPercent  settings.
            //cout<< gasPercent[j] <<"    "<<vacuumPressureEachGas[k][j]<<endl;
        }
    }   


    
}



void LatticeInterActionPoint::GetIonNumberPerInterAction(double electronNumPerBunch, int k)
{

    for(int p=0;p<gasSpec;p++)
    {
        ionLineDensity[k][p] = corssSectionEI[p] * vacuumPressureEachGas[k][p] / temperature[k] /Boltzmann * electronNumPerBunch;
        ionNumber[k][p]      = ionLineDensity[k][p] * interactionLength[k];
    }
}



void LatticeInterActionPoint::IonGenerator(double rmsRx, double rmsRy, double xAver,double yAver, int k)
{
    double tempx;
    double tempy;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    
    std::normal_distribution<> dx{xAver,rmsRx};
    std::normal_distribution<> dy{yAver,rmsRy};

    for(int p=0;p<gasSpec;p++)
    {
        macroIonCharge[k][p] = ionNumber[k][p]/macroIonNumber[k][p]; 
                    
        int i=0;
        while(i<macroIonNumber[k][p])
        {
            tempx = dx(gen);
            tempy = dy(gen);

            if( pow( (tempx-xAver)/rmsRx, 2)  + pow( (tempy-yAver)/rmsRy, 2) > 9.E0  ) // ions generated within 3 rms beam size.
            {
               continue;
            }
                           
            ionPositionX[k][p][i]=tempx;
            ionPositionY[k][p][i]=tempy;
            ionVelocityX[k][p][i]=0.E0;
            ionVelocityY[k][p][i]=0.E0;

            i=i+1;
        }   
   }     
}

void LatticeInterActionPoint::IonsUpdate(int k)
{
    // add the macroIonNumber[k][p] new generated ions to accumulated ions
    
    for(int p=0;p<gasSpec;p++)
    {
        for(int i=0;i<macroIonNumber[k][p];i++)
        {               
            ionAccumuPositionX[k][p].push_back(ionPositionX[k][p][i]);
            ionAccumuPositionY[k][p].push_back(ionPositionY[k][p][i]);
            ionAccumuVelocityX[k][p].push_back(ionVelocityX[k][p][i]);
            ionAccumuVelocityY[k][p].push_back(ionVelocityY[k][p][i]);
            ionAccumuFx[k][p].push_back(0.E0);
            ionAccumuFy[k][p].push_back(0.E0);  
        }
    
        ionAccumuNumber[k][p]  =  ionAccumuPositionX[k][p].size(); 
    }

    
    for(int p=0;p<gasSpec;p++)
    {
        if(ionAccumuNumber[k][p]>ionMaxNumberOneInterPoint)
        {
            cerr<<"the accumulated ions at "<<p<<"th interaction point is larger than limit "<<ionMaxNumberOneInterPoint<<endl;
            cerr<<"try reduce ions generate per interaction or enlarge the ionCalIonMaxNumber "<<endl;
            exit(0);
        }
    }
}


void LatticeInterActionPoint::IonRMSCal(int k)
{
    for(int p=0; p<gasSpec;p++)
    {
        IonRMSCal(k, p);
    }
}

void LatticeInterActionPoint::IonRMSCal(int k, int p)
{
    double x2Aver=0;
    double y2Aver=0;

    ionAccumuNumber[k][p]=ionAccumuPositionX[k][p].size();
            
    if(ionAccumuNumber[k][p]==0)
    {
        x2Aver  =  0;
        y2Aver  =  0;
    }
    else if(ionAccumuNumber[k][p]==1)
    {
        x2Aver  =  pow(ionAccumuPositionX[k][p][0],2);
        y2Aver  =  pow(ionAccumuPositionY[k][p][0],2);
    }
    else
    {
        ionAccumuAverX[k][p]   = accumulate(begin(ionAccumuPositionX[k][p]), end(ionAccumuPositionX[k][p]), 0.0) / ionAccumuNumber[k][p];
        ionAccumuAverY[k][p]   = accumulate(begin(ionAccumuPositionY[k][p]), end(ionAccumuPositionY[k][p]), 0.0) / ionAccumuNumber[k][p]; 

        for(int i=0;i<ionAccumuNumber[k][p];i++)
        {
            x2Aver  +=  pow(ionAccumuPositionX[k][p][i]-ionAccumuAverX[k][p] ,2);
            y2Aver  +=  pow(ionAccumuPositionY[k][p][i]-ionAccumuAverY[k][p] ,2);
        }
    }
    
    if(ionAccumuNumber[k][p]!=0)
    {
        x2Aver  = x2Aver  /ionAccumuNumber[k][p];
        y2Aver  = y2Aver  /ionAccumuNumber[k][p];
        
        ionAccumuRMSX[k][p] = sqrt(x2Aver);
        ionAccumuRMSY[k][p] = sqrt(y2Aver);
    }
    else
    {
        ionAccumuRMSX[k][p] = 0;
        ionAccumuRMSY[k][p] = 0;
    }
}

                                                 
void LatticeInterActionPoint::IonTransferDueToBunch(int bunchGap,int k, double bunchEffectiveSizeXMax, double bunchEffectiveSizeYMax)
{

    for(int p=0;p<gasSpec;p++)
    {
        for(int i=0;i<ionAccumuNumber[k][p];i++)
        {
            // kick----drift model to update ion velocity and position, time step is the bunchGap
            ionAccumuVelocityX[k][p][i] +=  ionAccumuFx[k][p][i];             // [m/s]
            ionAccumuVelocityY[k][p][i] +=  ionAccumuFy[k][p][i];

            ionAccumuPositionX[k][p][i] +=  ionAccumuVelocityX[k][p][i] * circRing/harmonics*bunchGap/CLight;
            ionAccumuPositionY[k][p][i] +=  ionAccumuVelocityY[k][p][i] * circRing/harmonics*bunchGap/CLight;
        }
    
    }

    // after tracking -- remove the lost ions here. 
    double ionLossXBoundary;
    double ionLossYBoundary;
    
    pipeAperatureX[k]>ionLossBoundary*bunchEffectiveSizeXMax ? (ionLossXBoundary = ionLossBoundary*bunchEffectiveSizeXMax) : (ionLossXBoundary = pipeAperatureX[k]);
    pipeAperatureY[k]>ionLossBoundary*bunchEffectiveSizeYMax ? (ionLossYBoundary = ionLossBoundary*bunchEffectiveSizeYMax) : (ionLossYBoundary = pipeAperatureY[k]);     
    
    totIonChargeAtInterPoint[k]=0;
    totMacroIonsAtInterPoint[k]=0;

    for(int p=0;p<gasSpec;p++)
    {    
        int i=0;
        while(i<ionAccumuPositionX[k][p].size())       
        {
            if(abs(ionAccumuPositionX[k][p][i])>ionLossXBoundary ||  abs(ionAccumuPositionY[k][p][i])>ionLossYBoundary )  //ion loss ceriteria
            {
                ionAccumuPositionX[k][p].erase(ionAccumuPositionX[k][p].begin()+i);
                ionAccumuPositionY[k][p].erase(ionAccumuPositionY[k][p].begin()+i);
                ionAccumuVelocityX[k][p].erase(ionAccumuVelocityX[k][p].begin()+i);
                ionAccumuVelocityY[k][p].erase(ionAccumuVelocityY[k][p].begin()+i);
                ionAccumuFx[k][p].erase(ionAccumuFx[k][p].begin()+i);
                ionAccumuFy[k][p].erase(ionAccumuFy[k][p].begin()+i);
                i--;
            }
            i++;
        }
        ionAccumuNumber[k][p]  =  ionAccumuPositionX[k][p].size();         

        totMacroIonsAtInterPoint[k] +=  ionAccumuNumber[k][p];
        totIonChargeAtInterPoint[k] +=  ionAccumuNumber[k][p] * macroIonCharge[k][p]; 
    }     
    //   cout<<ionLossXBoundary <<" "<<ionLossXBoundary<<" "<<totIonChargeAtInterPoint[k]<<"   "<<totMacroIonsAtInterPoint[k]
    //     <<" "<<ionLossBoundary<< "  "<< bunchEffectiveSizeXMax<<"   "<<bunchEffectiveSizeYMax<<endl;

    GetTotIonCharge();

}


void LatticeInterActionPoint:: GetTotIonCharge()
{

    totMacroIons = 0;
    totIonCharge = 0;

    for(int k=0;k<numberOfInteraction;k++)
    {
        for(int p=0;p<gasSpec;p++)
        {
            totIonCharge +=  ionAccumuNumber[k][p] * macroIonCharge[k][p];
            totMacroIons +=  ionAccumuNumber[k][p];
        }
    }

}

