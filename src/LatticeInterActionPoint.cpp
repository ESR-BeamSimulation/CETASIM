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
#include <gsl/gsl_matrix.h>


LatticeInterActionPoint::LatticeInterActionPoint()
{
}


LatticeInterActionPoint::~LatticeInterActionPoint()
{
    delete resDrivingTerms;
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
    averBetaX.resize(numberOfInteraction);
    averBetaY.resize(numberOfInteraction);
    
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
    

	// set the symplectic tracking map with dispersion, Twiss and rotation matrix
	symplecticMapB1H1.resize(numberOfInteraction);
	symplecticMapInvH1InvB1.resize(numberOfInteraction);
    symplecticMapInvH2InvB2.resize(numberOfInteraction);
    phaseAdvX12.resize(numberOfInteraction);
    phaseAdvY12.resize(numberOfInteraction);
    phaseAdvZ12.resize(numberOfInteraction);
   
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
    allIonAccumuRMSX.resize(numberOfInteraction);
    allIonAccumuRMSY.resize(numberOfInteraction);
    allIonAccumuAverX.resize(numberOfInteraction);
    allIonAccumuAverY.resize(numberOfInteraction);

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
    InitialLatticeIonInfo(inputParameter);
    InitialLatticeSympMat(inputParameter); 
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


    vector<double> intLengthTemp = interactionLength;    
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
}


void LatticeInterActionPoint::InitialLatticeIonInfo(const ReadInputSettings &inputParameter)
{
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

void LatticeInterActionPoint::InitialLatticeSympMat(const ReadInputSettings &inputParameter)
{
	double alphaX1,betaX1,phiX1,etaX1,etaXp1,alphaY1,betaY1,phiY1,etaY1,etaYp1,betaZ1; 
	double alphaX2,betaX2,phiX2,etaX2,etaXp2,alphaY2,betaY2,phiY2,etaY2,etaYp2,betaZ2;
 	double phiX21,phiY21,phiZ21;
 		
 	for(int i=0;i<numberOfInteraction;i++)
 	{
 		symplecticMapB1H1[i].Mat2DCreate(6,6);
 		symplecticMapInvH1InvB1[i].Mat2DCreate(6,6);
 		symplecticMapInvH2InvB2[i].Mat2DCreate(6,6);	
 	}
 	
	for(int i=0;i<numberOfInteraction;i++)
	{           
		alphaX1 = twissAlphaX[i];
		betaX1  = twissBetaX[i];
		alphaY1 = twissAlphaY[i];
		betaY1  = twissBetaY[i];
		phiX1   = xPhaseAdv[i];
		phiY1   = yPhaseAdv[i];
		etaX1   = twissDispX[i];
		etaXp1  = twissDispPX[i];
		etaY1   = twissDispY[i];
		etaYp1  = twissDispPY[i];
		betaZ1  = twissBetaZ[i];
		

		if(i<numberOfInteraction-1)
		{
			alphaX2 = twissAlphaX[i+1];
			betaX2  = twissBetaX[i+1];
			alphaY2 = twissAlphaY[i+1];
			betaY2  = twissBetaY[i+1];
			phiX2 	= xPhaseAdv[i+1]; 			
			phiY2   = yPhaseAdv[i+1];
			etaX2   = twissDispX[i+1];
			etaXp2  = twissDispPX[i+1];
			etaY2   = twissDispY[i+1];
			etaYp2  = twissDispPY[i+1];
			betaZ2  = twissBetaZ[i+1];           
		}
		else
		{
			alphaX2 = twissAlphaX[0];
			betaX2  = twissBetaX[0];
			alphaY2 = twissAlphaY[0];
			betaY2  = twissBetaY[0];
			phiX2   = 2 * PI * inputParameter.ringParBasic->workQx; 			
			phiY2   = 2 * PI * inputParameter.ringParBasic->workQy;
			etaX2   = twissDispX[0];
			etaXp2  = twissDispPX[0];
			etaY2   = twissDispY[0];
			etaYp2  = twissDispPY[0];
			betaZ2  = twissBetaZ[0];           
		}
				
		phaseAdvX12[i] = phiX2 - phiX1;
		phaseAdvY12[i] = phiY2 - phiY1;
		
		gsl_matrix *matH1  	  = gsl_matrix_alloc (6, 6); 	gsl_matrix_set_identity(matH1);
		gsl_matrix *matInvH1  = gsl_matrix_alloc (6, 6); 	gsl_matrix_set_identity(matInvH1);
		gsl_matrix *matInvH2  = gsl_matrix_alloc (6, 6); 	gsl_matrix_set_identity(matInvH2);
		gsl_matrix *matB1  	  = gsl_matrix_alloc (6, 6); 	gsl_matrix_set_zero(matB1);
		gsl_matrix *matInvB1  = gsl_matrix_alloc (6, 6); 	gsl_matrix_set_zero(matInvB1);
		gsl_matrix *matInvB2  = gsl_matrix_alloc (6, 6); 	gsl_matrix_set_zero(matInvB2);

		
		// set H1 and InvH1 and InvH2		
		// set H1
		gsl_matrix_set(matH1,0,5,-etaX1);   gsl_matrix_set(matH1,1,5,-etaXp1); gsl_matrix_set(matH1,2,5,-etaY1);   gsl_matrix_set(matH1,3,5,-etaYp1); 	
		gsl_matrix_set(matH1,4,0, etaXp1);  gsl_matrix_set(matH1,4,1,-etaX1);  gsl_matrix_set(matH1,4,2, etaYp1);  gsl_matrix_set(matH1,4,3,-etaY1); 	
		
		// set InvH1 		
		gsl_matrix_set(matInvH1,0,5, etaX1);   gsl_matrix_set(matInvH1,1,5, etaXp1); gsl_matrix_set(matInvH1,2,5, etaY1);    gsl_matrix_set(matInvH1,3,5, etaYp1); 	
		gsl_matrix_set(matInvH1,4,0,-etaXp1);  gsl_matrix_set(matInvH1,4,1, etaX1);  gsl_matrix_set(matInvH1,4,2, -etaYp1);  gsl_matrix_set(matInvH1,4,3, etaY1); 
		
		// set InvH2 	
		gsl_matrix_set(matInvH2,0,5, etaX2);   gsl_matrix_set(matInvH2,1,5, etaXp2); gsl_matrix_set(matInvH2,2,5, etaY2);    gsl_matrix_set(matInvH2,3,5, etaYp2); 	
		gsl_matrix_set(matInvH2,4,0,-etaXp2);  gsl_matrix_set(matInvH2,4,1, etaX2);  gsl_matrix_set(matInvH2,4,2,-etaYp2);   gsl_matrix_set(matInvH2,4,3, etaY2); 	
		
		
		// set B1 and InvB1 and InvB2
		// set B1
		gsl_matrix_set(matB1,0,0, 1 / sqrt(betaX1)); gsl_matrix_set(matB1,1,0,alphaX1 / sqrt(betaX1) );  gsl_matrix_set(matB1,1,1,sqrt(betaX1));
		gsl_matrix_set(matB1,2,2, 1 / sqrt(betaY1)); gsl_matrix_set(matB1,3,2,alphaY1 / sqrt(betaY1) );  gsl_matrix_set(matB1,3,3,sqrt(betaY1));
		gsl_matrix_set(matB1,4,4, 1 / sqrt(betaZ1)); 													 gsl_matrix_set(matB1,5,5,sqrt(betaZ1));
		
		// set invB1
		gsl_matrix_set(matInvB1,0,0, sqrt(betaX1)); gsl_matrix_set(matInvB1,1,0,- alphaX1 / sqrt(betaX1)); gsl_matrix_set(matInvB1,1,1,1 / sqrt(betaX1));
		gsl_matrix_set(matInvB1,2,2, sqrt(betaY1)); gsl_matrix_set(matInvB1,3,2,- alphaY1 / sqrt(betaY1)); gsl_matrix_set(matInvB1,3,3,1 / sqrt(betaY1));
		gsl_matrix_set(matInvB1,4,4, sqrt(betaZ1)); 													   gsl_matrix_set(matInvB1,5,5,1 / sqrt(betaZ1));		
		
		// set invB2
		gsl_matrix_set(matInvB2,0,0, sqrt(betaX2)); gsl_matrix_set(matInvB2,1,0,- alphaX2 / sqrt(betaX2)); gsl_matrix_set(matInvB2,1,1,1 / sqrt(betaX2));
		gsl_matrix_set(matInvB2,2,2, sqrt(betaY2)); gsl_matrix_set(matInvB2,3,2,- alphaY2 / sqrt(betaY2)); gsl_matrix_set(matInvB2,3,3,1 / sqrt(betaY2));
		gsl_matrix_set(matInvB2,4,4, sqrt(betaZ2)); 													   gsl_matrix_set(matInvB2,5,5,1 / sqrt(betaZ2));	

		
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,matB1,   matH1,   0.0,symplecticMapB1H1[i].mat2D);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,matInvH1,matInvB1,0.0,symplecticMapInvH1InvB1[i].mat2D);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,matInvH2,matInvB2,0.0,symplecticMapInvH2InvB2[i].mat2D);
		
		

      	//PrintGSLMatrix(symplecticMapB1H1[i].mat2D);
      	//PrintGSLMatrix(symplecticMapInvH2InvB2[i].mat2D);
          	
      	gsl_matrix_free(matH1);
      	gsl_matrix_free(matInvH2);
      	gsl_matrix_free(matB1);
      	gsl_matrix_free(matInvB2);
      	gsl_matrix_free(matInvH1);
      	gsl_matrix_free(matInvB1);
      	    			                   
     }   
		
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

    gsl_matrix * R2  = gsl_matrix_alloc (2, 2);     // coupling matrix Eq.(6)
    gsl_matrix_set_zero(R2);                        // no coupling at the point in SR calculation
    // gsl_matrix_set(R2,0,0,tengRMatR2[0]);
    // gsl_matrix_set(R2,0,1,tengRMatR2[1]);
    // gsl_matrix_set(R2,1,0,tengRMatR2[2]);
    // gsl_matrix_set(R2,1,1,tengRMatR2[3]);

    double det = get_det(R2);
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

void LatticeInterActionPoint::GetTransLinearCouplingCoef(const ReadInputSettings &inputParameter)
{
    double alphax,betax,alphay,betay,gammax,gammay;
    alphax = twissAlphaX[0];
    alphay = twissAlphaY[0];
    betax  = twissBetaX[0];
    betay  = twissBetaY[0];
    gammax = (1 + pow(alphax,2))/betax;
    gammay = (1 + pow(alphay,2))/betay;
    
    double nux = inputParameter.ringParBasic->workQx;
    double nuy = inputParameter.ringParBasic->workQy; 

    double phix = 2 * PI * nux ;
    double phiy = 2 * PI * nuy ;

    gsl_matrix *skewQuad44  = gsl_matrix_alloc (4, 4); 
    gsl_matrix_set_identity(skewQuad44);
    gsl_matrix_set(skewQuad44,1,2,inputParameter.ringParBasic->skewQuadK);
    gsl_matrix_set(skewQuad44,3,0,inputParameter.ringParBasic->skewQuadK);

    gsl_matrix *ILmatrix44  = gsl_matrix_alloc (4, 4); 
    gsl_matrix_set_zero(ILmatrix44);


    double tmp;
    tmp = cos(phix) + alphax * sin(phix);   gsl_matrix_set(ILmatrix44,0,0,tmp);   //R11
    tmp =             betax  * sin(phix);   gsl_matrix_set(ILmatrix44,0,1,tmp);   //R12
    tmp =           - gammax * sin(phix);   gsl_matrix_set(ILmatrix44,1,0,tmp);   //R21
    tmp = cos(phix) - alphax * sin(phix);   gsl_matrix_set(ILmatrix44,1,1,tmp);   //R22

    tmp = cos(phiy) + alphay * sin(phiy);   gsl_matrix_set(ILmatrix44,2,2,tmp);   //R33
    tmp =             betay  * sin(phiy);   gsl_matrix_set(ILmatrix44,2,3,tmp);   //R34
    tmp =           - gammay * sin(phiy);   gsl_matrix_set(ILmatrix44,3,2,tmp);   //R43
    tmp = cos(phiy) - alphay * sin(phiy);   gsl_matrix_set(ILmatrix44,3,3,tmp);   //R44
    

    gsl_matrix *oneTurn44  = gsl_matrix_alloc (4, 4);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,ILmatrix44,skewQuad44,0.0,oneTurn44);

    gsl_matrix *MM  = gsl_matrix_alloc (2, 2);
    gsl_matrix *mm  = gsl_matrix_alloc (2, 2); 
    gsl_matrix *NN  = gsl_matrix_alloc (2, 2);
    gsl_matrix *nn  = gsl_matrix_alloc (2, 2); 

    gsl_matrix *J2 = gsl_matrix_alloc (2, 2);
    gsl_matrix_set_zero(J2);
    gsl_matrix_set (J2, 0, 1, 1);
    gsl_matrix_set (J2, 1, 0,-1);


    // Ref. PRAB,2,074001, 1999. initialize the coupling terms with matrix formation
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            gsl_matrix_set(MM,i,j,gsl_matrix_get(oneTurn44,i,  j  ));
            gsl_matrix_set(NN,i,j,gsl_matrix_get(oneTurn44,i+2,j+2));
            gsl_matrix_set(mm,i,j,gsl_matrix_get(oneTurn44,i,  j+2));
            gsl_matrix_set(nn,i,j,gsl_matrix_get(oneTurn44,i+2,j  ));
        }
    }

 
    // Ref. Eq. (10)
    gsl_matrix *nnp     = gsl_matrix_alloc (2, 2);
    gsl_matrix *temp22  = gsl_matrix_alloc (2, 2);
    gsl_matrix_set_zero(temp22);
    gsl_matrix_transpose_memcpy(nnp,nn);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1.0, J2,     nnp,  0.0, temp22);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, temp22,  J2,  0.0, nnp); 

    gsl_matrix *HH = gsl_matrix_alloc (2, 2);
    gsl_matrix_set_zero(HH);
    gsl_matrix_memcpy(HH,nnp);
    gsl_matrix_add(HH,mm);

    // calculate equation (8) and (9) and (17) (18)
    gsl_matrix *MMsubNN = gsl_matrix_alloc (2, 2);
    gsl_matrix_set_zero(MMsubNN);
    gsl_matrix_memcpy(MMsubNN,MM);
    gsl_matrix_sub(MMsubNN,NN);

    double trace = gsl_get_trace(MMsubNN);
    double det = get_det(HH);
    int sign = (trace < 0) ? -1 : 1;
   
    double temp1  = trace==0? 0 :  pow(trace,2) / (pow(trace,2) + 4 * det);
    
    double gamma0[2],coeff0[2];
    for(int i=0;i<2;i++) coeff0[i] = 0;  
    gamma0[0] = sqrt(0.5 + 0.5*sqrt(temp1));
    gamma0[1] = sqrt(0.5 - 0.5*sqrt(temp1));
    
    if(gamma0[0]!=0) coeff0[0] = - sign / gamma0[0] / sqrt(pow(trace,2) + 4 * det); 
    if(gamma0[1]!=0) coeff0[1] =   sign / gamma0[1] / sqrt(pow(trace,2) + 4 * det); 


    gsl_matrix *CC  = gsl_matrix_alloc (2, 2);
    gsl_matrix *CCp = gsl_matrix_alloc (2, 2);
  
    gsl_matrix *VV        = gsl_matrix_alloc (4, 4);
    gsl_matrix *VVI       = gsl_matrix_alloc (4, 4);
    gsl_matrix *UU        = gsl_matrix_alloc (4, 4);
    gsl_matrix *temp44    = gsl_matrix_alloc (4, 4);
    gsl_matrix *Ga        = gsl_matrix_alloc (2, 2);
    gsl_matrix *Gb        = gsl_matrix_alloc (2, 2);
    gsl_matrix *GbI       = gsl_matrix_alloc (2, 2);
    gsl_matrix *CCB       = gsl_matrix_alloc (2, 2);
    
    gsl_matrix_set_zero(CC);
    gsl_matrix_set_zero(CCp);
    gsl_matrix_set_zero(VV);
    gsl_matrix_set_zero(VVI);
    gsl_matrix_set_zero(UU);
    gsl_matrix_set_zero(temp44); 
    gsl_matrix_set_zero(Ga);
    gsl_matrix_set_zero(Gb); 
    gsl_matrix_set_zero(GbI);
    gsl_matrix_set_zero(CCB);
    
    
    double gamma;
    double  coef;  
    if(inputParameter.ringParBasic->skewQuadK==0)
    {
        gamma = gamma0[0];
        coef  = coeff0[0];
    }
    else
    {
        gamma = (nuy - floor(nuy)) -  (nux - floor(nux)) <=0? gamma0[0]:gamma0[0];
        coef  = (nuy - floor(nuy)) -  (nux - floor(nux)) <=0? coeff0[0]:coeff0[0];
    }
    
    // set C matrix
    gsl_matrix_memcpy(CC,HH);
    gsl_matrix_scale(CC,coef);
    gsl_matrix_transpose_memcpy(CCp,CC);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1.0, J2    , CCp, 0.0,temp22);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, temp22, J2,  0.0,CCp);
    // set V matrix
    gsl_matrix_set_identity(VV);
    gsl_matrix_scale(VV,gamma);

            
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
        {
            gsl_matrix_set(VV,i,   j+2,  gsl_matrix_get(CC, i,j));
            gsl_matrix_set(VV,i+2, j,   -gsl_matrix_get(CCp,i,j));
        }     
    }

    gsl_matrix_memcpy(VVI,VV);
    gsl_matrix_inv(VVI);
 
    // eq(2) to get UU
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, VVI,    oneTurn44, 0.0, temp44);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, temp44, VV,        0.0, UU);

    // CCB is Eq (50),   
    gsl_matrix_set(Ga,0,0,    1./sqrt(betax));
    gsl_matrix_set(Ga,1,0,alphax/sqrt(betax));
    gsl_matrix_set(Ga,1,1,       sqrt(betax));
    gsl_matrix_set(Gb,0,0,    1./sqrt(betay));
    gsl_matrix_set(Gb,1,0,alphay/sqrt(betay));
    gsl_matrix_set(Gb,1,1,       sqrt(betay));

    gsl_matrix_memcpy(GbI,Gb);
    gsl_matrix_inv(GbI);

    gsl_matrix_set_zero(temp22);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, Ga,     CC,        0.0, temp22);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, temp22, GbI,       0.0, CCB);

    
    // get the RDT terms from the CCB matrix notation. Ref. PRAB.8.034001
    resDrivingTerms-> f1001 = ( gsl_matrix_get(CCB,0,1) - gsl_matrix_get(CCB,1,0) + li * gsl_matrix_get(CCB,0,0) + li * gsl_matrix_get(CCB,1,1) ) / (4 * gamma);
    resDrivingTerms-> f1010 = (-gsl_matrix_get(CCB,0,1) - gsl_matrix_get(CCB,1,0) + li * gsl_matrix_get(CCB,0,0) - li * gsl_matrix_get(CCB,1,1) ) / (4 * gamma);
    resDrivingTerms-> f0110 = conj(resDrivingTerms-> f1001);
    resDrivingTerms-> f0101 = conj(resDrivingTerms-> f1010);
    resDrivingTerms->linearCouplingFactor = sqrt(betax * betay / 4 / pow(PI,2) ) * inputParameter.ringParBasic->skewQuadK;  // []
    
    traceAB[0] = gsl_matrix_get(UU,0,0) + gsl_matrix_get(UU,1,1); // it gives the eigen tune of the eigen mode
    traceAB[1] = gsl_matrix_get(UU,2,2) + gsl_matrix_get(UU,3,3); // it gives the eigen tune of the eigen mode
    traceAB[2] = trace;
    gammaC[0]  = gamma0[0];
    gammaC[1]  = gamma0[1]; 
    detH =  det;  

    // get xy-alpha PRAB 10 064003, Eq.(15)
    double deltaNu = nuy - floor(nuy) - (nux - floor(nux) );
    if(deltaNu==0)
    {
        resDrivingTerms->xyAlpha = resDrivingTerms->linearCouplingFactor>0 ? PI / 4: -PI / 4   ;
    }
    else
    {
        resDrivingTerms->xyAlpha = atan2(resDrivingTerms->linearCouplingFactor,deltaNu) / 2.0 ;
    }

    tengRMatR2[0] = gsl_matrix_get(CC,0,0);
    tengRMatR2[1] = gsl_matrix_get(CC,0,1);
    tengRMatR2[2] = gsl_matrix_get(CC,1,0);
    tengRMatR2[3] = gsl_matrix_get(CC,1,1);



    // cout<<resDrivingTerms->linearCouplingFactor<<"  "<<deltaNu<<endl;

    // equaiton(3)  PRAB 8, 034001, 2005 / get f1001 and f1010 according the defintion to compare
    // equaiton(41) PRAB 10, 064003, 2007, the same as above
    // double phase0 = 2* PI * (nux - nuy);
    // double phase1 = 2* PI * (nux + nuy); 
    // complex<double> coef0 = 1.0 / 4.0 / (1.0 -  exp(li * phase0));
    // complex<double> coef1 = 1.0 / 4.0 / (1.0 -  exp(li * phase1));
    // double kernel = sqrt(betax * betay) * inputParameter.ringParBasic->skewQuadK;
    // complex<double> f1001 = - coef0 * kernel;
    // complex<double> f1010 = - coef1 * kernel;
    // cout<<abs(f1001)<<" "<<arg(f1001)<<endl;
    // cout<<abs(f1010)<<" "<<arg(f1010)<<endl;
    // getchar();

    gsl_matrix_free(MM);
    gsl_matrix_free(mm);
    gsl_matrix_free(NN);
    gsl_matrix_free(nn);
    gsl_matrix_free(nnp); 
    gsl_matrix_free(temp22);
    gsl_matrix_free(J2); 
    gsl_matrix_free(HH);
    gsl_matrix_free(MMsubNN);
    gsl_matrix_free(CCp);
    gsl_matrix_free(CC);
    gsl_matrix_free(VV);
    gsl_matrix_free(UU);  
    gsl_matrix_free(VVI);  
    gsl_matrix_free(temp44); 
    gsl_matrix_free(Ga); 
    gsl_matrix_free(Gb);
    gsl_matrix_free(GbI);
    gsl_matrix_free(CCB);        

    gsl_matrix_free(oneTurn44);
    gsl_matrix_free(skewQuad44);
    gsl_matrix_free(ILmatrix44);
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

            if( pow( (tempx-xAver)/rmsRx, 2)  + pow( (tempy-yAver)/rmsRy, 2) > 4.E0  ) // ions generated within 3 rms beam size.
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

    int totIonNum = 0;
    double xSum = 0.E0;
    double ySum = 0.E0;

    for(int p=0; p<gasSpec;p++)
    {
        totIonNum += ionAccumuNumber[k][p] ;
        
        for(int i=0;i<ionAccumuNumber[k][p];i++)
        {
            xSum += ionAccumuPositionX[k][p][i];
            ySum += ionAccumuPositionY[k][p][i];
        }
    }

    if(totIonNum!=0)
    {
        allIonAccumuAverX[k] = xSum / totIonNum;
        allIonAccumuAverY[k] = ySum / totIonNum;
    }
    else
    {
        allIonAccumuAverX[k] = 0 ;
        allIonAccumuAverY[k] = 0 ;
    }

    double x2Sum = 0;
    double y2Sum = 0;

    for(int p=0; p<gasSpec;p++)
    {
        totIonNum += ionAccumuNumber[k][p] ;
        
        for(int i=0;i<ionAccumuNumber[k][p];i++)
        {
            x2Sum += pow(ionAccumuPositionX[k][p][i] - allIonAccumuAverX[k],2) ;
            y2Sum += pow(ionAccumuPositionY[k][p][i] - allIonAccumuAverY[k],2) ;
        }
    }

    if(totIonNum!=0)
    {
        allIonAccumuRMSX[k]  = sqrt( x2Sum / totIonNum ) ;
        allIonAccumuRMSY[k]  = sqrt( y2Sum / totIonNum ) ;
    }
    else
    {
        allIonAccumuRMSX[k] = 0;
        allIonAccumuRMSY[k] = 0;
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

