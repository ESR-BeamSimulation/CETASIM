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

