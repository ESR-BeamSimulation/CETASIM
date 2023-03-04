//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             

#include"ReadInputSettings.h"
#include"Global.h"


#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstdlib>
#include<algorithm>
#include<cmath>

using namespace std;
using std::vector;

ReadInputSettings::ReadInputSettings()
{
  
}
ReadInputSettings::~ReadInputSettings()
{
   
    delete ringParBasic;         
    delete ringParRf;
    delete ringFillPatt;
    delete ringBunchPara;
    delete ringIonEffPara;
    delete ringBBFB;
    delete ringLRWake;
    delete ringSRWake;
    delete ringBBImp;
    delete ringRun;
    delete driveMode; 
    
}

int ReadInputSettings::ParamRead(int argc, char *argv[])
{

    string inputFile;
    if(argc>=2)
    {
      inputFile=argv[1];
    }
    else
    {
      inputFile="input.dat";
    }

    ifstream fin(inputFile);
    vector<string> strVec;
    string         str;
    
    if (! fin.is_open())
    {
        cerr<< "Error opening file input.dat"<<endl;
        exit(0);
    }
    

    while (!fin.eof())
    {
        getline(fin,str);
        
        if(str.length()==0 || str.find("=")==-1 )  continue;
        strVec.clear();
        StringVecSplit(str, strVec);
        
        string inputKey=strVec[0];

        // 1) set the vlaue for struct RingParSet 		
        if(strVec[0]=="ringcircring")
        {
            ringParBasic->circRing = stod(strVec[1]);
        }                      
        if(strVec[0]=="ringworkq")
        {
            ringParBasic->workQx = stod(strVec[1]);
            ringParBasic->workQy = stod(strVec[2]);
            ringParBasic->workQz = stod(strVec[3]);
        }
        if(strVec[0]=="ringchrom")
        {
            ringParBasic->chrom[0] = stod(strVec[1]);
            ringParBasic->chrom[1] = stod(strVec[2]);
        }
           
        if(strVec[0]=="ringsynchraddamptime")
        {
            ringParBasic->synchRadDampTime[0]=stod(strVec[1]);
            ringParBasic->synchRadDampTime[1]=stod(strVec[2]); 
            ringParBasic->synchRadDampTime[2]=stod(strVec[3]);
        }         
        
        if(strVec[0]=="ringpipeaperaturex")
        {
            ringParBasic->pipeAperature[1] = stod(strVec[1]);            
        }        
        if(strVec[0]=="ringpipeaperaturey")
        {
            ringParBasic->pipeAperature[2] = stod(strVec[1]);
        }
        if(strVec[0]=="ringpipeaperaturer")
        {
          ringParBasic->pipeAperature[0] = stod(strVec[1]);
        }   
        if(strVec[0]=="ringelectronbeamenergy")
        {
          ringParBasic->electronBeamEnergy = stod(strVec[1]);
        }
        if(strVec[0]=="ringalphac")
        {
          ringParBasic->alphac[0] = stod(strVec[1]);
          ringParBasic->alphac[1] = stod(strVec[2]);
          ringParBasic->alphac[2] = stod(strVec[3]);
        }

        if(strVec[0]=="ringadtx")
        {
          ringParBasic->aDTX[0] = stod(strVec[1]);
          ringParBasic->aDTX[1] = stod(strVec[2]);
        }
        if(strVec[0]=="ringadty")
        {
          ringParBasic->aDTY[0] = stod(strVec[1]);
          ringParBasic->aDTY[1] = stod(strVec[2]);
        }
        if(strVec[0]=="ringadtxy")
        {
          ringParBasic->aDTXY[0] = stod(strVec[1]);
          ringParBasic->aDTXY[1] = stod(strVec[2]);
        }

        if(strVec[0]=="ringu0")
        {
            ringParBasic->u0 = stod(strVec[1]);
        }
        if(strVec[0]=="ringsdelta0")
        {
            ringParBasic->sdelta0 = stod(strVec[1]);
        }  
        if(strVec[0]=="ringcurrent")
        {
            ringParBasic->ringCurrent = stod(strVec[1]);
        }
        if(strVec[0]=="ringnaturalemit")
        {
            ringParBasic->naturalEmit = stod(strVec[1]);
        }
        if(strVec[0]=="ringcoupling")
        {
            ringParBasic->couplingfactorXY = stod(strVec[1]);
        }     
              
        //----------------------------------------------------------------
         
      
        // 2) set the vlaue for struct rf set 	          
        
        if(strVec[0]=="rfringharm")
        {
            ringParRf->ringHarm = stoi(strVec[1]);         
        }

        
        if (strVec[0]=="rfresnum") 
        {
            ringParRf->resNum = stoi(strVec[1]);
            
            ringParRf->resType        .resize(ringParRf->resNum);
            ringParRf->resHarm        .resize(ringParRf->resNum);
            ringParRf->resVol         .resize(ringParRf->resNum);
            ringParRf->resPhase       .resize(ringParRf->resNum);
            ringParRf->resShuntImpRs  .resize(ringParRf->resNum);
            ringParRf->resQualityQ0   .resize(ringParRf->resNum);
            ringParRf->resCouplingBeta.resize(ringParRf->resNum);
            ringParRf->resDetuneFre   .resize(ringParRf->resNum); 
            ringParRf->resCold        .resize(ringParRf->resNum);
            ringParRf->resExciteIntability.resize(ringParRf->resNum);
            ringParRf->resDirFB.resize(ringParRf->resNum);
            ringParRf->resDirFBGain.resize(ringParRf->resNum);
            ringParRf->resDirFBPhase.resize(ringParRf->resNum);
            ringParRf->resDirFBDelay.resize(ringParRf->resNum);               
            ringParRf->resAmpFBRatioForTotSelfLoss.resize(ringParRf->resNum);                                 
        }

        if (strVec[0]=="rfresexciteintability") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resExciteIntability[i] = stod(strVec[i+1]);               
            }                    
        }

        if (strVec[0]=="rfresdirfb") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resDirFB[i] = stod(strVec[i+1]);               
            }                    
        }

        if (strVec[0]=="rfresdirfbgain") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resDirFBGain[i] = stod(strVec[i+1]);               
            }                    
        }

        if (strVec[0]=="rfresdirfbphaseshift") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resDirFBPhase[i] = stod(strVec[i+1]);               
            }                    
        }
        
        if (strVec[0]=="rfresdirfbdelay") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resDirFBDelay[i] = stod(strVec[i+1]);               
            }                    
        }


        if (strVec[0]=="rfresampfbratiofortotselfloss") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resAmpFBRatioForTotSelfLoss[i] = stod(strVec[i+1]);               
            }
                       
        }      
        if (strVec[0]=="rfrestype") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resType[i] = stoi(strVec[i+1]);
                if (ringParRf->resType[0]==0)
                {
                  cerr<<"set main cavity as active, rfrestype=1 "<<endl;
                }
            }
        }
        if (strVec[0]=="rfresharm") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resHarm[i] = stoi(strVec[i+1]);
            }
        }
        
        if (strVec[0]=="rfresshuntimprs") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resShuntImpRs[i] = stod(strVec[i+1]);
            }
        }
        
        if (strVec[0]=="rfresqualityq0") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resQualityQ0[i] = stod(strVec[i+1]);
            }
        }       
     
        if (strVec[0]=="rfrescouplingbeta") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resCouplingBeta[i] = stod(strVec[i+1]);
            }
        }
                         
        if (strVec[0]=="rfresdetunefre") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resDetuneFre[i] = stod(strVec[i+1]);
            }
        }                
               
        if (strVec[0]=="rfresvol") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resVol[i] = stod(strVec[i+1]);
                
            }
        }
        
        if (strVec[0]=="rfresphase") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resPhase[i] = stod(strVec[i+1]);
            }
        }
        
                    
        if (strVec[0]=="rftransientcavparwriteto") 
        {
            ringParRf->transResonParWriteTo = strVec[1];
        }

        if (strVec[0]=="rfbunchbinnum") 
        {
            ringParRf->rfBunchBinNum = stod(strVec[1]);
        }
        if (strVec[0]=="rfmethodforvb") 
        {
            ringParRf->methodForVb = strVec[1];
            ringParRf->methodForVb = "soft";
        }    
                 
        if (strVec[0]=="rfrescold") 
        {
            for(int i=0;i<ringParRf->resNum;i++)
            {
                ringParRf->resCold[i] = stoi(strVec[i+1]);
            }        
        } 
        //---------------------------------------------------------------------- 

        
        // 3) set the vlaue of different filling pattern 
         
        if(strVec[0]=="fptotbunchnumber")
        {
          ringFillPatt->totBunchNumber = stoi(strVec[1]);
        }        
        if(strVec[0]=="fptrainnumber")
        {
          ringFillPatt->trainNumber = stoi(strVec[1]);
          ringFillPatt->bunchNumberPerTrain.resize(ringFillPatt->trainNumber);
          ringFillPatt->trainGaps          .resize(ringFillPatt->trainNumber);
          ringFillPatt->bunchGaps          .resize(ringFillPatt->trainNumber);          
        }
        
        if(strVec[0]=="fpbunchnumberpertrain")
        {
          for(int i=0;i<ringFillPatt->trainNumber;i++)  
          {
            ringFillPatt->bunchNumberPerTrain[i] = stoi(strVec[i+1]);
          }  
        }               
        
        if(strVec[0]=="fptraingaps")
        {
          for(int i=0;i<ringFillPatt->trainNumber;i++)
          { 
            ringFillPatt->trainGaps[i] = stoi(strVec[i+1]);
          }
        }        
        if(strVec[0]=="fpbunchgaps")
        {
          for(int i=0;i<ringFillPatt->trainNumber;i++)
          { 
            ringFillPatt->bunchGaps[i] = stoi(strVec[i+1]);
          }
        }
        if(strVec[0]=="fpsetbunchchargenum")
        {
          ringFillPatt->bunchChargeNum =  stoi(strVec[1]);
          ringFillPatt->bunchChargeIndex.resize(ringFillPatt->bunchChargeNum );
          ringFillPatt->bunchCharge.resize(ringFillPatt->bunchChargeNum);
        }
        if(strVec[0]=="fpsetbunchchargeindex")
        {
          for(int i=0;i<ringFillPatt->bunchChargeNum;i++)
          { 
            ringFillPatt->bunchChargeIndex[i] = stoi(strVec[i+1]);
          }
        }
        if(strVec[0]=="fpsetbunchcharge")
        {
          for(int i=0;i<ringFillPatt->bunchChargeNum;i++)
          { 
            ringFillPatt->bunchCharge[i] = stoi(strVec[i+1]);
          }
        }

        //---------------------------------------------------------------------- 


        //4) initial bunch parameters
        if(strVec[0]=="bunchcurrent")
        {
          ringBunchPara->current = stod(strVec[1]);
        }
        if(strVec[0]=="bunchemittance")
        {
          ringBunchPara->emittanceX = stod(strVec[1]);
        }
        if(strVec[0]=="bunchkappa")
        {
          ringBunchPara->kappa = stod(strVec[1]);
        }
        
        if(strVec[0]=="bunchrmsbunchlength")
        {
          ringBunchPara->rmsBunchLength = stod(strVec[1]);
        }
        
        if(strVec[0]=="bunchrmsenergyspread")
        {
         ringBunchPara-> rmsEnergySpread = stod(strVec[1]);
        }
        
        if(strVec[0]=="bunchdistributiontype")
        {
          ringBunchPara->distributionType = stod(strVec[1]);
        }
        if(strVec[0]=="bunchmacroelenumperbunch")
        {
          ringBunchPara->macroEleNumPerBunch = stoi(strVec[1]);
        } 
        
             
        if(strVec[0]=="bunchinitialstaticoffset")
        {
          for(int i=0;i<=5;i++)  
          {
            ringBunchPara->initialStaticOffSet[i] = stod(strVec[i+1]);                 
          }            
        }      
        
        if(strVec[0]=="bunchinitialdynamicoffset")
        {
          for(int i=0;i<=5;i++)  
          {
            ringBunchPara->initialDynamicOffSet[i] = stod(strVec[i+1]);               
          }            
        }      
         
        if(strVec[0]=="bunchdiswriteto")
        {
          ringBunchPara->bunchDisWriteTo = strVec[1];
        }  
        
        //----------------------------------------------------------------------  
                
        //5 initial ion efect para        
        if(strVec[0]=="ioncalgasspec")
        {         
          ringIonEffPara->gasSpec = stod(strVec[1]);
          ringIonEffPara->ionMassNumber.resize(ringIonEffPara->gasSpec);
          ringIonEffPara->corssSectionEI.resize(ringIonEffPara->gasSpec);
          ringIonEffPara->gasPercent.resize(ringIonEffPara->gasSpec);          
        }
        if(strVec[0]=="ioncalionmassnumber")
        {                   
          for(int i=0;i<ringIonEffPara->gasSpec;i++)
          {
            ringIonEffPara->ionMassNumber[i] = stod(strVec[i+1]);
          }
        }            
        if(strVec[0]=="ioncalcorsssectionei")
        {
          for(int i=0;i<ringIonEffPara->gasSpec;i++)
          {          
            ringIonEffPara->corssSectionEI[i] = stod(strVec[i+1]) * 1.e-22;
          }  
        }
        if(strVec[0]=="ioncaliongaspercent")
        {
          for(int i=0;i<ringIonEffPara->gasSpec;i++)
          {          
            ringIonEffPara->gasPercent[i] = stod(strVec[i+1]);
          }  
        }


        
        if(strVec[0]=="ioncalionmaxnumber")
        {
          ringIonEffPara->ionMaxNumber = stod(strVec[1]);
        }   
        if(strVec[0]=="ioncalionlossboundary")
        {
          ringIonEffPara->ionLossBoundary = stod(strVec[1]);
        }    
        
        if(strVec[0]=="ioncalnumberofionbeaminterpoint")
        {
          ringIonEffPara->numberofIonBeamInterPoint = stoi(strVec[1]);
        } 
        
        if(strVec[0]=="ioncalmacroionnumbergeneratedperip")
        {   
           ringIonEffPara->macroIonNumberGeneratedPerIP = stod(strVec[1]);
        }
        if(strVec[0]=="ioncalioninfoprintinterval")
        {
          ringIonEffPara->ionInfoPrintInterval = stod(strVec[1]);
        }   
        if(strVec[0]=="ioncaliondiswriteto")
        {
          ringIonEffPara->ionDisWriteTo = strVec[1];
        }   
        if(strVec[0]=="ioncaltwissinput")
        {
            ringIonEffPara->twissInput = strVec[1];
        } 
        //----------------------------------------------------------------------  

            
        //6 initial bunch-by-bunch FB  para 

        if(strVec[0]=="fbkickstrengthkx")
        {
          ringBBFB->kickStrengthK[0] = stod(strVec[1]);
        }  
        if(strVec[0]=="fbkickstrengthky")
        {
          ringBBFB->kickStrengthK[1] = stod(strVec[1]);
        } 
        if(strVec[0]=="fbkickstrengthf")
        {
          ringBBFB->kickStrengthK[2] = stod(strVec[1]);
        }           
        if(strVec[0]=="fbgain")
        {
          ringBBFB->gain = stod(strVec[1]);
        }  
        if(strVec[0]=="fbdelay")
        {
          ringBBFB->delay = stoi(strVec[1]);
        }  
        if(strVec[0]=="fbtaps")
        {
          ringBBFB->taps = stoi(strVec[1]);
          ringBBFB->firCoeffx.resize(ringBBFB->taps + ringBBFB->delay);
          ringBBFB->firCoeffy.resize(ringBBFB->taps + ringBBFB->delay);
          ringBBFB->firCoeffz.resize(ringBBFB->taps + ringBBFB->delay);
          ringBBFB->firCoeffxy.resize(ringBBFB->taps + ringBBFB->delay);          
        }  
        if(strVec[0]=="fbkickerdispp")
        {
          ringBBFB->kickerDispP = stod(strVec[1]);
        }  
        if(strVec[0]=="fbkickerdisp")
        {
          ringBBFB->kickerDisp = stod(strVec[1]);
        }  
        if(strVec[0]=="fbfirbunchbybunchfeedbackpowerlimit")
        {
          ringBBFB->fIRBunchByBunchFeedbackPowerLimit = stod(strVec[1]);
        }  
        if(strVec[0]=="fbfirbunchbybunchfeedbackkickerimped")
        {
          ringBBFB->fIRBunchByBunchFeedbackKickerImped = stod(strVec[1]);
        }  
        
        if(strVec[0]=="fbfircoeffx")
        {
          for(int i=0;i<(ringBBFB->delay+ringBBFB->taps); i++)
          {
            if(i<ringBBFB->delay)
            {
                ringBBFB->firCoeffx[i] = 0;
            }
            else
            {
                ringBBFB->firCoeffx[i] = stod(strVec[i-ringBBFB->delay+1]);
            }
          }
        }  
        if(strVec[0]=="fbfircoeffy")
        {
          for(int i=0;i<(ringBBFB->delay+ringBBFB->taps); i++)
          {
            if(i<ringBBFB->delay)
            {
                ringBBFB->firCoeffy[i] = 0;
            }
            else
            {
                ringBBFB->firCoeffy[i] = stod(strVec[i-ringBBFB->delay+1]);
            }
          }  
        }          
        if(strVec[0]=="fbfircoeffz")
        {
          for(int i=0;i<(ringBBFB->delay+ringBBFB->taps); i++)
          {
            if(i<ringBBFB->delay)
            {
                ringBBFB->firCoeffz[i] = 0;
            }
            else
            {
                ringBBFB->firCoeffz[i] = stod(strVec[i-ringBBFB->delay+1]);
            }
            
          }  
        }          
        if(strVec[0]=="fbfircoeffxy")
        {
          for(int i=0;i<(ringBBFB->delay+ringBBFB->taps); i++)
          {
            if(i<ringBBFB->delay)
            {
                ringBBFB->firCoeffxy[i] = 0;
            }
            else
            {
                ringBBFB->firCoeffxy[i] = stod(strVec[i-ringBBFB->delay+1]);
            }
          }  
        }                 

        // 7.1)  initial long range wake 
        if(strVec[0]=="lrwpipegeoinput")
        {
          ringLRWake->pipeGeoInput = strVec[1];          
        }
        if(strVec[0]=="lrwbbrinput")
        {
          ringLRWake->bbrInput = strVec[1];          
        }
        if(strVec[0]=="lrwnturnswaketrunction")
        {
          ringLRWake->nTurnswakeTrunction = stoi(strVec[1]) + 1;
        }
        
        // 7.1)  initial short range wake 
        if(strVec[0]=="srwpipegeoinput")
        {
          ringSRWake->pipeGeoInput = strVec[1];      
        }
        if(strVec[0]=="srwbbrinput")
        {
          ringSRWake->bbrInput = strVec[1];       
        }
        if(strVec[0]=="srwbunchbinnum")
        {
          ringSRWake->SRWBunchBinNum = stod(strVec[1]);      
        }
        
        if(strVec[0]=="srwwakepotenwriteto")
        {
          ringSRWake->SRWWakePotenWriteTo = strVec[1];       
        }
        
                                         
        // 8) impedance 
        if(strVec[0]=="bbiimpedinput")
        {
          ringBBImp->impedInput = strVec[1];
        } 
        if(strVec[0]=="bbibunchbinnumberz")
        {
          ringBBImp->bunchBinNumberZ = stoi(strVec[1]);
        } 
        
        // 9) driveMode
        if(strVec[0]=="drivemodeon")
        {
          driveMode->driveModeOn = stoi(strVec[1]);
        } 
        if(strVec[0]=="drivecbmindex")
        {
          driveMode->driveCBMIndex = stoi(strVec[1]);
        }
        if(strVec[0]=="driveend")
        {
          driveMode->driveEnd = stoi(strVec[1]);
        }
        if(strVec[0]=="drivestart")
        {
          driveMode->driveStart = stoi(strVec[1]);
        }
        if(strVec[0]=="driveamp")
        {
          driveMode->driveAmp = stod(strVec[1]);
        }
        if(strVec[0]=="drivefre")
        {
          driveMode->driveFre = stod(strVec[1]);
        }
        if(strVec[0]=="driveplane")
        {
          driveMode->drivePlane = stoi(strVec[1]);
        }
        if(strVec[0]=="drivehw")
        {
          driveMode->driveHW = stoi(strVec[1]);
        }
        
        // 10) run       
        
        if(strVec[0]=="runnturns")
        {
          ringRun->nTurns = stoi(strVec[1]);
        }
        
        if(strVec[0]=="runcalsetting")
        {
           ringRun->calSetting = stoi(strVec[1]);
        }
        if(strVec[0]=="runsynraddampingflag")
        {
           ringRun->synRadDampingFlag[0] = stoi(strVec[1]);
           ringRun->synRadDampingFlag[1] = stoi(strVec[2]);
        }

        if(strVec[0]=="runfirbunchbybunchfeedbackflag")
        {
           ringRun->fIRBunchByBunchFeedbackFlag = stoi(strVec[1]);
        }
              
        if(strVec[0]=="runbbiflag")
        {
           ringRun->bBImpFlag = stoi(strVec[1]);
        }          
        if(strVec[0]=="runlongrangewake")
        {
           ringRun->lRWakeFlag = stoi(strVec[1]);
        }
        if(strVec[0]=="runshortrangewake")
        {
           ringRun->sRWakeFlag = stoi(strVec[1]);
        }
        if(strVec[0]=="runtbtbunchaverdata")
        {
           ringRun->TBTBunchAverData = strVec[1];
        }
        if(strVec[0]=="runtbtbunchdisdata")
        {
           ringRun->TBTBunchDisData = strVec[1];
        }
        if(strVec[0]=="runtbtbunchcavvoldata")
        {
           ringRun->TBTBunchCavVolData = strVec[1];
        }
        if(strVec[0]=="runtbtbunchhaissinski")
        {
           ringRun->TBTBunchHaissinski = strVec[1];
        }
        if(strVec[0]=="runtbtbunchlongtraj")
        {
           ringRun->TBTBunchLongTraj = strVec[1];
        }

        if(strVec[0]=="runbeamionflag")
        {
           ringRun->beamIonFlag = stoi(strVec[1]);
        }
        if(strVec[0]=="runbunchinfoprintinterval")
        {
          ringRun->bunchInfoPrintInterval = stoi(strVec[1]);
        }
        if(strVec[0]=="rungrowthratefittingstart")
        {
          ringRun->growthRateFittingStart = stoi(strVec[1]);
        }
        if(strVec[0]=="rungrowthratefittingend")
        {
          ringRun->growthRateFittingEnd = stoi(strVec[1]);
        }
        
        
        if(strVec[0]=="runtbtbunchprintnum")
        {                       
          ringRun->TBTBunchPrintNum = stoi(strVec[1]) ;
          ringRun->TBTBunchDisDataBunchIndex.resize(ringRun->TBTBunchPrintNum);
        }

        if(strVec[0]=="runcbmgr")
        {                       
          ringRun->runCBMGR = strVec[1];
        }
            
        if(strVec[0]=="runtbtbunchdisdatabunchindex")
        {                                 
          for(int i=0;i<ringRun->TBTBunchDisDataBunchIndex.size();i++)
          {
            ringRun->TBTBunchDisDataBunchIndex[i] = stoi(strVec[i+1]);           
          }          
        }
                                                               
    }


  if (driveMode->driveStart > driveMode->driveStart)
  {
    cerr<<"wrong setting in DRIVEMode, driveMode->driveStart have to be less than driveMode->driveEnd"<<endl;
  }


  if(ringRun->growthRateFittingStart >= ringRun->growthRateFittingEnd  )
  {
    cerr<<"wrong settings: growthRateFittingEnd is smaller growthRateFittingStart"<<endl;
    exit(0);
  }

  if(ringRun->growthRateFittingEnd >  ringRun->nTurns )
  {
    cerr<<"wrong settings: nTurns small than growthRateFittingEnd"<<endl;
    exit(0);
  } 



    // debug -- print all bunch data
    // ringRun->TBTBunchPrintNum = ringFillPatt->totBunchNumber;
    // ringRun->TBTBunchDisDataBunchIndex.resize(ringRun->TBTBunchPrintNum );
    // for(int i=0;i<ringRun->TBTBunchPrintNum;i++)
    // {
    //    ringRun->TBTBunchDisDataBunchIndex[i] = i; 
    // }
    


    if(ringRun->TBTBunchPrintNum > ringFillPatt->totBunchNumber )
    {
      cerr<<"wrong settings: Total Bunch Number is "<<ringFillPatt->totBunchNumber<<endl;
      cerr<<ringRun->TBTBunchPrintNum  <<" bunches are to be printed out"<<endl;
      exit(0);
    }
    else
    {
      for(int i=0;i<ringRun->TBTBunchPrintNum;i++)
      {
        if(ringRun->TBTBunchDisDataBunchIndex[i] >= ringFillPatt->totBunchNumber)
        {
          cerr<<"wrong settings: runTBTBunchDisDataBunchIndex, the printed bunch index is larger than total bunch number"<<endl;  
          exit(0);
        }
      }
    }
    



    //1) data from ringPaBbasic
    double  circRing           = ringParBasic->circRing;
    double  sdelta0            = ringParBasic->sdelta0;
    double  *alphac            = ringParBasic->alphac;  
    double  electronBeamEnergy = ringParBasic->electronBeamEnergy;
    double  u0                 = ringParBasic->u0;
    double  rGamma             = electronBeamEnergy/ElectronMassEV;
    double  eta                = alphac[0] - 1./pow(rGamma,2);
    double  rBeta              = sqrt(1.E0 -pow(rGamma,-2) );
    double  t0                 = circRing / CLight / rBeta;
    double  f0                 = 1 / t0;
    double  omega0             = 2 * PI * f0;

 

    
    //2) data from ringParRF    
    int harmonics = ringParRf->ringHarm;

    
    //4.1) data from ringBunchPara
    double kappa               = ringBunchPara->kappa;
    double emittanceX          = ringBunchPara->emittanceX  / (1 + kappa);
    double emittanceY          = ringBunchPara->emittanceX  / (1 + kappa) * kappa;    
    ringBunchPara->emittanceY = emittanceY;
    ringBunchPara->emittanceX = emittanceX;

    double emittanceZ          = ringBunchPara->rmsEnergySpread * ringBunchPara->rmsBunchLength;
    ringBunchPara-> emittanceZ = emittanceZ;
 
 
 
   //1) set ringParBasic
    double vRF = ringParRf->resVol[0];
    double synPhase = PI - asin(u0/vRF);
        
    
    double workQz= sqrt(-1 * eta * harmonics * vRF * cos(synPhase) / (2 * PI * rBeta * electronBeamEnergy) );
    double sigmaT0 = eta * sdelta0 / (2 * PI * workQz * f0);
    
        
    ringParBasic->rGamma = rGamma;
    ringParBasic->rBeta = rBeta;
 
    ringParBasic->t0 = t0;
    ringParBasic->f0 = f0;
    
    ringParBasic->omega0 = omega0;  
    ringParBasic->kappa = kappa;
    ringParBasic->harmonics = harmonics;
    ringParBasic->sigmaT0 = sigmaT0;
    ringParBasic->workQz = workQz;
 

    ringParBasic->eta = eta;
    // ringParBasic->ringCurrent = ringBunchPara->current * ringFillPatt->totBunchNumber;
    ringParBasic->betaFunAver[0] = circRing / 2/ PI / ringParBasic->workQx ;
    ringParBasic->betaFunAver[1] = circRing / 2/ PI / ringParBasic->workQy ;
    ringParBasic->betaFunAver[2] = circRing / 2/ PI / ringParBasic->workQz ;
                
 
    //set the damping partition number, natural bunche length here
      
    double coupling = ringParBasic->couplingfactorXY;
  
    double alpha0 = u0 / (2 * t0 * electronBeamEnergy);

    ringParBasic->dampingPartJ[0] = 2 * electronBeamEnergy / u0 / ringParBasic->synchRadDampTime[0];
    ringParBasic->dampingPartJ[1] = 1;
    ringParBasic->dampingPartJ[2] = 3 - ringParBasic->dampingPartJ[0];
    
    double alpha[3];
    for (int i=0;i<3;i++) alpha[i] = ringParBasic->dampingPartJ[i] * alpha0;


    ringParBasic->radIntegral[0] = circRing * alphac[0];
    ringParBasic->radIntegral[1] = circRing * (alpha[0] +     alpha[2]) / (ElecClassicRadius * CLight * pow(rGamma,3));
    ringParBasic->radIntegral[3] = circRing * (alpha[2] - 2 * alpha[0]) / (ElecClassicRadius * CLight * pow(rGamma,3));
    ringParBasic->radIntegral[2] = (2 * ringParBasic->radIntegral[1] + ringParBasic->radIntegral[3] ) * pow(sdelta0,2) / (Cq * pow(rGamma,2));
    ringParBasic->radIntegral[4] = ringParBasic->emitNat[0] *  (ringParBasic->radIntegral[1] - ringParBasic->radIntegral[3] ) / (Cq * pow(rGamma,2));


    ringParBasic->naturalBunchLength =  eta * CLight * sdelta0 / (2 * PI * workQz * f0 );
    ringParBasic->emitNat[0] =  ringParBasic->naturalEmit / (1 + coupling);
    ringParBasic->emitNat[1] =  ringParBasic->naturalEmit / (1 + coupling) * coupling;
    ringParBasic->emitNat[2] =  ringParBasic->naturalBunchLength * sdelta0;


    fin.close();    

    return 0;
}



