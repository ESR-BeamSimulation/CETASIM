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
    delete ringImpedance;
    delete ringRun; 
    
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
        exit (1);
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
          ringParBasic->alphac = stod(strVec[1]);
        }
        if(strVec[0]=="ringu0")
        {
            ringParBasic->u0 = stod(strVec[1]);
        }
        if(strVec[0]=="ringsdelta0")
        {
            ringParBasic->sdelta0 = stod(strVec[1]);
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
            ringParRf->resAmpFBRatioForTotSelfLoss.resize(ringParRf->resNum);                                 
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
        //---------------------------------------------------------------------- 


        //4) initial bunch parameters
        if(strVec[0]=="bunchcurrent")
        {
          ringBunchPara->current = stod(strVec[1]);
        }
        if(strVec[0]=="bunchemittancex")
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

            
        //6 initial FB efect para 

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
          ringImpedance->impedInput = strVec[1];
        } 
        if(strVec[0]=="bbibunchbinnumberz")
        {
          ringImpedance->bunchBinNumberZ = stoi(strVec[1]);
        } 
        

        // 9) run       
        
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
           ringRun->synRadDampingFlag = stoi(strVec[1]);
        }

        if(strVec[0]=="runfirbunchbybunchfeedbackflag")
        {
           ringRun->fIRBunchByBunchFeedbackFlag = stoi(strVec[1]);
        }
              
        if(strVec[0]=="runbbiflag")
        {
           ringRun->impedanceFlag = stoi(strVec[1]);
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

        if(strVec[0]=="runtbtbunchprintnum")
        {                       
          ringRun->TBTBunchPrintNum = stoi(strVec[1]) ;
          ringRun->TBTBunchDisDataBunchIndex.resize(ringRun->TBTBunchPrintNum);
        }
    
        if(strVec[0]=="runtbtbunchdisdatabunchindex")
        {                                 
          for(int i=0;i<ringRun->TBTBunchDisDataBunchIndex.size();i++)
          {
            ringRun->TBTBunchDisDataBunchIndex[i] = stoi(strVec[i+1]);           
          }          
        }
                                                               
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
      cerr<<"wrong settings: Tottal Bunch Number is "<<ringFillPatt->totBunchNumber<<endl;
      cerr<<ringFillPatt->totBunchNumber <<" bunches are to be printed out"<<endl;
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
    double  alphac             = ringParBasic->alphac;
    double  electronBeamEnergy = ringParBasic->electronBeamEnergy;
    double  u0                 = ringParBasic->u0;
    double  rGamma             = electronBeamEnergy/ElectronMassEV;
    double  eta                = alphac - 1./pow(rGamma,2);
    double  rBeta              = sqrt(1.E0 -pow(rGamma,-2) );
    double  t0                 = circRing /  rBeta /  CLight;
    double  f0                 = 1 / t0;
    double  omega0             = 2 * PI * f0;

    
    //2) data from ringParRF    
    int harmonics = ringParRf->ringHarm;

    
    //4.1) data from ringBunchPara
    double kappa               = ringBunchPara->kappa;
    double emittanceX          = ringBunchPara->emittanceX;
    double emittanceZ          = ringBunchPara->rmsEnergySpread * ringBunchPara->rmsBunchLength;
    ringBunchPara-> emittanceZ = emittanceZ;
 
 
 
   //1) set ringParBasic
    double vRF = ringParRf->resVol[0];
    double synPhase = PI - asin(vRF/electronBeamEnergy);
        
    
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
    ringParBasic->current = ringBunchPara->current * ringFillPatt->totBunchNumber;
    ringParBasic->betaFunAver[0] = circRing / 2/ PI / ringParBasic->workQx ;
    ringParBasic->betaFunAver[1] = circRing / 2/ PI / ringParBasic->workQy ;
    ringParBasic->betaFunAver[2] = circRing / 2/ PI / ringParBasic->workQz ;
        
         
    double emittanceY = emittanceX * kappa;    
    ringBunchPara->emittanceY = emittanceY;
                  
          
    fin.close();    

    return 0;
}



