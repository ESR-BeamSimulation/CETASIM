#pragma once

#include"ReadInputSettings.h"

#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstdlib>
#include"Global.h"
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
    delete ringWake;
    delete ringImpedance;
    delete ringRun; 
    
}

int ReadInputSettings::ParamRead()
{

    /*
    ringParBasic = (RingParBasic*)      (malloc(sizeof(ringParBasic)));
    ringParRf    = (RingParRf*)         (malloc(sizeof(RingParRf)));
    ringFillPatt = (RingFillPatt*)      (malloc(sizeof(RingFillPatt)));
    ringBunchPara= (RingBunchPara*)     (malloc(sizeof(RingBunchPara)));
    ringIonEffPara=(RingIonEffPara*)    (malloc(sizeof(RingIonEffPara)));
    */



    ifstream fin("input.dat");
  
    vector<string> strVec;
    string         str;
    
    if (! fin.is_open())
    {
        cerr<< "Error opening file input.dat"; exit (1);
    }


    while (!fin.eof())
    {
        getline(fin,str);
        
        if(str.length()==0 || str.find("=")==-1 )  continue;
        strVec.clear();
        StringVecSplit(str, strVec);
        
        // 1) set the vlaue for struct RingParSet 		
        if(strVec[0]=="circring")
        {
            ringParBasic->circRing = stod(strVec[1]);
        }                      
        if(strVec[0]=="workq")
        {
            ringParBasic->workQx = stod(strVec[1]);
            ringParBasic->workQy = stod(strVec[2]);
            ringParBasic->workQz = stod(strVec[3]);       
        }
            
        if(strVec[0]=="synchraddamptime")
        {
            ringParBasic->synchRadDampTime[0]=stod(strVec[1]);
            ringParBasic->synchRadDampTime[1]=stod(strVec[2]); 
            ringParBasic->synchRadDampTime[2]=stod(strVec[3]);           
        }         
        
        if(strVec[0]=="pipeaperaturex")
        {
            ringParBasic->pipeAperature[1] = stod(strVec[1]);
        }        
        if(strVec[0]=="pipeaperaturey")
        {
            ringParBasic->pipeAperature[2] = stod(strVec[1]);
        }
        if(strVec[0]=="pipeaperaturer")
        {
          ringParBasic->pipeAperature[0] = stod(strVec[1]);
        }   
        if(strVec[0]=="electronbeamenergy")
        {
          ringParBasic->electronBeamEnergy = stod(strVec[1]);
        }
        if(strVec[0]=="alphac")
        {
          ringParBasic->alphac = stod(strVec[1]);
        }
        if(strVec[0]=="u0")
        {
            ringParBasic->u0 = stod(strVec[1]);
        }
        if(strVec[0]=="sdelta0")
        {
            ringParBasic->sdelta0 = stod(strVec[1]);
        }         

         
      
        // 2) set the vlaue for struct rf set 	          
        if(strVec[0]=="rfbasefrequency")
        {
            ringParRf->rfBaseFreq = stod(strVec[1]);           
        }
        if(strVec[0]=="vrf")
        {
            ringParRf->vRF = stod(strVec[1]);           
        }
        
        if(strVec[0]=="harmrf")
        {
            ringParRf->harmRF = stoi(strVec[1]);           
        }
        
        
        
        // 3) set the vlaue for simple filling pattern -- update it in the future....
         
        if(strVec[0]=="totbunchnumber")
        {
          ringFillPatt->totBunchNumber = stoi(strVec[1]);
        }        
        if(strVec[0]=="trainnumber")
        {
          ringFillPatt->trainNumber = stoi(strVec[1]);
          ringFillPatt->bunchNumberPerTrain.resize(ringFillPatt->trainNumber);
          ringFillPatt->trainGaps          .resize(ringFillPatt->trainNumber);
          ringFillPatt->bunchGaps          .resize(ringFillPatt->trainNumber);          
        }
        
        if(strVec[0]=="bunchnumberpertrain")
        {
          for(int i=1;i<=ringFillPatt->trainNumber;i++)  
          {
            ringFillPatt->bunchNumberPerTrain[i-1] = stoi(strVec[i]);         
          }  
        }               
        
        if(strVec[0]=="traingaps")
        {
          for(int i=1;i<=ringFillPatt->trainNumber;i++)
          { 
            ringFillPatt->trainGaps[i-1] = stoi(strVec[i]);
          }
        }        
        if(strVec[0]=="bunchgaps")
        {
          for(int i=1;i<=ringFillPatt->trainNumber;i++)
          { 
            ringFillPatt->bunchGaps[i-1] = stoi(strVec[i]);
          }
        }
        

        //4) initial bunch parameters
        if(strVec[0]=="current")
        {
          ringBunchPara->current = stod(strVec[1]);
        }
        if(strVec[0]=="emittancex")
        {
          ringBunchPara->emittanceX = stod(strVec[1]);
        }
        if(strVec[0]=="kappa")
        {
          ringBunchPara->kappa = stod(strVec[1]);
        }
        
        if(strVec[0]=="rmsbunchlength")
        {
          ringBunchPara->rmsBunchLength = stod(strVec[1]);
        }
        
        if(strVec[0]=="rmsenergyspread")
        {
         ringBunchPara-> rmsEnergySpread = stod(strVec[1]);
        }
        
        if(strVec[0]=="distributiontype")
        {
          ringBunchPara->distributionType = stod(strVec[1]);
        }
        if(strVec[0]=="macroelenumperbunch")
        {
          ringBunchPara->macroEleNumPerBunch = stoi(strVec[1]);
        } 
        
             
        if(strVec[0]=="initialstaticoffset")
        {
          for(int i=0;i<=5;i++)  
          {
            ringBunchPara->initialStaticOffSet[i] = stod(strVec[i+1]);                 
          }            
        }      
        
        if(strVec[0]=="initialdynamicoffset")
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
                        
        //5 initial ion efect para        
        if(strVec[0]=="gasspec")
        {         
          ringIonEffPara->gasSpec = stod(strVec[1]);
          ringIonEffPara->ionMassNumber.resize(ringIonEffPara->gasSpec);
          ringIonEffPara->corssSectionEI.resize(ringIonEffPara->gasSpec);
          ringIonEffPara->gasPercent.resize(ringIonEffPara->gasSpec);          
        }
        if(strVec[0]=="ionmassnumber")
        {                   
          for(int i=0;i<ringIonEffPara->gasSpec;i++)
          {
            ringIonEffPara->ionMassNumber[i] = stod(strVec[i+1]);
          }
        }            
        if(strVec[0]=="corsssectionei")
        {
          for(int i=0;i<ringIonEffPara->gasSpec;i++)
          {          
            ringIonEffPara->corssSectionEI[i] = stod(strVec[i+1]) * 1.e-22;
          }  
        }
        if(strVec[0]=="gaspercent")
        {
          for(int i=0;i<ringIonEffPara->gasSpec;i++)
          {          
            ringIonEffPara->gasPercent[i] = stod(strVec[i+1]);
          }  
        }


        
        if(strVec[0]=="ionmaxnumber")
        {
          ringIonEffPara->ionMaxNumber = stod(strVec[1]);
        }   
        if(strVec[0]=="ionlossboundary")
        {
          ringIonEffPara->ionLossBoundary = stod(strVec[1]);
        }    
        
        if(strVec[0]=="numberofionbeaminterpoint")
        {
          ringIonEffPara->numberofIonBeamInterPoint = stoi(strVec[1]);
        } 
        
        if(strVec[0]=="macroionnumbergeneratedperip")
        {   
           ringIonEffPara->macroIonNumberGeneratedPerIP = stod(strVec[1]);
        }
        if(strVec[0]=="ioninfoprintinterval")
        {
          ringIonEffPara->ionInfoPrintInterval = stod(strVec[1]);
        }   
        if(strVec[0]=="iondiswriteto")
        {
          ringIonEffPara->ionDisWriteTo = strVec[1];
        }   
        if(strVec[0]=="twissinput")
        {
            ringIonEffPara->twissInput = strVec[1];
        } 
       
            
        //6 initial FB efect para 

        if(strVec[0]=="kickstrengthkx")
        {
          ringBBFB->kickStrengthK[0] = stod(strVec[1]);
        }  
        if(strVec[0]=="kickstrengthky")
        {
          ringBBFB->kickStrengthK[1] = stod(strVec[1]);
        } 
        if(strVec[0]=="kickstrengthf")
        {
          ringBBFB->kickStrengthK[2] = stod(strVec[1]);
        }           
        if(strVec[0]=="gain")
        {
          ringBBFB->gain = stod(strVec[1]);
        }  
        if(strVec[0]=="delay")
        {
          ringBBFB->delay = stoi(strVec[1]);
        }  
        if(strVec[0]=="taps")
        {
          ringBBFB->taps = stoi(strVec[1]);
          ringBBFB->firCoeffx.resize(ringBBFB->taps + ringBBFB->delay);
          ringBBFB->firCoeffy.resize(ringBBFB->taps + ringBBFB->delay);
          ringBBFB->firCoeffz.resize(ringBBFB->taps + ringBBFB->delay);
          ringBBFB->firCoeffxy.resize(ringBBFB->taps + ringBBFB->delay);          
        }  
        if(strVec[0]=="kickerdispp")
        {
          ringBBFB->kickerDispP = stod(strVec[1]);
        }  
        if(strVec[0]=="kickerdisp")
        {
          ringBBFB->kickerDisp = stod(strVec[1]);
        }  
        if(strVec[0]=="firbunchbybunchfeedbackpowerlimit")
        {
          ringBBFB->fIRBunchByBunchFeedbackPowerLimit = stod(strVec[1]);
        }  
        if(strVec[0]=="firbunchbybunchfeedbackkickerimped")
        {
          ringBBFB->fIRBunchByBunchFeedbackKickerImped = stod(strVec[1]);
        }  
        
        if(strVec[0]=="fircoeffx")
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
        if(strVec[0]=="fircoeffy")
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
        if(strVec[0]=="fircoeffz")
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
        if(strVec[0]=="fircoeffxy")
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

        // 7)  initial wake 
        if(strVec[0]=="pipegeoinput")
        {
          ringWake->pipeGeoInput = strVec[1];          
        }
        if(strVec[0]=="bbrinput")
        {
          ringWake->bbrInput = strVec[1];          
        }
        if(strVec[0]=="rwflag")
        {
          ringWake->rwFlag = stoi(strVec[1]);          
        }
        if(strVec[0]=="bbrflag")
        {
          ringWake->bbrFlag = stoi(strVec[1]);          
        }
        if(strVec[0]=="nturnswaketrunction")
        {
          ringWake->nTurnswakeTrunction = stoi(strVec[1]);
        } 
                                          
        // 8) impedance 
        
        if(strVec[0]=="impedinput")
        {
          ringImpedance->impedInput = strVec[1];
        } 
        if(strVec[0]=="bunchbinnumberz")
        {
          ringImpedance->bunchBinNumberZ = stoi(strVec[1]);
        } 
        
 
        // 9) run       
        
        if(strVec[0]=="nturns")
        {
          ringRun->nTurns = stoi(strVec[1]);
        }
        
        if(strVec[0]=="calsetting")
        {
           ringRun->calSetting = stoi(strVec[1]);
        }
        if(strVec[0]=="synraddampingflag")
        {
           ringRun->synRadDampingFlag = stoi(strVec[1]);
        }

        if(strVec[0]=="firbunchbybunchfeedbackflag")
        {
           ringRun->fIRBunchByBunchFeedbackFlag = stoi(strVec[1]);
        }
              
        if(strVec[0]=="impedanceflag")
        {
           ringRun->impedanceFlag = stoi(strVec[1]);
        }          
        if(strVec[0]=="wakeflag")
        {
           ringRun->wakeFlag = stoi(strVec[1]);
        }
        if(strVec[0]=="tbtbunchdata")
        {
           ringRun->TBTBunchData = strVec[1];
        }
        if(strVec[0]=="beamionflag")
        {
           ringRun->beamIonFlag = stoi(strVec[1]);
        }
        if(strVec[0]=="bunchinfoprintinterval")
        {
          ringRun->bunchInfoPrintInterval = stoi(strVec[1]);
        }
                                                       
    }
    



    

    //1) data from ringPaBbasic
    double  circRing           = ringParBasic->circRing;
    double  sdelta0            = ringParBasic->sdelta0;
    double  alphac             = ringParBasic->alphac;
    double  electronBeamEnergy = ringParBasic->electronBeamEnergy;
    double  u0                 = ringParBasic->u0;
    
    //2) data from ringParRF
    double  rfBaseFreq = ringParRf->rfBaseFreq;
    double  vRF        = ringParRf->vRF;       
    double synPhase = PI - asin(u0/vRF);     
    ringParRf->synPhase = synPhase;

    double harmRF=ringParRf->harmRF; 
	if (harmRF!=0 && harmRF!=1)
	{
		double phaseTemp = pow(harmRF,2)/(pow(harmRF,2)-1) * u0 / vRF;			
		double synPhase = PI -  asin(phaseTemp);	
		
		
		phaseTemp = - harmRF * u0 / vRF / sqrt( pow(pow(harmRF,2)-1,2) - pow(harmRF,4) * pow(u0 / vRF,2)  );
		double phiHC = atan(phaseTemp);

		
		double kHC = sqrt( 1./pow(harmRF,2) - pow( u0 / vRF ,2) / ( pow(harmRF,2)-1)  ); 
		double vHC = kHC * vRF;	
        
        ringParRf->synPhase =  synPhase;
        ringParRf->synPhaseH = phiHC;
        ringParRf->vRFH      = vHC;                       
	}
        






    //4) data from ringBunchPara
    double kappa               = ringBunchPara->kappa;
    double emittanceX          = ringBunchPara->emittanceX;
    double emittanceZ          = ringBunchPara->rmsEnergySpread * ringBunchPara->rmsBunchLength;
    ringBunchPara-> emittanceZ = emittanceZ;
    
    double  rGamma = electronBeamEnergy/ElectronMassEV;
    double  rBeta  = sqrt(1.E0 -pow(rGamma,-2) );
    
    
    double  t0     = circRing /  rBeta /  CLight;
    double  f0     = 1/t0;
    double  omega0 = 2* PI * f0;
    int harmonics = floor(rfBaseFreq/f0);
 

 

    double eta    = alphac - 1./pow(rGamma,2);

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



