#include "Beam.h"
//#include "LongImpSingalBunch.h"
#include <fstream>
#include <stdlib.h> 
#include <stdio.h> 
#include <iostream>
#include "Faddeeva.h"
#include <time.h>
#include <string>
#include <cmath>
#include <gsl/gsl_fft_complex.h>
#include <complex>
#include<iomanip>
#include<cstring>
#include "FIRFeedBack.h"
#include "WakeFunction.h"




using namespace std;
using std::vector;
using std::complex;


Beam::Beam()
{

}

Beam::~Beam()
{
    delete weakStrongBunchInfo;
    delete strongStrongBunchInfo;
}

void Beam::Initial(Train &train, LatticeInterActionPoint &latticeInterActionPoint , ReadInputSettings &inputParameter)
{
    int totBunchNum = inputParameter.ringFillPatt->totBunchNumber;
    
    harmonics = inputParameter.ringParBasic->harmonics;
    ionLossBoundary = inputParameter.ringIonEffPara->ionLossBoundary;
    beamVec.resize(totBunchNum);



       
    for(int i=0;i<totBunchNum;i++)
    {
        beamVec[i].Initial(latticeInterActionPoint, inputParameter);
        beamVec[i].DistriGenerator(latticeInterActionPoint,inputParameter,i);
    }
   
  
        
    if (! inputParameter.ringBunchPara->bunchDisWriteTo.empty())
    {    
        PrintInitialBunchDistri(inputParameter);         
    }

         
           
    int counter=0;   
    for(int i=0;i<train.trainNumber;i++)
    {            
        for(int j=0;j<train.bunchNumberPerTrain[i];j++)
        {
			beamVec[counter].bunchHarmNum  = train.trainStart[i] + j * (train.bunchGaps[i] + 1);		
            counter++;             
        }
    }    
    
    for(int i=0;i<totBunchNum-1;i++)
    {
        beamVec[i].bunchGap = beamVec[i+1].bunchHarmNum - beamVec[i].bunchHarmNum;
    }
    beamVec[totBunchNum-1].bunchGap = harmonics - beamVec[totBunchNum-1].bunchHarmNum;
    	
	RMOutPutFiles(inputParameter);	
     
    /* 	   
    for(int i=0;i<totBunchNum;i++)
    {
        cout<<beamVec[i].bunchGap<<"    "<<i<<" "<<beamVec[i].bunchHarmNum<<endl;      
    }    
    getchar();
    */

} 


void Beam::RMOutPutFiles(ReadInputSettings &inputParameter)
{
    /*
    string filePrefix = inputParameter.ringIonEffPara->ionDisWriteTo;
    string fname = filePrefix + ".sdds";
    */
    string command ="rm -rf *.sdds";
    char *cstr = new char[command.length()+1];
    strcpy(cstr,command.c_str());    
    system(cstr);
    delete [] cstr;
}


  
void Beam::PrintInitialBunchDistri(ReadInputSettings &inputParameter)
{
    ofstream fout(inputParameter.ringBunchPara->bunchDisWriteTo+".sdds");

	fout<<"SDDS1"<<endl;
	fout<<"&parameter name=partNum, units=A, type=float,  &end"<<endl;
	fout<<"&column name=x,  units=m,   type=float,  &end"<<endl;
	fout<<"&column name=xp, units=rad, type=float,  &end"<<endl;
	fout<<"&column name=y,  units=m,   type=float,  &end"<<endl;
	fout<<"&column name=yp, units=rad, type=float,  &end"<<endl;
	fout<<"&column name=z,  units=m,   type=float,  &end"<<endl;	
	fout<<"&column name=dp,   units=rad,   type=float,  &end"<<endl;
	fout<<"&data mode=ascii, &end"<<endl;	

   
    for(int i=0;i<beamVec.size();i++)
    {
        fout<<"! page number "<<i+1<<endl;
        fout<<beamVec[i].electronNumPerBunch<<endl;
        fout<<beamVec[i].macroEleNumPerBunch<<endl;
        
        for(int j=0; j<beamVec[i].macroEleNumPerBunch; j++)
        {
            fout<<setw(15)<<left<<beamVec[i].ePositionX[j]
                <<setw(15)<<left<<beamVec[i].eMomentumX[j]
                <<setw(15)<<left<<beamVec[i].ePositionY[j]
                <<setw(15)<<left<<beamVec[i].eMomentumY[j]
                <<setw(15)<<left<<beamVec[i].ePositionZ[j]
                <<setw(15)<<left<<beamVec[i].eMomentumZ[j]
                <<endl;
        }        	
    }
    
    fout.close();   
}





void Beam::Run(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter)
{

//   calSettings =1 :  electron bunch is weak
//   calSettings =2 :  single_bunch_Longi_imped  --- signal-bunch
//   calSettings =3 :  single_bunch_Trans_imped  --- signal-bunch
//   calSettings =4 :  Long Range wakefield      --- multi-bunch-mulit-turn


 
    int nTurns                      = inputParameter.ringRun->nTurns;
    int calSettings                 = inputParameter.ringRun->calSetting;
    int synRadDampingFlag           = inputParameter.ringRun->synRadDampingFlag; 
    int fIRBunchByBunchFeedbackFlag = inputParameter.ringRun->fIRBunchByBunchFeedbackFlag;
    int impedanceFlag               = inputParameter.ringRun->impedanceFlag;
    int beamIonFlag                 = inputParameter.ringRun->beamIonFlag;
    int wakeFlag                    = inputParameter.ringRun->wakeFlag;             
    int totBunchNum                 = inputParameter.ringFillPatt->totBunchNumber;
    int ionInfoPrintInterval        = inputParameter.ringIonEffPara->ionInfoPrintInterval;
    int bunchInfoPrintInterval      = inputParameter.ringRun->bunchInfoPrintInterval;
 
 

 
                 
    // preapre the data for bunch-by-bunch system --------------- 
    FIRFeedBack firFeedBack;
    if(fIRBunchByBunchFeedbackFlag)
    {        
        firFeedBack.Initial(inputParameter);
    }
                   	
    // -----------------longRange wake function ---------------              
    WakeFunction wakeFunction;
    if(wakeFlag)
    { 
        wakeFunction.Initial(inputParameter);
    }



    ofstream fout ("result.sdds",ios::out);
    fout<<"SDDS1"<<endl;
        fout<<"&column name=Turns,                type=long,   &end"<<endl;	  
	    fout<<"&column name=IonCharge, units=e,   type=float,  &end"<<endl;
	    fout<<"&column name=MaxJx,     units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=MaxJy,     units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=MaxAverX,  units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=MaxAverY,  units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=rmsBunchX, units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=rmsBunchY, units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=rmsBunchZ, units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=rmsBunchZP, units=m,   type=float,  &end"<<endl;
	    fout<<"&data mode=ascii, &end"<<endl;	
	    fout<<"! page number "<<1<<endl;
	    fout<<nTurns<<endl;
    
    

        
    switch(calSettings)
    { 
        
        case 1 :    //   calSettings =1 :  electron bunch is weak--one macro-particle  
        {                         
            if(beamVec[0].macroEleNumPerBunch!=1)
            {
                cerr<<"reset the macroEleNumPerBunch=1 in Bunch namelist"<<endl;
                exit(0);                                     
            }
            // run loop starts, for nTrns and each trun for k interaction-points                
            for(int n=0;n<nTurns;n++)   
            { 					                                                           
                cout<<n<<"  turns"<<endl;
                
                if(beamIonFlag)
                {
                    for (int k=0;k<inputParameter.ringIonEffPara->numberofIonBeamInterPoint;k++)
                    {
                        WSBeamRMSCal(latticeInterActionPoint, k);
                        WSGetMaxBunchInfo();
                        WSBeamIonEffectOneInteractionPoint(inputParameter,latticeInterActionPoint, n, k); 
                        BeamTransferPerInteractionPointDueToLatticeT(latticeInterActionPoint,k);             //the longitudianl direciton is matrix notation                    
                        WSBeamRMSCal(latticeInterActionPoint, k);
                    }                                                                                       //--modify to data with rf setting Vrf * Sin[phases+deltaPhi] notation...          

                    BeamTransferPerTurnDueToLatticeL(inputParameter);
                    
                      
                    if(ionInfoPrintInterval && (n%ionInfoPrintInterval==0))
                    {
                        WSIonDataPrint(inputParameter,latticeInterActionPoint, n);
                    }                     				                         
                }
                
                                                 
                                                 
                if(wakeFlag)
                {                    
                    WSBeamRMSCal(latticeInterActionPoint, 0);
                    WSGetMaxBunchInfo();  
                    LRWakeBeamIntaction(inputParameter,wakeFunction);                
                    BeamTransferPerTurnDueWake();
                    BeamTransferPerTurnDueToLattice(latticeInterActionPoint,inputParameter); 
                }
 
                
                if(synRadDampingFlag)
                {
                     BeamSynRadDamping(inputParameter,latticeInterActionPoint);                    
                }
                if(fIRBunchByBunchFeedbackFlag)
                {
                    FIRBunchByBunchFeedback(firFeedBack,n);  
                }    
                
                if(bunchInfoPrintInterval && (n%bunchInfoPrintInterval==0))
                {   
                    WSBeamRMSCal(latticeInterActionPoint, 0);
                    WSGetMaxBunchInfo();
                    WSBeamDataPrintPerTurn(n,latticeInterActionPoint,inputParameter);
                }  
                                                                                                                                                   
                

                fout<<n<<"  "
                    <<setw(15)<<left<< latticeInterActionPoint.totIonCharge
                    <<setw(15)<<left<< log10(sqrt(weakStrongBunchInfo->actionJxMax))
                    <<setw(15)<<left<< log10(sqrt(weakStrongBunchInfo->actionJyMax))
                    <<setw(15)<<left<< weakStrongBunchInfo->bunchAverXMax
                    <<setw(15)<<left<< weakStrongBunchInfo->bunchAverYMax
                    <<setw(15)<<left<< weakStrongBunchInfo->bunchRmsSizeX
                    <<setw(15)<<left<< weakStrongBunchInfo->bunchRmsSizeY
                    <<setw(15)<<left<<weakStrongBunchInfo->bunchZ
                    <<setw(15)<<left<<weakStrongBunchInfo->bunchZP
                    <<endl; 
                                       
             }                                 
        }
                     
        default:
            cerr<<"Do Nothing...  "<<endl;
        
   
    }
      
    fout.close();
    
    
}



void Beam:: LRWakeBeamIntaction(ReadInputSettings &inputParameter, WakeFunction &wakefunction)
{
 

    int nTurnswakeTrunction     = inputParameter.ringWake->nTurnswakeTrunction;
    int harmonics               = inputParameter.ringParBasic->harmonics;
    double electronBeamEnergy   = inputParameter.ringParBasic->electronBeamEnergy;
   
    
    // prepare the previous turns data 
    vector<double> posxDataTemp(beamVec.size());
    vector<double> posyDataTemp(beamVec.size());
    vector<double> poszDataTemp(beamVec.size());

    for(int i=0;i<beamVec.size();i++)
    {   
        posxDataTemp[i] = beamVec[i].xAver;
        posyDataTemp[i] = beamVec[i].yAver;
        poszDataTemp[i] = beamVec[i].zAver; 
    }

    
    wakefunction.posxData.erase(wakefunction.posxData.begin());
    wakefunction.posyData.erase(wakefunction.posyData.begin());
    wakefunction.poszData.erase(wakefunction.poszData.begin());
    
    wakefunction.posxData.push_back(posxDataTemp);
    wakefunction.posyData.push_back(posyDataTemp);
    wakefunction.poszData.push_back(poszDataTemp);
    //--------------------


    vector<double> wakeForceTemp(3,0);

    double tauij=0.e0;
    double tRF  = inputParameter.ringParBasic->t0 / double(harmonics);
    
    
    
    for (int i=0;i<beamVec.size();i++)
    {    
	    beamVec[i].wakeForceAver[0] =0.e0;      // x
	    beamVec[i].wakeForceAver[1] =0.e0;      // y
	    beamVec[i].wakeForceAver[2] =0.e0;      // z
   
	    for(int j=0;j<beamVec.size();j++)
	    {
		    for(int n=0;n<nTurnswakeTrunction;n++)
		    {
			    if(j<i)
			    {
			        tauij = ( (n + 1) * harmonics + beamVec[j].bunchHarmNum - beamVec[i].bunchHarmNum) * tRF;
			    }
			    else
			    {
			        tauij = (  n      * harmonics + beamVec[j].bunchHarmNum - beamVec[i].bunchHarmNum) * tRF;
			    }
                
                // resistive wall wake calculation			    			    
			    if(inputParameter.ringWake->rwFlag)
			    {
			        if(tauij==0) 	// assume bunch does not affect itself
                    {
                        wakeForceTemp={0,0,0};
                    }
                    else
                    {
                        wakeForceTemp=wakefunction.GetRWWakeForce(tauij );
       		    	}
       		    			    
			        beamVec[i].wakeForceAver[0] += wakeForceTemp[0] * beamVec[j].electronNumPerBunch * wakefunction.posxData[nTurnswakeTrunction-1-n][j];    //[V/C/m] * [m] ->  [V/C]
			        beamVec[i].wakeForceAver[1] += wakeForceTemp[1] * beamVec[j].electronNumPerBunch * wakefunction.posyData[nTurnswakeTrunction-1-n][j];    //[V/C/m] * [m] ->  [V/C]
			        beamVec[i].wakeForceAver[2] += wakeForceTemp[2] * beamVec[j].electronNumPerBunch;                                                        //[V/C]         ->  [V/C]		    

		        }   
		        //---------------------------------------------------------------------------------
		        		        
		        // BBR long range wake function calculation
		        if(inputParameter.ringWake->bbrFlag)
			    {		        		            
                    wakeForceTemp=wakefunction.GetBBRWakeForce(tauij);
			        beamVec[i].wakeForceAver[0] += wakeForceTemp[0] * beamVec[j].electronNumPerBunch * wakefunction.posxData[nTurnswakeTrunction-1-n][j];    //[V/C/m] * [m] ->  [V/C]
			        beamVec[i].wakeForceAver[1] += wakeForceTemp[1] * beamVec[j].electronNumPerBunch * wakefunction.posyData[nTurnswakeTrunction-1-n][j];    //[V/C/m] * [m] ->  [V/C]
			        beamVec[i].wakeForceAver[2] += wakeForceTemp[2] * beamVec[j].electronNumPerBunch;                                                        //[V/C]         ->  [V/C]
		        }
		        //---------------------------------------------------------------------------------      		        		    
		    }		          
	    }
        
    

        beamVec[i].wakeForceAver[0] *= ElectronCharge / electronBeamEnergy ;  // [V/C] * [C] * [1e] / [eV] ->rad   	    
        beamVec[i].wakeForceAver[1] *= ElectronCharge / electronBeamEnergy ;  // [V/C] * [C] * [1e] / [eV] ->rad  
        beamVec[i].wakeForceAver[2] *= ElectronCharge / electronBeamEnergy ;  // [V/C] * [C] * [1e] / [eV] ->rad
        

    }
        
    // Reference to "simulation of transverse multi-bunch instabilities of proton beam in LHC, PHD thesis, Alexander Koshik P. 32, Eq. (3.22)"
} 

void Beam::BeamTransferPerTurnDueWake()
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToWake();
    }
}


void Beam::BeamTransferPerInteractionPointDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToLatticeT(latticeInterActionPoint,k);
    }
}


void Beam::BeamTransferPerTurnDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint)
{
    for(int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)
    {
        BeamTransferPerInteractionPointDueToLatticeT(latticeInterActionPoint,k); 
    }
    
}

void Beam::BeamTransferPerTurnDueToLatticeL(ReadInputSettings &inputParameter) 
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToLatticeL(inputParameter); 
    }
}

void Beam::BeamTransferPerTurnDueToLattice(LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter)
{
    BeamTransferPerTurnDueToLatticeT(latticeInterActionPoint);
    BeamTransferPerTurnDueToLatticeL(inputParameter);
}




void Beam::BeamSynRadDamping(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint)
{
    for(int j=0;j<beamVec.size();j++)
    {    
        beamVec[j].BunchSynRadDamping(inputParameter,latticeInterActionPoint);        
    }
}

void Beam::WSBeamRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].WSRMSCal(latticeInterActionPoint,k);    
    }
}





void Beam::WSBeamIonEffectOneInteractionPoint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k)
{
    int totBunchNum = beamVec.size();    
    
    
           
    for(int j=0;j<totBunchNum;j++)
    {	
        latticeInterActionPoint.GetIonNumberPerInterAction(beamVec[j].electronNumPerBunch, k);
        latticeInterActionPoint.IonGenerator(beamVec[j].rmsRx,beamVec[j].rmsRy,beamVec[j].xAver,beamVec[j].yAver,k);
        latticeInterActionPoint.WSIonsUpdate(k);
        //latticeInterActionPoint.IonRMSCal(k);  
        beamVec[j].WSIonBunchInteraction(latticeInterActionPoint,k);
        beamVec[j].BunchTransferDueToIon(latticeInterActionPoint,k); 
        latticeInterActionPoint.IonTransferDueToBunch(beamVec[j].bunchGap,k,weakStrongBunchInfo->bunchEffectiveSizeXMax,weakStrongBunchInfo->bunchEffectiveSizeYMax);     
    }
        
      /*
        ofstream fout ("test.sdds",ios_base::app);   
        for(int i=0;i<latticeInterActionPoint.ionAccumuPositionX[k].size();i++)
        {
            fout<<setw(15)<<left<<latticeInterActionPoint.ionAccumuPositionX[k][i]
                <<setw(15)<<left<<latticeInterActionPoint.ionAccumuPositionY[k][i]
                <<setw(15)<<left<<latticeInterActionPoint.ionAccumuFx[k][i]
                <<setw(15)<<left<<latticeInterActionPoint.ionAccumuFy[k][i]
                <<endl; 
        }
        fout.close(); 
        */       
               
}




//// //void Beam::SSBeamRMSCal()
//// //{
//// //   for(int j=0;j<beamVec.size();j++)
//// //    {
//// //        beamVec[j].SSRMSCal(latticeInterActionPoint,k);
//// //    }
//// //}
//// 
//// 
//// 
//// 
// void Beam::SSBeamIonEffectOneTurn( LatticeInterActionPoint &latticeInterActionPoint, ReadInputSettings &inputParameter, int nTurns, int intevalofTurnsIonDataPrint, int printInterval)
//// {
//// 
////     double corssSectionEI = inputParameter.corssSectionEI;
////     int totBunchNum = inputParameter.totBunchNumber;
//// 
////     int counter=nTurns*latticeInterActionPoint.numberOfInteraction*totBunchNum;
////     
////     for(int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)  
////     {
////         for(int j=0;j<totBunchNum;j++)
////         {
//// 
////             latticeInterActionPoint.ionLineDensity[k] = corssSectionEI * latticeInterActionPoint.vacuumPressure[k]
////                                     /latticeInterActionPoint.temperature[k]/Boltzmann * beamVec[j].electronNumPerBunch;
////             //(1) get the ion beam line density at the kth interaction points due to the jth beam bunch
//// 
////             latticeInterActionPoint.ionNumber[k] = latticeInterActionPoint.ionLineDensity[k] * 
////                                                    latticeInterActionPoint.interactionLength[k];
////             //(2) get the ion beam number to be generated at the kth interaction point due to the jth beam bunch
//// 
////             beamVec[j].RMSCal(latticeInterActionPoint,k);
//// 
////             //(3) get the rms size of jth beam at kth interaction point, rmsRx rmsRy
////             latticeInterActionPoint.IonGenerator(beamVec[j].rmsRx,beamVec[j].rmsRy,beamVec[j].xAver,beamVec[j].yAver,k);
////             //(4) generate the ions randomly (Gausion distribution truncated n*rms) at the kth interaction point 
//// 
//// 
////             latticeInterActionPoint.SSIonsUpdate(bunchEffectiveSizeXMax, bunchEffectiveSizeYMax,beamVec[j].xAver,beamVec[j].yAver,k);
////             
////             latticeInterActionPoint.IonRMSCal(k);
////             //(5) get the accumulated ions information at the kth point (new ions plus old ions)
////             //    get the rms vaues of the accumulated ion distribution.
//// 
////             beamVec[j].SSIonBunchInteraction(latticeInterActionPoint,k);
////             //(6) get the interactionforces between accumulated ions and the jth beam bunch at the kth interaction point 
//// 
////             //if(intevalofTurnsIonDataPrint && (nTurns%intevalofTurnsIonDataPrint==0) && j==beamVec.size()-1)
////             //{
////             //    SSBunchDataPrint(beamVec[j],nTurns);
////             //}
//// 
////             //if(printInterval && (nTurns%printInterval==0) && j==beamVec.size()-1)
////             //{
////             //   SSIonDataPrint1(latticeInterActionPoint, nTurns, k);
////             //}
//// 
//// 
////             beamVec[j].BunchTransferDueToIon(latticeInterActionPoint,k); 
////             //(7) kick the jth electron beam bunch due to accumulated ions at the kth interactoin point 
//// 
////             latticeInterActionPoint.IonTransferDueToBunch(beamVec[j].bunchGap,k,bunchEffectiveSizeXMax, bunchEffectiveSizeYMax, beamVec[j].macroEleNumPerBunch);         
////             //(8) kick the accumulated ions due to jth electron bunch 
////             
////             
////             beamVec[j].BunchTransferDueToLattice(latticeInterActionPoint,k); 
////             //(9) electron beam transfer to next interaction point 
//// 
//// 
////             if(nTurns%printInterval==0 && k==0)
////             {
////                 bunchInfoOneTurn[j][0] = nTurns;                                           //1
////                 bunchInfoOneTurn[j][1] = j;                                                //2
////                 bunchInfoOneTurn[j][2] = latticeInterActionPoint.ionAccumuNumber[k];       //3
////                 bunchInfoOneTurn[j][3] = beamVec[j].rmsEffectiveRx;                        //4
////                 bunchInfoOneTurn[j][4] = beamVec[j].rmsEffectiveRy;                        //5
////                 bunchInfoOneTurn[j][5] = beamVec[j].rmsEffectiveRingEmitX;                 //6
////                 bunchInfoOneTurn[j][6] = beamVec[j].rmsEffectiveRingEmitX;                 //7
////                 bunchInfoOneTurn[j][7] = beamVec[j].xAver;                                 //8
////                 bunchInfoOneTurn[j][8] = beamVec[j].yAver;                                 //9    
////                 bunchInfoOneTurn[j][9] = beamVec[j].pxAver;                                //10
////                 bunchInfoOneTurn[j][10] = beamVec[j].pyAver;                               //11
////                 bunchInfoOneTurn[j][11] = latticeInterActionPoint.ionAccumuAverX[k];       //12
////                 bunchInfoOneTurn[j][12] = latticeInterActionPoint.ionAccumuAverY[k];       //13
////                 bunchInfoOneTurn[j][13] = latticeInterActionPoint.ionAccumuRMSX[k];        //14
////                 bunchInfoOneTurn[j][14] = latticeInterActionPoint.ionAccumuRMSY[k];        //15
////             }
////             counter++;
////         }
////     }
//// }
//// 
//// 
//// 
//// 
//// 
//// 
//// 
//// void Beam::SingleBunchLongiImpedInterAction(LongImpSingalBunch &longImpSingalBunch, int nTurns, ofstream &fout)
//// {
//// 
////     cout<<"do the beam_impedance in fre domain"<<endl;
////     cout<<"only the beamVec[0] is used in calculation" <<endl;
//// 
////     beamVec[0].InitialBeamCurDenZProf(); 
//// 
////     for(int i=0;i<longImpSingalBunch.longImpedR.size();i++)
////     {
//// //        beamVec[0];
////     }
//// 
//// 
//// }
//// 
//// void Beam::SSIonDataPrint(LatticeInterActionPoint &latticeInterActionPoint, int nTurns)
//// {
////     for(int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)
////     {
//// //        string fname;
//// //        fname = "SSionAccumudis_"+std::to_string(nTurns)+"_"+to_string(k)+".dat";
////         
////         char fname[256];
////         sprintf(fname, "SSionAccumudis_%d_%d.dat", nTurns, k);
////         
////         ofstream fout(fname,ios::out);
////         for(int i=0;i<latticeInterActionPoint.ionAccumuNumber[k];i++)
////         {
////             fout<<i                     <<"      "
////                 <<latticeInterActionPoint.ionAccumuPositionX[k][i]   <<"      "       
////                 <<latticeInterActionPoint.ionAccumuPositionY[k][i]   <<"      "       
////                 <<latticeInterActionPoint.ionAccumuVelocityX[k][i]   <<"      "       
////                 <<latticeInterActionPoint.ionAccumuVelocityY[k][i]   <<"      "       
////                 <<latticeInterActionPoint.ionAccumuFx[k][i]   <<"      "       
////                 <<latticeInterActionPoint.ionAccumuFy[k][i]   <<"      "       
////                 <<endl;
////         }
////         
////         fout.close();
////     }
//// 
//// 
////     vector<double> histogrameX;
////     vector<double> histogrameY;
////     
////     double x_max = 0;
////     double x_min = 0;
////     
////     double y_max = 0;
////     double y_min = 0;
////     int bins = 100;
////     histogrameX.resize(bins);
////     histogrameY.resize(bins);
////     
////     
////     double binSizeX =0;
////     double binSizeY =0;
////     
//// 
////     double temp;
////     int indexX;
////     int indexY;
////     
////     for(int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)
////     {
//// 
////         x_max =  ionLossBoundary*bunchEffectiveSizeXMax;
////         x_min = -ionLossBoundary*bunchEffectiveSizeXMax;
////         y_max =  ionLossBoundary*bunchEffectiveSizeYMax;
////         y_min = -ionLossBoundary*bunchEffectiveSizeYMax;
//// 
//// //        x_max = latticeInterActionPoint.ionAccumuAverX[k]  + latticeInterActionPoint.pipeAperatureX;
//// //        x_min = latticeInterActionPoint.ionAccumuAverX[k]  - latticeInterActionPoint.pipeAperatureX;
//// //        y_max = latticeInterActionPoint.ionAccumuAverY[k]  + latticeInterActionPoint.pipeAperatureY;
//// //        y_min = latticeInterActionPoint.ionAccumuAverY[k]  - latticeInterActionPoint.pipeAperatureY;
//// 
//// 
//// 
////         binSizeX  = (x_max - x_min)/bins;
////         binSizeY  = (y_max - y_min)/bins;
//// 
////         for(int i=0;i<latticeInterActionPoint.ionAccumuNumber[k];i++)
////         {
////             temp   = (latticeInterActionPoint.ionAccumuPositionX[k][i] - x_min)/binSizeX;
////             indexX = floor(temp);
////             
////             temp   = (latticeInterActionPoint.ionAccumuPositionY[k][i] - y_min)/binSizeY;
////             indexY = floor(temp);
////             if(indexX<0 || indexY<0 || indexX>bins || indexY>bins) continue;
////             
////             histogrameX[indexX] = histogrameX[indexX] +1;
////             histogrameY[indexY] = histogrameY[indexY] +1;
////         }
////         
////         
////         
//// //        string fname;
//// //        fname = "SSionAccumuHis_"+to_string(nTurns)+"_"+to_string(k)+".dat";
////         
////         char fname[256];
////         sprintf(fname, "SSionAccumuHis_%d_%d.dat", nTurns, k);
////         
////         
////         ofstream foutHis(fname,ios::out);
//// 
////         for(int i=0;i<bins;i++)
////         {
////             foutHis<<i<<"    "
////                    <<(x_min+i*binSizeX)/bunchEffectiveSizeXMax<<"  "
////                    <<(y_min+i*binSizeY)/bunchEffectiveSizeYMax<<"  "
////                    <<(x_min+i*binSizeX)<<"  "
////                    <<(y_min+i*binSizeY)<<"  "
////                    <<histogrameX[i]<<"  "
////                    <<histogrameY[i]<<"  "
////                    <<endl;
////         }
////      
////         foutHis.close();
////     }
//// 
//// }
//// 
//// 
//// 
//// 
//// void Beam::SSIonDataPrint1(LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k)
//// {
//// 
//// 
//// 
//// //    string fname;
//// //    fname = "SSionAccumudis_"+to_string(nTurns)+"_"+to_string(k)+".dat";
////     
////     char fname[256];
////     sprintf(fname, "SSionAccumudis_%d_%d.dat", nTurns, k);
////     
////     ofstream fout(fname,ios::out);
////     for(int i=0;i<latticeInterActionPoint.ionAccumuNumber[k];i++)
////     {
////         fout<<i                     <<"      "
////             <<latticeInterActionPoint.ionAccumuPositionX[k][i]   <<"      "       
////             <<latticeInterActionPoint.ionAccumuPositionY[k][i]   <<"      "       
////             <<latticeInterActionPoint.ionAccumuVelocityX[k][i]   <<"      "       
////             <<latticeInterActionPoint.ionAccumuVelocityY[k][i]   <<"      "       
////             <<latticeInterActionPoint.ionAccumuFx[k][i]   <<"      "       
////             <<latticeInterActionPoint.ionAccumuFy[k][i]   <<"      "       
////             <<endl;
////     }
////     
////     fout.close();
//// 
//// 
//// 
////     vector<double> histogrameX;
////     vector<double> histogrameY;
////     
////     double x_max = 0;
////     double x_min = 0;
////     
////     double y_max = 0;
////     double y_min = 0;
////     int bins = 100;
////     histogrameX.resize(bins);
////     histogrameY.resize(bins);
////     
////     
////     double binSizeX =0;
////     double binSizeY =0;
////     
//// 
////     double temp;
////     int indexX;
////     int indexY;
////     
//// 
//// 
////     x_max =  ionLossBoundary*bunchEffectiveSizeXMax;
////     x_min = -ionLossBoundary*bunchEffectiveSizeXMax;
////     y_max =  ionLossBoundary*bunchEffectiveSizeYMax;
////     y_min = -ionLossBoundary*bunchEffectiveSizeYMax;
//// 
//// //  x_max = latticeInterActionPoint.ionAccumuAverX[k]  + latticeInterActionPoint.pipeAperatureX;
//// //  x_min = latticeInterActionPoint.ionAccumuAverX[k]  - latticeInterActionPoint.pipeAperatureX;
//// //  y_max = latticeInterActionPoint.ionAccumuAverY[k]  + latticeInterActionPoint.pipeAperatureY;
//// //  y_min = latticeInterActionPoint.ionAccumuAverY[k]  - latticeInterActionPoint.pipeAperatureY;
//// 
//// 
//// 
////     binSizeX  = (x_max - x_min)/bins;
////     binSizeY  = (y_max - y_min)/bins;
//// 
////     for(int i=0;i<latticeInterActionPoint.ionAccumuNumber[k];i++)
////     {
////         temp   = (latticeInterActionPoint.ionAccumuPositionX[k][i] - x_min)/binSizeX;
////         indexX = floor(temp);
////         
////         temp   = (latticeInterActionPoint.ionAccumuPositionY[k][i] - y_min)/binSizeY;
////         indexY = floor(temp);
////         if(indexX<0 || indexY<0 || indexX>bins || indexY>bins) continue;
////         
////         histogrameX[indexX] = histogrameX[indexX] +1;
////         histogrameY[indexY] = histogrameY[indexY] +1;
////     }
////     
////     
////     
//// 
//// 
//// 
//// //    fname = "SSionAccumuHis_"+to_string(nTurns)+"_"+to_string(k)+".dat";
////     sprintf(fname, "SSionAccumuHis_%d_%d.dat", nTurns, k);
////     ofstream foutHis(fname,ios::out);
//// 
////     for(int i=0;i<bins;i++)
////     {
////         foutHis<<i<<"    "
////                <<(x_min+i*binSizeX)/bunchEffectiveSizeXMax<<"  "
////                <<(y_min+i*binSizeY)/bunchEffectiveSizeYMax<<"  "
////                <<(x_min+i*binSizeX)<<"  "
////                <<(y_min+i*binSizeY)<<"  "
////                <<histogrameX[i]<<"  "
////                <<histogrameY[i]<<"  "
////                <<endl;
////     }
////  
////     foutHis.close();
//// 
//// }
//// 
//// 
//// 
//// 
//// void Beam::SSBunchDataPrint( Bunch &bunch, int nTurn) 
//// {
//// 
//// 
//// //    string fname;
//// //    fname  = "SSbeamdis_"+to_string(nTurn)+".dat";
////     
////     char fname[256];
////     sprintf(fname, "SSbeamdis_%d.dat", nTurn);
////     
////     ofstream fout(fname);
//// 
////     for(int i=0;i<bunch.macroEleNumPerBunch;i++)
////     {
//// 
////         fout<<i                     <<"      "
////             <<bunch.ePositionX[i]   <<"      "       
////             <<bunch.ePositionY[i]   <<"      "       
////             <<bunch.eMomentumX[i]   <<"      "       
////             <<bunch.eMomentumY[i]   <<"      "       
////             <<bunch.eFx[i]          <<"      "       
////             <<bunch.eFy[i]          <<"      "       
////             <<endl;
////     }
////     fout.close();
//// }
//// 
//// 
//// 
void Beam::WSBeamDataPrintPerTurn(int nTurns, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter)
{   
     // print out the ion distribution at different interaction points
    int numberOfInteraction =  latticeInterActionPoint.numberOfInteraction;
       
    string filePrefix = inputParameter.ringRun->TBTBunchData;
    string fname = filePrefix + ".sdds";    
    ofstream fout(fname,ios_base::app);
    if(nTurns==0)
	{
	    fout<<"SDDS1"<<endl;
	    fout<<"&column name=Turns,             type=long,  &end"<<endl;	  
	    fout<<"&column name=HarmIndex,         type=long,  &end"<<endl;
	    fout<<"&column name=AverX,  units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=AverY,  units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=AverZ,  units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=AverXP, units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=AverYP, units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=AverZP, units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=Jx,     units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=Jy,     units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=dpx,    units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=dpy,    units=rad, type=float,  &end"<<endl;	    	    
	    fout<<"&data mode=ascii, &end"<<endl;	
	}
	        
    fout<<"! page number "<<nTurns + 1<<endl;
    fout<<beamVec.size()<<endl;
    //fout<<1<<endl;
    for(int i=0;i<beamVec.size();i++)
    {
        fout<<setw(15)<<left<<nTurns
            <<setw(15)<<left<<beamVec[i].bunchHarmNum
            <<setw(15)<<left<<beamVec[i].xAver       
            <<setw(15)<<left<<beamVec[i].yAver
            <<setw(15)<<left<<beamVec[i].zAver
            <<setw(15)<<left<<beamVec[i].pxAver
            <<setw(15)<<left<<beamVec[i].pyAver
            <<setw(15)<<left<<beamVec[i].pzAver
            <<setw(15)<<left<<beamVec[i].actionJx
            <<setw(15)<<left<<beamVec[i].actionJy
            <<setw(15)<<left<<beamVec[i].eFx[0]
            <<setw(15)<<left<<beamVec[i].eFy[0]                    
            <<endl;
    }        
    fout.close();     
}


void Beam::WSIonDataPrint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int nTurns)
{
     
    // print out the ion distribution at different interaction points
    int numberOfInteraction =  latticeInterActionPoint.numberOfInteraction;
       
    string filePrefix = inputParameter.ringIonEffPara->ionDisWriteTo;
    string fname = filePrefix + ".sdds";    
    ofstream fout(fname,ios_base::app);
    if(nTurns==0)
	{
	    fout<<"SDDS1"<<endl;
	    fout<<"&parameter name=totIonCharge, units=e, type=float,  &end"<<endl;
	    fout<<"&column name=ionMass,                   type=float,  &end"<<endl;
	    fout<<"&column name=macroIonCharge,  units=e,   type=float,  &end"<<endl;
	    fout<<"&column name=x,  units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=y,  units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=xp, units=m/s, type=float,  &end"<<endl;
	    fout<<"&column name=yp, units=m/s, type=float,  &end"<<endl;
	    fout<<"&column name=Fx, units=m/s, type=float,  &end"<<endl;
	    fout<<"&column name=Fy, units=m/s, type=float,  &end"<<endl;
	    fout<<"&data mode=ascii, &end"<<endl;	
	}
	        
    for(int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)
    {
        fout<<"! page number "<<nTurns * numberOfInteraction + k + 1 <<endl;
        fout<<latticeInterActionPoint.totIonCharge<<endl;
        fout<<latticeInterActionPoint.totMacroIonsAtInterPoint[k]<<endl;               
        
        for(int p=0;p<latticeInterActionPoint.gasSpec;p++)
        {
            for(int i=0;i<latticeInterActionPoint.ionAccumuNumber[k][p];i++)
            {
                fout<<setw(15)<<left<<latticeInterActionPoint.ionMassNumber[p]
                    <<setw(15)<<left<<latticeInterActionPoint.macroIonCharge[k][p]
                    <<setw(15)<<left<<latticeInterActionPoint.ionAccumuPositionX[k][p][i]       
                    <<setw(15)<<left<<latticeInterActionPoint.ionAccumuPositionY[k][p][i]   
                    <<setw(15)<<left<<latticeInterActionPoint.ionAccumuVelocityX[k][p][i]       
                    <<setw(15)<<left<<latticeInterActionPoint.ionAccumuVelocityY[k][p][i]       
                    <<setw(15)<<left<<latticeInterActionPoint.ionAccumuFx[k][p][i]      
                    <<setw(15)<<left<<latticeInterActionPoint.ionAccumuFy[k][p][i]       
                    <<endl;
            }
        }
    }        
    fout.close();


/*    
 
//    fname = "WSionCenter_"+to_string(nTurns)+"_"+to_string(k)+".dat";
    sprintf(fname, "WSionCenter_%d_%d.dat", nTurns, k);
    ofstream fionCenter(fname);


    fionCenter<<latticeInterActionPoint.ionAccumuAverX[k]<<"      "       
              <<latticeInterActionPoint.ionAccumuAverY[k]<<"      "
              <<latticeInterActionPoint.ionAccumuRMSX[k]<<"      "
              <<latticeInterActionPoint.ionAccumuRMSY[k]<<"      "
              <<latticeInterActionPoint.ionAccumuAverVelX[k]<<"      "
              <<latticeInterActionPoint.ionAccumuAverVelY[k]<<"      "
              <<endl;

    fout.close();
    
    


    vector<double> histogrameX;
    vector<double> histogrameY;
    
    double x_max = 0;
    double x_min = 0;
    
    double y_max = 0;
    double y_min = 0;
    int bins = 100;
    histogrameX.resize(bins);
    histogrameY.resize(bins);
    
    
    double binSizeX =0;
    double binSizeY =0;
    

    double temp;
    int indexX;
    int indexY;
    
    for(int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)
    {

        if(ionLossBoundary*bunchEffectiveSizeXMax>latticeInterActionPoint.pipeAperatureX)
        {
            x_max =  latticeInterActionPoint.pipeAperatureX;
            x_min = -latticeInterActionPoint.pipeAperatureX;
        }
        else
        {
            x_max =  ionLossBoundary*bunchEffectiveSizeXMax;
            x_min = -ionLossBoundary*bunchEffectiveSizeXMax;
        }
        
        if(ionLossBoundary*bunchEffectiveSizeYMax>latticeInterActionPoint.pipeAperatureY)
        {
            y_max =  latticeInterActionPoint.pipeAperatureY;
            y_min = -latticeInterActionPoint.pipeAperatureY;
        }
        else
        {
            y_max =  ionLossBoundary*bunchEffectiveSizeYMax;
            y_min = -ionLossBoundary*bunchEffectiveSizeYMax;
        }




//        x_max = latticeInterActionPoint.ionAccumuAverX[k]  + latticeInterActionPoint.pipeAperatureX;
//        x_min = latticeInterActionPoint.ionAccumuAverX[k]  - latticeInterActionPoint.pipeAperatureX;
//        y_max = latticeInterActionPoint.ionAccumuAverY[k]  + latticeInterActionPoint.pipeAperatureY;
//        y_min = latticeInterActionPoint.ionAccumuAverY[k]  - latticeInterActionPoint.pipeAperatureY;



        binSizeX  = (x_max - x_min)/bins;
        binSizeY  = (y_max - y_min)/bins;

        for(int i=0;i<latticeInterActionPoint.ionAccumuNumber[k];i++)
        {
            temp   = (latticeInterActionPoint.ionAccumuPositionX[k][i] - x_min)/binSizeX;
            indexX = floor(temp);
            
            temp   = (latticeInterActionPoint.ionAccumuPositionY[k][i] - y_min)/binSizeY;
            indexY = floor(temp);
            if(indexX<0 || indexY<0 || indexX>bins || indexY>bins) continue;
            
            histogrameX[indexX] = histogrameX[indexX] +1;
            histogrameY[indexY] = histogrameY[indexY] +1;
        }
        

        sprintf(fname, "WSionAccumuHis_%d_%d.dat", nTurns, k);
        ofstream foutHis(fname,ios::out);

        for(int i=0;i<bins;i++)
        {
            foutHis<<i<<"    "
                   <<(x_min+i*binSizeX)/bunchEffectiveSizeXMax<<"  "
                   <<(y_min+i*binSizeY)/bunchEffectiveSizeYMax<<"  "
                   <<(x_min+i*binSizeX)<<"  "
                   <<(y_min+i*binSizeY)<<"  "
                   <<histogrameX[i]<<"  "
                   <<histogrameY[i]<<"  "
                   <<endl;
        }

        foutHis.close();
    }
    
    
  */  
    
}



// void Beam::SSGetMaxBunchInfo()
// {
//     double totBunchNum = beamVec.size();
//     double tempEmitx;
//     double tempEmity;
//     double tempEffectiveEmitx;
//     double tempEffectiveEmity;
//     double tempRMSSizeX;
//     double tempRMSSizeY;
//     double tempEffectiveRMSSizeX;
//     double tempEffectiveRMSSizeY;
//     double tempSizeX;
//     double tempSizeY;
//     
//     double tempAverX;
//     double tempAverY;
//     
//     tempEmitx   = beamVec[0].rmsEmitX;
//     tempEmity   = beamVec[0].rmsEmitY;
// 
//     tempEffectiveEmitx   = beamVec[0].rmsEffectiveRingEmitX;
//     tempEffectiveEmity   = beamVec[0].rmsEffectiveRingEmitY;
// 
//     tempRMSSizeX   = beamVec[0].rmsRx;
//     tempRMSSizeY   = beamVec[0].rmsRy;
// 
//     tempEffectiveRMSSizeX = beamVec[0].rmsEffectiveRx;
//     tempEffectiveRMSSizeY = beamVec[0].rmsEffectiveRy;
// 
// 
//     tempAverX   = beamVec[0].xAver;
//     tempAverY   = beamVec[0].yAver;
//     
//     
//     for(int i=1;i<totBunchNum;i++)
//     {
//         if(beamVec[i].rmsEmitX>tempEmitx)
//         {
//             tempEmitx = beamVec[i].rmsEmitX;
//         }
//         
//         if(beamVec[i].rmsEmitY>tempEmity)
//         {
//             tempEmity = beamVec[i].rmsEmitY;
//         }
// 
// 
//         if(beamVec[i].rmsEffectiveRingEmitX>tempEffectiveEmitx)
//         {
//             tempEffectiveEmitx = beamVec[i].rmsEffectiveRingEmitX;
//         }
//         
//         if(beamVec[i].rmsEffectiveRingEmitY>tempEffectiveEmity)
//         {
//             tempEffectiveEmity = beamVec[i].rmsEffectiveRingEmitY;
//         }
//         
// 
//         
//         if(beamVec[i].rmsRx>tempRMSSizeX)
//         {
//             tempRMSSizeX = beamVec[i].rmsRx;
//         }
//         
//         if(beamVec[i].rmsRy>tempRMSSizeY)
//         {
//             tempRMSSizeY = beamVec[i].rmsRy;
//         }
//         
//         
//          
//         if(beamVec[i].rmsEffectiveRx>tempEffectiveRMSSizeX)
//         {
//             tempEffectiveRMSSizeX = beamVec[i].rmsEffectiveRx;
//         }
//         
//         if(beamVec[i].rmsEffectiveRy>tempEffectiveRMSSizeY)
//         {
//             tempEffectiveRMSSizeY = beamVec[i].rmsEffectiveRy;
//         }
// 
// 
//         if(beamVec[i].xAver>tempAverX)
//         {
//             tempAverX = beamVec[i].xAver;
//         }
//         
//         if(beamVec[i].yAver>tempAverY)
//         {
//             tempAverY = beamVec[i].yAver;
//         }
// 
//     }
//     
//     emitXMax = tempEmitx;
//     emitYMax = tempEmity;
//     effectivEemitXMax = tempEffectiveEmitx;
//     effectivEemitYMax = tempEffectiveEmity;
//     
//     bunchEffectiveSizeXMax = tempEffectiveRMSSizeX;
//     bunchEffectiveSizeYMax = tempEffectiveRMSSizeY;
//     
//     
//     bunchSizeXMax = tempRMSSizeX;
//     bunchSizeYMax = tempRMSSizeY;
//     
//     bunchAverXMax  = tempAverX;
//     bunchAverYMax  = tempAverY;
//     
// }

void Beam::WSGetMaxBunchInfo()
{ 
    double totBunchNum = beamVec.size();
    double tempJx;
    double tempJy;

    double tempMaxAverX;
    double tempMaxAverY;
    
    tempJx      = beamVec[0].actionJx;
    tempJy      = beamVec[0].actionJy;

    tempMaxAverX = beamVec[0].xAver;
    tempMaxAverY = beamVec[0].yAver;

    
    for(int i=1;i<totBunchNum;i++)
    {
        if(beamVec[i].actionJx>tempJx)
        {
            tempJx = beamVec[i].actionJx;
        }
        
        if(beamVec[i].actionJy>tempJy)
        {
            tempJy = beamVec[i].actionJy;
        }

        
        if(beamVec[i].xAver>tempMaxAverX)
        {
            tempMaxAverX = beamVec[i].xAver;
        }
        
        if(beamVec[i].yAver>tempMaxAverY)
        {
            tempMaxAverY = beamVec[i].yAver;
        }
    }
    
    
    weakStrongBunchInfo->actionJyMax = tempJy;
    weakStrongBunchInfo->actionJxMax = tempJx;
    weakStrongBunchInfo->bunchAverXMax = tempMaxAverX;
    weakStrongBunchInfo->bunchAverYMax = tempMaxAverY;
    
    weakStrongBunchInfo->bunchEffectiveSizeXMax = beamVec[0].rmsRx + tempMaxAverX;
    weakStrongBunchInfo->bunchEffectiveSizeYMax = beamVec[0].rmsRy + tempMaxAverY;
        
    double bunchRmsSizeXTemp=0.E0;
    double bunchRmsSizeYTemp=0.E0;
    
    for(int i=0;i<totBunchNum;i++)
    {
        bunchRmsSizeXTemp += pow(beamVec[i].xAver,2);
        bunchRmsSizeYTemp += pow(beamVec[i].yAver,2);
    }
    
    weakStrongBunchInfo->bunchRmsSizeX = sqrt(bunchRmsSizeXTemp)/totBunchNum;
    weakStrongBunchInfo->bunchRmsSizeY = sqrt(bunchRmsSizeYTemp)/totBunchNum;   
    
    weakStrongBunchInfo->bunchZ  = beamVec[0].ePositionZ[0];
    weakStrongBunchInfo->bunchZP = beamVec[0].eMomentumZ[0];   

}


void Beam::FIRBunchByBunchFeedback(FIRFeedBack &firFeedBack,int nTurns)
{

    //y[0] = \sum_0^{N} a_k x[-k] 
    //Ref. Nakamura's paper spring 8 notation here used is the same with Nakamura's paper

    double eta  = firFeedBack.kickerDisp;
    double etaP = firFeedBack.kickerDispP;


    vector<double> tranAngleKicky;
    vector<double> tranAngleKickx;
    vector<double> energyKickU;

    tranAngleKicky.resize(beamVec.size());
    tranAngleKickx.resize(beamVec.size());
    energyKickU.resize(beamVec.size());
    
    int nTaps =  firFeedBack.firCoeffx.size();    
    
    vector<double> posxDataTemp(beamVec.size());
    vector<double> posyDataTemp(beamVec.size());
    vector<double> poszDataTemp(beamVec.size());

    for(int i=0;i<beamVec.size();i++)
    {   
        posxDataTemp[i] = beamVec[i].xAver;
        posyDataTemp[i] = beamVec[i].yAver;
        poszDataTemp[i] = beamVec[i].zAver;        
    }
    
    firFeedBack.posxData.erase(firFeedBack.posxData.begin());
    firFeedBack.posyData.erase(firFeedBack.posyData.begin());
    firFeedBack.poszData.erase(firFeedBack.poszData.begin());
    
    firFeedBack.posxData.push_back(posxDataTemp);
    firFeedBack.posyData.push_back(posyDataTemp);    
    firFeedBack.poszData.push_back(poszDataTemp);


    
            


    for(int i=0;i<beamVec.size();i++)
    {
        tranAngleKickx[i] = 0.E0;
        tranAngleKicky[i] = 0.E0;
        energyKickU[i]    = 0.E0;
    }

    
    for(int i=0;i<beamVec.size();i++)
    {
        for(int k=0;k<nTaps;k++)
        {
            tranAngleKickx[i] +=  firFeedBack.firCoeffx[nTaps-k-1] * firFeedBack.posxData[k][i];
            tranAngleKicky[i] +=  firFeedBack.firCoeffy[nTaps-k-1] * firFeedBack.posyData[k][i];
            energyKickU[i]    +=  firFeedBack.firCoeffz[nTaps-k-1] * firFeedBack.poszData[k][i]; 
        }        
        
        tranAngleKickx[i] =  tranAngleKickx[i] *  firFeedBack.gain * firFeedBack.kickStrengthKx;
        tranAngleKicky[i] =  tranAngleKicky[i] *  firFeedBack.gain * firFeedBack.kickStrengthKy;
        energyKickU[i]    =  energyKickU[i]    *  firFeedBack.gain * firFeedBack.kickStrengthF;
        
        
        if(tranAngleKickx[i]  > firFeedBack.fIRBunchByBunchFeedbackKickLimit)
        {
            tranAngleKickx[i] = firFeedBack.fIRBunchByBunchFeedbackKickLimit;
        }
        if(tranAngleKicky[i]  > firFeedBack.fIRBunchByBunchFeedbackKickLimit)
        {
            tranAngleKicky[i] = firFeedBack.fIRBunchByBunchFeedbackKickLimit;
        }
        if(energyKickU[i]    > firFeedBack.fIRBunchByBunchFeedbackKickLimit)
        {
            energyKickU[i]   = firFeedBack.fIRBunchByBunchFeedbackKickLimit;
        }
        

    }

    

    for(int i=0;i<beamVec.size();i++)
    {
        for(int j=0;j<beamVec[i].macroEleNumPerBunch;j++)
        {
            beamVec[i].eMomentumX[j] = beamVec[i].eMomentumX[j] + tranAngleKickx[i];
            beamVec[i].eMomentumY[j] = beamVec[i].eMomentumY[j] + tranAngleKicky[i];
            //beamVec[i].eMomentumZ[j] = beamVec[i].eMomentumZ[j] + energyKickU[i];
        }

    }
	

	// here to add functions for a simple feedback damping ---
	/*
	for(int i=0;i<beamVec.size();i++)
	{	
		for(int j=0;j<beamVec[i].macroEleNumPerBunch;j++)
		{	
			beamVec[i].ePositionX[j] = beamVec[i].ePositionX[j] * 0.98;
			beamVec[i].ePositionY[j] = beamVec[i].ePositionY[j] * 0.98; 
		}
	}
	*/







}
