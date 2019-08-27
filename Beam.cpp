#include "Beam.h"
#include "LongImpSingalBunch.h"
#include <fstream>
#include <stdlib.h> 
#include <iostream>
#include "Faddeeva.h"
#include <time.h>
#include <string>
#include <cmath>
#include <gsl/gsl_fft_complex.h>
#include <complex>
#include "FIRFeedBack.h"




using namespace std;
using std::vector;
using std::complex;


Beam::Beam()
{

}

Beam::~Beam()
{
    
}

void Beam::Initial(Train &train, LatticeInterActionPoint &latticeInterActionPoint , ReadInputSettings &inputParameter)
{
    int totBunchNum = inputParameter.totBunchNumber;
    beamVec.resize(totBunchNum);
    

    printInterval=inputParameter.printInterval;

    bunchInfoOneTurn.resize(totBunchNum);
    for(int i=0;i<totBunchNum;i++)
    {
        bunchInfoOneTurn[i].resize(15);
    }
    
    vector<int> bunchGapIndexTemp;
    bunchGapIndexTemp.resize(totBunchNum);
    
    
    actionJxMax   = 0.E0;
    actionJyMax   = 0.E0;
    bunchSizeXMax = 0.E0;
    bunchSizeYMax = 0.E0;
 
    synchRadDampTime.resize(3);
  // in the unit of number of truns for synchRadDampTime setting.
  
    for(int i=0;i<synchRadDampTime.size();i++)
    {
        synchRadDampTime[i] = inputParameter.synchRadDampTime[i];
    }



    for(int i=0;i<totBunchNum;i++)
    {
        beamVec[i].Initial(latticeInterActionPoint, inputParameter);
        beamVec[i].DistriGenerator(latticeInterActionPoint,i);
    }



    int counter=0;
    for(int i=0;i<train.trainNumber;i++)
    {
        for(int j=0;j<train.bunchNumberPerTrain[i];j++)
        {
            bunchGapIndexTemp[counter] = train.trainStart[i] + j;
//            cout<<counter<<"    "<< bunchGapIndexTemp[counter]<<endl;
            counter++; 
        }
    }    
    
    int harmonics = inputParameter.harmonics;

    for(int i=0;i<totBunchNum-1;i++)
    {
        beamVec[i].bunchGap = bunchGapIndexTemp[i+1]-bunchGapIndexTemp[i];
    }
    beamVec[totBunchNum-1].bunchGap = harmonics - bunchGapIndexTemp[totBunchNum-1];
    
//    for(int i=0;i<totBunchNum;i++)
//    {
//        cout<<beamVec[i].bunchGap<<"    "<<i<<endl;
//    }
//    getchar();
}   




void Beam::Run(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter)
{
    int nTurns = inputParameter.nTurns;
    int totBunchNum = inputParameter.totBunchNumber;
    int calSettings=inputParameter.calSettings;
    int synRadDampingFlag = inputParameter.synRadDampingFlag;
    int intevalofTurnsIonDataPrint = inputParameter.intevalofTurnsIonDataPrint;
    int fIRBunchByBunchFeedbackFlag = inputParameter.fIRBunchByBunchFeedbackFlag;
    
    ofstream fout ("result.dat",ios::out);


//   calSettings =1 :  beam_ion_interaction      --- multi-bunch-multi-turn
//   calSettings =2 :  single_bunch_Longi_imped  --- signal-bunch-signal-turn  ---board band
//   calSettings =3 :  single_bunch_Trans_imped  --- signal-bunch-signal-turn  ---board band



// preapre the data for bunch-by-bunch system --------------- 

    FIRFeedBack firFeedBack;
    if(fIRBunchByBunchFeedbackFlag)
    {        
        firFeedBack.Initial(totBunchNum);
    }
// ***************end ********************************



    switch(calSettings)
    {

        case 1 :    //   calSettings =1 :  beam_ion_interaction      --- multi-bunch-multi-turn
        {
            if(beamVec[0].macroEleNumPerBunch==1)
            {
                for(int i=0;i<nTurns;i++)   
                { 
                    for (int k=0;k<inputParameter.numberofIonBeamInterPoint;k++)
                    {
                        WSBeamRMSCal(latticeInterActionPoint,k);
                        WSBeamIonEffectOneInteractionPoint(inputParameter,latticeInterActionPoint, nTurns,  k);
                        BeamTransferPerInteractionPointDueToLattice(latticeInterActionPoint,k); 
                        
                        if(intevalofTurnsIonDataPrint && (i%intevalofTurnsIonDataPrint==0))
                        {
                            WSIonDataPrint(latticeInterActionPoint, i, k);
                        }
                    }

                    WSGetMaxBunchInfo();

                    if(synRadDampingFlag)
                    {
                         BeamSynRadDamping(synchRadDampTime,latticeInterActionPoint);
                    }
                    if(fIRBunchByBunchFeedbackFlag)
                    {
                        FIRBunchByBunchFeedback(firFeedBack,i);  
                        
                    }
  
                    IonBeamDataPrintPerTurn(i,latticeInterActionPoint,fout);

                    cout<<i<<" turns  "<<"  "
                        <<latticeInterActionPoint.ionAccumuNumber[0]<<"     "
                        <<sqrt(actionJxMax)<<"    "
                        <<sqrt(actionJyMax)<<"  "
                        <<beamVec[0].xAver<<"   "
                        <<beamVec[0].yAver<<"   "
                        <<beamVec[0].pxAver<<"  "
                        <<beamVec[0].pyAver<<"  "
                        <<log10(sqrt(actionJxMax))<<"    "
                        <<log10(sqrt(actionJyMax))
                        <<endl;
                }
            }
            else   //strong-strong beam model;
            {
                
                    
                for(int i=0;i<nTurns;i++)   
                { 

                    SSBeamIonEffectOneTurn( latticeInterActionPoint, inputParameter,i);
                
                    if(intevalofTurnsIonDataPrint && (i%intevalofTurnsIonDataPrint==0))
                    {
                        SSIonDataPrint(latticeInterActionPoint, i);
                    }
                    
                    SSGetMaxBunchInfo();

                    if(synRadDampingFlag)
                    {
                         BeamSynRadDamping(synchRadDampTime,latticeInterActionPoint);
                    }
                    if(fIRBunchByBunchFeedbackFlag)
                    {
    
                        FIRBunchByBunchFeedback(firFeedBack,i);  
                    }
                
                    
                    IonBeamDataPrintPerTurn(i,latticeInterActionPoint,fout);
                    
                    cout<<i<<" turns  "<<"  "
                        <<latticeInterActionPoint.ionAccumuNumber[0]<<"     "
                        <<emitXMax<<"    "
                        <<emitYMax<<"  "
                        <<beamVec[0].xAver<<"   "
                        <<beamVec[0].yAver<<"   "
                        <<beamVec[0].pxAver<<"  "
                        <<beamVec[0].pyAver<<"  "
                        <<endl;
                
                }
            }

            break;
        }
        case 2 :        //   calSettings =2 :  single_bunch_Longi_imped  --- signal-bunch-signal-turn
        {
            if(train.totBunchNum!=1)
            {

                cerr<<"The total buncher number is largaer than 1 in the single_bunch_Longi_imped calcualtion"<<endl;
                exit(0);
            }
            
            LongImpSingalBunch longImpSingalBunch;
            longImpSingalBunch.Initial();
            
            for(int i=0;i<nTurns;i++)   
            { 
		        BeamTransferPerTurnDueToLattice(latticeInterActionPoint);
                SingleBunchLongiImpedInterAction(longImpSingalBunch,i,fout);
		        BeamTransferPerTurnDueToLattice(latticeInterActionPoint);
            }
            break;
        }

        case 3 :
        {
            if(train.totBunchNum!=1)
            {
                cerr<<"The total bucher nummber is largaer than 1"<<endl;
                cerr<<" in the single_bunch_trans_imped calcualtion "<<endl;
                exit(0);
            }
            
            for(int i=0;i<nTurns;i++)
            { 
                cout<<"still ongoing..."<<endl;
            }
            break;
        }
        default:
            cerr<<"Do Nothing...  "<<endl;
    }



    fout.close();
}


void Beam::BeamTransferPerInteractionPointDueToLattice(LatticeInterActionPoint &latticeInterActionPoint, int k)
{
	for(int j=0;j<beamVec.size();j++)
	{
		beamVec[j].BunchTransferDueToLattice(latticeInterActionPoint,k);
	}
	
}

void Beam::BeamSynRadDamping(vector<double> &synchRadDampTime,LatticeInterActionPoint &latticeInterActionPoint)
{
	for(int j=0;j<beamVec.size();j++)
	{	
		beamVec[j].BunchSynRadDamping(synchRadDampTime,latticeInterActionPoint);		
	}
}

void Beam::WSBeamRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].WSRMSCal(latticeInterActionPoint,k);
    }
}


void Beam::BeamTransferPerTurnDueToLattice(LatticeInterActionPoint &latticeInterActionPoint)
{
    for(int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)
    {
        for(int j=0;j<beamVec.size();j++)
        {
            beamVec[j].BunchTransferDueToLattice(latticeInterActionPoint,k); 
        }
    }
}




void Beam::WSBeamIonEffectOneInteractionPoint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k)
{
    int totBunchNum = beamVec.size();
    
    double corssSectionEI = inputParameter.corssSectionEI;

    for(int j=0;j<totBunchNum;j++)
    {

        latticeInterActionPoint.ionLineDensity[k] = corssSectionEI * latticeInterActionPoint.vacuumPressure[k]
                                /latticeInterActionPoint.temperature[k]/Boltzmann * beamVec[j].electronNumPerBunch;
                                
        latticeInterActionPoint.ionNumber[k] = latticeInterActionPoint.ionLineDensity[k] * 
                                               latticeInterActionPoint.interactionLength[k];



        beamVec[j].WSRMSCal(latticeInterActionPoint,k);
       
        latticeInterActionPoint.IonGenerator(beamVec[j].rmsRx,beamVec[j].rmsRy,beamVec[j].xAver,beamVec[j].yAver,k);

        latticeInterActionPoint.IonsUpdate(k);
        latticeInterActionPoint.IonRMSCal(k);


        beamVec[j].WSIonBunchInteraction(latticeInterActionPoint,k);
        
        
        latticeInterActionPoint.IonTransferDueToBunch(beamVec[j].bunchGap,k);     
      
        beamVec[j].BunchTransferDueToIon(latticeInterActionPoint,k); 

        if(nTurns%printInterval==0 && k==0)
        {
            bunchInfoOneTurn[j][0] = nTurns;									    //1
            bunchInfoOneTurn[j][1] = j;												//2
            bunchInfoOneTurn[j][2] = latticeInterActionPoint.ionAccumuNumber[k];   //3
            bunchInfoOneTurn[j][3] = beamVec[j].rmsRx;								//4
            bunchInfoOneTurn[j][4] = beamVec[j].rmsRy;								//5
            bunchInfoOneTurn[j][5] = log10(sqrt(beamVec[j].actionJx));				//6
            bunchInfoOneTurn[j][6] = log10(sqrt(beamVec[j].actionJy));				//7
            bunchInfoOneTurn[j][7] = beamVec[j].xAver;								//8
            bunchInfoOneTurn[j][8] = beamVec[j].yAver;								//9    
            bunchInfoOneTurn[j][9] = beamVec[j].pxAver;								//10
            bunchInfoOneTurn[j][10] = beamVec[j].pyAver;							//11
            bunchInfoOneTurn[j][11] = latticeInterActionPoint.ionAccumuAverX[k];	//12
            bunchInfoOneTurn[j][12] = latticeInterActionPoint.ionAccumuAverY[k];	//13
            bunchInfoOneTurn[j][13] = latticeInterActionPoint.ionAccumuRMSX[k];     //14
            bunchInfoOneTurn[j][14] = latticeInterActionPoint.ionAccumuRMSY[k];		//15
        }           
    }     
}






void Beam::SSBeamIonEffectOneTurn( LatticeInterActionPoint &latticeInterActionPoint, ReadInputSettings &inputParameter, int nTurns)
{

    double corssSectionEI = inputParameter.corssSectionEI;
    int totBunchNum = inputParameter.totBunchNumber;


    int counter=nTurns*latticeInterActionPoint.numberOfInteraction*totBunchNum;
    
    for(int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)  
    {
        for(int j=0;j<totBunchNum;j++)
        {

            latticeInterActionPoint.ionLineDensity[k] = corssSectionEI * latticeInterActionPoint.vacuumPressure[k]
                                    /latticeInterActionPoint.temperature[k]/Boltzmann * beamVec[j].electronNumPerBunch;
            //(1) get the ion beam line density at the kth interaction points due to the jth beam bunch


            latticeInterActionPoint.ionNumber[k] = latticeInterActionPoint.ionLineDensity[k] * 
                                                   latticeInterActionPoint.interactionLength[k];
            //(2) get the ion beam number to be generated at the kth interaction point due to the jth beam bunch


            beamVec[j].RMSCal(latticeInterActionPoint,k);

            //(3) get the rms size of jth beam at kth interaction point, rmsRx rmsRy
            latticeInterActionPoint.IonGenerator(beamVec[j].rmsRx,beamVec[j].rmsRy,beamVec[j].xAver,beamVec[j].yAver,k);
            //(4) generate the ions randomly (Gausion distribution truncated 3*rmsRx) at the kth interaction point 

            latticeInterActionPoint.IonsUpdate(k);
            latticeInterActionPoint.IonRMSCal(k);
            //(5) get the accumulated ions information at the kth point (new ions plus old ions)
            //    get the rms vaues of the accumulated ion distribution.

            beamVec[j].SSIonBunchInteraction(latticeInterActionPoint,k);
            //(6) get the interactionforces between accumulated ions and the jth beam bunch at the kth interaction point 


            beamVec[j].BunchTransferDueToIon(latticeInterActionPoint,k); 
            //(7) kick the jth electron beam bunch due to accumulated ions at the kth interactoin point 


            latticeInterActionPoint.IonTransferDueToBunch(beamVec[j].bunchGap,k);         
            //(8) kick the accumulated ions due to jth electron bunch 
            
            beamVec[j].BunchTransferDueToLattice(latticeInterActionPoint,k); 
            //(9) electron beam transfer to next interaction point 


            if(nTurns%printInterval==0 && k==0)
            {
                bunchInfoOneTurn[j][0] = nTurns;									    //1
                bunchInfoOneTurn[j][1] = j;												//2
                bunchInfoOneTurn[j][2] = latticeInterActionPoint.ionAccumuNumber[k];   //3
                bunchInfoOneTurn[j][3] = beamVec[j].rmsRx;								//4
                bunchInfoOneTurn[j][4] = beamVec[j].rmsRy;								//5
                bunchInfoOneTurn[j][5] = beamVec[j].rmsEmitX;	            			//6
                bunchInfoOneTurn[j][6] = beamVec[j].rmsEmitY;				            //7
                bunchInfoOneTurn[j][7] = beamVec[j].xAver;								//8
                bunchInfoOneTurn[j][8] = beamVec[j].yAver;								//9    
                bunchInfoOneTurn[j][9] = beamVec[j].pxAver;								//10
                bunchInfoOneTurn[j][10] = beamVec[j].pyAver;							//11
                bunchInfoOneTurn[j][11] = latticeInterActionPoint.ionAccumuAverX[k];	//12
                bunchInfoOneTurn[j][12] = latticeInterActionPoint.ionAccumuAverY[k];	//13
                bunchInfoOneTurn[j][13] = latticeInterActionPoint.ionAccumuRMSX[k];     //14
                bunchInfoOneTurn[j][14] = latticeInterActionPoint.ionAccumuRMSY[k];		//15
            }           


            counter++;
        }

    }
    
    
}





void Beam::SingleBunchLongiImpedInterAction(LongImpSingalBunch &longImpSingalBunch, int nTurns, ofstream &fout)
{

    cout<<"do the beam_impedance in fre domain"<<endl;
    cout<<"only the beamVec[0] is used in calculation" <<endl;

    
    beamVec[0].InitialBeamCurDenZProf(); 




    for(int i=0;i<longImpSingalBunch.longImpedR.size();i++)
    {
//        beamVec[0];

    }


}





void Beam::SSBunchDataPrint( Bunch &bunch, int counter)
{

    char fname[256];
    
    sprintf(fname, "SSbeamdis_%d.dat", counter);

    ofstream fout(fname);


    for(int i=0;i<bunch.macroEleNumPerBunch;i++)
    {

        fout<<i                     <<"      "
            <<bunch.ePositionX[i]   <<"      "       
            <<bunch.ePositionY[i]   <<"      "       
            <<bunch.eMomentumX[i]   <<"      "       
            <<bunch.eMomentumX[i]   <<"      "       
            <<bunch.eFx[i]          <<"      "       
            <<bunch.eFy[i]          <<"      "       
            <<endl;
    }
    fout.close();
}



void Beam::IonBeamDataPrintPerTurn(int turns,  LatticeInterActionPoint &latticeInterActionPoint, ofstream &fout)
{


    if(beamVec[0].macroEleNumPerBunch==1)
    {
        for (int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)
        {
            fout<<turns<<"  "																					//1
                <<k<<"	"																						//2
                <<latticeInterActionPoint.ionAccumuNumber[k]*latticeInterActionPoint.macroIonCharge[k]<<" " 	//3
                <<log10(sqrt(actionJxMax))<<"    "																//4
                <<log10(sqrt(actionJyMax))<<"    "																//5
                <<log10(sqrt(2*beamVec[0].emittanceX))<<"	"  // constant representing the beam size in x 		//6
                <<log10(sqrt(2*beamVec[0].emittanceY))<<"	"  // constant representing the beam size in y 		//7
                <<bunchSizeXMax<<"    "																			//8
                <<bunchSizeYMax<<"    "																			//9
                <<beamVec[beamVec.size()-1].xAver<<"   "														//10
                <<beamVec[beamVec.size()-1].yAver<<"   "														//11
                <<beamVec[0].xAver<<"   "														                //12
                <<beamVec[0].yAver<<"   "														                //13
                <<latticeInterActionPoint.ionAccumuAverX[k]<<"  "												//14
                <<latticeInterActionPoint.ionAccumuAverY[k]<<"  "												//15
                <<latticeInterActionPoint.ionAccumuRMSX[k]<<"  "												//16
                <<latticeInterActionPoint.ionAccumuRMSY[k]<<"  "												//17
                <<endl;
        }
    }
    else
    {
        for (int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)
        {
            fout<<turns<<"  "																					//1
                <<k<<"	"																						//2
                <<latticeInterActionPoint.ionAccumuNumber[k]*latticeInterActionPoint.macroIonCharge[k]<<" " 	//3
                <<emitXMax<<"    "													                			//4
                <<emitXMax<<"    "												                				//5
                <<log10(sqrt(2*beamVec[0].emittanceX))<<"	"  // constant representing the beam size in x 		//6
                <<log10(sqrt(2*beamVec[0].emittanceY))<<"	"  // constant representing the beam size in y 		//7
                <<bunchSizeXMax<<"    "																			//8
                <<bunchSizeYMax<<"    "																			//9
                <<beamVec[beamVec.size()-1].xAver<<"   "														//10
                <<beamVec[beamVec.size()-1].yAver<<"   "														//11
                <<beamVec[0].xAver<<"   "														                //12
                <<beamVec[0].yAver<<"   "														                //13
                <<latticeInterActionPoint.ionAccumuAverX[k]<<"  "												//14
                <<latticeInterActionPoint.ionAccumuAverY[k]<<"  "												//15
                <<latticeInterActionPoint.ionAccumuRMSX[k]<<"  "												//16
                <<latticeInterActionPoint.ionAccumuRMSY[k]<<"  "												//17
                <<endl;
        }
    }
    

    


	// print the data on certain turns  at the first interaction point  
    if(turns%printInterval==0)
    {
        char fname[256];
        sprintf(fname, "Aver_%d.dat", turns);
        ofstream fAverOut(fname,ios::out);
        
        double dataX[2*bunchInfoOneTurn.size()];
        double dataY[2*bunchInfoOneTurn.size()];
          
        
        for(int i=0;i<bunchInfoOneTurn.size();i++)
        {
        	dataX[2*i  ] = bunchInfoOneTurn[i][10];	 
        	dataX[2*i+1] = 0.E0;
        	dataY[2*i  ] = bunchInfoOneTurn[i][11];	 
        	dataY[2*i+1] = 0.E0;
        }
        
 		gsl_fft_complex_wavetable * wavetable;
  		gsl_fft_complex_workspace * workspace;
 		
 		wavetable = gsl_fft_complex_wavetable_alloc (bunchInfoOneTurn.size());
  		workspace = gsl_fft_complex_workspace_alloc (bunchInfoOneTurn.size());
        
  		gsl_fft_complex_forward (dataX, 1, bunchInfoOneTurn.size(), wavetable, workspace);   // FFT for bunch[i].xAver in one turn;
  		        
        vector<complex<double> > beamCentoidXFFT;
        beamCentoidXFFT.resize(bunchInfoOneTurn.size());
        
        for(int i=0;i<bunchInfoOneTurn.size();i++)
        {
        	beamCentoidXFFT[i] = complex< double > ( dataX[2*i], dataX[2*i+1]  );
        }
        
        
        vector<complex<double> > beamCentoidYFFT;
        beamCentoidYFFT.resize(bunchInfoOneTurn.size());
        
        gsl_fft_complex_forward (dataY, 1, bunchInfoOneTurn.size(), wavetable, workspace);  // FFT for bunch[i].yAver in one turn;
        
        for(int i=0;i<bunchInfoOneTurn.size();i++)
        {
            beamCentoidYFFT[i] = complex< double > ( dataY[2*i], dataY[2*i+1]  );
        }
       
        
        
        for(int i=0;i<bunchInfoOneTurn.size();i++)
        {
            for(int j=0;j<bunchInfoOneTurn[0].size();j++)
            {
                fAverOut<< bunchInfoOneTurn[i][j]<<"    ";
            }
            fAverOut<<abs(beamCentoidXFFT[i])<<"    "<<abs(beamCentoidXFFT[i]);
            fAverOut<<endl;
        }

        fAverOut.close();
    }

}


void Beam::SSIonDataPrint(LatticeInterActionPoint &latticeInterActionPoint, int nTurns)
{
    for(int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)
    {
        char fname[256];
        sprintf(fname, "SSionAccumudis_%d_%d.dat", nTurns,k);
        ofstream fout(fname);

        for(int i=0;i<latticeInterActionPoint.ionAccumuNumber[k];i++)
        {
            fout<<i                     <<"      "
                <<latticeInterActionPoint.ionAccumuPositionX[k][i]   <<"      "       
                <<latticeInterActionPoint.ionAccumuPositionY[k][i]   <<"      "       
                <<latticeInterActionPoint.ionAccumuVelocityX[k][i]   <<"      "       
                <<latticeInterActionPoint.ionAccumuVelocityY[k][i]   <<"      "       
                <<latticeInterActionPoint.ionAccumuFx[k][i]   <<"      "       
                <<latticeInterActionPoint.ionAccumuFy[k][i]   <<"      "       
                <<endl;
        }
        
        fout.close();
        
    }
}


void Beam::WSIonDataPrint(LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k)
{

    char fname[256];
    
    sprintf(fname, "WSionAccumudis_%d_%d.dat", nTurns, k);

    ofstream fout(fname);

    for(int i=0;i<latticeInterActionPoint.ionAccumuNumber[k];i++)
    {
        fout<<i                     <<"      "
            <<latticeInterActionPoint.ionAccumuPositionX[k][i]   <<"      "       
            <<latticeInterActionPoint.ionAccumuPositionY[k][i]   <<"      "       
            <<latticeInterActionPoint.ionAccumuVelocityX[k][i]   <<"      "       
            <<latticeInterActionPoint.ionAccumuVelocityY[k][i]   <<"      "       
            <<latticeInterActionPoint.ionAccumuFx[k][i]   <<"      "       
            <<latticeInterActionPoint.ionAccumuFy[k][i]   <<"      "       
            <<endl;
    }
    
    fout.close();
}

void Beam::SSGetMaxBunchInfo()
{
    double totBunchNum = beamVec.size();
    double tempEmitx;
    double tempEmity;
    double tempSizeX;
    double tempSizeY;
    
    tempEmitx      = beamVec[0].rmsEmitX;
    tempEmity      = beamVec[0].rmsEmitY;
    tempSizeX      = beamVec[0].rmsRx;
    tempSizeY      = beamVec[0].rmsRy;
    

    
    for(int i=1;i<totBunchNum;i++)
    {
        if(beamVec[i].rmsEmitX>tempEmitx)
        {
            tempEmitx = beamVec[i].rmsEmitX;
        }
        
        if(beamVec[i].rmsEmitY>tempEmity)
        {
            tempEmity = beamVec[i].rmsEmitY;
        }
        
        if(beamVec[i].rmsRx>tempSizeX)
        {
            tempSizeX = beamVec[i].rmsRx;
        }
        
        if(beamVec[i].rmsRy>tempSizeY)
        {
            tempSizeY = beamVec[i].rmsRy;
        }
        
    }
    
    

    
    emitXMax = tempEmitx;
    emitYMax = tempEmity;
    bunchSizeXMax = tempSizeX;
    bunchSizeYMax = tempSizeY;
}

void Beam::WSGetMaxBunchInfo()
{

    double totBunchNum = beamVec.size();
    double tempJx;
    double tempJy;
    double tempSizeX;
    double tempSizeY;
    
    tempJx      = beamVec[0].actionJx;
    tempJy      = beamVec[0].actionJy;
    tempSizeX   = beamVec[0].rmsRx;
    tempSizeY   = beamVec[0].rmsRy;
    

    
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
        
        if(beamVec[i].rmsRx>tempSizeX)
        {
            tempSizeX = beamVec[i].rmsRx;
        }
        
        if(beamVec[i].rmsRy>tempSizeY)
        {
            tempSizeY = beamVec[i].rmsRy;
        }
        
    }
    
    
    actionJyMax = tempJy;
    actionJxMax = tempJx;
    bunchSizeXMax = tempSizeX;
    bunchSizeYMax = tempSizeY;
    

}


void Beam::FIRBunchByBunchFeedback(FIRFeedBack &firFeedBack,int nTurns)
{

    //y[0] = \sum_0^{N} a_k x[k] 
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

    int tapsIndex = nTurns%nTaps;


    // started from a simpliest case, as Section IX A. Eq. 116 only transverse kick in y direction for test calculation  
    for(int i=0;i<beamVec.size();i++)
    {
        firFeedBack.posxData[tapsIndex][i] = beamVec[i].xAver;
        firFeedBack.posyData[tapsIndex][i] = beamVec[i].yAver;
        firFeedBack.poszData[tapsIndex][i] = beamVec[i].zAver;
    }


    for(int i=0;i<beamVec.size();i++)
    {
        tranAngleKicky[i] = 0.E0;
        tranAngleKickx[i] = 0.E0;
        energyKickU[i]    = 0.E0;
    }

    for(int i=0;i<beamVec.size();i++)
    {
        for(int k=0;k<nTaps;k++)
        {
            tranAngleKickx[i] +=  firFeedBack.firCoeffx[nTaps-k-1] * firFeedBack.posxData[k][i];
            tranAngleKicky[i] +=  firFeedBack.firCoeffy[nTaps-k-1] * firFeedBack.posyData[k][i];
            energyKickU[i]    +=  firFeedBack.firCoeffz[nTaps-k-1] * firFeedBack.poszData[k][i]; 
//            if(i=beamVec.size()-cout)
//                1<<k<<"  "<<nTaps-k-1<<"   "<<firFeedBack.firCoeffx[nTaps-k-1] <<"	"<<firFeedBack.posxData[k][i]<<endl;
            
        }
        
//        getchar();

        tranAngleKickx[i] =  tranAngleKickx[i] * firFeedBack.kickStrengthKx;
        tranAngleKicky[i] =  tranAngleKicky[i] * firFeedBack.kickStrengthKy;	
        energyKickU[i]    =  energyKickU[i]    * firFeedBack.kickStrengthF;
    }



    for(int i=0;i<beamVec.size();i++)
    {
        for(int j=0;j<beamVec[i].macroEleNumPerBunch;j++)
        {
            beamVec[i].eMomentumX[j] = beamVec[i].eMomentumX[j] + tranAngleKickx[i];
            beamVec[i].eMomentumY[j] = beamVec[i].eMomentumY[j] + tranAngleKicky[i];
//          beamVec[i].eMomentumZ[j] = beamVec[i].eMomentumZ[j] + energyKickU[i];
        }

    }


//    char fname[256];
    

//	if(nTurns%100==0)
//	{
//	    sprintf(fname, "bunch_feed%d.dat", nTurns);
//    	ofstream fFBTestOut(fname);
//	
//		for(int i=0;i<beamVec.size();i++)
//		{
//			fFBTestOut<<i<<"	"<< beamVec[i].eMomentumX[i] <<"	"<<beamVec[i].eMomentumY[i]<<"	"<<beamVec[i].eMomentumZ[i]<<endl;
//		}
//		fFBTestOut.close();
//	}
	
	
	
}



