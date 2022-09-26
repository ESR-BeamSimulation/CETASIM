//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             


#include "FittingGSL.h"
#include "FIRFeedBack.h"
#include "MPBeam.h"
#include "Faddeeva.h"
#include "WakeFunction.h"
//#include "LongImpSingalBunch.h"
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <string>
#include <cmath>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <complex>
#include <iomanip>
#include <cstring>




using namespace std;
using std::vector;
using std::complex;


MPBeam::MPBeam()
{
}

MPBeam::~MPBeam()
{
      delete strongStrongBunchInfo;
}


void MPBeam::Initial(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter)
{
    
    int totBunchNum = inputParameter.ringFillPatt->totBunchNumber; 
    int harmonics = inputParameter.ringParBasic->harmonics;              
    beamVec.resize(totBunchNum);
    timeBetweenBunch.resize(totBunchNum);

    // set bunch harmonic number    
    int counter=0;
    for(int i=0;i<train.trainNumber;i++)
    {
        for(int j=0;j<train.bunchNumberPerTrain[i];j++)
        {
			beamVec[counter].bunchHarmNum  = train.trainStart[i] + j * (train.bunchGaps[i] + 1);
            counter++;
        }
    }
    
    // set the bunch initial distribution
    for(int i=0;i<totBunchNum;i++)
    {
        beamVec[i].InitialMPBunch(inputParameter);
        beamVec[i].DistriGenerator(latticeInterActionPoint,inputParameter,i);
    }
 
    
    for(int i=0;i<totBunchNum-1;i++)
    {
        beamVec[i].bunchGap = beamVec[i+1].bunchHarmNum - beamVec[i].bunchHarmNum;
    }
    beamVec[totBunchNum-1].bunchGap = harmonics - beamVec[totBunchNum-1].bunchHarmNum;

	RMOutPutFiles();
	     
}

void MPBeam::InitialcavityResonator(ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    // According the beam filling pattern, tracking from transient to get steady state Vb. Ref. Bill Chapter 7.4.2
    // The notatoion of Bill's book is the same in P.B.Wilson's Section 3. Slac Pub 6062  except li -> - li, the phase rotation oppsiste.

    int resNum      = inputParameter.ringParRf->resNum;
    int ringHarmH   = inputParameter.ringParBasic->harmonics;
    double beamCurr = inputParameter.ringBunchPara->current * beamVec.size();
    double f0       = inputParameter.ringParBasic->f0;

    double deltaL;
    double cPsi;
    double tB;
    double tF;
    complex<double> vb0 = (0.E0,0.E0);
    complex<double> vb0temp = (0.E0,0.E0);
    complex<double> vbAccum = (0.E0,0.E0);
    complex<double> vbAccumAver = (0.E0,0.E0);
    complex<double> vbAccumAver1 = (0.E0,0.E0);
    complex<double> vbAccumAver2 = (0.E0,0.E0);
    complex<double> vbKickAccum = (0.E0,0.E0);
    complex<double> vbKickAver  = (0.E0,0.E0);

    int nTurns;
    int bunchHarmIndex;


    
    for(int i=0; i<resNum;i++)
    {

        tF     = cavityResonator.resonatorVec[i].tF  ;      // [s]
        nTurns = ceil(100 * tF * f0);                       // 100 cavity field filling/damping time defautl

        vbAccum      =   complex<double>(0,0);
        vbAccumAver1 =   complex<double>(0,0);
        vbAccumAver2 =   complex<double>(0,0);

        vbKickAccum  =   complex<double>(0,0);  // get the accumme of  vb0/2
        vbKickAver   =   complex<double>(0,0);

        for(int n=0;n<nTurns;n++)
        {
            int k=0;
            for(int j=0;j<ringHarmH;j++)
            {

                if( k<beamVec.size() && beamVec[k].bunchHarmNum == j)
                {
                    vb0  = complex<double>(-1 * cavityResonator.resonatorVec[i].resFre * 2 * PI * cavityResonator.resonatorVec[i].resShuntImpRs
                                              /  cavityResonator.resonatorVec[i].resQualityQ0, 0.E0)
                                              *  beamVec[k].electronNumPerBunch * ElectronCharge;       
                                              // wilson's equation p.6, The same as Bill equation (7.76), cos convention, real part reresents the energy gain. 
                    vb0temp = vb0;
                    k++;

                    if( n==nTurns-1)
                    {
                        vbKickAccum += vb0temp/2.0 ;
                    }

                }
                else
                {
                    vb0 = complex<double> (0.E0,0.E0);
                }

                vbAccum += vb0;                                     //[V]

                tB = 1.0 / f0 / ringHarmH ;                         //[s]
                deltaL = tB / tF ;
                cPsi   = deltaL * tan(cavityResonator.resonatorVec[i].resDeTunePsi );

                if( n==nTurns-1)
                {
                    vbAccumAver1 +=  vbAccum;
                }

                vbAccum = vbAccum * exp(- deltaL ) * exp (li * cPsi);    // decay and rotate...  [V] P.B.  Eq. (3.12) -- only get the beam induced voltage at the cavity

                if( n==nTurns-1)
                {
                    vbAccumAver2 +=  vbAccum;
                }
            }
        }
        vbAccumAver =  (vbAccumAver1 + vbAccumAver2) /2.0/ double(ringHarmH);   // (1) average of bunch induced voltage
        vbAccumAver =  vbAccumAver2/double(ringHarmH);                          // (2) after decay and rotate, just before next bunch.
        
        cavityResonator.resonatorVec[i].vbAccum   =  vbAccumAver;               // set cavity beam induced voltage
        vbKickAver  =  vbKickAccum /double(beamVec.size());                     // average energy that bunch is kicked by self-loss Vb0/2. 

        // extra feedback condition to compensate the self-loss term vb0/2 
        //double genAddvbArg = arg(cavityResonator.resonatorVec[i].resCavVolReq);       
        //if(genAddvbArg==PI/2.0 or genAddvbArg==-PI/2.0)
        //{
        //    cerr<<"Required Cavity Phase is PI/2 or -PI/2"<<endl;
        //    cerr<<"Does not work right now when cavity feedback is included, because of the scheme of feedback by amplitude"<<endl;
        //    exit(0);
        //}
        //double genAddvbAbs = (cavityResonator.resonatorVec[i].resCavVolReq.real() - vbKickAver.real()) / cos( genAddvbArg );
        //genAddvbAbs = abs(genAddvbAbs);
        //cavityResonator.resonatorVec[i].resGenVol = genAddvbAbs * exp(li * genAddvbArg ) - vbAccumAver;
        
        cavityResonator.resonatorVec[i].resGenVol = cavityResonator.resonatorVec[i].resCavVolReq - vbAccumAver;

        //set the cold and warm cavity condition. Vb is set up before particle tracking in the "warm" condition
        if(cavityResonator.resonatorVec[i].resCold)
        {
            cavityResonator.resonatorVec[i].vbAccum = complex<double>(0.E0,0.E0);
        }

        // set the active and passive condition.
        if (cavityResonator.resonatorVec[i].resType==0)
        {
            cavityResonator.resonatorVec[i].resGenVol=complex<double>(0.E0,0.E0);
        }

        // set initial cavity voltage, beam indcued volage for each bunch         
        for(int j=0;j<beamVec.size();j++)
        {
            beamVec[j].cavFBCenInfo->cavVolBunchCen[i]    =     cavityResonator.resonatorVec[i].resCavVolReq;     // set vol as the requried voltage.
            beamVec[j].cavFBCenInfo->cavAmpBunchCen[i]    = abs(cavityResonator.resonatorVec[i].resCavVolReq);
            beamVec[j].cavFBCenInfo->cavPhaseBunchCen[i]  = arg(cavityResonator.resonatorVec[i].resCavVolReq);
            beamVec[j].cavFBCenInfo->induceVolBunchCen[i] =     cavityResonator.resonatorVec[i].vbAccum;
            beamVec[j].cavFBCenInfo->genVolBunchAver[i]   =     cavityResonator.resonatorVec[i].resGenVol;
            beamVec[j].cavFBCenInfo->selfLossVolBunchCen[i] = vbKickAver;
        }
        
    }

    // print out the data to show the cavity voltage buildup process. Generated voltage compansate the Vb0/2 already by amplitude feedback. 
    vb0=(0,0);
    vbAccum=(0,0);

    ofstream fout(inputParameter.ringParRf->transResonParWriteTo+".sdds");
	fout<<"SDDS1"<<endl;
    fout<<"&parameter name=CavAmpIdeal,      units=V,   type=float,  &end"<<endl;
	fout<<"&parameter name=CavPhaseIdeal,    units=rad, type=float,  &end"<<endl;
	fout<<"&parameter name=CavFreq,          units=Hz,  type=float,  &end"<<endl;
	fout<<"&parameter name=GenAmp,           units=V,   type=float,  &end"<<endl;
	fout<<"&parameter name=GenPhase,         units=rad, type=float,  &end"<<endl;
    fout<<"&parameter name=beamIndAmp,       units=V,   type=float,  &end"<<endl;
	fout<<"&parameter name=beamIndPhase,     units=rad, type=float,  &end"<<endl;
    fout<<"&parameter name=CavDetunPsi,      units=rad, type=float,  &end"<<endl;

	fout<<"&column name=Turns,              units=s       type=float,  &end"<<endl;
	fout<<"&column name=CavAmpReq,          units=V,      type=float,  &end"<<endl;
    fout<<"&column name=CavPhaseReq,        units=rad,    type=float,  &end"<<endl;
    fout<<"&column name=CavAmp,             units=V,      type=float,  &end"<<endl;
    fout<<"&column name=CavPhase,           units=rad,    type=float,  &end"<<endl;
	fout<<"&column name=BeamIndAmp,         units=V,      type=float,  &end"<<endl;
    fout<<"&column name=BeamIndPhase,       units=rad,    type=float,  &end"<<endl;
    fout<<"&column name=GenAmp,             units=V,      type=float,  &end"<<endl;
    fout<<"&column name=GenPhase,           units=rad,    type=float,  &end"<<endl;
	fout<<"&data mode=ascii, &end"<<endl;


    double time=0;
    for(int i=0; i<resNum;i++)
    {
        time =0;
        tF     = cavityResonator.resonatorVec[i].tF  ;      // [s]
        nTurns = ceil(100 * tF * f0);
        vbAccum     =   complex<double>(0,0);
         
        fout<<abs(cavityResonator.resonatorVec[i].resCavVolReq)           <<endl;
        fout<<arg(cavityResonator.resonatorVec[i].resCavVolReq)           <<endl;
        fout<<    cavityResonator.resonatorVec[i].resFre                  <<endl;
        fout<<abs(cavityResonator.resonatorVec[i].resGenVol)              <<endl;
        fout<<arg(cavityResonator.resonatorVec[i].resGenVol)              <<endl;
        fout<<abs(cavityResonator.resonatorVec[i].vbAccum  )              <<endl;
        fout<<arg(cavityResonator.resonatorVec[i].vbAccum  )              <<endl;
        fout<<    cavityResonator.resonatorVec[i].resDeTunePsi            <<endl;

        fout<<"! page number "<<i+1<<endl;
        fout<<nTurns*beamVec.size()<<endl;


        for(int n=0;n<nTurns;n++)
        {
            for(int j=0;j<beamVec.size();j++)
            {
                vb0  = complex<double>(-1 * cavityResonator.resonatorVec[i].resFre * 2 * PI * cavityResonator.resonatorVec[i].resShuntImpRs /
                                            cavityResonator.resonatorVec[i].resQualityQ0,0.E0) * beamVec[j].electronNumPerBunch * ElectronCharge;
                vbAccum += vb0;                                       //[V]

                tB = beamVec[j].bunchGap * 1 / f0 / ringHarmH ;  //[s]

                time +=  tB;
                deltaL = tB / tF ;
                cPsi   = deltaL * tan(cavityResonator.resonatorVec[i].resDeTunePsi);


                fout<<setw(15)<<left<<time*f0 
                    <<setw(15)<<left<<abs(cavityResonator.resonatorVec[i].resCavVolReq)   // required CavVolAmp
                    <<setw(15)<<left<<arg(cavityResonator.resonatorVec[i].resCavVolReq)   // required CavVolPhase
                    <<setw(15)<<left<<abs(vbAccum + cavityResonator.resonatorVec[i].resGenVol)
                    <<setw(15)<<left<<arg(vbAccum + cavityResonator.resonatorVec[i].resGenVol)
                    <<setw(15)<<left<<abs(vbAccum)
                    <<setw(15)<<left<<arg(vbAccum)
                    <<setw(15)<<left<<abs(cavityResonator.resonatorVec[i].resGenVol)
                    <<setw(15)<<left<<arg(cavityResonator.resonatorVec[i].resGenVol)
                    <<endl;

                vbAccum = vbAccum * exp(- deltaL ) * exp (li * cPsi);    // decay and rotate...  [V]-- use the voltage after decay and feed this into tracking.
            }
        }


    }

    fout.close();       
}

void MPBeam::Run(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter, CavityResonator &cavityResonator)
{
    int nTurns                      = inputParameter.ringRun->nTurns;
    int *synRadDampingFlag           = inputParameter.ringRun->synRadDampingFlag;
    int fIRBunchByBunchFeedbackFlag = inputParameter.ringRun->fIRBunchByBunchFeedbackFlag;
    int impedanceFlag               = inputParameter.ringRun->impedanceFlag;
    int beamIonFlag                 = inputParameter.ringRun->beamIonFlag;
    int lRWakeFlag                  = inputParameter.ringRun->lRWakeFlag;
    int sRWakeFlag                  = inputParameter.ringRun->sRWakeFlag;
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
    WakeFunction lRWakeFunction;
    WakeFunction sRWakeFunction;
    if(lRWakeFlag)
    {            
        lRWakeFunction.InitialLRWake(inputParameter,latticeInterActionPoint);
    }
    if(sRWakeFlag)
    {            
        sRWakeFunction.InitialSRWake(inputParameter,latticeInterActionPoint);
    }
     
    //-----------------------------------------------------------              
    // turn by turn data-- average of bunches 
    ofstream fout ("result.sdds",ios::out);
    fout<<"SDDS1"<<endl;
    fout<<"&column name=turns,          type=long,              &end"<<endl;
    fout<<"&column name=ionCharge,      units=e,   type=float,  &end"<<endl;
    fout<<"&column name=maxAverX,       units=m,   type=float,  &end"<<endl;
    fout<<"&column name=maxAverY,       units=m,   type=float,  &end"<<endl;
    fout<<"&column name=averAllBunchX,  units=m,   type=float,  &end"<<endl;
    fout<<"&column name=averAllBunchY,  units=m,   type=float,  &end"<<endl;
    fout<<"&column name=averAllBunchZ,  units=m,   type=float,  &end"<<endl;
    fout<<"&column name=averAllBunchPX, units=rad, type=float,  &end"<<endl;
    fout<<"&column name=averAllBunchPY, units=rad, type=float,  &end"<<endl;
    fout<<"&column name=averAllBunchPZ, units=rad, type=float,  &end"<<endl;
    fout<<"&column name=rmsAllBunchX,   units=m,   type=float,  &end"<<endl;
    fout<<"&column name=rmsAllBunchY,   units=m,   type=float,  &end"<<endl;
    fout<<"&column name=rmsAllBunchZ,   units=m,   type=float,  &end"<<endl;
    fout<<"&column name=rmsAllBunchPX,  units=rad, type=float,  &end"<<endl;
    fout<<"&column name=rmsAllBunchPY,  units=rad, type=float,  &end"<<endl;
    fout<<"&column name=rmsAllBunchPZ,  units=rad, type=float,  &end"<<endl;

    for(int j=0;j<inputParameter.ringRun->TBTBunchPrintNum;j++)
    {
        string colname = string("&column name=") + string("averX_")     + to_string(j) + string(", units=m,   type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("averPx_")            + to_string(j) + string(", units=rad, type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("averY_")             + to_string(j) + string(", units=m,   type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("averPy_")            + to_string(j) + string(", units=rad, type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("averZ_")             + to_string(j) + string(", units=m,   type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("averPz_")            + to_string(j) + string(", units=rad, type=float,  &end");
        fout<<colname<<endl;
        
        colname = string("&column name=") + string("rmsBunchZ_")     + to_string(j) + string(", units=m,   type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("rmsEnergySpread_") + to_string(j) + string(", units=rad, type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("rmsRx_")         + to_string(j) + string(", units=m,   type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("rmsRy_")         + to_string(j) + string(", units=m,   type=float,  &end");            
        fout<<colname<<endl;
        colname = string("&column name=") + string("rmsEmitx_")      + to_string(j) + string(", units=m*rad,   type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("rmsEmity_")      + to_string(j) + string(", units=m*rad,   type=float,  &end");            
        fout<<colname<<endl;
    }
    fout<<"&data mode=ascii, &end"<<endl;
    fout<<"! page number "<<1<<endl;
    fout<<nTurns<<endl;



    // run loop starts, for nTrns and each trun for k interaction-points
    for(int n=0;n<nTurns;n++)
    {
        if(n%10==0) cout<<n<<"  turns"<<endl;

        MPBeamRMSCal(latticeInterActionPoint, 0);
        MPGetBeamInfo();
                
        if(beamIonFlag) 
        {
            for (int k=0;k<inputParameter.ringIonEffPara->numberofIonBeamInterPoint;k++)
            {
                MPBeamRMSCal(latticeInterActionPoint, k);
                MPGetBeamInfo();
                SSBeamIonEffectOneInteractionPoint(inputParameter,latticeInterActionPoint, n, k);
                BeamTransferPerInteractionPointDueToLatticeT(latticeInterActionPoint,k);             //transverse transfor per interaction point                
            }

            if(ionInfoPrintInterval && (n%ionInfoPrintInterval==0))
            {
                SSIonDataPrint(inputParameter,latticeInterActionPoint, n);
            }
        }
        else
        {
            MPBeamRMSCal(latticeInterActionPoint, 0);
            BeamTransferPerTurnDueToLatticeT(latticeInterActionPoint);
        }

        BeamTransferPerTurnDueToLatticeL(inputParameter,cavityResonator); 

        if(lRWakeFlag)
        {
            MPBeamRMSCal(latticeInterActionPoint, 0);
            LRWakeBeamIntaction(inputParameter,lRWakeFunction,latticeInterActionPoint);
            BeamTransferPerTurnDueWake();        
        }
        
        if(sRWakeFlag) 
        {
            MPBeamRMSCal(latticeInterActionPoint, 0);
            SRWakeBeamIntaction(inputParameter,sRWakeFunction,latticeInterActionPoint,n);         
        }
        
        
        if(fIRBunchByBunchFeedbackFlag)
        {
            MPBeamRMSCal(latticeInterActionPoint,0);
            FIRBunchByBunchFeedback(firFeedBack,n);
        }
        
        if(synRadDampingFlag)  // only in transverse direction
        {
            MPBeamRMSCal(latticeInterActionPoint,0);
            BeamSynRadDamping(inputParameter,latticeInterActionPoint);
        }

        fout<<n<<"  "
            <<setw(15)<<left<< latticeInterActionPoint.totIonCharge
            <<setw(15)<<left<< strongStrongBunchInfo->bunchAverXMax
            <<setw(15)<<left<< strongStrongBunchInfo->bunchAverYMax
            <<setw(15)<<left<< strongStrongBunchInfo->bunchAverX        // over all bunches, \sum_i beamVec[i].xAver / totBunchNumber
            <<setw(15)<<left<< strongStrongBunchInfo->bunchAverY        // over all bunches
            <<setw(15)<<left<< strongStrongBunchInfo->bunchAverZ        // over all bunches
            <<setw(15)<<left<< strongStrongBunchInfo->bunchAverPX       // over all bunches
            <<setw(15)<<left<< strongStrongBunchInfo->bunchAverPY       // over all bunches
            <<setw(15)<<left<< strongStrongBunchInfo->bunchAverPZ       // over all bunches
            <<setw(15)<<left<< strongStrongBunchInfo->bunchRmsSizeX     // over all bunches \sum_i pow(beamVec[i].xAver,2)/totBunchNumber 
            <<setw(15)<<left<< strongStrongBunchInfo->bunchRmsSizeY     // over all bunches
            <<setw(15)<<left<< strongStrongBunchInfo->bunchRmsSizeZ     // over all bunches
            <<setw(15)<<left<< strongStrongBunchInfo->bunchRmsSizePX    // over all bunches
            <<setw(15)<<left<< strongStrongBunchInfo->bunchRmsSizePY    // over all bunches
            <<setw(15)<<left<< strongStrongBunchInfo->bunchRmsSizePZ;   // over all bunches

        for(int i=0;i<inputParameter.ringRun->TBTBunchPrintNum;i++)
        {
            int index = inputParameter.ringRun->TBTBunchDisDataBunchIndex[i];
            fout<<setw(15)<<left<<beamVec[index].xAver
                <<setw(15)<<left<<beamVec[index].pxAver
                <<setw(15)<<left<<beamVec[index].yAver
                <<setw(15)<<left<<beamVec[index].pyAver
                <<setw(15)<<left<<beamVec[index].zAver
                <<setw(15)<<left<<beamVec[index].pzAver
                <<setw(15)<<left<<beamVec[index].rmsBunchLength
                <<setw(15)<<left<<beamVec[index].rmsEnergySpread
                <<setw(15)<<left<<beamVec[index].rmsRx
                <<setw(15)<<left<<beamVec[index].rmsRy
                <<setw(15)<<left<<beamVec[index].emittanceX
                <<setw(15)<<left<<beamVec[index].emittanceY;
        }
        fout<<endl;        
        
       
        if(  bunchInfoPrintInterval && (n%bunchInfoPrintInterval==0)  )
        {             
            MPBeamDataPrintPerTurn(n,latticeInterActionPoint,inputParameter);           
        }
                                           
    }
    fout.close();
    
    cout<<"End of Tracking "<<nTurns<< "Turns"<<endl;
        
}


//void MPBeam::SetBeamPosHistoryDataWithinWindow()
//{
//    for(int j=0;j<beamVec.size();j++)
//    {
//        beamVec[j].SetBunchPosHistoryDataWithinWindow();
//    }
//}

void MPBeam::SRWakeBeamIntaction(const  ReadInputSettings &inputParameter, WakeFunction &sRWakeFunction, const  LatticeInterActionPoint &latticeInterActionPoint, int turns)
{
    for(int j=0;j<beamVec.size();j++)
    {
         beamVec[j].BunchTransferDueToSRWake(inputParameter,sRWakeFunction,latticeInterActionPoint,turns);
    }
}   


void MPBeam::MPBeamRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].GetMPBunchRMS(latticeInterActionPoint,k);
    } 
}


void MPBeam::MPGetBeamInfo()
{
    double totBunchNum = beamVec.size();
    double tempEmitx;
    double tempEmity;
    double tempEmitz;
    double tempEffectiveEmitx;
    double tempEffectiveEmity;
    double tempRMSSizeX;
    double tempRMSSizeY;
    double tempEffectiveRMSSizeX;
    double tempEffectiveRMSSizeY;
    double tempBunchLength;
    double tempBunchEnergySpread;

    double tempAverX;
    double tempAverY;
    double tempAverZ;

    tempEmitx   = beamVec[0].emittanceX;
    tempEmity   = beamVec[0].emittanceY;
    tempEmitz   = beamVec[0].emittanceZ;

    tempEffectiveEmitx   = beamVec[0].rmsEffectiveRingEmitX;
    tempEffectiveEmity   = beamVec[0].rmsEffectiveRingEmitY;

    tempRMSSizeX   = beamVec[0].rmsRx;
    tempRMSSizeY   = beamVec[0].rmsRy;
    tempBunchLength = beamVec[0].rmsBunchLength;
    tempBunchEnergySpread = beamVec[0].rmsEnergySpread;

    tempEffectiveRMSSizeX = beamVec[0].rmsEffectiveRx;
    tempEffectiveRMSSizeY = beamVec[0].rmsEffectiveRy;


    tempAverX   = beamVec[0].xAver;
    tempAverY   = beamVec[0].yAver;
    tempAverZ   = beamVec[0].zAver;

    for(int i=1;i<totBunchNum;i++)
    {
        
        if(beamVec[i].emittanceX>tempEmitx)
        {
            tempEmitx = beamVec[i].emittanceX;
        }

        if(beamVec[i].emittanceY>tempEmity)
        {
            tempEmity = beamVec[i].emittanceY;
        }

        if(beamVec[i].emittanceZ>tempEmitz)
        {
            tempEmitz = beamVec[i].emittanceZ;
        }


        if(beamVec[i].rmsEffectiveRingEmitX>tempEffectiveEmitx)
        {
            tempEffectiveEmitx = beamVec[i].rmsEffectiveRingEmitX;
        }

        if(beamVec[i].rmsEffectiveRingEmitY>tempEffectiveEmity)
        {
            tempEffectiveEmity = beamVec[i].rmsEffectiveRingEmitY;
        }

        if(beamVec[i].rmsRx>tempRMSSizeX)
        {
            tempRMSSizeX = beamVec[i].rmsRx;
        }
        if(beamVec[i].rmsRy>tempRMSSizeY)
        {
            tempRMSSizeY = beamVec[i].rmsRy;
        }
        if(beamVec[i].rmsBunchLength>tempBunchLength)
        {
            tempBunchLength = beamVec[i].rmsBunchLength;
        }
        if(beamVec[i].rmsEnergySpread>tempBunchEnergySpread)
        {
            tempBunchEnergySpread = beamVec[i].rmsEnergySpread;
        }


        if(beamVec[i].rmsEffectiveRx>tempEffectiveRMSSizeX)
        {
            tempEffectiveRMSSizeX = beamVec[i].rmsEffectiveRx;
        }

        if(beamVec[i].rmsEffectiveRy>tempEffectiveRMSSizeY)
        {
            tempEffectiveRMSSizeY = beamVec[i].rmsEffectiveRy;
        }


        if(beamVec[i].xAver>tempAverX)
        {
            tempAverX = beamVec[i].xAver;
        }

        if(beamVec[i].yAver>tempAverY)
        {
            tempAverY = beamVec[i].yAver;
        }
        if(beamVec[i].zAver>tempAverZ)
        {
            tempAverZ = beamVec[i].zAver;
        }

    }

    strongStrongBunchInfo->emitXMax = tempEmitx;
    strongStrongBunchInfo->emitYMax = tempEmity;
    strongStrongBunchInfo->effectivEemitXMax = tempEffectiveEmitx;
    strongStrongBunchInfo->effectivEemitYMax = tempEffectiveEmity;

    strongStrongBunchInfo->bunchEffectiveSizeXMax = tempEffectiveRMSSizeX;
    strongStrongBunchInfo->bunchEffectiveSizeYMax = tempEffectiveRMSSizeY;


    strongStrongBunchInfo->bunchSizeXMax = tempRMSSizeX;
    strongStrongBunchInfo->bunchSizeYMax = tempRMSSizeY;
    strongStrongBunchInfo->bunchLengthMax = tempBunchLength;
    strongStrongBunchInfo->bunchEnergySpreadMax = tempBunchEnergySpread;


    strongStrongBunchInfo->bunchAverXMax  = tempAverX;
    strongStrongBunchInfo->bunchAverYMax  = tempAverY;
    strongStrongBunchInfo->bunchAverZMax  = tempAverZ;


    double bunchRmsSizeXTemp=0.e0;
    double bunchRmsSizeYTemp=0.e0;
    double bunchRmsSizeZTemp=0.e0;
    double bunchRmsSizePXTemp=0.e0;
    double bunchRmsSizePYTemp=0.e0;
    double bunchRmsSizePZTemp=0.e0;


 
    for(int i=0;i<totBunchNum;i++)
    {
        bunchRmsSizeXTemp += pow(beamVec[i].xAver,2);
        bunchRmsSizeYTemp += pow(beamVec[i].yAver,2);
        bunchRmsSizeZTemp += pow(beamVec[i].zAver,2);
	    bunchRmsSizePXTemp+= pow(beamVec[i].pxAver,2);
	    bunchRmsSizePYTemp+= pow(beamVec[i].pyAver,2);
	    bunchRmsSizePZTemp+= pow(beamVec[i].pzAver,2);
    }

    strongStrongBunchInfo->bunchRmsSizeX  = sqrt(bunchRmsSizeXTemp  /totBunchNum);
    strongStrongBunchInfo->bunchRmsSizeY  = sqrt(bunchRmsSizeYTemp  /totBunchNum);
    strongStrongBunchInfo->bunchRmsSizeZ  = sqrt(bunchRmsSizeZTemp  /totBunchNum);
    strongStrongBunchInfo->bunchRmsSizePX = sqrt(bunchRmsSizePXTemp /totBunchNum);
    strongStrongBunchInfo->bunchRmsSizePY = sqrt(bunchRmsSizePYTemp /totBunchNum);
    strongStrongBunchInfo->bunchRmsSizePZ = sqrt(bunchRmsSizePZTemp /totBunchNum);

    double bunchAverX=0;
    double bunchAverY=0;
    double bunchAverZ=0;
    double bunchAverPX=0;
    double bunchAverPY=0;
    double bunchAverPZ=0;

    for(int i=0;i<totBunchNum;i++)
    {
        bunchAverX += beamVec[i].xAver;
        bunchAverY += beamVec[i].yAver;
        bunchAverZ += beamVec[i].zAver;
	    bunchAverPX+= beamVec[i].pxAver;
	    bunchAverPY+= beamVec[i].pyAver;
	    bunchAverPZ+= beamVec[i].pzAver;
    }

    strongStrongBunchInfo->bunchAverX  = bunchAverX  / totBunchNum ;
    strongStrongBunchInfo->bunchAverY  = bunchAverY  / totBunchNum ;
    strongStrongBunchInfo->bunchAverZ  = bunchAverZ  / totBunchNum ;
    strongStrongBunchInfo->bunchAverPX = bunchAverPX / totBunchNum;
    strongStrongBunchInfo->bunchAverPY = bunchAverPY / totBunchNum;
    strongStrongBunchInfo->bunchAverPZ = bunchAverPZ / totBunchNum;
}

void MPBeam::MPBeamDataPrintPerTurn(int nTurns, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter)
{

    gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (beamVec.size());
	workspace = gsl_fft_complex_workspace_alloc (beamVec.size());


    double tempX[2*beamVec.size()];
    double tempY[2*beamVec.size()];
    double tempZ[2*beamVec.size()];
    
    for(int i=0;i<beamVec.size();i++)
	{
        tempX[2*i]  = beamVec[i].xAver;
        tempY[2*i]  = beamVec[i].yAver;
        tempZ[2*i]  = beamVec[i].zAver;
        tempX[2*i+1]= beamVec[i].pxAver;
        tempY[2*i+1]= beamVec[i].pyAver;
        tempZ[2*i+1]= beamVec[i].pzAver;
        //tempX[2*i+1]= 0;
        //tempY[2*i+1]= 0;
        //tempZ[2*i+1]= 0;
	}

    gsl_fft_complex_forward (tempX,1, beamVec.size(), wavetable, workspace);
    gsl_fft_complex_forward (tempY,1, beamVec.size(), wavetable, workspace);
    gsl_fft_complex_forward (tempZ,1, beamVec.size(), wavetable, workspace);
    gsl_fft_complex_workspace_free (workspace);
    gsl_fft_complex_wavetable_free (wavetable);



    string filePrefix = inputParameter.ringRun->TBTBunchAverData;
    string fname = filePrefix + ".sdds";
    ofstream fout(fname,ios_base::app);
   
    if(!inputParameter.ringRun->TBTBunchAverData.empty())
    {
        string filePrefix = inputParameter.ringRun->TBTBunchAverData;
        string fname = filePrefix + ".sdds";
        ofstream fout(fname,ios_base::app);
        
        if(nTurns==0)
	    {
	        fout<<"SDDS1"<<endl;
	        fout<<"&column name=Turns,                           type=long,  &end"<<endl;
	        fout<<"&column name=HarmIndex,                       type=long,  &end"<<endl;
            fout<<"&column name=TotIonCharge,       units=e,     type=float,  &end"<<endl;
	        fout<<"&column name=AverX,              units=m,     type=float,  &end"<<endl;
	        fout<<"&column name=AverY,              units=m,     type=float,  &end"<<endl;
	        fout<<"&column name=AverZ,              units=m,     type=float,  &end"<<endl;
	        fout<<"&column name=AverXP,             units=rad,   type=float,  &end"<<endl;
	        fout<<"&column name=AverYP,             units=rad,   type=float,  &end"<<endl;
	        fout<<"&column name=AverZP,             units=rad,   type=float,  &end"<<endl;
	        fout<<"&column name=AverZAna,           units=m,     type=float,  &end"<<endl;
	        fout<<"&column name=rmsBunchLengthAna,  units=m,     type=float,  &end"<<endl; 
	        fout<<"&column name=RmsEmitX,           units=m*rad, type=float,  &end"<<endl;
	        fout<<"&column name=RmsEmitY,           units=m*rad, type=float,  &end"<<endl;
	        fout<<"&column name=RmsEmitZ,           units=m*rad, type=float,  &end"<<endl;
	        fout<<"&column name=RmsRX,              units=m*rad, type=float,  &end"<<endl;
	        fout<<"&column name=RmsRY,              units=m*rad, type=float,  &end"<<endl;
	        fout<<"&column name=RmsBunchLength,     units=m,     type=float,  &end"<<endl;
     	    fout<<"&column name=RmsBunchEnergySpread, units=rad, type=float,  &end"<<endl;
            fout<<"&column name=ModeNum,                               type=long,   &end"<<endl;
            fout<<"&column name=AbsFFTAverX,                           type=float,  &end"<<endl;
            fout<<"&column name=AbsFFTAverY,                           type=float,  &end"<<endl;
            fout<<"&column name=AbsFFTAverZ,                           type=float,  &end"<<endl; 

     	    for(int j=0; j<inputParameter.ringParRf->resNum;j++)
            {
    	        string colname = string("&column name=") + string("cavAmp_") + to_string(j) + string(", units=V,   type=float,  &end");
                fout<<colname<<endl;
                colname = string("&column name=") + string("cavPhase_") + to_string(j) + string(", units=rad,   type=float,  &end");
                fout<<colname<<endl;
                colname = string("&column name=") + string("cavReal_") + to_string(j) + string(", units=rad,   type=float,  &end");
                fout<<colname<<endl;

                colname = string("&column name=") + string("genAmp_") + to_string(j) + string(", units=V,   type=float,  &end");
                fout<<colname<<endl;
                colname = string("&column name=") + string("genPhase_") + to_string(j) + string(", units=rad,   type=float,  &end");
                fout<<colname<<endl;
                colname = string("&column name=") + string("beamIndAmp_") + to_string(j) + string(", units=rad,   type=float,  &end");
                fout<<colname<<endl;
                colname = string("&column name=") + string("beamIndPhase_") + to_string(j) + string(", units=rad,   type=float,  &end");
                fout<<colname<<endl;
 	        }
	        fout<<"&data mode=ascii, &end"<<endl;
	    }

        fout<<"! page number "<<nTurns + 1<<endl;
        fout<<beamVec.size()<<endl;

        for(int i=0;i<beamVec.size();i++)
        {
            fout<<setw(24)<<left<<setprecision(16)<<nTurns
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].bunchHarmNum
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].totIonCharge
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].xAver
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].yAver
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].zAver
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].pxAver
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].pyAver
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].pzAver
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].zAverAnalytical
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].bunchLengthAnalytical
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].emittanceX
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].emittanceY
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].emittanceZ
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].rmsRx
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].rmsRy
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].rmsBunchLength
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].rmsEnergySpread
                <<setw(24)<<left<<setprecision(16)<<i
                <<setw(24)<<left<<setprecision(16)<<sqrt( pow(tempX[2*i],2) + pow(tempX[2*i+1],2) )
                <<setw(24)<<left<<setprecision(16)<<sqrt( pow(tempY[2*i],2) + pow(tempY[2*i+1],2) )
                <<setw(24)<<left<<setprecision(16)<<sqrt( pow(tempZ[2*i],2) + pow(tempZ[2*i+1],2) );

            for(int j=0; j<inputParameter.ringParRf->resNum;j++)
            {
                fout<<setw(24)<<left<<setprecision(16)<<    beamVec[i].cavFBCenInfo->cavAmpBunchCen[j]
                    <<setw(24)<<left<<setprecision(16)<<    beamVec[i].cavFBCenInfo->cavPhaseBunchCen[j]
                    <<setw(24)<<left<<setprecision(16)<<    beamVec[i].cavFBCenInfo->cavAmpBunchCen[j] * cos(beamVec[i].cavFBCenInfo->cavPhaseBunchCen[j])
                    <<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].cavFBCenInfo->genVolBunchAver[j])
                    <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].cavFBCenInfo->genVolBunchAver[j])
                    <<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].cavFBCenInfo->induceVolBunchCen[j])
                    <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].cavFBCenInfo->induceVolBunchCen[j]);;
            }
            fout<<endl;

        }
        fout.close();
    }
    

    // print out the speciflied bunch distribution....
    int TBTBunchPrintNum                  = inputParameter.ringRun->TBTBunchPrintNum;
    vector<int> TBTBunchDisDataBunchIndex = inputParameter.ringRun->TBTBunchDisDataBunchIndex;

    if(!inputParameter.ringRun->TBTBunchDisData.empty()  &&  (TBTBunchPrintNum !=0) )
    {
        for(int i=0;i<TBTBunchPrintNum;i++)
        {
            string filePrefix = inputParameter.ringRun->TBTBunchDisData;
            string fname = filePrefix +"_" +to_string(TBTBunchDisDataBunchIndex[i])+".sdds";
            ofstream fout1(fname,ios_base::app);

            if(nTurns==0)
	        {
	            fout1<<"SDDS1"<<endl;
	            fout1<<"&parameter name=partNum,               type=float,  &end"<<endl;
	            fout1<<"&parameter name=AverX,    units=m,     type=float,  &end"<<endl;
	            fout1<<"&parameter name=AverY,    units=m,     type=float,  &end"<<endl;
	            fout1<<"&parameter name=AverZ,    units=m,     type=float,  &end"<<endl;
                fout1<<"&parameter name=AverZAnalytical,    units=m,     type=float,  &end"<<endl;

	            fout1<<"&parameter name=rmsEmitX, units=m*rad, type=float,  &end"<<endl;
	            fout1<<"&parameter name=rmsEmitY, units=m*rad, type=float,  &end"<<endl;
	            fout1<<"&parameter name=rmsEmitZ, units=m*rad, type=float,  &end"<<endl;
	            fout1<<"&parameter name=rmsRx,    units=m,     type=float,  &end"<<endl;
	            fout1<<"&parameter name=rmsRy,    units=m,     type=float,  &end"<<endl;
	            fout1<<"&parameter name=rmsBunchLen,    units=m,     type=float,  &end"<<endl;
                fout1<<"&parameter name=rmsBunchLenAnalytical,    units=m,     type=float,  &end"<<endl;
	            fout1<<"&parameter name=rmsBunchEnergySpread,    units=rad,     type=float,  &end"<<endl;


	            fout1<<"&column name=x,              units=m,     type=float,  &end"<<endl;
	            fout1<<"&column name=y,              units=m,     type=float,  &end"<<endl;
	            fout1<<"&column name=z,              units=m,     type=float,  &end"<<endl;
	            fout1<<"&column name=xp,             units=rad,   type=float,  &end"<<endl;
	            fout1<<"&column name=yp,             units=rad,   type=float,  &end"<<endl;
	            fout1<<"&column name=zp,             units=rad,   type=float,  &end"<<endl;
	            fout1<<"&data mode=ascii, &end"<<endl;
	        }

            int bunchIndex  = TBTBunchDisDataBunchIndex[i];

            fout1<<"! page number "<<nTurns * TBTBunchPrintNum + i +1 <<endl;

            fout1<<beamVec[bunchIndex].electronNumPerBunch<<endl;
            fout1<<beamVec[bunchIndex].xAver<<endl;
            fout1<<beamVec[bunchIndex].yAver<<endl;
            fout1<<beamVec[bunchIndex].zAver<<endl;
            fout1<<beamVec[bunchIndex].zAverAnalytical<<endl;
            fout1<<beamVec[bunchIndex].emittanceX<<endl;
            fout1<<beamVec[bunchIndex].emittanceY<<endl;
            fout1<<beamVec[bunchIndex].emittanceZ<<endl;
            fout1<<beamVec[bunchIndex].rmsRx<<endl;
            fout1<<beamVec[bunchIndex].rmsRy<<endl;
            fout1<<beamVec[bunchIndex].rmsBunchLength<<endl;
            fout1<<beamVec[bunchIndex].bunchLengthAnalytical<<endl;
            fout1<<beamVec[bunchIndex].rmsEnergySpread<<endl;

            fout1<<beamVec[bunchIndex].macroEleNumPerBunch<<endl;

            for(int j=0; j<beamVec[i].macroEleNumPerBunch;j++)
            {
                fout1<<setw(20)<<left<<beamVec[bunchIndex].ePositionX[j]
                     <<setw(20)<<left<<beamVec[bunchIndex].ePositionY[j]
                     <<setw(20)<<left<<beamVec[bunchIndex].ePositionZ[j]
                     <<setw(20)<<left<<beamVec[bunchIndex].eMomentumX[j]
                     <<setw(20)<<left<<beamVec[bunchIndex].eMomentumY[j]
                     <<setw(20)<<left<<beamVec[bunchIndex].eMomentumZ[j]
                     <<endl;
            }

            fout1.close();
        }

    }
    
    // print out the rf voltage and phase along the specified bunch
    if(!inputParameter.ringRun->TBTBunchCavVolData.empty()  &&  (TBTBunchPrintNum !=0) )
    {
        for(int i=0;i<TBTBunchPrintNum;i++)
        {
            string filePrefix = inputParameter.ringRun->TBTBunchCavVolData;
            string fname = filePrefix +"_" +to_string(i)+".sdds";
            ofstream fout2(fname,ios_base::app);

            if(nTurns==0)
	        {
	            fout2<<"SDDS1"<<endl;
	            fout2<<"&column name=z,              units=m,     type=float,  &end"<<endl;
	            fout2<<"&column name=rho,            units=C/m,   type=float,  &end"<<endl;
	            fout2<<"&column name=potenWell,                   type=float,  &end"<<endl;
	            fout2<<"&column name=rhoAnalytical,               type=float,  &end"<<endl;

    	        for(int j=0; j<inputParameter.ringParRf->resNum;j++)
	            {

                    string colname = string("&column name=") + string("cavReal_") + to_string(j) + string(", units=V,   type=float,  &end");
                    fout2<<colname<<endl;
                    colname = string("&column name=") + string("cavImag_") + to_string(j) + string(", units=V,   type=float,  &end");
                    fout2<<colname<<endl;

                    colname = string("&column name=") + string("cavAbs_") + to_string(j) + string(", units=V,   type=float,  &end");
                    fout2<<colname<<endl;
                    colname = string("&column name=") + string("cavArg_") + to_string(j) + string(", units=rad,   type=float,  &end");
                    fout2<<colname<<endl;
                    colname = string("&column name=") + string("cavDpz_") + to_string(j) + string(", units=rad,   type=float,  &end");
                    fout2<<colname<<endl;
	            }
	            string colname = string("&column name=") + string("cavDpzTot")  + string(", units=rad,   type=float,  &end");
	            fout2<<colname<<endl;
	            fout2<<"&data mode=ascii, &end"<<endl;
            }

            int bunchIndex  = TBTBunchDisDataBunchIndex[i];

            fout2<<"! page number "<<nTurns * TBTBunchPrintNum + i +1 <<endl;
            fout2<<inputParameter.ringImpedance->bunchBinNumberZ<<endl;
            
            for(int k=0;k<inputParameter.ringImpedance->bunchBinNumberZ;k++)
            {
                fout2<<setw(20)<<left<<beamVec[bunchIndex].posZBins[k]
                     <<setw(20)<<left<<beamVec[bunchIndex].densProfVsBin[k]
                     <<setw(20)<<left<<beamVec[bunchIndex].hamiltonPotenWell[k]
                     <<setw(20)<<left<<beamVec[bunchIndex].densProfVsBinAnalytical[k];

                double dpzTemp=0;
                for(int j=0; j<inputParameter.ringParRf->resNum;j++)
                {
                    complex<double> temp = beamVec[bunchIndex].cavVolInfoVsLongBins[j][k];
                    double temp1 =  beamVec[bunchIndex].cavForceInfoVsLongBins[j][k];
                    fout2<<setw(20)<<left<<temp.real()
                         <<setw(20)<<left<<temp.imag()
                         <<setw(20)<<left<<abs(temp)
                         <<setw(20)<<left<<arg(temp)
                         <<setw(20)<<left<<temp1;

                    dpzTemp +=temp1;
                }
                fout2<<setw(20)<<left<<dpzTemp<<endl;

            }


            fout2.close();
        }


    }
    


}

void MPBeam::SSBeamIonEffectOneInteractionPoint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k)
{
    int totBunchNum = beamVec.size();

    for(int j=0;j<totBunchNum;j++)
    {
        latticeInterActionPoint.GetIonNumberPerInterAction(beamVec[j].electronNumPerBunch, k);
        latticeInterActionPoint.IonGenerator(beamVec[j].rmsRx,beamVec[j].rmsRy,beamVec[j].xAver,beamVec[j].yAver,k);        
        latticeInterActionPoint.IonsUpdate(k);
        latticeInterActionPoint.IonRMSCal(k);
        beamVec[j].SSIonBunchInteraction(latticeInterActionPoint,k);
        beamVec[j].BunchTransferDueToIon(latticeInterActionPoint,k);
        latticeInterActionPoint.IonTransferDueToBunch(beamVec[j].bunchGap,k,strongStrongBunchInfo->bunchSizeXMax,strongStrongBunchInfo->bunchSizeYMax);
    }

}

void MPBeam::BeamTransferPerInteractionPointDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToLatticeT(latticeInterActionPoint,k);
    }
}

void MPBeam::BeamTransferPerTurnDueToLatticeL(ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    GetBinDistBetweenBunch(inputParameter);
    
    vector<double> resAmpFBRatioForTotSelfLoss = inputParameter.ringParRf->resAmpFBRatioForTotSelfLoss;

    vector<complex<double> > vbKickAver;
    vector<complex<double> > vbAccumAver;
    vector<complex<double> > selfLossToCompensate;
    complex<double> totSelfLoss (0.e0,0.e0);
    vbKickAver.resize(inputParameter.ringParRf->resNum);
    vbAccumAver.resize(inputParameter.ringParRf->resNum);
    selfLossToCompensate.resize(inputParameter.ringParRf->resNum);

    
    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        complex<double> vbKickAverSum=(0.e0,0.e0);
        complex<double> vbAccumAverSum=(0.e0,0.e0);

        for(int j=0;j<beamVec.size();j++)
        {
            vbKickAverSum      += beamVec[j].cavFBCenInfo->selfLossVolBunchCen[i];
            vbAccumAverSum     += beamVec[j].cavFBCenInfo->induceVolBunchCen[i];
        }
        vbKickAver[i]      = vbKickAverSum   / double(beamVec.size());
        vbAccumAver[i]     = vbAccumAverSum  / double(beamVec.size());
    }

    
    
    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        totSelfLoss += vbKickAver[i];
    }


    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        selfLossToCompensate[i] = resAmpFBRatioForTotSelfLoss[i] * totSelfLoss;
    }
    
    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        //each cavity generator compensated the related self-loss by it self
        //double genAddvbArg = arg(cavityResonator.resonatorVec[i].resCavVolReq);
        //double genAddvbAbs = (cavityResonator.resonatorVec[i].resCavVolReq.real() - selfLossToCompensate[i].real()) / cos( genAddvbArg );
        //cavityResonator.resonatorVec[i].resGenVol = genAddvbAbs * exp( li * genAddvbArg ) - vbAccumAver[i];

        cavityResonator.resonatorVec[i].resGenVol = cavityResonator.resonatorVec[i].resCavVolReq - vbAccumAver[i];

        for(int j=0;j<beamVec.size();j++)
        {
            beamVec[j].cavFBCenInfo->genVolBunchAver[i] =  cavityResonator.resonatorVec[i].resGenVol;
        }

    }
        
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToLatticeL(inputParameter,cavityResonator);
    }

}

void MPBeam::GetBinDistBetweenBunch(ReadInputSettings &inputParameter)
{
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double circRing   = inputParameter.ringParBasic->circRing ;
    double rfLen      = circRing  / ringHarmH;

    for(int i=0;i<beamVec.size();i++)
    {

        if(i<beamVec.size()-1)
        {
            timeBetweenBunch[i] = beamVec[i].bunchGap * rfLen - (beamVec[i].zAver - beamVec[i].rmsBunchLength/2.0) + (beamVec[i+1].zAver + beamVec[i+1].rmsBunchLength/2.0);
        }
        else
        {
            timeBetweenBunch[i] = beamVec[i].bunchGap * rfLen + (beamVec[i].zAver - beamVec[i].rmsBunchLength/2.0) + (beamVec[i+1].zAver + beamVec[i+1].rmsBunchLength/2.0);
        }

        timeBetweenBunch[i] /= CLight;
        beamVec[i].timeToLastBunch = timeBetweenBunch[i];
    }

}

void MPBeam::SSIonDataPrint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,int nTurns)
{
    // print out the ion distribution at different interaction points
    // probably just print out the total ions speices and numbers. print ions different location may not a good idea.  
    int numberOfInteraction =  latticeInterActionPoint.numberOfInteraction;

    string filePrefix = inputParameter.ringIonEffPara->ionDisWriteTo;
    string fname = filePrefix + ".sdds";
    ofstream fout(fname,ios_base::app);

    for(int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)
    {
        
        if(nTurns==0)
        {
            fout<<"SDDS1"<<endl;
            fout<<"&parameter name=totIonCharge, units=e, type=float,  &end"<<endl;            
            for(int p=0;p<latticeInterActionPoint.ionAccumuPositionX[k].size();p++)
            {
                string col= string("parameter name=") + string("numOfMacroIon_") + to_string(p) + string(",  type=long,  &end");   
                fout<<col<<endl;
            }
        
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

        fout<<"! page number "<<nTurns * numberOfInteraction + k + 1 <<endl;
        fout<<latticeInterActionPoint.totIonCharge<<endl;
        for(int p=0;p<latticeInterActionPoint.ionAccumuPositionX[k].size();p++)
        {
            fout<<latticeInterActionPoint.ionAccumuPositionX[k][p].size()<<endl;            
        }

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

} 

void MPBeam::BeamTransferPerTurnDueToLattice(LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    BeamTransferPerTurnDueToLatticeT(latticeInterActionPoint);
    BeamTransferPerTurnDueToLatticeL(inputParameter,cavityResonator); 
}

void MPBeam::BeamTransferPerTurnDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint)
{
    for(int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)
    {
        BeamTransferPerInteractionPointDueToLatticeT(latticeInterActionPoint,k);
    }
}



void MPBeam::FIRBunchByBunchFeedback(FIRFeedBack &firFeedBack,int nTurns)
{
    //y[0] = \sum_0^{N} a_k x[-k]. Ref. Nakamura's paper spring 8 notation here used is the same with Nakamura's paper

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
        if(energyKickU[i]     > firFeedBack.fIRBunchByBunchFeedbackKickLimit)
        {
            energyKickU[i]    = firFeedBack.fIRBunchByBunchFeedbackKickLimit;
        }
    }


    for(int i=0;i<beamVec.size();i++)
    {
        for(int j=0;j<beamVec[i].macroEleNumPerBunch;j++)
        {
            beamVec[i].eMomentumX[j] = beamVec[i].eMomentumX[j] + tranAngleKickx[i];
            beamVec[i].eMomentumY[j] = beamVec[i].eMomentumY[j] + tranAngleKicky[i];
            beamVec[i].eMomentumZ[j] = beamVec[i].eMomentumZ[j] + energyKickU[i];
        }

    }

}

void MPBeam::BeamSynRadDamping(const ReadInputSettings &inputParameter, const LatticeInterActionPoint &latticeInterActionPoint)
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchSynRadDamping(inputParameter,latticeInterActionPoint);
    }
}

void MPBeam::LRWakeBeamIntaction(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction,const  LatticeInterActionPoint &latticeInterActionPoint)
{
    
    int nTurnswakeTrunction     = inputParameter.ringLRWake->nTurnswakeTrunction;
    int harmonics               = inputParameter.ringParBasic->harmonics;
    double electronBeamEnergy   = inputParameter.ringParBasic->electronBeamEnergy;
    double rBeta                = inputParameter.ringParBasic->rBeta;

    // prepare bunch center position of the previous turns
    vector<double> posxDataTemp(beamVec.size());
    vector<double> posyDataTemp(beamVec.size());
    vector<double> poszDataTemp(beamVec.size());
    // current turn position
    for(int i=0;i<beamVec.size();i++)
    {
        posxDataTemp[i] = beamVec[i].xAver;
        posyDataTemp[i] = beamVec[i].yAver;
        poszDataTemp[i] = beamVec[i].zAver;
    }

    // [0,nTurnswakeTrunction] stores [pos[-nTurnswakeTrunction],pos[0]] 
    wakefunction.posxData.erase(wakefunction.posxData.begin());
    wakefunction.posyData.erase(wakefunction.posyData.begin());
    wakefunction.poszData.erase(wakefunction.poszData.begin());

    wakefunction.posxData.push_back(posxDataTemp);
    wakefunction.posyData.push_back(posyDataTemp);
    wakefunction.poszData.push_back(poszDataTemp);

    vector<double> wakeForceTemp(3,0);

    double tauij=0.e0;
    int nTauij=0;
    double tRF  = inputParameter.ringParBasic->t0 / double(harmonics);
    double deltaZij=0;

    int tempIndex0,tempIndex1;

    for (int j=0;j<beamVec.size();j++)
    {
	    beamVec[j].lRWakeForceAver[0] =0.E0;      // x
	    beamVec[j].lRWakeForceAver[1] =0.E0;      // y
	    beamVec[j].lRWakeForceAver[2] =0.E0;      // z
        
	    for(int n=0;n<nTurnswakeTrunction;n++)
	    {          
            // self-interation of bunch in current turn is excluded if tempIndex1 = j - 1, when n=0.             
            if(n==0)
            {
                tempIndex0   = 0;
                tempIndex1   = j ;
            }
            else if (n==nTurnswakeTrunction-1)
            {
                tempIndex0   = j;
                tempIndex1   = beamVec.size()-1;
            }
            else
            {
                tempIndex0   = 0;
                tempIndex1   = beamVec.size()-1;
            }
                        
            for(int i=tempIndex0;i<=tempIndex1;i++)
            {
                deltaZij = wakefunction.poszData[nTurnswakeTrunction-1-n][i] - beamVec[j].zAver ;
                nTauij   = beamVec[i].bunchHarmNum - beamVec[j].bunchHarmNum - n * harmonics;
                tauij    = nTauij * tRF + deltaZij / CLight;

	            if(!inputParameter.ringLRWake->pipeGeoInput.empty())
	            {                                       
                    wakeForceTemp = wakefunction.GetRWLRWakeFun(tauij); 
                    beamVec[j].lRWakeForceAver[0] -= wakeForceTemp[0] * beamVec[i].electronNumPerBunch * wakefunction.posxData[nTurnswakeTrunction-1-n][i];    //[V/C/m] * [m] ->  [V/C]
                    beamVec[j].lRWakeForceAver[1] -= wakeForceTemp[1] * beamVec[i].electronNumPerBunch * wakefunction.posyData[nTurnswakeTrunction-1-n][i];    //[V/C/m] * [m] ->  [V/C]
                    beamVec[j].lRWakeForceAver[2] -= wakeForceTemp[2] * beamVec[i].electronNumPerBunch ;    //    [V/C]
                }

                if(!inputParameter.ringLRWake->bbrInput.empty())
	            {
                    wakeForceTemp = wakefunction.GetBBRWakeFun(tauij);                 
                    beamVec[j].lRWakeForceAver[0] -= wakeForceTemp[0] * beamVec[i].electronNumPerBunch * wakefunction.posxData[nTurnswakeTrunction-1-n][i];    //[V/C/m] * [m] ->  [V/C]
                    beamVec[j].lRWakeForceAver[1] -= wakeForceTemp[1] * beamVec[i].electronNumPerBunch * wakefunction.posyData[nTurnswakeTrunction-1-n][i];    //[V/C/m] * [m] ->  [V/C]
                    beamVec[j].lRWakeForceAver[2] -= wakeForceTemp[2] * beamVec[i].electronNumPerBunch ;                                                       //[V/C]         ->  [V/C]
                }
            }
                   
        }
        
        
        beamVec[j].lRWakeForceAver[0] *=  ElectronCharge / electronBeamEnergy / pow(rBeta,2);   // [V/C] * [C] * [1e] / [eV] ->rad
        beamVec[j].lRWakeForceAver[1] *=  ElectronCharge / electronBeamEnergy / pow(rBeta,2);   // [V/C] * [C] * [1e] / [eV] ->rad
        beamVec[j].lRWakeForceAver[2] *=  ElectronCharge / electronBeamEnergy / pow(rBeta,2);   // [V/C] * [C] * [1e] / [eV] ->rad

    }

    // Reference to "simulation of transverse multi-bunch instabilities of proton beam in LHC, PHD thesis, Alexander Koshik P. 32, Eq. (3.22)"
    // For this subroutine, all bunch feels the same kick strength from long range wake kciks.
} 
void MPBeam::BeamTransferPerTurnDueWake()
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToWake();
    }
}