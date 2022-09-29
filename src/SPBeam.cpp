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
#include "SPBeam.h"
//#include "LongImpSingalBunch.h"
#include "Faddeeva.h"
#include "WakeFunction.h"
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
#include <gsl/gsl_multifit_nlin.h>
#include <complex>
#include <iomanip>
#include <cstring>




using namespace std;
using std::vector;
using std::complex;


SPBeam::SPBeam()
{
}

SPBeam::~SPBeam()
{
    delete weakStrongBeamInfo;
}


void SPBeam::Initial(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter)
{
    
    int totBunchNum = inputParameter.ringFillPatt->totBunchNumber; 
    int harmonics = inputParameter.ringParBasic->harmonics;              
    beamVec.resize(totBunchNum);

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
    // in one bunch train, leading and tail bunch covention [head, tail] = [0,1,2,...i,]
    for(int i=0;i<totBunchNum;i++)
    {
        beamVec[i].InitialSPBunch(inputParameter);
        beamVec[i].DistriGenerator(latticeInterActionPoint,inputParameter,i);

    }

    for(int i=0;i<totBunchNum-1;i++)
    {
        beamVec[i].bunchGap = beamVec[i+1].bunchHarmNum - beamVec[i].bunchHarmNum;
    }
    beamVec[totBunchNum-1].bunchGap = harmonics - beamVec[totBunchNum-1].bunchHarmNum;

    // get the initial bunch distance
    SPGetBeamInfo();
    GetTimeDisToNextBunchIntial(inputParameter);

    RMOutPutFiles();
    /*
    for(int i=0;i<totBunchNum;i++)
    {
        cout<<beamVec[i].bunchGap<<"    "<<i<<" "<<beamVec[i].bunchHarmNum<<endl;
    }
    getchar();
    */    

//----------------------------------------------------------------------------------------
//  set section is used to generate the beam filling pattern data from mb track.

    ofstream fout ("mb_track_filling.sdds",ios::out);
    counter=0;
    for(int i=0;i<harmonics;i++)
    {
        if(counter<beamVec.size())
        {
            if(i<beamVec[counter].bunchHarmNum)
            {
                fout<<0<<endl;
            }
            else if(i==beamVec[counter].bunchHarmNum)
            {
                fout<<1<<endl;
                counter += 1;
            }
        }
        else
        {
            fout<<0<<endl;
        }
    }
    fout.close();
    //------------------------------------end ------------------------------------------

    
    // set section is used to generate the beam filling pattern data for elegant .
    ofstream fout1 ("elegant_filling_train_para.sdds",ios::out);

    fout1<<"SDDS1"<<endl;
    fout1<<"&parameter name=BucketNumber, type=long, &end"<<endl;
    fout1<<"&parameter name=Intensity, type=long, &end"<<endl;
    fout1<<"&data mode=ascii, &end"<<endl;

    counter=0;

    for(int i=0;i<harmonics;i++)
    {
        if(counter<beamVec.size())
        {
            if(i==beamVec[counter].bunchHarmNum)
            {
                fout1<<"!page number "<<counter+1<<endl;
                fout1<<i<<endl;
                fout1<<1<<endl;
                counter += 1;
            }
        }

    }
    fout1.close();
    //------------------------------------end ------------------------------------------
     
}

void SPBeam::InitialcavityResonator(ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    // According the beam filling pattern, tracking from transient to get steady state Vb. Ref. Bill Chapter 7.4.2
    // The notatoion of Bill's book is the same in P.B.Wilson's Section 3. Slac Pub 6062  except li -> - li, the phase rotation oppsiste.

    int resNum      = inputParameter.ringParRf->resNum;
    int ringHarmH   = inputParameter.ringParBasic->harmonics;
    double beamCurr = inputParameter.ringBunchPara->current * beamVec.size();
    double f0       = inputParameter.ringParBasic->f0;
    double t0       = inputParameter.ringParBasic->t0;

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

        vbKickAccum  =   complex<double>(0,0);              // get the accumme of  vb0/2
        vbKickAver   =   complex<double>(0,0);

        for(int n=0;n<nTurns;n++)
        {
            
            for(int k=0;k<beamVec.size(); k++)
            {
                vb0  = complex<double>(-1 * cavityResonator.resonatorVec[i].resFre * 2 * PI * cavityResonator.resonatorVec[i].resShuntImpRs
                                              /  cavityResonator.resonatorVec[i].resQualityQ0, 0.E0)
                                              *  beamVec[k].electronNumPerBunch * ElectronCharge;       
                                              // wilson's equation p.6, The same as Bill equation (7.76), cos convention, real part reresents the energy gain.
                if( n==nTurns-1)
                {
                    vbKickAccum += vb0/2.0 ;
                }

                vbAccum += vb0;                
                tB = beamVec[k].bunchGap * t0 / ringHarmH ;           //[s]
                deltaL = tB / tF ;                
                
                // cPsi = 2 * PI * (cavityResonator.resonatorVec[i].resFre - ringHarmH * f0 * cavityResonator.resonatorVec[i].resHarm) * tB;
                // cPsi   = deltaL * tan(cavityResonator.resonatorVec[i].resDeTunePsi ); 
                cPsi = 2 * PI * cavityResonator.resonatorVec[i].resFre * tB;
                
                if( n==nTurns-1)
                {
                    vbAccumAver1 +=  vbAccum;
                }

                vbAccum = vbAccum * exp(- deltaL ) * exp (li * cPsi);    // decay and rotate [V] P.B. Eq.(3.12) get the beam induced voltage at the cavity

                if( n==nTurns-1)
                {
                    vbAccumAver2 +=  vbAccum;
                }
            }
        }
        
        vbAccumAver =  (vbAccumAver1 + vbAccumAver2) /2.0/ double(beamVec.size());   // (1)
        vbAccumAver =  vbAccumAver1/double(beamVec.size());                          // before decay and rotate, just  after bunch left.
        // vbAccumAver =  vbAccumAver2/double(beamVec.size());                       // after decay and rotate.
 
        cavityResonator.resonatorVec[i].vbAccum   =  vbAccumAver;               // set cavity beam induced voltage
        cavityResonator.resonatorVec[i].vbAccum0  =  vbAccumAver;
        cavityResonator.resonatorVec[i].vbAccumRFFrame = vbAccumAver;
        
        vbKickAver  =  vbKickAccum /double(beamVec.size());                     // average energy that bunch is kicked by self-loss Vb0/2. 


        // extra feedback condition to compensate the self-loss term vb0/2          
        // double genAddvbArg = arg(cavityResonator.resonatorVec[i].resCavVolReq);
        // if(genAddvbArg==PI/2.0 or genAddvbArg==-PI/2.0)
        // {
        //    cerr<<"Required Cavity Phase is PI/2 or -PI/2"<<endl;
        //    cerr<<"Does not work when cavity feedback is included"<<endl;
        //    exit(0);
        // }
        // double genAddvbAbs = (cavityResonator.resonatorVec[i].resCavVolReq.real() - vbKickAver.real()) / cos( genAddvbArg );
        // genAddvbAbs = abs(genAddvbAbs);
        // cavityResonator.resonatorVec[i].resGenVol = genAddvbAbs * exp(li * genAddvbArg ) - vbAccumAver;
        //end of the extra feedback 


        // DC solution -- which is exactly the same as above equaiton (1)'s reuslts, very good agreement.
        //vbAccumAver = 2.0 * beamCurr * cavityResonator.resonatorVec[i].resShuntImpRs / (1.0 + cavityResonator.resonatorVec[i].resCouplingBeta)
        //                             * cos(cavityResonator.resonatorVec[i].resDeTunePsi) * exp(li * (PI + cavityResonator.resonatorVec[i].resDeTunePsi)  );

        // simply condition. resGenVol can not compensate the self-loss term, bunch will have a center shift due to self-loss.  
        cavityResonator.resonatorVec[i].resGenVol = cavityResonator.resonatorVec[i].resCavVolReq - vbAccumAver;



         // set the cold and warm cavity condition. Vb is set up before particle tracking in the "warm" condition
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
            // beamVec[j].cavFBCenInfo->cavAmpBunchCen[i]    = abs(cavityResonator.resonatorVec[i].resCavVolReq);
            // beamVec[j].cavFBCenInfo->cavPhaseBunchCen[i]  = arg(cavityResonator.resonatorVec[i].resCavVolReq);
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
    fout<<"&parameter name=CavAmpIdeal,   units=V,   type=float,  &end"<<endl;
	fout<<"&parameter name=CavPhaseIdeal,    units=rad, type=float,  &end"<<endl;
	fout<<"&parameter name=CavFreq,          units=Hz,  type=float,  &end"<<endl;
	fout<<"&parameter name=GenAmp,           units=V,   type=float,  &end"<<endl;
	fout<<"&parameter name=GenPhase,         units=rad, type=float,  &end"<<endl;
    fout<<"&parameter name=beamIndAmp,       units=V,   type=float,  &end"<<endl;
	fout<<"&parameter name=beamIndPhase,     units=rad, type=float,  &end"<<endl;
    fout<<"&parameter name=CavDetunPsi,      units=rad, type=float,  &end"<<endl;

	fout<<"&column name=Turns,                            type=float,  &end"<<endl;
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
        tF      = cavityResonator.resonatorVec[i].tF  ;      // [s]
        nTurns  = ceil(100 * tF * f0);
        vbAccum = complex<double>(0,0);
         
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

                tB = beamVec[j].bunchGap * t0 / ringHarmH ;  //[s]
                
                time +=  tB;
                deltaL = tB / tF ;

                // cPsi = 2 * PI * (cavityResonator.resonatorVec[i].resFre - ringHarmH * f0 * cavityResonator.resonatorVec[i].resHarm) * tB; 
                cPsi = 2 * PI * cavityResonator.resonatorVec[i].resFre * tB;

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

                vbAccum = vbAccum * exp(- deltaL ) * exp ( li * cPsi);    // decay and rotate...  [V]-- use the voltage after decay and feed this into tracking.
            }
        }
    }

    fout.close(); 

    // beam induced voltate, generator voltage, cavity voltage are all printout in RF frame the initial valye      
}

void SPBeam::Run(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter, CavityResonator &cavityResonator)
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
    fout<<"&column name=Turns,          type=long,              &end"<<endl;
    fout<<"&column name=IonCharge,      units=e,   type=float,  &end"<<endl;
    fout<<"&column name=MaxAverX,       units=m,   type=float,  &end"<<endl;
    fout<<"&column name=MaxAverY,       units=m,   type=float,  &end"<<endl;
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
    fout<<"&column name=MaxJx,          units=m,   type=float,  &end"<<endl;
    fout<<"&column name=MaxJy,          units=m,   type=float,  &end"<<endl;

    for(int j=0;j<inputParameter.ringRun->TBTBunchPrintNum;j++)
    {
        int bunchPrinted=inputParameter.ringRun->TBTBunchDisDataBunchIndex[j];
        string colname = string("&column name=") + string("x_")     + to_string(bunchPrinted) + string(", units=m,   type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("px_")            + to_string(bunchPrinted) + string(", units=rad, type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("y_")             + to_string(bunchPrinted) + string(", units=m,   type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("py_")            + to_string(bunchPrinted) + string(", units=rad, type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("z_")             + to_string(bunchPrinted) + string(", units=m,   type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("pz_")            + to_string(bunchPrinted) + string(", units=rad, type=float,  &end");
        fout<<colname<<endl;
    }
    fout<<"&data mode=ascii, &end"<<endl;
    fout<<"! page number "<<1<<endl;
    fout<<nTurns<<endl;

    
    
    // run loop starts, for nTrns and each trun for k interaction-points
    for(int n=0;n<nTurns;n++)
    {
        if(n%100==0) 
        {
            cout<<n<<"  turns"<<endl;
        }
        
        if(beamIonFlag) 
        {
            for (int k=0;k<inputParameter.ringIonEffPara->numberofIonBeamInterPoint;k++)
            {
                SPBeamRMSCal(latticeInterActionPoint, k);
                WSBeamIonEffectOneInteractionPoint(inputParameter,latticeInterActionPoint, n, k);
                BeamTransferPerInteractionPointDueToLatticeT(latticeInterActionPoint,k);             //transverse transfor per interaction point
            }

            if(ionInfoPrintInterval && (n%ionInfoPrintInterval==0))
            {                
                SPBeamRMSCal(latticeInterActionPoint, 0);
                WSIonDataPrint(inputParameter,latticeInterActionPoint, n);
            } 
        }
        else
        {
            SPBeamRMSCal(latticeInterActionPoint, 0);
            BeamTransferPerTurnDueToLatticeT(latticeInterActionPoint);            
        }
        
        SPBeamRMSCal(latticeInterActionPoint, 0);
        BeamTransferPerTurnDueToLatticeL(inputParameter,latticeInterActionPoint,cavityResonator,n); 
   
        if(lRWakeFlag)
        {
            SPBeamRMSCal(latticeInterActionPoint, 0);
            LRWakeBeamIntaction(inputParameter,lRWakeFunction,latticeInterActionPoint,n);
            BeamTransferPerTurnDueWake();   
        }
        
        if(fIRBunchByBunchFeedbackFlag)  
        {     
            SPBeamRMSCal(latticeInterActionPoint, 0);
            FIRBunchByBunchFeedback(firFeedBack,n);
        }

        if(synRadDampingFlag[0]==1) // transverse SR effect
        {
            SPBeamRMSCal(latticeInterActionPoint, 0);
            BeamSynRadDamping(inputParameter,latticeInterActionPoint);
        }
        
        SPBeamRMSCal(latticeInterActionPoint, 0);
        SPGetBeamInfo();
       
        if(bunchInfoPrintInterval && (n%bunchInfoPrintInterval==0)  )
        {      
            SPBeamDataPrintPerTurn(n,latticeInterActionPoint,inputParameter);           
        }

        fout<<n<<"  "
            <<setw(15)<<left<< latticeInterActionPoint.totIonCharge
            <<setw(15)<<left<< weakStrongBeamInfo->bunchAverXMax     //
            <<setw(15)<<left<< weakStrongBeamInfo->bunchAverYMax     //
            <<setw(15)<<left<< weakStrongBeamInfo->bunchAverX        // over all bunches with the bunch centor <x> value.
            <<setw(15)<<left<< weakStrongBeamInfo->bunchAverY        // over all bunches
            <<setw(15)<<left<< weakStrongBeamInfo->bunchAverZ        // over all bunches
            <<setw(15)<<left<< weakStrongBeamInfo->bunchAverPX       // over all bunches
            <<setw(15)<<left<< weakStrongBeamInfo->bunchAverPY       // over all bunches
            <<setw(15)<<left<< weakStrongBeamInfo->bunchAverPZ       // over all bunches
            <<setw(15)<<left<< weakStrongBeamInfo->bunchRmsSizeX     // over all bunches -- with bunch center <x> of each bunch to get this rms value.
            <<setw(15)<<left<< weakStrongBeamInfo->bunchRmsSizeY     // over all bunches
            <<setw(15)<<left<< weakStrongBeamInfo->bunchRmsSizeZ     // over all bunches
            <<setw(15)<<left<< weakStrongBeamInfo->bunchRmsSizePX    // over all bunches
            <<setw(15)<<left<< weakStrongBeamInfo->bunchRmsSizePY    // over all bunches
            <<setw(15)<<left<< weakStrongBeamInfo->bunchRmsSizePZ    // over all bunches
            <<setw(15)<<left<< log10(sqrt(weakStrongBeamInfo->actionJxMax))    
            <<setw(15)<<left<< log10(sqrt(weakStrongBeamInfo->actionJyMax));
        
        for(int i=0;i<inputParameter.ringRun->TBTBunchPrintNum;i++)
        {
            int index = inputParameter.ringRun->TBTBunchDisDataBunchIndex[i];
        
            fout<<setw(15)<<left<<beamVec[index].xAver
                <<setw(15)<<left<<beamVec[index].pxAver
                <<setw(15)<<left<<beamVec[index].yAver
                <<setw(15)<<left<<beamVec[index].pyAver
                <<setw(15)<<left<<beamVec[index].zAver
                <<setw(15)<<left<<beamVec[index].pzAver;
        
        }
        fout<<endl;                                           
    
        // getchar();
    }
    fout.close();
    cout<<"End of Tracking "<<nTurns<< " Turns"<<endl;
       
    // after single particle tracking- with the RF data to get the Haissinski solution -- not a self-consistent process, since
    // the bunch centor shift due to the potenwell distortation also affect the cavity voltage and phase.
    // It must be treated by numerical simulaiton...  
     
    if( !inputParameter.ringRun->TBTBunchHaissinski.empty() &&  !cavityResonator.resonatorVec.empty())
    {
        GetHaissinski(inputParameter,cavityResonator,sRWakeFunction); 
        if(!inputParameter.ringRun->TBTBunchLongTraj.empty())
        {
            GetAnalyticalLongitudinalPhaseSpace(inputParameter,cavityResonator,sRWakeFunction);
        }
    }    
    //sddsplot();    
}


void SPBeam::SetBeamPosHistoryDataWithinWindow()
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].SetBunchPosHistoryDataWithinWindow();
    }
}

void SPBeam::SPBeamRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].GetSPBunchRMS(latticeInterActionPoint,k);
    } 
}

void SPBeam::SPGetBeamInfo()
{
    double totBunchNum = beamVec.size();
    vector<double> tempJx;
    vector<double> tempJy;
    vector<double> tempAverX;
    vector<double> tempAverY;
    
    for(int i=0;i<totBunchNum;i++)
    {
        tempJx.push_back(beamVec[i].actionJx);
        tempJy.push_back(beamVec[i].actionJy);        
        tempAverX.push_back(abs(beamVec[i].xAver));
        tempAverY.push_back(abs(beamVec[i].yAver));
    }
    
    double maxJx    = *max_element(tempJx.begin(),tempJx.end() );
    double maxJy    = *max_element(tempJy.begin(),tempJy.end() );
    
    double maxAverX = *max_element(tempAverX.begin(),tempAverX.end() );
    double maxAverY = *max_element(tempAverY.begin(),tempAverY.end() );
    
    

    weakStrongBeamInfo->actionJxMax    = maxJx;
    weakStrongBeamInfo->actionJyMax    = maxJy;    
    weakStrongBeamInfo->bunchAverXMax  = maxAverX;
    weakStrongBeamInfo->bunchAverYMax  = maxAverY;


    weakStrongBeamInfo->bunchEffectiveSizeXMax = beamVec[0].rmsRx + maxAverX;
    weakStrongBeamInfo->bunchEffectiveSizeYMax = beamVec[0].rmsRy + maxAverY;


    double bunchRmsSizeXTemp=0.E0;
    double bunchRmsSizeYTemp=0.E0;
    double bunchRmsSizeZTemp=0.E0;
    double bunchRmsSizePXTemp=0.E0;
    double bunchRmsSizePYTemp=0.E0;
    double bunchRmsSizePZTemp=0.E0;


    for(int i=0;i<totBunchNum;i++)
    {
        bunchRmsSizeXTemp += pow(beamVec[i].xAver,2);
        bunchRmsSizeYTemp += pow(beamVec[i].yAver,2);
        bunchRmsSizeZTemp += pow(beamVec[i].zAver,2);
	    bunchRmsSizePXTemp+= pow(beamVec[i].pxAver,2);
	    bunchRmsSizePYTemp+= pow(beamVec[i].pyAver,2);
	    bunchRmsSizePZTemp+= pow(beamVec[i].pzAver,2);
    }

    weakStrongBeamInfo->bunchRmsSizeX  = sqrt(bunchRmsSizeXTemp  /totBunchNum);
    weakStrongBeamInfo->bunchRmsSizeY  = sqrt(bunchRmsSizeYTemp  /totBunchNum);
    weakStrongBeamInfo->bunchRmsSizeZ  = sqrt(bunchRmsSizeZTemp  /totBunchNum);
    weakStrongBeamInfo->bunchRmsSizePX = sqrt(bunchRmsSizePXTemp /totBunchNum);
    weakStrongBeamInfo->bunchRmsSizePY = sqrt(bunchRmsSizePYTemp /totBunchNum);
    weakStrongBeamInfo->bunchRmsSizePZ = sqrt(bunchRmsSizePZTemp /totBunchNum);

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


    weakStrongBeamInfo->bunchAverX  = bunchAverX  / totBunchNum;
    weakStrongBeamInfo->bunchAverY  = bunchAverY  / totBunchNum;
    weakStrongBeamInfo->bunchAverZ  = bunchAverZ  / totBunchNum;
    weakStrongBeamInfo->bunchAverPX = bunchAverPX / totBunchNum;
    weakStrongBeamInfo->bunchAverPY = bunchAverPY / totBunchNum;
    weakStrongBeamInfo->bunchAverPZ = bunchAverPZ / totBunchNum;
}

void SPBeam::GetHilbertAnalyticalInOneTurn()
{    
    vector<double> xSignal(beamVec.size());
    vector<double> ySignal(beamVec.size());
    vector<double> zSignal(beamVec.size());

    for(int i=0;i<beamVec.size();i++)
	{
        xSignal[i] = beamVec[i].xAver;
        ySignal[i] = beamVec[i].yAver;
        zSignal[i] = beamVec[i].zAver;

    }
    vector<complex<double> > xAnalytical = GetHilbertAnalytical(xSignal);
    vector<complex<double> > yAnalytical = GetHilbertAnalytical(ySignal);
    vector<complex<double> > zAnalytical = GetHilbertAnalytical(zSignal);    

    for(int i=0;i<beamVec.size();i++)
	{
        beamVec[i].xAverAnalytical = xAnalytical[i];
        beamVec[i].yAverAnalytical = yAnalytical[i];
        beamVec[i].zAverAnalytical = zAnalytical[i];
    }

    // Alex Chao--Eq.(2.94)--verfication
    // with a RLC model, given only the real part of the impedacne and get the analytical signal from hilbert transform
    // it can build the image part of the impedacne.  

    // zSignal.resize(1000);
    // ofstream fout("test.dat");
    // for(int i=0;i<zSignal.size();i++)
    // {
    //     double denominator = 10.0 * (50. / i  - i /50.);
    //     complex<double> zSignalComplex = 1.0 / (1.0 + li * denominator);                
    //     zSignal[i]  = zSignalComplex.real();
    //     // test of a damping signla
    //     // zSignal[i]  = exp(0.001 * i ) * cos( 23.12 * i / zSignal.size() * 2 * PI  );  
    // }

    // vector<complex<double> > zSignalAnalytical =  GetHilbertAnalytical(zSignal);

    // for(int i=0;i<zSignal.size();i++)
    // {
    //     fout<<setw(15)<<left<<i
    //         <<setw(15)<<left<<zSignalAnalytical[i].real()
    //         <<setw(15)<<left<<zSignalAnalytical[i].imag()
    //         <<setw(15)<<left<<abs(zSignalAnalytical[i])
    //         <<setw(15)<<left<<arg(zSignalAnalytical[i])
    //         <<endl;
    // }
    // fout.close();
    // cout<<__LINE__<<endl;
    // getchar();
} 

vector<complex<double> > SPBeam::GetHilbertAnalytical(vector<double> signal)
{   
    // with hilbert transform to get analytical signal Bessy in VSR
    // Ref. Calculation of coupled bunch effects in the synchrotorn loight source bessy VSR
    gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (signal.size());
	workspace = gsl_fft_complex_workspace_alloc (signal.size());
    
    double xFFT[2*signal.size()];

    for(int i=0;i<signal.size();i++)
	{
        xFFT[2* i    ] = signal[i];
        xFFT[2* i + 1] = 0;
    }
    gsl_fft_complex_forward (xFFT,1, signal.size(), wavetable, workspace);

    // keep the DC term here, some reference set this term to zero...
    // for(int i=0;i<1;i++)
	// {        
    //     xFFT[2*i    ] = 0;
    //     xFFT[2*i + 1] = 0;
    // }

    for(int i=1;i<=signal.size()/2;i++)
	{        
        xFFT[2*i    ] = 2 * xFFT[2*i   ];
        xFFT[2*i + 1] = 2 * xFFT[2*i +1];
    }

    if(signal.size() % 2 ==0 )
    {
        int index = signal.size() / 2;
        xFFT[2 * index    ] =   xFFT[2 * index   ];
        xFFT[2 * index + 1] =   xFFT[2 * index + 1];
    }

    for (int i=signal.size()/2 +1;i<signal.size();i++)
    {
        xFFT[2*i    ] = 0.E0;
        xFFT[2*i + 1] = 0.E0;
    }

    gsl_fft_complex_inverse (xFFT,1, signal.size(), wavetable, workspace);

    vector<complex<double> > analyticalSingal(signal.size());

    for(int i=0;i<signal.size();i++)
    {
        analyticalSingal[i] = complex<double> (xFFT[2*i], xFFT[2*i +1]);
    }

    gsl_fft_complex_workspace_free (workspace);
    gsl_fft_complex_wavetable_free (wavetable);


    return analyticalSingal;
}


void SPBeam::SPBeamDataPrintPerTurn(int nTurns, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter)
{
    // get the spectrum of the center motion...
    
    GetHilbertAnalyticalInOneTurn(); 
 
    gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (beamVec.size());
	workspace = gsl_fft_complex_workspace_alloc (beamVec.size());

    double alphax = latticeInterActionPoint.twissAlphaX[0];
    double betax  = latticeInterActionPoint.twissBetaX [0];
    double alphay = latticeInterActionPoint.twissAlphaY[0];
    double betay  = latticeInterActionPoint.twissBetaY[0] ;
    double alphaz = 0;
    double betaz  = inputParameter.ringBunchPara->rmsBunchLength / inputParameter.ringBunchPara->rmsEnergySpread;
    double workQx = inputParameter.ringParBasic->workQx;
    double workQy = inputParameter.ringParBasic->workQy;
    double workQz = inputParameter.ringParBasic->workQz;

    double tempX[2*beamVec.size()];
    double tempY[2*beamVec.size()];
    double tempZ[2*beamVec.size()];
    double tempXFFT[2*beamVec.size()];
    double tempYFFT[2*beamVec.size()];
    double tempZFFT[2*beamVec.size()];

    // IPAC 2022 WEPOMS010 -- Diamond-II -- WangsiWei's IPAC paper. 
    for(int i=0;i<beamVec.size();i++)
	{
        double x  =  beamVec[i].xAver;
        double y  =  beamVec[i].yAver;
        double z  =  beamVec[i].zAver; 
        
        double px =  beamVec[i].pxAver;
        double py =  beamVec[i].pyAver;
        double pz =  beamVec[i].pzAver;

        // transverse is clockwise rotation from ith bunch to  bpm position
        // re-distribution bunch position to 
        double phasex = - 2.0 * PI * workQx * beamVec[i].bunchHarmNum / inputParameter.ringParBasic->harmonics;
        double phasey = - 2.0 * PI * workQy * beamVec[i].bunchHarmNum / inputParameter.ringParBasic->harmonics;
        double phasez =   2.0 * PI * workQz * beamVec[i].bunchHarmNum / inputParameter.ringParBasic->harmonics;

        complex<double> zx = ( x / sqrt(betax) - li * ( sqrt(betax) * px + alphax / sqrt(betax) * x ) ) * exp (li * phasex);
        complex<double> zy = ( y / sqrt(betay) - li * ( sqrt(betay) * py + alphay / sqrt(betay) * y ) ) * exp (li * phasey); 
        complex<double> zz = ( z / sqrt(betaz) - li * ( sqrt(betaz) * pz + alphaz / sqrt(betaz) * z ) ) * exp (li * phasez);  

        tempX[2*i]  = zx.real();
        tempY[2*i]  = zy.real();
        tempZ[2*i]  = zz.real();
        tempX[2*i+1]= zx.imag();
        tempY[2*i+1]= zy.imag();
        tempZ[2*i+1]= zz.imag();

        tempXFFT[2*i    ] = tempX[2*i    ];
        tempXFFT[2*i + 1] = tempX[2*i + 1];
        tempYFFT[2*i    ] = tempY[2*i    ];
        tempYFFT[2*i + 1] = tempY[2*i + 1];
        tempZFFT[2*i    ] = tempZ[2*i    ];
        tempZFFT[2*i + 1] = tempZ[2*i + 1];        
    }
        
    gsl_fft_complex_forward (tempXFFT,1, beamVec.size(), wavetable, workspace);
    gsl_fft_complex_forward (tempYFFT,1, beamVec.size(), wavetable, workspace);
    gsl_fft_complex_forward (tempZFFT,1, beamVec.size(), wavetable, workspace);
  
     
    string filePrefix = inputParameter.ringRun->TBTBunchAverData;
    string fname = filePrefix + ".sdds";
    ofstream fout(fname,ios_base::app);
   
    if(nTurns==0)
	{
	    fout<<"SDDS1"<<endl;
	    fout<<"&column name=Turns,                                 type=long,   &end"<<endl;
	    fout<<"&column name=HarmIndex,                             type=long,   &end"<<endl;
	    fout<<"&column name=AverX,                      units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=AverY,                      units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=AverZ,                      units=m,   type=float,  &end"<<endl;
	    fout<<"&column name=AverXP,                     units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=AverYP,                     units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=AverZP,                     units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=Jx,                         units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=Jy,                         units=rad, type=float,  &end"<<endl;
        fout<<"&column name=TotIonCharge,               units=e,   type=float,  &end"<<endl;
	    fout<<"&column name=dpxIon,                     units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=dpyIon,                     units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=dpxLRW,                     units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=dpyLRW,                     units=rad, type=float,  &end"<<endl;
	    fout<<"&column name=dpzLRW,                     units=rad, type=float,  &end"<<endl;
        fout<<"&column name=ModeNum,                               type=long,   &end"<<endl;
	    fout<<"&column name=AbsFFTAverX,                           type=float,  &end"<<endl;
	    fout<<"&column name=AbsFFTAverY,                           type=float,  &end"<<endl;
	    fout<<"&column name=AbsFFTAverZ,                           type=float,  &end"<<endl;
        fout<<"&column name=ArgFFTAverX,                           type=float,  &end"<<endl;
	    fout<<"&column name=ArgFFTAverY,                           type=float,  &end"<<endl;
	    fout<<"&column name=ArgFFTAverZ,                           type=float,  &end"<<endl;
        fout<<"&column name=AbsHilbertAverX,                       type=float,  &end"<<endl;
        fout<<"&column name=AbsHilbertAverY,                       type=float,  &end"<<endl;
        fout<<"&column name=AbsHilbertAverZ,                       type=float,  &end"<<endl;
        fout<<"&column name=ArgHilbertAverX,                       type=float,  &end"<<endl;
        fout<<"&column name=ArgHilbertAverY,                       type=float,  &end"<<endl;
        fout<<"&column name=ArgHilbertAverZ,                       type=float,  &end"<<endl;
        fout<<"&column name=HilbertAverX,                          type=float,  &end"<<endl;
        fout<<"&column name=HilbertAverY,                          type=float,  &end"<<endl;
        fout<<"&column name=HilbertAverZ,                          type=float,  &end"<<endl;
        fout<<"&column name=HilbertAverXP,                         type=float,  &end"<<endl;
        fout<<"&column name=HilbertAverYP,                         type=float,  &end"<<endl;
        fout<<"&column name=HilbertAverZP,                         type=float,  &end"<<endl;
        fout<<"&column name=timeToNextBunch,             units=s,  type=float,  &end"<<endl;

	    for(int j=0; j<inputParameter.ringParRf->resNum;j++)
	    {
	        string colname = string("&column name=") + string("cavAmp_") + to_string(j) + string(", units=V,   type=float,  &end");
            fout<<colname<<endl;
            colname = string("&column name=") + string("cavPhase_")      + to_string(j) + string(", units=rad, type=float,  &end");
            fout<<colname<<endl;
            colname = string("&column name=") + string("genAmp_")        + to_string(j) + string(", units=V,   type=float,  &end");
            fout<<colname<<endl;
            colname = string("&column name=") + string("genPhase_")      + to_string(j) + string(", units=rad, type=float,  &end");
            fout<<colname<<endl;
            colname = string("&column name=") + string("beamIndAmp_")    + to_string(j) + string(", units=rad, type=float,  &end");
            fout<<colname<<endl;
            colname = string("&column name=") + string("beamIndPhase_")  + to_string(j) + string(", units=rad, type=float,  &end");
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
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].xAver
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].yAver
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].zAver
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].pxAver
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].pyAver
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].pzAver
            // <<setw(24)<<left<<setprecision(16)<<tempX[2*i]
            // <<setw(24)<<left<<setprecision(16)<<tempY[2*i]
            // <<setw(24)<<left<<setprecision(16)<<tempZ[2*i] 
            // <<setw(24)<<left<<setprecision(16)<<tempX[2*i+1]
            // <<setw(24)<<left<<setprecision(16)<<tempY[2*i+1]
            // <<setw(24)<<left<<setprecision(16)<<tempZ[2*i+1]
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].actionJx
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].actionJy
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].totIonCharge
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].eFxDueToIon[0]
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].eFyDueToIon[0]
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].lRWakeForceAver[0]
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].lRWakeForceAver[1]
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].lRWakeForceAver[2]
            <<setw(24)<<left<<setprecision(16)<<i
            <<setw(24)<<left<<setprecision(16)<<sqrt( pow(tempXFFT[2*i],2) + pow(tempXFFT[2*i+1],2) )
            <<setw(24)<<left<<setprecision(16)<<sqrt( pow(tempYFFT[2*i],2) + pow(tempYFFT[2*i+1],2) )
            <<setw(24)<<left<<setprecision(16)<<sqrt( pow(tempZFFT[2*i],2) + pow(tempZFFT[2*i+1],2) )
            <<setw(24)<<left<<setprecision(16)<<atan2(    tempXFFT[2*i],         tempXFFT[2*i+1] )
            <<setw(24)<<left<<setprecision(16)<<atan2(    tempYFFT[2*i],         tempYFFT[2*i+1] )
            <<setw(24)<<left<<setprecision(16)<<atan2(    tempZFFT[2*i],         tempZFFT[2*i+1] )
            <<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].xAverAnalytical)
            <<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].yAverAnalytical)
            <<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].zAverAnalytical)
            <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].xAverAnalytical)
            <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].yAverAnalytical)
            <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].zAverAnalytical)
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].xAverAnalytical.real()
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].yAverAnalytical.real()
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].zAverAnalytical.real()
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].xAverAnalytical.imag()
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].yAverAnalytical.imag()
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].zAverAnalytical.imag()
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].timeFromCurrnetBunchToNextBunch; 


            for(int j=0; j<inputParameter.ringParRf->resNum;j++)
            {
                fout<<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].cavFBCenInfo->cavVolBunchCen[j])
                    <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].cavFBCenInfo->cavVolBunchCen[j])
                    <<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].cavFBCenInfo->genVolBunchAver[j])
                    <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].cavFBCenInfo->genVolBunchAver[j])
                    <<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].cavFBCenInfo->induceVolBunchCen[j])
                    <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].cavFBCenInfo->induceVolBunchCen[j]);
            }
            fout<<endl;
    }

    fout.close();
    
    
    // ofstream foutX("coupledBunchedModeX.sdds",ios_base::app);
    // ofstream foutY("coupledBunchedModeY.sdds",ios_base::app);

    // foutX<<setw(15)<<left<<nTurns;
    // foutY<<setw(15)<<left<<nTurns;
    // for(int i=0;i<beamVec.size();i++)
    // {    
    //     foutX<<setw(15)<<left<< tempXFFTAmp[i] ;
    //     foutY<<setw(15)<<left<< tempYFFTAmp[i] ;
    // }
    // foutX<<endl;
    // foutY<<endl;
    // foutX.close();
    // foutY.close();

    // ofstream foutZ("coupledBunchedModeZ.sdds",ios_base::app);
    // foutZ<<setw(15)<<left<<nTurns;
    // for(int i=0;i<beamVec.size();i++)
    // {    
    //     foutZ<<setw(15)<<left<< tempZFFTAmp[i] ;
    // }
    // foutZ<<endl;
    // foutZ.close();


    // get the growth rate of the coupled bunched mode
    
    int turns = nTurns + inputParameter.ringRun->bunchInfoPrintInterval;
   
    vector<double> tempXFFTAmp(beamVec.size());
    vector<double> tempYFFTAmp(beamVec.size());
    vector<double> tempZFFTAmp(beamVec.size());
    
    vector<double> tempXHilbertAmp(beamVec.size());
    vector<double> tempYHilbertAmp(beamVec.size());
    vector<double> tempZHilbertAmp(beamVec.size());

    vector<double> historyTempAverX(beamVec.size());
    vector<double> historyTempAverY(beamVec.size());
    vector<double> historyTempAverZ(beamVec.size());
   
    for (int i=0;i<beamVec.size();i++)
    {
        tempXFFTAmp[i] = log(sqrt( pow(tempXFFT[2*i],2) + pow(tempXFFT[2*i+1],2) ) );
        tempYFFTAmp[i] = log(sqrt( pow(tempYFFT[2*i],2) + pow(tempYFFT[2*i+1],2) ) );
        tempZFFTAmp[i] = log(sqrt( pow(tempZFFT[2*i],2) + pow(tempZFFT[2*i+1],2) ) );

        tempXHilbertAmp[i]  = log( abs(beamVec[i].xAverAnalytical) );
        tempYHilbertAmp[i]  = log( abs(beamVec[i].yAverAnalytical) );
        tempZHilbertAmp[i]  = log( abs(beamVec[i].zAverAnalytical) );

        historyTempAverX[i] =  beamVec[i].xAver;
        historyTempAverY[i] =  beamVec[i].yAver;
        historyTempAverZ[i] =  beamVec[i].zAver;
    }

    coupledBunchModeAmpX.push_back(tempXFFTAmp);
    coupledBunchModeAmpY.push_back(tempYFFTAmp);
    coupledBunchModeAmpZ.push_back(tempZFFTAmp);

    hilbertCoupledBunchModeAmpX.push_back(tempXHilbertAmp);
    hilbertCoupledBunchModeAmpY.push_back(tempYHilbertAmp);
    hilbertCoupledBunchModeAmpZ.push_back(tempZHilbertAmp);

    historyAverX.push_back(historyTempAverX);
    historyAverY.push_back(historyTempAverY);
    historyAverZ.push_back(historyTempAverZ);

    
    if( turns == inputParameter.ringRun->nTurns )
    {
        
        int dim =  coupledBunchModeAmpX.size();
        ofstream foutCoupMode("coupledBunchedModeGrowthRate.sdds");
        foutCoupMode<<"SDDS1"<<endl;
	    foutCoupMode<<"&column name=ModeIndex,                                type=long,   &end"<<endl;
	    foutCoupMode<<"&column name=CBMGRX,                      units=1/s,   type=float,   &end"<<endl;
	    foutCoupMode<<"&column name=CBMGRY,                      units=1/s,   type=float,  &end"<<endl;
        foutCoupMode<<"&column name=CBMGRZ,                      units=1/s,   type=float,  &end"<<endl;
        foutCoupMode<<"&column name=CBMHilbertGRX,               units=1/s,   type=float,   &end"<<endl;
	    foutCoupMode<<"&column name=CBMHilbertGRY,               units=1/s,   type=float,  &end"<<endl;
        foutCoupMode<<"&column name=CBMHilbertGRZ,               units=1/s,   type=float,  &end"<<endl;
        foutCoupMode<<"&column name=BunchHilbertGRX,             units=1/s,   type=float,   &end"<<endl;
	    foutCoupMode<<"&column name=BunchHilbertGRY,             units=1/s,   type=float,  &end"<<endl;
        foutCoupMode<<"&column name=BunchHilbertGRZ,             units=1/s,   type=float,  &end"<<endl;
        foutCoupMode<<"&data mode=ascii, &end"<<endl;
        foutCoupMode<<"! page number "<<1<<endl;
        foutCoupMode<<beamVec.size()<<endl;

        FittingGSL fittingGSL;        
        double fitWeight[dim];
        double fitX[dim];
        double tempXFFTAmpOneMode[dim];
        double tempYFFTAmpOneMode[dim];
        double tempZFFTAmpOneMode[dim]; 

        vector<double> xAverHis(dim);
        vector<double> yAverHis(dim);
        vector<double> zAverHis(dim);


        for(int i=0; i<beamVec.size(); i++)
        {
            // coupled bunch mode gorwth rate calculation--- mode-by-mode process -- FFT of averZ  in one turn.
            vector<double> resFitX(2,0.E0);
            vector<double> resFitY(2,0.E0);
            vector<double> resFitZ(2,0.E0);
                       
            for(int n=0;n<dim;n++)
            {                
                tempXFFTAmpOneMode[n]  = coupledBunchModeAmpX[n][i];
                tempYFFTAmpOneMode[n]  = coupledBunchModeAmpY[n][i];
                tempZFFTAmpOneMode[n]  = coupledBunchModeAmpZ[n][i];
                fitWeight[n]          = 1.E0;
                fitX[n]               = n * inputParameter.ringRun->bunchInfoPrintInterval;              
            } 

            resFitX = fittingGSL.FitALinear(fitX, fitWeight, tempXFFTAmpOneMode,dim);
            resFitY = fittingGSL.FitALinear(fitX, fitWeight, tempYFFTAmpOneMode,dim);
            resFitZ = fittingGSL.FitALinear(fitX, fitWeight, tempZFFTAmpOneMode,dim);
            ////////////////////////////////////////////////////////////////////////////////////////////

            // coupled bunch mode gorwth rate calculation--- mode-by-mode process -- FFT analytical signal in one turn.
            vector<double> resHilbertCBMFitX(2,0.E0);
            vector<double> resHilbertCBMFitY(2,0.E0);
            vector<double> resHilbertCBMFitZ(2,0.E0);

            for(int n=0;n<dim;n++)
            {    
                tempXFFTAmpOneMode[n]  = hilbertCoupledBunchModeAmpX[n][i];
                tempYFFTAmpOneMode[n]  = hilbertCoupledBunchModeAmpY[n][i];
                tempZFFTAmpOneMode[n]  = hilbertCoupledBunchModeAmpZ[n][i];
                fitWeight[n]          = 1.E0;
                fitX[n]               = n * inputParameter.ringRun->bunchInfoPrintInterval;    
            }
            
            resHilbertCBMFitX = fittingGSL.FitALinear(fitX, fitWeight, tempXFFTAmpOneMode,dim);
            resHilbertCBMFitY = fittingGSL.FitALinear(fitX, fitWeight, tempYFFTAmpOneMode,dim);
            resHilbertCBMFitZ = fittingGSL.FitALinear(fitX, fitWeight, tempZFFTAmpOneMode,dim);
            ////////////////////////////////////////////////////////////////////////////////////////////

            // bunch-by-bunch growth rate calculation --- bunch-by-bunch process 
            // With Hilber transform  get the envelope and do fitting.
            vector<double> resFitHilbertBunchX(2,0.E0);
            vector<double> resFitHilbertBunchY(2,0.E0);
            vector<double> resFitHilbertBunchZ(2,0.E0);

            for(int n=0;n<dim;n++)
            {                
                xAverHis[n]  = historyAverX[n][i];
                yAverHis[n]  = historyAverY[n][i];
                zAverHis[n]  = historyAverZ[n][i];
            } 

            vector<complex<double> > xAnalytical = GetHilbertAnalytical(xAverHis);
            vector<complex<double> > yAnalytical = GetHilbertAnalytical(yAverHis);
            vector<complex<double> > zAnalytical = GetHilbertAnalytical(zAverHis);    

            double tempXHilberOneBunch[dim], tempYHilberOneBunch[dim], tempZHilberOneBunch[dim];
            
            for(int n=0;n<dim;n++)
            {                
                tempXHilberOneBunch[n] = log(abs(xAnalytical[n]));
                tempYHilberOneBunch[n] = log(abs(yAnalytical[n]));
                tempZHilberOneBunch[n] = log(abs(zAnalytical[n]));   
                fitWeight[n]          = 1.E0;
                fitX[n]               = n * inputParameter.ringRun->bunchInfoPrintInterval; 
             
            } 

            resFitHilbertBunchX = fittingGSL.FitALinear(fitX, fitWeight, tempXHilberOneBunch,dim);
            resFitHilbertBunchY = fittingGSL.FitALinear(fitX, fitWeight, tempYHilberOneBunch,dim);
            resFitHilbertBunchZ = fittingGSL.FitALinear(fitX, fitWeight, tempZHilberOneBunch,dim);
        
            // print out results
            foutCoupMode<<setw(15)<<left<<i
                        <<setw(15)<<left<<resFitX[1]             / inputParameter.ringParBasic->t0
                        <<setw(15)<<left<<resFitY[1]             / inputParameter.ringParBasic->t0
                        <<setw(15)<<left<<resFitZ[1]             / inputParameter.ringParBasic->t0
                        <<setw(15)<<left<<resHilbertCBMFitX[1]   / inputParameter.ringParBasic->t0
                        <<setw(15)<<left<<resHilbertCBMFitY[1]   / inputParameter.ringParBasic->t0
                        <<setw(15)<<left<<resHilbertCBMFitX[1]   / inputParameter.ringParBasic->t0
                        <<setw(15)<<left<<resFitHilbertBunchX[1] / inputParameter.ringParBasic->t0
                        <<setw(15)<<left<<resFitHilbertBunchY[1] / inputParameter.ringParBasic->t0
                        <<setw(15)<<left<<resFitHilbertBunchZ[1] / inputParameter.ringParBasic->t0 
                        <<endl;              
        }

        foutCoupMode.close();
    }


}

void SPBeam::WSBeamIonEffectOneInteractionPoint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k)
{

    int totBunchNum = beamVec.size();


    // ofstream fout("ioncharge.sdds",ios_base::app);
    for(int j=0;j<totBunchNum;j++)
    {
        latticeInterActionPoint.GetIonNumberPerInterAction(beamVec[j].electronNumPerBunch, k);
        latticeInterActionPoint.IonGenerator(beamVec[j].rmsRx,beamVec[j].rmsRy,beamVec[j].xAver,beamVec[j].yAver,k);
        latticeInterActionPoint.IonsUpdate(k);
        beamVec[j].WSIonBunchInteraction(latticeInterActionPoint,k);
        beamVec[j].BunchTransferDueToIon(latticeInterActionPoint,k);
        latticeInterActionPoint.IonTransferDueToBunch(beamVec[j].bunchGap,k,weakStrongBeamInfo->bunchEffectiveSizeXMax,weakStrongBeamInfo->bunchEffectiveSizeYMax);
        
        // fout<<setw(15)<<left<<j + nTurns * totBunchNum     
        //     <<setw(15)<<left<<latticeInterActionPoint.totIonCharge
        //     <<setw(15)<<left<<weakStrongBeamInfo->bunchEffectiveSizeXMax
        //     <<setw(15)<<left<<weakStrongBeamInfo->bunchEffectiveSizeYMax
        //     <<endl;

    }
    // cout<<__LINE__<<endl;
    // getchar();

}

void SPBeam::BeamTransferPerTurnDueToLatticeL(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,CavityResonator &cavityResonator, int turns)
{ 
       
    vector<double> resAmpFBRatioForTotSelfLoss = inputParameter.ringParRf->resAmpFBRatioForTotSelfLoss;

    vector<complex<double> > vbKickAver;
    vector<complex<double> > vbAccumAver;
    vector<complex<double> > selfLossToCompensate;
    complex<double> totSelfLoss;
    vbKickAver.resize(inputParameter.ringParRf->resNum);
    vbAccumAver.resize(inputParameter.ringParRf->resNum);
    selfLossToCompensate.resize(inputParameter.ringParRf->resNum);
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double circRing   = inputParameter.ringParBasic->circRing;
    double rfLen      = circRing  / ringHarmH;
    double tRF        = inputParameter.ringParBasic->t0 / ringHarmH; 
    
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
        
    // for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    // {
    //     totSelfLoss += vbKickAver[i];
    // }


    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        // selfLossToCompensate[i] = resAmpFBRatioForTotSelfLoss[i] * totSelfLoss;
        selfLossToCompensate[i] = vbKickAver[i];
    }
    

    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        //each cavity generator compensated the related self-loss by it self
        //double genAddvbArg = arg(cavityResonator.resonatorVec[i].resCavVolReq);
        //double genAddvbAbs = (cavityResonator.resonatorVec[i].resCavVolReq.real() - selfLossToCompensate[i].real()) / cos( genAddvbArg );
        //cavityResonator.resonatorVec[i].resGenVol = genAddvbAbs * exp( li * genAddvbArg ) - vbAccumAver[i];
        
        // active cavity, and cavity feedback is set as 1 in input file, HeFei's method.
        if(cavityResonator.resonatorVec[i].rfResCavVolFB==1  && cavityResonator.resonatorVec[i].resType==1)
        {    
            cavityResonator.resonatorVec[i].resGenVol = cavityResonator.resonatorVec[i].resCavVolReq - vbAccumAver[i];
        }
        else if(cavityResonator.resonatorVec[i].rfResCavVolFB==0  && cavityResonator.resonatorVec[i].resType==1)
        {
            cavityResonator.resonatorVec[i].resGenVol = cavityResonator.resonatorVec[i].resGenVol;
        }
        else if(cavityResonator.resonatorVec[i].rfResCavVolFB==1  && cavityResonator.resonatorVec[i].resType==0)
        {
            cerr<<"Wrong setting: passive cavity can set cavitiy FB"<<endl;
            exit(0);
        }
        else if (cavityResonator.resonatorVec[i].rfResCavVolFB==0  && cavityResonator.resonatorVec[i].resType==0)
        {
            cavityResonator.resonatorVec[i].resGenVol = complex<double>(0.E0,0.E0);
        }
        
        for(int j=0;j<beamVec.size();j++)
        {
            beamVec[j].cavFBCenInfo->genVolBunchAver[i] =  cavityResonator.resonatorVec[i].resGenVol;
        }
        
        // set voltage and phase for hainssinki solver
        for(int j=0;j<beamVec.size();j++)
        {              
            beamVec[j].haissinski->cavAmp[i]   = abs(beamVec[j].cavFBCenInfo->induceVolBunchCen[i] + cavityResonator.resonatorVec[i].resGenVol);
            beamVec[j].haissinski->cavPhase[i] = arg(beamVec[j].cavFBCenInfo->induceVolBunchCen[i] + cavityResonator.resonatorVec[i].resGenVol);            
        }

    }


    //longitudinal tracking ---
    if(beamVec.size()==1)
    {
        beamVec[0].timeFromCurrnetBunchToNextBunch  = beamVec[0].bunchGap * tRF - (beamVec[0].zAver - beamVec[0].zAverLastTurn) / CLight;
        beamVec[0].zAverLastTurn                    = beamVec[0].zAver;        
        // beamVec[0].BunchTransferDueToLatticeLMatarix(inputParameter,cavityResonator);
        beamVec[0].BunchTransferDueToLatticeLTest(inputParameter,cavityResonator);
        beamVec[0].GetSPBunchRMS(latticeInterActionPoint, 0);   
    }
    else
    {
        // save the last turn position info    
        for(int i=0;i<beamVec.size();i++)
        {
            beamVec[i].zAverLastTurn = beamVec[i].zAver;
        }
        
        for(int i=0;i<beamVec.size();i++)
        {     
            // get the time from current bunch to next bunch 
            if(i<beamVec.size()-1)
            {
                beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF - (beamVec[i].zAverLastTurn - beamVec[i+1].zAverLastTurn)/ CLight; 
            }
            else
            {
                beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF - (beamVec[i].zAverLastTurn - beamVec[0].zAver ) / CLight;
            }

            // get the time from last bunch to current bunch not used in simulation 
            if(i==0)
            {
                beamVec[i].timeFromLastBunchToCurrentBunch  = beamVec[beamVec.size()-1].bunchGap * tRF 
                                                            -(beamVec[beamVec.size()-1].zAverLastTurn - beamVec[i].zAverLastTurn) / CLight; 
            }
            else
            {
                beamVec[i].timeFromLastBunchToCurrentBunch  = beamVec[i-1             ].bunchGap * tRF 
                                                            -(beamVec[i -1            ].zAverLastTurn - beamVec[i].zAverLastTurn) / CLight;
            }

            // beamVec[i].BunchTransferDueToLatticeL(inputParameter,cavityResonator,turns);          
            // beamVec[i].BunchTransferDueToLatticeLMatarix(inputParameter,cavityResonator);
            beamVec[i].BunchTransferDueToLatticeLTest(inputParameter,cavityResonator);             
            beamVec[i].GetSPBunchRMS(latticeInterActionPoint, 0);   
        }
    }
    
        
}

void SPBeam::GetTimeDisToNextBunchIntial(ReadInputSettings &inputParameter)
{
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double circRing   = inputParameter.ringParBasic->circRing;
    double rfLen      = circRing  / ringHarmH;
    double temp=0;
    // get the bunch distance
    if(beamVec.size()>1)
    {
        for(int i=0;i<beamVec.size();i++)
        {
            // get the time from current bunch to next bunch
            if(i<beamVec.size()-1)
            {
                beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * rfLen - (beamVec[i].zAver - beamVec[i+1].zAver); 
            }
            else
            {
                beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * rfLen - (beamVec[i].zAver - beamVec[0  ].zAver ) ;
            }
            beamVec[i].timeFromCurrnetBunchToNextBunch /= CLight;     

            // get the time from last bunch to current bunch
            if(i==0)
            {
                beamVec[i].timeFromLastBunchToCurrentBunch  = beamVec[beamVec.size()-1].bunchGap * rfLen - (beamVec[beamVec.size()-1].zAver - beamVec[0].zAver ); 
            }
            else
            {
                beamVec[i].timeFromLastBunchToCurrentBunch  = beamVec[i-1             ].bunchGap * rfLen - (beamVec[i             -1].zAver - beamVec[i].zAver ) ;
            }

            beamVec[i].timeFromLastBunchToCurrentBunch /= CLight;  
        }
    }
    else
    {
        // if one bunch the inital bunch distance is set as ring circumference.
        beamVec[0].timeFromCurrnetBunchToNextBunch   = beamVec[0].bunchGap * rfLen ; 
        beamVec[0].timeFromCurrnetBunchToNextBunch  /= CLight;          
        beamVec[0].timeFromLastBunchToCurrentBunch   = beamVec[0].bunchGap * rfLen ; 
        beamVec[0].timeFromLastBunchToCurrentBunch  /= CLight;    
    }
}



void SPBeam::WSIonDataPrint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,int nTurns)
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
                string col= string("&parameter name=") + string("numOfMacroIon_") + to_string(p) + string(",  type=long,  &end");   
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


void SPBeam::BeamTransferPerInteractionPointDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToLatticeT(latticeInterActionPoint,k);
    }
}

void SPBeam::BeamTransferPerTurnDueToLattice(LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter,CavityResonator &cavityResonator,int turns)
{
    BeamTransferPerTurnDueToLatticeT(latticeInterActionPoint);
    BeamTransferPerTurnDueToLatticeL(inputParameter,latticeInterActionPoint,cavityResonator,turns); 
}

void SPBeam::BeamTransferPerTurnDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint)
{
    
    for(int k=0;k<latticeInterActionPoint.numberOfInteraction;k++)
    {
        BeamTransferPerInteractionPointDueToLatticeT(latticeInterActionPoint,k);
    }
}

void SPBeam::GetHaissinski(ReadInputSettings &inputParameter,CavityResonator &cavityResonator,WakeFunction &sRWakeFunction)
{
    vector<int> TBTBunchDisDataBunchIndex = inputParameter.ringRun->TBTBunchDisDataBunchIndex;
    int TBTBunchPrintNum  = inputParameter.ringRun->TBTBunchPrintNum;        
    int bunchIndex; 

    for(int i=0;i<TBTBunchPrintNum;i++)
    {        
        bunchIndex  = TBTBunchDisDataBunchIndex[i];
        beamVec[bunchIndex].GetBunchHaissinski(inputParameter,cavityResonator,sRWakeFunction);                     
    }
    string filePrefix = inputParameter.ringRun->TBTBunchHaissinski;
    string fname = filePrefix + ".sdds";
    ofstream fout(fname,ios_base::app);        
    
    fout<<"SDDS1"<<endl;
    for(int j=0; j<inputParameter.ringParRf->resNum;j++)
    {
        string parname = string("&parameter name=") + string("cavAmp_") + to_string(j) + string(", units=V,   type=float,  &end");
        fout<<parname<<endl;
        parname = string("&parameter name=") + string("cavPhase_") + to_string(j) + string(", units=rad,   type=float,  &end");
        fout<<parname<<endl;
    }   
    
    fout<<"&parameter name=bunchHarm,                       type=long,  &end"<<endl;
    fout<<"&parameter name=averZ,       units=m,            type=float,  &end"<<endl;
    fout<<"&parameter name=rmsZ,        units=m,            type=float,  &end"<<endl;
   
    fout<<"&column name=z,              units=m,            type=float,  &end"<<endl;
    fout<<"&column name=rfHamilton,                         type=float,  &end"<<endl;
    fout<<"&column name=wakeHamilton,                       type=float,  &end"<<endl;
    fout<<"&column name=totHamilton,                        type=float,  &end"<<endl;
    fout<<"&column name=rwWakePoten,                        type=float,  &end"<<endl;
    fout<<"&column name=bbrWakePoten,                        type=float,  &end"<<endl;
    fout<<"&column name=totWakePoten,                        type=float,  &end"<<endl;   
    fout<<"&column name=profile,       units=arb.units      type=float,  &end"<<endl;
    fout<<"&data mode=ascii, &end"<<endl;
    
    if(!inputParameter.ringRun->TBTBunchHaissinski.empty()  &&  (TBTBunchPrintNum !=0) )
    {        
        for(int i=0;i<TBTBunchPrintNum;i++)
        {
            bunchIndex  = TBTBunchDisDataBunchIndex[i];
            fout<<"! page number "<<i + 1<<endl;
        
            for(int j=0; j<inputParameter.ringParRf->resNum;j++)
            {
                fout<<beamVec[bunchIndex].haissinski->cavAmp[j]<<endl;
                fout<<beamVec[bunchIndex].haissinski->cavPhase[j]<<endl;              
            }             
            fout<<beamVec[bunchIndex].bunchHarmNum<<endl; 
            fout<<beamVec[bunchIndex].haissinski->averZ<<endl;                         
            fout<<beamVec[bunchIndex].haissinski->rmsZ<<endl;
            fout<<beamVec[bunchIndex].haissinski->nz<<endl;

            double norm=0;
            for(int k=0; k<beamVec[bunchIndex].haissinski->nz;k++)
            {
                norm += beamVec[bunchIndex].haissinski->bunchProfile[k];
            }       
            
            for(int k=0;k<beamVec[bunchIndex].haissinski->nz;k++)
            {
                fout<<setw(15)<<left<<beamVec[bunchIndex].haissinski->bunchPosZ[k]
                    <<setw(15)<<left<<beamVec[bunchIndex].haissinski->rfHamiltonian[k]
                    <<setw(15)<<left<<beamVec[bunchIndex].haissinski->wakeHamiltonian[k]
                    <<setw(15)<<left<<beamVec[bunchIndex].haissinski->totHamiltonian[k]
                    <<setw(15)<<left<<beamVec[bunchIndex].haissinski->rwWakePoten[k]
                    <<setw(15)<<left<<beamVec[bunchIndex].haissinski->bbrWakePoten[k]
                    <<setw(15)<<left<<beamVec[bunchIndex].haissinski->totWakePoten[k]
                    <<setw(15)<<left<<beamVec[bunchIndex].haissinski->bunchProfile[k]/norm             
                    <<endl;    
            }                                
        }
    }

    fout.close();    
}


void SPBeam::GetAnalyticalLongitudinalPhaseSpace(ReadInputSettings &inputParameter,CavityResonator &cavityResonator,WakeFunction &sRWakeFunction)
{
    
    vector<int> TBTBunchDisDataBunchIndex = inputParameter.ringRun->TBTBunchDisDataBunchIndex;
    int TBTBunchPrintNum  = inputParameter.ringRun->TBTBunchPrintNum;        
    int bunchIndex; 
    
    for(int i=0;i<TBTBunchPrintNum;i++)
    {        
        bunchIndex  = TBTBunchDisDataBunchIndex[i];
        beamVec[bunchIndex].GetParticleLongitudinalPhaseSpace(inputParameter,cavityResonator,i);                           
    }
}


void SPBeam::FIRBunchByBunchFeedback(FIRFeedBack &firFeedBack,int nTurns)
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

void SPBeam::BeamSynRadDamping(const ReadInputSettings &inputParameter, const LatticeInterActionPoint &latticeInterActionPoint)
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchSynRadDamping(inputParameter,latticeInterActionPoint);
    }
}

void SPBeam::LRWakeBeamIntaction(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction,const  LatticeInterActionPoint &latticeInterActionPoint,int turns)
{

    int nTurnswakeTrunction     = inputParameter.ringLRWake->nTurnswakeTrunction;
    int harmonics               = inputParameter.ringParBasic->harmonics;
    double electronBeamEnergy   = inputParameter.ringParBasic->electronBeamEnergy;
    double rBeta                = inputParameter.ringParBasic->rBeta;
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double circRing   = inputParameter.ringParBasic->circRing;


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
    vector<double> wakeStasticForceTemp(3,0);

    double tauij=0.e0;
    double tauijStastic=0.e0;
    int nTauij=0;
    double deltaZij=0;
    double tRF   = inputParameter.ringParBasic->t0 / double(harmonics);
    double rfLen = circRing  / ringHarmH;

    int tempIndex0,tempIndex1;

    // bunch j in the wittness particle during the simulation
    for (int j=0;j<beamVec.size();j++)
    {
        beamVec[j].lRWakeForceAver[0] =0.E0;      // x rad
        beamVec[j].lRWakeForceAver[1] =0.E0;      // y rad
        beamVec[j].lRWakeForceAver[2] =0.E0;      // z rad
        
        for(int n=0;n<nTurnswakeTrunction;n++)
        {
            // self-interation of bunch in current turn is excluded if tempIndex1 = j - 1, when n=0.             
            if(n==0)
            {
                tempIndex0   = 0;
                tempIndex1   = j - 1;
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
            // bunch i in the leadng particle in previous turn during the simulation --Alex eq.4.3
            // dLij =  -nkL_rf + z_j - z_i
            //      =  -nkT_rf - dz_i /c + dz_j / c 
            //      =  -nkT_rf - (dz_i - dz_j ) / c
            for(int i=tempIndex0;i<=tempIndex1;i++)
            {                 
                // Alex Chao Eq.(4.3) -- W'(k C - n C + z_n -z_k) i(k): leading bunch; j(n): wittness bunch
                deltaZij = wakefunction.poszData[nTurnswakeTrunction-1-n][i] - beamVec[j].zAver;
                nTauij   = beamVec[i].bunchHarmNum - beamVec[j].bunchHarmNum - n * harmonics;
                tauijStastic = nTauij * tRF;
                if(beamVec.size()==1)
                {
                    tauij        = nTauij * tRF  - deltaZij / CLight;   
                }
                else
                {
                    tauij        = nTauij * tRF  + deltaZij / CLight;   
                }
                
                // notification: 
                // ensure the wakefucntion return the focusing strength in transverse and energy loss in longitudinal. 
                // then: beamVec[j].lRWakeForceAver[?] -=  mins here.   

                if(!inputParameter.ringLRWake->pipeGeoInput.empty())
                {                                       
                    wakeForceTemp = wakefunction.GetRWLRWakeFun(tauij); 
                    beamVec[j].lRWakeForceAver[0] -= wakeForceTemp[0] * beamVec[i].electronNumPerBunch * wakefunction.posxData[nTurnswakeTrunction-1-n][i];    //[V/C/m] * [m] ->  [V/C]
                    beamVec[j].lRWakeForceAver[1] -= wakeForceTemp[1] * beamVec[i].electronNumPerBunch * wakefunction.posyData[nTurnswakeTrunction-1-n][i];    //[V/C/m] * [m] ->  [V/C]
                    beamVec[j].lRWakeForceAver[2] -= wakeForceTemp[2] * beamVec[i].electronNumPerBunch;                                                       //[V/C]         ->  [V/C]                    
                }

                if(!inputParameter.ringLRWake->bbrInput.empty())
                {
                    wakeForceTemp = wakefunction.GetBBRWakeFun(tauij);                                       
                    beamVec[j].lRWakeForceAver[0] -= wakeForceTemp[0] * beamVec[i].electronNumPerBunch * wakefunction.posxData[nTurnswakeTrunction-1-n][i];    //[V/C/m] * [m] ->  [V/C]
                    beamVec[j].lRWakeForceAver[1] -= wakeForceTemp[1] * beamVec[i].electronNumPerBunch * wakefunction.posyData[nTurnswakeTrunction-1-n][i];    //[V/C/m] * [m] ->  [V/C]                   
                    beamVec[j].lRWakeForceAver[2] -= wakeForceTemp[2] * beamVec[i].electronNumPerBunch ; 
                    // to subsctract the "stastic term"
                    // wakeStasticForceTemp = wakefunction.GetBBRWakeFun(tauijStastic);
                    // beamVec[j].lRWakeForceAver[2] -= (wakeForceTemp[2] - wakeStasticForceTemp[2]) * beamVec[i].electronNumPerBunch;                                                        //[V/C]         ->  [V/C]                
                }             
            }   
              
        }
        
        
        beamVec[j].lRWakeForceAver[0] *=  ElectronCharge / electronBeamEnergy / pow(rBeta,2);   // [V/C] * [C] * [1e] / [eV] ->rad
        beamVec[j].lRWakeForceAver[1] *=  ElectronCharge / electronBeamEnergy / pow(rBeta,2);   // [V/C] * [C] * [1e] / [eV] ->rad
        beamVec[j].lRWakeForceAver[2] *=  ElectronCharge / electronBeamEnergy / pow(rBeta,2);   // [V/C] * [C] * [1e] / [eV] ->rad

    }


    // Reference to "simulation of transverse multi-bunch instabilities of proton beam in LHC, PHD thesis, Alexander Koshik P. 32, Eq. (3.22)"
    // with BBR cavity model, to get the coupled bunch mode growth rate, the potential well distortation has to be substrcted.
    // tauijStastic = nTauij * tRF;
    // wakeStasticForceTemp = wakefunction.GetBBRWakeFun(tauijStastic); 
    // beamVec[j].lRWakeForceAver[2] -= (wakeForceTemp[2] - wakeStasticForceTemp[2] )* beamVec[i].electronNumPerBunch;
} 
void SPBeam::BeamTransferPerTurnDueWake()
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToWake();
    }
}