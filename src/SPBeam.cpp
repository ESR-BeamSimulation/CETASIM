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
    freXIQDecompScan.resize(totBunchNum);
    freYIQDecompScan.resize(totBunchNum);
    freZIQDecompScan.resize(totBunchNum);

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
    
    // set the different bunch charge according to index read in
    for(int i=0;i<inputParameter.ringFillPatt->bunchChargeNum;i++)
    {
        int index = inputParameter.ringFillPatt->bunchChargeIndex[i];
        beamVec[index].normCurrent = inputParameter.ringFillPatt->bunchCharge[i];
    }
    int totNormCurrent=0;
    for(int i=0;i<totBunchNum;i++)
    {
        totNormCurrent += beamVec[i].normCurrent; 
    }

    for(int i=0;i<totBunchNum;i++)
    {
        beamVec[i].current = beamVec[i].normCurrent * inputParameter.ringParBasic->ringCurrent / totNormCurrent;
    }
    //---------------------------------------------------------------------------------------

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
    GetTimeDisToNextBunch(inputParameter);

    // set the IQDecompose scan frequency
    double workQx = inputParameter.ringParBasic->workQx;
    double workQy = inputParameter.ringParBasic->workQy;
    double workQz = inputParameter.ringParBasic->workQz;
    double f0       = inputParameter.ringParBasic->f0;

    // workQx = workQx - floor(workQx);
    // workQy = workQy - floor(workQy);
    // workQz = workQz - floor(workQz);

    // only works for unifrom filling patterns setings.
    for(int i=0;i<totBunchNum;i++)
    {
        freXIQDecompScan[i] = (i + workQx) * f0;
        freYIQDecompScan[i] = (i + workQy) * f0;
        freZIQDecompScan[i] = (i + workQz) * f0;
    }
    
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
    // double beamCurr = inputParameter.ringBunchPara->current * beamVec.size();
    double beamCurr = inputParameter.ringParBasic->ringCurrent ;
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

                vbAccum = vbAccum * exp(- deltaL ) * exp ( li * cPsi);    // decay and rotate [V] P.B. Eq.(3.12) get the beam induced voltage at the cavity

                if( n==nTurns-1)
                {
                    vbAccumAver2 +=  vbAccum;
                }
            }
        }
     
        
        vbAccumAver =  (vbAccumAver1 + vbAccumAver2) /2.0/ double(beamVec.size());   // (1)
        vbAccumAver =  vbAccumAver1/double(beamVec.size());                          // before decay and rotate, just  after bunch left.
        // vbAccumAver =  vbAccumAver2/double(beamVec.size());                       // after decay and rotate.
        // DC solution -- which is exactly the same as above equaiton (1)'s reuslts, very good agreement.
        vbAccumAver = 2.0 * beamCurr * cavityResonator.resonatorVec[i].resShuntImpRs / (1.0 + cavityResonator.resonatorVec[i].resCouplingBeta)
                                    * cos(cavityResonator.resonatorVec[i].resDeTunePsi) * exp(li * (PI + cavityResonator.resonatorVec[i].resDeTunePsi)  );
        

        cavityResonator.resonatorVec[i].vbAccum   =  vbAccumAver;               // set cavity beam induced voltage. It is the value particle feel
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
            beamVec[j].bunchRFModeInfo->cavVolBunchCen[i]    =     cavityResonator.resonatorVec[i].resCavVolReq;     // set vol as the requried voltage.
            beamVec[j].bunchRFModeInfo->induceVolBunchCen[i] =     cavityResonator.resonatorVec[i].vbAccum;
            beamVec[j].bunchRFModeInfo->genVolBunchAver[i]   =     cavityResonator.resonatorVec[i].resGenVol;
            beamVec[j].bunchRFModeInfo->selfLossVolBunchCen[i] =    vbKickAver;
        }
    }
    
    cavityResonator.GetIntialGenIg();
    cavityResonator.GetInitalCavityPowerInfo(inputParameter);

    // print out the data to show the cavity voltage buildup process.  
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
    fout<<"&parameter name=GenAmpIg,         units=Ampere, type=float,  &end"<<endl;
    fout<<"&parameter name=GenArgIg,         units=rad, type=float,  &end"<<endl;
    fout<<"&parameter name=LoadAnlge,        units=rad, type=float,  &end"<<endl;
    fout<<"&parameter name=BeamPower,        units=watt,type=float,  &end"<<endl;
    fout<<"&parameter name=GenPower,         units=watt,type=float,  &end"<<endl;
    fout<<"&parameter name=CavPower,         units=watt,type=float,  &end"<<endl;
    fout<<"&parameter name=GenReflectPower,  units=watt,type=float,  &end"<<endl;

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
    for(int i=0;i<resNum;i++)
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
        fout<<abs(cavityResonator.resonatorVec[i].resGenIg)               <<endl;
        fout<<arg(cavityResonator.resonatorVec[i].resGenIg)               <<endl;
        fout<<arg(cavityResonator.resonatorVec[i].resGenIg)  - arg(cavityResonator.resonatorVec[i].resCavVolReq)   <<endl;    
        fout<<    cavityResonator.resonatorVec[i].resBeamPower            <<endl;
        fout<<    cavityResonator.resonatorVec[i].resGenPower             <<endl;   
        fout<<    cavityResonator.resonatorVec[i].resCavPower             <<endl;      
        fout<<    cavityResonator.resonatorVec[i].resGenPowerReflect      <<endl; 

        fout<<"! page number "<<i+1<<endl;
        fout<<nTurns*beamVec.size()<<endl;

        
        for(int n=0;n<nTurns;n++)
        {
            for(int j=0;j<beamVec.size();j++)
            {
                vb0  = complex<double>(-1 * cavityResonator.resonatorVec[i].resFre * 2 * PI * cavityResonator.resonatorVec[i].resShuntImpRs /
                                            cavityResonator.resonatorVec[i].resQualityQ0,0.E0) * beamVec[j].electronNumPerBunch * ElectronCharge;
                vbAccum += vb0;                                       //[V]

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

                tB = beamVec[j].bunchGap * t0 / ringHarmH ;  //[s]
                time +=  tB;
                deltaL = tB / tF ;
                // cPsi = 2 * PI * (cavityResonator.resonatorVec[i].resFre - ringHarmH * f0 * cavityResonator.resonatorVec[i].resHarm) * tB; 
                cPsi = 2 * PI * cavityResonator.resonatorVec[i].resFre * tB;
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
            // BeamTransferPerTurnDueToLatticeT(latticeInterActionPoint);       // here use the one turn matrix with elegant apprpach     
            BeamTransferPerTurnDueToLatticeTOneTurnR66(inputParameter,latticeInterActionPoint);
        }
        
        SPBeamRMSCal(latticeInterActionPoint, 0);
        // here tracking partilce in longitudinal due to the RF filed from cavity
        // BeamTransferPerTurnDueToLatticeL(inputParameter,latticeInterActionPoint,cavityResonator,n); 
        BeamTransferPerTurnDueToLatticeLTest(inputParameter,latticeInterActionPoint,cavityResonator,n); 
        
        if(lRWakeFlag)
        {
            SPBeamRMSCal(latticeInterActionPoint, 0);
            LRWakeBeamIntaction(inputParameter,lRWakeFunction,latticeInterActionPoint,n);  
        }
        
        if(fIRBunchByBunchFeedbackFlag)  
        {     
            SPBeamRMSCal(latticeInterActionPoint, 0);
            FIRBunchByBunchFeedback(inputParameter,firFeedBack,n);
        }

        if(synRadDampingFlag[0]==1) // transverse SR effect
        {
            SPBeamRMSCal(latticeInterActionPoint, 0);
            BeamSynRadDamping(inputParameter,latticeInterActionPoint);
        }

        if((inputParameter.driveMode->driveModeOn!=0) && (inputParameter.driveMode->driveStart <n )  &&  (inputParameter.driveMode->driveEnd >n) )
        {
            SPBeamRMSCal(latticeInterActionPoint, 0);
            BeamTransferDuetoDriveMode(inputParameter,n);
        }   

        SPBeamRMSCal(latticeInterActionPoint, 0);
        SPGetBeamInfo();
                    
        if(bunchInfoPrintInterval && (n%bunchInfoPrintInterval==0) )
        {      
            SPBeamDataPrintPerTurn(n,latticeInterActionPoint,inputParameter,cavityResonator);           
            
            if(inputParameter.driveMode->driveModeOn!=0)
            {
                GetDriveModeGrowthRate(n,inputParameter);
            }
            if(!inputParameter.ringRun->runCBMGR.empty())
            {
                GetCBMGR(n,latticeInterActionPoint,inputParameter);
            }
        }

        MarkParticleLostInBunch(inputParameter,latticeInterActionPoint);

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

void SPBeam::BeamTransferPerTurnDueToLatticeTOneTurnR66(const ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint)
{
    for(int i=0;i<beamVec.size();i++)
    {
        beamVec[i].BunchTransferDueToLatticeOneTurnT66(inputParameter,latticeInterActionPoint);
    }
}


void SPBeam::MarkParticleLostInBunch(const ReadInputSettings &inputParameter, const LatticeInterActionPoint &latticeInterActionPoint)
{
    for(int i=0;i<beamVec.size();i++)
    {
        beamVec[i].MarkLostParticle(inputParameter,latticeInterActionPoint);
    }
}


void SPBeam::BeamTransferDuetoDriveMode(const ReadInputSettings &inputParameter, const int n)
{     
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToDriveMode(inputParameter,n);     
    }
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

void SPBeam::GetHilbertAnalyticalInOneTurn(const ReadInputSettings &inputParameter)
{    
    double workQx = inputParameter.ringParBasic->workQx;
    double workQy = inputParameter.ringParBasic->workQy;
    double workQz = inputParameter.ringParBasic->workQz;
    int harmonics = inputParameter.ringParBasic->harmonics;
    
    vector<double> xSignal(beamVec.size());
    vector<double> ySignal(beamVec.size());
    vector<double> zSignal(beamVec.size());

    double phasex, phasey, phasez;

    for(int i=0;i<beamVec.size();i++)
	{
        phasex = - 2.0 * PI * workQx * beamVec[i].bunchHarmNum / harmonics;
        phasey = - 2.0 * PI * workQy * beamVec[i].bunchHarmNum / harmonics;
        phasez = - 2.0 * PI * workQz * beamVec[i].bunchHarmNum / harmonics;

        xSignal[i] = beamVec[i].xAver; //* exp (li * phasex).real();
        ySignal[i] = beamVec[i].yAver; //* exp (li * phasey).real();
        zSignal[i] = beamVec[i].zAver; //* exp (li * phasez).real();

    }
    vector<complex<double> > xAnalytical = GetHilbertAnalytical(xSignal,0.01,workQx);
    vector<complex<double> > yAnalytical = GetHilbertAnalytical(ySignal,0.01,workQy);
    vector<complex<double> > zAnalytical = GetHilbertAnalytical(zSignal,0.01,workQz);    

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

vector<complex<double> > SPBeam::GetHilbertAnalytical(vector<double> signal, const double filterBandWithdNu,  double workQ)
{   
    // with hilbert transform to get analytical signal Bessy in VSR
    // Ref. Calculation of coupled bunch effects in the synchrotorn loight source bessy VSR
    
    workQ -= floor(workQ);
    if(workQ>0.5) workQ = 1 - workQ;  

    int nuCenter = workQ * signal.size();
    double sigmaNu  = filterBandWithdNu  * signal.size();

    // generate a gaussian filter according to the working point

    vector<double> gaussFilter(signal.size(),0.E0);   
    for(int i=0;i<signal.size();i++)
    {
        gaussFilter[i] =  exp(- 1.0 / 2.0 * pow( (i - nuCenter) / sigmaNu ,2) );
    }
 
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
        xFFT[2 * index    ] =   xFFT[2 * index   ]  / 2.0;
        xFFT[2 * index + 1] =   xFFT[2 * index + 1] / 2.0;
    }

    for (int i=signal.size()/2 +1;i<signal.size();i++)
    {
        xFFT[2*i    ] = 0.E0;
        xFFT[2*i + 1] = 0.E0;
    }

    // for(int i=0;i<signal.size();i++)
    // {
    //     xFFT[2*i    ] *= gaussFilter[i];       
    //     xFFT[2*i + 1] *= gaussFilter[i];
    // }

    // ofstream fout("test.sdds");   
    // for(int i=0;i<signal.size();i++)
    // {
    //     fout<<setw(15)<<left<<i
    //         <<setw(15)<<left<<xFFT[2*i]
    //         <<setw(15)<<left<<xFFT[2*i+1]
    //         <<setw(15)<<left<<sqrt(pow(xFFT[2*i],2) + pow(xFFT[2*i+1],2) )
    //         <<setw(15)<<left<<gaussFilter[i]
    //         <<endl;
    // }    
    // fout.close();
    // cout<<__LINE__<<endl;
    // getchar();

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


void SPBeam::GetDriveModeGrowthRate(const int turns, const ReadInputSettings &inputParameter)
{
   
    double time = 0.E0;
    int ringHarm        = inputParameter.ringParRf->ringHarm;
    double f0           = inputParameter.ringParBasic->f0;
    double t0           = inputParameter.ringParBasic->t0;
    double rBeta        = inputParameter.ringParBasic->rBeta;
    double tRF          = t0 / ringHarm;
    double driveAmp     = inputParameter.driveMode->driveAmp;
    double driveFre     = inputParameter.driveMode->driveFre;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    double workQ;
    vector<complex<double> > signalIQ(beamVec.size());

    if(inputParameter.driveMode->drivePlane == 0)   // IQ signale decomposition...
    {
        workQ = inputParameter.ringParBasic->workQz;
        for(int i=0;i<beamVec.size();i++)
        {
            time =  beamVec[i].bunchHarmNum * tRF ;
            signalIQ[i] = beamVec[i].zAver  * exp( - li * 2.0 * PI * driveFre * time) ;
        }
    }
    else if (inputParameter.driveMode->drivePlane == 1)
    {
        workQ = inputParameter.ringParBasic->workQx;
        for(int i=0;i<beamVec.size();i++)
        {
            time =  beamVec[i].bunchHarmNum * tRF ;
            signalIQ[i] = beamVec[i].xAver  * exp( - li * 2.0 * PI * driveFre * time) ;
        }
    }
    else if (inputParameter.driveMode->drivePlane == 2)
    {
        workQ = inputParameter.ringParBasic->workQy;
        for(int i=0;i<beamVec.size();i++)
        {
            time =  beamVec[i].bunchHarmNum * tRF ;
            signalIQ[i] = beamVec[i].yAver  * exp( - li * 2.0 * PI * driveFre * time);
        }
    }

    complex<double> signalIQAver =  complex<double>(0.0,0.0);

    if (inputParameter.driveMode->driveHW !=0)
    {
        signalIQAver = accumulate(signalIQ.begin(), signalIQ.end(), complex<double>(0.0,0.0) ) / double(beamVec.size());
    }
    else
    {            
        double weigh=0.E0;
        double weighSum=0.E0;

        for(int i=0;i<beamVec.size();i++)
        {
            weigh =  pow( sin( i * PI / beamVec.size() ), 2);
            signalIQAver += weigh * signalIQ[i];
            weighSum     += weigh;
        }        
        signalIQAver /= weighSum;
    }
    
    ampIQ.push_back(abs(signalIQAver)); 
    phaseIQ.push_back(arg(signalIQAver));  


    if( (turns + inputParameter.ringRun->bunchInfoPrintInterval)  == inputParameter.ringRun->nTurns   )
    {        
        int indexStart =  int(inputParameter.ringRun->growthRateFittingStart / inputParameter.ringRun->bunchInfoPrintInterval) ;
        int indexEnd   =  int(inputParameter.ringRun->growthRateFittingEnd   / inputParameter.ringRun->bunchInfoPrintInterval);  
        int dim        =  indexEnd - indexStart ;


        double fitWeight[dim];
        double fitX[dim];
        double fitY[dim];
        vector<double> resFit(2,0.E0);
        double modePhase = 0.E0;
        

        for(int n=0;n<dim;n++)
        {                
            fitWeight[n]        = 1.E0;
            fitX[n]             = n * inputParameter.ringRun->bunchInfoPrintInterval;
            fitY[n]             = log(ampIQ[n+indexStart]);
            modePhase          += phaseIQ[n+indexStart];
        }
        modePhase /= dim ; 
    
        FittingGSL fittingGSL;  
        resFit = fittingGSL.FitALinear(fitX,fitWeight,fitY,dim);
                    
        double decimal = driveFre / f0 / beamVec.size() - floor( driveFre / f0 / beamVec.size() ) ;

        ofstream fout("driveModeGrowthRate.sdds");
        fout<<"SDDS1"                                                                    <<endl;
        fout<<"&parameter name=ModeIndex                              type=float    &end"<<endl;
        fout<<"&parameter name=ModeDriveFre,            units=Hz,     type=float    &end"<<endl;
        fout<<"&parameter name=ModeGrowthRate,          units=1/s,    type=float    &end"<<endl;
        fout<<"&parameter name=ModePhase,               units=rad,    type=float    &end"<<endl;
        fout<<"&parameter name=ModeFre,                 units=1/s,    type=float    &end"<<endl;

        fout<<"&column name=Turns,                                    type=long     &end"<<endl;
        fout<<"&column name=ampIQ,                                    type=float     &end"<<endl;
        fout<<"&column name=phaseIQ,                                  type=float     &end"<<endl;
        fout<<"&data mode=ascii &end"<<endl;          
        fout<<"! page number 1"<<endl;
        
        fout<<setw(15)<<left<<floor(decimal * beamVec.size())<<endl;
        fout<<setw(15)<<left<<driveFre<<endl; 
        fout<<setw(15)<<left<<resFit[1] / inputParameter.ringParBasic->t0<<endl; 
        fout<<setw(15)<<left<<modePhase<<endl;  
        fout<<setw(15)<<left<<modePhase / (2 * PI * tRF * ringHarm / beamVec.size()  ) <<endl;         
        fout<<ampIQ.size()<<endl;

        for(int i=0;i<ampIQ.size();i++)
        {
            fout<<setw(15)<<left<<i * inputParameter.ringRun->bunchInfoPrintInterval
                <<setw(15)<<left<<log(ampIQ[i])
                <<setw(15)<<left<<phaseIQ[i]<<endl;
        }    
        fout.close();
    }
   
}

void SPBeam::GetCBMGR(const int turns, const LatticeInterActionPoint &latticeInterActionPoint, const ReadInputSettings &inputParameter)
{
    double alphax = latticeInterActionPoint.twissAlphaX[0];
    double betax  = latticeInterActionPoint.twissBetaX [0];
    double alphay = latticeInterActionPoint.twissAlphaY[0];
    double betay  = latticeInterActionPoint.twissBetaY[0] ;
    double alphaz = 0;
    double betaz  = inputParameter.ringBunchPara->rmsBunchLength / inputParameter.ringBunchPara->rmsEnergySpread;
    double workQx = inputParameter.ringParBasic->workQx;
    double workQy = inputParameter.ringParBasic->workQy;
    double workQz = inputParameter.ringParBasic->workQz;
    int harmonics = inputParameter.ringParBasic->harmonics;
    double tRF    = inputParameter.ringParBasic->t0 / harmonics; 

    gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (beamVec.size());
	workspace = gsl_fft_complex_workspace_alloc (beamVec.size());

    double tempXFFT[2*beamVec.size()];
    double tempYFFT[2*beamVec.size()];
    double tempZFFT[2*beamVec.size()];

    // (1) IPAC 2022 WEPOMS010 -- Diamond-II -- WangSiWei's IPAC paper. 
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
        double phasex = - 2.0 * PI * workQx * beamVec[i].bunchHarmNum / harmonics;
        double phasey = - 2.0 * PI * workQy * beamVec[i].bunchHarmNum / harmonics;
        double phasez = - 2.0 * PI * workQz * beamVec[i].bunchHarmNum / harmonics;

        complex<double> zx = ( x / sqrt(betax) - li * ( sqrt(betax) * px + alphax / sqrt(betax) * x ) ) * exp (li * phasex);
        complex<double> zy = ( y / sqrt(betay) - li * ( sqrt(betay) * py + alphay / sqrt(betay) * y ) ) * exp (li * phasey); 
        complex<double> zz = ( z / sqrt(betaz) + li * ( sqrt(betaz) * pz + alphaz / sqrt(betaz) * z ) ) * exp (li * phasez);  

        tempXFFT[2*i]  = zx.real();
        tempYFFT[2*i]  = zy.real();
        tempZFFT[2*i]  = zz.real();
        tempXFFT[2*i+1]= zx.imag();
        tempYFFT[2*i+1]= zy.imag();
        tempZFFT[2*i+1]= zz.imag();
    }

    gsl_fft_complex_forward (tempXFFT,1, beamVec.size(), wavetable, workspace);
    gsl_fft_complex_forward (tempYFFT,1, beamVec.size(), wavetable, workspace);
    gsl_fft_complex_forward (tempZFFT,1, beamVec.size(), wavetable, workspace);
    

    // store the data after fft turn by turn for fitting
    vector<double> tempXFFTAmp(beamVec.size());
    vector<double> tempYFFTAmp(beamVec.size());
    vector<double> tempZFFTAmp(beamVec.size());
    vector<double> tempXFFTArg(beamVec.size());
    vector<double> tempYFFTArg(beamVec.size());
    vector<double> tempZFFTArg(beamVec.size());
        
    for (int i=0;i<beamVec.size();i++)
    {
        tempXFFTAmp[i] = sqrt( pow(tempXFFT[2*i],2) + pow(tempXFFT[2*i+1],2) );
        tempYFFTAmp[i] = sqrt( pow(tempYFFT[2*i],2) + pow(tempYFFT[2*i+1],2) );
        tempZFFTAmp[i] = sqrt( pow(tempZFFT[2*i],2) + pow(tempZFFT[2*i+1],2) );
        tempXFFTArg[i] = atan2(tempXFFT[2*i], tempXFFT[2*i +1]);
        tempYFFTArg[i] = atan2(tempYFFT[2*i], tempYFFT[2*i +1]);
        tempZFFTArg[i] = atan2(tempZFFT[2*i], tempZFFT[2*i +1]);
    }

    coupledBunchModeAmpX.push_back(tempXFFTAmp);
    coupledBunchModeAmpY.push_back(tempYFFTAmp);
    coupledBunchModeAmpZ.push_back(tempZFFTAmp);
    coupledBunchModeArgX.push_back(tempXFFTArg);
    coupledBunchModeArgY.push_back(tempYFFTArg);
    coupledBunchModeArgZ.push_back(tempZFFTArg);

    //(2) Store IQ decomposition without excitation-- only show the unstable mode   
    vector<double> xIQAmpOneTurn(beamVec.size());
    vector<double> yIQAmpOneTurn(beamVec.size());
    vector<double> zIQAmpOneTurn(beamVec.size());
    vector<double> xIQArgOneTurn(beamVec.size());
    vector<double> yIQArgOneTurn(beamVec.size());
    vector<double> zIQArgOneTurn(beamVec.size());
                
    for(int j=0;j<freXIQDecompScan.size();j++)
    {
        vector<complex<double> > xIQ(beamVec.size());
        vector<complex<double> > yIQ(beamVec.size());
        vector<complex<double> > zIQ(beamVec.size());
        complex<double> xIQAver,yIQAver,zIQAver;
        
        double time = 0.0;
        double phaseX,phaseY,phaseZ;
                
        for(int i=0;i<beamVec.size();i++)
        {
            time =  beamVec[i].bunchHarmNum * tRF ;

            phaseX = - 2.0 * PI * workQx * beamVec[i].bunchHarmNum / harmonics;
            phaseY = - 2.0 * PI * workQy * beamVec[i].bunchHarmNum / harmonics;
            phaseZ = - 2.0 * PI * workQz * beamVec[i].bunchHarmNum / harmonics;

            xIQ[i] = beamVec[i].xAver  * exp( - li * 2.0 * PI * freXIQDecompScan[j] * time) ;
            yIQ[i] = beamVec[i].yAver  * exp( - li * 2.0 * PI * freYIQDecompScan[j] * time) ;
            zIQ[i] = beamVec[i].zAver  * exp( - li * 2.0 * PI * freZIQDecompScan[j] * time) ;            
        }

        xIQAver =   accumulate(xIQ.begin(), xIQ.end(), complex<double>(0.0,0.0) ) / double(beamVec.size());
        yIQAver =   accumulate(yIQ.begin(), yIQ.end(), complex<double>(0.0,0.0) ) / double(beamVec.size());
        zIQAver =   accumulate(zIQ.begin(), zIQ.end(), complex<double>(0.0,0.0) ) / double(beamVec.size());       

        xIQAmpOneTurn[j] = abs(xIQAver);
        yIQAmpOneTurn[j] = abs(yIQAver);
        zIQAmpOneTurn[j] = abs(zIQAver);

        xIQArgOneTurn[j] = arg(xIQAver);
        yIQArgOneTurn[j] = arg(yIQAver);
        zIQArgOneTurn[j] = arg(zIQAver);
    }

    ampXIQ.push_back(xIQAmpOneTurn);
    ampYIQ.push_back(yIQAmpOneTurn);
    ampZIQ.push_back(zIQAmpOneTurn);
    argXIQ.push_back(xIQArgOneTurn);
    argYIQ.push_back(yIQArgOneTurn); 
    argZIQ.push_back(zIQArgOneTurn);      
    //IQ decomposition --------------------------------


    //(3)  store beam pos data for bunch-by-bunch Growth rate calculation
    vector<double > xOneTurn(beamVec.size());
    vector<double > yOneTurn(beamVec.size());
    vector<double > zOneTurn(beamVec.size());
    for(int i=0;i<beamVec.size();i++)
    {
        xOneTurn[i] =  beamVec[i].xAver;
        yOneTurn[i] =  beamVec[i].yAver;
        zOneTurn[i] =  beamVec[i].zAver;
    }
    historyAverX.push_back(xOneTurn);
    historyAverY.push_back(yOneTurn);
    historyAverZ.push_back(zOneTurn);

    // Ideal method agrees with Analytical method. 
    // To get stable and unstbale coupled bunch mode grwoth-- have to rebuild the (x-px) along the ring.
    // the px info can be re-build from analyticla single from Hilbert transform--togethe with a gaussian filter. 
    // The IO approach with exictiation only give the unstbale growthrate. 
    // if no px info is generated.--- cannnot get the grwoth rate if there is no excitation. 
    // 
    if( (turns + inputParameter.ringRun->bunchInfoPrintInterval)  == inputParameter.ringRun->nTurns   )
    {               
        GetAnalyticalWithFilter(inputParameter);
    
        int indexStart =  int(inputParameter.ringRun->growthRateFittingStart / inputParameter.ringRun->bunchInfoPrintInterval);
        int indexEnd   =  int(inputParameter.ringRun->growthRateFittingEnd / inputParameter.ringRun->bunchInfoPrintInterval);  
        int dim        =  indexEnd - indexStart ;

        double fitWeight[dim];
        double fitX[dim];
        double fitFFTX[dim],    fitFFTY[dim],  fitFFTZ[dim];
        double fitAmpIQX[dim], fitAmpIQY[dim], fitAmpIQZ[dim];
        double fitAbsAverXAna[dim],fitAbsAverYAna[dim],fitAbsAverZAna[dim];
        double fitHilbertCBMX[dim],fitHilbertCBMY[dim],fitHilbertCBMZ[dim];
        vector<double> resFitFFTX(2,0.E0), resFitFFTY(2,0.E0), resFitFFTZ(2,0.E0);
        vector<double> resFitIQX(2,0.E0),  resFitIQY(2,0.E0),  resFitIQZ(2,0.E0);
        vector<double> resFitAbsAverXAna(2,0.E0),  resFitAbsAverYAna(2,0.E0),  resFitAbsAverZAna(2,0.E0);
        vector<double> resfitHilbertCBMX(2,0.E0),  resfitHilbertCBMY(2,0.E0),  resfitHilbertCBMZ(2,0.E0);
        
        // get the CBMAmp from hilbert transform....
        for(int n=0;n<ampXIQ.size();n++)
        {
            double time = 0.0;
            double phasex,phasey,phasez;
            double x,y,z,px,py,pz;
            
            for(int i=0; i<beamVec.size(); i++)
            {
                // tempXFFT[2*i  ] = anaSignalAverX[n][i].real();
                // tempXFFT[2*i+1] = anaSignalAverX[n][i].imag();
                // tempYFFT[2*i  ] = anaSignalAverY[n][i].real();
                // tempYFFT[2*i+1] = anaSignalAverY[n][i].imag();
                // tempZFFT[2*i  ] = anaSignalAverY[n][i].real();
                // tempZFFT[2*i+1] = anaSignalAverZ[n][i].imag();

                x  = anaSignalAverX[n][i].real();
                y  = anaSignalAverY[n][i].real();
                z  = anaSignalAverZ[n][i].real();
                px = -anaSignalAverX[n][i].imag();
                py = -anaSignalAverY[n][i].imag();
                pz = -anaSignalAverZ[n][i].imag();

                phasex = - 2.0 * PI * workQx * beamVec[i].bunchHarmNum / harmonics;
                phasey = - 2.0 * PI * workQy * beamVec[i].bunchHarmNum / harmonics;
                phasez = - 2.0 * PI * workQz * beamVec[i].bunchHarmNum / harmonics;
                
                complex<double> zx = ( x / sqrt(betax) - li * ( sqrt(betax) * px + alphax / sqrt(betax) * x ) ) * exp (li * phasex);
                complex<double> zy = ( y / sqrt(betay) - li * ( sqrt(betay) * py + alphay / sqrt(betay) * y ) ) * exp (li * phasey); 
                complex<double> zz = ( z / sqrt(betaz) + li * ( sqrt(betaz) * pz + alphaz / sqrt(betaz) * z ) ) * exp (li * phasez);

                tempXFFT[2*i]  = zx.real();
                tempYFFT[2*i]  = zy.real();
                tempZFFT[2*i]  = zz.real();
                tempXFFT[2*i+1]= zx.imag();
                tempYFFT[2*i+1]= zy.imag();
                tempZFFT[2*i+1]= zz.imag();
            }
            gsl_fft_complex_forward (tempXFFT,1, beamVec.size(), wavetable, workspace);
            gsl_fft_complex_forward (tempYFFT,1, beamVec.size(), wavetable, workspace);
            gsl_fft_complex_forward (tempZFFT,1, beamVec.size(), wavetable, workspace);

            for (int i=0;i<beamVec.size();i++)
            {
                tempXFFTAmp[i] = sqrt( pow(tempXFFT[2*i],2) + pow(tempXFFT[2*i+1],2) );
                tempYFFTAmp[i] = sqrt( pow(tempYFFT[2*i],2) + pow(tempYFFT[2*i+1],2) );
                tempZFFTAmp[i] = sqrt( pow(tempZFFT[2*i],2) + pow(tempZFFT[2*i+1],2) );
                tempXFFTArg[i] = atan2(tempXFFT[2*i], tempXFFT[2*i +1]);
                tempYFFTArg[i] = atan2(tempYFFT[2*i], tempYFFT[2*i +1]);
                tempZFFTArg[i] = atan2(tempZFFT[2*i], tempZFFT[2*i +1]);
            }

            hilbertCoupledBunchModeAmpX.push_back(tempXFFTAmp);
            hilbertCoupledBunchModeAmpY.push_back(tempYFFTAmp);
            hilbertCoupledBunchModeAmpZ.push_back(tempZFFTAmp);
            hilbertCoupledBunchModeArgX.push_back(tempXFFTArg);
            hilbertCoupledBunchModeArgY.push_back(tempYFFTArg);
            hilbertCoupledBunchModeArgZ.push_back(tempZFFTArg);

        }
        
        //////////////////////////////////////////////////////////////////////////
        
        FittingGSL fittingGSL; 
        ofstream fout(inputParameter.ringRun->runCBMGR+".sdds");
        fout<<"SDDS1"                                                                    <<endl;
        fout<<"&parameter name=dataPointsForFit                       type=long     &end"<<endl;
        fout<<"&parameter name=ModeIndex                              type=long     &end"<<endl;
        fout<<"&parameter name=CBMIdealGRX,                  units=1/s,    type=float    &end"<<endl;
        fout<<"&parameter name=CBMIdealGRY,                  units=1/s,    type=float    &end"<<endl;
        fout<<"&parameter name=CBMIdealGRZ,                  units=1/s,    type=float    &end"<<endl;
        fout<<"&parameter name=CBMIQGRX,                units=1/s,    type=float    &end"<<endl;
        fout<<"&parameter name=CBMIQGRY,                units=1/s,    type=float    &end"<<endl;
        fout<<"&parameter name=CBMIQGRZ,                units=1/s,    type=float    &end"<<endl;
        fout<<"&parameter name=CBMHTGRX,                units=1/s,    type=float    &end"<<endl;
        fout<<"&parameter name=CBMHTGRY,                units=1/s,    type=float    &end"<<endl;
        fout<<"&parameter name=CBMHTGRZ,                units=1/s,    type=float    &end"<<endl;

        fout<<"&parameter name=BunchHarmIndex                         type=long     &end"<<endl;
        fout<<"&parameter name=BGRX,                    units=1/s,    type=float    &end"<<endl;
        fout<<"&parameter name=BGRY,                    units=1/s,    type=float    &end"<<endl;
        fout<<"&parameter name=BGRZ,                    units=1/s,    type=float    &end"<<endl;


        fout<<"&column name=Turns,                                    type=long      &end"<<endl;
        fout<<"&column name=AbsCBMIdealX,                              type=float     &end"<<endl;
        fout<<"&column name=AbsCBMIdealY,                              type=float     &end"<<endl;
        fout<<"&column name=AbsCBMIdealZ,                              type=float     &end"<<endl;
        fout<<"&column name=ArgCBMIdealX,                units=rad,    type=float     &end"<<endl;
        fout<<"&column name=ArgCBMIdealY,                units=rad,    type=float     &end"<<endl;
        fout<<"&column name=ArgCBMIdealZ,                units=rad,    type=float     &end"<<endl;
        
        fout<<"&column name=AbsIQX,                                   type=float     &end"<<endl;
        fout<<"&column name=AbsIQY,                                   type=float     &end"<<endl;
        fout<<"&column name=AbsIQZ,                                   type=float     &end"<<endl;
        fout<<"&column name=ArgIQX,                     units=rad,    type=float     &end"<<endl;
        fout<<"&column name=ArgIQY,                     units=rad,    type=float     &end"<<endl;
        fout<<"&column name=ArgIQZ,                     units=rad,    type=float     &end"<<endl;

        fout<<"&column name=AbsCBMHilbertAnaX,         units=m,      type=float     &end"<<endl;
        fout<<"&column name=AbsCBMHilbertAnaY,         units=m,      type=float     &end"<<endl;
        fout<<"&column name=AbsCBMHilbertAnaZ,         units=m,      type=float     &end"<<endl;
        fout<<"&column name=ArgCBMHilbertAnaX,         units=rad,    type=float     &end"<<endl;
        fout<<"&column name=ArgCBMHilbertAnaY,         units=rad,    type=float     &end"<<endl;
        fout<<"&column name=ArgCBMHilbertAnaZ,         units=rad,    type=float     &end"<<endl;

        fout<<"&column name=AverX,                      units=m,      type=float     &end"<<endl;
        fout<<"&column name=AverY,                      units=m,      type=float     &end"<<endl;
        fout<<"&column name=AverZ,                      units=m,      type=float     &end"<<endl;
        fout<<"&column name=AnalyticalAverX,            units=m,      type=float     &end"<<endl;
        fout<<"&column name=AnalyticalAverY,            units=m,      type=float     &end"<<endl;
        fout<<"&column name=AnalyticalAverZ,            units=m,      type=float     &end"<<endl;
        // fout<<"&column name=AnalyticalAverPX,           units=rad,    type=float     &end"<<endl;
        // fout<<"&column name=AnalyticalAverPY,           units=rad,    type=float     &end"<<endl;
        // fout<<"&column name=AnalyticalAverPZ,           units=rad,    type=float     &end"<<endl;      
        fout<<"&column name=AbsHilbertAnaX,                 units=m,      type=float     &end"<<endl;
        fout<<"&column name=AbsHilbertAnaY,                 units=m,      type=float     &end"<<endl;
        fout<<"&column name=AbsHilbertAnaZ,                 units=m,      type=float     &end"<<endl;        
        fout<<"&data mode=ascii &end"<<endl;          
            
        for(int i=0; i<beamVec.size(); i++)
        {
            for(int n=0;n<dim;n++)
            {                
                fitWeight[n]            = 1.E0;
                fitX[n]                 = n * inputParameter.ringRun->bunchInfoPrintInterval;
                fitFFTX[n]              = log(coupledBunchModeAmpX[n+indexStart][i]);
                fitFFTY[n]              = log(coupledBunchModeAmpY[n+indexStart][i]);
                fitFFTZ[n]              = log(coupledBunchModeAmpZ[n+indexStart][i]);
                fitAmpIQX[n]            = log(              ampXIQ[n+indexStart][i]);
                fitAmpIQY[n]            = log(              ampYIQ[n+indexStart][i]);
                fitAmpIQZ[n]            = log(              ampZIQ[n+indexStart][i]);
                fitHilbertCBMX[n]       = log(hilbertCoupledBunchModeAmpX[n+indexStart][i]);
                fitHilbertCBMY[n]       = log(hilbertCoupledBunchModeAmpY[n+indexStart][i]);
                fitHilbertCBMZ[n]       = log(hilbertCoupledBunchModeAmpZ[n+indexStart][i]);

                fitAbsAverXAna[n]       = log(abs(          anaSignalAverX[n+indexStart][i]));
                fitAbsAverYAna[n]       = log(abs(          anaSignalAverY[n+indexStart][i]));
                fitAbsAverZAna[n]       = log(abs(          anaSignalAverZ[n+indexStart][i]));
                

            }
                     
            resFitFFTX        = fittingGSL.FitALinear(fitX,fitWeight,fitFFTX,dim);
            resFitFFTY        = fittingGSL.FitALinear(fitX,fitWeight,fitFFTY,dim);   
            resFitFFTZ        = fittingGSL.FitALinear(fitX,fitWeight,fitFFTZ,dim); 
            resFitIQX         = fittingGSL.FitALinear(fitX,fitWeight,fitAmpIQX,dim);
            resFitIQY         = fittingGSL.FitALinear(fitX,fitWeight,fitAmpIQY,dim);
            resFitIQZ         = fittingGSL.FitALinear(fitX,fitWeight,fitAmpIQZ,dim);

            resfitHilbertCBMX = fittingGSL.FitALinear(fitX,fitWeight,fitHilbertCBMX,dim);
            resfitHilbertCBMY = fittingGSL.FitALinear(fitX,fitWeight,fitHilbertCBMY,dim);
            resfitHilbertCBMZ = fittingGSL.FitALinear(fitX,fitWeight,fitHilbertCBMZ,dim);

            resFitAbsAverXAna = fittingGSL.FitALinear(fitX,fitWeight,fitAbsAverXAna,dim);
            resFitAbsAverYAna = fittingGSL.FitALinear(fitX,fitWeight,fitAbsAverYAna,dim);
            resFitAbsAverZAna = fittingGSL.FitALinear(fitX,fitWeight,fitAbsAverZAna,dim);


            fout<<"! page number "<<i<<endl;
            fout<<setw(15)<<left<<dim<<endl;           
            fout<<setw(15)<<left<<i<<endl;
            fout<<setw(15)<<left<<resFitFFTX[1] / inputParameter.ringParBasic->t0<<endl;
            fout<<setw(15)<<left<<resFitFFTY[1] / inputParameter.ringParBasic->t0<<endl; 
            fout<<setw(15)<<left<<resFitFFTZ[1] / inputParameter.ringParBasic->t0<<endl;
            fout<<setw(15)<<left<<resFitIQX[1]  / inputParameter.ringParBasic->t0<<endl;
            fout<<setw(15)<<left<<resFitIQY[1]  / inputParameter.ringParBasic->t0<<endl; 
            fout<<setw(15)<<left<<resFitIQZ[1]  / inputParameter.ringParBasic->t0<<endl;
            fout<<setw(15)<<left<<resfitHilbertCBMX[1] / inputParameter.ringParBasic->t0<<endl;
            fout<<setw(15)<<left<<resfitHilbertCBMY[1] / inputParameter.ringParBasic->t0<<endl; 
            fout<<setw(15)<<left<<resfitHilbertCBMZ[1] / inputParameter.ringParBasic->t0<<endl; 

            fout<<setw(15)<<left<<beamVec[i].bunchHarmNum<<endl;
            fout<<setw(15)<<left<<resFitAbsAverXAna[1] / inputParameter.ringParBasic->t0<<endl;
            fout<<setw(15)<<left<<resFitAbsAverYAna[1] / inputParameter.ringParBasic->t0<<endl; 
            fout<<setw(15)<<left<<resFitAbsAverZAna[1] / inputParameter.ringParBasic->t0<<endl;


            fout<<ampXIQ.size()<<endl;

            for(int n=0;n<ampXIQ.size();n++)
            {
                fout<<setw(15)<<left<<n * inputParameter.ringRun->bunchInfoPrintInterval
                    <<setw(15)<<left<<coupledBunchModeAmpX[n][i]
                    <<setw(15)<<left<<coupledBunchModeAmpY[n][i]
                    <<setw(15)<<left<<coupledBunchModeAmpZ[n][i]
                    <<setw(15)<<left<<coupledBunchModeArgX[n][i]
                    <<setw(15)<<left<<coupledBunchModeArgY[n][i]
                    <<setw(15)<<left<<coupledBunchModeArgZ[n][i]
                    <<setw(15)<<left<<ampXIQ[n][i]
                    <<setw(15)<<left<<ampYIQ[n][i]
                    <<setw(15)<<left<<ampZIQ[n][i]
                    <<setw(15)<<left<<argXIQ[n][i]
                    <<setw(15)<<left<<argXIQ[n][i]
                    <<setw(15)<<left<<argXIQ[n][i]
                    <<setw(15)<<left<<hilbertCoupledBunchModeAmpX[n][i]
                    <<setw(15)<<left<<hilbertCoupledBunchModeAmpX[n][i]
                    <<setw(15)<<left<<hilbertCoupledBunchModeAmpZ[n][i]
                    <<setw(15)<<left<<hilbertCoupledBunchModeArgX[n][i]
                    <<setw(15)<<left<<hilbertCoupledBunchModeArgX[n][i]
                    <<setw(15)<<left<<hilbertCoupledBunchModeArgZ[n][i]
                    <<setw(15)<<left<<historyAverX[n][i]
                    <<setw(15)<<left<<historyAverY[n][i]
                    <<setw(15)<<left<<historyAverZ[n][i]
                    <<setw(15)<<left<<anaSignalAverX[n][i].real()
                    <<setw(15)<<left<<anaSignalAverY[n][i].real()
                    <<setw(15)<<left<<anaSignalAverZ[n][i].real()
                    <<setw(15)<<left<<abs(anaSignalAverX[n][i])
                    <<setw(15)<<left<<abs(anaSignalAverY[n][i])
                    <<setw(15)<<left<<abs(anaSignalAverZ[n][i])
                    <<endl;
            }    
            
        }

        fout.close();
    }

    gsl_fft_complex_workspace_free (workspace);
    gsl_fft_complex_wavetable_free (wavetable);
}

void SPBeam::GetAnalyticalWithFilter(const ReadInputSettings &inputParameter)
{
    int turns = historyAverX.size();
    double workQx = inputParameter.ringParBasic->workQx;
    double workQy = inputParameter.ringParBasic->workQy;
    double workQz = inputParameter.ringParBasic->workQz;
    int harmonics = inputParameter.ringParBasic->harmonics;
    vector<double> xSignal(turns);
    vector<double> ySignal(turns);
    vector<double> zSignal(turns);

    anaSignalAverX.resize(turns);
    anaSignalAverY.resize(turns);
    anaSignalAverZ.resize(turns);
    
    for(int n=0;n<turns;n++)
    {
        anaSignalAverX[n].resize(beamVec.size());
        anaSignalAverY[n].resize(beamVec.size());
        anaSignalAverZ[n].resize(beamVec.size());
    }


    for(int i=0;i<beamVec.size();i++)
    {
        for(int n=0;n<turns;n++)
        {
            xSignal[n] = historyAverX[n][i];
            ySignal[n] = historyAverY[n][i];
            zSignal[n] = historyAverZ[n][i];
        }
        vector<complex<double> > xAnalytical = GetHilbertAnalytical(xSignal,0.01,workQx);
        vector<complex<double> > yAnalytical = GetHilbertAnalytical(ySignal,0.01,workQy);
        vector<complex<double> > zAnalytical = GetHilbertAnalytical(zSignal,0.01,workQz);  

        for(int n=0;n<turns;n++)
        {
            anaSignalAverX[n][i] = xAnalytical[n];
            anaSignalAverY[n][i] = yAnalytical[n];
            anaSignalAverZ[n][i] = zAnalytical[n];
        }
    }

}


void SPBeam::SPBeamDataPrintPerTurn(int nTurns, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    // get the spectrum of the center motion...

    double alphax = latticeInterActionPoint.twissAlphaX[0];
    double betax  = latticeInterActionPoint.twissBetaX [0];
    double alphay = latticeInterActionPoint.twissAlphaY[0];
    double betay  = latticeInterActionPoint.twissBetaY[0] ;
    double alphaz = 0;
    double betaz  = inputParameter.ringBunchPara->rmsBunchLength / inputParameter.ringBunchPara->rmsEnergySpread;
    double workQx = inputParameter.ringParBasic->workQx;
    double workQy = inputParameter.ringParBasic->workQy;
    double workQz = inputParameter.ringParBasic->workQz;
    int harmonics = inputParameter.ringParBasic->harmonics;
    double tRF    = inputParameter.ringParBasic->t0 / harmonics; 
     
    string filePrefix = inputParameter.ringRun->TBTBunchAverData;
    string fname = filePrefix + ".sdds";
    ofstream fout(fname,ios_base::app);
    
    if(nTurns == 0)
	{
	    fout<<"SDDS1"<<endl;
	    fout<<"&column name=Turns,                                 type=long,   &end"<<endl;
        fout<<"&column name=BunchIndex,                            type=long,   &end"<<endl;
	    fout<<"&column name=HarmIndex,                             type=long,   &end"<<endl;
        fout<<"&column name=BunchCurrent,               units=mA   type=float,  &end"<<endl;
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
            colname = string("&column name=") + string("beamIndAmp_")    + to_string(j) + string(", units=V, type=float,  &end");
            fout<<colname<<endl;
            colname = string("&column name=") + string("beamIndPhase_")  + to_string(j) + string(", units=rad, type=float,  &end");
            fout<<colname<<endl;

            // 2023-02-15  get the cavity Voltage Sampled m*tRf used in feedback
            colname = string("&column name=") + string("cavAmpSamp_")    + to_string(j) + string(", units=V,   type=float,  &end");
            fout<<colname<<endl; 
            colname = string("&column name=") + string("cavPhaseSamp_")  + to_string(j) + string(", units=rad, type=float,  &end");
            fout<<colname<<endl;
            colname = string("&column name=") + string("genAmpSamp_")    + to_string(j) + string(", units=V,   type=float,  &end");
            fout<<colname<<endl; 
            colname = string("&column name=") + string("genPhaseSamp_")  + to_string(j) + string(", units=rad, type=float,  &end");
            fout<<colname<<endl;
            colname = string("&column name=") + string("beamIndAmpSamp_")    + to_string(j) + string(", units=V,   type=float,  &end");
            fout<<colname<<endl; 
            colname = string("&column name=") + string("beamIndPhaseSamp_")  + to_string(j) + string(", units=rad, type=float,  &end");
            fout<<colname<<endl;
            colname = string("&column name=") + string("deltaCavAmpSamp_")    + to_string(j) + string(", units=V,   type=float,  &end");
            fout<<colname<<endl; 
            colname = string("&column name=") + string("deltaCavPhaseSamp_")  + to_string(j) + string(", units=rad, type=float,  &end");
            fout<<colname<<endl;
            colname = string("&column name=") + string("genIgAmp_")      + to_string(j) + string(", units=ampere, type=float,  &end");
            fout<<colname<<endl;
            colname = string("&column name=") + string("genIgPhase_")    + to_string(j) + string(", units=rad, type=float,  &end");
            fout<<colname<<endl;
        }

	    fout<<"&data mode=ascii, &end"<<endl;
	}

    fout<<"! page number "<<nTurns + 1<<endl;
    fout<<beamVec.size()<<endl;

    for(int i=0;i<beamVec.size();i++)
    {
        fout<<setw(24)<<left<<setprecision(16)<<nTurns
            <<setw(24)<<left<<setprecision(16)<<i
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].bunchHarmNum
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].current
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].xAver
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].yAver
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].zAver
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].pxAver
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].pyAver
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].pzAver
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].actionJx
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].actionJy
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].totIonCharge
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].eFxDueToIon[0]
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].eFyDueToIon[0]
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].lRWakeForceAver[0]
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].lRWakeForceAver[1]
            <<setw(24)<<left<<setprecision(16)<<beamVec[i].lRWakeForceAver[2];

            for(int j=0; j<inputParameter.ringParRf->resNum;j++)
            {
                double absCavVolBunchCen =abs(beamVec[i].bunchRFModeInfo->cavVolBunchCen[j]);
                double argCavVolBunchCen =arg(beamVec[i].bunchRFModeInfo->cavVolBunchCen[j]);
                if(absCavVolBunchCen==0) argCavVolBunchCen=0;                
                double absGenVolBunchAver =abs(beamVec[i].bunchRFModeInfo->genVolBunchAver[j]);
                double argGenVolBunchAver =arg(beamVec[i].bunchRFModeInfo->genVolBunchAver[j]);  
                if(absGenVolBunchAver==0) argGenVolBunchAver=0;
                double absInduceVolBunchCen = abs(beamVec[i].bunchRFModeInfo->induceVolBunchCen[j]);
                double argInduceVolBunchCen = arg(beamVec[i].bunchRFModeInfo->induceVolBunchCen[j]);
                if(absInduceVolBunchCen==0) argInduceVolBunchCen=0;

                double argCavPhaseFB = arg(cavityResonator.resonatorVec[j].vCavSample[beamVec[i].bunchHarmNum]);
                double absCavAmpFB   = abs(cavityResonator.resonatorVec[j].vCavSample[beamVec[i].bunchHarmNum]);
                double argGenPhaseFB = arg(cavityResonator.resonatorVec[j].vGenSample[beamVec[i].bunchHarmNum]);
                double absGenPhaseFB = abs(cavityResonator.resonatorVec[j].vGenSample[beamVec[i].bunchHarmNum]);
                double argVBPhaseFB  = arg(cavityResonator.resonatorVec[j].vBSample[beamVec[i].bunchHarmNum]);
                double absVBPhaseFB  = abs(cavityResonator.resonatorVec[j].vBSample[beamVec[i].bunchHarmNum]); 
                double deltaCavAmp   = abs(cavityResonator.resonatorVec[j].deltaVCavSample[beamVec[i].bunchHarmNum]);
                double deltaCavPhase = arg(cavityResonator.resonatorVec[j].deltaVCavSample[beamVec[i].bunchHarmNum]);

                double absIg             = abs(beamVec[i].bunchRFModeInfo->genIgBunchAver[j]);
                double argIg             = arg(beamVec[i].bunchRFModeInfo->genIgBunchAver[j]);
                if(absIg==0)  argIg=0;



                fout<<setw(24)<<left<<setprecision(16)<<absCavVolBunchCen
                    <<setw(24)<<left<<setprecision(16)<<argCavVolBunchCen
                    <<setw(24)<<left<<setprecision(16)<<absGenVolBunchAver
                    <<setw(24)<<left<<setprecision(16)<<argGenVolBunchAver
                    <<setw(24)<<left<<setprecision(16)<<absInduceVolBunchCen
                    <<setw(24)<<left<<setprecision(16)<<argInduceVolBunchCen
                    <<setw(24)<<left<<setprecision(16)<<absCavAmpFB
                    <<setw(24)<<left<<setprecision(16)<<argCavPhaseFB
                    <<setw(24)<<left<<setprecision(16)<<absGenPhaseFB
                    <<setw(24)<<left<<setprecision(16)<<argGenPhaseFB
                    <<setw(24)<<left<<setprecision(16)<<absVBPhaseFB
                    <<setw(24)<<left<<setprecision(16)<<argVBPhaseFB
                    <<setw(24)<<left<<setprecision(16)<<deltaCavAmp
                    <<setw(24)<<left<<setprecision(16)<<deltaCavPhase
                    <<setw(24)<<left<<setprecision(16)<<absIg
                    <<setw(24)<<left<<setprecision(16)<<argIg;
            }
            fout<<endl;
    }

    fout.close();
    
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
         
    // vector<double> resAmpFBRatioForTotSelfLoss = inputParameter.ringParRf->resAmpFBRatioForTotSelfLoss;
    vector<complex<double> > vbKickAver;
    vector<complex<double> > vbAccumAver;
    vector<complex<double> > selfLossToCompensate;
    vector<complex<double> > genVolAver;
    complex<double> totSelfLoss;

    vbKickAver.resize(inputParameter.ringParRf->resNum);
    vbAccumAver.resize(inputParameter.ringParRf->resNum);
    genVolAver.resize(inputParameter.ringParRf->resNum);
    selfLossToCompensate.resize(inputParameter.ringParRf->resNum);

    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double circRing   = inputParameter.ringParBasic->circRing;
    double rfLen      = circRing  / ringHarmH;
    double tRF        = inputParameter.ringParBasic->t0 / ringHarmH; 
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double t0         = inputParameter.ringParBasic->t0;

    
    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        complex<double> vbKickAverSum=(0.e0,0.e0);
        complex<double> vbAccumAverSum=(0.e0,0.e0);
        complex<double> genVolAverSum=(0.e0,0.e0);
        for(int j=0;j<beamVec.size();j++)
        {
            vbKickAverSum      += beamVec[j].bunchRFModeInfo->selfLossVolBunchCen[i];
            vbAccumAverSum     += beamVec[j].bunchRFModeInfo->induceVolBunchCen[i];
            genVolAverSum      += beamVec[j].bunchRFModeInfo->genVolBunchAver[i];
        }
        vbKickAver[i]      = vbKickAverSum   / double(beamVec.size());
        vbAccumAver[i]     = vbAccumAverSum  / double(beamVec.size());
        genVolAver[i]      = genVolAverSum   / double(beamVec.size());
    }



    

    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        //each cavity generator compensated the related self-loss by it self
        //double genAddvbArg = arg(cavityResonator.resonatorVec[i].resCavVolReq);
        //double genAddvbAbs = (cavityResonator.resonatorVec[i].resCavVolReq.real() - selfLossToCompensate[i].real()) / cos( genAddvbArg );
        //cavityResonator.resonatorVec[i].resGenVol = genAddvbAbs * exp( li * genAddvbArg ) - vbAccumAver[i];
        
        // active cavity, and cavity feedback is set as 1 in input file, HeFei's method.
        // if(cavityResonator.resonatorVec[i].rfResCavVolFB==1  && cavityResonator.resonatorVec[i].resType==1)
        // {    
        //     cavityResonator.resonatorVec[i].resGenVol = cavityResonator.resonatorVec[i].resCavVolReq - vbAccumAver[i];
        // }
        // else if(cavityResonator.resonatorVec[i].rfResCavVolFB==0  && cavityResonator.resonatorVec[i].resType==1)
        // {
        //     cavityResonator.resonatorVec[i].resGenVol = cavityResonator.resonatorVec[i].resGenVol;
        // }
        // else if(cavityResonator.resonatorVec[i].rfResCavVolFB==1  && cavityResonator.resonatorVec[i].resType==0)
        // {
        //     cerr<<"Wrong setting: passive cavity can set cavitiy FB"<<endl;
        //     exit(0);
        // }
        // else if (cavityResonator.resonatorVec[i].rfResCavVolFB==0  && cavityResonator.resonatorVec[i].resType==0)
        // {
        //     cavityResonator.resonatorVec[i].resGenVol = complex<double>(0.E0,0.E0);
        // }
                
        // set voltage and phase for hainssinki solver
        for(int j=0;j<beamVec.size();j++)
        {              
            beamVec[j].haissinski->cavAmp[i]   = abs(beamVec[j].bunchRFModeInfo->cavVolBunchCen[i]);
            beamVec[j].haissinski->cavPhase[i] = arg(beamVec[j].bunchRFModeInfo->cavVolBunchCen[i]);            
        }

    }   


    //longitudinal tracking ---
    if(beamVec.size()==1)   
    {
        //get the time distance from current bunch to next bunch - in SP model 
        beamVec[0].timeFromCurrnetBunchToNextBunch  = beamVec[0].bunchGap * tRF + (beamVec[0].zAverLastTurn - beamVec[0].zAver ) / CLight / rBeta;
        beamVec[0].zAverLastTurn                    = beamVec[0].zAver; 
        // beamVec[i].BunchTransferDueToLatticeLMatarix(inputParameter,cavityResonator);            
        // beamVec[0].BunchTransferDueToLatticeLYamamoto(inputParameter,cavityResonator);
        beamVec[0].BunchTransferDueToLatticeL(inputParameter,cavityResonator); 
        beamVec[0].GetSPBunchRMS(latticeInterActionPoint, 0);                 
    }
    else    // when multi-bunches in tracking
    {
        // save the last turn position info    
        for(int i=0;i<beamVec.size();i++)
        {
            beamVec[i].zAverLastTurn = beamVec[i].zAver;
        }
        
        for(int i=0;i<beamVec.size();i++)
        {     
            // get the time for next bunch move to RF interaction point:  Trf + deltaT[i+1] - deltaT[i] 
            if(i<beamVec.size()-1)
            {
                beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zAver - beamVec[i+1].zAver)/ CLight / rBeta; 
            }
            else
            {
                beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zAver - beamVec[0].zAver ) / CLight / rBeta;
            }
            // beamVec[i].BunchTransferDueToLatticeLMatarix(inputParameter,cavityResonator);            
            // beamVec[i].BunchTransferDueToLatticeLYamamoto(inputParameter,cavityResonator);
             beamVec[i].BunchTransferDueToLatticeL(inputParameter,cavityResonator);      
            beamVec[i].GetSPBunchRMS(latticeInterActionPoint,0);          
        }
    }       
}

void SPBeam::BeamTransferPerTurnDueToLatticeLTest(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,CavityResonator &cavityResonator, int turns)
{ 
         
    // vector<double> resAmpFBRatioForTotSelfLoss = inputParameter.ringParRf->resAmpFBRatioForTotSelfLoss;
    // vector<complex<double> > vbKickAver;
    // vector<complex<double> > vbAccumAver;
    // vector<complex<double> > selfLossToCompensate;
    // vector<complex<double> > genVolAver;
    // complex<double> totSelfLoss;

    // vbKickAver.resize(inputParameter.ringParRf->resNum);
    // vbAccumAver.resize(inputParameter.ringParRf->resNum);
    // genVolAver.resize(inputParameter.ringParRf->resNum);
    // selfLossToCompensate.resize(inputParameter.ringParRf->resNum);

    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double circRing   = inputParameter.ringParBasic->circRing;
    double rfLen      = circRing  / ringHarmH;
    double tRF        = inputParameter.ringParBasic->t0 / ringHarmH; 
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double t0         = inputParameter.ringParBasic->t0;

    
    // for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    // {
    //     complex<double> vbKickAverSum=(0.e0,0.e0);
    //     complex<double> vbAccumAverSum=(0.e0,0.e0);
    //     complex<double> genVolAverSum=(0.e0,0.e0);
    //     for(int j=0;j<beamVec.size();j++)
    //     {
    //         vbKickAverSum      += beamVec[j].bunchRFModeInfo->selfLossVolBunchCen[i];
    //         vbAccumAverSum     += beamVec[j].bunchRFModeInfo->induceVolBunchCen[i];
    //         genVolAverSum      += beamVec[j].bunchRFModeInfo->genVolBunchAver[i];
    //     }
    //     vbKickAver[i]      = vbKickAverSum   / double(beamVec.size());
    //     vbAccumAver[i]     = vbAccumAverSum  / double(beamVec.size());
    //     genVolAver[i]      = genVolAverSum   / double(beamVec.size());
    // }



    

    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        //each cavity generator compensated the related self-loss by it self
        //double genAddvbArg = arg(cavityResonator.resonatorVec[i].resCavVolReq);
        //double genAddvbAbs = (cavityResonator.resonatorVec[i].resCavVolReq.real() - selfLossToCompensate[i].real()) / cos( genAddvbArg );
        //cavityResonator.resonatorVec[i].resGenVol = genAddvbAbs * exp( li * genAddvbArg ) - vbAccumAver[i];
        
        // active cavity, and cavity feedback is set as 1 in input file, HeFei's method.
        // if(cavityResonator.resonatorVec[i].rfResCavVolFB==1  && cavityResonator.resonatorVec[i].resType==1)
        // {    
        //     cavityResonator.resonatorVec[i].resGenVol = cavityResonator.resonatorVec[i].resCavVolReq - vbAccumAver[i];
        // }
        // else if(cavityResonator.resonatorVec[i].rfResCavVolFB==0  && cavityResonator.resonatorVec[i].resType==1)
        // {
        //     cavityResonator.resonatorVec[i].resGenVol = cavityResonator.resonatorVec[i].resGenVol;
        // }
        // else if(cavityResonator.resonatorVec[i].rfResCavVolFB==1  && cavityResonator.resonatorVec[i].resType==0)
        // {
        //     cerr<<"Wrong setting: passive cavity can set cavitiy FB"<<endl;
        //     exit(0);
        // }
        // else if (cavityResonator.resonatorVec[i].rfResCavVolFB==0  && cavityResonator.resonatorVec[i].resType==0)
        // {
        //     cavityResonator.resonatorVec[i].resGenVol = complex<double>(0.E0,0.E0);
        // }
                
        // set voltage and phase for hainssinki solver
        for(int j=0;j<beamVec.size();j++)
        {              
            beamVec[j].haissinski->cavAmp[i]   = abs(beamVec[j].bunchRFModeInfo->cavVolBunchCen[i]);
            beamVec[j].haissinski->cavPhase[i] = arg(beamVec[j].bunchRFModeInfo->cavVolBunchCen[i]);            
        }

    }   


    //longitudinal tracking ---
    if(beamVec.size()==1)   
    {
        //get the time distance from current bunch to next bunch - in SP model 
        beamVec[0].timeFromCurrnetBunchToNextBunch  = beamVec[0].bunchGap * tRF + (beamVec[0].zAverLastTurn - beamVec[0].zAver ) / CLight / rBeta;
        beamVec[0].zAverLastTurn                    = beamVec[0].zAver; 
        // beamVec[i].BunchTransferDueToLatticeLMatarix(inputParameter,cavityResonator);            
        // beamVec[0].BunchTransferDueToLatticeLYamamoto(inputParameter,cavityResonator);
        beamVec[0].BunchTransferDueToLatticeL(inputParameter,cavityResonator); 
        beamVec[0].GetSPBunchRMS(latticeInterActionPoint, 0);                 
    }
    else    // when multi-bunches in tracking
    {
        // save the last turn position info    
        for(int i=0;i<beamVec.size();i++)
        {
            beamVec[i].zAverLastTurn = beamVec[i].zAver;
        }
        
        for(int i=0;i<beamVec.size();i++)
        {     
            // get the time for next bunch move to RF interaction point:  Trf + deltaT[i+1] - deltaT[i] 
            if(i<beamVec.size()-1)
            {
                beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zAver - beamVec[i+1].zAver)/ CLight / rBeta; 
            }
            else
            {
                beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zAver - beamVec[0].zAver ) / CLight / rBeta;
            }
            // beamVec[i].BunchTransferDueToLatticeLMatarix(inputParameter,cavityResonator);            
            // beamVec[i].BunchTransferDueToLatticeLYamamoto(inputParameter,cavityResonator);
            beamVec[i].BunchTransferDueToLatticeL(inputParameter,cavityResonator);      
            beamVec[i].GetLongiKickDueToCavFB(inputParameter,cavityResonator);            
            beamVec[i].GetSPBunchRMS(latticeInterActionPoint,0);          

        }
    }       
}

void SPBeam::GetTimeDisToNextBunch(ReadInputSettings &inputParameter)
{
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double circRing   = inputParameter.ringParBasic->circRing;
    double rfLen      = circRing  / ringHarmH;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double tRF        = inputParameter.ringParBasic->t0 / ringHarmH; 

    double temp=0;
    // get the bunch distance

    if(beamVec.size()>1)
    {
        for(int i=0;i<beamVec.size();i++)
        {
            // get the time from current bunch to next bunch
            if(i<beamVec.size()-1)
            {
                beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zAver - beamVec[i+1].zAver) / CLight / rBeta; 
            }
            else
            {
                beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zAver - beamVec[0  ].zAver ) / CLight / rBeta;
            }
   
            // // get the time from last bunch to current bunch
            // if(i==0)
            // {
            //     beamVec[i].timeFromLastBunchToCurrentBunch  = beamVec[beamVec.size()-1].bunchGap * tRF + (beamVec[beamVec.size()-1].zAver - beamVec[0].zAver ) / CLight / rBeta; 
            // }
            // else
            // {
            //     beamVec[i].timeFromLastBunchToCurrentBunch  = beamVec[i-1             ].bunchGap * tRF + (beamVec[i             -1].zAver - beamVec[i].zAver ) / CLight / rBeta;
            // }
        }
    }
    else
    {
        // if one bunch the inital bunch distance is set as ring circumference.
        beamVec[0].timeFromCurrnetBunchToNextBunch   = beamVec[0].bunchGap * tRF;        
        beamVec[0].timeFromLastBunchToCurrentBunch   = beamVec[0].bunchGap * tRF; 
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


void SPBeam::FIRBunchByBunchFeedback(const ReadInputSettings &inputParameter,FIRFeedBack &firFeedBack,int nTurns)
{
    //y[0] = \sum_0^{N} a_k x[-k]
    //Ref. Nakamura's paper spring 8 notation here used is the same with Nakamura's paper
    double rBeta = inputParameter.ringParBasic->rBeta;
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
            beamVec[i].eMomentumZ[j] = beamVec[i].eMomentumZ[j] + energyKickU[i] / pow(rBeta,2);
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
    double deltaTij=0;
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
            //  Tij =  nij * Trf + (deltaT_j - deltaT_i)
            //      =  nij * Trf - (deltaZ_j - deltaZ_i) / c 
            //      =  nij * Trf - (deltaZ_j - deltaZ_i) / c
            
            // finnally  -Tij  is applied in the wake-force subroutine.
            // -Tij =  - nij * Trf + (deltaZ_j - deltaZ_i) / c

            for(int i=tempIndex0;i<=tempIndex1;i++)
            {                 
                nTauij   = beamVec[i].bunchHarmNum - beamVec[j].bunchHarmNum - n * harmonics;
                tauijStastic = nTauij * tRF;
                
                deltaTij = (beamVec[j].zAver -  wakefunction.poszData[nTurnswakeTrunction-1-n][i]) / CLight / rBeta;
                tauij    = tauijStastic  + deltaTij;   
                
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
                    beamVec[j].lRWakeForceAver[2] -= wakeForceTemp[2] * beamVec[i].electronNumPerBunch; 
                    // to subsctract the "stastic term -- only when coupled bunch modes"
                    // wakeStasticForceTemp = wakefunction.GetBBRWakeFun(tauijStastic);
                    // beamVec[j].lRWakeForceAver[2] -= (wakeForceTemp[2] - wakeStasticForceTemp[2]) * beamVec[i].electronNumPerBunch;                                                        //[V/C]         ->  [V/C]                
                }
            }   
        }
        

        beamVec[j].lRWakeForceAver[0] *=  ElectronCharge / electronBeamEnergy / pow(rBeta,2);   // [V/C] * [C] * [1e] / [eV] ->rad
        beamVec[j].lRWakeForceAver[1] *=  ElectronCharge / electronBeamEnergy / pow(rBeta,2);   // [V/C] * [C] * [1e] / [eV] ->rad
        beamVec[j].lRWakeForceAver[2] *=  ElectronCharge / electronBeamEnergy / pow(rBeta,2);   // [V/C] * [C] * [1e] / [eV] ->rad
    }

    BeamTransferPerTurnDueWake(); 

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