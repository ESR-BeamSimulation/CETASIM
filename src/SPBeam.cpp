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
        beamVec[i].InitialSPBunch(inputParameter);
        beamVec[i].DistriGenerator(latticeInterActionPoint,inputParameter,i);
        if (i==0)
        {
            beamVec[0].ePositionX[0] = 1.E-4;
            beamVec[0].ePositionY[0] = 1.E-4;
            beamVec[0].ePositionZ[0] = 0.E0;
            beamVec[0].eMomentumX[0] = 0.E0;
            beamVec[0].eMomentumY[0] = 0.E0;
            beamVec[0].eMomentumZ[0] = 0.E0;
        }
    }

    SPBeamRMSCal(latticeInterActionPoint, 0);
    
    for(int i=0;i<totBunchNum-1;i++)
    {
        beamVec[i].bunchGap = beamVec[i+1].bunchHarmNum - beamVec[i].bunchHarmNum;
    }
    beamVec[totBunchNum-1].bunchGap = harmonics - beamVec[totBunchNum-1].bunchHarmNum;

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
            
            /* do the average over the whole rf bucket.
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
            */
            // do the averarage over the occupied bucket 

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
                tB = beamVec[k].bunchGap * 1.0 / f0 / ringHarmH ;           //[s]
                deltaL = tB / tF ;
                cPsi   = deltaL * tan(cavityResonator.resonatorVec[i].resDeTunePsi );

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
        vbAccumAver =  vbAccumAver2/double(beamVec.size());                          // after decay and rotate, just before next bunch comes in.

        cavityResonator.resonatorVec[i].vbAccum   =  vbAccumAver;               // set cavity beam induced voltage
        vbKickAver  =  vbKickAccum /double(beamVec.size());                     // average energy that bunch is kicked by self-loss Vb0/2. 

        // extra feedback condition to compensate the self-loss term vb0/2          
        //double genAddvbArg = arg(cavityResonator.resonatorVec[i].resCavVolReq);
        //if(genAddvbArg==PI/2.0 or genAddvbArg==-PI/2.0)
        //{
        //    cerr<<"Required Cavity Phase is PI/2 or -PI/2"<<endl;
        //    cerr<<"Does not work when cavity feedback is included"<<endl;
        //    exit(0);
        //}
        //double genAddvbAbs = (cavityResonator.resonatorVec[i].resCavVolReq.real() - vbKickAver.real()) / cos( genAddvbArg );
        //genAddvbAbs = abs(genAddvbAbs);
        //cavityResonator.resonatorVec[i].resGenVol = genAddvbAbs * exp(li * genAddvbArg ) - vbAccumAver;
        // end of the extra feedback 

        // simply condition. resGenVol can not compensate the self-loss term, bunch will have a center shift due to self-loss.  
        cavityResonator.resonatorVec[i].resGenVol = cavityResonator.resonatorVec[i].resCavVolReq - vbAccumAver;

        //cout<<abs(cavityResonator.resonatorVec[i].resCavVolReq)<<" "<<arg(cavityResonator.resonatorVec[i].resCavVolReq)<<endl;
        //cout<<abs(cavityResonator.resonatorVec[i].resGenVol)<<" "<<arg(cavityResonator.resonatorVec[i].resGenVol)<<endl;

        // DC solution -- which is exactly the same as above equaiton (1)'s reuslts, very good agreement.
        //vbAccumAver = 2.0 * beamCurr * cavityResonator.resonatorVec[i].resShuntImpRs / (1.0 + cavityResonator.resonatorVec[i].resCouplingBeta)
        //                             * cos(cavityResonator.resonatorVec[i].resDeTunePsi) * exp(li * (PI + cavityResonator.resonatorVec[i].resDeTunePsi)  );


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
    fout<<"&parameter name=CavAmpIdeal,   units=V,   type=float,  &end"<<endl;
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

void SPBeam::Run(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter, CavityResonator &cavityResonator)
{
    int nTurns                      = inputParameter.ringRun->nTurns;
    int synRadDampingFlag           = inputParameter.ringRun->synRadDampingFlag;
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
        
        vector<double> posxDataTemp(beamVec.size());
        vector<double> posyDataTemp(beamVec.size());
        vector<double> poszDataTemp(beamVec.size());
        
        for(int i=0;i<beamVec.size();i++)
        {
            posxDataTemp[i] = beamVec[i].xAver;
            posyDataTemp[i] = beamVec[i].yAver;
            poszDataTemp[i] = beamVec[i].zAver;
        }

        // set positions of bunches at the first turn 
        lRWakeFunction.posxData.push_back(posxDataTemp);
        lRWakeFunction.posyData.push_back(posyDataTemp);
        lRWakeFunction.poszData.push_back(poszDataTemp);

        lRWakeFunction.posxData.erase(lRWakeFunction.posxData.begin());
        lRWakeFunction.posyData.erase(lRWakeFunction.posyData.begin());
        lRWakeFunction.poszData.erase(lRWakeFunction.poszData.begin());

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
        string colname = string("&column name=") + string("x_")     + to_string(j) + string(", units=m,   type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("px_")            + to_string(j) + string(", units=rad, type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("y_")             + to_string(j) + string(", units=m,   type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("py_")            + to_string(j) + string(", units=rad, type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("z_")             + to_string(j) + string(", units=m,   type=float,  &end");
        fout<<colname<<endl;
        colname = string("&column name=") + string("pz_")            + to_string(j) + string(", units=rad, type=float,  &end");
        fout<<colname<<endl;
    }
    fout<<"&data mode=ascii, &end"<<endl;
    fout<<"! page number "<<1<<endl;
    fout<<nTurns<<endl;
             
    // run loop starts, for nTrns and each trun for k interaction-points
    for(int n=0;n<nTurns;n++)
    {
        if(n%100==0) cout<<n<<"  turns"<<endl;

        SPBeamRMSCal(latticeInterActionPoint, 0);
        SPGetBeamInfo();
        //SetBeamPosHistoryDataWithinWindow(); -- Maro's approaches to get the coupled bunch growthrate.. PRAB
                
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
        
        BeamTransferPerTurnDueToLatticeL(inputParameter,cavityResonator); 

        if(lRWakeFlag)
        {
            SPBeamRMSCal(latticeInterActionPoint, 0);
            LRWakeBeamIntaction(inputParameter,lRWakeFunction,latticeInterActionPoint);
            BeamTransferPerTurnDueWake();   
        }
        
        if(synRadDampingFlag)
        {
            SPBeamRMSCal(latticeInterActionPoint, 0);
            BeamSynRadDamping(inputParameter,latticeInterActionPoint);
        }
        
        if(fIRBunchByBunchFeedbackFlag)
        {            
            SPBeamRMSCal(latticeInterActionPoint, 0);
            FIRBunchByBunchFeedback(firFeedBack,n);
        }

        SPBeamRMSCal(latticeInterActionPoint, 0);
        SPGetBeamInfo(); 

        if(  bunchInfoPrintInterval && (n%bunchInfoPrintInterval==0)  )
        {             
            SPBeamRMSCal(latticeInterActionPoint, 0);
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
    }
    fout.close();
    cout<<"End of Tracking "<<nTurns<< "Turns"<<endl;
       
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
        tempAverX.push_back(beamVec[i].xAver);
        tempAverY.push_back(beamVec[i].yAver);
    }
    
    auto maxJx    = max_element(tempJx.begin(),tempJx.end() );
    auto maxJy    = max_element(tempJy.begin(),tempJy.end() );
    auto maxAverX = max_element(tempAverX.begin(),tempAverX.end() );
    auto maxAverY = max_element(tempAverY.begin(),tempAverY.end() );
    

    weakStrongBeamInfo->actionJxMax    = *maxJx;
    weakStrongBeamInfo->actionJyMax    = *maxJy;    
    weakStrongBeamInfo->bunchAverXMax  = *maxAverX;
    weakStrongBeamInfo->bunchAverYMax  = *maxAverY;

    weakStrongBeamInfo->bunchEffectiveSizeXMax = beamVec[0].rmsRx + *maxAverX;
    weakStrongBeamInfo->bunchEffectiveSizeYMax = beamVec[0].rmsRy + *maxAverY;

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

void SPBeam::SPBeamDataPrintPerTurn(int nTurns, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter)
{
    // get the spectrum of the center motion...
    
    gsl_fft_complex_wavetable * wavetable;
	gsl_fft_complex_workspace * workspace;
	wavetable = gsl_fft_complex_wavetable_alloc (beamVec.size());
	workspace = gsl_fft_complex_workspace_alloc (beamVec.size());

    double alphax = latticeInterActionPoint.twissAlphaX[0];
    double betax  = latticeInterActionPoint.twissBetaX [0];
    double alphay = latticeInterActionPoint.twissAlphaY[0];
    double betay  = latticeInterActionPoint.twissBetaY[0] ;
    double alphaz = 0;
    double betaz  = 1;
    double workQx = inputParameter.ringParBasic->workQx;
    double workQy = inputParameter.ringParBasic->workQy;
    double workQz = inputParameter.ringParBasic->workQz;


    double tempX[2*beamVec.size()];
    double tempY[2*beamVec.size()];
    double tempZ[2*beamVec.size()];
    double tempXFFT[2*beamVec.size()];
    double tempYFFT[2*beamVec.size()];
    double tempZFFT[2*beamVec.size()];
 
    // IPAC 2022 WEPOMS010 -- Diamond-II
    for(int i=0;i<beamVec.size();i++)
	{
        double x  =  beamVec[i].xAver;
        double y  =  beamVec[i].yAver;
        double z  =  beamVec[i].zAver;
        
        double px =  beamVec[i].pxAver;
        double py =  beamVec[i].pyAver;
        double pz =  beamVec[i].pzAver;


        double phasex = - 2.0 * PI * workQx * i / beamVec.size();
        double phasey = - 2.0 * PI * workQy * i / beamVec.size();
        double phasez = - 2.0 * PI * workQz * i / beamVec.size();
        complex<double> zx = x / sqrt(betax) - li * ( sqrt(betax) * px + alphax / sqrt(betax) * x ) * exp (li * phasex);
        complex<double> zy = y / sqrt(betay) - li * ( sqrt(betay) * py + alphay / sqrt(betay) * y ) * exp (li * phasey); 
        complex<double> zz = z / sqrt(betaz) - li * ( sqrt(betaz) * pz + alphaz / sqrt(betaz) * z ) * exp (li * phasez);  

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
    gsl_fft_complex_workspace_free (workspace);
    gsl_fft_complex_wavetable_free (wavetable);
     
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
            <<setw(24)<<left<<setprecision(16)<<tempX[2*i]
            <<setw(24)<<left<<setprecision(16)<<tempY[2*i]
            <<setw(24)<<left<<setprecision(16)<<tempZ[2*i]
            <<setw(24)<<left<<setprecision(16)<<tempX[2*i+1]
            <<setw(24)<<left<<setprecision(16)<<tempY[2*i+1]
            <<setw(24)<<left<<setprecision(16)<<tempZ[2*i+1]
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
            <<setw(24)<<left<<setprecision(16)<<sqrt( pow(tempZFFT[2*i],2) + pow(tempZFFT[2*i+1],2) );

            for(int j=0; j<inputParameter.ringParRf->resNum;j++)
            {
                fout<<setw(24)<<left<<setprecision(16)<<    beamVec[i].cavFBCenInfo->cavAmpBunchCen[j]
                    <<setw(24)<<left<<setprecision(16)<<    beamVec[i].cavFBCenInfo->cavPhaseBunchCen[j]
                    <<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].cavFBCenInfo->genVolBunchAver[j])
                    <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].cavFBCenInfo->genVolBunchAver[j])
                    <<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].cavFBCenInfo->induceVolBunchCen[j])
                    <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].cavFBCenInfo->induceVolBunchCen[j]);
            }
            fout<<endl;
    }

    fout.close();

    ofstream foutX("coupledBunchedModeX.sdds",ios_base::app);
    ofstream foutY("coupledBunchedModeY.sdds",ios_base::app);
    
    foutX<<setw(15)<<left<<nTurns;
    foutY<<setw(15)<<left<<nTurns;
    for(int i=0;i<beamVec.size();i++)
    {    
        foutX<<setw(15)<<left<< log( sqrt( pow(tempXFFT[2*i],2) + pow(tempXFFT[2*i+1],2) ) );
        foutY<<setw(15)<<left<< log( sqrt( pow(tempYFFT[2*i],2) + pow(tempYFFT[2*i+1],2) ) );
    }
    foutX<<endl;
    foutY<<endl;
    foutX.close();
    foutY.close();


    /*
    double nux =  inputParameter.ringParBasic->workQx;
    vector<double> xHistoryData ;
    vector<double> sinParamFit(2,0);
    for(int j=0;j<beamVec.size();j++)
    {
        xHistoryData =  beamVec[j].xyzHistoryDataToFit[0];
        FittingGSL fittingGSL;
        sinParamFit = fittingGSL.FitASin(xHistoryData,nux);
    }
    */
}

void SPBeam::WSBeamIonEffectOneInteractionPoint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k)
{
    int totBunchNum = beamVec.size();

    for(int j=0;j<totBunchNum;j++)
    {
        latticeInterActionPoint.GetIonNumberPerInterAction(beamVec[j].electronNumPerBunch, k);
        latticeInterActionPoint.IonGenerator(beamVec[j].rmsRx,beamVec[j].rmsRy,beamVec[j].xAver,beamVec[j].yAver,k);
        latticeInterActionPoint.IonsUpdate(k);
        beamVec[j].WSIonBunchInteraction(latticeInterActionPoint,k);
        beamVec[j].BunchTransferDueToIon(latticeInterActionPoint,k);
        latticeInterActionPoint.IonTransferDueToBunch(beamVec[j].bunchGap,k,weakStrongBeamInfo->bunchEffectiveSizeXMax,weakStrongBeamInfo->bunchEffectiveSizeYMax);

    }

}

void SPBeam::BeamTransferPerTurnDueToLatticeL(ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{ 
    GetBinDistBetweenBunch(inputParameter);
    
    vector<double> resAmpFBRatioForTotSelfLoss = inputParameter.ringParRf->resAmpFBRatioForTotSelfLoss;


    vector<complex<double> > vbKickAver;
    vector<complex<double> > vbAccumAver;
    vector<complex<double> > selfLossToCompensate;
    complex<double> totSelfLoss;
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
        // extra feedback condition to compensate the self-loss term vb0/2 
        //double genAddvbArg = arg(cavityResonator.resonatorVec[i].resCavVolReq);
        //double genAddvbAbs = (cavityResonator.resonatorVec[i].resCavVolReq.real() - vbKickAver[i].real()) / cos( genAddvbArg );
        //cavityResonator.resonatorVec[i].resGenVol = genAddvbAbs * exp( li * genAddvbArg ) - vbAccumAver[i];
        
        cavityResonator.resonatorVec[i].resGenVol = cavityResonator.resonatorVec[i].resCavVolReq - vbAccumAver[i];

        for(int j=0;j<beamVec.size();j++)
        {
            beamVec[j].cavFBCenInfo->genVolBunchAver[i] =  cavityResonator.resonatorVec[i].resGenVol;
        }

    }
    
     // set the phase and vol used in haissinski solver
     /*
    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {    
        complex<double> gen = cavityResonator.resonatorVec[i].resCavVolReq - vbAccumAver[i];                    
        for(int j=0;j<beamVec.size();j++)
        {              
            beamVec[j].haissinski->cavAmp[i]   = abs(beamVec[j].cavFBCenInfo->induceVolBunchCen[i] + gen);
            beamVec[j].haissinski->cavPhase[i] = arg(beamVec[j].cavFBCenInfo->induceVolBunchCen[i] + gen);            
            //with this cavity amp and pahse, self-loss are included 
            //beamVec[j].haissinski->cavAmp[i]   = abs(beamVec[j].cavFBCenInfo->cavVolBunchCen[i]);
            //beamVec[j].haissinski->cavPhase[i] = arg(beamVec[j].cavFBCenInfo->cavVolBunchCen[i]);
        }        
    }
    */
    // longitudinal tracking due to rf
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToLatticeL(inputParameter,cavityResonator);
    }



}

void SPBeam::GetBinDistBetweenBunch(ReadInputSettings &inputParameter)
{

    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double circRing   = inputParameter.ringParBasic->circRing ;
    double rfLen      = circRing  / ringHarmH;

    for(int i=0;i<beamVec.size();i++)
    {
        if(beamVec[i].macroEleNumPerBunch==1) 
        {
            if(i<beamVec.size()-1)
            {
                timeBetweenBunch[i] = beamVec[i].bunchGap * rfLen + (beamVec[i].zAver - beamVec[i+1].zAver) ;
            }
            else
            {
                timeBetweenBunch[i] = beamVec[i].bunchGap * rfLen + (beamVec[i].zAver - beamVec[0].zAver) ;
            }
            timeBetweenBunch[i] /= CLight;

            beamVec[i].timeToFirstBinOfNextBunch = timeBetweenBunch[i];
        }
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

void SPBeam::BeamTransferPerTurnDueToLattice(LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter,CavityResonator &cavityResonator)
{
    BeamTransferPerTurnDueToLatticeT(latticeInterActionPoint);
    BeamTransferPerTurnDueToLatticeL(inputParameter,cavityResonator); 
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

void SPBeam::LRWakeBeamIntaction(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction,const  LatticeInterActionPoint &latticeInterActionPoint)
{

    int nTurnswakeTrunction     = inputParameter.ringLRWake->nTurnswakeTrunction;
    int harmonics               = inputParameter.ringParBasic->harmonics;
    double electronBeamEnergy   = inputParameter.ringParBasic->electronBeamEnergy;

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

    int tempIndex0,tempIndex1;
    
    //ofstream fout ("test.dat",ios_base::app);
    for (int j=0;j<beamVec.size();j++)
    {
	    beamVec[j].lRWakeForceAver[0] =0.E0;      // x rad
	    beamVec[j].lRWakeForceAver[1] =0.E0;      // y rad
	    beamVec[j].lRWakeForceAver[2] =0.E0;      // z rad
        
	    for(int n=0;n<nTurnswakeTrunction;n++)
	    {
            // self-interation of bunch in current turn is excluded. 
            if(n==0)
            {
                tempIndex0   = j + 1;
                tempIndex1   = beamVec.size() - 1;
            }
            else if (n==nTurnswakeTrunction-1)
            {
                tempIndex0   = 0;
                tempIndex1   = j;
            }
            else
            {
                tempIndex0   = 0;
                tempIndex1   = beamVec.size()-1;
            }

            for(int i=tempIndex0;i<=tempIndex1;i++)
            {
                
                nTauij = beamVec[j].bunchHarmNum - beamVec[i].bunchHarmNum - n * harmonics;
                tauij  = nTauij * tRF;
	            if(!inputParameter.ringLRWake->pipeGeoInput.empty())
	            {                                       
                    wakeForceTemp = wakefunction.GetRWLRWakeFun(tauij); // V/C/m
                    beamVec[j].lRWakeForceAver[0] -= wakeForceTemp[0] * beamVec[i].electronNumPerBunch * wakefunction.posxData[nTurnswakeTrunction-1-n][i];    //[V/C/m] * [m] ->  [V/C]
                    beamVec[j].lRWakeForceAver[1] -= wakeForceTemp[1] * beamVec[i].electronNumPerBunch * wakefunction.posyData[nTurnswakeTrunction-1-n][i];    //[V/C/m] * [m] ->  [V/C]
                    beamVec[j].lRWakeForceAver[2] -= wakeForceTemp[2] * beamVec[i].electronNumPerBunch ;                                                       //[V/C]         ->  [V/C]
                }

                if(!inputParameter.ringLRWake->bbrInput.empty())
	            {
                    wakeForceTemp = wakefunction.GetBBRWakeFun(tauij);  // V/C/m
	                beamVec[j].lRWakeForceAver[0] -= wakeForceTemp[0] * beamVec[i].electronNumPerBunch * wakefunction.posxData[nTurnswakeTrunction-1-n][i];    //[V/C/m] * [m] ->  [V/C]
	                beamVec[j].lRWakeForceAver[1] -= wakeForceTemp[1] * beamVec[i].electronNumPerBunch * wakefunction.posyData[nTurnswakeTrunction-1-n][i];    //[V/C/m] * [m] ->  [V/C]
	                beamVec[j].lRWakeForceAver[2] -= wakeForceTemp[2] * beamVec[i].electronNumPerBunch;                                                        //[V/C]         ->  [V/C]
                }

                //fout<<setw(15)<<left<<n
                //    <<setw(15)<<left<<beamVec[i].lRWakeForceAver[0]
                //    <<setw(15)<<left<<wakeForceTemp[0]
                //    <<setw(15)<<left<<wakefunction.posxData[nTurnswakeTrunction-1-n][j]
                //    <<setw(15)<<left<<beamVec[j].electronNumPerBunch
                //    <<endl;

                //fout<<setw(15)<<left<<nTurnswakeTrunction-1-n
                //    <<setw(15)<<left<<j
                //    <<setw(15)<<left<<i
                //    <<setw(15)<<left<<tauij*CLight
                //    <<setw(15)<<left<<left<<nTauij
                //    <<setw(15)<<left<<wakeForceTemp[0]
                //    <<setw(15)<<left<<wakeForceTemp[1]
                //    <<setw(15)<<left<<wakeForceTemp[2]
                //    <<endl;

                //cout<<n<<" Turns  i="<<beamVec[i].bunchHarmNum<<" j="<<beamVec[j].bunchHarmNum<<" Trf="<<nTauij<<"  "<<nTurnswakeTrunction-1-n
                //       <<"   "<<wakefunction.posxData[nTurnswakeTrunction-1-n][j]<<"   "<<wakefunction.posxData[nTurnswakeTrunction-1-n][j]<<endl;
            }
                   
        }
        
        
        //fout<<i<<"  "<<beamVec[i].lRWakeForceAver[0]<<endl;
        //fout<<endl;
        //cout<<__LINE__<<endl;
        //getchar();
        
        beamVec[j].lRWakeForceAver[0] *=  ElectronCharge / electronBeamEnergy ;  // [V/C] * [C] * [1e] / [eV] ->rad
        beamVec[j].lRWakeForceAver[1] *=  ElectronCharge / electronBeamEnergy ;  // [V/C] * [C] * [1e] / [eV] ->rad
        beamVec[j].lRWakeForceAver[2] *=  ElectronCharge / electronBeamEnergy ;  // [V/C] * [C] * [1e] / [eV] ->rad

        //cout<<setw(15)<<left<<beamVec[i].xAver
        //    <<setw(15)<<left<<beamVec[i].lRWakeForceAver[0]
        //    <<setw(15)<<left<<beamVec[i].lRWakeForceAver[1]
        //    <<setw(15)<<left<<beamVec[i].lRWakeForceAver[2]
        //    <<setw(15)<<left<<__LINE__<<__FILE__<<endl;
    }
    //getchar();
    //fout.close();

    // Reference to "simulation of transverse multi-bunch instabilities of proton beam in LHC, PHD thesis, Alexander Koshik P. 32, Eq. (3.22)"
} 
void SPBeam::BeamTransferPerTurnDueWake()
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToWake();
    }
}