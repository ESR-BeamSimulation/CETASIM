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
#include "BoardBandImp.h"
#include "Ramping.h"
#include "CUDAFunction.cuh"
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
#include <fftw3.h>
#include <complex.h>




using namespace std;
using std::vector;
using std::complex;


MPBeam::MPBeam()
{
}

MPBeam::~MPBeam()
{
      delete strongStrongBunchInfo;
      delete quasiWakePoten;
    //   delete partCord;
}


void MPBeam::Initial(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter)
{
    
    int totBunchNum = inputParameter.ringFillPatt->totBunchNumber; 
    int harmonics = inputParameter.ringParBasic->harmonics;              
    beamVec.resize(totBunchNum);
    freXIQDecompScan.resize(totBunchNum);
    freYIQDecompScan.resize(totBunchNum);
    freZIQDecompScan.resize(totBunchNum);
    bunchZMinZMax.resize(totBunchNum);
    for(int i=0;i<totBunchNum;i++)
    {
        bunchZMinZMax[i].resize(2);
    }

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
    // -----------------------------------------------------------------------

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

    // set the bunch initial distribution and prepare partCord for GPU
    int totMacroPartNum = 0;
    for(int i=0;i<totBunchNum;i++)
    {
        beamVec[i].InitialMPBunch(inputParameter);
        beamVec[i].DistriGenerator(latticeInterActionPoint,inputParameter,i);
        totMacroPartNum =beamVec[i].macroEleNumPerBunch; 
    }
    // partCord = new double(6*totMacroPartNum);

    // -----------------------------------------------------------------------
    
    for(int i=0;i<totBunchNum-1;i++)
    {
        beamVec[i].bunchGap = beamVec[i+1].bunchHarmNum - beamVec[i].bunchHarmNum;
    }
    beamVec[totBunchNum-1].bunchGap = harmonics - beamVec[totBunchNum-1].bunchHarmNum;
    // -----------------------------------------------------------------------

    MPGetBeamInfo();
	GetTimeDisToNextBunchIntial(inputParameter); 
    
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

void MPBeam::GetTimeDisToNextBunchIntial(ReadInputSettings &inputParameter)
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
    complex<double> vbKickAver[2];

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
        vbKickAver[i] =   complex<double>(0,0);

        for(int n=0;n<nTurns;n++)
        {
            int k=0;    
            for(int j=0;j<ringHarmH;j++)
            {

                if( k<beamVec.size() && beamVec[k].bunchHarmNum == j)
                {
                    vb0  = complex<double>(-1 * cavityResonator.resonatorVec[i].resFre * 2 * PI * cavityResonator.resonatorVec[i].resShuntImpRs
                                              /  cavityResonator.resonatorVec[i].resQualityQ0, 0.E0) *  beamVec[k].electronNumPerBunch * ElectronCharge;       
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
        vbKickAver[i]  =  vbKickAccum /double(beamVec.size());                  // average energy that bunch is kicked by self-loss Vb0/2. 
    }

    for(int i=0; i<resNum;i++)
    {

        if (cavityResonator.resonatorVec[i].resType==0) // passive cavity
        {
            cavityResonator.resonatorVec[i].resGenVol=complex<double>(0.E0,0.E0);
        }
        else                                           // active cavity
        {
            // double cavArg = arg(cavityResonator.resonatorVec[i].resCavVolReq);
            // double cavAbs = (cavityResonator.resonatorVec[i].resCavVolReq.real() - vbKickAver[i].real()) / cos( cavArg );
            // genAddvbAbs = abs(cavAbs);
            // cavityResonator.resonatorVec[i].resGenVol = cavAbs * exp(li * cavArg ) - cavityResonator.resonatorVec[i].vbAccum;
            
            cavityResonator.resonatorVec[i].resGenVol = cavityResonator.resonatorVec[i].resCavVolReq - cavityResonator.resonatorVec[i].vbAccum;
        }

        if(cavityResonator.resonatorVec[i].resCold || cavityResonator.resonatorVec[i].resRfMode==0)   // cavity is cold or ideal cavity
        {
            cavityResonator.resonatorVec[i].vbAccum = complex<double>(0.E0,0.E0);
        }
        cavityResonator.resonatorVec[i].GetInitialResonatorGenIg();
        cavityResonator.resonatorVec[i].GetInitialResonatorPower(inputParameter);

        // set initial cavity voltage, beam indcued volage for each bunch         
        for(int j=0;j<beamVec.size();j++)
        {
            beamVec[j].bunchRFModeInfo->cavVolBunchCen[i]    =     cavityResonator.resonatorVec[i].resCavVolReq;     // requried voltage.
            beamVec[j].bunchRFModeInfo->induceVolBunchCen[i] =     cavityResonator.resonatorVec[i].vbAccum;
            beamVec[j].bunchRFModeInfo->genVolBunchAver[i]   =     cavityResonator.resonatorVec[i].resGenVol;
            beamVec[j].bunchRFModeInfo->selfLossVolBunchCen[i] =    vbKickAver[i];
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
    fout<<"&parameter name=GenAmpIg,         units=Ampere, type=float,  &end"<<endl;
    fout<<"&parameter name=GenArgIg,         units=rad, type=float,  &end"<<endl;
    fout<<"&parameter name=LoadAnlge,        units=rad, type=float,  &end"<<endl;
    fout<<"&parameter name=BeamPower,        units=watt,type=float,  &end"<<endl;
    fout<<"&parameter name=GenPower,         units=watt,type=float,  &end"<<endl;
    fout<<"&parameter name=CavPower,         units=watt,type=float,  &end"<<endl;
    fout<<"&parameter name=GenReflectPower,  units=watt,type=float,  &end"<<endl;

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

void MPBeam::GetQuasiWakePoten(const ReadInputSettings &inputParameter,const BoardBandImp &boardBandImp)
{
    double rBeta              = inputParameter.ringParBasic->rBeta;
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    double zMinBin = -boardBandImp.zMax;
    double dzBin   =  boardBandImp.dz;
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double rfLen   = inputParameter.ringParBasic->t0 / ringHarmH * CLight;
    int    nBins   =  boardBandImp.binPosZ.size(); // bins for bunch profile, that is 2*boardBandImp.zZImp.size() - 1
    double rmsBunchLength = inputParameter.ringBBImp->quasiWakeBunchLen;
    quasiWakePoten->bunchLength = rmsBunchLength;

    int indexStart =  int ( ( -rfLen / 2. - zMinBin + dzBin / 2 ) / dzBin);
    int indexEnd   =  int ( (  rfLen / 2. - zMinBin + dzBin / 2 ) / dzBin);
    int N = indexEnd - indexStart;
    double wz[N],wDx[N],wDy[N],wQx[N],wQy[N];


    double temp[nBins];
    double wakePoten[nBins];
    
    for(int i=0;i<nBins;i++)
    {
        temp[i] = 1.0/sqrt(2*PI)/rmsBunchLength * exp(-pow(boardBandImp.binPosZ[i],2)/2/pow(rmsBunchLength,2)) * CLight ;   // [C/s] 
    }
    
    fftw_complex *r2cout  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nBins /2 +1) );  // the same size as boardBandImp.zZImp.size() // bunch specturm
    fftw_complex *c2rin   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nBins /2 +1) );  // the same size as boardBandImp.zZImp.size()
    
    fftw_plan p = fftw_plan_dft_r2c_1d(nBins, temp, r2cout, FFTW_ESTIMATE);    // ffw, according to g++/nvcc complile flag is FFTW_MEASURE/FFTW_ESTIMATE
    fftw_execute(p);


    // (0) - get longitudianl quasi-wake-poten r2cout stores the pspectrum info         
    for(int i=0;i<boardBandImp.zZImp.size();i++)
    {
        complex<double> indVF = complex<double> (r2cout[i][0], r2cout[i][1] ) * boardBandImp.zZImp[i];
        // complex<double> indVF = boardBandImp.zZImp[i];
        c2rin[i][0] = indVF.real();
        c2rin[i][1] = indVF.imag();                                    // [C/s] * [Ohm] = [V]   
    }
    
    p = fftw_plan_dft_c2r_1d(nBins, c2rin, wakePoten, FFTW_ESTIMATE);
    fftw_execute(p);
    for(int i=0;i<N;i++)   wz[i] = wakePoten[i+indexStart] / nBins;

    //(1) transverse dipole x 
    for(int i=0;i<boardBandImp.zZImp.size();i++)
    {
        complex<double> indVF = complex<double> (r2cout[i][0], r2cout[i][1] ) * boardBandImp.zDxImp[i] * li;
        c2rin[i][0] = indVF.real();
        c2rin[i][1] = indVF.imag();                              // [C/s] * [Ohm] = [V]   
    }
    p = fftw_plan_dft_c2r_1d(nBins, c2rin, wakePoten, FFTW_ESTIMATE);
    fftw_execute(p); 
    for(int i=0;i<N;i++)  wDx[i] = wakePoten[i+indexStart] / nBins;


    //(2) transverse dipole y 
    for(int i=0;i<boardBandImp.zZImp.size();i++)
    {
        complex<double> indVF = complex<double> (r2cout[i][0], r2cout[i][1] ) * boardBandImp.zDyImp[i] * li;
        c2rin[i][0] = indVF.real();
        c2rin[i][1] = indVF.imag();                              // [C/s] * [Ohm] = [V]   
    }
    p = fftw_plan_dft_c2r_1d(nBins, c2rin, wakePoten, FFTW_ESTIMATE);
    fftw_execute(p); 
    for(int i=0;i<N;i++)  wDy[i] = wakePoten[i+indexStart] / nBins;

    //(3) transverse quad x 
    for(int i=0;i<boardBandImp.zZImp.size();i++)
    {
        complex<double> indVF = complex<double> (r2cout[i][0], r2cout[i][1] ) * boardBandImp.zQxImp[i] * li;
        c2rin[i][0] = indVF.real();
        c2rin[i][1] = indVF.imag();                              // [C/s] * [Ohm] = [V]   
    }
    p = fftw_plan_dft_c2r_1d(nBins, c2rin, wakePoten, FFTW_ESTIMATE);
    fftw_execute(p); 
    for(int i=0;i<N;i++)  wQx[i] = wakePoten[i+indexStart] / nBins;

    //(4) transverse quad y 
    for(int i=0;i<boardBandImp.zZImp.size();i++)
    {
        complex<double> indVF = complex<double> (r2cout[i][0], r2cout[i][1] ) * boardBandImp.zQyImp[i] * li;
        c2rin[i][0] = indVF.real();
        c2rin[i][1] = indVF.imag();                              // [C/s] * [Ohm] = [V]   
    }
    p = fftw_plan_dft_c2r_1d(nBins, c2rin, wakePoten, FFTW_ESTIMATE);
    fftw_execute(p); 
    for(int i=0;i<N;i++)  wQy[i] = wakePoten[i+indexStart] / nBins;


    fftw_destroy_plan(p);
    fftw_free(r2cout);
    fftw_free(c2rin);


    // re-generate a quasi-wake from impedance here
    double poszHead =   6 * rmsBunchLength; // wake starts point at poz=-6 rms bunch length. 
    int indexTail   =   0; 
    int indexHead   =   int ( ( poszHead  + rfLen / 2  ) / dzBin);
    int indexCen    =   int ( ( 0.0       + rfLen / 2  ) / dzBin);
    int flipNum     =   indexHead - indexCen - 1;  
    int dim         =   indexHead - indexTail;
    

    quasiWakePoten->wz.resize(dim);
    quasiWakePoten->wDx.resize(dim);
    quasiWakePoten->wDy.resize(dim);
    quasiWakePoten->wQx.resize(dim);
    quasiWakePoten->wQy.resize(dim);
    quasiWakePoten->binPosZ.resize(dim);

    for(int i=0;i<dim;i++)
    {
        quasiWakePoten->binPosZ[i]  = i * dzBin - dim * dzBin;
        quasiWakePoten->wz[i]  = wz[i];
        quasiWakePoten->wDx[i] = wDx[i];
        quasiWakePoten->wDy[i] = wDy[i];
        quasiWakePoten->wQx[i] = wQx[i];
        quasiWakePoten->wQy[i] = wQy[i];
    }
   
    // flip the longitudianl wake poten
    // quasiWakePoten->wz[dim-1] = wz[indexCen];
    for(int i=0;i<flipNum;i++)
    {
        quasiWakePoten->wz[i]  = 0.E0;
    }
    for(int i=flipNum;i<indexCen;i++)
    {
        quasiWakePoten->wz[i] = wz[i-flipNum];
    }
    
    for(int i=indexCen;i<dim-1;i++)
    {
        quasiWakePoten->wz[i] = wz[i-flipNum] + wz[2* (dim-1) - i - flipNum ];
    }    
    quasiWakePoten->wz[dim-1] = wz[indexCen];
    // flip end 

    // reverse here ensure that  
    // quasiWakePoten->wz ad quasiWakePoten->wDx... alin from head to tail as show in Fig.2.6, here head -> tail whne posZ->+infinit
    reverse(quasiWakePoten->binPosZ.begin(),quasiWakePoten->binPosZ.end());
    reverse(quasiWakePoten->wz.begin(),quasiWakePoten->wz.end());
    reverse(quasiWakePoten->wDx.begin(),quasiWakePoten->wDx.end());
    reverse(quasiWakePoten->wDy.begin(),quasiWakePoten->wDy.end());
    reverse(quasiWakePoten->wQx.begin(),quasiWakePoten->wQx.end());
    reverse(quasiWakePoten->wDy.begin(),quasiWakePoten->wQx.end());
    

    // get wakepoten from wakes
    // benckmark test of wake poten of differetn approach

    // int nBinBunchDen = rfLen / dzBin;
    // int nBinWake = quasiWakePoten->wz.size();
    
    // ofstream fout("wakePotenL_fliped.dat");
    // fout<<"SDDS1"<<endl;
    // fout<<"&column name=z,              units=m,              type=float,  &end" <<endl;
    // fout<<"&column name=profile,        units=m,              type=float,  &end" <<endl;
    // fout<<"&column name=indVZ,                                type=float,  &end" <<endl;
    // fout<<"&data mode=ascii, &end"                                               <<endl;

    // double densityProfile[nBinBunchDen];
    // for(int i=0;i<nBinBunchDen;i++) densityProfile[i]=0.E0;
    
    // for(int k=0;k<10;k++)
    // {
    //     fout<<"! page number "<<k + 1<<endl;
    //     fout<<nBinBunchDen<<endl;

    //     rmsBunchLength = (k+1) * 1.0E-3;
    //     for(int i=0;i<nBinBunchDen;i++)
    //     {
    //         double x          = i * dzBin - rfLen / 2.0;
    //         densityProfile[i] = 1.0/sqrt(2*PI)/rmsBunchLength * exp(-pow(x,2)/2/pow(rmsBunchLength,2)) ;   // [C/m] 
    //     }

    //     double wakePotenL = 0;

    //     for(int i=0;i<nBinBunchDen;i++)
    //     {               
    //         wakePotenL = 0;
    //         for(int j=i;j<nBinBunchDen;j++)    
    //         {
    //             int tij = j - i;
    //             if(tij>nBinWake) break;
    //             wakePotenL -= quasiWakePoten->wz[tij] * densityProfile[j] * dzBin;   //[v/c]*[c/m] * dz
    //         }
            
    //         fout<<setw(15)<<left<<i*dzBin-rfLen/ 2.0 
    //             <<setw(15)<<left<<densityProfile[i]
    //             <<setw(15)<<left<<wakePotenL
    //             <<endl;
    //     }
    // }

    // fout.close();


    // // get wakepoten from impedacen
    // ofstream fout1("wakePotenFre.dat");
    // fout1<<"SDDS1"<<endl;
    // fout1<<"&column name=z,              units=m,              type=float,  &end" <<endl;
    // fout1<<"&column name=profile,        units=m,              type=float,  &end" <<endl;
    // fout1<<"&column name=indVZ,                                type=float,  &end" <<endl;
    // fout1<<"&data mode=ascii, &end"                                               <<endl;

    // for(int k=0;k<10;k++)
    // {
    //     fout1<<"! page number "<<k + 1<<endl;
    //     fout1<<N<<endl;

    //     rmsBunchLength = (k+1) * 1.0E-3;
    //     for(int i=0;i<nBins;i++)
    //     {
    //         temp[i] = 1.0/sqrt(2*PI)/rmsBunchLength * exp(-pow(boardBandImp.binPosZ[i],2)/2/pow(rmsBunchLength,2)) * CLight ;   // [1/s] 
    //     }
        
    //     fftw_complex *r2cout  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nBins /2 +1) );  // the same size as boardBandImp.zZImp.size() // bunch specturm
    //     fftw_complex *c2rin   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (nBins /2 +1) );  // the same size as boardBandImp.zZImp.size()

    //     fftw_plan p = fftw_plan_dft_r2c_1d(nBins, temp, r2cout, FFTW_ESTIMATE);    // ffw, according to g++/nvcc complile flag is FFTW_MEASURE/FFTW_ESTIMATE
    //     fftw_execute(p);

    //     // (0) - get longitudianl quasi-wake-poten r2cout stores the pspectrum info     
    //     for(int i=0;i<boardBandImp.zZImp.size();i++)
    //     {
    //         complex<double> indVF = complex<double> (r2cout[i][0], r2cout[i][1] ) * boardBandImp.zZImp[i];
    //         c2rin[i][0] = indVF.real();
    //         c2rin[i][1] = indVF.imag();                              // [C/s] * [Ohm] = [V]   
    //     }
    //     p = fftw_plan_dft_c2r_1d(nBins, c2rin, wakePoten, FFTW_ESTIMATE);
    //     fftw_execute(p);
    //     for(int i=0;i<N;i++)   wz[i] = wakePoten[i+indexStart] / nBins;

    //     for(int i=0;i<N;i++)
    //     {
    //         fout1<<setw(15)<<left<<boardBandImp.binPosZ[i+indexStart]
    //             <<setw(15)<<left<<temp[i+indexStart]
    //             <<setw(15)<<left<<-wz[i]
    //             <<endl;
    //     }
    // }
    // fout1.close();

    // for(int i=0;i<N;i++) quasiWakePoten->binPosZ[i] = -rfLen / 2 + i * dzBin;

    // ofstream fout("wake_green_function1.dat");
    // fout<<"SDDS1"<<endl;
    // fout<<"&column name=z,              units=m,              type=float,  &end" <<endl;
    // fout<<"&column name=profile,        units=m,              type=float,  &end" <<endl;
    // fout<<"&column name=indVZ,                                type=float,  &end" <<endl;
    // fout<<"&column name=indVDX,                               type=float,  &end" <<endl;
    // fout<<"&column name=indVDY,                               type=float,  &end" <<endl;
    // fout<<"&column name=indVQX,                                type=float,  &end" <<endl;
    // fout<<"&column name=indVQY,                                type=float,  &end" <<endl;
    // fout<<"&data mode=ascii, &end"                                               <<endl;
    // fout<<"! page number " << 1 <<endl;
    // fout<<N<<endl;

    // for(int i=0;i<dim;i++) 
    // {
    //     fout<<setw(15)<<left<<quasiWakePoten->binPosZ[i]
    //         <<setw(15)<<left<<temp[i+indexStart]
    //         <<setw(15)<<left<<quasiWakePoten->wz[i]
    //         <<setw(15)<<left<<quasiWakePoten->wDx[i]
    //         <<setw(15)<<left<<quasiWakePoten->wDy[i]
    //         <<setw(15)<<left<<quasiWakePoten->wQx[i]
    //         <<setw(15)<<left<<quasiWakePoten->wQy[i]
    //         <<endl;
    // }
    // cout<<__LINE__<<__FILE__<<endl;
    // getchar();
}
void MPBeam::LongWakePotenBenchmark(const ReadInputSettings &inputParameter,const BoardBandImp &boardBandImp)
{

}


void MPBeam::GetQuasiWakePoten(const ReadInputSettings &inputParameter)
{
    string fileName = inputParameter.ringBBImp->wakeInput;

    int lineNumber = 0;
    string str;
    vector<string> strVec;
    int index = 0;
    int cenIndex=0;

    // file is from sddsprintout command and impedance is with ele
    ifstream fin(fileName);
    while (!fin.eof())
    {
        if(lineNumber<5)
        {
            getline(fin,str);		 
            lineNumber++;
            continue;
        }
        getline(fin,str);
        if(str.length()==0)  continue;
        StringSplit2(str,strVec);
        
        quasiWakePoten->binPosZ.push_back(stod(strVec[0]));
        quasiWakePoten->wz.push_back( -1*stod(strVec[1]));   // -1 change to Alex Chao definition
        quasiWakePoten->wDx.push_back(-1*stod(strVec[2]));
        quasiWakePoten->wDy.push_back(-1*stod(strVec[3]));
        quasiWakePoten->wQx.push_back(-1*stod(strVec[4]));
        quasiWakePoten->wQy.push_back(-1*stod(strVec[5]));
        if(stod(strVec[0])==0) cenIndex = index;
        index ++;
    }
    
    double z0 = quasiWakePoten->binPosZ[0];
    for(int i=0 ;i<quasiWakePoten->binPosZ.size(); i++) quasiWakePoten->binPosZ[i] -= z0; 

    // flip the longitudinal wakes
    int wakeBins = quasiWakePoten->binPosZ.size();
    vector<double> temp(index,0.E0);
    for(int i=0;i<temp.size();i++)  temp[i] = quasiWakePoten->wz[i];
    quasiWakePoten->wz[0] = temp[cenIndex] ; 
    for(int i=1;i<cenIndex;i++)                    
    {
        quasiWakePoten->wz[i] = temp[cenIndex + i] + temp[cenIndex - i];
    }
    for(int i=cenIndex;i<wakeBins - cenIndex;i++)  
    {
        quasiWakePoten->wz[i] = temp[cenIndex + i];
    }
    for(int i=wakeBins - cenIndex ;i<wakeBins;i++) 
    {
        quasiWakePoten->wz[i] =  0.E0;
    }
    // end of flips

    // quasiWakePoten->wz ad quasiWakePoten->wDx... alin from head to tail as show in Fig.2.6
    // quasiWakePoten->wz[0] is the head particle. 


    // ofstream fout("test.dat");
    // for(int i=0;i<quasiWakePoten->wz.size();i++)
    // {
    //     fout<<setw(15)<<left<<quasiWakePoten->binPosZ[i]
    //         <<setw(15)<<left<<quasiWakePoten->wz[i]
    //         <<setw(15)<<left<<quasiWakePoten->wDx[i]
    //         <<setw(15)<<left<<quasiWakePoten->wDy[i]
    //         <<setw(15)<<left<<quasiWakePoten->wQx[i]
    //         <<setw(15)<<left<<quasiWakePoten->wQy[i]    
    //         <<endl;
    // }

    // fout.close();
    // cout<<__LINE__<<__FILE__<<endl;
    // getchar();
}


void MPBeam::Run(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter, CavityResonator &cavityResonator)
{
    int nTurns                      = inputParameter.ringRun->nTurns;
    int synRadDampingFlag          = inputParameter.ringRun->synRadDampingFlag;
    int fIRBunchByBunchFeedbackFlag = inputParameter.ringRun->fIRBunchByBunchFeedbackFlag;
    int bBImpFlag                   = inputParameter.ringRun->bBImpFlag;
    int beamIonFlag                 = inputParameter.ringRun->beamIonFlag;
    int lRWakeFlag                  = inputParameter.ringRun->lRWakeFlag;
    int sRWakeFlag                  = inputParameter.ringRun->sRWakeFlag;
    int totBunchNum                 = inputParameter.ringFillPatt->totBunchNumber;
    int ionInfoPrintInterval        = inputParameter.ringIonEffPara->ionInfoPrintInterval;
    int bunchInfoPrintInterval      = inputParameter.ringRun->bunchInfoPrintInterval;
    int rampFlag                    = inputParameter.ringRun->rampFlag;

    // prepare the ramping class
    Ramping ramping;

    // preapre the data for bunch-by-bunch system ---------------   
    FIRFeedBack firFeedBack;
    if(fIRBunchByBunchFeedbackFlag)
    {
        firFeedBack.Initial(inputParameter);
    }
    
    // preapre the data for board band impedance ---------------   
    BoardBandImp boardBandImp;
    
    if(bBImpFlag)
    {
        if(inputParameter.ringBBImp->timeDomain==0)  // allocate the vector used to store the impedance data
        {
            boardBandImp.ReadInImp(inputParameter);
            for(int i=0;i<beamVec.size();i++)
            {
                beamVec[i].profileForBunchBBImp.resize(boardBandImp.nBins*2-1,0E0);
                beamVec[i].wakePotenFromBBI->wakePotenZ.resize(boardBandImp.nBins*2-1,0E0);
                beamVec[i].wakePotenFromBBI->wakePotenDx.resize(boardBandImp.nBins*2-1,0E0);
                beamVec[i].wakePotenFromBBI->wakePotenDy.resize(boardBandImp.nBins*2-1,0E0);
                beamVec[i].wakePotenFromBBI->wakePotenQx.resize(boardBandImp.nBins*2-1,0E0);
                beamVec[i].wakePotenFromBBI->wakePotenQy.resize(boardBandImp.nBins*2-1,0E0);         
            }
        }
        else
        {
            // GetQuasiWakePoten(inputParameter,boardBandImp);
            GetQuasiWakePoten(inputParameter);
            
        }      
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
    fout<<"&parameter name=beamCurr,    units=mA,   type=float,  &end"<<endl;
    fout<<"&parameter name=I1,          units=m,   type=float,  &end"<<endl;
    fout<<"&parameter name=I2,          units=1/m, type=float,  &end"<<endl;
    fout<<"&parameter name=I3,          units=1/m^2,   type=float,  &end"<<endl;
    fout<<"&parameter name=I4,          units=1/m, type=float,  &end"<<endl;
    fout<<"&parameter name=I5,          units=1/m, type=float,  &end"<<endl;
    fout<<"&parameter name=Jx,                     type=float,  &end"<<endl;
    fout<<"&parameter name=Jy,                     type=float,  &end"<<endl;
    fout<<"&parameter name=Jz,                     type=float,  &end"<<endl;
    fout<<"&parameter name=f1001Abs,               type=float,  &end"<<endl;
    fout<<"&parameter name=f1010Abs,               type=float,  &end"<<endl;
    fout<<"&parameter name=f0110Abs,               type=float,  &end"<<endl;
    fout<<"&parameter name=f0101Abs,               type=float,  &end"<<endl;
    fout<<"&parameter name=f1001Arg,  units=rad,   type=float,  &end"<<endl;
    fout<<"&parameter name=f1010Arg,  units=rad,   type=float,  &end"<<endl;
    fout<<"&parameter name=f0110Arg,  units=rad,   type=float,  &end"<<endl;
    fout<<"&parameter name=f0101Arg,  units=rad,   type=float,  &end"<<endl;
    fout<<"&parameter name=LCFactor,               type=float,  &end"<<endl;
    fout<<"&column name=Turns,                     type=long,   &end"<<endl;    
    fout<<"&column name=IonCharge,      units=e,   type=float,  &end"<<endl;
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
    fout<<"&column name=Nux,                       type=float,  &end"<<endl;
    fout<<"&column name=Nuy,                       type=float,  &end"<<endl;
    fout<<"&column name=deltaNu,                   type=float,  &end"<<endl;
    fout<<"&column name=f1001Abs,               type=float,  &end"<<endl;
    fout<<"&column name=f1010Abs,               type=float,  &end"<<endl;
    fout<<"&column name=f1001Arg,  units=rad,   type=float,  &end"<<endl;
    fout<<"&column name=f1010Arg,  units=rad,   type=float,  &end"<<endl;
    fout<<"&column name=LCFactor,               type=float,  &end"<<endl;
    fout<<"&column name=xyAlpha,   units=rad,   type=float,  &end"<<endl;
    fout<<"&column name=traceA,                 type=float,  &end"<<endl;
    fout<<"&column name=traceB,                 type=float,  &end"<<endl;
    fout<<"&column name=traceMSubN,             type=float,  &end"<<endl;
    fout<<"&column name=gamma0,                 type=float,  &end"<<endl;
    fout<<"&column name=gamma1,                 type=float,  &end"<<endl;
    fout<<"&column name=detH,                   type=float,  &end"<<endl;

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
        colname = string("&column name=") + string("eigenEmitx_")    + to_string(j) + string(", units=m*rad,   type=float,  &end");            
        fout<<colname<<endl;
        colname = string("&column name=") + string("eigenEmity_")    + to_string(j) + string(", units=m*rad,    type=float,  &end");            
        fout<<colname<<endl;
        colname = string("&column name=") + string("transmission_")  + to_string(j) + string(",                type=float,  &end");            
        fout<<colname<<endl;
        colname = string("&column name=") + string("xyCouplingAlpha_")  + to_string(j) + string(", units=rad,  type=float,  &end");            
        fout<<colname<<endl;
       
    }
    fout<<"&data mode=ascii, &end"<<endl;
    fout<<inputParameter.ringParBasic->ringCurrent <<endl;
    for(int i=0;i<5;i++) fout<<inputParameter.ringParBasic->radIntegral[i]<<endl;
    for(int i=0;i<3;i++) fout<<inputParameter.ringParBasic->dampingPartJ[i]<<endl;
    fout<<abs(latticeInterActionPoint.resDrivingTerms->f1001)<<endl;
    fout<<abs(latticeInterActionPoint.resDrivingTerms->f1010)<<endl;
    fout<<abs(latticeInterActionPoint.resDrivingTerms->f0110)<<endl;
    fout<<abs(latticeInterActionPoint.resDrivingTerms->f0101)<<endl;
    fout<<arg(latticeInterActionPoint.resDrivingTerms->f1001)<<endl;
    fout<<arg(latticeInterActionPoint.resDrivingTerms->f1010)<<endl;
    fout<<arg(latticeInterActionPoint.resDrivingTerms->f0110)<<endl;
    fout<<arg(latticeInterActionPoint.resDrivingTerms->f0101)<<endl;
    fout<<latticeInterActionPoint.resDrivingTerms->linearCouplingFactor<<endl;

    fout<<"! page number "<<1<<endl;
    fout<<nTurns<<endl;


    // run loop starts, for nTrns and each trun for k interaction-points
    for(int n=0;n<nTurns;n++)
    {
        if(n%10==0) cout<<n<<"  turns, bunch_0 transmission: "<<beamVec[0].transmission <<endl;

        MPBeamRMSCal(latticeInterActionPoint, 0);
        MPGetBeamInfo();
                
        if(beamIonFlag) 
        {
            for (int k=0;k<inputParameter.ringIonEffPara->numberofIonBeamInterPoint;k++)
            {
                MPBeamRMSCal(latticeInterActionPoint, k);
                MPGetBeamInfo();
                SSBeamIonEffectOneInteractionPoint(inputParameter,latticeInterActionPoint, n, k);
                BeamTransferPerInteractionPointDueToLatticeT(inputParameter,latticeInterActionPoint,k);             //transverse transfor per interaction point                
            }

            MPBeamRMSCal(latticeInterActionPoint,0);
             
            if(ionInfoPrintInterval && (n%ionInfoPrintInterval==0))
            {
                MPBeamRMSCal(latticeInterActionPoint, 0);
                SSIonDataPrint(inputParameter,latticeInterActionPoint, n);
            }
        }
        else
        {
            // BeamTransferPerTurnDueToLatticeT(inputParameter,latticeInterActionPoint);
            // BeamTransferPerTurnR66AndSynRadGPU(inputParameter,latticeInterActionPoint);
            BeamTransferPerTurnDueToLatticeTOneTurnR66(inputParameter,latticeInterActionPoint);
            MPBeamRMSCal(latticeInterActionPoint, 0);
        }

        // BeamMomtumUpdateDueToRF(inputParameter,latticeInterActionPoint,cavityResonator);
        BeamMomtumUpdateDueToRFTest(inputParameter,latticeInterActionPoint,cavityResonator);
        BeamLongiPosTransferOneTurn(inputParameter);
        BeamEnergyLossOneTurn(inputParameter);

        // Subroutine in below only change the momentum 
        if(bBImpFlag)  BBImpBeamInteraction(inputParameter,boardBandImp,latticeInterActionPoint);

        if(lRWakeFlag) LRWakeBeamIntaction(inputParameter,lRWakeFunction,latticeInterActionPoint);
           
        if(sRWakeFlag) SRWakeBeamIntaction(inputParameter,sRWakeFunction,latticeInterActionPoint,n);
             
        if(fIRBunchByBunchFeedbackFlag) FIRBunchByBunchFeedback(inputParameter,firFeedBack,n);

        if(rampFlag) ramping.RampingPara(inputParameter,latticeInterActionPoint,n);


        if((inputParameter.driveMode->driveModeOn!=0) && (inputParameter.driveMode->driveStart <n )  &&  (inputParameter.driveMode->driveEnd >n) )
        {
            BeamTransferDuetoDriveMode(inputParameter,n);
        }

        // update the eign emittnace for SR simulation
       
        if(synRadDampingFlag==1) BeamSynRadDamping(inputParameter,latticeInterActionPoint);
        
        MarkParticleLostInBunch(inputParameter,latticeInterActionPoint);   

        MPBeamRMSCal(latticeInterActionPoint, 0);
        MPGetBeamInfo();

        double nux = inputParameter.ringParBasic->workQx - floor(inputParameter.ringParBasic->workQx); 
        double nuy = inputParameter.ringParBasic->workQy - floor(inputParameter.ringParBasic->workQy);
        fout<<setw(15)<<left<<n
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
            <<setw(15)<<left<< strongStrongBunchInfo->bunchRmsSizePZ   // over all bunches
            <<setw(15)<<left<< nux
            <<setw(15)<<left<< nuy
            <<setw(15)<<left<< nuy-nux
            <<setw(15)<<left<< abs(latticeInterActionPoint.resDrivingTerms->f1001)
            <<setw(15)<<left<< abs(latticeInterActionPoint.resDrivingTerms->f1010)
            <<setw(15)<<left<< arg(latticeInterActionPoint.resDrivingTerms->f1001)
            <<setw(15)<<left<< arg(latticeInterActionPoint.resDrivingTerms->f1010)
            <<setw(15)<<left<< latticeInterActionPoint.resDrivingTerms->linearCouplingFactor
            <<setw(15)<<left<< latticeInterActionPoint.resDrivingTerms->xyAlpha
            <<setw(15)<<left<< latticeInterActionPoint.traceAB[0]
            <<setw(15)<<left<< latticeInterActionPoint.traceAB[1]
            <<setw(15)<<left<< latticeInterActionPoint.traceAB[2]
            <<setw(15)<<left<< latticeInterActionPoint.gammaC[0]
            <<setw(15)<<left<< latticeInterActionPoint.gammaC[1]
            <<setw(15)<<left<< latticeInterActionPoint.detH;

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
                <<setw(15)<<left<<beamVec[index].emittanceY
                <<setw(15)<<left<<beamVec[index].eigenEmitX
                <<setw(15)<<left<<beamVec[index].eigenEmitY
                <<setw(15)<<left<<beamVec[index].transmission
                <<setw(15)<<left<<beamVec[index].xyCouplingAlpha;

        }
        fout<<endl;        
        
        if(bunchInfoPrintInterval && (n%bunchInfoPrintInterval==0)  )
        {             
            MPBeamDataPrintPerTurn(n,latticeInterActionPoint,inputParameter); 
            if(inputParameter.driveMode->driveModeOn!=0)
            {
                GetDriveModeGrowthRate(n,inputParameter);
            }
            if(!inputParameter.ringRun->runCBMGR.empty())
            {
                GetCBMGR(n,latticeInterActionPoint,inputParameter);
            }          
        }                                       
    }
    fout.close();

    cout<<"End of Tracking "<<nTurns<< "Turns"<<endl;


    // if( !inputParameter.ringRun->TBTBunchHaissinski.empty() &&  !cavityResonator.resonatorVec.empty())
    // {
    //     GetHaissinski(inputParameter,cavityResonator,sRWakeFunction); 
    //     if(!inputParameter.ringRun->TBTBunchLongTraj.empty())
    //     {
    //         GetAnalyticalLongitudinalPhaseSpace(inputParameter,cavityResonator,sRWakeFunction);
    //     }
    //     cout<<"Haissinski and longitudinal phase space trajectory"<<endl;
    // } 

}


void MPBeam::TuneRamping(ReadInputSettings &inputParameter, double turns)
{
    int nTurns       = inputParameter.ringRun->nTurns;       
    if(inputParameter.ramping->rampingNu[0]==1 ) inputParameter.ringParBasic->workQx +=  inputParameter.ramping->deltaNuPerTurn[0];
    if(inputParameter.ramping->rampingNu[1]==1 ) inputParameter.ringParBasic->workQy +=  inputParameter.ramping->deltaNuPerTurn[1];
}

void MPBeam::BeamLongiPosTransferOneTurn(const ReadInputSettings &inputParameter)
{
    for(int i=0;i<beamVec.size();i++) 
        beamVec[i].BunchLongPosTransferOneTurn(inputParameter);
}


void MPBeam::BeamSynRadDamping(const ReadInputSettings &inputParameter, LatticeInterActionPoint &latticeInterActionPoint)
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchSynRadDamping(inputParameter,latticeInterActionPoint);
    }

    // GPU version of synRadDamping simulation
    // int totalPartiNum=0;
    // for(int i=0;i<beamVec.size();i++)
    // {
    //     totalPartiNum += beamVec[i].macroEleNumPerBunch;
    // }
    // double partCord[6*totalPartiNum];
    
    // int  paraNum = sizeof(latticeInterActionPoint.latticeSynRadBRH)/sizeof(latticeInterActionPoint.latticeSynRadBRH[0]);

    // CopyPartCordToGPU(partCord,totalPartiNum); 
    // GPU_PartiSynRad(totalPartiNum,partCord,paraNum,latticeInterActionPoint.latticeSynRadBRH,latticeInterActionPoint.latticeSynRadCoff);
    // CopyPartCordFromGPU(partCord);
}


void MPBeam::BeamTransferPerTurnR66AndSynRadGPU(const ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint)
{
    int totalPartiNum=0;
    for(int i=0;i<beamVec.size();i++)
    {
        totalPartiNum += beamVec[i].macroEleNumPerBunch;
    }
    double partCord[6*totalPartiNum];
    int  oneTurnMatrixParaNum = sizeof(latticeInterActionPoint.latticeParaForOneTurnMap)/sizeof(latticeInterActionPoint.latticeParaForOneTurnMap[0]);
    
    CopyPartCordToGPU(partCord,totalPartiNum);
    GPU_PartiOneTurnTransferAndSynRad(totalPartiNum,partCord,oneTurnMatrixParaNum,latticeInterActionPoint.latticeParaForOneTurnMap,
                                                                                  latticeInterActionPoint.latticeSynRadBRH,
                                                                                  latticeInterActionPoint.latticeSynRadCoff);
    CopyPartCordFromGPU(partCord);

}

void MPBeam::BeamTransferPerTurnDueToLatticeTOneTurnR66(const ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint)
{
    for(int i=0;i<beamVec.size();i++)
    {
        //  beamVec[i].BunchTransferDueToLatticeOneTurnT66GPU(inputParameter,latticeInterActionPoint);
        beamVec[i].BunchTransferDueToLatticeOneTurnT66(inputParameter,latticeInterActionPoint);
    }

    // GPU version of the beam transfer
    // int totalPartiNum=0;
    // for(int i=0;i<beamVec.size();i++)
    // {
    //     totalPartiNum += beamVec[i].macroEleNumPerBunch;
    // }
    // double partCord[6*totalPartiNum];
    // int  paraNum = sizeof(latticeInterActionPoint.latticeParaForOneTurnMap)/sizeof(latticeInterActionPoint.latticeParaForOneTurnMap[0]);
    
    // CopyPartCordToGPU(partCord,totalPartiNum);
    // GPU_PartiOneTurnTransfer(totalPartiNum,partCord,paraNum,latticeInterActionPoint.latticeParaForOneTurnMap);
    // CopyPartCordFromGPU(partCord);
   
    
}
void MPBeam::CopyPartCordFromGPU(double *partCord)
{
    int indStart = 0;
    for(int i=0;i<beamVec.size();i++)
    {
        for(int j=0;j<beamVec[i].macroEleNumPerBunch;j++)
        {
            beamVec[i].ePositionX[j] = partCord[indStart + 6*j  ] ;
            beamVec[i].eMomentumX[j] = partCord[indStart + 6*j+1] ;
            beamVec[i].ePositionY[j] = partCord[indStart + 6*j+2] ;
            beamVec[i].eMomentumY[j] = partCord[indStart + 6*j+3] ;
            beamVec[i].ePositionZ[j] = partCord[indStart + 6*j+4] ;
            beamVec[i].eMomentumZ[j] = partCord[indStart + 6*j+5] ;        
        }
        indStart += 6 * beamVec[i].macroEleNumPerBunch; 
    }
}

void MPBeam::CopyPartCordToGPU(double *partCord, int totalPartiNum)
{
    int indStart = 0;
    for(int i=0;i<beamVec.size();i++)
    {
        for(int j=0;j<beamVec[i].macroEleNumPerBunch;j++)
        {
            partCord[indStart + 6*j  ] =  beamVec[i].ePositionX[j];
            partCord[indStart + 6*j+1] =  beamVec[i].eMomentumX[j];
            partCord[indStart + 6*j+2] =  beamVec[i].ePositionY[j];
            partCord[indStart + 6*j+3] =  beamVec[i].eMomentumY[j];
            partCord[indStart + 6*j+4] =  beamVec[i].ePositionZ[j];
            partCord[indStart + 6*j+5] =  beamVec[i].eMomentumZ[j];        
        
        }
        indStart += 6 * beamVec[i].macroEleNumPerBunch;
    }
    
}

void MPBeam::BBImpBeamInteraction(const ReadInputSettings &inputParameter, const BoardBandImp &boardBandImp, const LatticeInterActionPoint &latticeInterActionPoint)
{
    
    int wakeBins = quasiWakePoten->wz.size();
    vector<vector<double>> wakePoten(6,vector<double>(wakeBins,0));
    wakePoten[1] = quasiWakePoten->wz;
    wakePoten[2] = quasiWakePoten->wDx;
    wakePoten[3] = quasiWakePoten->wDy;
    wakePoten[4] = quasiWakePoten->wQx;
    wakePoten[5] = quasiWakePoten->wQy;
    wakePoten[0] = quasiWakePoten->binPosZ;


    for(int j=0;j<beamVec.size();j++)
    {
        if(beamVec[j].macroEleCharge==0) continue;
        if(inputParameter.ringBBImp->timeDomain==1)
        {
            beamVec[j].BBImpBunchInteractionTD(inputParameter,boardBandImp,latticeInterActionPoint,wakePoten);
        }
        else if(inputParameter.ringBBImp->timeDomain==0)
        {
            beamVec[j].BBImpBunchInteraction(inputParameter,boardBandImp,latticeInterActionPoint);
        }
        else cerr<<"wrong setting in BoardBandImpedance namelist"<<endl;    
    }

}

void MPBeam::BeamTransferDuetoDriveMode(const ReadInputSettings &inputParameter, const int n)
{     
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToDriveMode(inputParameter,n);     
    }
}

void MPBeam::MarkParticleLostInBunch(const ReadInputSettings &inputParameter, const LatticeInterActionPoint &latticeInterActionPoint)
{
    for(int i=0;i<beamVec.size();i++)
    {
        beamVec[i].MarkLostParticle(inputParameter,latticeInterActionPoint);
    }
}




void MPBeam::GetDriveModeGrowthRate(const int turns, const ReadInputSettings &inputParameter)
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


void MPBeam::GetCBMGR(const int turns, const LatticeInterActionPoint &latticeInterActionPoint, const ReadInputSettings &inputParameter)
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
void MPBeam::GetAnalyticalWithFilter(const ReadInputSettings &inputParameter)
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

vector<complex<double> > MPBeam::GetHilbertAnalytical(vector<double> signal, const double filterBandWithdNu,  double workQ)
{
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

void MPBeam::GetHaissinski(ReadInputSettings &inputParameter,CavityResonator &cavityResonator,WakeFunction &sRWakeFunction)
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


void MPBeam::GetAnalyticalLongitudinalPhaseSpace(ReadInputSettings &inputParameter,CavityResonator &cavityResonator,WakeFunction &sRWakeFunction)
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
            fout<<"&column name=BunchIndex,                       type=long,  &end"<<endl;
            fout<<"&column name=TotIonCharge,       units=e,     type=float,  &end"<<endl;
	        fout<<"&column name=AverX,              units=m,     type=float,  &end"<<endl;
	        fout<<"&column name=AverY,              units=m,     type=float,  &end"<<endl;
	        fout<<"&column name=AverZ,              units=m,     type=float,  &end"<<endl;
	        fout<<"&column name=AverXP,             units=rad,   type=float,  &end"<<endl;
	        fout<<"&column name=AverYP,             units=rad,   type=float,  &end"<<endl;
	        fout<<"&column name=AverZP,             units=rad,   type=float,  &end"<<endl;
	        fout<<"&column name=RmsEmitX,           units=m*rad, type=float,  &end"<<endl;
	        fout<<"&column name=RmsEmitY,           units=m*rad, type=float,  &end"<<endl;
	        fout<<"&column name=RmsEmitZ,           units=m*rad, type=float,  &end"<<endl;
	        fout<<"&column name=RmsRX,              units=m*rad, type=float,  &end"<<endl;
	        fout<<"&column name=RmsRY,              units=m*rad, type=float,  &end"<<endl;
	        fout<<"&column name=RmsBunchLength,     units=m,     type=float,  &end"<<endl;
     	    fout<<"&column name=RmsBunchEnergySpread, units=rad, type=float,  &end"<<endl;

     	    for(int j=0; j<inputParameter.ringParRf->resNum;j++)
            {
    	        string colname = string("&column name=") + string("cavAmp_") + to_string(j) + string(", units=V,   type=float,  &end");
                fout<<colname<<endl;
                colname = string("&column name=") + string("cavPhase_") + to_string(j) + string(", units=rad,   type=float,  &end");
                fout<<colname<<endl;
                colname = string("&column name=") + string("cavReal_") + to_string(j) + string(", units=V,   type=float,  &end");
                fout<<colname<<endl;

                colname = string("&column name=") + string("genAmp_") + to_string(j) + string(", units=V,   type=float,  &end");
                fout<<colname<<endl;
                colname = string("&column name=") + string("genPhase_") + to_string(j) + string(", units=rad,   type=float,  &end");
                fout<<colname<<endl;
                colname = string("&column name=") + string("beamIndAmp_") + to_string(j) + string(", units=V,   type=float,  &end");
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
                <<setw(24)<<left<<setprecision(16)<<i
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].totIonCharge
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].xAver
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].yAver
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].zAver
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].pxAver
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].pyAver
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].pzAver
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].emittanceX
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].emittanceY
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].emittanceZ
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].rmsRx
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].rmsRy
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].rmsBunchLength
                <<setw(24)<<left<<setprecision(16)<<beamVec[i].rmsEnergySpread;

            for(int j=0; j<inputParameter.ringParRf->resNum;j++)
            {
                fout<<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].bunchRFModeInfo->cavVolBunchCen[j])
                    <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].bunchRFModeInfo->cavVolBunchCen[j])
                    <<setw(24)<<left<<setprecision(16)<<    beamVec[i].bunchRFModeInfo->cavVolBunchCen[j].real()
                    <<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].bunchRFModeInfo->genVolBunchAver[j])
                    <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].bunchRFModeInfo->genVolBunchAver[j])
                    <<setw(24)<<left<<setprecision(16)<<abs(beamVec[i].bunchRFModeInfo->induceVolBunchCen[j])
                    <<setw(24)<<left<<setprecision(16)<<arg(beamVec[i].bunchRFModeInfo->induceVolBunchCen[j]);
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
    if(!inputParameter.ringRun->TBTBunchPro.empty()  &&  (TBTBunchPrintNum !=0) )
    {
        for(int i=0;i<TBTBunchPrintNum;i++)
        {
            string filePrefix = inputParameter.ringRun->TBTBunchPro;
            string fname = filePrefix +"_" +to_string(i)+".sdds";
            ofstream fout2(fname,ios_base::app);

            if(nTurns==0)
	        {
	            fout2<<"SDDS1"<<endl;
	            fout2<<"&column name=z,              units=m,     type=float,  &end"<<endl;
	            fout2<<"&column name=rho,            units=C/m,   type=float,  &end"<<endl;
	            // fout2<<"&column name=potenWell,                   type=float,  &end"<<endl;
	            // fout2<<"&column name=rhoAnalytical,               type=float,  &end"<<endl;

    	        // for(int j=0; j<inputParameter.ringParRf->resNum;j++)
	            // {

                //     string colname = string("&column name=") + string("cavReal_") + to_string(j) + string(", units=V,   type=float,  &end");
                //     fout2<<colname<<endl;
                //     colname = string("&column name=") + string("cavImag_") + to_string(j) + string(", units=V,   type=float,  &end");
                //     fout2<<colname<<endl;

                //     colname = string("&column name=") + string("cavAbs_") + to_string(j) + string(", units=V,   type=float,  &end");
                //     fout2<<colname<<endl;
                //     colname = string("&column name=") + string("cavArg_") + to_string(j) + string(", units=rad,   type=float,  &end");
                //     fout2<<colname<<endl;
                //     colname = string("&column name=") + string("cavDpz_") + to_string(j) + string(", units=rad,   type=float,  &end");
                //     fout2<<colname<<endl;
	            // }
	            // string colname = string("&column name=") + string("cavDpzTot")  + string(", units=rad,   type=float,  &end");
	            // fout2<<colname<<endl;
	            fout2<<"&data mode=ascii, &end"<<endl;
            }

            int bunchIndex  = TBTBunchDisDataBunchIndex[i];

            fout2<<"! page number "<< nTurns + 1<<endl;
            fout2<<inputParameter.ringParRf->rfBunchBinNum<<endl;
            
            for(int k=0;k<inputParameter.ringParRf->rfBunchBinNum;k++)
            {
                fout2<<setw(20)<<left<<beamVec[bunchIndex].posZBins[k]
                     <<setw(20)<<left<<beamVec[bunchIndex].densProfVsBin[k]<<endl;
                    //  <<setw(20)<<left<<beamVec[bunchIndex].hamiltonPotenWell[k]
                    //  <<setw(20)<<left<<beamVec[bunchIndex].densProfVsBinAnalytical[k];

                // double dpzTemp=0;
                // for(int j=0; j<inputParameter.ringParRf->resNum;j++)
                // {
                //     complex<double> temp = beamVec[bunchIndex].cavVolInfoVsLongBins[j][k];
                //     double temp1 =  beamVec[bunchIndex].cavForceInfoVsLongBins[j][k];
                //     fout2<<setw(20)<<left<<temp.real()
                //          <<setw(20)<<left<<temp.imag()
                //          <<setw(20)<<left<<abs(temp)
                //          <<setw(20)<<left<<arg(temp)
                //          <<setw(20)<<left<<temp1;

                //     dpzTemp +=temp1;
                // }
                // fout2<<setw(20)<<left<<dpzTemp<<endl;

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

void MPBeam::BeamTransferPerInteractionPointDueToLatticeT(const ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    for(int j=0;j<beamVec.size();j++)
    {
        beamVec[j].BunchTransferDueToLatticeT(inputParameter,latticeInterActionPoint,k);
    }
}


void MPBeam::BeamMomtumUpdateDueToRFTest(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,CavityResonator &cavityResonator)
{
    
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double circRing   = inputParameter.ringParBasic->circRing;
    double rfLen      = circRing  / ringHarmH;
    double tRF        = inputParameter.ringParBasic->t0 / ringHarmH; 
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double t0         = inputParameter.ringParBasic->t0;
    
    
    GetTimeDisToNextBunch(inputParameter);
    for (int j=0;j<inputParameter.ringParRf->resNum;j++)
    {     
        double oneTurnTime=0;
        if (cavityResonator.resonatorVec[j].resRfMode ==0) // rfca element ideal cavity 
        {
            for(int i=0;i<beamVec.size();i++)
            {
                beamVec[i].BunchMomentumUpdateDueToRFCA(inputParameter,cavityResonator.resonatorVec[j]);
            }
        }
        else  // rfmode element  -- beam loading in included in the tracking
        {
            for (int i=0;i<beamVec.size();i++)
            {
                // particle momentum update
                if(cavityResonator.resonatorVec[j].resExciteIntability==0)
                {
                    beamVec[i].BunchMomentumUpdateDueToRFModeStable(inputParameter,cavityResonator.resonatorVec[j],j);
                }
                else
                {
                    beamVec[i].BunchMomentumUpdateDueToRFMode(inputParameter,cavityResonator.resonatorVec[j],j);
                }

                // store the info sampled by the cavnty at n*tRF, stored the data for cavity feedbacks 
                for (int k=0;k<beamVec[i].bunchGap;k++)
                {
                    int harmonicIndex = beamVec[i].bunchHarmNum + k;
                    double dt = k * tRF  - beamVec[i].zMinCurrentTurn / CLight; 
                    cavityResonator.resonatorVec[j].GetResonatorInfoAtNTrf(harmonicIndex,dt);
                    cavityResonator.resonatorVec[j].ResonatorDynamics(tRF);  // cavity always is updated by tRf
                }

                // resonator beamInduced voltage updated till next bunch -- also control condition for instability excitation
                double timeTemp;

                if(cavityResonator.resonatorVec[j].resExciteIntability==0)
                {
                    timeTemp = beamVec[i].bunchGap * tRF;
                    // timeTemp = beamVec[i].timeFromCurrnetBunchToNextBunch;
                }
                else 
                {
                    timeTemp = beamVec[i].timeFromCurrnetBunchToNextBunch;
                }     
                cavityResonator.resonatorVec[j].GetBeamInducedVol(timeTemp);
     
            }           
        }    
    }

    // here the cavity feedback can be set here, the infomation is stored at cavityResonator.resonatorVec.deltaVCavSample[harmonicNum]. one turn at n*tRF 
    // direct feedback is treated in below. 
    for (int j=0;j<inputParameter.ringParRf->resNum;j++)
    {
        if(cavityResonator.resonatorVec[j].resRfMode ==1 && cavityResonator.resonatorVec[j].resDirFB==1)
        {
            for (int i=0;i<beamVec.size();i++)
            {
                beamVec[i].GetLongiKickDueToCavFB(inputParameter,cavityResonator.resonatorVec[j]);
            }
        }
    }

    
} 

void MPBeam::BeamEnergyLossOneTurn(const ReadInputSettings &inputParameter)
{
    for(int i=0;i<beamVec.size();i++)
    {
        beamVec[i].BunchEnergyLossOneTurn(inputParameter);
    }
} 

void MPBeam::GetTimeDisToNextBunch(const ReadInputSettings &inputParameter)
{
    // bunch length is included in the the calcualton
    // time to next bunch: t = bunchGap * tRf - (  beamVec[i].tMax     - beamVec[i+1].tMin)
    //                     t = bunchGap * tRf - (- beamVec[i].zMin / c + beamVec[i+1].zMax / c)
    //                     t = bunchGap * tRf + (  beamVec[i].zMin     - beamVec[i+1].zMax ) / c

    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double rBeta      = inputParameter.ringParBasic->rBeta;
    double tRF        = inputParameter.ringParBasic->t0 / ringHarmH; 

    GetBunchMinZMaxZ();

    if(beamVec.size()==1)
    {        
        beamVec[0].timeFromCurrnetBunchToNextBunch  =  beamVec[0].bunchGap * tRF + (beamVec[0].zMinCurrentTurn - beamVec[0].zMaxCurrentTurn) / CLight / rBeta;
        beamVec[0].zAverLastTurn                    =  beamVec[0].zAver;                                                   
    }
    else 
    {
        for(int i=0;i<beamVec.size();i++)
        {
            if(i<beamVec.size()-1)
            {
                beamVec[i].GetZMinMax();
                beamVec[i+1].GetZMinMax();
                beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zMinCurrentTurn - beamVec[i+1].zMaxCurrentTurn) / CLight / rBeta;
            }
            else
            {
                beamVec[i].GetZMinMax();
                beamVec[0].GetZMinMax();
                beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zMinCurrentTurn - beamVec[0  ].zMaxCurrentTurn) / CLight / rBeta;
            }           
        }        
    }
}





void MPBeam::BeamMomtumUpdateDueToRF(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,CavityResonator &cavityResonator)
{
   
    // cout<<__LINE__<<__FILE__<<endl;
    vector<double> resAmpFBRatioForTotSelfLoss = inputParameter.ringParRf->resAmpFBRatioForTotSelfLoss;
    
    vector<complex<double> > vbKickAver;
    vector<complex<double> > vbAccumAver;
    vector<complex<double> > selfLossToCompensate;
    
   
    complex<double> totSelfLoss (0.e0,0.e0);
    vbKickAver.resize(inputParameter.ringParRf->resNum);
    vbAccumAver.resize(inputParameter.ringParRf->resNum);
    selfLossToCompensate.resize(inputParameter.ringParRf->resNum);
    
    int ringHarmH     = inputParameter.ringParRf->ringHarm;
    double circRing   = inputParameter.ringParBasic->circRing;
    double rfLen      = circRing  / ringHarmH;
    double tRF        = inputParameter.ringParBasic->t0 / ringHarmH; 
    double rBeta      = inputParameter.ringParBasic->rBeta; 
    

    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        complex<double> vbKickAverSum=(0.e0,0.e0);
        complex<double> vbAccumAverSum=(0.e0,0.e0);

        for(int j=0;j<beamVec.size();j++)
        {
            vbKickAverSum      += beamVec[j].bunchRFModeInfo->selfLossVolBunchCen[i];
            vbAccumAverSum     += beamVec[j].bunchRFModeInfo->induceVolBunchCen[i];
        }
        vbKickAver[i]      = vbKickAverSum   / double(beamVec.size());
        vbAccumAver[i]     = vbAccumAverSum  / double(beamVec.size());
    }

     
    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {
        selfLossToCompensate[i]  =   vbKickAver[i];  
    }
    
    for(int i=0;i<inputParameter.ringParRf->resNum;i++)
    {        
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
   
   
    if(inputParameter.ringParRf->methodForVb=="rigid")
    {
        //  He Tianlong code approaches..    
        if(beamVec.size()==1)
        {        
            beamVec[0].timeFromCurrnetBunchToNextBunch  =  beamVec[0].bunchGap * tRF + (beamVec[0].zAverLastTurn  - beamVec[0].zAver) / CLight / rBeta;
            beamVec[0].zAverLastTurn                    =  beamVec[0].zAver;                                                   
            beamVec[0].BunchMomentumUpdateDuetoRFRigid(inputParameter,cavityResonator);    
            beamVec[0].GetMPBunchRMS(latticeInterActionPoint, 0);
        }
        else 
        {
            for(int i=0;i<beamVec.size();i++)
            {
                // prepare the condition for longitudinal tracking -- bunch-by-bunch, inside bunch bin-by-bin
                if(i<beamVec.size()-1)
                {
                    beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zAver - beamVec[i+1].zAver) / CLight / rBeta;
                }
                else
                {
                    beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zAver - beamVec[0  ].zAver) / CLight / rBeta;
                }
                beamVec[i].BunchMomentumUpdateDuetoRFRigid(inputParameter,cavityResonator);
                beamVec[i].GetMPBunchRMS(latticeInterActionPoint, 0);
            }        
        }   
    }
    else if (inputParameter.ringParRf->methodForVb=="soft")
    {
        // elegnat and MBtrack's approach; It requires initial beam induced voltage well calculated. 
        GetBunchMinZMaxZ();
        if(beamVec.size()==1)
        {        
            beamVec[0].timeFromCurrnetBunchToNextBunch  =  beamVec[0].bunchGap * tRF + (beamVec[0].zMinCurrentTurn - beamVec[0].zMaxCurrentTurn) / CLight / rBeta;
            beamVec[0].zAverLastTurn                    =  beamVec[0].zAver;                                                   
            beamVec[0].BunchMomentumUpdateDuetoRFBinByBin(inputParameter,cavityResonator);    
            beamVec[0].GetMPBunchRMS(latticeInterActionPoint, 0);
        }
        else 
        {
            for(int i=0;i<beamVec.size();i++)
            {
                // prepare the condition for longitudinal tracking -- bunch-by-bunch, inside bunch bin-by-bin
                if(i<beamVec.size()-1)
                {
                    beamVec[i].GetZMinMax();
                    beamVec[i+1].GetZMinMax();
                    beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zMinCurrentTurn - beamVec[i+1].zMaxCurrentTurn) / CLight / rBeta;
                    // beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zAver - beamVec[i+1].zAver) / CLight / rBeta;  
                }
                else
                {
                    beamVec[i].GetZMinMax();
                    beamVec[0].GetZMinMax();
                    beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zMinCurrentTurn - beamVec[0  ].zMaxCurrentTurn) / CLight / rBeta;
                    // beamVec[i].timeFromCurrnetBunchToNextBunch  = beamVec[i].bunchGap * tRF + (beamVec[i].zAver - beamVec[0  ].zAver) / CLight / rBeta;
                }
                
                beamVec[i].BunchMomentumUpdateDuetoRFBinByBin(inputParameter,cavityResonator);
                beamVec[i].GetMPBunchRMS(latticeInterActionPoint, 0);
            }        
        }  
    }
    else if (inputParameter.ringParRf->methodForVb=="noinstability")
    {
        
    }

    // cavity dynamics and cavity feedback is not included here yet.  It is alwasy "soft" methods here at this moment..  

}

void MPBeam::GetBunchMinZMaxZ()
{
    for(int i=0;i<beamVec.size();i++)
    {
        beamVec[i].GetZMinMax();
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



void MPBeam::FIRBunchByBunchFeedback(const ReadInputSettings &inputParameter,FIRFeedBack &firFeedBack,int nTurns)
{
    //y[0] = \sum_0^{N} a_k x[-k]. Ref. Nakamura's paper spring 8 notation here used is the same with Nakamura's paper
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


void MPBeam::LRWakeBeamIntaction(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction,const  LatticeInterActionPoint &latticeInterActionPoint)
{
    
    int nTurnswakeTrunction     = inputParameter.ringLRWake->nTurnswakeTrunction;
    int harmonics               = inputParameter.ringParBasic->harmonics;
    double electronBeamEnergy   = inputParameter.ringParBasic->electronBeamEnergy;
    double rBeta                = inputParameter.ringParBasic->rBeta;
    double tRF                  = inputParameter.ringParBasic->t0 / double(harmonics);

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
    double deltaTij=0;
    double tauijStastic;

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
                tempIndex1   = j - 1 ;
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

     BeamTransferPerTurnDueWake(); 
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