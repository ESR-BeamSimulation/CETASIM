//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#ifndef READINPUTSETTINGS_H
#define READINPUTSETTINGS_H

#include <string>
#include <iostream>
#include <vector>

using namespace std;
using std::vector;


class ReadInputSettings
{
public:
    ReadInputSettings();
    ~ReadInputSettings();

    int ParamRead(int argc, char *argv[]);

// read from initial file;
   //1) ring basic parameter 
    struct RingParBasic{
        // read from input
        double circRing;
        double workQx;
        double workQy;
        double  u0;
        double electronBeamEnergy;
        double synchRadDampTime[3];
        double pipeAperature[3];
        double chrom[2]; 
        double alphac[3];
        double aDTX[5];
        double aDTY[5];
        double aDTXY[5];
        double sdelta0;
        double couplingfactorXY;
        double naturalEmit;
        double naturalBunchLength= 0.E0;
        double tengR21 =0.E0;
        double skewQuadK=0;
        // calculated
        double  rGamma; 
        double  rBeta;  
        double  t0;    
        double  f0;    
        double  omega0;
        int     harmonics; 
        double  eta; 
        double workQz;
        double kappa;
        double sigmaT0;
        double ringCurrent;
        double betaFunAver[3];      // x y z
        double radIntegral[5];

        double emitNat[3];        // x y z
        double dampingPartJ[3];   // x y z
        int    ringSectNum;        // how many section the ring is cutted  
        string twissInput = "twiss.dat"; 
    };
    RingParBasic *ringParBasic = new RingParBasic;
   
    //2) rf basic parameter 
    struct RingParRf
    {
        // read from input

        int    ringHarm;
        
        // 2021-11-29 -- set the cavity togethre with generator
        int resNum;
        int rfBunchBinNum;        
        vector<int>    resHarm;
        vector<int>    rfMode;
        vector<int>    resType;
        vector<double> resVol;
        vector<double> resShuntImpRs;
        vector<double> resQualityQ0;
        vector<double> resCouplingBeta;
        vector<double> resDetuneFre;
        vector<int>    resCold;
        double synPhaseOffset;
        vector<double> resAmpFBRatioForTotSelfLoss;
        vector<double> resPhase;
        vector <int>    resExciteIntability;
        vector <int>    resDirFB;
        vector <double> resDirFBGain;
        vector <double> resDirFBPhase;
        vector <double> resDirFBDelay;
        string transResonParWriteTo;
        string methodForVb="soft";

        
    };     
    RingParRf *ringParRf = new RingParRf;
    
   
    //3) filling pattern set
    struct RingFillPatt
    {
        int trainNumber;
        int totBunchNumber;
        vector<int> bunchNumberPerTrain;
        vector<int> trainGaps;
        vector<int> bunchGaps; 
        int bunchChargeNum;
        vector<int> bunchChargeIndex;
        vector<int> bunchCharge;
    };
    RingFillPatt *ringFillPatt = new RingFillPatt;


   
    //4) initial bunch parameter 
    struct RingBunchPara
    {       
        double current;
        double emittanceX;
        double emittanceY;
        double kappa;
        double rmsBunchLength;
        double rmsEnergySpread;
        double initialStaticOffSet[6];
        double initialDynamicOffSet[6];     
        int distributionType;
        int macroEleNumPerBunch;
        double emittanceZ;  
        string bunchDisWriteTo;   
    };
    RingBunchPara *ringBunchPara = new RingBunchPara;
    

   // 5) ion effect settings
   
    struct RingIonEffPara
    {
        vector<double> ionMassNumber;
        vector<double> corssSectionEI;
        vector<double> gasPercent;
        int gasSpec;
        
        double ionLossBoundary;
        double ionMaxNumber;
        int numberofIonBeamInterPoint;
        int macroIonNumberGeneratedPerIP;
        int   ionInfoPrintInterval;
        int ionCalSCMethod; 
        string ionDisWriteTo;
        string twissInput="twiss.dat"; 
    };
    RingIonEffPara *ringIonEffPara = new RingIonEffPara;
    
   
    // 6) bunch-by-bunch feedbcak
    struct RingBBFB
    {    
        int mode = 0;
        double kickStrengthK[3];
        double gain;
        int delay;
        int taps;
        int nSections = 0;
        vector<int> start;
        vector<int> end;
        double kickerDispP;
        double kickerDisp;
        double fIRBunchByBunchFeedbackPowerLimit;                  // power wat limit on feedback   
        double fIRBunchByBunchFeedbackKickerImped;              // ohm 
        vector<double >  firCoeffx;
        vector<double >  firCoeffy;
        vector<double >  firCoeffz;
        vector<double >  firCoeffxy;   //Nakamura's TDLSF approaches, one group pickup and kicker for 2D                      
    };    
    RingBBFB * ringBBFB =  new RingBBFB;
        
    // 7.1 LRWake
    struct RingLRWake
    {
        string pipeGeoInput;
        string bbrInput;
        int nTurnswakeTrunction=100;
    };    
    RingLRWake * ringLRWake =  new RingLRWake;    

    // 7.2 SRWake
    struct RingSRWake
    {
        string pipeGeoInput;
        string bbrInput;
        string SRWWakePotenWriteTo;
        int SRWBunchBinNum=100;

    };    
    RingSRWake * ringSRWake =  new RingSRWake;    
    

    // 8 impedance 
    struct RingImpedance
    {
        int bunchBinNumberZ=100;
        double quasiWakeBunchLen=0.E0;
        string impedInput;
        string wakeInput;
        int timeDomain = 1;            // tracking in time domain, if not, tracking in freqeuncy domain 
        int impedSimFlag[5]={1,1,1,1,1};
    };    
    RingImpedance *ringBBImp =  new RingImpedance;    
    
    // 9 DRIVEDAMP 
    struct DriveMode
    {
        int driveModeOn = 0;
        int driveCBMIndex = 0;
        int driveEnd = 0;       // in unit of turns
        int driveStart = 0;     // in unit of turns
        double driveAmp = 0; 
        double driveFre = 0.E0; // Hz 
        int drivePlane = 0;
        int driveHW = 0;       
    };    
    DriveMode * driveMode =  new DriveMode;  
        
    // 10 run setting
    struct RingRun
    {
        int calSetting;
        int synRadDampingFlag;
        int fIRBunchByBunchFeedbackFlag;
        int beamIonFlag;
        int nTurns;
        int growthRateFittingStart;
        int growthRateFittingEnd;
        int bBImpFlag;
        int lRWakeFlag;
        int sRWakeFlag;
        int TBTBunchPrintNum;
        int rampFlag;
        int spaceChargeFlag = 0;
        int scMeshNum[3] = {32,32,33};
        vector<int> TBTBunchDisDataBunchIndex;
        string TBTBunchAverData;
        string TBTBunchDisData;
        string TBTBunchPro;
        string TBTBunchHaissinski;
        string TBTBunchLongTraj;
        string runCBMGR;
        

        
        int bunchInfoPrintInterval;
    };       
    RingRun * ringRun =  new RingRun; 

    // 11 run ramping
    struct Ramping
    {
        int rampingNu[2];       //x and y
        double deltaNuPerTurn[2]; //dnuyx and dnuy
        int rampingSKQ;
        double deltaSKQKPerTurn;
        int deltaTurns;
        int rampingTurns[2];
    };       
    Ramping * ramping =  new Ramping; 
};

#endif
