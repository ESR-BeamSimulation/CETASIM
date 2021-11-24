#ifndef READINPUTSETTINGS_H
#define READINPUTSETTINGS_H

#include<string>
#include<iostream>
#include<vector>


using namespace std;
using std::vector;

class ReadInputSettings
{
public:
    ReadInputSettings();
    ~ReadInputSettings();

    int ParamRead();

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
        double alphac;
        double sdelta0;
 
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
        double current;
        double betaFunAver[3];  
                            
     };
     RingParBasic *ringParBasic = new RingParBasic;
   
    //2) rf basic parameter 
    struct RingParRf
    {
        // read from input
        double rfBaseFreq;
        double vRF;
        int   harmRF;
        // default is the ideal value        
        double synPhase =0.E0;
        double synPhaseH=0.E0;
        double vRFH=0.E0;
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
        string ionDisWriteTo;
        string twissInput; 
    };
    RingIonEffPara *ringIonEffPara = new RingIonEffPara;
    
   
    // 6) bunch-by-bunch feedbcak
    struct RingBBFB
    {    
        double kickStrengthK[3];
        double gain;
        int delay;
        int taps;
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
    
    
    // 7 wake
    struct RingWake
    {
        string wakeInput;
        string pipeGeoInput;
        string bbrInput;
        int bbrFlag;
        int rwFlag;
        int nTurnswakeTrunction;
    };    
    RingWake * ringWake =  new RingWake;    
    
    // 8 impedance 
    struct RingImpedance
    {
        int bunchBinNumberZ;
        string impedInput;
    };    
    RingImpedance * ringImpedance =  new RingImpedance;    
        
    // 9 run
    struct RingRun
    {
        int calSetting;
        int synRadDampingFlag;
        int fIRBunchByBunchFeedbackFlag;
        int beamIonFlag;
        int nTurns;
        int impedanceFlag;
        int wakeFlag;
        string TBTBunchData;
        int bunchInfoPrintInterval;
    };       
    RingRun * ringRun =  new RingRun; 
    
};

#endif
