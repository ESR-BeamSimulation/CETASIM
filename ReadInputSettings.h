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

    int ParamRead();


// read from initial file;

    double circRing               = 1360.4E0   ;                 
    double workQx                 = 114.05E0 ;
    double workQy                 = 106.231E0;
    double workQz                 = 0.00521 ;
    double rfBaseFrequency        = 166.6E6 ;


    double corssSectionEI         = 2.0e-22;                     //  [m^2] Collision crossing section of CO
    double pipeAperatureX         = 1.1E-3;                      // [m] 
    double pipeAperatureY         = 1.1E-3;                      // [m] 
    double pipeAperatureR         = 1.1E-3;                      // [m] 
    double electronBeamEnergy     = 6.0E+9;
    
    int macroIonNumberGeneratedPerIP = 0;

    int trainNumber;
    int totBunchNumber;
    int printInterval;
    vector<double> synchRadDampTime;                            // in the unit of turns
    int macroEleNumPerBunch;
    int bunchBinNumberZ;
    double current;
    
    double emittanceX      = 34.2E-10;
    double kappa           = 0.1;
    double emittanceY      = emittanceX * kappa;
    double rmsBunchLength  = 30.E-3;        
    double rmsEnergySpread = 1.E-3;
    int nTurns;
    int calSettings;
    int    ionMaxNumber           = 1E+5;                        // max limit number of ion generated
    int    numberofIonBeamInterPoint = 1;                           // number of interaction point considered
    int    distributionType;
    double initialDisDx = 0;
    double initialDisDy = 0;

//  file read end    

    double  omegas=0.E0;                // [Hz] sychronous revolution frequency 
    double  t0=0.E0;                    // [s]  EvolutionPeriod  T0
    int     harmonics=0;                // []   number of harmonics 
    int     synRadDampingFlag=0;         // flag to turn on the SynRadDamping
    int     intevalofTurnsIonDataPrint = 0;
    int     fIRBunchByBunchFeedbackFlag = 0;
    
    

};

#endif
