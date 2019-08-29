#ifndef BUNCH_H
#define BUNCH_H

#include <vector>
#include <complex>
#include "Global.h"
#include "LatticeInterActionPoint.h"
#include "Train.h"
#include "ReadInputSettings.h"

using namespace std;
using std::vector;
using std::complex;


class Bunch
{


public:
    Bunch();
    ~Bunch();

    double CorssSectionEI;
    double pipeAperatureR;
    double PipeAperatureX;
    double PipeAperatureY;

    int    distributionType;
    int    macroEleNumPerBunch;
    int    macroEleNumSurivePerBunch;
    double macroEleCharge;
    double electronNumPerBunch;
    double current;             // [A]
    double electronEnergy;              // [eV]
    double rGamma;
    double rBeta;
 
    vector<double > ePositionX; // [m]
    vector<double > ePositionY; // [m]
    vector<double > ePositionZ; // [m]
    vector<double > eMomentumX; // [rad]  -- BeamTransferDueToIon(), px(i)/pz 
    vector<double > eMomentumY; // [rad]  -- BeamTransferDueToIon(), py(i)/pz
    vector<double > eMomentumZ; // [rad]  -- BeamTransferDueToIon(), pz(i)/pz
    vector<double > eFx;        // [rad]
    vector<double > eFy;        // [rad]
    vector<int    >  eSurive;
    
    double actionJx;    // used for weak-strong simulation
    double actionJy;    // used for weak-strong simulation

    

    // use in single-bunch-board-band-impedance calculation in fre domain---------
    vector<double > beamCurDenZProf;   // store the   rho(z) as funciton of z.   <I(z)>
    vector<double > beamCurDenTProf;   // store the y*rho(z) as function of z.   <y I(z)>
    double bunchBinNumberZ;
    double bunchBinNumberT;
    void InitialBeamCurDenZProf();
    void InitialBeamCurDenTProf();
    //-----------------------------------------------------------------------------

    
    void Initial(LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);


    double emittanceX;          // [m rad] non-normalized rms emittance of beam in one bunch
    double emittanceY;          // [m rad] non-normalized rms emittance of beam in one bunch
    double emittanceZ; 			// [m rad] non-normalized rms emittance of beam in one bunch
    double kappa;               // emittance coupling factor   emittanceYNorm/emittanceXNorm <1 in light source
    void   DistriGenerator(LatticeInterActionPoint &latticeInterActionPoint, int i);
    double GSSlover(double if0);
    
    void InducedIonDensity(LatticeInterActionPoint &latticeInterActionPoint); // ion density induced by per bunch
    int    bunchGap;            // the number of rf period needed for the coming bunch 
    
    double rmsRx;
    double rmsRy;
    double rmsEmitX=0.E0;
    double rmsEmitY=0.E0;
    double xAver=0.E0;
    double yAver=0.E0;
    double zAver=0.E0;
    double pxAver=0.E0;
    double pyAver=0.E0;
    double pzAver=0.E0;
    double rmsRingEmitX=0.E0;       
    double rmsRingEmitY=0.E0;       
    double rmsRingEmitZ=0.E0;       
    double initialDisDx;
    double initialDisDy;
    
    double rmsBunchLength=0.E0;     
    double rmsEnergySpread=0.E0;    //rad
    double omegaE=0.E0;				// Ohmi's paper Eq.10
    double omegaI=0.E0;				// Ohmi's paper Eq.10
    void RMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void WSRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k);
    
    void SSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void WSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BassettiErskine(double posx,double posy,double rmsRxTemp, double rmsRyTemp,double &tempFx,double &tempFy);
    
    void BunchTransferDueToLattice(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BunchTransferDueToIon(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BunchSynRadDamping(vector<double> &synchRadDampTime,LatticeInterActionPoint &latticeInterActionPoint);

    
private:

};




#endif
