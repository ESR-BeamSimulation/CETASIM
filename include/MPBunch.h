//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#ifndef MPBUNCH_H
#define MPBUNCH_H

#include <vector>
#include <complex>
#include "Global.h"
#include "PIC3D.h"
#include "LatticeInterActionPoint.h"
#include "Train.h"
#include "ReadInputSettings.h"
#include "CavityResonator.h"
#include "Bunch.h"
#include "WakeFunction.h"
#include "Spline.h"
#include "BoardBandImp.h"
#include "BeamIon2DPIC.h"


using namespace std;
using std::vector;
using std::complex;

using v1d= vector<double> ;
using v2d= vector<vector<double> > ;
using v3d= vector<vector<vector<double> > > ;


class MPBunch : public Bunch
{


public:
    MPBunch();
    ~MPBunch();

    int macroEleNumSurivePerBunch;
    double rmsEffectiveRingEmitX,rmsEffectiveRingEmitY,rmsEffectiveRingEmitZ;
    double rmsEffectiveRx,rmsEffectiveRy;
    
    // used in the soft bunch transient beam loading simulation-- bin-by-bin process.
    vector<double > beamCurDenZProf;   // store the   rho(z) as funciton of z.   <I(z)>
    vector<vector<complex<double> > >  cavVolInfoVsLongBins;
    vector<vector<double > > cavForceInfoVsLongBins;
    vector<double> posZBins;
    vector<double> densProfVsBin;
    vector<double> hamiltonPotenWell;
    vector<double> densProfVsBinAnalytical;
    double zAverAnalytical;
    double bunchLengthAnalytical;
    vector<vector<double>> srWakePoten;
    

    // the wakepoten here decleared for solver in frequecy domain
    struct WakePotenFromBBI
    {
        vector<double> wakePotenZ;
        vector<double> wakePotenDx;
        vector<double> wakePotenDy;
        vector<double> wakePotenQx;
        vector<double> wakePotenQy;     
    };    
    WakePotenFromBBI *wakePotenFromBBI =  new WakePotenFromBBI;
        
    vector<double> profileForBunchBBImp;
    

    // 2.5d PIC model for space charge simulation
    vector<vector<vector<double> > > slicedBunchDisInfo ;  // instore beam along for 2.5D PIC model (-5 simgaz, 5*simgaz) 11 slices as default.   
    vector<vector<double> > slicedBunchChargeInfo;
    vector<vector<int   > > slicedBunchPartiIndex;
         
    void InitialMPBunch(const ReadInputSettings &inputParameter);        
    void DistriGenerator(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter, int i);
    double GSSlover(const double if0);
    void GetMPBunchRMS(const LatticeInterActionPoint &latticeInterActionPoint, int k);    
    void SSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void SSIonBunchInteractionPIC(BeamIon2DPIC &beamIon2DPIC,LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BunchTransferDueToSRWake(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction, const LatticeInterActionPoint &latticeInterActionPoint, int turns);
    void GetZMinMax();
    void BBImpBunchInteraction(const ReadInputSettings &inputParameter, const BoardBandImp &boardBandImp, const LatticeInterActionPoint &latticeInterActionPoint);
    void BBImpBunchInteractionTD(const ReadInputSettings &inputParameter,const BoardBandImp &boardBandImp,const LatticeInterActionPoint &latticeInterActionPoint,vector<vector<double>> wakePoten);
    void GetSmoothedBunchProfileGassionFilter(const ReadInputSettings &inputParameter, const BoardBandImp &boardBandImp );
    void GetSmoothedBunchProfileGassionFilter(double *profile,int nBins);
    void GetEigenEmit(const LatticeInterActionPoint &latticeInterActionPoint);
    void BunchMomentumUpdateDueToSpaceChargePIC(PIC3D &picSCBeam3D,LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BunchMomentumUpdateDueToSpaceChargeAnalytical(LatticeInterActionPoint &latticeInterActionPoint, int k, const ReadInputSettings &inputParameter);
    void SetSlicedBunchInfo(const int nz);
    void UpdateTransverseMomentumBEModel(vector<vector<double>> &particles, vector<double> eCharge, const double ds);
    void UpdateTransverseMomentumLinear(vector<vector<double>> &particles, vector<double> eCharge, const double ds);

    void BunchMomentumUpdateDueToRFMode(const ReadInputSettings &inputParameter,Resonator &resonator, int resIndex);
    void BunchMomentumUpdateDueToRFModeStable(const ReadInputSettings &inputParameter,Resonator &resonator, int resIndex);      
    void BunchMomentumUpdateDuetoRFRigid(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator);

    // BeamInduced Voltage calculated once per bunch. Bunch distance is from center to center.
    // beaminduced voltage rotate and decay once per bunch
    

    void BunchMomentumUpdateDuetoRFBinByBin(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator);
    // BeamInduced Voltage calculated bin-by-bin per bunch. Bunch distance is from zMax of ith bunch to zMin of (i+1)th bunch.
    
    void BunchTransferDueToLatticeLNoInstability(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator);
    // very first version, now used now. 
private:

};


#endif
