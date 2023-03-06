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
#include "LatticeInterActionPoint.h"
#include "Train.h"
#include "ReadInputSettings.h"
#include "CavityResonator.h"
#include "Bunch.h"
#include "WakeFunction.h"
#include "Spline.h"
#include "BoardBandImp.h"
using namespace std;
using std::vector;
using std::complex;

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
    vector<double> profileForBunchBBImp;
    vector<double> beamIndVFromBBImpZ;
    vector<double> beamIndVFromBBImpX;
    vector<double> beamIndVFromBBImpY;

         
    void InitialMPBunch(const ReadInputSettings &inputParameter);        
    void DistriGenerator(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter, int i);
    double GSSlover(const double if0);
    void GetMPBunchRMS(const LatticeInterActionPoint &latticeInterActionPoint, int k);    
    void SSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BunchLongiMomentumUpdateDuetoRF(const ReadInputSettings &inputParameter);
    void BunchTransferDueToSRWake(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction, const LatticeInterActionPoint &latticeInterActionPoint, int turns);
    void GetZMinMax();
    void BBImpBunchInteraction(const ReadInputSettings &inputParameter, const BoardBandImp &boardBandImp );
    void GetBunchProfileForBeamBBImpEffect(const ReadInputSettings &inputParameter, const BoardBandImp &boardBandImp );
    

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
