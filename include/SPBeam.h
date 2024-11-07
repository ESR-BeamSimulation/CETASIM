//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#ifndef SPBEAM_H
#define SPBEAM_H

#include "Global.h"
#include "SPBunch.h"
#include "Train.h"
#include "LatticeInterActionPoint.h"
#include "ReadInputSettings.h"
#include "WakeFunction.h"
#include "FIRFeedBack.h"
#include <fstream>
#include <vector>
#include <complex>
#include "CavityResonator.h" 
#include <iostream>
#include <iomanip>


using namespace std;
using std::vector;
using std::complex;


class SPBeam 
{
public:
    SPBeam();
    ~SPBeam();
    int currentTurnNum = 0;
    // rms of all bunches...
    struct WeakStrongBeamInfo{
        
        double actionJxMax=0;
        double actionJyMax=0;        
        double bunchAverXMax=0;
        double bunchAverYMax=0;	    
        double bunchEffectiveSizeXMax=0;
        double bunchEffectiveSizeYMax=0;
	    
	    double bunchAverX=0;
	    double bunchAverY=0;
	    double bunchAverZ=0;
	    double bunchAverPX=0;
	    double bunchAverPY=0;
	    double bunchAverPZ=0;
	    double bunchRmsSizeX=0;
	    double bunchRmsSizeY=0;
	    double bunchRmsSizeZ=0;
	    double bunchRmsSizePX=0;
	    double bunchRmsSizePY=0;
	    double bunchRmsSizePZ=0;

    	vector<vector<double> > beamInfoOneTurn;   
    };
    WeakStrongBeamInfo *weakStrongBeamInfo = new WeakStrongBeamInfo;    
    
    vector<SPBunch> beamVec;
    


    // for coupled bunch mode or bunch-by-bunch grwoth rate calculation
    // nominal method to get the coupled bunch grwothe rate
    vector<vector<double> > coupledBunchModeAmpX;
    vector<vector<double> > coupledBunchModeAmpY;
    vector<vector<double> > coupledBunchModeAmpZ;
    vector<vector<double> > coupledBunchModeArgX;
    vector<vector<double> > coupledBunchModeArgY;
    vector<vector<double> > coupledBunchModeArgZ;
    vector<double>CBMIdealGRX;
    vector<double>CBMIdealGRY;
    vector<double>CBMIdealGRZ;
    vector<double>CBMHTGRX;
    vector<double>CBMHTGRY;
    vector<double>CBMHTGRZ;
    vector<double>CBMIQGRX;
    vector<double>CBMIQGRY;
    vector<double>CBMIQGRZ;
    
    // IQ decomposition to get the coupled bunch grwothe rate  -- only the unstable mode can be identified
    vector<vector<double > > ampXIQ;
    vector<vector<double > > ampYIQ;
    vector<vector<double > > ampZIQ;
    vector<vector<double > > argXIQ;
    vector<vector<double > > argYIQ;
    vector<vector<double > > argZIQ;

    // analytical signal along the tracking turns.
    vector<vector<complex<double> > > anaSignalAverX;
    vector<vector<complex<double> > > anaSignalAverY;
    vector<vector<complex<double> > > anaSignalAverZ;
    //--------------------------------------------------------------------------------

    vector<vector<double> > hilbertCoupledBunchModeAmpX;
    vector<vector<double> > hilbertCoupledBunchModeAmpY;
    vector<vector<double> > hilbertCoupledBunchModeAmpZ; 
    vector<vector<double> > hilbertCoupledBunchModeArgX;
    vector<vector<double> > hilbertCoupledBunchModeArgY;
    vector<vector<double> > hilbertCoupledBunchModeArgZ;



    vector<vector<double> > hilbertAmpX;
    vector<vector<double> > hilbertAmpY;
    vector<vector<double> > hilbertAmpZ; 

    vector<vector<double> > historyAverX;
    vector<vector<double> > historyAverY;
    vector<vector<double> > historyAverZ;
    vector<vector<double> > historyAverPX;
    vector<vector<double> > historyAverPY;
    vector<vector<double> > historyAverPZ;
    
    // for excitation print 
    vector<double > freXIQDecompScan;
    vector<double > freYIQDecompScan;
    vector<double > freZIQDecompScan;
    vector<double> ampIQ;
    vector<double> phaseIQ;



    void Run(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter, CavityResonator &cavityResonator);
    void SPBeamRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void SPGetBeamInfo();
    void SPBeamDataPrintPerTurn(int turns, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter,CavityResonator &cavityResonator);
    void WSBeamIonEffectOneInteractionPoint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k);
    void GetHaissinski(ReadInputSettings &inputParameter,CavityResonator &cavityResonator,WakeFunction &sRWakeFunction);
    void GetAnalyticalLongitudinalPhaseSpace(ReadInputSettings &inputParameter,CavityResonator &cavityResonator,WakeFunction &sRWakeFunction);
    void GetTimeDisToNextBunch(ReadInputSettings &inputParameter);
    void GetDriveModeGrowthRate(const int n, const ReadInputSettings &inputParameter);
    void GetCBMGR(const int turns, const LatticeInterActionPoint &latticeInterActionPoint, const ReadInputSettings &inputParameter);
    
    void GetCBMGR1(const int turns, const LatticeInterActionPoint &latticeInterActionPoint, const ReadInputSettings &inputParameter);
    void SetIdealCoupledBunchModeData(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter);
    void SetHilbertCoupledBunchModeData(const LatticeInterActionPoint &latticeInterActionPoint, const ReadInputSettings &inputParameter);
    void SetIQCoupledBunchModeData(const LatticeInterActionPoint &latticeInterActionPoint, const ReadInputSettings &inputParameter);

    // shared funcitons by MP and SP cases.         
    void TuneRamping(ReadInputSettings &inputParameter,double n);
    void Initial(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);
    void InitialcavityResonator(ReadInputSettings &inputParameter,CavityResonator &cavityResonator);    
    void BeamTransferPerInteractionPointDueToLatticeT(const ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BeamMomentumUpdateDueToRF(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,CavityResonator &cavityResonator,int turns);
    void BeamMomentumUpdateDueToRFTest(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,CavityResonator &cavityResonator,int turns);
    void BeamLongiPosTransferOneTurn(const ReadInputSettings &inputParameter); 
    void BeamEnergyLossOneTurn(const ReadInputSettings &inputParameter);
    void PrintHaissinski(ReadInputSettings &inputParameter);
    void BeamTransferDueToSkewQuad(const ReadInputSettings &inputParameter);
    

    // void BeamTransferPerTurnDueToLattice(LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter,CavityResonator &cavityResonator,int turns);
    // void BeamTransferPerTurnDueToLatticeT(const ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint);
    void WSIonDataPrint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,int count);   
    void SetBeamPosHistoryDataWithinWindow();
    void GetAnalyticalWithFilter(const ReadInputSettings &inputParameter);
    vector<complex<double> > GetHilbertAnalytical(vector<double> signal, const double filterBandWithdNu,  double workQ);
    vector<complex<double> > GetHilbertAnalytical(vector<complex<double> >  signal, const double filterBandWithdNu,  double workQ);
    
    void GetHilbertAnalyticalInOneTurn(const ReadInputSettings &inputParameter);
    void BeamTransferDuetoDriveMode(const ReadInputSettings &inputParameter,const int n);                 
    void MarkParticleLostInBunch(const ReadInputSettings &inputParameter, const LatticeInterActionPoint &latticeInterActionPoint);

    void BeamSynRadDamping(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint);
    void FIRBunchByBunchFeedback(const ReadInputSettings &inputParameter,FIRFeedBack &firFeedBack,int nTurns);
    void BeamTransferPerTurnDueWake();
    void BeamTransferPerTurnDueToLatticeTOneTurnR66(const ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint);
    

    // for long range RW wake function
	void LRWakeBeamIntaction(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction, const  LatticeInterActionPoint &latticeInterActionPoint,int turns);  
    

      		
private:
};




#endif
