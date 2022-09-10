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



using namespace std;
using std::vector;
using std::complex;


class SPBeam 
{
public:
    SPBeam();
    ~SPBeam();
    
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
    vector<vector<double> > coupledBunchModeAmpX;
    vector<vector<double> > coupledBunchModeAmpY;
    vector<vector<double> > coupledBunchModeAmpZ; 
    vector<vector<double> > hilbertAmpX;
    vector<vector<double> > hilbertAmpY;
    vector<vector<double> > hilbertAmpZ; 
    vector<vector<double> > historyAverX;
    vector<vector<double> > historyAverY;
    vector<vector<double> > historyAverZ;



    void Run(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter, CavityResonator &cavityResonator);
    void SPBeamRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void SPGetBeamInfo();
    void SPBeamDataPrintPerTurn(int turns, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);
    void WSBeamIonEffectOneInteractionPoint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k);
    void GetBinDistBetweenBunch(ReadInputSettings &inputParameter);   
    void GetHaissinski(ReadInputSettings &inputParameter,CavityResonator &cavityResonator,WakeFunction &sRWakeFunction);
    void GetAnalyticalLongitudinalPhaseSpace(ReadInputSettings &inputParameter,CavityResonator &cavityResonator,WakeFunction &sRWakeFunction);
    void GetTimeDisToNextBunchIntial(ReadInputSettings &inputParameter);

    // shared funcitons by MP and SP cases.         
    void Initial(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);
    void InitialcavityResonator(ReadInputSettings &inputParameter,CavityResonator &cavityResonator);    
    void BeamTransferPerInteractionPointDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BeamTransferPerTurnDueToLatticeL(ReadInputSettings &inputParameter,CavityResonator &cavityResonator,int turns);
    void BeamTransferPerTurnDueToLattice(LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter,CavityResonator &cavityResonator,int turns);
    void BeamTransferPerTurnDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint);
    void WSIonDataPrint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,int count);   
    void SetBeamPosHistoryDataWithinWindow();
    vector<complex<double> > GetHilbertAnalytical(vector<double> signal);
    // vector<complex<double> > GetHilbertAnalytical();
    void GetHilbertAnalyticalInOneTurn();                  

    void BeamSynRadDamping(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint);
    void FIRBunchByBunchFeedback(FIRFeedBack &firFeedBack,int nTurns);
    void BeamTransferPerTurnDueWake();
    // for long range RW wake function
	void LRWakeBeamIntaction(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction, const  LatticeInterActionPoint &latticeInterActionPoint,int turns);  
    

      		
private:
};




#endif
