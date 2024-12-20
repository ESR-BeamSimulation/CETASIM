//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#ifndef MPBEAM_H
#define MPBEAM_H

#include "Global.h"
#include "MPBunch.h"
#include "Train.h"
#include "LatticeInterActionPoint.h"
#include "ReadInputSettings.h"
#include "WakeFunction.h"
#include "FIRFeedBack.h"
#include "BoardBandImp.h"
#include "BeamIon2DPIC.h"
#include <fstream>
#include <vector>
#include <complex>
#include "CavityResonator.h" 
#include "PIC3D.h"

using namespace std;
using std::vector;
using std::complex;



class MPBeam 
{
public:
    MPBeam();
    ~MPBeam();
    
    int currentTurnNum = 0;
    struct StrongStrongBunchInfo{

        double emitYMax;
        double emitXMax;
        double emitZMax;
        double effectivEemitXMax;
        double effectivEemitYMax;
        double bunchSizeXMax;
        double bunchSizeYMax;
        double bunchLengthMax;
        double bunchEnergySpreadMax;
        double bunchEffectiveSizeXMax;
        double bunchEffectiveSizeYMax;
    	vector<vector<double> > bunchInfoOneTurn;

        double bunchAverXMax;
        double bunchAverYMax;
        double bunchAverZMax;
        double bunchRmsSizeX=0;
	    double bunchRmsSizeY=0;
	    double bunchRmsSizeZ=0;
	    double bunchRmsSizePX=0;
	    double bunchRmsSizePY=0;
	    double bunchRmsSizePZ=0;
	    double bunchAverX=0;
	    double bunchAverY=0;
	    double bunchAverZ=0;
	    double bunchAverPX=0;
	    double bunchAverPY=0;
	    double bunchAverPZ=0;
        
    };
    StrongStrongBunchInfo *strongStrongBunchInfo = new StrongStrongBunchInfo;      
    vector<MPBunch> beamVec;
    vector<vector<double> > bunchZMinZMax; // bunchTMaxTMinTAver[i][0,1] -> [mic,max]

    struct QuasiWakePoten{
        double bunchLength;
        vector<double> wz;              // V/C
        vector<double> wDx;             // V/C/m
        vector<double> wDy;             // V/C/m
        vector<double> wQx;             // V/C/m    
        vector<double> wQy;             // V/C/m
        vector<double> binPosZ;
    };
    QuasiWakePoten *quasiWakePoten = new QuasiWakePoten;

    // for GPU, all particle cord 
    // double  *partCord; 
    //
    
    // for coupled bunch mode or bunch-by-bunch growth rate calculation
    // nominal method to get the coupled bunch grwothe rate
    vector<vector<double> > coupledBunchModeAmpX;
    vector<vector<double> > coupledBunchModeAmpY;
    vector<vector<double> > coupledBunchModeAmpZ;
    vector<vector<double> > coupledBunchModeArgX;
    vector<vector<double> > coupledBunchModeArgY;
    vector<vector<double> > coupledBunchModeArgZ;
    // IQ decomposition to get the coupled bunch grwothe rate  -- only the unstable mode can be identified
    vector<vector<double > > ampXIQ;
    vector<vector<double > > ampYIQ;
    vector<vector<double > > ampZIQ;
    vector<vector<double > > argXIQ;
    vector<vector<double > > argYIQ;
    vector<vector<double > > argZIQ;
    // history bunch position info
    vector<vector<double> > historyAverX;
    vector<vector<double> > historyAverY;
    vector<vector<double> > historyAverZ;
      // analytical signal along the tracking turns.
    vector<vector<complex<double> > > anaSignalAverX;
    vector<vector<complex<double> > > anaSignalAverY;
    vector<vector<complex<double> > > anaSignalAverZ;
    // hilbert transfer for compled bunch mode calculation
    vector<vector<double> > hilbertCoupledBunchModeAmpX;
    vector<vector<double> > hilbertCoupledBunchModeAmpY;
    vector<vector<double> > hilbertCoupledBunchModeAmpZ; 
    vector<vector<double> > hilbertCoupledBunchModeArgX;
    vector<vector<double> > hilbertCoupledBunchModeArgY;
    vector<vector<double> > hilbertCoupledBunchModeArgZ;


    vector<double > freXIQDecompScan;
    vector<double > freYIQDecompScan;
    vector<double > freZIQDecompScan;
    //--------------------------------------------------------------------------------


    vector<double> ampIQ;
    vector<double> phaseIQ;

    void Run(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter, CavityResonator &cavityResonator);
    void MPBeamRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void MPGetBeamInfo();
    void MPBeamDataPrintPerTurn(int turns, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);
    void SSBeamIonEffectOneInteractionPoint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k, BeamIon2DPIC &beamIon2DPIC);
    void SRWakeBeamIntaction(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction, const  LatticeInterActionPoint &latticeInterActionPoint,int turns);    
    void GetTimeDisToNextBunchIntial(ReadInputSettings &inputParameter);
    void GetBunchMinZMaxZ();
    void GetQuasiWakePoten(const ReadInputSettings &inputParameter, const BoardBandImp &boardBandImp);
    void GetQuasiWakePoten(const ReadInputSettings &inputParameter);
    void LongWakePotenBenchmark(const ReadInputSettings &inputParameter,const BoardBandImp &boardBandImp);
    void TuneRamping(ReadInputSettings &inputParameter, double turns);
    //// shared funcitons by MP and SP--just copy of each other.         
    void Initial(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);
    void InitialcavityResonator(ReadInputSettings &inputParameter,CavityResonator &cavityResonator);    
    void BeamTransferPerInteractionPointDueToLatticeT(const ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BeamMomtumUpdateDueToRF(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,CavityResonator &cavityResonator);
    void BeamMomtumUpdateDueToRFTest(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,CavityResonator &cavityResonator);
    void BeamLongiPosTransferOneTurn(const ReadInputSettings &inputParameter);
    void GetTimeDisToNextBunch(const ReadInputSettings &inputParameter);
    void BeamEnergyLossOneTurn(const ReadInputSettings &inputParameter); 
    void BeamTransferPerTurnDueToLatticeTOneTurnR66(const ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint);
    void BeamTransferDueToSkewQuad(const ReadInputSettings &inputParameter);
    void BeamTransferDueToSpaceChargePIC(PIC3D &scBeam3D,LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BeamTransferDueToSpaceChargeAnalytical(LatticeInterActionPoint &latticeInterActionPoint, int k, ReadInputSettings &inputParameter);
    void BeamTransferPerTurnR66AndSynRadGPU(const ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint);
    void SSIonDataPrint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,int count);  
    void GetAnalyticalLongitudinalPhaseSpace(ReadInputSettings &inputParameter,CavityResonator &cavityResonator,WakeFunction &sRWakeFunction);
    void GetHaissinski(ReadInputSettings &inputParameter,CavityResonator &cavityResonator,WakeFunction &sRWakeFunction);
    void BeamTransferDuetoDriveMode(const ReadInputSettings &inputParameter, const int n);
    void MarkParticleLostInBunch(const ReadInputSettings &inputParameter, const LatticeInterActionPoint &latticeInterActionPoint);
    void GetDriveModeGrowthRate(const int turns, const ReadInputSettings &inputParameter);
    void GetCBMGR(const int turns, const LatticeInterActionPoint &latticeInterActionPoint, const ReadInputSettings &inputParameter);
    void GetAnalyticalWithFilter(const ReadInputSettings &inputParameter); 
    vector<complex<double> > GetHilbertAnalytical(vector<double> signal, const double filterBandWithdNu, double workQ);  
    //void SetBeamPosHistoryDataWithinWindow();

    void MPGetAccumuPhaseAdv(const LatticeInterActionPoint &latticeInterActionPoint,const  ReadInputSettings &inputParameter);    // not used          
                               	                       
    void BeamSynRadDamping(const ReadInputSettings &inputParameter, LatticeInterActionPoint &latticeInterActionPoint);
    void FIRBunchByBunchFeedback(const ReadInputSettings &inputParameter,FIRFeedBack &firFeedBack,int nTurns);
    void BeamTransferPerTurnDueWake();
    //// for long range RW wake function
	void LRWakeBeamIntaction(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction, const  LatticeInterActionPoint &latticeInterActionPoint);  
    void BBImpBeamInteraction(const ReadInputSettings &inputParameter, const BoardBandImp &boardBandImp,const LatticeInterActionPoint &latticeInterActionPoint );
    void CopyPartCordToGPU(double *partCord, int totalPartiNum);
    void CopyPartCordFromGPU(double *partCord);

      		
private:
};




#endif
