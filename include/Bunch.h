//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************

#ifndef BUNCH_H
#define BUNCH_H

#include <vector>
#include <complex>
#include "Global.h"
#include "LatticeInterActionPoint.h"
#include "Train.h"
#include "ReadInputSettings.h"
#include "CavityResonator.h"
#include "WakeFunction.h"
#include "Resonator.h"
#include "Spline.h"


using std::vector;
using std::complex;
using v1i= vector<int> ;
using v1d= vector<double> ;
using v2d= vector<vector<double> > ;
using v3d= vector<vector<vector<double> > > ;

class Bunch
{


public:
    Bunch();
    ~Bunch();

    int currentTurnNum = 0;    
    int    bunchGap;                // the number of rf period needed for the coming bunch
    int    bunchHarmNum;    

    vector<double> ePositionX;      // m
    vector<double> ePositionY;
    vector<double> ePositionZ;    
    vector<double> eMomentumX;      // rad
    vector<double> eMomentumY;
    vector<double> eMomentumZ;
    vector<double> eFxDueToIon;     // rad  
    vector<double> eFyDueToIon;
    vector<double> eFzDueToIon;
    vector<int> eSurive;         // if not Surive--throw out loss infomation
    vector<vector<double> > accPhaseAdvX;    // phaseAdvX[np][3], 0 1
    vector<vector<double> > accPhaseAdvY;    // accumulated phase advance of each particle for tune-spread simulation
    vector<vector<double> > accPhaseAdvZ;    // accumulated phase advance of each particle for tune-spread simulation       

    double xAver=0.E0;
    double yAver=0.E0;
    double zAver=0.E0;
    double pxAver=0.E0;
    double pyAver=0.E0;
    double pzAver=0.E0;
    double zAverLastTurn=0.E0;
    double pzAverLastTurn=0.E0;
    double transmission=1.E0;

    // for beam loading to get the time distance between bunches
    double zMinCurrentTurn =0.E0;
    double zMaxCurrentTurn =0.E0;
    
    complex<double> xAverAnalytical;
    complex<double> yAverAnalytical;
    complex<double> zAverAnalytical;
    // double xAverHilbertAnalytical=0.E0;
    // double yAverHilbertAnalytical=0.E0;
    // double zAverHilbertAnalytical=0.E0;
    // double pxAverHilbertAnalytical=0.E0;
    // double pyAverHilbertAnalytical=0.E0;
    // double pzAverHilbertAnalytical=0.E0;
    
    int normCurrent = 1;
    int macroEleNumPerBunch=1;
    double macroEleCharge;               // [C]
    double electronNumPerBunch;
    double current;                      // [A]
    double electronEnergy;               // [eV]
    double rGamma;
    double rBeta;
    double timeFromCurrnetBunchToNextBunch;    //[s]
    double timeFromLastBunchToCurrentBunch;
    double timeToNextBunch;
    double timeToLastBunch;
    
    double rmsBunchLength;
    double rmsBunchLengthLastTurn = 0.E0;
    double rmsEnergySpread;
    double rmsRx;
    double rmsRy;
    double emittanceX;      // rms emittance projected
    double emittanceY;      // rms emittance
    double emittanceZ;
    double totIonCharge;
    double eigenEmitX,eigenEmitY;
    double emitXY4D=0;
    double xyCouplingAlpha=0;  //rms value in x-y space
    int latticeSetionPassedCount = 0;


    struct BunchRFModeInfo{
        vector<complex<double> > induceVolBunchCen;
        vector<complex<double> > genVolBunchAver;
        vector<complex<double> > selfLossVolBunchCen;         
        vector<complex<double> > cavVolBunchCen;
        vector<complex<double> > genIgBunchAver;
        vector<double>           genPower;
        vector<double>           beamPower;
        vector<double>           cavPower;
        vector<double>           refPower;
    };
    BunchRFModeInfo *bunchRFModeInfo = new BunchRFModeInfo;
        
    struct Haissinski{
    double dt = 1.e-12 ;      // set as a default value.. 1 ps 
    double dz = dt * CLight;      
    int    nz  ;   
    double zMax;
    double zMin;
    double averZ;
    double rmsZ;
    double averNus = 0;
    double rmsNus = 0;
    vector<double> nus;
    vector<double> actionJ;
    vector<double> bunchProfile;                // head to tail -> [+,-]--ensure to have a positive hamiltonian--SYL Eq. 3.36
    vector<double> bunchPosZ;
    vector<double> totWakePoten;
    vector<double> rwWakePoten;
    vector<double> bbrWakePoten;            
    vector<double> wakeHamiltonian;          
    vector<double> rfHamiltonian;           
    vector<double> totHamiltonian;
    vector<double> vRF;

    vector<double> cavAmp;
    vector<double> cavPhase;                                  
    };
    Haissinski *haissinski = new Haissinski;
    
    
    vector<double> lRWakeForceAver;        // used in the long range wakefunction simulation
    
    vector<vector<double>> xyzHistoryDataToFit;   

    
    void Initial(const ReadInputSettings &inputParameter);
    void BassettiErskine1(double posx,double posy,double rmsRxTemp, double rmsRyTemp,double &tempFx,double &tempFy);
    void GaussianField(double posx,double posy,double rmsRxTemp, double rmsRyTemp,double &tempFx,double &tempFy);
    void BunchTransferDueToIon(const LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BunchTransferDueToLatticeT(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BunchTransferDueToLatticeTSymplectic(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BunchTransferDueToLatticeOneTurnT66(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint);
	void BunchTransferDuetoSkewQuad(const ReadInputSettings &inputParameter);
    void InitialAccumPhaseAdV(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter );
    
    // void BunchTransferDueToLatticeOneTurnT66GPU(const ReadInputSettings &inputParameter, LatticeInterActionPoint &latticeInterActionPoint);
    void BunchLongPosTransferOneTurn(const ReadInputSettings &inputParameter);
    void SetBunchPosHistoryDataWithinWindow();
    void MarkLostParticle(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint);
    void BunchSynRadDamping(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint);
    void BunchTransferDueToWake();
    void BunchTransferDueToDriveMode(const ReadInputSettings &inputParameter, const int n);
    void GetLongiKickDueToCavFB(const ReadInputSettings &inputParameter,Resonator &resonator);
    void BunchMomentumUpdateDueToRFCA(const ReadInputSettings &inputParameter,Resonator &resonator,int resIndex);   
    void BunchEnergyLossOneTurn(const ReadInputSettings &inputParameter);

    // bunch haissinski solution//deal the data in Haissinski structure. 
    void GetBunchHaissinski(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,WakeFunction &sRWakeFunction);
    void GetWakeHamiltonian(const ReadInputSettings &inputParameter,WakeFunction &sRWakeFunction);
    void GetWakeHamiltonianFromRW(const ReadInputSettings &inputParameter,WakeFunction &sRWakeFunction);
    void GetWakeHamiltonianFromBBR(const ReadInputSettings &inputParameter,WakeFunction &sRWakeFunction);   
    void GetRFHamiltonian(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator);
    void GetTotHamiltonian();
    vector<double> GetProfile(const ReadInputSettings &inputParameter);
    void GetBunchAverAndBunchLengthFromHaissinskiSolution();  

    // once got the hassinski--then get the particle longitudinal phase space 
    void GetParticleLongitudinalPhaseSpace(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,int bunchIndex);
    void GetParticleLongitudinalPhaseSpace1(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,int bunchIndex);
    vector<double> LeapFrog(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,vector<double> z,const tk::spline &wakePotenFit);

    void GetAccumuPhaseAdv(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter);

private:

};




#endif
