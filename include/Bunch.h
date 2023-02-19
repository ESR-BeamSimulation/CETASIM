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
#include "Spline.h"
using namespace std;
using std::vector;
using std::complex;

class Bunch
{


public:
    Bunch();
    ~Bunch();
        
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

    double xAver=0.E0;
    double yAver=0.E0;
    double zAver=0.E0;
    double pxAver=0.E0;
    double pyAver=0.E0;
    double pzAver=0.E0;
    double zAverLastTurn=0.E0;
    double pzAverLastTurn=0.E0;

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
    double macroEleCharge;
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
    double emittanceX;      // rms emittance
    double emittanceY;      // rms emittance
    double emittanceZ;
    double totIonCharge;

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
    vector<double> bunchProfile;                // head to tail -> [+,-]--ensure to have a positive hamiltonian--SYL Eq. 3.36
    vector<double> bunchPosZ;
    vector<double> totWakePoten;
    vector<double> rwWakePoten;
    vector<double> bbrWakePoten;            
    vector<double> wakeHamiltonian;          
    vector<double> rfHamiltonian;           
    vector<double> totHamiltonian;           
    
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
    void BunchTransferDueToLatticeT(const LatticeInterActionPoint &latticeInterActionPoint, int k);
    void SetBunchPosHistoryDataWithinWindow();
    void MarkLostParticle(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint);
    void BunchSynRadDamping(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint);
    void BunchTransferDueToWake();
    void BunchTransferDueToDriveMode(const ReadInputSettings &inputParameter, const int n);
    void GetLongiKickDueToCavFB(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator);
    void GetVBSampledDueToIthBunch(const int j,const int k,const ReadInputSettings &inputParameter, CavityResonator &cavityResonator);

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
    vector<double> LeapFrog(const ReadInputSettings &inputParameter,const CavityResonator &cavityResonator,vector<double> z,const tk::spline &wakePotenFit);


private:

};




#endif
