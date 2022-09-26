//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#ifndef SBBUNCH_H
#define SBBUNCH_H

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
using namespace std;
using std::vector;
using std::complex;

class SPBunch : public Bunch
{


public:
    SPBunch();
    ~SPBunch();

    
    double actionJx=0.E0;        // used for weak-strong simulation
    double actionJy=0.E0;        // used for weak-strong simulation
  
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
    
    

    void InitialSPBunch(const ReadInputSettings &inputParameter);        
    void DistriGenerator(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter, int i);
    void GetSPBunchRMS(const LatticeInterActionPoint &latticeInterActionPoint, int k);    
    void WSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BunchTransferDueToLatticeL(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator,int turns);
    void BunchTransferDueToLatticeLTest(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator,int turns);
    void BunchTransferDueToLatticeLMatarix(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator,int turns);
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
