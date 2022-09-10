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
    complex<double> xAverAnalytical;
    complex<double> yAverAnalytical;
    complex<double> zAverAnalytical;
    // double xAverHilbertAnalytical=0.E0;
    // double yAverHilbertAnalytical=0.E0;
    // double zAverHilbertAnalytical=0.E0;
    // double pxAverHilbertAnalytical=0.E0;
    // double pyAverHilbertAnalytical=0.E0;
    // double pzAverHilbertAnalytical=0.E0;
    

    int macroEleNumPerBunch=1;
    double macroEleCharge;
    double electronNumPerBunch;
    double current;                      // [A]
    double electronEnergy;               // [eV]
    double rGamma;
    double rBeta;
    double timeToNextBunch;    //[s]
    double timeToLastBunch;
    
    double rmsBunchLength;
    double rmsEnergySpread;
    double rmsRx;
    double rmsRy;
    double emittanceX;      // rms emittance
    double emittanceY;      // rms emittance
    double emittanceZ;
    double totIonCharge;

    struct CAVFBCenInfo{
        
        vector<complex<double> > induceVolBunchCen;
        vector<complex<double> > genVolBunchAver;
        vector<complex<double> > selfLossVolBunchCen;         
        vector<complex<double> > cavVolBunchCen;
        vector<double> cavAmpBunchCen;
        vector<double> cavPhaseBunchCen;
           
    };
    CAVFBCenInfo *cavFBCenInfo = new CAVFBCenInfo;
    
    vector<double> lRWakeForceAver;        // used in the long range wakefunction simulation
    
    vector<vector<double>> xyzHistoryDataToFit;   

    
    void Initial(const ReadInputSettings &inputParameter);
    void BassettiErskine1(double posx,double posy,double rmsRxTemp, double rmsRyTemp,double &tempFx,double &tempFy);
    void GaussianField(double posx,double posy,double rmsRxTemp, double rmsRyTemp,double &tempFx,double &tempFy);
    void BunchTransferDueToIon(const LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BunchTransferDueToLatticeT(const LatticeInterActionPoint &latticeInterActionPoint, int k);
    void SetBunchPosHistoryDataWithinWindow();

    void BunchSynRadDamping(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint);
    void BunchTransferDueToWake();
    


private:

};




#endif
