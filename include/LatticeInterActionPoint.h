//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#ifndef LATTICEINTERACTIONPOINT_H
#define LATTICEINTERACTIONPOINT_H


#include "Global.h"
#include <vector>
#include "ReadInputSettings.h"

using namespace std;
using std::vector;
using std::complex;


class LatticeInterActionPoint
{


public:
    LatticeInterActionPoint();
    ~LatticeInterActionPoint();

    int numberOfInteraction;
    double circRing;
    vector<double> ionMassNumber;
    vector<double> corssSectionEI;
    vector<double> gasPercent;
    int gasSpec;
     
    int harmonics;
    

    // (1) lattice function
    vector<double> twissAlphaX;
    vector<double> twissBetaX;              // m/rad  
    vector<double> twissAlphaY;
    vector<double> twissBetaY;              // m/rad
    vector<double> twissAlphaZ;
    vector<double> twissBetaZ;              // m/rad
    
    vector<double> twissDispX;				// m
    vector<double> twissDispPX;             // rad
    vector<double> twissDispY;				// m
    vector<double> twissDispPY;             // rad    
   
    vector<vector<double> > transferMatrix;
    vector<vector<double> > xTransferMatrix;
    vector<vector<double> > yTransferMatrix;
    vector<vector<double> > zTransferMatrix;
    vector<double> xPhaseAdv;
    vector<double> yPhaseAdv;
    vector<double> zPhaseAdv;
    
    vector<double>vacuumPressure;                          // Pa
    vector<vector<double> >vacuumPressureEachGas;          // Pa
    vector<double> temperature;                            // K
    vector<vector<double> >ionLineDensity;                 // [1/m] ion number per meter in rings. 
    vector<double> interactionLength;                      // [m]   effective length of interaction point \sum interactionLength[k] = circ

    vector<double> pipeAperatureX;
    vector<double> pipeAperatureY;


    //(2) ion information definition
    vector<vector<double>  >ionNumber;                       // pth type ion number to be generted at ith interaction point
    vector<vector<int   >  >macroIonNumber;                  // pth type macro ionnumber  at ith interaction point 
    vector<vector<double>  >macroIonCharge;                  // pth type macroion charge  at ith interaction point  [e]  
    vector<vector<vector<double> > >ionPositionX;            // m    pth type ion info at ith interction point
    vector<vector<vector<double> > >ionPositionY;            // m
    vector<vector<vector<double> > >ionVelocityX;            // m/s
    vector<vector<vector<double> > >ionVelocityY;            // m/s
    

    vector<vector<int> >ionAccumuNumber;
    vector<vector<vector<double> > >ionAccumuPositionX;      // m --  pth type accumulated ions info at ith interction point
    vector<vector<vector<double> > >ionAccumuPositionY;      // m
    vector<vector<vector<double> > >ionAccumuVelocityX;      // m/s
    vector<vector<vector<double> > >ionAccumuVelocityY;      // m/s
    
    vector<vector<double> >ionAccumuAverX;
    vector<vector<double> >ionAccumuAverY;
    vector<vector<double> >ionAccumuAverVelX;
    vector<vector<double> >ionAccumuAverVelY;

    vector<vector<double> >ionAccumuRMSX;                    // p type ions at ith interaction point 
    vector<vector<double> >ionAccumuRMSY; 
    
    vector<vector<vector<double> > >ionAccumuFx;             
    vector<vector<vector<double> > >ionAccumuFy;             
      
    int ionMaxNumberOneInterPoint;                          // pth type ions at ith interaction point.  Maxiuimum allowed macro ions number
          
    vector<double> totIonChargeAtInterPoint;                // info along the whole ring.
    vector<int>  totMacroIonsAtInterPoint;
    int    totMacroIons;
    double totIonCharge;
    
    double ionLossBoundary;
    double latticeParaForOneTurnMap[22];
    double latticeSynRadBRH[72];
    double latticeSynRadCoff[6];  //Ref ZhangYuan -- used in GPU
    double tengRMatR2[4];         //Ref ZhangYuan Eq.(6)  (0,0),(0,1),(1,0),(1,1)
    double traceAB[3];            //Sagan Eq.(3) trace of A B and M-N   
    double gammaC[2];             //Sagan Eq.(8)
    double detH;      
    struct ResDrivingTerms{
        // skew quadrupole 
        complex<double> f1001;
        complex<double> f1010;
        complex<double> f0101;
        complex<double> f0110;
        double linearCouplingFactor;             // according to the definition of Yong-Chul' note
        double xyAlpha;

    };
    ResDrivingTerms *resDrivingTerms = new ResDrivingTerms;

    void Initial(const ReadInputSettings &inputParameter);
    void InitialLattice(const ReadInputSettings &inputParameter);
    void IonGenerator(double rmsRx, double rmsRy, double xAver,double yAver, int k);
    void IonsUpdate(int k);
    void IonRMSCal(int k);
    void IonRMSCal(int k,int p);
    void IonTransferDueToBunch(int bunchGap, int k, double bunchSizeXMax, double bunchSizeYMax);    // ion loss certeria done at here. 
    void GetTotIonCharge();
    void GetIonNumberPerInterAction(double electronNumPerBunch, int k);
    void SetLatticeParaForOneTurnMap(const ReadInputSettings &inputParameter);
    void SetLatticeBRHForSynRad(const ReadInputSettings &inputParameter);  //Hirata 1997 paper
    void GetTransLinearCouplingCoef(const ReadInputSettings &inputParameter);
    
private:

};




#endif
