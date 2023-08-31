//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#ifndef WakeFunction_H
#define WakeFunction_H

#include <cmath>
#include <complex>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <vector>
#include <string>
#include <algorithm>
#include "Global.h"
#include "ReadInputSettings.h"
#include "LatticeInterActionPoint.h"
class WakeFunction
{

public:
    WakeFunction();
    ~WakeFunction();
    
    vector<vector<double> > posxData;
    vector<vector<double> > posyData;
    vector<vector<double> > poszData;
    

	double betaFunIntPoint[2];       // x y
    double betaFunAver[2];
    
	vector<double> sectorRadiusX;
	vector<double> sectorRadiusY;
	vector<double> sectorLength;
	vector<double> sectorBetaX; 
	vector<double> sectorBetaY;
	vector<double> sectorNum; 
	vector<double> sectormatSigma;        
    vector<vector<double> > yokoyaFactor;

	vector<double> GetRWSRWakeFun (double tau) ; // Alex Chao notation -- full interaction of quasi-green function (3.11) 
    vector<double> GetRWLRWakeFun (double tau) ; // Alex Chao notation -- Eq. 2.53 
    double GetRWLongWakeTerm2(double u, double tau0); 
    vector<double> GetRWPusdoWakeFun(double u) ;                            // Eq. 3.11 and 3.56
    vector<double> GetBBRPusdoWakeFun(double u, int i,double sigmat);       // Eq. 3.14 and 3.58

    // definition of the bbr wake and impedance can be also found in elegant mannul.   
    //Refer to: Ref. NIMA 221-230 806 (2016) Nagaoka for the long range wake wake function estimation. 
    
    vector<double> lRs;     //ohm
    vector<double> lQ;
    vector<double> lOmega;
    vector<double> tyRs;     //ohm /m
    vector<double> tyQ;
    vector<double> tyOmega;
    vector<double> txRs;     //ohm /m
    vector<double> txQ;
    vector<double> txOmega;
    

	vector<double> GetBBRWakeFun(double tau); // works for both long and short range wakefunction
    vector<double> GetBBRWakeFun1(double tau) ; // works for both long and short range wakefunction--1mm bunch pusedo-wake potential as wake function

    void InitialLRWake(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint);
    void InitialSRWake(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint);
    void BBRWakeParaReadIn(string inputfilename);
    void RWWakeParaReadIn(string inputfilename);
            
	void GetyokoyaFactor(double a,double b,vector<double> &yokoyaFactor);
	double Lfunction1 (double r,double t,double mu);
	double Lfunction2 (double r,double t,double mu);
	double Lfunction  (double r,double t,double mu);
	double Ldfunction1(double r,double t,double mu);
	double Ldfunction2(double r,double t,double mu);
	double Ldxfunction(double r,double t,double mu);
	double Ldyfunction(double r,double t,double mu);	

};


#endif
