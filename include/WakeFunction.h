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

class WakeFunction
{

public:
    WakeFunction();
    ~WakeFunction();
    
    vector<vector<double> > posxData;
    vector<vector<double> > posyData;
    vector<vector<double> > poszData;
    

	double betaFunAver[3]; 
    
	vector<double> sectorRadiusX;
	vector<double> sectorRadiusY;
	vector<double> sectorLength;
	vector<double> sectorBetaX; 
	vector<double> sectorBetaY;
	vector<double> sectorNum; 
	vector<double> sectormatSigma;        
    vector<vector<double> > yokoyaFactor;
	vector<double> GetRWWakeForce (double tau);


    // definition of the bbr wake and impedance can be found in elegant mannul paper followed Refer to: Ref. NIMA 221-230 806 (2016) Nagaoka
    vector<double> lRs;     //ohm
    vector<double> lQ;
    vector<double> lOmega;
    vector<double> tRs;     //ohm /m
    vector<double> tQ;
    vector<double> tOmega;
	vector<double> GetBBRWakeForce(double tau);

    void Initial(ReadInputSettings &inputParameter);
    
    
    
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
