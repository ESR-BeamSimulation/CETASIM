//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             


#include <vector>
#include <complex>
#include <iostream>
#include<fstream>	
#include <numeric>
#include <cmath>
#include "WakeFunction.h"
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_bessel.h>
#include "Global.h"
#include <iomanip>

using namespace std;
using std::vector;
using std::complex;



WakeFunction::WakeFunction()
{	    
    
}

WakeFunction::~WakeFunction()
{   
}

void WakeFunction::InitialLRWake(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint)
{
    int nTurnswakeTrunction = inputParameter.ringLRWake->nTurnswakeTrunction;   
    int totBunchNum         = inputParameter.ringFillPatt->totBunchNumber;    

    // wake interatin is applied at the first point in twiss.dat file
    betaFunIntPoint[0]  = latticeInterActionPoint.twissBetaX[0];
    betaFunIntPoint[1]  = latticeInterActionPoint.twissBetaY[0];

    betaFunAver[0]      = inputParameter.ringParBasic->betaFunAver[0];  // average beta funtion in x and y
    betaFunAver[1]      = inputParameter.ringParBasic->betaFunAver[1];

    posxData.resize(nTurnswakeTrunction);
    posyData.resize(nTurnswakeTrunction);
    poszData.resize(nTurnswakeTrunction);
    
    for(int i=0;i<nTurnswakeTrunction;i++)
    {
        posxData[i].resize(totBunchNum);
        posyData[i].resize(totBunchNum);
        poszData[i].resize(totBunchNum);
    }

    for(int i=0;i<nTurnswakeTrunction;i++)
    {
        for(int j=0;j<totBunchNum;j++)
        {
            posxData[i][j] = 0.E0;
            posyData[i][j] = 0.E0;
            poszData[i][j] = 0.E0;
        }
    }

    if(!inputParameter.ringLRWake->pipeGeoInput.empty())
    {    
        RWWakeParaReadIn(inputParameter.ringLRWake->pipeGeoInput);
    }
    if(!inputParameter.ringLRWake->bbrInput.empty())
    {    
        BBRWakeParaReadIn(inputParameter.ringLRWake->bbrInput);    
    }
}

void WakeFunction::InitialSRWake(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint)
{
    betaFunIntPoint[0]  = latticeInterActionPoint.twissBetaX[0];
    betaFunIntPoint[1]  = latticeInterActionPoint.twissBetaY[0];
    betaFunAver[0]      = inputParameter.ringParBasic->betaFunAver[0];
    betaFunAver[1]      = inputParameter.ringParBasic->betaFunAver[1];


    if(!inputParameter.ringSRWake->pipeGeoInput.empty())
    {    
        RWWakeParaReadIn(inputParameter.ringSRWake->pipeGeoInput);
    }
    if(!inputParameter.ringSRWake->bbrInput.empty())
    {    
        BBRWakeParaReadIn(inputParameter.ringSRWake->bbrInput);   
    } 
}


void WakeFunction::BBRWakeParaReadIn(string inputfilename)
{
     ifstream fin1(inputfilename);
     
    if (! fin1.is_open())
    {
        cerr<< "Error opening file: "<< inputfilename <<endl; 
        exit (1);
    }
    else
    {        
        vector<string> strVec;
        string         str;
        getline(fin1,str);
        
        int i=0;
        while (!fin1.eof())
        {
            getline(fin1,str);
            if(str.length()==0)  continue;
		                
	        StringSplit2(str,strVec);

	        lRs     .push_back(stod(strVec[0]));
	        lQ      .push_back(stod(strVec[1]));
	        lOmega  .push_back(stod(strVec[2]) * 2 * PI );
	        tRs     .push_back(stod(strVec[3]));
	        tQ      .push_back(stod(strVec[4]));
	        tOmega  .push_back(stod(strVec[5])* 2 * PI);

	        i++;
        }			
    }
    fin1.close();    
 
    
}

void WakeFunction::RWWakeParaReadIn(string filename)
{     
    ifstream fin(filename);        
    if (! fin.is_open())
    {
        cerr<< "Error opening file: "<< filename <<endl; 
        exit (1);
    }
    else
    {   // read in the geo-parameters     
        vector<string> strVec;
        string         str;
        getline(fin,str);  // skip the first line -- the input have to follow the format.
        
        int i=0;
        while (!fin.eof())
        {
            getline(fin,str);
            if(str.length()==0)  continue;
		                
	        StringSplit2(str,strVec);

	        sectorRadiusX .push_back(stod(strVec[0]));
	        sectorRadiusY .push_back(stod(strVec[1]));
	        sectorLength  .push_back(stod(strVec[2]));
	        sectorBetaX   .push_back(stod(strVec[3]));
	        sectorBetaY   .push_back(stod(strVec[4]));
	        sectorNum     .push_back(stod(strVec[5]));
	        sectormatSigma.push_back(stod(strVec[6]));                 
	        i++;           
        }			

        yokoyaFactor.resize(sectorRadiusX.size());
        
        
        for(int i=0;i<yokoyaFactor.size();i++)
        {
            yokoyaFactor[i].resize(5);
        }

        vector<double> yokoyaFactorTemp(5,1.0);
            
        double temp=0;
        double radius;
        
        // get the Yokoya factor of according to the input data
        for(int i=0;i<sectorRadiusX.size();i++)
        {    
            radius = sectorRadiusX[i] >= sectorRadiusY[i] ? sectorRadiusY[i] : sectorRadiusX[i]; 
            
            if(sectorRadiusX[i]>sectorRadiusY[i])
            {
                GetyokoyaFactor(sectorRadiusX[i],sectorRadiusY[i],yokoyaFactorTemp);
            }
            else if(sectorRadiusX[i]<sectorRadiusY[i])
            {
                GetyokoyaFactor(sectorRadiusY[i],sectorRadiusX[i],yokoyaFactorTemp);
                temp=yokoyaFactorTemp[1];
                yokoyaFactorTemp[1] = yokoyaFactorTemp[2];
                yokoyaFactorTemp[2] = temp;
                temp=yokoyaFactorTemp[3];
                yokoyaFactorTemp[3] = yokoyaFactorTemp[4];
                yokoyaFactorTemp[4] = temp;
            }
                
            for(int j=0;j<yokoyaFactorTemp.size();j++)
            {
                yokoyaFactor[i][j] = yokoyaFactorTemp[j];
            }
        }

    }
    fin.close();

}

double WakeFunction::GetRWLongWakeTerm2(double u, double tau0)
{
    if(u>0)
    {
        cerr<<"intergation is wrong, u can not larger than 0"<<endl;
        exit(0);
    }    
    
    double x = 0;
    double dx=0.01;
    int nz = floor(10./dx); 
    double core;
    double sum=0.0;
        
    for(int i=0;i<nz;i++)
    {
        x = i * dx;
        core  = pow(x,2) * exp(u * pow(x,2) ) / ( pow(x,6) + 8 );
        sum  += core * dx;
    }
    
    return sum;
}



vector<double> WakeFunction::GetRWPusdoWakeFun(double u) 
{
  // u is delta/sigmaT0, where sigmaT0 is the bunch length for pusedo wakefunction
  // Alex Chao Eq. 3.11 as the pusde-green fucntion in longitudinal.  Eq. 3.11 with a unit [V/C]/L
  // Alex Chao Eq. 3.56 as the pusde-green fucntion in transverse.    Eq. 3.56 with a unit rad unit. here y0 is not included in this function  
  // only the dimensionless term fu is calculated here. 
  // longitudnal wake -- head is positive -- loss energy
  // transverse  wake -- head is negative -- defocusing 
   vector<double> fu(3,0);
   if(abs(u)<1.e-2 || abs(u)>40)    // assume the pusdo wake only in the range of (-10*sigma=1cm, 10*sigma=1cm) 
   {
       return fu;
   }
   
   double x = pow(u,2)/4.0; 
       
   double a =  0.25;
   double fa    =  gsl_sf_bessel_Inu(a,   x);
   double fap1  =  gsl_sf_bessel_Inu(1+a, x);
   double fam1  =  2.0 * a / x * fa +  fap1; 
   
   double besselP025  = fa;
   double besselM075  = fam1;
   

   a =  0.75;
   fa    =  gsl_sf_bessel_Inu(a,   x);
   fap1  =  gsl_sf_bessel_Inu(1+a, x);
   fam1  =  2.0 * a / x * fa +  fap1; 

   double besselP075  = fa;
   double besselM025  = fam1;  


   double coeff1 =  pow( abs(u), 1.5) * exp (- x );   
   double coeff2 =  u > 0 ? 1 : -1;
   double coeff3 =  besselM025 - besselP075; 
   double coeff4 = -besselP025 + besselM075;
   
   fu[2] = coeff1 * ( coeff2 * coeff3 +  coeff4  );   

   coeff1 =  pow( abs(u), 0.5) * exp (- x );
   coeff2 =  u > 0 ? 1 : -1;
   coeff3 =  besselP025;
   coeff4 = -besselM025; 
   
   fu[0] = coeff1 * ( coeff2 * coeff3 +  coeff4  ); 
   fu[1] = fu[0]; 

   return fu;
}

vector<double> WakeFunction::GetRWSRWakeFun(double tau) 
{
    vector<double> wakeFun(3,0.E0);  
    double radius, pipeMatSigma,chi,s0,tau0;


    // check my mathematica notebook and USSR-CDR-20220-0106 Eq.(7-19) Eq. (7-20), which is applied here     
    /*  With resistive wakefunction to get the wake potential --  not a good idea -- a much longer time is requried, compared with a peusdo-green-wake approach. 
    for(int i=0;i<sectorRadiusX.size();i++)
    {
        radius = sectorRadiusX[i] >= sectorRadiusY[i] ? sectorRadiusY[i] : sectorRadiusX[i]; 
        pipeMatSigma = sectormatSigma[i] * FactorGaussSI;
        chi          = CLight / (4 * PI) / pipeMatSigma / radius;   //[1]
        s0           = pow(2*chi,1.0/3) * radius ;                  //[m]
        tau0         = s0/CLight ;                                  //[s]
        
        double u = tau / tau0;                                          
        double coff1 = 16 /  pow(radius,2);                         // 1/[m^2]
        double coeff2 = 1 / 3.0 * exp(u) * cos (sqrt(3.0)*u);
        double coeff3 = GetRWLongWakeTerm2(u,tau0);

        wakeFun[2] += coff1 * (coeff2 - sqrt(2.0)/PI * coeff3 ) * sectorLength[i] *  sectorNum[i] * yokoyaFactor[i][0] * FactorGaussSI;    // 1/[m] -> V/C ;       
    }
    */
      
    // us 1mm wake poten as pusdo wake function as a default value here, both longitudinal and transvese wake poten are given      
    double sigmat=1.e-3 / CLight;
    double u,coeffL,coeffT;
    vector<double> fu(3,0);          
    for(int i=0;i<sectorRadiusX.size();i++)
    {
        radius = sectorRadiusX[i] >= sectorRadiusY[i] ? sectorRadiusY[i] : sectorRadiusX[i]; 
        pipeMatSigma =  sectormatSigma[i] * FactorGaussSI;
        u            =  tau / sigmat;        
        coeffL       =  1.0 / (4 * radius) / pow(sigmat* CLight,1.5) * sqrt( CLight/(2*PI*pipeMatSigma              ));       // 1/m^2
        coeffT       =  ElecClassicRadius  / (2.0 * pow(radius,3))   * sqrt( CLight/(  PI*pipeMatSigma*sigmat*CLight));       // 1/m^2
        fu           =  GetRWPusdoWakeFun(u);                                                                                 // [dimenless] or [rad] unit                                  
        
        
        wakeFun[0]  +=  coeffT * fu[0] * sectorLength[i] *  sectorNum[i] * yokoyaFactor[i][1] * FactorGaussSI * sectorBetaX[i] / betaFunIntPoint[0];        // [1/m^2] * [m] -> [1/m] 
        wakeFun[1]  +=  coeffT * fu[1] * sectorLength[i] *  sectorNum[i] * yokoyaFactor[i][2] * FactorGaussSI * sectorBetaY[i] / betaFunIntPoint[1];        // [1/m^2] * [m] -> [1/m] 
        wakeFun[2]  +=  coeffL * fu[2] * sectorLength[i] *  sectorNum[i] * yokoyaFactor[i][0] * FactorGaussSI;                                              // [1/m^2] * [m] -> [1/m]-> V/C         
    }   
   
    // wakeFun[0] * dN(particle nuumer in bin) * y0 / gamma goes into the rad unit. momentum shift due to transverse wakefunction in this (bin). --Eq. 3.49 
	
    // longitudnal wake -- head is positive -- loss energy
    // transverse  wake -- head is negative -- defocusing

    return wakeFun;
}


vector<double> WakeFunction::GetRWLRWakeFun(double tau)
{
    // wake sign follows Alex. Chao's Fig. 2.6 Notation. 
    // return the energy loss and focusing strength. 

    // notification: 
    // longitudnal wake -- head is positive -- loss energy
    // transverse  wake -- head is negative -- defocusing

    if(tau>0)
    {
        cerr<<"Resistive wall wake calculation, tau is larger than zero"<<endl;
        exit(0);
    }

    vector<double> wakeFun(3,0.E0); // equation diverge due to approximation
    if (tau==0)                     // set wakefunction to zero when tau==0
    {
        return wakeFun;
    }
    
    // Alex Chao 2.53 and Eq.(2.76) are wake function impedance paris. The result as Nagao's Eq.(24) and Eq.(26) -- have to change to c*tau=z<0 .
	
    //in Nagaka's paper
	//wakeFun[0] =   1./(PI*pow(radius,3)) * sqrt(VaccumZ0/pipeMatSigma*CLight/PI)/pow(tau,0.5) * ringCirc;     // [V/C/m/m] positive, defocusing --per unit lengh
	//wakeFun[1] =   wakeFun[0];		
	//wakeFun[2] = - 1./(4*PI*radius)      * sqrt(VaccumZ0/pipeMatSigma/CLight/PI)/pow(tau,1.5) * ringCirc;     // [V/C/m ]  negative, beam energy loss
              
    //in Alex Chao 
    //wakeFun[0] = 2./PI/pow(radius,3) * sqrt(1/pipeMatSigma) / pow(tau,0.5)  * ringCirc * FactorGaussSI;
    //wakeFun[2] = 0;
    //wakeFun[2] = - 1./(2*PI*radius)  * sqrt(1/pipeMatSigma/CLight) / pow(tau,1.5)  * ringCirc * FactorGaussSI;  

    // in below we use the Alex Chao's notation, and it is proved, both notation gives the same results. 
     
    double radius, pipeMatSigma;
    
    for(int i=0;i<sectorRadiusX.size();i++)
    {
        radius = sectorRadiusX[i] >= sectorRadiusY[i] ? sectorRadiusY[i] : sectorRadiusX[i]; 
        pipeMatSigma = sectormatSigma[i] * FactorGaussSI;
                
        wakeFun[0] += 1./ pow(radius,3) * sqrt(1.E0/pipeMatSigma) *  sectorLength[i] *  sectorNum[i] * yokoyaFactor[i][1] * sectorBetaX[i];  // [m^-3] * [s^{1/2}] * [m] --> [m^-2] * [s^{1/2}]
        wakeFun[1] += 1./ pow(radius,3) * sqrt(1.E0/pipeMatSigma) *  sectorLength[i] *  sectorNum[i] * yokoyaFactor[i][2] * sectorBetaY[i];   
        wakeFun[2] += 1./     radius    * sqrt(1.E0/pipeMatSigma) *  sectorLength[i] *  sectorNum[i] * yokoyaFactor[i][0]; // [m^-1] * [s^{1/2}] * [m] -->          [s^{1/2}]    
    }
    // cout<< wakeFun[0] *  ( -2.) / PI * FactorGaussSI <<endl;
    // getchar();

    wakeFun[0] = wakeFun[0] *  ( -2.) / PI            /  pow( abs(tau),0.5) * FactorGaussSI / betaFunIntPoint[0] ;  // [m^-2] * [s^{1/2}] * [s^{-1/2}]  -> [m^-2] -> [s/m 1/m 1/s] ->[ohm]/m/s ->V/I/m/s->V/C/m  
    wakeFun[1] = wakeFun[1] *  ( -2.) / PI            /  pow( abs(tau),0.5) * FactorGaussSI / betaFunIntPoint[1]; 
    wakeFun[2] = wakeFun[2] /  (  2.) / PI / CLight   /  pow( abs(tau),1.5) * FactorGaussSI;                        //          [s^{1/2}] * [s^{-3/2}] [s]/[m]    -> 1/[m] -> [ohm/s] -> V/C

	return wakeFun;
    // longitudnal wake -- head is positive -- loss energy
    // transverse  wake -- head is negative -- defocusing 

}

vector<double> WakeFunction:: GetBBRPusdoWakeFun(double u, int i,double sigmat)
{
    vector<double> wakeFun(3,0.E0); 
    double lRsTemp    = lRs[i];
    double lQTemp     = lQ[i];
    double lOmegaTemp = lOmega[i];
    double tRsTemp    = tRs[i];
    double tQTemp     = tQ[i];
    double tOmegaTemp = tOmega[i];


    double sigmaz = sigmat * CLight; 
    
    // longitudinal pusdo wake poten   ---Eq. 3.14
    double vz = lOmegaTemp * sigmaz / CLight;
    int    nz   = 1000;
    double zMax = 3 * vz / (2 * lQTemp);
    double dz   = zMax / nz;

    double z=0;
    double term0,term1,term2,term3;
    double sum=0;

    for(int j=0;j<nz;j++)
    {
        z = j * dz;
        
        term0 = -1 / 2.0 * pow( u + 2 * lQTemp / vz * z ,2) - z;
        term1 = cos(z * sqrt(4 * lQTemp * lQTemp - 1));
        term2 = sin(z * sqrt(4 * lQTemp * lQTemp - 1)) / sqrt(4 * lQTemp * lQTemp - 1);
        sum  += exp(term0) * (term1 - term2) * dz;
    }
    double coeffL = sqrt(2. / PI ) * lRsTemp * CLight / sigmaz;
    wakeFun[2] = coeffL * sum;


    // tranverse pusdo wake poten  -- 3.58
    vz   = tOmegaTemp * sigmaz / CLight;
    zMax = 3 * vz / (2 * tQTemp);
    dz   = zMax /  nz;

    sum = 0;
    for(int j=0;j<nz;j++)
    {
        z = j * dz;
        
        term0 = -1 / 2.0 * pow( u + 2 * tQTemp / vz * z ,2) - z;
        term2 = sin(z * sqrt(4 * tQTemp * tQTemp - 1)) / sqrt(4 * tQTemp * tQTemp - 1);
        sum  += 2 * tQTemp / vz * exp(term0) * term2 * dz;
    }
    double coeffT = sqrt(2. / PI ) * tRsTemp * CLight;
    wakeFun[0] = sum * coeffT;
    wakeFun[1] = sum * coeffT;

    
    return wakeFun;
}

vector<double> WakeFunction::GetBBRWakeFun1(double tau) 
{
    vector<double> wakeFun(3,0);
    vector<double> wakeFunTemp;

    double sigmat=1.e-3 / CLight;  //1mm pusedo wake poten
    double u,coeffL,coeffT;

    for(int i=0;i<tRs.size();i++)
    {
        u           =  tau / sigmat;  
        wakeFunTemp = GetBBRPusdoWakeFun(u, i, sigmat);
        wakeFun[0] += wakeFunTemp[0];
        wakeFun[1] += wakeFunTemp[1];
        wakeFun[2] += wakeFunTemp[2];
    }


    return wakeFun;
}


vector<double> WakeFunction::GetBBRWakeFun(double tau)     // requires tau < 0;
{
    // longitudinal wake funciton, Refer to Alex 2.82 and 2.84 and 2.87 and 2.88    
    // wake sign follows Alex. Chao's Fig. 2.6 Notation. 
    // return the energy loss and focusing strength. 
    
    // notification: 
    // longitudnal wake -- head is positive -- loss energy
    // transverse  wake -- head is negative -- defocusing
     
    
    if(tau>0)
    {   
        cerr<<"BBR model, tau is larger than zero"<<endl;
        exit(0);
    }
    
    vector<double> wakeFun(3,0);
    double coeff,coeff1,coeff2;
    double temp;
    double omegab;
    double alpha;

    // transverse BBR wake function  Alex Chao's notation Eq. 2.88 head is negative -- defocusing
        
    for(int i=0;i<tRs.size();i++)
    {                               
        alpha  =  tOmega[i] / 2.0 / tQ[i];
        omegab =  sqrt( pow(tOmega[i],2) - pow(alpha,2) );
        coeff  =  CLight * tRs[i] * tOmega[i] / tQ[i] / omegab  * exp( alpha * tau ) * sin ( omegab * tau  ) ; //[m/s] * [ohm]/[m^2] ->  [ohm]/[m s] -> [V/C/m]
    
        wakeFun[0] +=  coeff * betaFunAver[0] / betaFunIntPoint[0];                                                    // [V/C/m]            
        wakeFun[1] +=  coeff * betaFunAver[1] / betaFunIntPoint[1];                                                    // [V/C/m]           
    }

    
    // longitudinal BBR wake function Alex Chao's Eq. 2.84--- head is postive -- loss energy
    for(int i=0;i<lRs.size();i++)
    {  
        alpha  =  lOmega[i] / 2.0 / lQ[i];
        omegab =  sqrt( pow(lOmega[i],2) - pow(alpha,2) ); 
               
        coeff  =  2 * alpha * lRs[i] * exp( alpha * tau);
        coeff1 =  cos( omegab * tau );
        coeff2 =  sin( omegab * tau ) * alpha / omegab;
                        
        if(tau==0)
        {
            wakeFun[2] +=  coeff / 2;                                                       // [V/C]   
        }
        else
        {
            wakeFun[2] +=  coeff * (coeff1 + coeff2);                                       // [V/C]  
        }
    }

    return wakeFun;
}


void WakeFunction:: GetyokoyaFactor(double a, double b, vector<double> &yokoyaFactor)
{

	double f=sqrt(abs(pow(a,2)-pow(b,2)));
	double mu= acosh(a/f);

	vector<double> coef(5,0);
	
	coef[0] = 2*sqrt(2)/PI*b/f;
	coef[1] = 	sqrt(2)/PI*pow(b/f,3);
	coef[2] = 	coef[1]; 
	coef[3] = - coef[1];
	coef[4] =   coef[1];



	vector<double> lFunPL(5,0);

	int kronekDeltaP=0;
	int kronekDeltaL=0;


	vector<double> temp(5,0);

	

	for (int L=0;L<50;L++)
	{
		if(L==0) 
		{	
			kronekDeltaL=2;
		}
		else
		{
			kronekDeltaL=1;
		}			
	
		for (int P=0;P<50;P++)
		{
						
			if(P==0)
			{	
				kronekDeltaP=2;
			}
			else
			{
				kronekDeltaP=1;
			}	

			lFunPL[0] = Lfunction(P,L,mu);
			lFunPL[3] = lFunPL[0];
			lFunPL[4] = lFunPL[0];
			lFunPL[1] = Ldxfunction(P,L,mu);
			lFunPL[2] = Ldyfunction(P,L,mu);


			temp[0] = temp[0] + lFunPL[0] *	pow(-1,P+L) * 1. / cosh(2* P   *mu) / cosh(2 *L   *mu) / kronekDeltaP / kronekDeltaL ;
			temp[1] = temp[1] +	lFunPL[1] * pow(-1,P+L) * 1. / cosh((2*P+1)*mu) / cosh((2*L+1)*mu) * (2*P+1) * (2*L+1);
			temp[2] = temp[2] +	lFunPL[2] * pow(-1,P+L) * 1. / sinh((2*P+1)*mu) / sinh((2*L+1)*mu) * (2*P+1) * (2*L+1);
			temp[3] = temp[3] +	lFunPL[3] * pow(-1,P+L) * 1. / cosh(2* P   *mu) / cosh(2* L   *mu) * pow(2*P,2) / kronekDeltaL ;
			temp[4] = temp[4] +	lFunPL[4] * pow(-1,P+L) * 1. / cosh(2* P   *mu) / cosh(2* L   *mu) * pow(2*P,2) / kronekDeltaL; 
	
		}	
	}	


	for(int i=0;i<5;i++)
	{
		yokoyaFactor[i]=coef[i]*temp[i];
	}
	

	

}

double WakeFunction:: Lfunction1 (double r,double t,double mu)
{
	double lfun1;

	int rmt = r-t;
	int rpt = r+t;
	int absrmt = abs(rmt);	

	if(r!=t)
	{
		double gamma1 = gsl_sf_gamma(0.5 + absrmt);
		double gamma2 = gsl_sf_gamma(0.5);


		gsl_sf_result  factorial;
		int stauts1 = gsl_sf_fact_e(absrmt,&factorial);	
		

		gsl_sf_result  hyper2f1;
		int  stauts2= gsl_sf_hyperg_2F1_e(0.5,absrmt+0.5,absrmt+1,exp(-4*mu), &hyper2f1);	


		lfun1 = sqrt(2.) * PI * exp(- (2*absrmt+1) * mu) * gamma1 / gamma2 / factorial.val * hyper2f1.val;
		
	}
	else
	{
		lfun1 =2 * sqrt(2) * exp(-mu) * gsl_sf_ellint_Kcomp(sqrt(exp(-4*mu)),GSL_PREC_DOUBLE);

	}
	
	return lfun1;
}

double WakeFunction:: Lfunction2 (double r,double t,double mu)
{
	
	double lfun2;

	int rmt = r-t;
	int rpt = r+t;


	double gamma1 = gsl_sf_gamma(0.5 + rpt);
	double gamma2 = gsl_sf_gamma(0.5);


	gsl_sf_result  factorial;
	int stauts1 = gsl_sf_fact_e(rpt,&factorial);	
	

	gsl_sf_result  hyper2f1;
	int  stauts2= gsl_sf_hyperg_2F1_e(0.5,rpt+0.5,rpt+1,exp(-4*mu), &hyper2f1);	

	lfun2 = sqrt(2.) * PI * exp(- (2*rpt+1) * mu) * gamma1 / gamma2 / factorial.val * hyper2f1.val;


	return lfun2;

}

double WakeFunction:: Lfunction (double r,double t, double mu)
{
	double lf1,lf2;
	
	lf1 = Lfunction1(r,t,mu);	
	lf2 = Lfunction2(r,t,mu);



	return lf1+lf2;
}

double WakeFunction:: Ldfunction1 (double r,double t,double mu)
{
	
	double ldfun1;

	int rmt = r-t;
	int rpt = r+t;
	int absrmt = abs(rmt);	


	double gamma1 = gsl_sf_gamma(0.5 + absrmt);
	double gamma2 = gsl_sf_gamma(0.5);


	gsl_sf_result  factorial;
	int stauts1 = gsl_sf_fact_e(absrmt,&factorial);	
	


	gsl_sf_result  hyper2f1;
	int  stauts2= gsl_sf_hyperg_2F1_e(0.5,absrmt+0.5,absrmt+1,exp(-4*mu), &hyper2f1);	

	ldfun1 = sqrt(2.) * PI * exp(- (2*absrmt+1) * mu) * gamma1 / gamma2 / factorial.val * hyper2f1.val;

	return ldfun1;	
}

double WakeFunction:: Ldfunction2 (double r,double t,double mu)
{
	
	double ldfun2;

	int rmt = r-t;
	int rpt = r+t;
	int absrmt = abs(rmt);	


	double gamma1 = gsl_sf_gamma(1.5 + rpt);
	double gamma2 = gsl_sf_gamma(0.5);

	gsl_sf_result  factorial;
	int stauts1 = gsl_sf_fact_e(rpt+1,&factorial);	
	
	gsl_sf_result  hyper2f1;
	int  stauts2= gsl_sf_hyperg_2F1_e(0.5,rpt+1.5,rpt+2,exp(-4*mu), &hyper2f1);	

	ldfun2 = sqrt(2.) * PI * exp(- (2*rpt+3) * mu) * gamma1 / gamma2 / factorial.val * hyper2f1.val;

	return ldfun2;	
}


double WakeFunction:: Ldxfunction (double r,double t,double mu)
{
	double ld1,ld2;
	
	ld1 = Ldfunction1(r,t,mu);
	ld2 = Ldfunction2(r,t,mu);

	return ld1+ld2;
	
}

double WakeFunction:: Ldyfunction (double r,double t,double mu)
{
	double ld1,ld2;	
	ld1 = Ldfunction1(r,t,mu);
	ld2 = Ldfunction2(r,t,mu);

	return ld1-ld2;
}





