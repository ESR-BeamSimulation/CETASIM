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
#include "CavityResonator.h"
#include <iostream>
#include<fstream>
#include<iomanip>
	
#include <numeric>
#include <cmath>

using namespace std;
using std::vector;
using std::complex;


CavityResonator::CavityResonator()
{

}

CavityResonator::~CavityResonator()
{
    
}

void CavityResonator::Initial(ReadInputSettings &inputParameter)
{
    int resNum     = inputParameter.ringParRf->resNum;
    int ringHarmH  = inputParameter.ringParBasic->harmonics;
    double f0      = inputParameter.ringParBasic->f0;
    
           
    if(resNum==0) 
    {
        cerr<<"reset the resNum number, it has to be bigger than one"<<endl;
        exit(0);
    }
    
    resonatorVec.resize(resNum);
    
     //initialize the resonator paramters     
    for(int i=0;i<resNum;i++)
    { 
        // read from input file
        resonatorVec[i].resHarm            = inputParameter.ringParRf-> resHarm[i];
        resonatorVec[i].resType            = inputParameter.ringParRf-> resType[i];
        resonatorVec[i].resShuntImpRs      = inputParameter.ringParRf-> resShuntImpRs[i];
        resonatorVec[i].resQualityQ0       = inputParameter.ringParRf-> resQualityQ0[i];
        resonatorVec[i].resCouplingBeta    = inputParameter.ringParRf-> resCouplingBeta[i];
        resonatorVec[i].resDetuneFre       = inputParameter.ringParRf-> resDetuneFre[i];        
        resonatorVec[i].resVolAbsReq       = inputParameter.ringParRf-> resVol[i];
        resonatorVec[i].resPhaseReq        = inputParameter.ringParRf-> resPhase[i];                 
        resonatorVec[i].resCold            = inputParameter.ringParRf-> resCold[i];
             
                                                
        // calcuated....
        resonatorVec[i].resQualityQL       = resonatorVec[i].resQualityQ0 / (1 + resonatorVec[i].resCouplingBeta );                
        resonatorVec[i].resFre             = resonatorVec[i].resHarm * ringHarmH * f0 + resonatorVec[i].resDetuneFre;                
        double tanPsi                      = resonatorVec[i].resQualityQL * 
                                            (
                                              resonatorVec[i].resFre / (resonatorVec[i].resHarm * ringHarmH * f0 ) 
                                            - resonatorVec[i].resHarm * ringHarmH * f0 / resonatorVec[i].resFre     
                                            );  // Eq. (7.23)           
                                                    
        resonatorVec[i].resDeTunePsi      = atan(tanPsi); 
        resonatorVec[i].tF                = 2 * resonatorVec[i].resQualityQL / (2 * PI * resonatorVec[i].resFre );  //Eq.(7.25) [s] same notation as P.B. Wilson [s]                         
        resonatorVec[i].resGenVolFB       = complex<double>(0.e0,0.e0);
    }

        
    double  u0 = inputParameter.ringParBasic->u0;

    if(resNum==1)   // single cavity case, the default syn phase setting ensure to compensate the one turn loss
    {                
        resonatorVec[0].resPhaseReq = PI - asin(u0 /  resonatorVec[0].resVolAbsReq) - PI / 2. ;            
        resonatorVec[0].resCavVolReq   = complex<double>(resonatorVec[0].resVolAbsReq * cos(resonatorVec[0].resPhaseReq ), resonatorVec[0].resVolAbsReq * sin(resonatorVec[0].resPhaseReq ) );         
    }

    double tempU0=0.E0;    
    for (int i=0; i<resNum;i++)
    {
	    resonatorVec[i].resCavVolReq = complex<double>(resonatorVec[i].resVolAbsReq * cos(resonatorVec[i].resPhaseReq ), resonatorVec[i].resVolAbsReq * sin(resonatorVec[i].resPhaseReq ) );
        tempU0 += resonatorVec[i].resCavVolReq.real();
    }

    
    if ((tempU0-u0)>1.e5)
        cerr<<"initial settins: \sum resCavVolReq.real() - U0 > 0.1MeV " <<endl;
    
             
    /*
    // assume that the first cavity is active -- change it into cos covention. 

     
	if (resNum==2 && resonatorVec[1].resType==1)  // the second cavity is active -- RF settings supply an ideal bunch lengthing effect  
	{
		int n = resonatorVec[1].resHarm;
		double vRF = resonatorVec[0].resVolAbsReq;
		
		double phaseTemp = pow(n,2)/(pow(n,2)-1) * u0 / vRF;			
		resonatorVec[0].resPhaseReq = PI -  asin(phaseTemp) - PI/2. ;	
        
            						
		phaseTemp = - n * u0 / vRF / sqrt( pow(pow(n,2)-1,2) - pow(n,4) * pow(u0 / vRF,2)  );
		resonatorVec[1].resPhaseReq = atan(phaseTemp) - PI/2;
	
		// set the the main cavity as single RF value--HC as  (-PI/2) 
		//resonatorVec[0].resPhaseReq = PI - asin(u0 /  resonatorVec[0].resVolAbsReq) - PI/2.;        
		//resonatorVec[1].resPhaseReq = - PI/2.0;

	    double kHC = sqrt( 1./pow(n,2) - pow( u0 / vRF ,2) / ( pow(n,2)-1)  ); 
	    resonatorVec[1].resVolAbsReq = kHC * resonatorVec[0].resVolAbsReq;        
              
	    resonatorVec[0].resCavVolReq = complex<double>(resonatorVec[0].resVolAbsReq * cos(resonatorVec[0].resPhaseReq ), resonatorVec[0].resVolAbsReq * sin(resonatorVec[0].resPhaseReq ) );
	    resonatorVec[1].resCavVolReq = complex<double>(resonatorVec[1].resVolAbsReq * cos(resonatorVec[1].resPhaseReq ), resonatorVec[1].resVolAbsReq * sin(resonatorVec[1].resPhaseReq ) ); 
	    
	   	    
	}
	else if(resNum==2 && resonatorVec[1].resType==0) // the second cavity is passive 
	{
	    // to do Ref. Yamamoto's Eq. (6) and (7). 
	}
    */
 
}



double CavityResonator::GetWakefunction(double tau)  
{
    //BBR model as natural  inlcluding the self-loss when tauij=0
    //longitudinal BBR wake function Alex Chao's Eq. 2.84
    //with the same notation, tau has to be smaller that zero;
    
    
    double wakefun = 0;
    double coeff,coeff1,coeff2;
    double temp;
    double omegab,alpha,omegaRes;    
    
    if(tau>0)
    {
        cerr<<" subroutine in wake function is not correct, tau must be <=0 "<<endl;
        exit(0);
    }
    
    for(int i=0;i<resonatorVec.size();i++)
    {  
        if(resonatorVec[i].resQualityQL<0.5)
        {
            cerr<<"cavity loaded quality factor has to be larger than 0.5, Alex Chao' Eq. 2.84"<<endl;
            exit(0);
        }        
        
        omegaRes = 2 * PI * resonatorVec[i].resFre;
        alpha    = omegaRes  / 2.0 / resonatorVec[i].resQualityQL;        
        omegab   = sqrt( pow(omegaRes,2) - pow(alpha,2) ); 
                                       
        coeff  =  2 * alpha * resonatorVec[i].resShuntImpRs * exp( alpha * tau);
        coeff1 =  cos( omegab * tau );
        coeff2 =  sin( omegab * tau ) * alpha / omegab;
        
                
        if(tau==0)
        {
            wakefun +=  coeff / 2;                                                       // [V/C]
        }
        else
        {
            wakefun +=  coeff * (coeff1 + coeff2);                                       // [V/C]
        }        
    }
       
    return wakefun;
               
}



