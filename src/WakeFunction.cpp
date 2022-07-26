

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
void WakeFunction::Initial(ReadInputSettings &inputParameter)
{
    int nTurnswakeTrunction = inputParameter.ringWake->nTurnswakeTrunction;   
    int totBunchNum         = inputParameter.ringFillPatt->totBunchNumber;    

 
    
    for(int i=0;i<3;i++)
    {
        betaFunAver[i]  = inputParameter.ringParBasic->betaFunAver[i];
    }
    
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
    
          
    ifstream fin(inputParameter.ringWake->pipeGeoInput);
    
    if (! fin.is_open())
    {
        cerr<< "Error opening file: "<< inputParameter.ringWake->pipeGeoInput <<endl; exit (1);
    }
    else
    {   // read in the geo-parameters     
        vector<string> strVec;
        string         str;
        getline(fin,str);
        
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
            
        double factor=CLight/(4*PI) * VaccumZ0;  
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

     
      
   // for bbr model input data reading in 
    
    ifstream fin1(inputParameter.ringWake->bbrInput);
 
    if (! fin1.is_open())
    {
        cerr<< "Error opening file: "<< inputParameter.ringWake->bbrInput <<endl; 
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
		    lOmega  .push_back(stod(strVec[2]));
		    tRs     .push_back(stod(strVec[3]));
		    tQ      .push_back(stod(strVec[4]));
		    tOmega  .push_back(stod(strVec[5]));

		    i++;
        }			
    }
	fin1.close();
              
}


vector<double> WakeFunction::GetRWWakeForce(double tau)
{
	// 1)Ref. NIMA 221-230 806 (2016) Nagaoka	- Eq.(24) and Eq.(26)
	// To compare wake without approximaiton, Eq (25) is used, whereas an addition factor s0/a have to added in Eq(25)

    // Alex Chao 2.53 and Eq.(2.76) are wake function impedance paris. The result as Nagao's Eq.(24) and Eq.(26) -- have to change to z>0 direction.
	
    // have to read into the ring geo-data here. To get the total wake force in total. here just for unit wake force. 

    /* in Nagaka's paper, with his notation
	wakeForce[0] =   1./(PI*pow(radius,3)) * sqrt(VaccumZ0/pipeMatSigma*CLight/PI)/pow(tau,0.5) * 1100.537;     // [V/C/m/m] positive, defocusing --per unit lengh
	wakeForce[1] =   wakeForce[0];		
	wakeForce[2] = - 1./(4*PI*radius)      * sqrt(VaccumZ0/pipeMatSigma/CLight/PI)/pow(tau,1.5) * 1100.537;     // [V/C/m ]  negative, decrease the energy
    */
            
    // from Alex Chao 2.53 wake function.   Eq.(2.76) is the related impedance used. 
    //wakeForce[0] = 2./PI/pow(radius,3) * sqrt(1/pipeMatSigma) / pow(tau,0.5)  * 1100.537 * factor;
    //wakeForce[2] =0;
    //wakeForce[2] = - 1./(2*PI*radius)  * sqrt(1/pipeMatSigma/CLight) / pow(tau,1.5)  * 1100.537 * factor;  
    // in below we use the Alex Chao's notation.

    vector<double> wakeForce(3,0.E0);
    
    double factor=CLight/(4*PI) * VaccumZ0;  
    double radius, pipeMatSigma;

    for(int i=0;i<sectorRadiusX.size();i++)
    {
        radius = sectorRadiusX[i] >= sectorRadiusY[i] ? sectorRadiusY[i] : sectorRadiusX[i]; 
        pipeMatSigma = sectormatSigma[i] * factor;

        wakeForce[0] += 1./ pow(radius,3) * sqrt(1/pipeMatSigma) *  sectorLength[i] *   sectorNum[i]  * sectorBetaX[i] * yokoyaFactor[i][1]; 
        wakeForce[1] += 1./ pow(radius,3) * sqrt(1/pipeMatSigma) *  sectorLength[i] *   sectorNum[i]  * sectorBetaY[i] * yokoyaFactor[i][2];
        //wakeForce[2] += 1/    radius    * sqrt(1/pipeMatSigma) *  sectorLength[i] *   sectorNum[i]  * yokoyaFactor[i][0];
        wakeForce[2] += 0;
    }

    wakeForce[0] = wakeForce[0] * 2. / PI           /  pow(tau,0.5) * factor / betaFunAver[0]; 
    wakeForce[1] = wakeForce[1] * 2. / PI           /  pow(tau,0.5) * factor / betaFunAver[1];
    wakeForce[2] = wakeForce[2] / 2. / sqrt(CLight) /  pow(tau,1.5) * factor; 

	return wakeForce;
}


vector<double> WakeFunction::GetBBRWakeForce(double tau)
{
    // longitudinal wake funciton, Refer to Alex 2.82 and 2.84 and 2.87 and 2.88 -- however have to change the physical reange to z>0 and unit.    
    // Here the wake are Eq. 12 and 20 to be used  Ref. NIMA 221-230 806 (2016) Nagaoka.

    
    vector<double> wakeForce(3,0);
    double coeff,coeff1,coeff2;
    double temp;
    double omegab,alpha;

    // transverse BBR wake function  Alex Chao's notation Eq. 2.88
        
    for(int i=0;i<tRs.size();i++)
    {                               
        alpha  =  tOmega[i] / 2 / tQ[i];
        omegab =  sqrt( pow(tOmega[i],2) - pow(alpha,2) );
                
        coeff   = CLight * tRs[i] * tOmega[i] / tQ[i] / omegab  * exp( - alpha * tau ) * sin (  omegab * tau  ) ; //[m/s] * [ohm]/[m^2] ->  [ohm]/[m s] -> [V/C/m]
              
        wakeForce[0] +=  coeff ;                                                    // [V/C/m]
        wakeForce[1] +=  coeff ;                                                    // [V/C/m]
        
    }



    // longitudinal BBR wake function Alex Chao's Eq. 2.84
    for(int i=0;i<lRs.size();i++)
    {  
        alpha  =  lOmega[i] / 2 / tQ[i];
        omegab =  sqrt( pow(lOmega[i],2) - pow(alpha,2) ); 
               
        coeff  =  2 * alpha * lRs[i] * exp(- alpha * tau);
        coeff1 =  cos( omegab * tau );
        coeff2 =  sin( omegab * tau ) * alpha / omegab;
        
        //coeff   = lOmega[i] * lRs[i] / lQ[i]  * exp(- lOmega[i] / 2 / lQ[i] * tau);       // [Ohm/s] ->[V/C]
        //coeff1  = cos(lOmega[i] * tau * temp); 
        //coeff2  = 1 / 2 / lQ[i] / temp * sin(- lOmega[i] * tau * temp); 
                
        if(tau==0)
        {
            wakeForce[2] +=  coeff / 2;                                                       // [V/C]
        }
        else
        {
            wakeForce[2] +=  coeff * (coeff1 + coeff2);                                       // [V/C]
        }        
    }
    wakeForce[2]=0;
    
    return wakeForce;
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





