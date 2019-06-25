
// Ref. Nakamura paper ""Single loop multi-dimensional digitla feedback by fir filter
// y[0] = K \sum_0^{N} a_k x[n-k] 

#include "FIRFeedBack.h"
#include <vector>
#include <iostream>
 #include<numeric>


using namespace std;
using std::vector;
using std::complex;


FIRFeedBack::FIRFeedBack()
{

}

FIRFeedBack::~FIRFeedBack()
{

}

void FIRFeedBack::Initial()
{
	delay = 1;
	taps  = 9;
	kickerDisp = 0;
	kickerDispP = 0;
	int firOrder = delay + taps;
	
	gain.resize(firOrder);
	firCoeffx.resize(firOrder);
	firCoeffy.resize(firOrder);
	firCoeffz.resize(firOrder);
	firCoeffxy.resize(firOrder);
	
	posxData.resize(firOrder);
	posyData.resize(firOrder);
	poszData.resize(firOrder);

	for(int i=0;i<delay;i++)
	{
		firCoeffx[i] = 0.E0;
		firCoeffy[i] = 0.E0;
	}
	
//	for(int i=0;i<firOrder;i++)
//	{
//		firCoeffy[i] = 0.E0;
//	}
	
	
	
	firCoeffx[1] = 0.3039;
	firCoeffx[2] = 0.6124;
	firCoeffx[3] = 0.3034;
	firCoeffx[4] =-0.1764;
	firCoeffx[5] =-0.3555;
	firCoeffx[6] =-0.1712;
	firCoeffx[7] = 0.0347;
	firCoeffx[8] =-0.0811;
	firCoeffx[9] =-0.4701;
	
	
	
	firCoeffy[1] = 0.6923;
	firCoeffy[2] = 0.1701;
	firCoeffy[3] =-0.5195;
	firCoeffy[4] =-0.2701;
	firCoeffy[5] = 0.0343;
	firCoeffy[6] =-0.1399;
	firCoeffy[7] =-0.0406;
	firCoeffy[8] = 0.2606;
	firCoeffy[9] =-0.1874;
	
	for(int i=0;i<firOrder;i++)
	{
		firCoeffxy[i] = 0.E0;
		firCoeffz[i] = 0.E0;
	}
	
	double tempCoefx=1;
	double tempCoefy=1;
	
	kickStrengthKx = tempCoefx * 1.0E-1;
	kickStrengthKy = tempCoefy * 1.0E-1;
	kickStrengthF = 0.E0;
	
	double sumx =   accumulate(begin(firCoeffx), end(firCoeffx), 0.0);
	double sumy =   accumulate(begin(firCoeffx), end(firCoeffx), 0.0);
//	cout<<sumx<<"the total is "<<sumy<<endl;
//	getchar();
	
	
	
	
}
