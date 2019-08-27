#ifndef FIRFeedBack_H
#define FIRFeedBack_H

#include <vector>
#include <complex>
#include "Global.h"



using namespace std;
using std::vector;
using std::complex;


class FIRFeedBack
{


public:
    FIRFeedBack();
    ~FIRFeedBack();

	int delay;
	int taps;
	double kickStrengthKx;		// Eq. 116
	double kickStrengthKy;		// Eq. 116
	double kickStrengthF;		// Eq. 117
	double kickerDisp;			// [m] disperson function at kicker
	double kickerDispP;			//     D(disperson)/Ds function at kicker
	vector<double >  gain;
    vector<double >  firCoeffx;
    vector<double >  firCoeffy;
    vector<double >  firCoeffz;
    vector<double >  firCoeffxy;   //Nakamura's TDLSF approaches, one group pickup and kicker for 2D  
    vector<vector<double> > posxData;
    vector<vector<double> > posyData;
    vector<vector<double> > poszData;
    
    void Initial(int totBunchNum);
      
private:

};




#endif
