
// Ref. Nakamura paper ""Single loop multi-dimensional digitla feedback by fir filter
// y[0] = K \sum_0^{N} a_k x[n-k] 

#include "FIRFeedBack.h"
#include <vector>
#include <iostream>
#include<numeric>
#include <fstream>
 


using namespace std;
using std::vector;
using std::complex;


FIRFeedBack::FIRFeedBack()
{

}

FIRFeedBack::~FIRFeedBack()
{

}

void FIRFeedBack::Initial(int totBunchNum, double beamEnergy)
{

    ifstream fin("FIR_input.dat",ios::in);
    if (! fin.is_open())
    {
        cerr<< "Error opening file of the FIR_input.dat"; exit (1);
    }
    
    cout<<"------------ The parameter of the FIR filter- --------"<<endl;

    vector<string> strVec;
    string         str;
    
    while (!fin.eof())
    {
        getline(fin,str);
        
        if(str.length()==0 || str.find("=")==-1 )  continue;
        StringVecSplit(str, strVec);
        
        if(strVec[0]=="delay")
        {
          delay = stod(strVec[1]);
          cout<<"the delay turns of FIR filter                    :"<< stod(strVec[1])<<"turns"<<endl;
        }
        if(strVec[0]=="taps")
        {
          taps = stod(strVec[1]);
          cout<<"the taps of the FIR filter is                    :"<< stod(strVec[1])<<"turns"<<endl;
        }
        
        if(strVec[0]=="kickerdisp")
        {
          kickerDisp = stod(strVec[1]);
          cout<<"the dispersion at kicker                           :"<< stod(strVec[1])<<"[m]"<<endl;
        }
        
        if(strVec[0]=="kickerdispp")
        {
          kickerDispP = stod(strVec[1]);
          cout<<"the Ddispersion/Ds at kicker                         :"<< stod(strVec[1])<<""<<endl;
        }
        
        if(strVec[0]=="firbunchbybunchfeedbackpowerlimit")
        {
          fIRBunchByBunchFeedbackPowerLimit = stod(strVec[1]);
          cout<<"the feed back power limit is                         :"<< stod(strVec[1])<<" W"<<endl;
        }

        if(strVec[0]=="firbunchbybunchfeedbackkickerimped")
        {
          fIRBunchByBunchFeedbackKickerImped = stod(strVec[1]);
          cout<<"the kicker impedance of feedback system               :"<< stod(strVec[1])<<" ohm"<<endl;
        }
        
        
    }
    strVec.clear();
    fin.clear();
    fin.close();

    
    fIRBunchByBunchFeedbackKickLimit = sqrt(2*fIRBunchByBunchFeedbackPowerLimit*fIRBunchByBunchFeedbackKickerImped)/beamEnergy;
    //    cout<<fIRBunchByBunchFeedbackKickLimit<<endl;
    //    getchar();


    int firOrder = delay + taps;
    
    gain.resize(3); 
    gain[0] =  1.0E+5;   //gain in x direction
    gain[1] =  1.1E+5;   //gain in y direction
    gain[2] =  1.0E+5;   //gain in z direction
    

//    gainX.resize(totBunchNum);
//    gainY.resize(totBunchNum);
//    gainZ.resize(totBunchNum);
//    
    firCoeffx.resize(firOrder);
    firCoeffy.resize(firOrder);
    firCoeffz.resize(firOrder);
    firCoeffxy.resize(firOrder);

    posxData.resize(firOrder);
    posyData.resize(firOrder);
    poszData.resize(firOrder);
    
    
    for(int i=0;i<firCoeffx.size();i++)
    {
        posxData[i].resize(totBunchNum);
        posyData[i].resize(totBunchNum);
        poszData[i].resize(totBunchNum);
    }

    for(int i=0;i<firOrder;i++)
    {
        for(int j=0;j<totBunchNum;j++)
        {
            posxData[i][j] = 0.E0;
            posyData[i][j] = 0.E0;
            poszData[i][j] = 0.E0;
        }
    }



    for(int i=0;i<delay;i++)
    {
        firCoeffx[i] = 0.E0;
        firCoeffy[i] = 0.E0;
        firCoeffz[i] = 0.E0;
        firCoeffxy[i] = 0.E0;
    }
    

    fin.open("FIR_input.dat");



    while (!fin.eof())
    {
        getline(fin,str);
        StringSplit(str, strVec);
        if(str.length()==0 )  continue;
        transform(strVec[0].begin(), strVec[0].end(), strVec[0].begin(), ::tolower);

        if(strVec[0]=="fircoeffx")
        {
            cout<<"the coefficients of of FIR for x direction  :";
            for(int i=1;i<strVec.size();i++)
            {
                firCoeffx[i-1+delay] = stod(strVec[i]);
                cout<<firCoeffx[i-1+delay]<<"       ";
            }
            cout<<endl;
        }
        
        if(strVec[0]=="fircoeffy")
        {
            cout<<"the coefficients of of FIR for y direction  :";
            for(int i=1;i<strVec.size();i++)
            {
                firCoeffy[i-1+delay] = stod(strVec[i]);
                cout<<firCoeffy[i-1+delay]<<"       ";
            }
            cout<<endl;
        }
        if(strVec[0]=="fircoeffz")
        {
            cout<<"the coefficients of of FIR for z direction  :";
            for(int i=1;i<strVec.size();i++)
            {
                firCoeffz[i-1+delay] = stod(strVec[i]);
                cout<<firCoeffz[i-1+delay]<<"       ";
            }
            cout<<endl;
        }
        if(strVec[0]=="fircoeffxy")
        {
            cout<<"the coefficients of of FIR for x-y coupling  :";
            for(int i=1;i<strVec.size();i++)
            {
                firCoeffxy[i-1+delay] = stod(strVec[i]);
                cout<<firCoeffxy[i-1+delay]<<"       ";
            }
            cout<<endl;
        }
        
    }
    
    strVec.clear();
    fin.clear();
    fin.close();



    kickStrengthKx = 1.0E-6; 
    kickStrengthKy = 1.0E-6;
    kickStrengthF = 0.E0;   


    

//    double tempCoefx=1;
//    double tempCoefy=1;

//    kickStrengthKx = tempCoefx * 1.0E-1; 
//    kickStrengthKy = tempCoefy * 1.0E-1;
//    kickStrengthF = 0.E0;

//    double sumx =   accumulate(begin(firCoeffx), end(firCoeffx), 0.0);
//    double sumy =   accumulate(begin(firCoeffy), end(firCoeffy), 0.0);
//    cout<<sumx<<"the total is "<<sumy<<endl;
//    getchar();


}
