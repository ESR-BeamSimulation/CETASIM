//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             

#include "Bunch.h"
#include "Global.h"
#include "BoardBandImp.h"
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <numeric>
#include <cmath>
#include <random>
#include <gsl/gsl_fft_complex.h>
#include <algorithm>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_histogram.h>
#include <fftw3.h>

using namespace std;
using std::vector;
using std::complex;

BoardBandImp::BoardBandImp()
{
}

BoardBandImp::~BoardBandImp()
{

}

void BoardBandImp::ReadInImp(const  ReadInputSettings &inputParameter)
{    
    string fileName = inputParameter.ringBBImp->impedInput;

    int lineNumber = 0;
    string str;
    vector<string> strVec;

    // file is from sddsprintout command and impedance is with ele
    ifstream fin(fileName);
    while (!fin.eof())
    {
        if(lineNumber<5)
        {
            getline(fin,str);		 
            lineNumber++;
            continue;
        }
        getline(fin,str);
        if(str.length()==0)  continue;
        StringSplit2(str,strVec);
        // change to Alex Chao' notation
        freq.push_back(stod(strVec[0]));
        zZImp.push_back ( complex<double>( stod(strVec[1]),  -stod(strVec[6])) );
        zDxImp.push_back( complex<double>(-stod(strVec[7]),  -stod(strVec[2])) );
        zDyImp.push_back( complex<double>(-stod(strVec[8]),  -stod(strVec[3])) );
        
        zQxImp.push_back( complex<double>(-stod(strVec[9]),  -stod(strVec[4])) );
        zQyImp.push_back( complex<double>(-stod(strVec[10]),  -stod(strVec[5])) );
    }

    // vector<double> freq0   = freq ;
    // vector<complex<double>> zZImp0  = zZImp ;
    // vector<complex<double>> zDxImp0 = zDxImp ;
    // vector<complex<double>> zDyImp0 = zDyImp ;

    // // to set the minus frequency part of the impedance
    // reverse(freq.begin(),freq.end());
    // reverse(zZImp.begin(),zZImp.end());
    // reverse(zDxImp.begin(),zDxImp.end());
    // reverse(zDyImp.begin(),zDyImp.end());

    // freq.pop_back();
    // zZImp.pop_back();
    // zDxImp.pop_back();
    // zDyImp.pop_back();

    // for(int i=0;i<freq.size();i++)
    // {
    //     freq[i]    *= (-1.E0);
    //     zZImp[i]    = conj(zZImp[i]);
    //     zDxImp[i]   = conj(zDxImp[i]) * (-1.E0);
    //     zDyImp[i]   = conj(zDyImp[i]) * (-1.E0);
    // }

    // freq.insert(    freq.end(),  freq0.begin(),  freq0.end());
    // zZImp.insert(  zZImp.end(), zZImp0.begin(), zZImp0.end());
    // zDxImp.insert(zDxImp.end(),zDxImp0.begin(),zDxImp0.end());
    // zDyImp.insert(zDyImp.end(),zDyImp0.begin(),zDyImp0.end());


    double rBeta        = inputParameter.ringParBasic->rBeta;

    nBins   = freq.size();                 
    freqMax = freq.back();
    dt      = 1.0 / (2 * freqMax);
    tMax    = dt * (nBins - 1) ; 
    dz      = dt * CLight * rBeta;      
    zMax    = dz * (nBins - 1);

    binPosZ.resize(2*nBins-1,0.E0);

    for(int i=0;i<binPosZ.size();i++ )
    {
        binPosZ[i] = (- tMax + i * dt ) * CLight * rBeta;
    }


    // test fft from impedance to wake directly.....
    // double Rs = 1;
    // double fR = 1.0E+9;
    // double Q = 10;
    // int nn=pow(2,16);
    // double df = 10 * fR / nn ;

    // // double alpha  = 2.0 * PI * fR / 2 / Q;
    // // double omegab = sqrt(pow(2.0*PI*fR,2) - pow(alpha,2) ); 
    
    // double freTemp=0;
    // complex<double>  zZImpTemp[nn];


    // ofstream fout("test.dat");
    // for(int i=0;i<nn;i++ )
    // {
    //     freTemp = i * df;
    //     if(freTemp==0)
    //     {
    //         zZImpTemp[i] = complex<double>(0.0,0.0);
    //     }
    //     else
    //     {
    //         zZImpTemp[i]   = Rs / (1.0 + li * Q * ( fR / freTemp -  freTemp / fR));    // [C/s*Ohm] = [V]
    //     }
    // }
    
    // fftw_complex *c2rin   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nn );
    
    // for(int i=0;i<nn ;i++)
    // {
    //     c2rin[i][0] = zZImpTemp[i].real();
    //     c2rin[i][1] = zZImpTemp[i].imag(); 
    // }

    // double temp[2*nn+1];
    
    // fftw_plan p = fftw_plan_dft_c2r_1d(nn, c2rin, temp, FFTW_ESTIMATE);    // ffw, according to g++/nvcc complile flag is FFTW_MEASURE/FFTW_ESTIMATE
    // fftw_execute(p);

    
    
    // fout<<"SDDS1"<<endl;
    // fout<<"&column name=z,              units=m,              type=float,  &end" <<endl;
    // fout<<"&column name=wakeZ,                                 type=float,  &end" <<endl;
    // fout<<"&data mode=ascii, &end"                                               <<endl;
    // fout<<"! page number " << 1 <<endl;
    // fout<<nn<<endl;


    // for(int i=0;i<2*nn;i++)
    // {
    //     double posz = i*dz;
    //     fout<<setw(15)<<left<<posz
    //         <<setw(15)<<left<<temp[i]/nn
    //         // <<2*Rs*alpha*exp(-alpha * posz / CLight) * (cos(-omegab * posz / CLight) + alpha / omegab * sin(-omegab * posz / CLight)  )
    //         <<endl;
    // }

    // cout<<__LINE__<<__FILE__<<endl;
    // getchar();

    
}


