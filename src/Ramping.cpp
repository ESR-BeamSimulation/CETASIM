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
#include<iomanip>
#include <numeric>
#include <cmath>
#include "Ramping.h"



using namespace std;
using std::vector;
using std::complex;



Ramping::Ramping()
{

}


Ramping::~Ramping()
{
    
}

void Ramping::RampingPara(ReadInputSettings &inputParameter, LatticeInterActionPoint &latticeInterActionPoint,int n)
{
    int nTurns       = inputParameter.ringRun->nTurns;  
    int dTurns       = inputParameter.ramping->deltaTurns;
    
    int flag = (n%dTurns==0) && (n >= inputParameter.ramping->rampingTurns[0]) && ( n <= inputParameter.ramping->rampingTurns[1]);
    if(flag==1)
    {
        // linear ramping the tunes
        if(inputParameter.ramping->rampingNu[0]!=0 ) inputParameter.ringParBasic->workQx +=  inputParameter.ramping->deltaNuPerTurn[0] * dTurns;
        if(inputParameter.ramping->rampingNu[1]!=0 ) inputParameter.ringParBasic->workQy +=  inputParameter.ramping->deltaNuPerTurn[1] * dTurns;
        if(inputParameter.ramping->rampingSKQ!=0)    inputParameter.ringParBasic->skewQuadK += inputParameter.ramping->deltaSKQKPerTurn* dTurns;

        // linear ramping the skew quad       
        // if(n<=nTurns/2)
        // {
        //     if(inputParameter.ramping->rampingSKQ!=0)    inputParameter.ringParBasic->skewQuadK += inputParameter.ramping->deltaSKQKPerTurn* dTurns ;
        // }
        // else
        // {
        //     if(inputParameter.ramping->rampingSKQ!=0)    inputParameter.ringParBasic->skewQuadK -= inputParameter.ramping->deltaSKQKPerTurn* dTurns ;
        // }
        latticeInterActionPoint.GetTransLinearCouplingCoef(inputParameter);
        latticeInterActionPoint.SetLatticeBRHForSynRad(inputParameter);
        // cout<<setw(15)<<left<<setprecision(9)<<flag
        //     <<setw(15)<<left<<setprecision(9)<<n
        //     <<setw(15)<<left<<setprecision(9)<<inputParameter.ringParBasic->workQy
        //     <<setw(15)<<left<<setprecision(9)<<inputParameter.ringParBasic->workQx
        //     <<endl;
    }
    
    
}


void Ramping::RampSKQ(LatticeInterActionPoint &latticeInterActionPoint,int n)
{
     
}














