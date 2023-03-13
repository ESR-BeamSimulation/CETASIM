//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             

#include "ReadInputSettings.h"
#include "LatticeInterActionPoint.h"
#include <vector>
#include <complex>
#include <iostream>
#include<fstream>	
#include <numeric>
#include <cmath>
#include "Global.h"
#include <iomanip>
#include <cuda_runtime.h>


using namespace std;
using std::vector;
using std::complex;



void GPU_PartiOneTurnTransfer(int macroEleNumPerBunch, double *partCord, int paraNum, double *ringPara);
void GPU_PartiSynRad(int totalPartiNum,double *partCord,int paraNum, double *radMatrix,double *coeff);
void GPU_PartiOneTurnTransferAndSynRad(int totalPartiNum, double *partCord, int oneTurnMatrixParaNum, double *ringPara, double *radMatrixBRH, double *synRadCoff);


