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
#include "Global.h"
#include <iomanip>
#include <cuda_runtime.h>
#include <cufftXt.h>
#include <cufft.h>
#include "CUDAFunction.cuh"

using namespace std;
using std::vector;
using std::complex;

__global__ void increment(double *partCord, double *oneTurnMap,  const int dim)
{
    // GPU only deal with matrix multiplication, and matrix is fixed.
}

void GPU_PartiOneTurnTransfer(double *partCord, double *oneTurnMap, const int dim)
{
    
    double *d_partCord,*d_oneTurnMap;
    cudaMalloc((void**) &d_partCord, dim * sizeof(double) );
    cudaMemcpy(d_partCord,partCord,dim * sizeof(double), cudaMemcpyHostToDevice);
    
    cudaMalloc((void**) &d_oneTurnMap, 36 * sizeof(double) );
    cudaMemcpy(d_oneTurnMap,oneTurnMap,36 * sizeof(double), cudaMemcpyHostToDevice);

    dim3 block(1024);
    dim3 grid((dim + block.x -1)/ block.x);

    increment<<<grid, block>>>(d_partCord,d_oneTurnMap,dim);
    
    cudaFree(d_partCord);
    cudaFree(d_oneTurnMap);

    cout<<"test"<<endl;
    getchar();
}



