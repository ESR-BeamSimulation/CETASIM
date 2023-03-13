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
#include "Global.h"
#include "CUDAFunction.cuh"
#include "math_constants.h"
#include <curand.h>
#include <cuda.h>
#include <curand_kernel.h>



using namespace std;
using std::vector;
using std::complex;

__global__ void GPU_OneTurnMap(int partNum, double *partCord, double *ringPara )
{
    // partCord {x,px,y,py,z,pz,...}

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // from ringPara input
    double alpha[2]  = {ringPara[0],ringPara[1]};
    double beta[2]   = {ringPara[2],ringPara[3]};
    double nu[2]     = {ringPara[4],ringPara[5]};
    double chrom[2]  = {ringPara[6],ringPara[7]};
    double eta[2]    = {ringPara[8],ringPara[9]};
    double etap[2]   = {ringPara[10],ringPara[11]};
    double aDTX[2]   = {ringPara[12],ringPara[13]};
    double aDTY[2]   = {ringPara[14],ringPara[15]};
    double aDTXY[2]  = {ringPara[16],ringPara[17]};
    double alphac[3] = {ringPara[18],ringPara[19],ringPara[20]};
    double circRing  = ringPara[21];
  

    // // temp during calculation
    double amp[2],nutmp[2],phi[2],gamma[2],sinPhi[2],cosPhi[2];  // used in x and y direction
    double oneTurnMap[6][6];
    memset(oneTurnMap,0,sizeof(oneTurnMap));
    double vectX[6]={0,0,0,0,0,0};
    double vectY[6]={0,0,0,0,0,0};
    int partIndex = idx;


    if(partIndex < partNum)
    {
        #pragma unroll 
        for(int j=0;j<6;j++) vectX[j]  = partCord[partIndex * 6 + j];   // set the vectX (x,px,y,py,z,pz)
                     
        vectX[0] -= eta[0]  * vectX[5];   //x 
        vectX[1] -= etap[0] * vectX[5];   //px
        vectX[2] -= eta[1]  * vectX[5];   //y
        vectX[3] -= etap[1] * vectX[5];   //py

    
        #pragma unroll 
        for(int j=0;j<2;j++)
        {           
            amp[j]  = ( vectX[2*j] * vectX[2*j] + (alpha[j]*vectX[2*j] + beta[j]*vectX[2*j+1]) * (alpha[j]*vectX[2*j] + beta[j]*vectX[2*j+1]) ) / beta[j];
        }
        
        nutmp[0] = nu[0] + chrom[0] * vectX[5] +  aDTX[0] * amp[0] + aDTX[1] * amp[0] * amp[0]  / 2 + aDTXY[0] * amp[0] * amp[1];
        nutmp[1] = nu[1] + chrom[1] * vectX[5] +  aDTY[0] * amp[1] + aDTY[1] * amp[1] * amp[1]  / 2 + aDTXY[1] * amp[0] * amp[1];
        
        #pragma unroll
        for(int j=0;j<2;j++)
        {
            phi[j]      = 2.0  * M_PI * nutmp[j];
            sinPhi[j]   = sin( phi[j] );
            cosPhi[j]   = cos( phi[j] );
            gamma[j]    = (1 + alpha[j] * alpha[j] ) / beta[j];    
        }
        

        oneTurnMap[0][0] = cosPhi[0] + alpha[0] * sinPhi[0];
        oneTurnMap[0][1] =              beta[0] * sinPhi[0];
        oneTurnMap[1][0] =           - gamma[0] * sinPhi[0];
        oneTurnMap[1][1] = cosPhi[0] - alpha[0] * sinPhi[0];


        oneTurnMap[2][2] = cosPhi[1] + alpha[1] * sinPhi[1];
        oneTurnMap[2][3] =              beta[1] * sinPhi[1];
        oneTurnMap[3][2] =           - gamma[1] * sinPhi[1];
        oneTurnMap[3][3] = cosPhi[1] - alpha[1] * sinPhi[1];

        
        oneTurnMap[0][5] =  eta[0]  -  eta[0] * cosPhi[0] - (alpha[0] * eta[0] + beta[0] * etap[0]) * sinPhi[0];
        oneTurnMap[2][5] =  eta[1]  -  eta[1] * cosPhi[1] - (alpha[1] * eta[1] + beta[1] * etap[1]) * sinPhi[1];
        oneTurnMap[1][5] = -etap[0] - etap[0] * cosPhi[0] + ( eta[0] + alpha[0] * alpha[0] * eta[0] + alpha[0] * beta[0] * etap[0] ) * sinPhi[0] / beta[0];
        oneTurnMap[3][5] = -etap[1] - etap[1] * cosPhi[1] + ( eta[1] + alpha[1] * alpha[1] * eta[1] + alpha[1] * beta[1] * etap[1] ) * sinPhi[1] / beta[1];

        oneTurnMap[4][0] = -etap[0] + etap[0] * cosPhi[0] + ( eta[0] + alpha[0] * alpha[0] * eta[0] + alpha[0] * beta[0] * etap[0] ) * sinPhi[0] / beta[0];
        oneTurnMap[4][2] = -etap[1] + etap[1] * cosPhi[1] + ( eta[1] + alpha[1] * alpha[1] * eta[1] + alpha[1] * beta[1] * etap[1] ) * sinPhi[1] / beta[1];

        oneTurnMap[4][1] = eta[0]   - eta[0] * cosPhi[0] + (alpha[0] * eta[0] + beta[0] * etap[0]) * sinPhi[0];
        oneTurnMap[4][3] = eta[1]   - eta[1] * cosPhi[1] + (alpha[1] * eta[1] + beta[1] * etap[1]) * sinPhi[1];


        oneTurnMap[4][4] = 1;
        oneTurnMap[5][5] = 1;

        #pragma unroll
        for(int j=0;j<6;j++)
        {
            #pragma unroll
            for(int k=0;k<6;k++)
            {
                vectY[j] += oneTurnMap[j][k] * vectX[k];
            }
        }
        vectY[4] -= circRing * (alphac[0] * vectY[5]  + alphac[1] * vectY[5] * vectY[5] + alphac[2] * vectY[5] * vectY[5] * vectY[5] );
        
        #pragma unroll
        for(int j=0;j<6;j++) partCord[partIndex * 6 + j] = vectY[j];
    }
    
}

void GPU_PartiOneTurnTransfer(int totalPartiNum, double *partCord, int paraNum,double *ringPara)
{
    
    int partCordMem = totalPartiNum * 6 * sizeof(double);
    int paraMem     = paraNum *  sizeof(double); 
    double *d_partCord,*d_ringPara;
    
    cudaMalloc((void**) &d_partCord, partCordMem);
    cudaMalloc((void**) &d_ringPara, paraMem);
    cudaMemcpy(d_partCord,partCord,  partCordMem, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ringPara,ringPara,  paraMem,     cudaMemcpyHostToDevice);
    
    // each particle is locate at each thread on GPU. 
    int blockDim;
    if(totalPartiNum>1000) 
    {
        blockDim =128;
    }
    else
    {
        int count = 0;
        int temp=totalPartiNum;
        while(temp/2!=0)
        {
            count++;
            temp /=2;    
        }
        blockDim = pow(2,count);
    }
  
    dim3 block(blockDim);
    dim3 grid((totalPartiNum + block.x -1)/ block.x);

    GPU_OneTurnMap<<<grid, block>>>(totalPartiNum,d_partCord,d_ringPara);

    cudaMemcpy(partCord, d_partCord, partCordMem, cudaMemcpyDeviceToHost);
    cudaFree(d_partCord);
    cudaFree(d_ringPara);

}


__global__ void GPU_OneTurnSynRad(int partNum, double *partCord, double *radMatrixBRH, double *radCoeff, double *gaussRand)
{
    // partCord {x,px,y,py,z,pz,...}curandStateXORWOW_t  curandStateXORWOW_t* states

    int pId = blockIdx.x * blockDim.x + threadIdx.x;

    if(pId < partNum)
    {
        double vectX[6]={0,0,0,0,0,0};
        double vectY[6]={0,0,0,0,0,0};

        #pragma unroll 
        for(int i=0;i<6;i++) vectX[i]  = partCord[pId * 6 + i];   // set the vectX (x,px,y,py,z,pz)
        
        // transfer to norm cord from acc cord
        #pragma unroll
        for(int i=0;i<6;i++)
        {
            #pragma unroll
            for(int j=0;j<6;j++)
            {
                vectY[i] += radMatrixBRH[6*i +j] * vectX[j];
            }
        }
        
        vectY[0] *=  radCoeff[0];
        vectY[1] *=  radCoeff[0];
        vectY[2] *=  radCoeff[1];
        vectY[3] *=  radCoeff[1];
        vectY[5] *=  radCoeff[2];

        
        vectY[0] += radCoeff[3] * gaussRand[pId * 6 + 0]; 
        vectY[1] += radCoeff[3] * gaussRand[pId * 6 + 1]; 
        vectY[2] += radCoeff[4] * gaussRand[pId * 6 + 2]; 
        vectY[3] += radCoeff[4] * gaussRand[pId * 6 + 3]; 
        // vectY[4] += radCoeff[5] * gaussRand[pId * 6 + 4]; 
        vectY[5] += radCoeff[5] * gaussRand[pId * 6 + 5]; 

        // double gaussRand = curand_normal_double(state);
        // printf("%15.8f\n",gaussRand[pId * 6]);

        #pragma unroll
        for(int i=0;i<6;i++)
        {
            vectX[i] = 0;
            #pragma unroll
            for(int j=0;j<6;j++)
            {
                vectX[i] += radMatrixBRH[6*i +j + 36] * vectY[j];
            }
        }
         
        // transfer to acc cord from norm cord
        #pragma unroll
        for(int j=0;j<6;j++) partCord[pId * 6 + j] = vectX[j];
    }
    
}


void GPU_PartiSynRad(int totalPartiNum, double *partCord, int paraNum, double *radMatrixBRH, double *radCoeff)
{
    int partCordMem      = totalPartiNum * 6 * sizeof(double);
    int radMatrixBRHMem  = paraNum *  sizeof(double);
    int radCeoffMem      = 6 * sizeof(double);
    
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> dx{0,1};

    double gaussRand[totalPartiNum*6];
    for(int i=0;i<totalPartiNum*6;i++) gaussRand[i    ] = dx(gen);
    for(int i=0;i<totalPartiNum;  i++) gaussRand[6*i+4] = 0;
    
    double *d_partCord, *d_radMatrixBRH,*d_radCoeff,*d_gaussRand;
    cudaMalloc((void**) &d_partCord,     partCordMem);
    cudaMalloc((void**) &d_gaussRand,    partCordMem);
    cudaMalloc((void**) &d_radMatrixBRH, radMatrixBRHMem);
    cudaMalloc((void**) &d_radCoeff,     radCeoffMem);
    

    cudaMemcpy(d_partCord,     partCord,      partCordMem,         cudaMemcpyHostToDevice);
    cudaMemcpy(d_gaussRand,    gaussRand,     partCordMem,         cudaMemcpyHostToDevice);
    cudaMemcpy(d_radMatrixBRH, radMatrixBRH,  radMatrixBRHMem,     cudaMemcpyHostToDevice);
    cudaMemcpy(d_radCoeff,     radCoeff,      radCeoffMem,         cudaMemcpyHostToDevice);  

    int blockDim;
    if(totalPartiNum>1000) 
    {
        blockDim =128;
    }
    else
    {
        int count = 0;
        int temp=totalPartiNum; 
        while(temp/2!=0)
        {
            count++;
            temp /=2;    
        }
        blockDim = pow(2,count);
    }

    dim3 block(blockDim);
    dim3 grid((totalPartiNum + block.x -1)/ block.x);

    GPU_OneTurnSynRad<<<grid, block>>>(totalPartiNum,d_partCord,d_radMatrixBRH,d_radCoeff,d_gaussRand);

    cudaMemcpy(partCord, d_partCord, partCordMem, cudaMemcpyDeviceToHost);
    // cout<<"host"<<endl;
    // getchar();
    cudaFree(d_partCord);
    cudaFree(d_radMatrixBRH);
    cudaFree(d_radCoeff);
}


__global__ void  GPU_OneTurnMapAndSynRad(int partNum,double *partCord,double *ringPara,double *radMatrixBRH,double *radCoeff,double *gaussRand)
{

    int pId = blockIdx.x * blockDim.x + threadIdx.x;

    if(pId < partNum)
    {
        // one turn map transfer 
        double alpha[2]  = {ringPara[0],ringPara[1]};
        double beta[2]   = {ringPara[2],ringPara[3]};
        double nu[2]     = {ringPara[4],ringPara[5]};
        double chrom[2]  = {ringPara[6],ringPara[7]};
        double eta[2]    = {ringPara[8],ringPara[9]};
        double etap[2]   = {ringPara[10],ringPara[11]};
        double aDTX[2]   = {ringPara[12],ringPara[13]};
        double aDTY[2]   = {ringPara[14],ringPara[15]};
        double aDTXY[2]  = {ringPara[16],ringPara[17]};
        double alphac[3] = {ringPara[18],ringPara[19],ringPara[20]};
        double circRing  = ringPara[21];


        // // temp during calculation
        double amp[2],nutmp[2],phi[2],gamma[2],sinPhi[2],cosPhi[2];  // used in x and y direction
        double oneTurnMap[6][6];
        memset(oneTurnMap,0,sizeof(oneTurnMap));
        double vectX[6]={0,0,0,0,0,0};
        double vectY[6]={0,0,0,0,0,0};
        
        #pragma unroll 
        for(int j=0;j<6;j++) vectX[j]  = partCord[pId * 6 + j];   // set the vectX (x,px,y,py,z,pz)
                     
        vectX[0] -= eta[0]  * vectX[5];   //x 
        vectX[1] -= etap[0] * vectX[5];   //px
        vectX[2] -= eta[1]  * vectX[5];   //y
        vectX[3] -= etap[1] * vectX[5];   //py

    
        #pragma unroll 
        for(int j=0;j<2;j++)
        {           
            amp[j]  = ( vectX[2*j] * vectX[2*j] + (alpha[j]*vectX[2*j] + beta[j]*vectX[2*j+1]) * (alpha[j]*vectX[2*j] + beta[j]*vectX[2*j+1]) ) / beta[j];
        }
        
        nutmp[0] = nu[0] + chrom[0] * vectX[5] +  aDTX[0] * amp[0] + aDTX[1] * amp[0] * amp[0]  / 2 + aDTXY[0] * amp[0] * amp[1];
        nutmp[1] = nu[1] + chrom[1] * vectX[5] +  aDTY[0] * amp[1] + aDTY[1] * amp[1] * amp[1]  / 2 + aDTXY[1] * amp[0] * amp[1];
        
        #pragma unroll
        for(int j=0;j<2;j++)
        {
            phi[j]      = 2.0  * M_PI * nutmp[j];
            sinPhi[j]   = sin( phi[j] );
            cosPhi[j]   = cos( phi[j] );
            gamma[j]    = (1 + alpha[j] * alpha[j] ) / beta[j];    
        }
        

        oneTurnMap[0][0] = cosPhi[0] + alpha[0] * sinPhi[0];
        oneTurnMap[0][1] =              beta[0] * sinPhi[0];
        oneTurnMap[1][0] =           - gamma[0] * sinPhi[0];
        oneTurnMap[1][1] = cosPhi[0] - alpha[0] * sinPhi[0];


        oneTurnMap[2][2] = cosPhi[1] + alpha[1] * sinPhi[1];
        oneTurnMap[2][3] =              beta[1] * sinPhi[1];
        oneTurnMap[3][2] =           - gamma[1] * sinPhi[1];
        oneTurnMap[3][3] = cosPhi[1] - alpha[1] * sinPhi[1];

        
        oneTurnMap[0][5] =  eta[0]  -  eta[0] * cosPhi[0] - (alpha[0] * eta[0] + beta[0] * etap[0]) * sinPhi[0];
        oneTurnMap[2][5] =  eta[1]  -  eta[1] * cosPhi[1] - (alpha[1] * eta[1] + beta[1] * etap[1]) * sinPhi[1];
        oneTurnMap[1][5] = -etap[0] - etap[0] * cosPhi[0] + ( eta[0] + alpha[0] * alpha[0] * eta[0] + alpha[0] * beta[0] * etap[0] ) * sinPhi[0] / beta[0];
        oneTurnMap[3][5] = -etap[1] - etap[1] * cosPhi[1] + ( eta[1] + alpha[1] * alpha[1] * eta[1] + alpha[1] * beta[1] * etap[1] ) * sinPhi[1] / beta[1];

        oneTurnMap[4][0] = -etap[0] + etap[0] * cosPhi[0] + ( eta[0] + alpha[0] * alpha[0] * eta[0] + alpha[0] * beta[0] * etap[0] ) * sinPhi[0] / beta[0];
        oneTurnMap[4][2] = -etap[1] + etap[1] * cosPhi[1] + ( eta[1] + alpha[1] * alpha[1] * eta[1] + alpha[1] * beta[1] * etap[1] ) * sinPhi[1] / beta[1];

        oneTurnMap[4][1] = eta[0]   - eta[0] * cosPhi[0] + (alpha[0] * eta[0] + beta[0] * etap[0]) * sinPhi[0];
        oneTurnMap[4][3] = eta[1]   - eta[1] * cosPhi[1] + (alpha[1] * eta[1] + beta[1] * etap[1]) * sinPhi[1];


        oneTurnMap[4][4] = 1;
        oneTurnMap[5][5] = 1;

        #pragma unroll
        for(int j=0;j<6;j++)
        {
            #pragma unroll
            for(int k=0;k<6;k++)
            {
                vectY[j] += oneTurnMap[j][k] * vectX[k];
            }
        }
        vectY[4] -= circRing * (alphac[0] * vectY[5]  + alphac[1] * vectY[5] * vectY[5] + alphac[2] * vectY[5] * vectY[5] * vectY[5] );
        
        #pragma unroll
        for(int j=0;j<6;j++) partCord[pId * 6 + j] = vectY[j];
        //------------------------end of one turn map--------------------------------------------
        
       
       // synchrontron radiation damping--------------------------------------------- 

        #pragma unroll 
        for(int i=0;i<6;i++) vectX[i]  = partCord[pId * 6 + i];   // set the vectX (x,px,y,py,z,pz)
        
        // transfer to norm cord from acc cord
        #pragma unroll
        for(int i=0;i<6;i++)
        {
            vectY[i] = 0.E0;
            #pragma unroll
            for(int j=0;j<6;j++)
            {
                vectY[i] += radMatrixBRH[6*i +j] * vectX[j];
            }
        }
        
        vectY[0] *=  radCoeff[0];
        vectY[1] *=  radCoeff[0];
        vectY[2] *=  radCoeff[1];
        vectY[3] *=  radCoeff[1];
        vectY[5] *=  radCoeff[2];

        
        vectY[0] += radCoeff[3] * gaussRand[pId * 6 + 0]; 
        vectY[1] += radCoeff[3] * gaussRand[pId * 6 + 1]; 
        vectY[2] += radCoeff[4] * gaussRand[pId * 6 + 2]; 
        vectY[3] += radCoeff[4] * gaussRand[pId * 6 + 3]; 
        // vectY[4] += radCoeff[5] * gaussRand[pId * 6 + 4]; 
        vectY[5] += radCoeff[5] * gaussRand[pId * 6 + 5]; 

        // double gaussRand = curand_normal_double(state);
        // printf("%15.8f\n",gaussRand[pId * 6]);

        #pragma unroll
        for(int i=0;i<6;i++)
        {
            vectX[i] = 0;
            #pragma unroll
            for(int j=0;j<6;j++)
            {
                vectX[i] += radMatrixBRH[6*i +j + 36] * vectY[j];
            }
        }
         
        // transfer to acc cord from norm cord
        #pragma unroll
        for(int j=0;j<6;j++) partCord[pId * 6 + j] = vectX[j];
        //----------syn end-----------------------------------------------------------------------------------------------
        
    }

}

void GPU_PartiOneTurnTransferAndSynRad(int totalPartiNum, double *partCord, int oneTurnMatrixParaNum, double *ringPara, double *radMatrixBRH, double *radCoeff)
{

    int partCordMem              = totalPartiNum * 6    * sizeof(double);
    int ringParaMem              = oneTurnMatrixParaNum * sizeof(double);
    int radMatrixBRHMem          = 72                   * sizeof(double);
    int radCeoffMem              = 6                    * sizeof(double); 
    
    double *d_partCord,*d_ringPara,*d_radMatrixBRH,*d_radCoeff,*d_gaussRand ;
    
    cudaMalloc((void**) &d_partCord,     partCordMem);
    cudaMalloc((void**) &d_gaussRand,    partCordMem);
    cudaMalloc((void**) &d_ringPara,     ringParaMem);
    cudaMalloc((void**) &d_radMatrixBRH, radMatrixBRHMem);
    cudaMalloc((void**) &d_radCoeff,     radCeoffMem);


    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> dx{0,1};

    double gaussRand[totalPartiNum*6];
    for(int i=0;i<totalPartiNum*6;i++) gaussRand[i    ] = dx(gen);
    for(int i=0;i<totalPartiNum;  i++) gaussRand[6*i+4] = 0;
    

    cudaMemcpy(d_partCord,     partCord,      partCordMem,              cudaMemcpyHostToDevice);
    cudaMemcpy(d_gaussRand,    gaussRand,     partCordMem,              cudaMemcpyHostToDevice);
    cudaMemcpy(d_ringPara,     ringPara,      ringParaMem,              cudaMemcpyHostToDevice);
    cudaMemcpy(d_radMatrixBRH, radMatrixBRH,  radMatrixBRHMem,          cudaMemcpyHostToDevice);
    cudaMemcpy(d_radCoeff,     radCoeff,      radCeoffMem,              cudaMemcpyHostToDevice);  
    
    // each particle is locate at each thread on GPU. 
    int blockDim;
    if(totalPartiNum>1000) 
    {
        blockDim =128;
    }
    else
    {
        int count = 0;
        int temp=totalPartiNum;
        while(temp/2!=0)
        {
            count++;
            temp /=2;    
        }
        blockDim = pow(2,count);
    }
    dim3 block(blockDim);
    dim3 grid((totalPartiNum + block.x -1)/ block.x);

    GPU_OneTurnMapAndSynRad<<<grid, block>>>(totalPartiNum,d_partCord,d_ringPara,d_radMatrixBRH,d_radCoeff,d_gaussRand);


    cudaMemcpy(partCord, d_partCord, partCordMem, cudaMemcpyDeviceToHost);
    cudaFree(d_partCord);
    cudaFree(d_ringPara);

}