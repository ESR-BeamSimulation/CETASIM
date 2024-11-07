//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once    

#include "Global.h"
#include "PIC2D.h"
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

using v1i= vector<int> ;
using v1d= vector<double> ;
using v2d= vector<vector<double> > ;
using v3d= vector<vector<vector<double> > > ;

PIC2D::PIC2D()
{

}

PIC2D::~PIC2D()
{
    fftw_free(in_xy);
    fftw_free(out_xy);
    fftw_destroy_plan(rho_X_Y_To_Rho_Kx_Ky);
    fftw_destroy_plan(phi_Kx_Ky_To_Phi_X_Y);

}

void PIC2D::InitialPIC2D(vector<int> meshNum)
{
    numberOfGrid = meshNum;
    int nx = numberOfGrid[0];
    int ny = numberOfGrid[1];
  
    meshx = v1d(nx);
    meshy = v1d(ny);

    rho = v2d(nx, v1d(ny) );
    phi = v2d(nx, v1d(ny) );
    ex  = v2d(nx, v1d(ny) );
    ey  = v2d(nx, v1d(ny) );

    in_xy  = (double*)             fftw_malloc(sizeof(double)      * nx * ny );
    out_xy = (double*)             fftw_malloc(sizeof(double)      * nx * ny );

    int numofgridx = nx;
    int numofgridy = ny;
    int numofgridz = 1;
    
    // 1
    int n_xy[]={numofgridx,numofgridy};
    rho_X_Y_To_Rho_Kx_Ky = fftw_plan_many_r2r(2, n_xy, numofgridz,
    in_xy  ,n_xy,
    numofgridz, 1,
    out_xy ,n_xy,
    numofgridz, 1,
    kind_xy_forward, FFTW_MEASURE);
    //2
    // fftw_plan_with_nthreads(omp_get_max_threads());
    phi_Kx_Ky_To_Phi_X_Y = fftw_plan_many_r2r(2, n_xy, numofgridz,
    out_xy,n_xy,
    numofgridz, 1,
    in_xy,n_xy ,
    numofgridz, 1,
    kind_xy_backward, FFTW_MEASURE);
}

void PIC2D::Set2DMesh(double rmsXY[2], double beamCen[2])
{
    int dims = 2;
    vector<double> meshLength(dims,0.E0);
    vector<double> meshStart(dims,0.E0);
    
    for(int plane=0; plane<dims;plane++)
    {
        beamrmssize[plane]  = rmsXY[plane];
        meshCen[plane]      = beamCen[plane];
        meshLength[plane]   = beamrmssize[plane] * 2 * range;
    }

    double meshLengthMax = (meshLength[0] > meshLength[1]) ?  meshLength[0] : meshLength[1];
    
    dv = 1.E0;
    for (int plane=0; plane<dims;plane++)
    {  
        meshLength[plane] = meshLengthMax;
        meshWidth[plane]  = meshLength[plane] / (numberOfGrid[plane] - 1);
        meshStart[plane]  = beamCen[plane] - meshLength[plane] / 2; 

        dv *= meshWidth[plane];

        for(int i=0; i<numberOfGrid[plane]; i++) 
        {
            if(plane==0) meshx[i] = meshStart[plane] + meshWidth[plane] * i;
            if(plane==1) meshy[i] = meshStart[plane] + meshWidth[plane] * i;
        } 
    }
}


void PIC2D::Set2DMesh(vector<vector<double>>  &particles)
{
    int dims = 2; 
    vector<double> beamCen(dims,0.E0);
    vector<double> beamRms(dims,0.E0);
    vector<double> meshLength(dims,0.E0);
    vector<double> meshStart(dims,0.E0);

    dv = 1.E0;

    for (int plane=0; plane<dims;plane++)
    {   
        vector<double> pos = particles[plane];
        beamCen[plane] = accumulate(pos.begin(), pos.end(), 0.E0) / pos.size();
        beamRms[plane] = 0.E0;
        for(int i = 0; i<pos.size();i++ )
        {
            beamRms[plane] += pow(pos[i] - beamCen[plane] ,2);
        }
        beamRms[plane] =  sqrt(beamRms[plane] / pos.size() );
        meshCen[plane] = beamCen[plane];
        meshLength[plane] = beamRms[plane] * 2 * range;
    }
    // dynamical mesh generates good result only when  meshLength[0]~meshLength[1]. 
    // if not Rx>>Ry, then dynamical meshlength meshLength[0] > meshLength[1], 
    // simply set meshLength[0] = meshLength[1] which gives right physics 
    if(beamRms[0]==0 || beamRms[1]==0)
    {
        cerr<<"2D PIC mesh length, x:"<< beamRms[0]<<", y:"<<beamRms[1]<<endl;
    }


    double meshLengthMax = (meshLength[0] > meshLength[1]) ?  meshLength[0] : meshLength[1]; 
    
    for (int plane=0; plane<dims;plane++)
    {  
        meshLength[plane] = meshLengthMax;
        meshWidth[plane]  = meshLength[plane] / (numberOfGrid[plane] - 1);
        meshStart[plane]  = beamCen[plane] - meshLength[plane] / 2; 

        dv *= meshWidth[plane];

        for(int i=0; i<numberOfGrid[plane]; i++) 
        {
            if(plane==0) meshx[i] = meshStart[plane] + meshWidth[plane] * i;
            if(plane==1) meshy[i] = meshStart[plane] + meshWidth[plane] * i;
        } 
        beamrmssize[plane] = beamRms[plane]; 
    }
}

void PIC2D::Set2DRho (vector<vector<double>>  &particles, vector<double> charge)
{
    // initialize vecotor size accroding to the input parameters
    partNum =  particles[0].size();
    
    partMeshIndX = v1i(partNum);
    partMeshIndY = v1i(partNum); 
    isPartOutMesh = v1i(partNum,0); 
    weigh = v2d(partNum,v1d(4,0.E0));
    partExEyField = v2d(2,v1d(partNum,0.E0));          //(Ex(v1d),Ey(v1d))
    
    // charge represents the line density  [C/m]

    for(int i=0;i<partNum;i++)
    {
        partMeshIndX[i] = floor( (particles[0][i] - meshx[0] ) / meshWidth[0] ) ;
        partMeshIndY[i] = floor( (particles[1][i] - meshy[0] ) / meshWidth[1] ) ;
    }
    
    // get the linear weighting factor from particle to mesh. 
    for (int i=0;i<partNum;i++)
    {
        bool outMeshFlagX, outMeshFlagY, outMeshFlagZ, outMeshFlag; 

        outMeshFlagX = partMeshIndX[i] >= numberOfGrid[0] - 1 || partMeshIndX[i] < 0;  
        outMeshFlagY = partMeshIndY[i] >= numberOfGrid[1] - 1 || partMeshIndY[i] < 0;
        outMeshFlag  = outMeshFlagX || outMeshFlagY ;

        if(outMeshFlag) 
        {
            isPartOutMesh[i] = 1;
            continue;
        }

        double sum = 0;
        for(int j=0;j<4; j++)
        {
            int idx, idy;
            idx = (j%2==0) ? 1 : 0;
            idy = (j%4< 2) ? 1 : 0;

            weigh[i][j] =   
            abs( particles[0][i] - meshx[ partMeshIndX[i] + idx ]) *
            abs( particles[1][i] - meshy[ partMeshIndY[i] + idy ])  / dv;

            sum += weigh[i][j];
        }
        if(abs(sum-1)>1.E-9)  
        {
            cerr<<"weighing is wrong:"<< sum - 1<<endl;
            exit(0);
        }
    }
    
    // get the charge density on mesh
    for(int x=0;x<numberOfGrid[0];x++)
    {
        for(int y=0;y<numberOfGrid[1];y++)
        {
            rho[x][y] =0.E0;
        }
    }   
        
    for (int i=0;i<partNum;i++)
    {
        if(isPartOutMesh[i]==1) continue; 
        int idx  =  partMeshIndX[i];
        int idy  =  partMeshIndY[i];
        rho[idx  ][idy  ] += weigh[i][0] * charge[i] / dv;    
        rho[idx+1][idy  ] += weigh[i][1] * charge[i] / dv;
        rho[idx  ][idy+1] += weigh[i][2] * charge[i] / dv;
        rho[idx+1][idy+1] += weigh[i][3] * charge[i] / dv; //   C/m^3
    }
   

    // comment: With un-symmetric beam (electron beam with non 1 coupling factor)
    // The PIC solver does not produce the the electric field well. 
    // In below, it is an example from a uniform profile in real space.
    // analytically with a 2d beam.  lineDensity = totCharge / circ [C/m]
    // density lambda =  lineDensity/ pi / (4 a b) [C/m^3]
    // within the profile: Ex[x] =  lambda / epsilon0 b /(a+b) x;  Ey[y] =  lambda / epsilon0 a /(a+b) y;  
    // KV beam distributin for benchmark 
    // double a = beamrmssize[0] * 2 ;
    // double b = beamrmssize[1] * 2 ;
    // double x,y;
    // for(int i=0;i<numberOfGrid[0];i++)
    // {
    //     for(int j=0;j<numberOfGrid[1];j++)
    //     {
    //         x = meshx[0] + i * meshWidth[0];
    //         y = meshy[0] + j * meshWidth[1];
    //         if ( pow(x,2) / pow(a,2)  + pow(y,2) / pow(b,2) < 1 ) 
    //         {   
    //             rho[i][j] = 3.33564E-12 / PI / a / b ;
    //         }
    //         else 
    //         {
    //             rho[i][j] = 0.E0;
    //         }
            
    //     }
    // }
}


void PIC2D::Set2DPhi()
{
    int nx = numberOfGrid[0];
    int ny = numberOfGrid[1];

    // (x,y) real space to (kx,ky)
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny;j++)
        {
            phi[i][j] =0.E0;
            in_xy[i * ny  + j ] = rho[i][j] ;                       //   C/m^3
        }
    }   

    fftw_execute(rho_X_Y_To_Rho_Kx_Ky);                             //  in_xy->out_xy

    double K_rho2phi=0, knx=0, kny=0;
    int coun445;

    for(int i=0; i<nx; ++i)
    {
        for(int j=0; j<ny; ++j)
        {  
            knx=PI / 2 * (i+1) / (nx + 1);
            kny=PI / 2 * (j+1) / (ny + 1);
            
            // The SK function is uniquely decided
            K_rho2phi = pow(2*sin(knx) / meshWidth[0], 2)
                      + pow(2*sin(kny) / meshWidth[1], 2);            // unit is 1/m^2   

            coun445 = i * ny  + j ;
            out_xy[coun445] = out_xy[coun445] / K_rho2phi /Epsilon;   // unit is  (C/m^3) / (1/m^2) /(C/m*V) = V             
        }
    }

    fftw_execute(phi_Kx_Ky_To_Phi_X_Y);                              //  out_xy -> in_xy

    for(int i=0; i<nx; ++i)
    {
        for(int j=0; j<ny; ++j)
        {
            coun445 = i * ny  + j ;
            phi[i][j] = in_xy[coun445] / (4 * (nx + 1) * (ny + 1)) ;        // [V]
        }
    }

}


void PIC2D::Set2DEField()
{
    int nx = numberOfGrid[0];
    int ny = numberOfGrid[1];
    double phix1 = 0.E0, phix2 = 0.E0;
    double phiy1 = 0.E0, phiy2 = 0.E0;

    for(int i=0;i<nx;++i)
    {
        for(int j=0;j<ny;++j)
        {                
            phix1 = (i==0   ) ? 0 :  phi[i-1][j  ];
            phix2 = (i==nx-1) ? 0 :  phi[i+1][j  ];
            
            phiy1 = (j==0   ) ? 0 :  phi[i  ][j-1];
            phiy2 = (j==ny-1) ? 0 :  phi[i  ][j+1];
            
            ex[i][j] =  -(phix2 - phix1)  / (2 * meshWidth[0]);      // [V/m]
            ey[i][j] =  -(phiy2 - phiy1)  / (2 * meshWidth[1]);      // [V/m]            
        }
    }

    // within the profile: Ex[x] =  lambda / epsilon0 b /(a+b) x;  Ey[y] =  lambda / epsilon0 a /(a+b) y;  
//     ofstream fout("test.sdds");
//     double a = beamrmssize[0] * 2 ;
//     double b = beamrmssize[1] * 2 ;
//     double x,y,tempx,tempy;
//     for(int i=0;i<numberOfGrid[0];i++)
//     {
//         for(int j=0;j<numberOfGrid[1];j++)
//         {
//             x = meshx[0] + i * meshWidth[0];
//             y = meshy[0] + j * meshWidth[1];
//             if ( pow(x,2) / pow(a,2)  + pow(y,2) / pow(b,2) < 1 ) 
//             {   
//                 tempx =  x /a / (a+b);
//                 tempy =  y /b / (a+b);
//             }
//             else 
//             {
//                 tempx = 0.E0;
//                 tempy = 0.E0;
//             }

//             fout<<setw(15)<< x 
//                 <<setw(15)<< y
//                 <<setw(15)<< rho[i][j]
//                 <<setw(15)<< phi[i][j]
//                 <<setw(15)<< ex[i][j]
//                 <<setw(15)<< ey[i][j]
//                 <<setw(15)<< tempx   
//                 <<setw(15)<< tempy   
//                 <<endl;
//         }
//     }  
//     cout<<__LINE__<<__FILE__<<endl;
//     getchar();
}


void PIC2D::SetPartSCField()
{
    int idx, idy;
    for(int i=0;i<partNum;++i)
    {
        if(isPartOutMesh[i]==1)
        {
            partExEyField[0][i] = 0.E0 ;
            partExEyField[1][i] = 0.E0 ;
        }
        else
        {
            idx  =  partMeshIndX[i];
            idy  =  partMeshIndY[i];
        
            partExEyField[0][i] = ex[idx  ][idy  ]     * weigh[i][0]
                                + ex[idx+1][idy  ]     * weigh[i][1]
                                + ex[idx  ][idy+1]     * weigh[i][2]
                                + ex[idx+1][idy+1]     * weigh[i][3];
                    
             partExEyField[1][i]= ey[idx  ][idy  ]     * weigh[i][0]
                                + ey[idx+1][idy  ]     * weigh[i][1]
                                + ey[idx  ][idy+1]     * weigh[i][2]
                                + ey[idx+1][idy+1]     * weigh[i][3];            
        }
    }
}

void PIC2D::UpdateTransverseMomentum(vector<vector<double>>  &particles, const double ds)
{
    double gamma0 = electronBeamEnergy / ElectronMassEV; 
    double beta0  = sqrt(1 - 1 / pow(gamma0, 2));
    double p0     = beta0 * gamma0;
    double px, py, pz, p, gamma, beta, factorT, factorL; 
    double deltaPx,deltaPy,deltaPz; 
    
    // ofstream fout("test.sdds");

    for(int i=0;i<partNum;++i)
    {
        pz = (1 + particles[5][i]) * p0;
        px = pz * particles[1][i];
        py = pz * particles[3][i];
        p = sqrt(pow(pz,2) + pow(px,2) + pow(py,2) );
        gamma = sqrt(1 + pow(p,2) );
        beta  = p / gamma;
  
        deltaPx =  partExEyField[0][i] / pow(gamma,2) * ElectronCharge * ds / (beta * CLight) / (ElectronMass * p * CLight );
        deltaPy =  partExEyField[1][i] / pow(gamma,2) * ElectronCharge * ds / (beta * CLight) / (ElectronMass * p * CLight );

        particles[3][i] += deltaPx;
        particles[4][i] += deltaPy;

        // fout<<setw(15)<< particles[0][i] 
        //     <<setw(15)<< particles[1][i]
        //     <<setw(15)<< deltaPx
        //     <<setw(15)<< deltaPy
        //     <<endl; 

    }
    
    // cout<<__LINE__<<__FILE__<<endl;
    // getchar();

}





vector<vector<double> > PIC2D::GetPartSCField(vector<vector<double>> &particles)
{    
    //(1) get the weighing localy according to the input particle dis...
    int np = particles[0].size();
    vector<vector<double> > partEField(2,vector<double>(np,0.E0));
    vector<int> partMeshIndX =  vector<int> (np, 0); 
    vector<int> partMeshIndY =  vector<int> (np, 0);
    vector<int> isPartOutMesh = vector<int> (np, 0);
    v2d weigh = v2d(np,v1d(4,0.E0));

    for(int i=0;i<np;i++)
    {
        partMeshIndX[i] = floor( (particles[0][i] - meshx[0] ) / meshWidth[0] ) ;
        partMeshIndY[i] = floor( (particles[1][i] - meshy[0] ) / meshWidth[1] ) ;
    }
    
    // get the linear weighting factor from particle to mesh. 
    for (int i=0;i<np;i++)
    {
        bool outMeshFlagX, outMeshFlagY, outMeshFlag; 

        outMeshFlagX = partMeshIndX[i] >= numberOfGrid[0] - 1 || partMeshIndX[i] < 0;  
        outMeshFlagY = partMeshIndY[i] >= numberOfGrid[1] - 1 || partMeshIndY[i] < 0;
        outMeshFlag  = outMeshFlagX || outMeshFlagY ;

        if(outMeshFlag) 
        {
            isPartOutMesh[i] = 1;
            continue;
        }

        double sum = 0;
        int idx, idy;
        for(int j=0;j<4; j++)
        {
            idx = (j%2==0) ? 1 : 0;
            idy = (j%4< 2) ? 1 : 0;

            weigh[i][j] =   
            abs( particles[0][i] - meshx[ partMeshIndX[i] + idx ]) *
            abs( particles[1][i] - meshy[ partMeshIndY[i] + idy ])  /dv;

            sum += weigh[i][j];
        }
        if(abs(sum-1)>1.E-9)  
        {
            cerr<<"weighing is wrong:"<< sum - 1<<endl;
            exit(0);
        }
    }

    //（2) accrording to the weighing obtained get the efield
    int idx, idy;
    for(int i=0;i<np;++i)
    {
        if(isPartOutMesh[i]==1)
        {
            partEField[0][i] = 0.E0 ;
            partEField[1][i] = 0.E0 ;
        }
        else
        {
            idx  =  partMeshIndX[i];
            idy  =  partMeshIndY[i];
        
            partEField[0][i]    = ex[idx  ][idy  ]     * weigh[i][0]
                                + ex[idx+1][idy  ]     * weigh[i][1]
                                + ex[idx  ][idy+1]     * weigh[i][2]
                                + ex[idx+1][idy+1]     * weigh[i][3];
                    
             partEField[1][i]   = ey[idx  ][idy  ]     * weigh[i][0]
                                + ey[idx+1][idy  ]     * weigh[i][1]
                                + ey[idx  ][idy+1]     * weigh[i][2]
                                + ey[idx+1][idy+1]     * weigh[i][3];
            
        }
    }    

    return partEField;

}