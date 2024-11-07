//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once    

#include "Global.h"
#include "PIC3D.h"
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

PIC3D::PIC3D()
{
};

PIC3D::~PIC3D()
{

    fftw_free(in_xy);
    fftw_free(out_xy);
    fftw_free(in);
    fftw_free(out);
    fftw_free(outf);
    fftw_free(in_z);
    fftw_free(out_z);
    fftw_free(outf_z);

    fftw_free(in_xyz);
    fftw_free(out_xyz);

    fftw_destroy_plan(rho_X_Y_Z_To_Rho_KX_KY_Z);
    fftw_destroy_plan(phi_KX_KY_Z_To_Phi_X_Y_Z);
    fftw_destroy_plan(rho_KX_KY_Z_To_Rho_KX_KY_KZ);
    fftw_destroy_plan(phi_KX_KY_KZ_To_Phi_KX_KY_Z);

    fftw_destroy_plan(rho_X_Y_Z_To_Rho_KX_KY_KZ);
    fftw_destroy_plan(phi_KX_KY_KZ_To_Phi_X_Y_Z);

    delete slicedBunchPIC2D;

}

void PIC3D::InitialSC3D(const ReadInputSettings &inputParameter)
{
    electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    for(int i=0;i<3;++i)
    {
        numberOfGrid[i] = inputParameter.ringRun->scMeshNum[i];  
    }
    
    int nx = numberOfGrid[0];
    int ny = numberOfGrid[1];
    int nz = numberOfGrid[2];
  
    meshx = v1d(nx);
    meshy = v1d(ny);
    meshz = v1d(nz);

    rho = v3d(nx, v2d(ny, v1d(nz) ) );
    phi = v3d(nx, v2d(ny, v1d(nz) ) );
    ex  = v3d(nx, v2d(ny, v1d(nz) ) );
    ey  = v3d(nx, v2d(ny, v1d(nz) ) );
    ez  = v3d(nx, v2d(ny, v1d(nz) ) );


    in     = (double*)             fftw_malloc(sizeof(double)      * nx * ny * nz);
    out    = (double*)             fftw_malloc(sizeof(double)      * nx * ny * nz);
    outf   = (double*)             fftw_malloc(sizeof(double)      * nx * ny * nz);
    in_xy  = (double*)             fftw_malloc(sizeof(double)      * nx * ny * nz);
    out_xy = (double*)             fftw_malloc(sizeof(double)      * nx * ny * nz);
    in_z   = (fftw_complex*)       fftw_malloc(sizeof(fftw_complex)* nx * ny * nz);
    out_z  = (fftw_complex*)       fftw_malloc(sizeof(fftw_complex)* nx * ny * nz);
    outf_z = (fftw_complex*)       fftw_malloc(sizeof(fftw_complex)* nx * ny * nz);

    in_xyz     = (double*)          fftw_malloc(sizeof(double)      * nx * ny * nz);
    out_xyz    = (double*)          fftw_malloc(sizeof(double)      * nx * ny * nz);


    int numofgridx = nx;
    int numofgridy = ny;
    int numofgridz = nz; 

    // 1
    int n_xy[]={numofgridx,numofgridy};
    rho_X_Y_Z_To_Rho_KX_KY_Z = fftw_plan_many_r2r(2, n_xy, numofgridz,
    in_xy  ,n_xy,
    numofgridz, 1,
    out_xy ,n_xy,
    numofgridz, 1,
    kind_xy_forward, FFTW_MEASURE);
    //2
    // fftw_plan_with_nthreads(omp_get_max_threads());
    int n_z[]={numofgridz};
    rho_KX_KY_Z_To_Rho_KX_KY_KZ = fftw_plan_many_dft(1, n_z, numofgridx*numofgridy,
    in_z , n_z ,
    1,numofgridz,
    out_z, n_z ,
    1,numofgridz,
    FFTW_FORWARD, FFTW_MEASURE);
    //3
    // fftw_plan_with_nthreads(omp_get_max_threads());
    phi_KX_KY_KZ_To_Phi_KX_KY_Z = fftw_plan_many_dft(1, n_z, numofgridx*numofgridy,
    out_z  , n_z ,
    1,numofgridz,
    outf_z, n_z ,
    1,numofgridz,
    FFTW_BACKWARD, FFTW_MEASURE);
    //4
    // fftw_plan_with_nthreads(omp_get_max_threads());
    phi_KX_KY_Z_To_Phi_X_Y_Z = fftw_plan_many_r2r(2, n_xy, numofgridz,
    out_xy,n_xy,
    numofgridz, 1,
    in_xy,n_xy ,
    numofgridz, 1,
    kind_xy_backward, FFTW_MEASURE);


    //test 3d sine fft  
    int n_xyz[]={numofgridx,numofgridy,numofgridz};
    rho_X_Y_Z_To_Rho_KX_KY_KZ = fftw_plan_r2r(3, n_xyz,  in_xyz, out_xyz, kind_xyz_forward,  FFTW_MEASURE);
    phi_KX_KY_KZ_To_Phi_X_Y_Z = fftw_plan_r2r(3, n_xyz, out_xyz,  in_xyz, kind_xyz_backward, FFTW_MEASURE);

    vector<int> meshNum = vector<int>{numofgridx,numofgridy};
    slicedBunchPIC2D =  new PIC2D();
    slicedBunchPIC2D->InitialPIC2D(meshNum);
    slicedBunchPIC2D->electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;
    slicedBunchPIC2D->range = range;
}


void PIC3D::Set3DMesh(vector<vector<double>> &particles)
{
    // particles -> [3,np] dimension, store the info of position (x, y, z)
    int dims = 3; 
    vector<double> beamCen(dims,0.E0);
    beamRms = vector<double> (dims,0.E0);
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


    // set mesh length as the same value -- maximum ones.
    
    double meshLengthMax = (meshLength[0] > meshLength[1]) ?  meshLength[0] : meshLength[1]; 
    meshLengthMax = (meshLengthMax > meshLength[2]) ?  meshLengthMax : meshLength[2]; 
    meshLength[0] = meshLengthMax; 
    meshLength[1] = meshLengthMax;
    meshLength[2] = meshLengthMax; 


    for (int plane=0; plane<dims;plane++)
    {
        meshStart[plane]  = beamCen[plane] - meshLength[plane] / 2; 
        meshWidth[plane]  = meshLength[plane] / (numberOfGrid[plane] - 1);

        dv *= meshWidth[plane];

        for(int i=0; i<numberOfGrid[plane]; i++) 
        {
            if(plane==0) meshx[i] = meshStart[plane] + meshWidth[plane] * i;
            if(plane==1) meshy[i] = meshStart[plane] + meshWidth[plane] * i;
            if(plane==2) meshz[i] = meshStart[plane] + meshWidth[plane] * i;
        } 
    }

    partNum =  particles[0].size();
    partMeshIndX = v1i(partNum);
    partMeshIndY = v1i(partNum); 
    partMeshIndZ = v1i(partNum);
    isPartOutMesh = v1i(partNum,0); 
    weigh = v2d(partNum,v1d(8,0.E0));
    partEField = v2d(3,v1d(partNum,0.E0));
    
}

void PIC3D::Set3DRho(vector<vector<double>> &particles, double macroCharge)
{
    // macroCharge is charge per macro-partilce has with a unit [C]
    for(int i=0;i<partNum;i++)
    {
        partMeshIndX[i] = floor( (particles[0][i] - meshx[0] ) / meshWidth[0] ) ;
        partMeshIndY[i] = floor( (particles[1][i] - meshy[0] ) / meshWidth[1] ) ;
        partMeshIndZ[i] = floor( (particles[2][i] - meshz[0] ) / meshWidth[2] ) ;
    }

   
    // get the linear weighting factor from particle to mesh. 
    for (int i=0;i<partNum;i++)
    {
        bool outMeshFlagX, outMeshFlagY, outMeshFlagZ, outMeshFlag; 

        outMeshFlagX = (partMeshIndX[i] >= numberOfGrid[0] - 1) || (partMeshIndX[i] < 0);  
        outMeshFlagY = (partMeshIndY[i] >= numberOfGrid[1] - 1) || (partMeshIndY[i] < 0);
        outMeshFlagZ = (partMeshIndZ[i] >= numberOfGrid[2] - 1) || (partMeshIndZ[i] < 0);  
        outMeshFlag  = outMeshFlagX || outMeshFlagY || outMeshFlagZ ;

        if(outMeshFlag) 
        {
            isPartOutMesh[i] = 1;
            continue;
        }

        double sum = 0;
        for(int j=0;j<8; j++)
        {
            int idx, idy, idz;
            idx = (j%2==0) ? 1 : 0;
            idy = (j%4< 2) ? 1 : 0;
            idz = (j <= 3) ? 1 : 0;

            weigh[i][j] =   
            abs( particles[0][i] - meshx[ partMeshIndX[i] + idx ]) *
            abs( particles[1][i] - meshy[ partMeshIndY[i] + idy ]) *
            abs( particles[2][i] - meshz[ partMeshIndZ[i] + idz ]) / dv;   
            
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
            for(int z=0;z<numberOfGrid[2];z++)
            {
                 rho[x][y][z] = 0.E0;
            }
        }
    }    
    
   
    for (int i=0;i<partNum;i++)
    {
        int idx  =  partMeshIndX[i];
        int idy  =  partMeshIndY[i];
        int idz  =  partMeshIndZ[i];
        
        rho[idx  ][idy  ][idz  ] += weigh[i][0] * macroCharge / dv;    
        rho[idx+1][idy  ][idz  ] += weigh[i][1] * macroCharge / dv;
        rho[idx  ][idy+1][idz  ] += weigh[i][2] * macroCharge / dv;
        rho[idx+1][idy+1][idz  ] += weigh[i][3] * macroCharge / dv;
        rho[idx  ][idy  ][idz+1] += weigh[i][4] * macroCharge / dv;
        rho[idx+1][idy  ][idz+1] += weigh[i][5] * macroCharge / dv;
        rho[idx  ][idy+1][idz+1] += weigh[i][6] * macroCharge / dv;
        rho[idx+1][idy+1][idz+1] += weigh[i][7] * macroCharge / dv;   //   [C/m^3]

    }
    
    // benechmark -- with a unifrom beam distribution
    // double sum = 0;
    // double a, b, c;
    // a = beamRms[0]; 
    // b = beamRms[1]; 
    // c = beamRms[2];
    // cout<<a<<endl;
    // cout<<b<<endl;
    // cout<<c<<endl;
    // cout<<meshWidth[0]<<endl;
    // cout<<meshWidth[1]<<endl;
    // cout<<meshWidth[2]<<endl;
    
    
    // for(int i=0;i<numberOfGrid[0];i++)
    // {
    //     for(int j=0;j<numberOfGrid[1];j++)
    //     {
    //         for(int k=0;k<numberOfGrid[2];k++)
    //         {
    //             double x = meshx[0] + i * meshWidth[0];
    //             double y = meshy[0] + j * meshWidth[1];
    //             double z = meshz[0] + k * meshWidth[2];

    //             if ( pow(x,2) / pow(a,2)  + pow(y,2) / pow(b,2) + pow(z,2) / pow(c,2) < 1 ) 
    //             {   
    //                 rho[i][j][k] = 1. / (4 / 3. * PI * a * b * c);
    //             }
    //             else 
    //             {
    //                 rho[i][j][k] = 0.E0;
    //             }

    //         }
    //     }
    // }    
}

void PIC3D::Set3DPhi()
{
    int nx = numberOfGrid[0];
    int ny = numberOfGrid[1];
    int nz = numberOfGrid[2];

    // (x,y) real space to (kx,ky)
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny;j++)
        {
            for(int k=0; k<nz; k++)
            {
                phi[i][j][k] = 0.E0;
                in_xy[i * ny * nz + j * nz + k] = rho[i][j][k] ;   //   C/m^3
            }
        }
    }    

    fftw_execute(rho_X_Y_Z_To_Rho_KX_KY_Z);                         //  in_xy->out_xy

    int coun365;
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nz;k++)
            {
                coun365 = i * ny * nz + j * nz + k;
                in_z[coun365][0]=out_xy[coun365];                  //   out_xy-> in_z 
                in_z[coun365][1]=0;
            }
        }
    }

    fftw_execute(rho_KX_KY_Z_To_Rho_KX_KY_KZ);                   //   in_z-> out_z

    double K_rho2phi=0, knx=0, kny=0, knz=0;
    int coun445;

    for(int i=0; i<nx; ++i)
    {
        for(int j=0; j<ny; ++j)
        {
            for(int k=0; k<nz; ++k)
            {
                knx=PI / 2 * (i+1) / (nx + 1);
                kny=PI / 2 * (j+1) / (ny + 1);
                knz=PI     * (k)   / (nz);              // chech the formular used in FFTW mannul in function FFT,  
                
                // The SK function is uniquely decided
                K_rho2phi=   pow(2*sin(knx) / meshWidth[0], 2)
                            +pow(2*sin(kny) / meshWidth[1], 2)
                            +pow(2*sin(knz) / meshWidth[2], 2);   // unit is 1/m^2

                coun445 = i * ny * nz + j * nz + k;
                out_z[coun445][0] = out_z[coun445][0] /K_rho2phi /Epsilon;   // unit is  (C/m^3) / (1/m^2) /(C/m*V)=V                
                out_z[coun445][1] = out_z[coun445][1] /K_rho2phi /Epsilon;
            }
        }
    }



    fftw_execute(phi_KX_KY_KZ_To_Phi_KX_KY_Z);                           // out_z  -> outf_z

    for(int i=0; i<nx; ++i)
    {
        for(int j=0; j<ny; ++j)
        {
            for(int k=0; k<nz; ++k)
            {
                coun365 = i * ny * nz + j * nz + k;
                out_xy[coun365] = outf_z[coun365][0];                     // outf_z -> out_xy
            }
        }
    }

    fftw_execute(phi_KX_KY_Z_To_Phi_X_Y_Z);                              //  out_xy -> in_xy

    for(int i=0; i<nx; ++i)
    {
        for(int j=0; j<ny; ++j)
        {
            for(int k=0; k<nz; ++k)
            {
                coun365 = i * ny * nz + j * nz + k;
                phi[i][j][k] = in_xy[coun365] / (4 * (nx + 1) * (ny + 1) * nz ) ;        // [V]
            }
        }
    }
}


void PIC3D::Set3DPhi1()
{
    int nx = numberOfGrid[0];
    int ny = numberOfGrid[1];
    int nz = numberOfGrid[2];

    // (x,y，z) real space to (kx,ky,kz)
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny;j++)
        {
            for(int k=0; k<nz; k++)
            {
                phi[i][j][k] = 0.E0;
                in_xyz[i * ny * nz + j * nz + k] = rho[i][j][k] ;   //   C/m^3
            }
        }
    }    

    fftw_execute(rho_X_Y_Z_To_Rho_KX_KY_KZ);                       //  in_xyz->out_xyz

    double K_rho2phi=0, knx=0, kny=0, knz=0;
    int coun445;

    for(int i=0; i<nx; ++i)
    {
        for(int j=0; j<ny; ++j)
        {
            for(int k=0; k<nz; ++k)
            {
                knx=PI / 2 * (i+1) / (nx + 1);
                kny=PI / 2 * (j+1) / (ny + 1);
                knz=PI / 2 * (k+1) / (nz + 1);               // chech the formular used in FFTW mannul in function FFT,  
                
                // The SK function is uniquely decided
                K_rho2phi=   pow(2*sin(knx) / meshWidth[0], 2)
                            +pow(2*sin(kny) / meshWidth[1], 2)
                            +pow(2*sin(knz) / meshWidth[2], 2);   // unit is 1/m^2

                coun445 = i * ny * nz + j * nz + k;
                out_xyz[coun445] = out_xyz[coun445] /K_rho2phi / Epsilon;   // unit is  (C/m^3) / (1/m^2) /(C/m*V)=V                
            }
        }
    }

    fftw_execute(phi_KX_KY_KZ_To_Phi_X_Y_Z);                              //  out_xyz -> in_xyz
    int coun365;
    for(int i=0; i<nx; ++i)
    {
        for(int j=0; j<ny; ++j)
        {
            for(int k=0; k<nz; ++k)
            {
                coun365 = i * ny * nz + j * nz + k;
                phi[i][j][k] = in_xyz[coun365] / (8 * (nx + 1) * (ny + 1) * (nz + 1) ) ;        // [V]
            }
        }
    }

}


void PIC3D::Set3DEField()
{
    int nx = numberOfGrid[0];
    int ny = numberOfGrid[1];
    int nz = numberOfGrid[2];
    double phix1 = 0.E0, phix2 = 0.E0;
    double phiy1 = 0.E0, phiy2 = 0.E0;
    double phiz1 = 0.E0, phiz2 = 0.E0;

    ofstream fout("test.sdds");
    for(int i=0;i<nx;++i)
    {
        for(int j=0;j<ny;++j)
        {
            for(int k=0;k<nz;++k)
            {
                
                phix1 = (i==0   ) ? 0 :  phi[i-1][j  ][k  ];
                phix2 = (i==nx-1) ? 0 :  phi[i+1][j  ][k  ];
                
                phiy1 = (j==0   ) ? 0 :  phi[i  ][j-1][k  ];
                phiy2 = (j==ny-1) ? 0 :  phi[i  ][j+1][k  ];

                phiz1 = (k==0   ) ? 0 :  phi[i  ][j  ][k-1];
                phiz2 = (k==nz-1) ? 0 :  phi[i  ][j  ][k+1];
                
                ex[i][j][k] =  -(phix2 - phix1)  / (2 * meshWidth[0]);  //* beamRms[2] * 2 / 2304;      // [V/m]
                ey[i][j][k] =  -(phiy2 - phiy1)  / (2 * meshWidth[1]);  //* beamRms[2] * 2 / 2304;      // [V/m]
                ez[i][j][k] =  -(phiz2 - phiz1)  / (2 * meshWidth[2]);  // * beamRms[2] * 2 / 2304;      // [V/m]

                fout<<setw(15)<< meshx[i] 
                    <<setw(15)<< meshy[j]
                    <<setw(15)<< meshz[k]
                    <<setw(15)<< ex[i][j][k]
                    <<setw(15)<< ey[i][j][k]
                    <<setw(15)<< ez[i][j][k]
                    <<setw(15)<< phi[i][j][k]
                    <<setw(15)<< rho[i][j][k]
                    <<endl; 

            }
        }
    }

    cout<<__LINE__<<__FILE__<<endl;
    getchar();
}


void PIC3D::SetPartSCField(vector<vector<double>>  &particles)
{
    double e1=0,e2=0,e3=0,b1=0,b2=0,b3=0;
    int idx, idy, idz;

    for(int i=0;i<partNum;++i)
    {
        if(isPartOutMesh[i]==1)
        {
            partEField[0][i] = 0.E0 ;
            partEField[1][i] = 0.E0 ;
            partEField[2][i] = 0.E0 ;

        }
        else
        {
            idx  =  partMeshIndX[i];
            idy  =  partMeshIndY[i];
            idz  =  partMeshIndZ[i];
        
            partEField[0][i]    = ex[idx  ][idy  ][idz  ]     * weigh[i][0]
                                + ex[idx+1][idy  ][idz  ]     * weigh[i][1]
                                + ex[idx  ][idy+1][idz  ]     * weigh[i][2]
                                + ex[idx+1][idy+1][idz  ]     * weigh[i][3]
                                + ex[idx  ][idy  ][idz+1]     * weigh[i][4]
                                + ex[idx+1][idy  ][idz+1]     * weigh[i][5]
                                + ex[idx  ][idy+1][idz+1]     * weigh[i][6]
                                + ex[idx+1][idy+1][idz+1]     * weigh[i][7];
                    
             partEField[1][i]   = ey[idx  ][idy  ][idz  ]     * weigh[i][0]
                                + ey[idx+1][idy  ][idz  ]     * weigh[i][1]
                                + ey[idx  ][idy+1][idz  ]     * weigh[i][2]
                                + ey[idx+1][idy+1][idz  ]     * weigh[i][3]
                                + ey[idx  ][idy  ][idz+1]     * weigh[i][4]
                                + ey[idx+1][idy  ][idz+1]     * weigh[i][5]
                                + ey[idx  ][idy+1][idz+1]     * weigh[i][6]
                                + ey[idx+1][idy+1][idz+1]     * weigh[i][7];
            
             partEField[2][i]   = ez[idx  ][idy  ][idz  ]     * weigh[i][0]
                                + ez[idx+1][idy  ][idz  ]     * weigh[i][1]
                                + ez[idx]  [idy+1][idz  ]     * weigh[i][2]
                                + ez[idx+1][idy+1][idz  ]     * weigh[i][3]
                                + ez[idx  ][idy  ][idz+1]     * weigh[i][4]
                                + ez[idx+1][idy  ][idz+1]     * weigh[i][5]
                                + ez[idx  ][idy+1][idz+1]     * weigh[i][6]
                                + ez[idx+1][idy+1][idz+1]     * weigh[i][7];   // [V/m]
        }
    }
}


void PIC3D::UpdatMomentum(vector<vector<double>> &particles, const double ds)
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
        
        
        factorT = 1 / electronBeamEnergy / pow(beta,2) / gamma0; 
        factorL = 1 / electronBeamEnergy / pow(beta,2);

        // particles[3][i] += partEField[0][i] * factorT;
        // particles[4][i] += partEField[1][i] * factorT;
        // particles[5][i] += partEField[2][i] * factorL;
        // deltaPx = partEField[0][i] * ds / electronBeamEnergy / beta / p;
        // deltaPy = partEField[1][i] * ds / electronBeamEnergy / beta / p;
        
        deltaPx = partEField[0][i] * factorT * ds / gamma0;
        deltaPy = partEField[1][i] * factorT * ds / gamma0;

        particles[3][i] += deltaPx;
        particles[4][i] += deltaPy;

        // fout<<setw(15)<< particles[0][i] 
        //     <<setw(15)<< particles[1][i]
        //     <<setw(15)<< particles[2][i]
        //     <<setw(15)<< deltaPx
        //     <<setw(15)<< deltaPy
        //     <<setw(15)<< deltaPz
        //     <<endl; 

    }       

    // cout<<__LINE__<<__FILE__<<endl;
    // getchar();
    


}



