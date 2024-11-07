//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************

#ifndef PIC2D_H
#define PIC2D_H

#include <vector>
#include <complex>
#include <fftw3.h>
#include "ReadInputSettings.h"
using namespace std;
using std::complex;


using std::vector;
using v1i= vector<int> ;
using v1d= vector<double> ;
using v2d= vector<vector<double> > ;
using v3d= vector<vector<vector<double> > > ;


class PIC2D
{


public:
    PIC2D();
    ~PIC2D();

    v1i numberOfGrid;          // nx; ny;
    double meshWidth[2] = {0.E0,0E0};           //
    double meshCen[2]   = {0.E0,0E0};
    void InitialPIC2D(vector<int> meshNum);

    double beamrmssize[2];

    double electronBeamEnergy;
    double dv = 1.0;
    int partNum;
    v1d meshx;    // meshx.size() = nx
    v1d meshy;    
    v2d partExEyField;


    v1i partMeshIndX;  // countx[partNum] 
    v1i partMeshIndY;  // county[partNum]
    v1i isPartOutMesh;
    v2d weigh;

    v2d rho;    // nx*ny
    v2d phi;
    v2d ex;
    v2d ey;
    double range = 10.E0;


    double              *in_xy,   *out_xy;
    fftw_plan           rho_X_Y_To_Rho_Kx_Ky, phi_Kx_Ky_To_Phi_X_Y;
    fftw_r2r_kind       kind_xy_forward[2]  = { FFTW_RODFT00,FFTW_RODFT00 };
    fftw_r2r_kind       kind_xy_backward[2] = { FFTW_RODFT00,FFTW_RODFT00 };

    void Set2DMesh(vector<vector<double>>  &particles);
    void Set2DMesh(double rmsXY[2], double averXY[2]);  // mesh is generated accoring to the rms beam size from the whole bunch
    void Set2DRho (vector<vector<double>>  &particles, vector<double> charge);    
    void Set2DPhi();
    void Set2DEField();
    void SetPartSCField ();
    void UpdateTransverseMomentum  (vector<vector<double>>  &particles, const double ds);

    vector<vector<double> > GetPartSCField(vector<vector<double>>  &particles);  // get the sc field accoriding the input particle distribution

};




#endif
