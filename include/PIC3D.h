//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************

#ifndef PIC3D_H
#define PIC3D_H


#include <vector>
#include <complex>
#include <fftw3.h>
#include "PIC2D.h"
#include "ReadInputSettings.h"
using namespace std;
using std::complex;

using std::vector;
using v1i= vector<int> ;
using v1d= vector<double> ;
using v2d= vector<vector<double> > ;
using v3d= vector<vector<vector<double> > > ;


class PIC3D 
{

public: 
    PIC3D();
    ~PIC3D();

    // for self space-charge
    void InitialSC3D(const ReadInputSettings &inputParameter);

    int numberOfGrid[3] = {128,128,64};               // nx; ny; nz;
    double meshWidth[3] = {0.E0,0E0,0E0};           //
    double meshCen[3]   = {0.E0,0E0,0E0};
    vector<double> beamRms;
    double electronBeamEnergy;
    double dv = 1.0;
    int partNum;
    v1d meshx;    // meshx.size() = nx
    v1d meshy;    
    v1d meshz;
    v2d partEField;


    v1i partMeshIndX;  // countx[partNum] 
    v1i partMeshIndY;  // county[partNum]
    v1i partMeshIndZ;  // countz[partNum]
    v1i isPartOutMesh;
    v2d weigh;

    v3d rho;    // nx*ny*nz
    v3d phi;
    v3d ex;
    v3d ey;
    v3d ez;
    double range = 10.E0;

    // 3D FFT
    double              *in_xy,   *out_xy;
    double              *in,      *out,      *outf;  
    fftw_complex        *in_z,    *out_z,    *outf_z;
    fftw_plan           rho_X_Y_Z_To_Rho_KX_KY_Z,  phi_KX_KY_Z_To_Phi_X_Y_Z;
    fftw_r2r_kind       kind_xy_forward[2]  = { FFTW_RODFT00,FFTW_RODFT00 };
    fftw_r2r_kind       kind_xy_backward[2] = { FFTW_RODFT00,FFTW_RODFT00 };

    fftw_plan           rho_KX_KY_Z_To_Rho_KX_KY_KZ,     phi_KX_KY_KZ_To_Phi_KX_KY_Z;
  
    void Set3DMesh      (vector<vector<double>>  &particles);        // particles[3 , np], 2d array particle distribution in real space  
    void Set3DRho       (vector<vector<double>>  &particles, double macroCharge);     
    void Set3DPhi       ();        //
    void Set3DEField    ();
    void SetPartSCField (vector<vector<double>>  &particles);
    void UpdatMomentum  (vector<vector<double>>  &particles, const double ds);

    void Set3DPhi1();
    double           *in_xyz,   *out_xyz;
    fftw_plan        rho_X_Y_Z_To_Rho_KX_KY_KZ,  phi_KX_KY_KZ_To_Phi_X_Y_Z;
    fftw_r2r_kind    kind_xyz_forward[3]  = { FFTW_RODFT00,FFTW_RODFT00, FFTW_RODFT00 };
    fftw_r2r_kind    kind_xyz_backward[3] = { FFTW_RODFT00,FFTW_RODFT00, FFTW_RODFT00 };

    // applied for 2.5D model, sliced in longitudinnal, each longitudianl slice track partilce with a 2D PIC   
    PIC2D *slicedBunchPIC2D; 
    // set the 1d space charge solver for SC Ez


};




#endif
