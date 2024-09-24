//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************

#ifndef MATRIXDEF2D_H
#define MATRIXDEF2D_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_matrix.h>
#include <vector>
#include <complex>
#include <iostream>
#include <iomanip>


using namespace std;
using std::vector;
using std::complex;

class MatrixDef2D
{

public:

    MatrixDef2D();
	gsl_matrix *mat2D;  
    void Mat2DCreate(int rows, int cols);
    void Mat2DFree();
    
    ~MatrixDef2D();

private:


};


#endif
