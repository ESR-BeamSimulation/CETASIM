//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************

#pragma once

#include "MatrixDef2D.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_matrix.h>

using namespace std;
using std::vector;
using std::complex;

MatrixDef2D::MatrixDef2D()
{
}

MatrixDef2D::~MatrixDef2D()
{
	Mat2DFree();
}


void MatrixDef2D::Mat2DCreate(int rows, int cols)
{
	mat2D = gsl_matrix_alloc (rows, cols);
	gsl_matrix_set_zero(mat2D);
}

void MatrixDef2D::Mat2DFree()
{
	gsl_matrix_free(mat2D);
}







