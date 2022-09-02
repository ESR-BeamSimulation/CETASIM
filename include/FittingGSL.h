//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#ifndef FittingGSL_H
#define FittingGSL_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
// #include <gsl/gsl_multifit_nlinear.h>

#include <vector>
using std::vector;




class FittingGSL
{
public:
    FittingGSL();
    ~FittingGSL();

    vector<double> FitASin(const vector<double> &xHistoryData, const double nux);   
    
    vector<double> FitALinear(double *x, double *w, double *y, int n);

private:

};




#endif
