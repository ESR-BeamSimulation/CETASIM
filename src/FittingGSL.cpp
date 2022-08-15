//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
                                                          
#pragma once

#include "FittingGSL.h"
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit_nlinear.h>
#include "Global.h"

using namespace std;

FittingGSL::FittingGSL()
{
}

FittingGSL::~FittingGSL()
{
}

struct dataToFit {
size_t n;
double * t;
double * y;
double nu;
};

int sin_f(const gsl_vector * x, void *data, gsl_vector * f)
{
  size_t  n = ((struct dataToFit *)data)->n;
  double *t = ((struct dataToFit *)data)->t;
  double *y = ((struct dataToFit *)data)->y;
  double nu = ((struct dataToFit *)data)->nu;

  double A     = gsl_vector_get (x, 0);
  double theta = gsl_vector_get (x, 1);
  double b     = gsl_vector_get (x, 2);

  //double nu=0.3722;
     
  for(int i = 0; i < n; i++)
  {
    /* Model Yi = A * exp(-lambda * t_i) + b */
    double Yi = A * sin (2 * PI * nu * t[i] + theta) + b;
    gsl_vector_set (f, i, Yi - y[i]);
  }
 
    return GSL_SUCCESS;
}


int sin_df (const gsl_vector * x, void *data, gsl_matrix * J)
{
  size_t  n = ((struct dataToFit *)data)->n;
  double *t = ((struct dataToFit *)data)->t;
  double nu = ((struct dataToFit *)data)->nu;

  double A     = gsl_vector_get (x, 0);
  double theta = gsl_vector_get (x, 1);
  double b     = gsl_vector_get (x, 2);
  
  //double nu=0.3722;

  for (size_t i = 0; i < n; i++)
  {
    /* Jacobian matrix J(i,j) = dfi / dxj, */
    /* where fi = (Yi - yi)/sigma[i],      */
    /*       Yi = A * exp(-lambda * t_i) + b  */
    /* and the xj are the parameters (A,lambda,b) */
    
    double d1 =     sin(2 * PI * nu * t[i] + theta);
    double d2 = A * cos(2 * PI * nu * t[i] + theta);
    double d3 = 1;
    gsl_matrix_set (J, i, 0, d1);
    gsl_matrix_set (J, i, 1, d2);
    gsl_matrix_set (J, i, 2, d3);
  }

  return GSL_SUCCESS;
}

void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double rcond;

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);

  //fprintf(stderr, "iter %2zu: A = %.8f, theta = %.8f, b = %.8f, cond(J) = %8.4f, |f(x)| = %.8f\n",
  //        iter,
  //        gsl_vector_get(x, 0),
  //        gsl_vector_get(x, 1),
  //        gsl_vector_get(x, 2),
  //        1.0 / rcond,
  //        gsl_blas_dnrm2(f));
}


//int (FittingGSL::* cos_f) (const gsl_vector * x, void *data, gsl_vector * f)
//{
//    return 0;
//}


vector<double> FittingGSL::FitASin(const vector<double> &xHistoryData,const double nux)
{
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
 
  const size_t n = xHistoryData.size();
  const size_t p = 3;

  gsl_vector *f;
  gsl_matrix *J;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  double t[n], y[n], weights[n];
  struct dataToFit d = { n, t, y, nux};
  double x_init[3] = { 1.0, 1.0, 0.0 }; /* starting values */
  gsl_vector_view x = gsl_vector_view_array (x_init, p);
  gsl_vector_view wts = gsl_vector_view_array(weights, n);
  gsl_rng * r;
  double chisq, chisq0;
  int status, info;
  size_t i;

  const double xtol = 1e-10;
  const double gtol = 1e-10;
  const double ftol = 0.0;

  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);

  /* define the function to be minimized */
  fdf.f = sin_f;
  fdf.df = sin_df;   /* set to NULL for finite-difference Jacobian */
  fdf.fvv = NULL;     /* not using geodesic acceleration */
  fdf.n = n;
  fdf.p = p;
  fdf.params = &d;

  /* this is the data to be fitted */
  for (i = 0; i < n; i++)
    {
      /*
      double ti = i * TMAX / (n - 1.0);
      double yi = 1.0 + 5 * exp (-1.5 * ti);
      double si = 0.1 * yi;
      double dy = gsl_ran_gaussian(r, si);
      */
      t[i] = i+1;
      y[i] = xHistoryData[i];
      double si   = 0.1 * y[i];
      weights[i] = 1.0 / (si * si);
    };

  /* allocate workspace with default parameters */
  w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

  /* compute initial cost function */
  f = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(f, f, &chisq0);

  /* solve the system with a maximum of 100 iterations */
  status = gsl_multifit_nlinear_driver(100, xtol, gtol, ftol,callback, NULL, &info, w);

  /* compute covariance of best fit parameters */
  J = gsl_multifit_nlinear_jac(w);
  gsl_multifit_nlinear_covar (J, 0.0, covar);

  /* compute final cost */
  gsl_blas_ddot(f, f, &chisq);


  /*
  fprintf(stderr, "summary from method '%s/%s'\n",  gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
  fprintf(stderr, "number of iterations: %zu\n",    gsl_multifit_nlinear_niter(w));
  fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
  fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
  fprintf(stderr, "reason for stopping: %s\n",   (info == 1) ? "small step size" : "small gradient");
  fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
  fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));
 
  double dof = n - p;
  double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

  fprintf(stderr, "chisq/dof = %g\n", chisq / dof);

  fprintf (stderr, "A      = %.8f +/- %.5f\n", gsl_vector_get(w->x, 0), c*sqrt(gsl_matrix_get(covar,0,0)));
  fprintf (stderr, "theta  = %.8f +/- %.5f\n", gsl_vector_get(w->x, 1), c*sqrt(gsl_matrix_get(covar,1,1)));
  fprintf (stderr, "b      = %.8f +/- %.5f\n", gsl_vector_get(w->x, 2), c*sqrt(gsl_matrix_get(covar,2,2)))
  

  fprintf (stderr, "status = %s\n", gsl_strerror (status));
 */

  std::vector<double> fitRes;
  fitRes.resize(2);

  fitRes[0] = gsl_vector_get(w->x, 0);
  fitRes[1] = gsl_vector_get(w->x, 1);
  
  gsl_multifit_nlinear_free (w);
  gsl_matrix_free (covar);
  gsl_rng_free (r);

    
  return fitRes;
    
}


vector<double> FittingGSL:: FitALinear(double *x, double *w, double *y, int n)
{
  double c0, c1, cov00, cov01, cov11, chisq;
  
  gsl_fit_wlinear (x, 1, w, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
  
  vector<double> fitRes(2,0.E0);
  fitRes[0] = c0;
  fitRes[1] = c1;

  return fitRes;

}
