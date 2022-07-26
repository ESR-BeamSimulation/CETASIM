#ifndef FittingGSL_H
#define FittingGSL_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#include <vector>
using std::vector;




class FittingGSL
{
public:
    FittingGSL();
    ~FittingGSL();

 

    vector<double> FitASin(const vector<double> &xHistoryData, const double nux);
    
private:

};




#endif
