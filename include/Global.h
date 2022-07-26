#ifndef GLOBAL_H
#define GLOBAL_H

#include <cmath>
#include <complex>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <vector>
#include <string>
#include <algorithm>

using std::complex;
using std::vector;
using namespace std;


const double Epsilon 				= 8.854187817E-12;
const double CLight 				= 299792458;                  // m/s
const double ElectronCharge         = 1.602017733e-19;            // C 
const double ElectronMass           = 9.10938215E-31;             // KG 
const double ElectronMassEV         = 0.51099906e6;               // eV
const double IonMass                = 1.672621777e-27;            // KG
const double IonMassEV              = 931.49432e+6;               // eV
const double PI                     = 4*atan(1.E0);
const double ElecClassicRadius      = 2.81794092E-15;             // m
const double ProtonClassicRadius    = 1.535E-18;                  // m
const double Boltzmann              = 1.3806505E-23;              // J/K
const double VaccumZ0				= 120 * PI;
const double FactorGaussSI        = VaccumZ0 * CLight / (4 * PI); //[ohm] [m/s]
const  complex<double> li(0,1);

extern int      numProcess;
extern int      myRank;


double Gaussrand(double rms, double aver,double randomIndex);
void gsl_matrix_mul(gsl_matrix *a,gsl_matrix *b,gsl_matrix *c);
void gsl_matrix_inv(gsl_matrix *a);
double get_det(gsl_matrix * A);
int StringVecSplit(string str, vector<string> &strVec);
int StringSplit(string str, vector<string> &strVec);
int StringSplit2(string str, vector<string> &strVec);
void RMOutPutFiles();
double vectorDisNorm2(vector<double> a,vector<double> b);


#endif
