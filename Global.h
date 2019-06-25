#ifndef GLOBAL_H
#define GLOBAL_H

#include <cmath>
#include <complex>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

using std::complex;



const double Epsilon 				= 8.854187817E-12;
const double CLight 				= 299792458;                  // m/s
const double ElectronCharge         = 1.602017733e-19;             // C 
const double ElectronMass           = 9.10938215E-31;              // KG 
const double ElectronMassEV         = 0.51099906e6;                // eV
const double IonMass                = 1.672621777e-27;             // KG
const double IonMassEV              = 931.49432e+6;                  // eV
const double PI                     = 4*atan(1.E0);
const double ElecClassicRadius      = 2.81794092E-15;             // m
const double Boltzmann              = 1.3806505E-23;             // J/K



const double CorssSectionEI         = 2.0e-22;                    //  [m^2] Collision crossing section of CO
const double PipeAperatureX         = 1.1E-3;                      // [m] 
const double PipeAperatureY         = 1.1E-3;                      // [m] 
const double PipeAperatureR         = 1.1E-3;                      // [m] 
const int    IonMaxNumber           = 1E+5;                 
const  complex<double> li(0,1);

extern int      numProcess;
extern int      myRank;
extern double   Omegas;                // [Hz] sychronous revolution frequency 
extern double   T0;                    // [s]  EvolutionPeriod  T0
extern int      Harmonics;             // []   number of harmonics 
double Gaussrand(double rms, double aver,int randomIndex);
void gsl_matrix_mul(gsl_matrix *a,gsl_matrix *b,gsl_matrix *c);
void gsl_matrix_inv(gsl_matrix *a);
double get_det(gsl_matrix * A);


// HEPS lattice parameters
const double CircRing               = 1360.4E0;                     // m
const double WorkQx                 = 114.05E0; //114.141E0; //114.19;
const double WorkQy                 = 106.231E0; //.1571;
const double WorkQz                 = 0.00521; //.1571;
const double RFBaseFrequency        = 166.6E6;                    // Hz
//************************************************************************

// SIAP lattice parameters
//const double CircRing               = 432;                     // m
//const double WorkQx                 = 22.22;
//const double WorkQy                 = 11.29;
//const double RFBaseFrequency        = 499.645E6;                    // Hz
//************************************************************************






#endif
