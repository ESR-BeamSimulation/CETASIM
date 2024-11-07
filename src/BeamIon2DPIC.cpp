//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             

#include "Bunch.h"
#include "Global.h"
#include "BeamIon2DPIC.h"
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

BeamIon2DPIC::BeamIon2DPIC()
{
}

BeamIon2DPIC::~BeamIon2DPIC()
{

}

void BeamIon2DPIC::InitialPIC2D()
{
    pic2DBeam.InitialPIC2D(meshNumBeam);
    pic2DIon.InitialPIC2D(meshNumIon);
}
