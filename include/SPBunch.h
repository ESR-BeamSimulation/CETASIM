//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#ifndef SBBUNCH_H
#define SBBUNCH_H

#include <vector>
#include <complex>
#include "Global.h"
#include "LatticeInterActionPoint.h"
#include "Train.h"
#include "ReadInputSettings.h"
#include "CavityResonator.h"
#include "Bunch.h"
#include "WakeFunction.h"
#include "Spline.h"
using namespace std;
using std::vector;
using std::complex;

class SPBunch : public Bunch
{


public:
    SPBunch();
    ~SPBunch();
  
    double actionJx=0.E0;        // used for weak-strong simulation
    double actionJy=0.E0;        // used for weak-strong simulation
      
    void InitialSPBunch(const ReadInputSettings &inputParameter);        
    void DistriGenerator(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter, int i);
    void GetSPBunchRMS(const LatticeInterActionPoint &latticeInterActionPoint, int k);    
    void WSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k);
    // beam loading simulation, the same approaches as elegant.
    // have to benchmark with cases when there exist instabilites.....(single and double cavities.....)
    void BunchMomentumUpdateDueToRF(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator);
    
    // beam loading simulation, yamamoto's paper, it is not 100% correct approaches. 
    // For certain partilce, Vg exp(li * phi_p) + Vb = Vc, phi_p is particle phase refer to generator. 
    void BunchMomentumUpdateDueToRFYamamoto(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator);    
    // beam loading simulation, used to benckmark the coupled bunch instablity growth rate as idea cases
    void BunchMomentumUpdateDueToRFMatrix(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator);
    
private:

};




#endif
