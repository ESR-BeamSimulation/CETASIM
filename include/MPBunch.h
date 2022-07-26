#ifndef MBBUNCH_H
#define MBBUNCH_H

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

class MPBunch : public Bunch
{


public:
    MPBunch();
    ~MPBunch();

    int macroEleNumSurivePerBunch;
    double rmsEffectiveRingEmitX,rmsEffectiveRingEmitY,rmsEffectiveRingEmitZ;
    double rmsEffectiveRx,rmsEffectiveRy;
    
    // used in the soft bunch transient beam loading simulation-- bin-by-bin process.
    vector<double > beamCurDenZProf;   // store the   rho(z) as funciton of z.   <I(z)>
    vector<vector<complex<double> > >  cavVolInfoVsLongBins;
    vector<vector<double > > cavForceInfoVsLongBins;
    vector<double> posZBins;
    vector<double> densProfVsBin;
    vector<double> hamiltonPotenWell;
    vector<double> densProfVsBinAnalytical;
    double zAverAnalytical;
    double bunchLengthAnalytical;
    vector<vector<double>> srWakePoten;

         
    void InitialMPBunch(const ReadInputSettings &inputParameter);        
    void DistriGenerator(const LatticeInterActionPoint &latticeInterActionPoint,const ReadInputSettings &inputParameter, int i);
    double GSSlover(const double if0);
    void GetMPBunchRMS(const LatticeInterActionPoint &latticeInterActionPoint, int k);    
    void SSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BunchTransferDueToLatticeL(const ReadInputSettings &inputParameter,CavityResonator &cavityResonator);
    void BunchTransferDueToSRWake(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction, const LatticeInterActionPoint &latticeInterActionPoint, int turns);
    vector<double> GetZMinMax();
private:

};




#endif
