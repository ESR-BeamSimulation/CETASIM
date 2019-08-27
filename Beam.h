#ifndef BEAM_H
#define BEAM_H


#include "Global.h"
#include "Bunch.h"
#include "Train.h"
#include "LatticeInterActionPoint.h"
#include "LongImpSingalBunch.h"
#include <fstream>
#include <vector>
#include <complex>
#include "FIRFeedBack.h"
#include "ReadInputSettings.h"



using namespace std;
using std::vector;
using std::complex;


class Beam
{

public:
    Beam();
    ~Beam();
    

    

    int printInterval;    
    double actionJxMax;
    double actionJyMax;
    double bunchSizeXMax;
    double bunchSizeYMax;
    vector<double> synchRadDampTime;   // in the unit of number of truns, Ref. ZhangYuan and Ohmi, PRAB 8 074402.

	vector<vector<double> > bunchInfoOneTurn;
    vector<Bunch> beamVec;
    void Initial(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);
    
    
    void Run(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);
    void SSBunchDataPrint(Bunch &bunch,int count);
    void WSIonDataPrint(LatticeInterActionPoint &latticeInterActionPoint,int count, int k);
    void WSBeamIonEffectOneTurn(Train &train, LatticeInterActionPoint &latticeInterActionPoint, int nTurns, ofstream &fout);
    void SSBeamIonEffectOneTurn(Train &train, LatticeInterActionPoint &latticeInterActionPoint, int nTurns, ofstream &fout);

    void SingleBunchLongiImpedInterAction(LongImpSingalBunch &longImpSingalBunch, int nTurns, ofstream &fout);
    void GetMaxBunchInfo();
    void FIRBunchByBunchFeedback(FIRFeedBack &firFeedBack,int nTurns);
	void IonBeamDataPrintPerTurn(int turns, LatticeInterActionPoint &latticeInterActionPoint, ofstream &fout);
	void BeamSynRadDamping(vector<double> &synchRadDampTime,LatticeInterActionPoint &latticeInterActionPoint);
	void WSBeamRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k);
	void BeamTransferPerTurnDueToLattice(LatticeInterActionPoint &latticeInterActionPoint);
	void WSBeamIonEffectOneInteractionPoint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k);
	void BeamTransferPerInteractionPointDueToLattice(LatticeInterActionPoint &latticeInterActionPoint, int k);
	
private:
};




#endif
