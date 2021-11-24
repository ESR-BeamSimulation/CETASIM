#ifndef BEAM_H
#define BEAM_H


#include "Global.h"
#include "Bunch.h"
#include "Train.h"
#include "LatticeInterActionPoint.h"
//#include "LongImpSingalBunch.h"
#include <fstream>
#include <vector>
#include <complex>
#include "FIRFeedBack.h"
#include "ReadInputSettings.h"
#include "WakeFunction.h"



using namespace std;
using std::vector;
using std::complex;


class Beam
{

public:
    Beam();
    ~Beam();
    

    struct WeakStrongBunchInfo{
        
        double actionJxMax=0;
        double actionJyMax=0;        
        double bunchAverXMax=0;
        double bunchAverYMax=0;
	    double bunchRmsSizeX=0;
	    double bunchRmsSizeY=0;
        double bunchEffectiveSizeXMax=0;
        double bunchEffectiveSizeYMax=0;
        double bunchZ=0;
        double bunchZP=0;
    	vector<vector<double> > bunchInfoOneTurn;   
    };
    WeakStrongBunchInfo *weakStrongBunchInfo = new WeakStrongBunchInfo;    
    
    struct StrongStrongBunchInfo{

        double emitYMax;
        double emitXMax;
        double effectivEemitXMax;
        double effectivEemitYMax;
        double bunchSizeXMax;
        double bunchSizeYMax;
        double bunchEffectiveSizeXMax;
        double bunchEffectiveSizeYMax;
    	vector<vector<double> > bunchInfoOneTurn;

            
        double bunchAverXMax;
        double bunchAverYMax;
	    double bunchRmsSizeX;
	    double bunchRmsSizeY;

    };
    StrongStrongBunchInfo *strongStrongBunchInfo = new StrongStrongBunchInfo;    

    
    int harmonics;
    
    double ionLossBoundary;
        
    vector<double> synchRadDampTime;   // in the unit of number of truns, Ref. ZhangYuan and Ohmi, PRAB 8 074402.
    vector<Bunch> beamVec;
	
	
    // function in weak-strong calculation
    void Initial(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);
    void PrintInitialBunchDistri(ReadInputSettings &inputParameter);
    void WSGetMaxBunchInfo();  
    void WSBeamRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void WSBeamIonEffectOneInteractionPoint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k);     
    void WSIonDataPrint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,int count);    
    void RMOutPutFiles(ReadInputSettings &inputParameter);
    void BeamSynRadDamping(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint);
    void FIRBunchByBunchFeedback(FIRFeedBack &firFeedBack,int nTurns);
 
    void Run(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);
    void BeamTransferPerInteractionPointDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void WSBeamDataPrintPerTurn(int turns, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);
    void BeamTransferPerTurnDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint);
    void BeamTransferPerTurnDueToLatticeL(ReadInputSettings &inputParameter);
    void BeamTransferPerTurnDueToLattice(LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);

    void BeamTransferPerTurnDueWake();
    // for long range RW wake function
	void LRWakeBeamIntaction(ReadInputSettings &inputParameter, WakeFunction &wakefunction);  


    //- not checked yet.........    
    void SSBunchDataPrint(Bunch &bunch,int count);
    void SSBunchDataPrint(Bunch &bunch,int count, int k);
    void SSIonDataPrint(LatticeInterActionPoint &latticeInterActionPoint, int nTurns);
    void SSIonDataPrint1(LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k);
    void SSBeamIonEffectOneTurn( LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter, int nTurns,int intevalofTurnsIonDataPrint, int printInterval);

//    void SingleBunchLongiImpedInterAction(LongImpSingalBunch &longImpSingalBunch, int nTurns, ofstream &fout);

    void SSGetMaxBunchInfo();


  
		
private:
};




#endif
