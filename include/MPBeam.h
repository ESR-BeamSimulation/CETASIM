#ifndef MPBEAM_H
#define MPBEAM_H

#include "Global.h"
#include "MPBunch.h"
#include "Train.h"
#include "LatticeInterActionPoint.h"
#include "ReadInputSettings.h"
#include "WakeFunction.h"
#include "FIRFeedBack.h"
#include <fstream>
#include <vector>
#include <complex>
#include "CavityResonator.h" 


using namespace std;
using std::vector;
using std::complex;


class MPBeam 
{
public:
    MPBeam();
    ~MPBeam();
   
    struct StrongStrongBunchInfo{

        double emitYMax;
        double emitXMax;
        double emitZMax;
        double effectivEemitXMax;
        double effectivEemitYMax;
        double bunchSizeXMax;
        double bunchSizeYMax;
        double bunchLengthMax;
        double bunchEnergySpreadMax;
        double bunchEffectiveSizeXMax;
        double bunchEffectiveSizeYMax;
    	vector<vector<double> > bunchInfoOneTurn;

        double bunchAverXMax;
        double bunchAverYMax;
        double bunchAverZMax;
        double bunchRmsSizeX=0;
	    double bunchRmsSizeY=0;
	    double bunchRmsSizeZ=0;
	    double bunchRmsSizePX=0;
	    double bunchRmsSizePY=0;
	    double bunchRmsSizePZ=0;
	    double bunchAverX=0;
	    double bunchAverY=0;
	    double bunchAverZ=0;
	    double bunchAverPX=0;
	    double bunchAverPY=0;
	    double bunchAverPZ=0;
        
    };
    StrongStrongBunchInfo *strongStrongBunchInfo = new StrongStrongBunchInfo;      
    vector<MPBunch> beamVec;
    vector<double> timeBetweenBunch; 

    void Run(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter, CavityResonator &cavityResonator);
    void MPBeamRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void MPGetBeamInfo();
    void MPBeamDataPrintPerTurn(int turns, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);
    void SSBeamIonEffectOneInteractionPoint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint, int nTurns, int k);
    void GetBinDistBetweenBunch(ReadInputSettings &inputParameter);
    void SRWakeBeamIntaction(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction, const  LatticeInterActionPoint &latticeInterActionPoint,int turns);    
    


    //// shared funcitons by MP and SP--just copy of each other.         
    void Initial(Train &train, LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter);
    void InitialcavityResonator(ReadInputSettings &inputParameter,CavityResonator &cavityResonator);    
    void BeamTransferPerInteractionPointDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint, int k);
    void BeamTransferPerTurnDueToLatticeL(ReadInputSettings &inputParameter,CavityResonator &cavityResonator);
    void BeamTransferPerTurnDueToLattice(LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter,CavityResonator &cavityResonator);
    void BeamTransferPerTurnDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint);
    void SSIonDataPrint(ReadInputSettings &inputParameter,LatticeInterActionPoint &latticeInterActionPoint,int count);  
    //void SetBeamPosHistoryDataWithinWindow();
                  
                               	                       
    void BeamSynRadDamping(const ReadInputSettings &inputParameter,const LatticeInterActionPoint &latticeInterActionPoint);
    void FIRBunchByBunchFeedback(FIRFeedBack &firFeedBack,int nTurns);
    void BeamTransferPerTurnDueWake();
    //// for long range RW wake function
	void LRWakeBeamIntaction(const  ReadInputSettings &inputParameter, WakeFunction &wakefunction, const  LatticeInterActionPoint &latticeInterActionPoint);  
    

      		
private:
};




#endif
