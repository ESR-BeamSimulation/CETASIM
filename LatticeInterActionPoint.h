#ifndef LATTICEINTERACTIONPOINT_H
#define LATTICEINTERACTIONPOINT_H


#include "Global.h"
#include <vector>
#include "ReadInputSettings.h"

using namespace std;
using std::vector;
using std::complex;


class LatticeInterActionPoint
{


public:
    LatticeInterActionPoint();
    ~LatticeInterActionPoint();

    int numberOfInteraction;
    double circRing;
    
    
    double pipeAperatureR;
    double pipeAperatureX;
    double pipeAperatureY;
    int harmonics;

    vector<double> twissAlphaX;
    vector<double> twissBetaX;              // m/rad   betatron function at center the high betatron straight section
    vector<double> twissAlphaY;
    vector<double> twissBetaY;              // m/rad
    vector<double> twissAlphaZ;
    vector<double> twissBetaZ;              // m/rad
    
    vector<double> twissDispX;				// m
    vector<double> twissDispPX;             // rad
    vector<double> twissDispY;				// m
    vector<double> twissDispPY;             // rad    
   
    

    vector<vector<double> > ionPositionX;            // m
    vector<vector<double> > ionPositionY;            // m
    vector<vector<double> > ionVelocityX;            // m/s
    vector<vector<double> > ionVelocityY;            // m/s
    

    
    vector<int> ionAccumuNumber;
    vector<vector<double> > ionAccumuPositionX;            // m -- max number accumulated is set to 1000;
    vector<vector<double> > ionAccumuPositionY;            // m
    vector<vector<double> > ionAccumuVelocityX;            // m/s
    vector<vector<double> > ionAccumuVelocityY;            // m/s
    
    vector<double> ionAccumuAverX;
    vector<double> ionAccumuAverY;
    vector<double> ionAccumuRMSX;
    vector<double> ionAccumuRMSY; 
    
    vector<vector<double> > ionAccumuFx;            // 
    vector<vector<double> > ionAccumuFy;            // 
    
    vector<double> vacuumPressure;          // Pa
    vector<double> temperature;             // K
    vector<double> ionLineDensity;          // [1/m] ion number per meter in rings. 
    vector<double> interactionLength;       // [m]   effective length of interaction point
                                            //       unifrom disributed with numberOfInteraction slices

    vector <int>   macroIonNumber;

    vector<double> ionNumber;
    int ionMaxNumberOneInterPoint;
    vector<double> macroIonCharge;
    
    vector<vector<double> > transferMatrix;
    vector<vector<double> > xTransferMatrix;
    vector<vector<double> > yTransferMatrix;
    vector<vector<double> > zTransferMatrix;
    vector<double> xPhaseAdv;
    vector<double> yPhaseAdv;
    vector<double> zPhaseAdv;
    

    void Initial(ReadInputSettings &inputParameter);
    void InitialLattice(ReadInputSettings &inputParameter);
    void IonGenerator(double rmsRx, double rmsRy, double xAver,double yAver, int k);
    void IonsUpdate(int k);
    void IonTransferDueToBunch(int bunchGap, int k);
    void IonRMSCal(int k);

private:

};




#endif
