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
    vector<double> ionMassNumber;
    vector<double> corssSectionEI;
    vector<double> gasPercent;
    int gasSpec;
     
    double pipeAperatureX;
    double pipeAperatureY;
    double pipeAperatureR;

    int harmonics;

    vector<double> twissAlphaX;
    vector<double> twissBetaX;              // m/rad  
    vector<double> twissAlphaY;
    vector<double> twissBetaY;              // m/rad
    vector<double> twissAlphaZ;
    vector<double> twissBetaZ;              // m/rad
    
    vector<double> twissDispX;				// m
    vector<double> twissDispPX;             // rad
    vector<double> twissDispY;				// m
    vector<double> twissDispPY;             // rad    
   
    

    vector<vector<vector<double> > >ionPositionX;            // m   //at the kth interaction point, pth ion species..
    vector<vector<vector<double> > >ionPositionY;            // m
    vector<vector<vector<double> > >ionVelocityX;            // m/s
    vector<vector<vector<double> > >ionVelocityY;            // m/s
    

    vector<vector<int> >ionAccumuNumber;
    vector<vector<vector<double> > >ionAccumuPositionX;            // m -- max number accumulated is set to 1000;
    vector<vector<vector<double> > >ionAccumuPositionY;            // m
    vector<vector<vector<double> > >ionAccumuVelocityX;            // m/s
    vector<vector<vector<double> > >ionAccumuVelocityY;            // m/s
    
    vector<vector<double> >ionAccumuAverX;
    vector<vector<double> >ionAccumuAverY;
    vector<vector<double> >ionAccumuAverVelX;
    vector<vector<double> >ionAccumuAverVelY;

    vector<vector<double> >ionAccumuRMSX;
    vector<vector<double> >ionAccumuRMSY; 
    
    vector<vector<vector<double> > >ionAccumuFx;            // 
    vector<vector<vector<double> > >ionAccumuFy;            // 
    
    vector<double>vacuumPressure;          // Pa
    vector<vector<double> >vacuumPressureEachGas;          // Pa
    vector<double> temperature;             // K
    vector<vector<double> >ionLineDensity;          // [1/m] ion number per meter in rings. 
    vector<double> interactionLength;       // [m]   effective length of interaction point
                                            //       unifrom disributed with numberOfInteraction slices
    vector<vector<double>  >ionNumber;
    vector<vector<int   >  >macroIonNumber;
    vector<vector<double>  >macroIonCharge;
    

    int ionMaxNumberOneInterPoint;
    
    vector<vector<double> > transferMatrix;
    vector<vector<double> > xTransferMatrix;
    vector<vector<double> > yTransferMatrix;
    vector<vector<double> > zTransferMatrix;
    vector<double> xPhaseAdv;
    vector<double> yPhaseAdv;
    vector<double> zPhaseAdv;
    
    double ionLossBoundary;
    vector<double> totIonChargeAtInterPoint;
    vector<int>  totMacroIonsAtInterPoint;
    int    totMacroIons;
    double totIonCharge;
    
    void Initial(ReadInputSettings &inputParameter);
    void InitialLattice(ReadInputSettings &inputParameter);
    void IonGenerator(double rmsRx, double rmsRy, double xAver,double yAver, int k);
    void WSIonsUpdate(int k);
    void IonRMSCal(int k, int p);
    
    void SSIonsUpdate(double bunchEffectiveSizeXMax, double bunchEffectiveSizeYMax, double xAver,double yAver, int k);

    void IonTransferDueToBunch(int bunchGap, int k, double bunchSizeXMax, double bunchSizeYMax);

    void GetTotIonCharge();
    void GetIonNumberPerInterAction(double electronNumPerBunch, int k);
    
private:

};




#endif
