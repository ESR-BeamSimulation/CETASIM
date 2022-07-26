#include "Bunch.h"
#include "Global.h"
#include "Faddeeva.h"
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <complex>
#include <iostream>
#include<iomanip>
#include <fstream>
#include <numeric>
#include <cmath>
#include <random>
#include <gsl/gsl_fft_complex.h>
#include <algorithm>
#include <gsl/gsl_matrix.h>



Bunch::Bunch()
{
}

Bunch::~Bunch()
{
}

void Bunch::Initial(LatticeInterActionPoint &latticeInterActionPoint, ReadInputSettings &inputParameter)
{


    macroEleNumPerBunch = inputParameter.ringBunchPara->macroEleNumPerBunch;      // macroEleNumPerBunch=1 -- electron beam weak-strong model
                                                                                  // macroEleNumPerBunch!=1 -- electron beam strong-strong model

    kappa =  inputParameter.ringBunchPara->kappa;
    pipeAperatureR =   inputParameter.ringParBasic->pipeAperature[0];
    pipeAperatureX =   inputParameter.ringParBasic->pipeAperature[1];
    pipeAperatureY =   inputParameter.ringParBasic->pipeAperature[2];                            



    ePositionX.resize(macroEleNumPerBunch);
    ePositionX.resize(macroEleNumPerBunch);
    ePositionY.resize(macroEleNumPerBunch);
    ePositionZ.resize(macroEleNumPerBunch);

    eMomentumX.resize(macroEleNumPerBunch);
    eMomentumY.resize(macroEleNumPerBunch);
    eMomentumZ.resize(macroEleNumPerBunch);

    eFx.resize(macroEleNumPerBunch);
    eFy.resize(macroEleNumPerBunch);
    eSurive.resize(macroEleNumPerBunch);


    
    bunchBinNumberZ = inputParameter.ringImpedance->bunchBinNumberZ;
    beamCurDenZProf.resize(bunchBinNumberZ);
    beamCurDenTProf.resize(bunchBinNumberT);
    
    current = inputParameter.ringBunchPara->current;

    electronEnergy  = inputParameter.ringParBasic->electronBeamEnergy;
    rmsEnergySpread = inputParameter.ringBunchPara-> rmsEnergySpread;
    rmsBunchLength  = inputParameter.ringBunchPara-> rmsBunchLength;

    emittanceZ      = inputParameter.ringBunchPara->emittanceZ;
    emittanceX      = inputParameter.ringBunchPara->emittanceX;
    emittanceY      = inputParameter.ringBunchPara->emittanceY;
    

    rGamma = inputParameter.ringParBasic->rGamma;
    rBeta  = inputParameter.ringParBasic->rBeta;


    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        ePositionX[i]=0.E0;
        ePositionY[i]=0.E0;
        ePositionZ[i]=0.E0;
        eMomentumX[i]=0.E0;
        eMomentumY[i]=0.E0;
        eMomentumZ[i]=0.E0;
        eFx[i]       =0.E0;
        eFy[i]       =0.E0;
        eSurive[i]   =1;
    }


    electronNumPerBunch = current / inputParameter.ringParBasic->f0  / ElectronCharge;
            
    macroEleCharge = electronNumPerBunch / macroEleNumPerBunch;
    distributionType = inputParameter.ringBunchPara->distributionType;    
}



void Bunch::DistriGenerator(LatticeInterActionPoint &latticeInterActionPoint,ReadInputSettings &inputParameter,int randomIndex)
{

// the first particle is set as referecen particle -- it is with (x,x',y,yp,z,z') = 0 

// longitudial bunch pahse space generatetion --simple rms in both z and z' phase space.

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> doffset{0,1}; 
    

    double initialStaticOffSet[6];
    double initialDynamicOffSet[6];


    for(int i=0;i<6;i++)
    {
       initialStaticOffSet[i]   = inputParameter.ringBunchPara->initialStaticOffSet[i];
       initialDynamicOffSet[i]  = inputParameter.ringBunchPara->initialDynamicOffSet[i];
    }    


    srand(time(0)+randomIndex);
    double disDx = doffset(gen) * initialDynamicOffSet[0] + initialStaticOffSet[0];
    double disDy = doffset(gen) * initialDynamicOffSet[1] + initialStaticOffSet[1];
    double disDz = doffset(gen) * initialDynamicOffSet[2] + initialStaticOffSet[2];
    double disMx = doffset(gen) * initialDynamicOffSet[3] + initialStaticOffSet[3];
    double disMy = doffset(gen) * initialDynamicOffSet[4] + initialStaticOffSet[4];
    double disMz = doffset(gen) * initialDynamicOffSet[5] + initialStaticOffSet[5];

    
   
 
    std::normal_distribution<> dx{0,rmsBunchLength};
    std::normal_distribution<> dy{0,rmsEnergySpread};

    ePositionX[0]=0.E0;
    ePositionY[0]=0.E0;
    ePositionZ[0]=0.E0;
    eMomentumX[0]=0.E0;
    eMomentumY[0]=0.E0;
    eMomentumZ[0]=0.E0;


    double tempx;
    double tempy;
    double temp;
	int i=1;

    while(i<macroEleNumPerBunch)
    {
        tempx = dx(gen);
        tempy = dy(gen);
        
        temp =  pow(tempx/rmsBunchLength,2) + pow(tempy/rmsEnergySpread,2);
        
        if( temp<pow(3,2))                      // longitudinal is truncted at 3 sigma both in z and dp direction 
        {
            ePositionZ[i]=tempx;			                    // m
            eMomentumZ[i]=tempy;			                    // rad			
        }
        else
        {
            continue;
        }
        i++;
    }


    double gammaX;
    double gammaY;
    double sigmaX;
    double sigmaY;
    double dSigmaX;
    double dSigmaY;
    double alphaX;
    double alphaY;
    double betaX;
    double betaY;
	double dispersionX;
	double dispersionY;
	double dispersionPX;
	double dispersionPY;


    // bunch distribution is generated accroding the twiss parameter at the first interaction points.
    
    alphaX = latticeInterActionPoint.twissAlphaX[0];   
    betaX  = latticeInterActionPoint.twissBetaX[0];
    
    alphaY = latticeInterActionPoint.twissAlphaY[0];
    betaY  = latticeInterActionPoint.twissBetaY[0];
	dispersionX =  latticeInterActionPoint.twissDispX[0];
	dispersionY =  latticeInterActionPoint.twissDispY[0];
	dispersionPX =  latticeInterActionPoint.twissDispPX[0];
	dispersionPY =  latticeInterActionPoint.twissDispPY[0];



    gammaX = (1+pow(alphaX,2))/betaX;
    gammaY = (1+pow(alphaY,2))/betaY;

    sigmaX  =  sqrt(betaX);
    sigmaY  =  sqrt(betaY);
    dSigmaX = -alphaX/sqrt(betaX);
    dSigmaY = -alphaY/sqrt(betaY);


    double rX;
    double rY;

    rX   = sqrt(emittanceX*betaX);
    rY   = sqrt(emittanceY*betaY);

 
    double f0;
    double if0;
    double fi;
    double axax;
    double ayay;
    double ax;
    double ay;
    double phaseX;
    double phaseY;
    double ix;


    i=1;   
         
    while(i<macroEleNumPerBunch)
    {    
      switch(distributionType)
        {
            case 1:
                f0  = 4*emittanceX; 
                if0 = double(std::rand())/RAND_MAX;
                fi = f0;
                break;
            case 2:
                f0  = 6*emittanceX;
                if0 = double(std::rand())/RAND_MAX;
                fi = f0 * sqrt(if0);
                break;
            case 3:
                f0  = 3*emittanceX;
                if0 = double(std::rand())/RAND_MAX;
                fi  = GSSlover(if0);
                fi  = fi * f0 ;
                break;
            default:
            cerr<<"Do Nothing, no distribution type in generation.  "<<endl;
        }
        

        ix   = double(rand())/RAND_MAX;
        
        axax =  fi* ix;
        ayay =  (fi - axax) * kappa;  
        
        ax  = sqrt(axax);
        ay  = sqrt(ayay);
 

        phaseX = 2 * PI* double(rand())/RAND_MAX;
        phaseY = 2 * PI* double(rand())/RAND_MAX;


        ePositionX[i] = ax *   sigmaX * cos( phaseX ) ;
        ePositionY[i] = ay *   sigmaY * cos( phaseY ) ;
        eMomentumX[i] = ax * (dSigmaX * cos( phaseX ) - sin(phaseX)/sigmaX );
        eMomentumY[i] = ay * (dSigmaY * cos( phaseY ) - sin(phaseY)/sigmaY );


        // to suppress the numerical noise on average calculation
        //ePositionX[i+1] = -ePositionX[i];
        //ePositionY[i+1] = -ePositionY[i];
        //eMomentumX[i+1] = -eMomentumX[i];
        //eMomentumY[i+1] = -eMomentumY[i];
        
    

        if(distributionType==3)
        {
            if(pow(ePositionX[i]/rX,2) + pow(ePositionY[i]/rY,2)>9)
            {
                continue;
            }
        }
        
        i++;
    }
    
 
	// take the dispersion and initial error into accout   
 	   
	for(int i=0;i<macroEleNumPerBunch;i++)
	{
        ePositionZ[i] +=  disDz; 
        eMomentumZ[i] +=  disMz;    
        
                                                                
		ePositionX[i] +=  disDx + dispersionX  * eMomentumZ[i];
		ePositionY[i] +=  disDy + dispersionY  * eMomentumZ[i];
				
		eMomentumX[i] +=  disMx + dispersionPX * eMomentumZ[i];
		eMomentumY[i] +=  disMy + dispersionPY * eMomentumZ[i];		   
	}

    //eMomentumZ[0] +=  0.0001; 

	
}

double Bunch::GSSlover(double if0)
{

    double fLow;
    double fUp;
    double fMid;
    double rangeLow;
    double rangeUP;
    double mid;
    
    rangeLow = 0.E0;
    rangeUP  = 8.E0;    // seems like 8 rms emittance truncated..
    
    fLow = 1 - if0 - (1+rangeLow) * exp(-rangeLow);    // compared with yuri' equation in NIMA BeamPath. 
    fUp  = 1 - if0 - (1+rangeUP ) * exp(-rangeUP );


    while(rangeUP - rangeLow > 1.0E-5 )
    {
        mid  = (rangeLow + rangeUP)/2.E0;
        fMid = 1 - if0 - (1+mid) * exp(-mid);

        if(abs(fMid)<1.0E-6)
        {
            return mid;
        }

        if(fMid<0)
        {
            rangeLow = mid; 
        }
        else if (fMid>0)
        {
            rangeUP = mid;
        }
    }
    
    mid  = (rangeLow + rangeUP)/2.E0;
    return mid;
}



//void Bunch::RMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k)
//{
//
//
//
//    pxAver=0.E0;
//    pyAver=0.E0;
//    xAver =0.E0;
//    yAver =0.E0;
//    
//    macroEleNumSurivePerBunch = 0;
//    
////    xAver   =   accumulate(begin(ePositionX), end(ePositionX), 0.0); 
////    yAver   =   accumulate(begin(ePositionY), end(ePositionY), 0.0); 
////    pxAver  =   accumulate(begin(eMomentumX), end(eMomentumX), 0.0); 
////    pyAver  =   accumulate(begin(eMomentumY), end(eMomentumY), 0.0); 
//    
//
//    for(int i=0;i<macroEleNumPerBunch;i++)
//    {
//        if(eSurive[i]==0)    
//        {
//            continue;
//        }
//        xAver   =  xAver    + ePositionX[i];
//        yAver   =  yAver    + ePositionY[i];
//        pxAver  = pxAver    + eMomentumX[i];
//        pyAver  = pyAver    + eMomentumY[i];
//        
//        macroEleNumSurivePerBunch++;
//    }
//
//    xAver   = xAver   /macroEleNumSurivePerBunch;
//    yAver   = yAver   /macroEleNumSurivePerBunch;
//    pxAver  = pxAver  /macroEleNumSurivePerBunch;
//    pyAver  = pyAver  /macroEleNumSurivePerBunch;
//
//
//
//
//    double x2Aver=0.E0;
//    double y2Aver=0.E0;
//    double px2Aver=0.E0;
//    double py2Aver=0.E0;
//    double xpxAver=0.E0;
//    double ypyAver=0.E0;
//
//// The below section is used to calculate the effective emittance
//
//
//    for(int i=0;i<macroEleNumPerBunch;i++)
//    {
//        if(eSurive[i]==0)    
//        {
//            continue;
//        }
//
//        x2Aver  = x2Aver    + pow(ePositionX[i],2);
//        y2Aver  = y2Aver    + pow(ePositionY[i],2);
//        px2Aver = px2Aver   + pow(eMomentumX[i],2);
//        py2Aver = py2Aver   + pow(eMomentumY[i],2);
//
//        xpxAver = xpxAver   + (ePositionX[i]) * (eMomentumX[i]);
//        ypyAver = ypyAver   + (ePositionY[i]) * (eMomentumY[i]);
//
//    }
//
//
//
//    x2Aver  = x2Aver  /macroEleNumSurivePerBunch;
//    y2Aver  = y2Aver  /macroEleNumSurivePerBunch;
//    px2Aver = px2Aver /macroEleNumSurivePerBunch;
//    py2Aver = py2Aver /macroEleNumSurivePerBunch;
//    xpxAver = xpxAver /macroEleNumSurivePerBunch;
//    ypyAver = ypyAver /macroEleNumSurivePerBunch;
//    
//    rmsEffectiveRingEmitX = sqrt(x2Aver * px2Aver - pow(xpxAver,2));
//    rmsEffectiveRingEmitY = sqrt(y2Aver * py2Aver - pow(ypyAver,2));
//    
//
//    rmsEffectiveRx = sqrt(rmsEffectiveRingEmitX * latticeInterActionPoint.twissBetaX[k]);
//    rmsEffectiveRy = sqrt(rmsEffectiveRingEmitY * latticeInterActionPoint.twissBetaY[k]);
//
//
//
//
//
//// The below section is used to calculate the  rms emittance
//// the obtained rms size is used to calculate the interaction between beam and ion
//
//    x2Aver=0.E0;
//    y2Aver=0.E0;
//    px2Aver=0.E0;
//    py2Aver=0.E0;
//    xpxAver=0.E0;
//    ypyAver=0.E0;
//
//
//    for(int i=0;i<macroEleNumPerBunch;i++)
//    {
//        if(eSurive[i]==0)    
//        {
//            continue;
//        }
//        
//        x2Aver  = x2Aver    + pow(ePositionX[i]-xAver ,2);
//        y2Aver  = y2Aver    + pow(ePositionY[i]-yAver ,2);
//        px2Aver = px2Aver   + pow(eMomentumX[i]-pxAver,2);
//        py2Aver = py2Aver   + pow(eMomentumY[i]-pyAver,2);
//
//        xpxAver = xpxAver   + (ePositionX[i]-xAver) * (eMomentumX[i]-pxAver);
//        ypyAver = ypyAver   + (ePositionY[i]-yAver) * (eMomentumY[i]-pyAver);
//    }
//
//
//
//    x2Aver  = x2Aver  /macroEleNumSurivePerBunch;
//    y2Aver  = y2Aver  /macroEleNumSurivePerBunch;
//    px2Aver = px2Aver /macroEleNumSurivePerBunch;
//    py2Aver = py2Aver /macroEleNumSurivePerBunch;
//    xpxAver = xpxAver /macroEleNumSurivePerBunch;
//    ypyAver = ypyAver /macroEleNumSurivePerBunch;
//
//    rmsEmitX = sqrt(x2Aver * px2Aver - pow(xpxAver,2));
//    rmsEmitY = sqrt(y2Aver * py2Aver - pow(ypyAver,2));
//
//
//    rmsRx = sqrt(rmsEmitX * latticeInterActionPoint.twissBetaX[k]);
//    rmsRy = sqrt(rmsEmitY * latticeInterActionPoint.twissBetaY[k]);
//
//
//
//}



void Bunch::WSRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k)
{

    rmsRx = sqrt(emittanceX * latticeInterActionPoint.twissBetaX[k] +  pow(rmsEnergySpread *latticeInterActionPoint.twissDispX[k],2) );
    rmsRy = sqrt(emittanceY * latticeInterActionPoint.twissBetaY[k] +  pow(rmsEnergySpread *latticeInterActionPoint.twissDispY[k],2) );
    
    xAver = ePositionX[0];
    yAver = ePositionY[0];
    zAver = ePositionZ[0];
    pxAver = eMomentumX[0];
    pyAver = eMomentumY[0];
    pzAver = eMomentumZ[0];

    actionJx = 0.E0;
    actionJy = 0.E0;


    double twissAlphaXTemp = latticeInterActionPoint.twissAlphaX[k];
    double twissAlphaYTemp = latticeInterActionPoint.twissAlphaY[k];
    double twissBetaXTemp  = latticeInterActionPoint.twissBetaX[k];
    double twissBetaYTemp  = latticeInterActionPoint.twissBetaY[k];
    

    actionJx = (1+pow(twissAlphaXTemp,2))/twissBetaXTemp * pow(xAver,2) 
             +  2*twissAlphaXTemp* xAver * pxAver 
             +   twissBetaXTemp * pow(pxAver,2);

    actionJy = (1+pow(twissAlphaYTemp,2))/twissBetaYTemp * pow(yAver,2) 
            + 2*twissAlphaYTemp* yAver * pyAver 
            +   twissBetaYTemp * pow(pyAver,2);
    
    actionJx = actionJx /2.;
    actionJy = actionJy /2.;
		
}


//void Bunch::SSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k)
//{
//
//// (1) get the force of accumulated ions due to the bunch electron beam --- BassettiErskine model.
//// The process to get the force from electron to accumulated ion is the same as a weak-strong model. 
//
//    double nE0;
//    double coeffI;
//
//    nE0 = electronNumPerBunch;
//    coeffI = 2.0*nE0*ElecClassicRadius*ElectronMassEV/IonMassEV/28*CLight;
//
//
//    double tempFx;
//    double tempFy;
//    double posx;
//    double posy;
//    double rmsRxTemp;
//    double rmsRyTemp;
//
//    rmsRxTemp = rmsRx;
//    rmsRyTemp = rmsRy;
//
//
//    for(int j=0;j<latticeInterActionPoint.ionAccumuNumber[k];j++)
//    {
//
//        latticeInterActionPoint.ionAccumuFx[k][j]=0.E0;
//        latticeInterActionPoint.ionAccumuFy[k][j]=0.E0;
//
//
//        posx    =  latticeInterActionPoint.ionAccumuPositionX[k][j] - xAver;
//        posy    =  latticeInterActionPoint.ionAccumuPositionY[k][j] - yAver;
//
//        BassettiErskine(posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy);
//
//        latticeInterActionPoint.ionAccumuFx[k][j]= coeffI*tempFx;
//        latticeInterActionPoint.ionAccumuFy[k][j]= coeffI*tempFy;
//
//    }
////    cout<< latticeInterActionPoint.ionAccumuFx[k][0]<<endl;
////    cout<< latticeInterActionPoint.ionAccumuFy[k][0]<<endl;
////    cout<<__LINE__<<endl;
////    getchar();
//
//
//// (2) get the force of certain bunched beam due to accumulated ions --- BassettiErskine model.
//// The process to get the force from accumulated ion beam is the similar to inverse "strong-weak" model. 
//
//
//    double coeffE;
//    double nI0;
//    nI0 = latticeInterActionPoint.macroIonCharge[k] * latticeInterActionPoint.ionAccumuNumber[k];
//    coeffE = 2.0*nI0*ElecClassicRadius/rGamma; 
//
//
//    rmsRxTemp = latticeInterActionPoint.ionAccumuRMSX[k];
//    rmsRyTemp = latticeInterActionPoint.ionAccumuRMSY[k];
//
//
//
//    for(int i=0;i<macroEleNumPerBunch;i++)
//    {
//
//        eFx[i] =0.E0;
//        eFy[i] =0.E0;
//        posx   =0.E0;
//        posy   =0.E0;
//        tempFx =0.E0;
//        tempFy =0.E0;
//
//
//        if(eSurive[i] == 0)
//        {
//            continue;
//        }
//
//        posx    =  ePositionX[i] - latticeInterActionPoint.ionAccumuAverX[k];
//        posy    =  ePositionY[i] - latticeInterActionPoint.ionAccumuAverY[k];
//
//        BassettiErskine(posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy);
//    
//        eFx[i] =    coeffE * tempFx;
//        eFy[i] =    coeffE * tempFy;
//
//
//    }
////    cout<<eFx[0]<<"  "<<eFx[0]<<" "<<__LINE__<<"  strong strong force of Electron "<<endl;
////    
////    getchar();   
//
//
//
//
//
//
//
//// test of the bassettin formular
///*    ofstream sctestfile("SC_test_1.dat");
//   double sc_x;
//   double sc_y;
//   tempFx=0;
//   tempFy=0;
//   
//   for (int i=-50;i<=50;i++)
//   {
//       for (int j=-25;j<=25;j++)
//           {
//               posx =i;
//               posy =j;
//               
//               BassettiErskine(posx,posy,4.,2.,tempFx,tempFy);
//               
//               sctestfile<<i<<"  "<<j<<"   "<<tempFx<<"   "<<tempFy<<endl;
//           }
//   }
//   
//   cout<<__LINE__<<endl;
//   getchar();
// */
//
////getchar();
//
//
//
//}
//
//
//
void Bunch::WSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k)
{

    double nE0 = electronNumPerBunch;
    double nI0;
    double ionMassNumber;
    double omegaE, omegaI;
    double coeffI, coeffE;
    double tempFx=0.E0;
    double tempFy=0.E0;
    double posx=0.E0;
    double posy=0.E0;
    double rmsRxTemp = rmsRx;
    double rmsRyTemp = rmsRy;
    double eFxTemp;
    double eFyTemp;
    
    eFx[0] =0.E0;
    eFy[0] =0.E0;
    eFxTemp = 0.E0;
    eFyTemp = 0.E0;

    

    for(int p=0;p<latticeInterActionPoint.gasSpec;p++)
    {
        nI0 = latticeInterActionPoint.macroIonCharge[k][p];
        ionMassNumber = latticeInterActionPoint.ionMassNumber[p];
        
        coeffI = 2.0*nE0*ElecClassicRadius*ElectronMassEV/IonMassEV/ionMassNumber * CLight;   // [m * m/s]
        coeffE = 2.0*ElecClassicRadius/rGamma;                                                // [m] 
       
        for(int j=0;j<latticeInterActionPoint.ionAccumuNumber[k][p];j++)
        {
            latticeInterActionPoint.ionAccumuFx[k][p][j]=0.E0;
            latticeInterActionPoint.ionAccumuFy[k][p][j]=0.E0;

            posx    =  latticeInterActionPoint.ionAccumuPositionX[k][p][j] - xAver;
            posy    =  latticeInterActionPoint.ionAccumuPositionY[k][p][j] - yAver;  
            
            if(abs(posx/rmsRxTemp) + abs(posy/rmsRyTemp)<1.0E-5)
            {
                tempFx=0;
                tempFy=0;
            }
            else
            {
                if( (rmsRxTemp - rmsRyTemp)/rmsRyTemp > 1.e-4 )
                {
                    BassettiErskine1(posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy);                   //tempFx ->[1/m]; // assume electron rms as gaussion-- get force at ion position
                }
                else if ( (rmsRxTemp - rmsRyTemp)/rmsRyTemp < -1.e-4  )
                {
                    BassettiErskine1(posy,posx,rmsRyTemp,rmsRxTemp,tempFy,tempFx); 
                }
                else
                {
                    GaussianField(posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy);
                }
            }
                
            latticeInterActionPoint.ionAccumuFx[k][p][j]= coeffI*tempFx;                  //[1/m] * [m * m/s] - > [m/s]; -> (13) integrate along dt gives ion velovity change
            latticeInterActionPoint.ionAccumuFy[k][p][j]= coeffI*tempFy;                  

            eFxTemp = eFxTemp + tempFx;
            eFyTemp = eFyTemp + tempFy;
        }
               
        eFx[0] += -1*eFxTemp * coeffE * nI0;   // since the e and ion with opposite charge state  [1/m * m ]-> [rad]  integrage (12) along ds -- beam dpx change. 
        eFy[0] += -1*eFyTemp * coeffE * nI0;  

    }

                                    
    

    /*
    cout<<rmsRx<<"   "<<rmsRy<<endl;

    ofstream fout ("test.sdds",ios_base::app);  
    for(int i=0;i<latticeInterActionPoint.ionAccumuPositionX[k].size();i++)
    {
        fout<<setw(15)<<left<<latticeInterActionPoint.ionPositionX[k][i]/rmsRx
            <<setw(15)<<left<<latticeInterActionPoint.ionPositionY[k][i]/rmsRy
            <<setw(15)<<left<<latticeInterActionPoint.ionAccumuFx[k][i]
            <<setw(15)<<left<<latticeInterActionPoint.ionAccumuFy[k][i]
            <<endl; 
    }
    cout<<__LINE__<<"   "<<__FILE__<<endl;
    getchar();
    */

    
//    omegaE = coeffE * nI0 * pow(CLight,2) /rmsRy/(rmsRx+rmsRy);
//    omegaE = sqrt(omegaE);
//    
//    omegaI = coeffI * CLight /rmsRy/(rmsRx+rmsRy);
//    omegaI = sqrt(omegaI);




 
 /*
    // to verify the BassettiErskine formula 2020-10-22 chaoli DESY     
   ofstream sctestfile1("SC_test_1.dat");
   double sc_x;
   double sc_y;
   tempFx=0;
   tempFy=0;
   
   for (int i=-30;i<=30;i++)
   {
       for (int j=-20;j<=20;j++)
           {
               posx =i;
               posy =j;               
               GaussianField(posx,posy,2.,2.,tempFx,tempFy);               
               sctestfile1<<posx<<"  "<<posy<<"   "<<tempFx<<"   "<<tempFy <<endl;
           }
   }
   sctestfile1.close();

 
 
   ofstream sctestfile2("SC_test_2.dat");
   
     for (int i=-30;i<=30;i++)
   {
       for (int j=-20;j<=20;j++)
           {
               posx =i;
               posy =j;
               
               BassettiErskine1(posx,posy,2.1,2.,tempFx,tempFy);
               
               sctestfile2<<posx<<"  "<<posy<<"   "<<tempFx<<"   "<<tempFy<<endl;
           }
   }
   sctestfile2.close();
  

      ofstream sctestfile3("SC_test_3.dat");
   
    for (int i=-50;i<=50;i++)
    {
       for (int j=-25;j<=25;j++)
           {
               posx =i;
               posy =j;
               
               BassettiErskine(posx,posy,16.,8.,tempFx,tempFy);
               
               sctestfile3<<i<<"  "<<j<<"   "<<tempFx<<"   "<<tempFy<<"	"<<sqrt(pow(tempFx,2)+pow(tempFy,2))<<endl;
           }
    }
    sctestfile3.close();
      
   cout<<__LINE__<<endl;
   getchar();
 */

}


void Bunch::GaussianField(double posx,double posy,double rmsRxTemp, double rmsRyTemp,double &tempFx,double &tempFy)
{
    double r2 = pow(posx,2) + pow(posy,2);
    double sigma=rmsRxTemp;
    double sigma2 = sigma * sigma; 
    
    tempFx = - posx / r2 * (1 - exp(- r2 /2/sigma2 ) );
    tempFy = - posy / r2 * (1 - exp(- r2 /2/sigma2 ) );

}


void Bunch::BassettiErskine0(double posx,double posy,double rmsRxTemp, double rmsRyTemp,double &tempFx,double &tempFy)
{
    complex<double> z1(0.E0,0.E0);
    complex<double> z2(0.E0,0.E0);
    complex<double> w1(0.E0,0.E0);
    complex<double> w2(0.E0,0.E0);
    double z3;
    double w3;

    complex<double> w0(0.E0,0.E0);
    
    double sigma;
    double ryOverRx;
    double temp;
    double tempPosix;
    double tempPosiy;
    tempFx=0;
    tempFy=0;

    tempPosix=0.E0;
    tempPosiy=0.E0;
    ryOverRx=0.E0;
    sigma   =0.E0;
    temp    =0.E0;
    
    

    if (pow(rmsRxTemp,2)>pow(rmsRyTemp,2))
    {
        tempPosix    =  abs(posx);
        tempPosiy    =  abs(posy);
        
        if(abs(tempPosix/rmsRxTemp) + abs(tempPosiy/rmsRyTemp)<1.0E-6)
        {
            tempFx=0.E0;
            tempFy=0.E0;
        }
        else
        {
            sigma    = sqrt(2*pow(rmsRxTemp,2)-2*pow(rmsRyTemp,2));
            ryOverRx = rmsRyTemp/rmsRxTemp;
        
            double coeffBE;
            coeffBE    = -1 * sqrt(PI)/sigma;


            z1           =  {tempPosix/sigma,tempPosiy/sigma};
            w1           =  Faddeeva::w(z1);

            z2           =  {ryOverRx*tempPosix/sigma,tempPosiy/sigma/ryOverRx};
            w2           =  Faddeeva::w(z2);

            z3           =  - pow(tempPosix/rmsRxTemp,2)/2 - pow(tempPosiy/rmsRyTemp,2)/2;
            w3           =  - exp(z3); 
            
            w0           =  coeffBE * (w1 + w3 * w2 ); 
            
            tempFx       = w0.imag();
            tempFy       = w0.real();
            

            if(posx<=0)
            {
                tempFx  = -tempFx;
            }

            if(posy<=0)
            {
                tempFy = -tempFy; 
            }

        }

    }
    else
    {
        temp  = posx;
        posx  = posy;
        posy  = -temp;

        
        temp        = rmsRyTemp;
        rmsRyTemp   = rmsRxTemp;
        rmsRxTemp   = temp;
        
        
        tempPosix    =  abs(posx);
        tempPosiy    =  abs(posy);
        
        if(abs(tempPosix/rmsRxTemp) + abs(tempPosiy/rmsRyTemp)<1.0E-6)
        {
            tempFx=0.E0;
            tempFy=0.E0;
        }
        else
        {
            sigma    = sqrt(2*pow(rmsRxTemp,2)-2*pow(rmsRyTemp,2));
            ryOverRx = rmsRyTemp/rmsRxTemp;
        


            double coeffBE;
            coeffBE    = -1 * sqrt(PI)/sigma;


            z1           =  {tempPosix/sigma,tempPosiy/sigma};
            w1           =  Faddeeva::w(z1);

            z2           =  {ryOverRx*tempPosix/sigma,tempPosiy/sigma/ryOverRx};
            w2           =  Faddeeva::w(z2);

            z3           =  - pow(tempPosix/rmsRxTemp,2)/2 - pow(tempPosiy/rmsRyTemp,2)/2;
            w3           =  - exp(z3); 
            
            w0           =  coeffBE * (w1 + w3 * w2 ); 
            
            tempFx       = w0.imag();
            tempFy       = w0.real();
            

            if(posx<=0)
            {
                tempFx  = -tempFx;
            }

            if(posy<=0)
            {
                tempFy = -tempFy; 
            }
        
            temp     = tempFy;
            tempFy   = tempFx;
            tempFx   = -temp;
                    
        }               
    }
    
 

}


void Bunch::BassettiErskine1(double posx,double posy,double rmsRxTemp, double rmsRyTemp,double &tempFx,double &tempFy)
{
    complex<double> z1(0.E0,0.E0);
    complex<double> z2(0.E0,0.E0);
    complex<double> w1(0.E0,0.E0);
    complex<double> w2(0.E0,0.E0);
    double z3;
    double w3;

    complex<double> w0(0.E0,0.E0);
    
    double sigma=0.E0;
    double ryOverRx=0.E0;
    double tempPosix=0.E0;
    double tempPosiy=0.E0;
    tempFx=0;
    tempFy=0;
    
    tempPosix    =  abs(posx);
    tempPosiy    =  abs(posy);
    
    sigma    = sqrt(2*pow(rmsRxTemp,2)-2*pow(rmsRyTemp,2));
    ryOverRx = rmsRyTemp/rmsRxTemp;

    double coeffBE;
    coeffBE    = -1 * sqrt(PI)/sigma;

    z1           =  {tempPosix/sigma,tempPosiy/sigma};
    w1           =  Faddeeva::w(z1);

    z2           =  {ryOverRx*tempPosix/sigma,tempPosiy/sigma/ryOverRx};
    w2           =  Faddeeva::w(z2);

    z3           =  - pow(tempPosix/rmsRxTemp,2)/2 - pow(tempPosiy/rmsRyTemp,2)/2;
    w3           =  - exp(z3); 
    
    w0           =  coeffBE * (w1 + w3 * w2 ); 
    
    tempFx       = w0.imag();
    tempFy       = w0.real();    
    

    if(posx<=0)
    {
        tempFx  = -tempFx;
    }

    if(posy<=0)
    {
        tempFy = -tempFy; 
    }

}



void Bunch::BunchTransferDueToIon(LatticeInterActionPoint &latticeInterActionPoint, int k)
{

    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        if(eSurive[i]==0)    
        {
            continue;
        }
      
        eMomentumX[i] +=  eFx[i]  ;    // rad
        eMomentumY[i] +=  eFy[i]  ;    // rad 
    }
}

void Bunch::BunchTransferDueToWake()
{
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        eMomentumX[i] += wakeForceAver[0];
        eMomentumY[i] += wakeForceAver[1];
        eMomentumZ[i] += wakeForceAver[2];
    }
}


void Bunch::BunchTransferDueToLatticeT(LatticeInterActionPoint &latticeInterActionPoint, int k)
{

    double xtemp=0.E0;
    double ytemp=0.E0;
    double xPtemp=0.E0;
    double yPtemp=0.E0;



    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        // refers to SY. Lee Eq. 2.67 

        xtemp  = latticeInterActionPoint.xTransferMatrix[k][0] * ePositionX[i]
               + latticeInterActionPoint.xTransferMatrix[k][1] * eMomentumX[i];
        
        xPtemp = latticeInterActionPoint.xTransferMatrix[k][2] * ePositionX[i]
               + latticeInterActionPoint.xTransferMatrix[k][3] * eMomentumX[i];
                       
        ytemp  = latticeInterActionPoint.yTransferMatrix[k][0] * ePositionY[i]
               + latticeInterActionPoint.yTransferMatrix[k][1] * eMomentumY[i];
        
        yPtemp = latticeInterActionPoint.yTransferMatrix[k][2] * ePositionY[i]
               + latticeInterActionPoint.yTransferMatrix[k][3] * eMomentumY[i];
               
                                               
        ePositionX[i] = xtemp;
        ePositionY[i] = ytemp;
        eMomentumX[i] = xPtemp;        
        eMomentumY[i] = yPtemp;
       
        if(pow(ePositionX[i]/pipeAperatureX,2) + pow(ePositionY[i]/pipeAperatureY,2)>1)
        {
            eSurive[i] = 0;
        }   

    }

}

void Bunch::BunchTransferDueToLatticeL(ReadInputSettings &inputParameter)
{
    double ztemp=0.E0;
    double ttemp=0.E0;
    double zPtemp=0.E0;

    double t0         = inputParameter.ringParBasic->t0;
    double rBeta      = inputParameter.ringParBasic->rBeta;  
    double eta        = inputParameter.ringParBasic->eta;
    double rfBaseFreq = inputParameter.ringParRf->rfBaseFreq;
    double vRF        = inputParameter.ringParRf->vRF;
    double synPhase   = inputParameter.ringParRf->synPhase;
    double u0         = inputParameter.ringParBasic->u0;  
    double electronBeamEnergy = inputParameter.ringParBasic->electronBeamEnergy;


    double zPtemp1 =  u0 / electronBeamEnergy;
    double zPtemp2 ;

    // higher order harmoinc cavities is no included yet. --- key is to 
    double synPhaseH  = inputParameter.ringParRf->synPhaseH;
    double vRFH       = inputParameter.ringParRf->vRFH;
    double harmRF     = inputParameter.ringParRf->harmRF;



    
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
               
        if (harmRF!=0 && harmRF!=1)
        {
            zPtemp2 =  vRF  / electronBeamEnergy * sin (2. * PI *          rfBaseFreq * ePositionZ[i] / ( rBeta * CLight)  +  synPhase) 
                    +  vRFH / electronBeamEnergy * sin (2. * PI * harmRF * rfBaseFreq * ePositionZ[i] / ( rBeta * CLight)  +  synPhaseH); 
                    // refer to Petra4 CDR Eq. 6.7 ~ Eq. 6.10    
        }
        else
        {
            zPtemp2 =  vRF / electronBeamEnergy * sin (2. * PI * rfBaseFreq * ePositionZ[i] / ( rBeta * CLight)  +  synPhase) ;         
        } 
          
        zPtemp  =  eMomentumZ[i]  -  zPtemp1 +  zPtemp2;            //Refer. Nagaoka Eq.(1~4); or S. Y. Lee 3.30 and 3.31        
        ztemp   =  ePositionZ[i]  +  eta * t0 * zPtemp * CLight;    // deltaZ = deltaT*CLight = eta * T0 * deltaPOverP *CLight

        ePositionZ[i] =  ztemp; 
        eMomentumZ[i] =  zPtemp;       
    }

}



void Bunch::BunchSynRadDamping(ReadInputSettings &inputParameter, LatticeInterActionPoint &latticeInterActionPoint)
{

    
    // in the unit of number of truns for synchRadDampTime setting.
    vector<double> synchRadDampTime;
    synchRadDampTime.resize(3);    
    for(int i=0;i<synchRadDampTime.size();i++)
    {
        synchRadDampTime[i] = inputParameter.ringParBasic->synchRadDampTime[i];
    }
        

    int k=0;    
    //	 set the J2_sympeletic matrix
    gsl_matrix * sympleMarixJ = gsl_matrix_alloc (2, 2); 
    gsl_matrix_set_zero(sympleMarixJ);
    gsl_matrix_set (sympleMarixJ, 0, 1, 1);
    gsl_matrix_set (sympleMarixJ, 1, 0,-1);

    //	// set the H matrix
    gsl_matrix * dispMatrixH  = gsl_matrix_alloc (6, 6);
    gsl_matrix * dispMatrixHX = gsl_matrix_alloc (2, 2);
    gsl_matrix * dispMatrixHY = gsl_matrix_alloc (2, 2);

    gsl_matrix * H31   = gsl_matrix_alloc (2, 2); 
    gsl_matrix * H32   = gsl_matrix_alloc (2, 2);  


    gsl_matrix_set_identity(dispMatrixH);
    gsl_matrix_set_zero(dispMatrixHX);
    gsl_matrix_set_zero(dispMatrixHY);	
    gsl_matrix_set_zero(H31);
    gsl_matrix_set_zero(H32);	
	    


    gsl_matrix_set(dispMatrixHX,0,1,latticeInterActionPoint.twissDispX[k]);
    gsl_matrix_set(dispMatrixHX,1,1,latticeInterActionPoint.twissDispPX[k]);

    gsl_matrix_set(dispMatrixHY,0,1,latticeInterActionPoint.twissDispY[k]);
    gsl_matrix_set(dispMatrixHY,1,1,latticeInterActionPoint.twissDispPY[k]);


    gsl_matrix * tempMatrix1 = gsl_matrix_alloc (2, 2);
    gsl_matrix * tempMatrix2 = gsl_matrix_alloc (2, 2);
    gsl_matrix_set_zero(tempMatrix1);
    gsl_matrix_set_zero(tempMatrix2);


    //	// get H31 sub_matrix;
    gsl_matrix_transpose(dispMatrixHX);

    gsl_matrix_mul(dispMatrixHX,sympleMarixJ,tempMatrix1);
    gsl_matrix_mul(sympleMarixJ,tempMatrix1,tempMatrix2);
    gsl_matrix_scale (tempMatrix2, -1);
    gsl_matrix_memcpy (H31, tempMatrix2);

    gsl_matrix_transpose(dispMatrixHX);



    //	// get H32 sub_matrix;
    gsl_matrix_transpose(dispMatrixHY);

    gsl_matrix_mul(dispMatrixHY,sympleMarixJ,tempMatrix1);
    gsl_matrix_mul(sympleMarixJ,tempMatrix1,tempMatrix2);
    gsl_matrix_scale (tempMatrix2, -1);

    gsl_matrix_memcpy (H32, tempMatrix2);	
    gsl_matrix_transpose(dispMatrixHY);


    //	// set H Marix
    gsl_matrix_set(dispMatrixH,0,4,-gsl_matrix_get(dispMatrixHX,0,0));
    gsl_matrix_set(dispMatrixH,0,5,-gsl_matrix_get(dispMatrixHX,0,1));
    gsl_matrix_set(dispMatrixH,1,4,-gsl_matrix_get(dispMatrixHX,1,0));
    gsl_matrix_set(dispMatrixH,1,5,-gsl_matrix_get(dispMatrixHX,1,1));
    
    gsl_matrix_set(dispMatrixH,2,4,-gsl_matrix_get(dispMatrixHY,0,0));
    gsl_matrix_set(dispMatrixH,2,5,-gsl_matrix_get(dispMatrixHY,0,1));
    gsl_matrix_set(dispMatrixH,3,4,-gsl_matrix_get(dispMatrixHY,1,0));
    gsl_matrix_set(dispMatrixH,3,5,-gsl_matrix_get(dispMatrixHY,1,1));

    gsl_matrix_set(dispMatrixH,4,0,gsl_matrix_get(H31,0,0));
    gsl_matrix_set(dispMatrixH,5,0,gsl_matrix_get(H31,1,0));
    gsl_matrix_set(dispMatrixH,4,1,gsl_matrix_get(H31,0,1));
    gsl_matrix_set(dispMatrixH,5,1,gsl_matrix_get(H31,1,1));

    gsl_matrix_set(dispMatrixH,4,2,gsl_matrix_get(H32,0,0));
    gsl_matrix_set(dispMatrixH,5,2,gsl_matrix_get(H32,1,0));
    gsl_matrix_set(dispMatrixH,4,3,gsl_matrix_get(H32,0,1));
    gsl_matrix_set(dispMatrixH,5,3,gsl_matrix_get(H32,1,1));


   
    //set Teng Marrix R

    gsl_matrix * tengMatrixR  = gsl_matrix_alloc (6, 6);
    gsl_matrix_set_identity(tengMatrixR);

    gsl_matrix * R2  = gsl_matrix_alloc (2, 2);     //coupling matrix Eq.(6)
//    gsl_matrix_set_identity(R2);                  // full coupling
    gsl_matrix_set_zero(R2);                        // no coupling at the point in SR calculation

    double b = get_det(R2);
    b = sqrt(1-b);

    gsl_matrix * R12  = gsl_matrix_alloc (2, 2);
    gsl_matrix * R21  = gsl_matrix_alloc (2, 2);


    //	//set R12
    gsl_matrix_transpose(R2);
    gsl_matrix_mul(R2,sympleMarixJ,tempMatrix1);
    gsl_matrix_mul(sympleMarixJ,tempMatrix1,tempMatrix2);
    gsl_matrix_memcpy (R12, tempMatrix2);
    gsl_matrix_transpose(R2);
    //set R21
    gsl_matrix_memcpy (R21, R2);


    gsl_matrix_set(tengMatrixR,0,0,b);
    gsl_matrix_set(tengMatrixR,1,1,b);
    gsl_matrix_set(tengMatrixR,2,2,b);
    gsl_matrix_set(tengMatrixR,3,3,b);

    gsl_matrix_set(tengMatrixR,0,2,gsl_matrix_get(R12,0,0));
    gsl_matrix_set(tengMatrixR,0,3,gsl_matrix_get(R12,0,1));
    gsl_matrix_set(tengMatrixR,1,2,gsl_matrix_get(R12,1,0));
    gsl_matrix_set(tengMatrixR,1,3,gsl_matrix_get(R12,1,1));

    gsl_matrix_set(tengMatrixR,2,0,gsl_matrix_get(R21,0,0));
    gsl_matrix_set(tengMatrixR,2,1,gsl_matrix_get(R21,0,1));
    gsl_matrix_set(tengMatrixR,3,0,gsl_matrix_get(R21,1,0));
    gsl_matrix_set(tengMatrixR,3,1,gsl_matrix_get(R21,1,1));

    
    
    //set Twiss B Matrix  
    gsl_matrix * twissMatrixB = gsl_matrix_alloc (6, 6);
    gsl_matrix_set_zero(twissMatrixB);

    
    double a00,a01,a10,a11;

    a00 = 1.0/sqrt(latticeInterActionPoint.twissBetaX[k]);
    a10 = latticeInterActionPoint.twissAlphaX[k] / sqrt(latticeInterActionPoint.twissBetaX[k]);
    a11 = sqrt(latticeInterActionPoint.twissBetaX[k]);


    gsl_matrix_set(twissMatrixB,0,0,a00);
    gsl_matrix_set(twissMatrixB,1,0,a10);
    gsl_matrix_set(twissMatrixB,1,1,a11);


    a00 = 1.0/sqrt(latticeInterActionPoint.twissBetaY[k]);
    a10 = latticeInterActionPoint.twissAlphaY[k] / sqrt(latticeInterActionPoint.twissBetaY[k]);
    a11 = sqrt(latticeInterActionPoint.twissBetaY[k]);


    gsl_matrix_set(twissMatrixB,2,2,a00);
    gsl_matrix_set(twissMatrixB,3,2,a10);
    gsl_matrix_set(twissMatrixB,3,3,a11);


    a00 = 1.0/sqrt(latticeInterActionPoint.twissBetaZ[k]);
    a10 = latticeInterActionPoint.twissAlphaZ[k] / sqrt(latticeInterActionPoint.twissBetaZ[k]);
    a11 = sqrt(latticeInterActionPoint.twissBetaZ[k]);


    gsl_matrix_set(twissMatrixB,4,4,a00);
    gsl_matrix_set(twissMatrixB,5,4,a10);
    gsl_matrix_set(twissMatrixB,5,5,a11);

   
    gsl_matrix * cordTransfer  = gsl_matrix_alloc (6, 6);
    gsl_matrix * invCordTransfer  = gsl_matrix_alloc (6, 6);
    gsl_matrix * cordTransferTemp  = gsl_matrix_alloc (6, 6);
    gsl_matrix_mul(tengMatrixR,dispMatrixH,cordTransferTemp);
    gsl_matrix_mul(twissMatrixB,cordTransferTemp,cordTransfer);
    

    //get the   invCordTransfer =  (cordTransfer)^-1 ---- approach 1 

    gsl_matrix *  sympleMarixJ6 = gsl_matrix_alloc (6, 6);
    gsl_matrix_set_zero(sympleMarixJ6);
    gsl_matrix_set(sympleMarixJ6,0,1, 1);
    gsl_matrix_set(sympleMarixJ6,1,0,-1);
    gsl_matrix_set(sympleMarixJ6,2,3, 1);
    gsl_matrix_set(sympleMarixJ6,3,2,-1);
    gsl_matrix_set(sympleMarixJ6,4,5, 1);
    gsl_matrix_set(sympleMarixJ6,5,4,-1);

    gsl_matrix_transpose(cordTransfer);

    gsl_matrix_mul(cordTransfer,sympleMarixJ6,cordTransferTemp);
    gsl_matrix_mul(sympleMarixJ6,cordTransferTemp,invCordTransfer);
    gsl_matrix_scale (invCordTransfer, -1);
    
    gsl_matrix_transpose(cordTransfer);
    //  end of approach 1 -------------------------------------------

    //  get the   invCordTransfer =  (cordTransfer)^-1 ---- approach 2  in a general sense but slower 
    //  gsl_matrix_memcpy(invCordTransfer,cordTransfer);
    //  gsl_matrix_inv(invCordTransfer);
    ////end of approach 2 -------------------------------------------

    // Ref. Zhang Yuan PRAB 8 074402

        
    double lambda[3];	
    double coeff[3];

    
    for(int i=0;i<3;i++)
    { 
        lambda[i] = exp(-1.0/synchRadDampTime[i]);
    }


    
    coeff[0] = sqrt(1 - pow(lambda[0],2)) * sqrt(emittanceX); 
    coeff[1] = sqrt(1 - pow(lambda[1],2)) * sqrt(emittanceY);
    coeff[2] = sqrt(1 - pow(lambda[2],2)) * sqrt(emittanceZ);
    
    gsl_matrix * vecX   = gsl_matrix_alloc (6, 1);
    gsl_matrix * vecNX  = gsl_matrix_alloc (6, 1);

    
    double tempX,tempPX,tempY,tempPY,tempZ,tempPZ;
    double randR[6];

    std::random_device rd{};
    std::mt19937 gen{rd()};

    std::normal_distribution<> dx{0,1};
    
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        //(0) transfer x to X;
        gsl_matrix_set(vecX,0,0,ePositionX[i]);
        gsl_matrix_set(vecX,1,0,eMomentumX[i]);
        gsl_matrix_set(vecX,2,0,ePositionY[i]);
        gsl_matrix_set(vecX,3,0,eMomentumY[i]);
        gsl_matrix_set(vecX,4,0,ePositionZ[i]);
        gsl_matrix_set(vecX,5,0,eMomentumZ[i]);
	    
//		for(int i=0;i<6;i++)
//		{
//			cout<<gsl_matrix_get(vecX,i,0)<<endl;
//		}
//		
//		cout<<endl;
//		
//		for(int i=0;i<3;i++)
//		{
//			cout<<lambda[i]<<endl;
//		}
//		cout<<endl;
		
		gsl_matrix_mul(cordTransfer,vecX,vecNX);

        tempX  = gsl_matrix_get(vecNX,0,0);
        tempPX = gsl_matrix_get(vecNX,1,0);
        tempY  = gsl_matrix_get(vecNX,2,0);
        tempPY = gsl_matrix_get(vecNX,3,0);
        tempZ  = gsl_matrix_get(vecNX,4,0);
        tempPZ = gsl_matrix_get(vecNX,5,0);	

		//(1)   Eq(8~10) dampling and excitation


//		for(int i=0;i<6;i++)
//		{
//			cout<<gsl_matrix_get(vecNX,i,0)<<endl;
//		}
//		cout<<endl;
		
        // synchron-radiation-damping
        tempX  = tempX  * lambda[0] ;
        tempPX = tempPX * lambda[0] ; 
        tempY  = tempY  * lambda[1] ;
        tempPY = tempPY * lambda[1] ; 
        tempZ  = tempZ  			;
        tempPZ = tempPZ * pow(lambda[2],2) ; 


        // synchron-radiation-excitiation only in multi-particle case
        if(macroEleNumPerBunch!=1)
        {
            for(int j=0;j<6;j++)
            {
                //randR[j] = dx(gen);	
                //if(j==4)
                //{
                //    randR[j] = 0.E0;
                //}
                randR[j]=Gaussrand(1,0,100.0*j);
                if(j==4)
                {
                    randR[j]=0.E0;
                }
                
                cout<<randR[j]<<endl;
            }

        tempX  = tempX   + coeff[0] * randR[0];
        tempPX = tempPX  + coeff[0] * randR[1]; 
        tempY  = tempY   + coeff[1] * randR[2];
        tempPY = tempPY  + coeff[1] * randR[3]; 
        tempZ  = tempZ   + coeff[2] * randR[4];
        tempPZ = tempPZ  + coeff[2] * randR[5]; 
        }

                 
        gsl_matrix_set(vecNX,0,0,tempX);
        gsl_matrix_set(vecNX,1,0,tempPX);
        gsl_matrix_set(vecNX,2,0,tempY);
        gsl_matrix_set(vecNX,3,0,tempPY);
        gsl_matrix_set(vecNX,4,0,tempZ);
        gsl_matrix_set(vecNX,5,0,tempPZ);

		//(2) transfer X to x,  Eq(11)
		
//		for(int i=0;i<6;i++)
//		{
//			cout<<gsl_matrix_get(vecNX,i,0)<<endl;
//		}
//		cout<<endl;

        gsl_matrix_mul(invCordTransfer,vecNX,vecX);         

//		for(int i=0;i<6;i++)
//		{
//			cout<<gsl_matrix_get(vecX,i,0)<<endl;
//		}
//			 
//		getchar(); 
         
         
        ePositionX[i] = gsl_matrix_get(vecX,0,0);
        eMomentumX[i] = gsl_matrix_get(vecX,1,0);
        ePositionY[i] = gsl_matrix_get(vecX,2,0);
        eMomentumY[i] = gsl_matrix_get(vecX,3,0);
        ePositionZ[i] = gsl_matrix_get(vecX,4,0);
        eMomentumZ[i] = gsl_matrix_get(vecX,5,0);

        if(pow(ePositionX[i]/pipeAperatureX,2) + pow(ePositionY[i]/pipeAperatureY,2)>1)
        {
            eSurive[i] = 0;
        }    

        
    }
    
    

    gsl_matrix_free (sympleMarixJ);
    gsl_matrix_free (dispMatrixH);
    gsl_matrix_free (dispMatrixHX);
    gsl_matrix_free (dispMatrixHY);
    gsl_matrix_free (H31);
    gsl_matrix_free (H32);
    gsl_matrix_free (tempMatrix1);
    gsl_matrix_free (tempMatrix2);
    
    gsl_matrix_free (tengMatrixR);
    gsl_matrix_free (R2);
    gsl_matrix_free (R21);
    gsl_matrix_free (R12);

   
    gsl_matrix_free (twissMatrixB);

    gsl_matrix_free (cordTransfer);
    gsl_matrix_free (invCordTransfer);
    gsl_matrix_free (cordTransferTemp);
    gsl_matrix_free (sympleMarixJ6);

     
    gsl_matrix_free (vecX);
    gsl_matrix_free (vecNX);
      
    
}




//void Bunch::InitialBeamCurDenZProf()
//{
//    
////    std::vector<double>::iterator zMax = std::max_element(std::begin(ePositionZ), std::end(ePositionZ));
////    std::cout << "Max element is " << *zMax<< " at position " << std::distance(std::begin(ePositionZ), zMax) << std::endl;
////    
////    
////    std::vector<double>::iterator zMin = std::min_element(std::begin(ePositionZ), std::end(ePositionZ));
////    std::cout << "Max element is " << *zMin<< " at position " << std::distance(std::begin(ePositionZ), zMin) << std::endl;
//
////    double rangeBeamSizeZ;
////    if(abs(*zMin)>abs(*zMax))
////    {
////        rangeBeamSizeZ = abs(*zMin);
////    }
////    else
////    {
////        rangeBeamSizeZ = abs(*zMax);
////    }
//
////    double binZSize;
////    
////    binZSize = 2 * rangeBeamSizeZ / bunchBinNumberZ;
////    
////    for(int i=0;i<beamCurDenZProf.size();i++)
////    {
////        beamCurDenZProf[i]=0.E0 ;
////    }
//
//
////    ofstream fout("beamCurDenZProf.dat");
//
//
////    double tempDeltaz;
////    double poszTemp;
////    int poszTempIndex;
////    for(int i=0;i<macroEleNumPerBunch;i++)
////    {
////        poszTemp         = abs(rangeBeamSizeZ-ePositionZ[i]) / binZSize;
////        poszTempIndex    = floor(poszTemp);
//////        tempDeltaz       = poszTemp - poszTempIndex;
//
////        beamCurDenZProf[poszTempIndex  ]  = beamCurDenZProf[poszTempIndex  ]+1-tempDeltaz;
////        beamCurDenZProf[poszTempIndex+1]  = beamCurDenZProf[poszTempIndex+1]+  tempDeltaz;
//
//////        beamCurDenZProf[poszTempIndex]  =   beamCurDenZProf[poszTempIndex] +1;
////    }
//
////    for(int i=0;i<beamCurDenZProf.size();i++)
////    {
////        fout<< i<<" "<<beamCurDenZProf[i]<<endl ;
////    }
////    
////    cout<<"TTTTTT"<<endl;
////    getchar();
////    fout.close();
//    
//}
//
//void Bunch::InitialBeamCurDenTProf()
//{
//
//}





