#include "Bunch.h"
#include "Global.h"
#include "Faddeeva.h"
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <complex>
#include <iostream>
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


    macroEleNumPerBunch = inputParameter.macroEleNumPerBunch;                    // macroEleNumPerBunch=1 -- electron beam weak-strong model
                                                                                  // macroEleNumPerBunch!=1 -- electron beam strong-strong model
    kappa =  inputParameter.kappa;
    pipeAperatureR =   inputParameter.pipeAperatureR;
    PipeAperatureX =   inputParameter.pipeAperatureX;
    PipeAperatureY =   inputParameter.pipeAperatureY;                            

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
    initialDisDx = inputParameter.initialDisDx;
    initialDisDy = inputParameter.initialDisDy;

    
    bunchBinNumberZ = inputParameter.bunchBinNumberZ;
    beamCurDenZProf.resize(bunchBinNumberZ);
    beamCurDenTProf.resize(bunchBinNumberT);
    
    current = inputParameter.current;
    electronEnergy  = inputParameter.electronBeamEnergy;
    
    
    
    emittanceX      = inputParameter.emittanceX;
    emittanceY      = inputParameter.emittanceY;

    rmsBunchLength  = inputParameter.rmsBunchLength;        
    rmsEnergySpread = inputParameter.rmsEnergySpread;
    

    rGamma = electronEnergy/ElectronMassEV;
    rBeta  = sqrt(1.E0-1.E0/pow(rGamma,2));


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


    electronNumPerBunch = current / inputParameter. omegas * 2*PI /inputParameter.totBunchNumber/ElectronCharge;

    macroEleCharge = electronNumPerBunch / macroEleNumPerBunch;
    distributionType = inputParameter.distributionType;

}



void Bunch::DistriGenerator(LatticeInterActionPoint &latticeInterActionPoint,int randomIndex)
{



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

    alphaX = latticeInterActionPoint.twissAlphaX[0];   // here only the data at the first points is used. 
    betaX  = latticeInterActionPoint.twissBetaX[0];
    
    alphaY = latticeInterActionPoint.twissAlphaY[0];
    betaY  = latticeInterActionPoint.twissBetaY[0];
    
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


//    cout<<rX<<" "<<rY<<endl;
//    cout<<emittanceX<<" "<<betaX<<endl;
//    cout<<emittanceY<<" "<<betaY<<endl;

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

    double disDx;
    double disDy;

    srand(time(0)+randomIndex);
    
    disDx = 2*(double(std::rand())/RAND_MAX - 0.5) * initialDisDx;
    disDy = 2*(double(std::rand())/RAND_MAX - 0.5) * initialDisDy;

    int i=1;   
               // the first is set as reference particle 
               // ePositionX[0]=0
               // ePositionY[0]=0
               // eMomentumX[0]=0
               // eMomentumY[0]=0
    
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
                f0  = 2*emittanceX;
                if0 = double(std::rand())/RAND_MAX;
                fi  = GSSlover(if0);
                fi  = fi * f0 ;
                break;
            default:
            cerr<<"Do Nothing, no distribution type in generation.  "<<endl;
        }
        

        ix   = double(rand())/RAND_MAX;
        
        axax =  fi* ix;
        ayay =  (fi - axax) *kappa;  
        
        ax  = sqrt(axax);
        ay  = sqrt(ayay);
 

        phaseX = 2 * PI* double(rand())/RAND_MAX;
        phaseY = 2 * PI* double(rand())/RAND_MAX;





        ePositionX[i] = ax *   sigmaX * cos( phaseX ) + disDx;
        ePositionY[i] = ay *   sigmaY * cos( phaseY ) + disDy;
        eMomentumX[i] = ax * (dSigmaX * cos( phaseX ) - sin(phaseX)/sigmaX );
        eMomentumY[i] = ay * (dSigmaY * cos( phaseY ) - sin(phaseY)/sigmaY );



        // to suppress the numerical noise on average calculation
        ePositionX[i+1] = -ePositionX[i];
        ePositionY[i+1] = -ePositionY[i];
        eMomentumX[i+1] = -eMomentumX[i];
        eMomentumY[i+1] = -eMomentumY[i];
        
    

    

        if(distributionType==3)
        {
            if(pow(ePositionX[i]/rX,2) + pow(ePositionY[i]/rY,2)>9)
            {
                continue;
            }
            
            if(pow(ePositionX[i+1]/rX,2) + pow(ePositionY[i+1]/rY,2)>9)
            {
                continue;
            }
        }

        i=i+2;
    }
    
    
    


//    double tempAverXTest=0;
//    double tempAverYTest=0;
//    double tempAverpXTest=0;
//    double tempAverpYTest=0;
//    tempAverXTest  =  accumulate(begin(ePositionX), end(ePositionX), 0.0); 
//    tempAverYTest  =  accumulate(begin(ePositionY), end(ePositionY), 0.0); 
//    tempAverpXTest =  accumulate(begin(eMomentumX), end(eMomentumX), 0.0); 
//    tempAverpYTest =  accumulate(begin(eMomentumY), end(eMomentumY), 0.0); 
//    tempAverXTest  =  tempAverXTest/1000;
//    tempAverYTest  =  tempAverYTest/1000;
//    tempAverpXTest =  tempAverpXTest/1000;
//    tempAverpYTest =  tempAverpYTest/1000;


//    cout<<__LINE__<<"   "
//        <<tempAverXTest<<"  "
//        <<tempAverYTest<<"  "
//        <<tempAverpXTest<<" "
//        <<tempAverpYTest<<endl;
//        
//        cout<<ePositionX[0]<<"  "<<ePositionX[1]<<endl;
//        getchar();

//    ofstream fout("dist_test.dat");
//    for(int i=0;i<macroEleNumPerBunch;i++)
//    {

//        fout<<i<<"     "
//            <<ePositionX[i]  <<"     "
//            <<ePositionY[i]  <<"     "
//            <<eMomentumX[i]<<"     "
//            <<eMomentumY[i]<<"     "
//            <<ax<<"     "
//            <<ay<<"     "
//            <<sigmaX<<"     "
//            <<sigmaY<<"     "
//            <<dSigmaX<<"     "
//            <<dSigmaY<<"     "
//            <<endl;
//    
//    }
//    
//    fout.close();
//    cout<<"test"<<endl;
//    getchar();




    std::random_device rd{};
    std::mt19937 gen{rd()};
 

    std::normal_distribution<> dx{0,rmsBunchLength};
    std::normal_distribution<> dy{0,rmsEnergySpread};


    double tempx;
    double tempy;

    i=0;
    while(i<macroEleNumPerBunch)
    {
        tempx = dx(gen);
        tempy = dy(gen);

        // Longitudinally, beam is 3 sigma truncted.
        if( abs((tempx-rmsBunchLength)/rmsBunchLength)<3 || abs((tempy-rmsEnergySpread)/rmsEnergySpread)<3)  
        {
            ePositionZ[i]=tempx;
            eMomentumZ[i]=tempy;
        }
        else
        {
            continue;
        }

        i++;
    }

    
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

void Bunch::InducedIonDensity(LatticeInterActionPoint &latticeInterActionPoint)
{


    for(int i=0;i<latticeInterActionPoint.numberOfInteraction;i++)
    {
        latticeInterActionPoint.ionLineDensity[i] = CorssSectionEI*latticeInterActionPoint.vacuumPressure[i]
                                    /latticeInterActionPoint.temperature[i]/Boltzmann*electronNumPerBunch;
    }
    
    
    
}

void Bunch::RMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k)
{



    pxAver=0.E0;
    pyAver=0.E0;
    xAver =0.E0;
    yAver =0.E0;
    
    macroEleNumSurivePerBunch = 0;
    
//    xAver   =   accumulate(begin(ePositionX), end(ePositionX), 0.0); 
//    yAver   =   accumulate(begin(ePositionY), end(ePositionY), 0.0); 
//    pxAver  =   accumulate(begin(eMomentumX), end(eMomentumX), 0.0); 
//    pyAver  =   accumulate(begin(eMomentumY), end(eMomentumY), 0.0); 
    

    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        if(eSurive[i]==0)    
        {
            continue;
        }
        xAver   =  xAver    + ePositionX[i];
        yAver   =  yAver    + ePositionY[i];
        pxAver  = pxAver    + eMomentumX[i];
        pyAver  = pyAver    + eMomentumY[i];
        
        macroEleNumSurivePerBunch++;
    }

    xAver   = xAver   /macroEleNumSurivePerBunch;
    yAver   = yAver   /macroEleNumSurivePerBunch;
    pxAver  = pxAver  /macroEleNumSurivePerBunch;
    pyAver  = pyAver  /macroEleNumSurivePerBunch;




    double x2Aver=0.E0;
    double y2Aver=0.E0;
    double px2Aver=0.E0;
    double py2Aver=0.E0;
    double xpxAver=0.E0;
    double ypyAver=0.E0;

// The below section is used to calculate the center rms emittance

//    for(int i=0;i<macroEleNumPerBunch;i++)
//    {
//        if(eSurive[i]==0)    
//        {
//            continue;
//        }
//        x2Aver  = x2Aver    + pow(ePositionX[i]-xAver ,2);
//        y2Aver  = y2Aver    + pow(ePositionY[i]-yAver ,2);
//        px2Aver = px2Aver   + pow(eMomentumX[i]-pxAver,2);
//        py2Aver = py2Aver   + pow(eMomentumY[i]-pyAver,2);

//        xpxAver = xpxAver   + (ePositionX[i]-xAver) * (eMomentumX[i]-pxAver);
//        ypyAver = ypyAver   + (ePositionY[i]-yAver) * (eMomentumY[i]-pyAver);
//    }



//    x2Aver  = x2Aver  /macroEleNumSurivePerBunch;
//    y2Aver  = y2Aver  /macroEleNumSurivePerBunch;
//    px2Aver = px2Aver /macroEleNumSurivePerBunch;
//    py2Aver = py2Aver /macroEleNumSurivePerBunch;
//    xpxAver = xpxAver /macroEleNumSurivePerBunch;
//    ypyAver = ypyAver /macroEleNumSurivePerBunch;
    

// The below section is used to calculate the effective emittance

    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        if(eSurive[i]==0)    
        {
            continue;
        }
        x2Aver  = x2Aver    + pow(ePositionX[i],2);
        y2Aver  = y2Aver    + pow(ePositionY[i],2);
        px2Aver = px2Aver   + pow(eMomentumX[i],2);
        py2Aver = py2Aver   + pow(eMomentumY[i],2);

        xpxAver = xpxAver   + (ePositionX[i]) * (eMomentumX[i]);
        ypyAver = ypyAver   + (ePositionY[i]) * (eMomentumY[i]);
    }



    x2Aver  = x2Aver  /macroEleNumSurivePerBunch;
    y2Aver  = y2Aver  /macroEleNumSurivePerBunch;
    px2Aver = px2Aver /macroEleNumSurivePerBunch;
    py2Aver = py2Aver /macroEleNumSurivePerBunch;
    xpxAver = xpxAver /macroEleNumSurivePerBunch;
    ypyAver = ypyAver /macroEleNumSurivePerBunch;

    rmsEmitX = sqrt(x2Aver * px2Aver - pow(xpxAver,2));
    rmsEmitY = sqrt(y2Aver * py2Aver - pow(ypyAver,2));


    rmsRx = sqrt(rmsEmitX * latticeInterActionPoint.twissBetaX[k]);
    rmsRy = sqrt(rmsEmitY * latticeInterActionPoint.twissBetaY[k]);

}




void Bunch::WSRMSCal(LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    // in the calculation In ZhangYuan's Equation (8,9,19) what is the value of emittancex? 
    // is it obtained from the rms calculation of beam distribution or what else? Seems do not affect too much 

    rmsRx = sqrt(emittanceX * latticeInterActionPoint.twissBetaX[k]);
    rmsRy = sqrt(emittanceY * latticeInterActionPoint.twissBetaY[k]);

    xAver = ePositionX[0];
    yAver = ePositionY[0];

    pxAver = eMomentumX[0];
    pyAver = eMomentumY[0];

    actionJx = 0.E0;
    actionJy = 0.E0;

    double twissAlphaXTemp;
    double twissAlphaYTemp;
    double twissBetaXTemp;
    double twissBetaYTemp;
    
    twissAlphaXTemp = latticeInterActionPoint.twissAlphaX[k];
    twissAlphaYTemp = latticeInterActionPoint.twissAlphaY[k];
    twissBetaXTemp  = latticeInterActionPoint.twissBetaX[k];
    twissBetaYTemp  = latticeInterActionPoint.twissBetaY[k];


    actionJx = (1+pow(twissAlphaXTemp,2))/twissBetaXTemp * pow(xAver,2) 
             + 2*twissAlphaXTemp* xAver * pxAver 
             +   twissBetaXTemp * pow(pxAver,2);


    actionJy = (1+pow(twissAlphaYTemp,2))/twissBetaYTemp * pow(yAver,2) 
            + 2*twissAlphaYTemp* yAver * pyAver 
            +   twissBetaYTemp * pow(pyAver,2);
    
    actionJx = actionJx /2.;
    actionJy = actionJy /2.;
    
}


void Bunch::SSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k)
{

// (1) get the force of accumulated ions due to the bunch electron beam --- BassettiErskine model.
// The process to get the force from electron to accumulated ion is the same as a weak-strong model. 

    double nE0;
    double coeffI;

    nE0 = electronNumPerBunch;
    coeffI = 2.0*nE0*ElecClassicRadius*ElectronMassEV/IonMassEV/28*CLight;


    double tempFx;
    double tempFy;
    double posx;
    double posy;
    double rmsRxTemp;
    double rmsRyTemp;

    rmsRxTemp = rmsRx;
    rmsRyTemp = rmsRy;

    
 

    for(int j=0;j<latticeInterActionPoint.ionAccumuNumber[k];j++)
    {

        latticeInterActionPoint.ionAccumuFx[k][j]=0.E0;
        latticeInterActionPoint.ionAccumuFy[k][j]=0.E0;



        posx    =  latticeInterActionPoint.ionAccumuPositionX[k][j] - xAver;
        posy    =  latticeInterActionPoint.ionAccumuPositionY[k][j] - yAver;

        BassettiErskine(posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy);

        latticeInterActionPoint.ionAccumuFx[k][j]= coeffI*tempFx;
        latticeInterActionPoint.ionAccumuFy[k][j]= coeffI*tempFy;

    }
//    cout<< latticeInterActionPoint.ionAccumuFx[k][0]<<endl;
//    cout<< latticeInterActionPoint.ionAccumuFy[k][0]<<endl;
//    cout<<__LINE__<<endl;
//    getchar();


// (2) get the force of certain bunched beam due to accumulated ions --- BassettiErskine model.
// The process to get the force from accumulated ion beam is the similar to inverse "strong-weak" model. 


    double coeffE;
    double nI0;
    nI0 = latticeInterActionPoint.macroIonCharge[k] * latticeInterActionPoint.ionAccumuNumber[k];
    coeffE = 2.0*nI0*ElecClassicRadius/rGamma; 


    rmsRxTemp = latticeInterActionPoint.ionAccumuRMSX[k];
    rmsRyTemp = latticeInterActionPoint.ionAccumuRMSY[k];



    for(int i=0;i<macroEleNumPerBunch;i++)
    {

        eFx[i] =0.E0;
        eFy[i] =0.E0;
        posx   =0.E0;
        posy   =0.E0;
        tempFx =0.E0;
        tempFy =0.E0;
        

        if(eSurive[i] == 0)
        {
            continue;
        }

        posx    =  ePositionX[i] - latticeInterActionPoint.ionAccumuAverX[k];
        posy    =  ePositionY[i] - latticeInterActionPoint.ionAccumuAverY[k];

        BassettiErskine(posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy);
    
        eFx[i] =    coeffE * tempFx;
        eFy[i] =    coeffE * tempFy;


    }
//    cout<<eFx[0]<<"  "<<eFx[0]<<" "<<__LINE__<<"  strong strong force of Electron "<<endl;
//    
//    getchar();   







// test of the bassettin formular
//    ofstream sctestfile("SC_test.dat");
//    double sc_x;
//    double sc_y;
//    tempFx=0;
//    tempFy=0;
//    
//    for (int i=-25;i<=25;i++)
//    {
//        for (int j=-50;j<=50;j++)
//            {
//                posx =i;
//                posy =j;
//                
//                BassettiErskine(posx,posy,2.,4.,tempFx,tempFy);
//                
//                sctestfile<<i<<"  "<<j<<"   "<<tempFx<<"   "<<tempFy<<endl;
//            }
//    }
//    
//    cout<<__LINE__<<endl;
//    getchar();






}




void Bunch::WSIonBunchInteraction(LatticeInterActionPoint &latticeInterActionPoint, int k)
{
    double nE0;
    double nI0;

    nE0 = electronNumPerBunch;
    nI0 = latticeInterActionPoint.macroIonCharge[k];
    

    
    omegaE =0.E0;
    omegaI =0.E0;

    double coeffI;
    double coeffE;
    coeffI = 2.0*nE0*ElecClassicRadius*ElectronMassEV/IonMassEV/28*CLight;
    coeffE = 2.0*ElecClassicRadius/rGamma;
    // refers to the Eq. 13 of Ohmi's PRE paper and EPAC 2008 Guoqing Xia, MOPPO68  F(x,y) = -2*N_0*r_e*me*c^2*f(x,y)   
    // there equations are wrong, need to be modifiled. 


    
    double tempFx=0.E0;
    double tempFy=0.E0;
    double posx=0.E0;
    double posy=0.E0;
    double rmsRxTemp;
    double rmsRyTemp;
    double eFxTemp;
    double eFyTemp;
    rmsRxTemp = rmsRx;
    rmsRyTemp = rmsRy;




    eFx[0] =0.E0;
    eFy[0] =0.E0;
    eFxTemp = 0.E0;
    eFyTemp = 0.E0;



    for(int j=0;j<latticeInterActionPoint.ionAccumuNumber[k];j++)
    {

        latticeInterActionPoint.ionAccumuFx[k][j]=0.E0;
        latticeInterActionPoint.ionAccumuFy[k][j]=0.E0;


        posx    =  latticeInterActionPoint.ionAccumuPositionX[k][j] - xAver;
        posy    =  latticeInterActionPoint.ionAccumuPositionY[k][j] - yAver;

        BassettiErskine(posx,posy,rmsRxTemp,rmsRyTemp,tempFx,tempFy);

        latticeInterActionPoint.ionAccumuFx[k][j]= coeffI*tempFx;
        latticeInterActionPoint.ionAccumuFy[k][j]= coeffI*tempFy;

        eFxTemp = eFxTemp + tempFx;
        eFyTemp = eFyTemp + tempFy;

    }

    eFx[0] = -1*eFxTemp * coeffE * nI0;    // since the e and ion with oppsite charge state
    eFy[0] = -1*eFyTemp * coeffE * nI0;    // since the e and ion with oppsite charge state
    
//    cout<<latticeInterActionPoint.ionAccumuFx[k][0]<<endl;
//    cout<<latticeInterActionPoint.ionAccumuFy[k][0]<<endl;
//    cout<<posx<<"   "<<posy<<"   "<<eFx[0]<<" "<<eFy[0]<<"   "<<__LINE__<<"  "<< "weak strong force of electron"<<endl;
    
    //    getchar();
    
    omegaE = coeffE * nI0 * pow(CLight,2) /rmsRy/(rmsRx+rmsRy);
    omegaE = sqrt(omegaE);
    
    omegaI = coeffI * CLight /rmsRy/(rmsRx+rmsRy);
    omegaI = sqrt(omegaI);
    
//    cout<<"ssss "<<omegaE/Omegas<<"  "<<omegaI/Omegas<<endl;
//    getchar();





}



void Bunch::BassettiErskine(double posx,double posy,double rmsRxTemp, double rmsRyTemp,double &tempFx,double &tempFy)
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
            coeffBE    = sqrt(PI)/sigma;


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
            coeffBE    = sqrt(PI)/sigma;


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










//    tempPosix   = abs(posx);
//    tempPosiy   = abs(posy);


//    if(pow(rmsRxTemp,2)<pow(rmsRyTemp,2))
//    {
//        temp        = rmsRyTemp;
//        rmsRyTemp   = rmsRxTemp;
//        rmsRxTemp   = temp;
////        tempPosix   = abs(posy);
////        tempPosiy   = abs(posx);

//    }


//    if(abs(tempPosix/rmsRxTemp) + abs(tempPosiy/rmsRyTemp)<1.0E-6)
//    {
//        tempFx=0.E0;
//        tempFy=0.E0;

//    }
//    else
//    {
//        sigma    = sqrt(2*pow(rmsRxTemp,2)-2*pow(rmsRyTemp,2));
//        ryOverRx = rmsRyTemp/rmsRxTemp;
//    


//        double coeffBE;
//        coeffBE    = sqrt(PI)/sigma;


//        z1           =  {tempPosix/sigma,tempPosiy/sigma};
//        w1           =  Faddeeva::w(z1);

//        z2           =  {ryOverRx*tempPosix/sigma,tempPosiy/sigma/ryOverRx};
//        w2           =  Faddeeva::w(z2);

//        z3           =  - pow(tempPosix/rmsRxTemp,2)/2 - pow(tempPosiy/rmsRyTemp,2)/2;
//        w3           =  - exp(z3); 
//        
//        w0           =  coeffBE * (w1 + w3 * w2 ); 
//        
//        tempFx       = w0.imag();
//        tempFy       = w0.real();
//        

//        if(posx<=0)
//        {
//            tempFx  = -tempFx;
//        }

//        if(posy<=0)
//        {
//            tempFy = -tempFy; 
//        }

//    }


//    if(pow(rmsRxTemp,2)<pow(rmsRyTemp,2))
//    {
//        temp     = tempFy;
//        tempFy   = tempFx;
//        tempFx   = temp;
//    }

        tempFx = -tempFx;
        tempFy = -tempFy;


//    cout<<"Bunch:: BassettiErskine  "<<tempFx<<" "<<tempFy<<"   "<<w0<<"    "<<z1<<"    "<<z2<<"    "<<z3<<"    "
//        <<w1<<"    "<<w2<<"    "<<w3<<"    "<<tempPosix<<"    "<<tempPosiy<<" "<<sigma <<"  "<<rmsRx<<" "<<rmsRy<<endl;

////////////////// this part is used to check the calculation of the complex error functions//////////
//    for(int i=0;i<=100;i++)
//    {
//        for(int j=0;j<=20;j++)
//        {
//            posx        = -10*rmsRx + 10*rmsRx/50*i;
//            posy        = -10*rmsRy + 10*rmsRy/10*j;
//            
//            tempPosix = abs(posx);
//            tempPosiy = abs(posy);
//            
//            z1           =  {tempPosix/sigma,tempPosiy/sigma};
//            w1           =  Faddeeva::w(z1);
//            z2           =  {ryOverRx*tempPosix/sigma,tempPosiy/sigma/ryOverRx};
//            w2           =  Faddeeva::w(z2);

//            z3           =  - pow(tempPosix/rmsRx,2)/2 - pow(tempPosiy/rmsRy,2)/2;
//            w3           =  - exp(z3); 
//            
//            w0           =  coeff * (w1 + w3 * w2 );


//            tempFx      =   w0.imag();          
//            tempFy      =   w0.real();

//            if(posx<=0.E0)
//            {
//                tempFx  = -tempFx;
//            }

//            if(posy<0.E0)
//            {
//                tempFy = -tempFy; 
//            }
//       
//          
//            fout<<i<<"     "
//                <<j<<"     "
//                <<posx<<"     "
//                <<posy<<"     "
//                <<tempFx<<"     "
//                <<tempFy<<"     "
//                <<endl;
//       
//        }
//    }
//  cout<<"over1111ssssssssss11111"<<endl;
//  getchar();
////////////////////////////////////test is over ///////////////////////


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




void Bunch::BunchTransferDueToLattice(LatticeInterActionPoint &latticeInterActionPoint, int k)
{

    double xtemp=0.E0;
    double ytemp=0.E0;
    double ztemp=0.E0;
    double xPtemp=0.E0;
    double yPtemp=0.E0;
	double zPtemp=0.E0;


    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        xtemp  = latticeInterActionPoint.xTransferMatrix[k][0] * ePositionX[i]
               + latticeInterActionPoint.xTransferMatrix[k][1] * eMomentumX[i];
        
        xPtemp = latticeInterActionPoint.xTransferMatrix[k][2] * ePositionX[i]
               + latticeInterActionPoint.xTransferMatrix[k][3] * eMomentumX[i];
                       
        ytemp  = latticeInterActionPoint.yTransferMatrix[k][0] * ePositionY[i]
               + latticeInterActionPoint.yTransferMatrix[k][1] * eMomentumY[i];
        
        yPtemp = latticeInterActionPoint.yTransferMatrix[k][2] * ePositionY[i]
               + latticeInterActionPoint.yTransferMatrix[k][3] * eMomentumY[i];
               
        ztemp  = latticeInterActionPoint.zTransferMatrix[k][0] * ePositionZ[i]
               + latticeInterActionPoint.zTransferMatrix[k][1] * eMomentumZ[i];
               
        zPtemp = latticeInterActionPoint.zTransferMatrix[k][2] * ePositionZ[i]
			   + latticeInterActionPoint.zTransferMatrix[k][3] * eMomentumZ[i];
                       
                       
                      
           
           
        ePositionX[i] = xtemp;
        ePositionY[i] = ytemp;
        eMomentumX[i] = xPtemp;
        eMomentumY[i] = yPtemp;
        eMomentumZ[i] = xPtemp;
        eMomentumZ[i] = yPtemp;
        
        if(pow(ePositionX[i],2) + pow(ePositionY[i],2)>pipeAperatureR)
        {
            eSurive[i] = 0;
        }   
    }


}


void Bunch::BunchSynRadDamping(vector<double> &synchRadDampTime, LatticeInterActionPoint &latticeInterActionPoint)
{

    int k=0;

    //	 set the J2_sympeletic matrix
    gsl_matrix * sympleMarixJ = gsl_matrix_alloc (2, 2); 
    gsl_matrix_set_zero(sympleMarixJ);
    gsl_matrix_set (sympleMarixJ, 0, 1,1);
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
    //	
    //	// get H31 sub_matrix;
    gsl_matrix_transpose(dispMatrixHX);

    gsl_matrix_mul(dispMatrixHX,sympleMarixJ,tempMatrix1);
    gsl_matrix_mul(sympleMarixJ,tempMatrix1,tempMatrix2);
    gsl_matrix_scale (tempMatrix2, -1);

    gsl_matrix_memcpy (H31, tempMatrix2);

    gsl_matrix_transpose(dispMatrixHX);

    gsl_matrix_set(dispMatrixH,4,0,gsl_matrix_get(H31,0,0));
    gsl_matrix_set(dispMatrixH,5,0,gsl_matrix_get(H31,1,0));
    gsl_matrix_set(dispMatrixH,4,1,gsl_matrix_get(H31,0,1));
    gsl_matrix_set(dispMatrixH,5,1,gsl_matrix_get(H31,1,1));


    //	// get H32 sub_matrix;
    gsl_matrix_transpose(dispMatrixHY);

    gsl_matrix_mul(dispMatrixHY,sympleMarixJ,tempMatrix1);
    gsl_matrix_mul(sympleMarixJ,tempMatrix1,tempMatrix2);
    gsl_matrix_scale (tempMatrix2, -1);

    gsl_matrix_memcpy (H32, tempMatrix2);	
    gsl_matrix_transpose(dispMatrixHY);

    //	// set H Marix
    gsl_matrix_set(dispMatrixH,0,5,gsl_matrix_get(dispMatrixHX,0,1));
    gsl_matrix_set(dispMatrixH,1,5,gsl_matrix_get(dispMatrixHX,1,1));
    gsl_matrix_set(dispMatrixH,2,5,gsl_matrix_get(dispMatrixHY,0,1));
    gsl_matrix_set(dispMatrixH,3,5,gsl_matrix_get(dispMatrixHY,1,1));

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

    gsl_matrix * R2  = gsl_matrix_alloc (2, 2);   //coupling matrix Eq.(6)
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

    emittanceZ = rmsBunchLength * rmsEnergySpread;

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
            randR[j] = dx(gen);	
            if(j==4)
            {
                randR[j] = 0.E0;
            }
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

        if(pow(ePositionX[i],2) + pow(ePositionY[i],2)>pipeAperatureR)
        {
            eSurive[i] = 0;
        }     

    }
    



    gsl_matrix_free (dispMatrixH);
    gsl_matrix_free (dispMatrixHX);
    gsl_matrix_free (dispMatrixHY);

    gsl_matrix_free (tengMatrixR);
    gsl_matrix_free (R2);
    gsl_matrix_free (R21);
    gsl_matrix_free (R12);

    gsl_matrix_free (twissMatrixB);
    gsl_matrix_free (sympleMarixJ);

    gsl_matrix_free (tempMatrix1);
    gsl_matrix_free (tempMatrix2);
    gsl_matrix_free (H31);
    gsl_matrix_free (H32);

    gsl_matrix_free (cordTransfer);
    gsl_matrix_free (invCordTransfer);
    gsl_matrix_free (cordTransferTemp);
    gsl_matrix_free (sympleMarixJ6);

    gsl_matrix_free (vecX);
    gsl_matrix_free (vecNX);
  
}




void Bunch::InitialBeamCurDenZProf()
{
    
    std::vector<double>::iterator zMax = std::max_element(std::begin(ePositionZ), std::end(ePositionZ));
    std::cout << "Max element is " << *zMax<< " at position " << std::distance(std::begin(ePositionZ), zMax) << std::endl;
    
    
    std::vector<double>::iterator zMin = std::min_element(std::begin(ePositionZ), std::end(ePositionZ));
    std::cout << "Max element is " << *zMin<< " at position " << std::distance(std::begin(ePositionZ), zMin) << std::endl;

    double rangeBeamSizeZ;
    if(abs(*zMin)>abs(*zMax))
    {
        rangeBeamSizeZ = abs(*zMin);
    }
    else
    {
        rangeBeamSizeZ = abs(*zMax);
    }

    double binZSize;
    
    binZSize = 2 * rangeBeamSizeZ / bunchBinNumberZ;
    
    for(int i=0;i<beamCurDenZProf.size();i++)
    {
        beamCurDenZProf[i]=0.E0 ;
    }


    ofstream fout("beamCurDenZProf.dat");


    double tempDeltaz;
    double poszTemp;
    int poszTempIndex;
    for(int i=0;i<macroEleNumPerBunch;i++)
    {
        poszTemp         = abs(rangeBeamSizeZ-ePositionZ[i]) / binZSize;
        poszTempIndex    = floor(poszTemp);
//        tempDeltaz       = poszTemp - poszTempIndex;

        beamCurDenZProf[poszTempIndex  ]  = beamCurDenZProf[poszTempIndex  ]+1-tempDeltaz;
        beamCurDenZProf[poszTempIndex+1]  = beamCurDenZProf[poszTempIndex+1]+  tempDeltaz;

//        beamCurDenZProf[poszTempIndex]  =   beamCurDenZProf[poszTempIndex] +1;
    }

    for(int i=0;i<beamCurDenZProf.size();i++)
    {
        fout<< i<<" "<<beamCurDenZProf[i]<<endl ;
    }
    
    cout<<"TTTTTT"<<endl;
    getchar();
    fout.close();
    
}

void Bunch::InitialBeamCurDenTProf()
{

}





