#include "LatticeInterActionPoint.h"
#include "Global.h"
#include <vector>
#include <stdlib.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <time.h>
#include <random>
#include<string>



LatticeInterActionPoint::LatticeInterActionPoint()
{

}





LatticeInterActionPoint::~LatticeInterActionPoint()
{

}

void LatticeInterActionPoint::Initial(ReadInputSettings &inputParameter)
{

    numberOfInteraction = inputParameter.numberofIonBeamInterPoint;                //HEPS lattice 24 super period
    ionMaxNumberOneInterPoint = inputParameter.ionMaxNumber;
    
    
    circRing =  inputParameter.circRing;
    harmonics = inputParameter.harmonics;
 
    pipeAperatureR = inputParameter.pipeAperatureR;
    pipeAperatureX = inputParameter.pipeAperatureX;
    pipeAperatureY = inputParameter.pipeAperatureY;
    
    

    twissAlphaX.resize(numberOfInteraction);
    twissBetaX.resize(numberOfInteraction);
    twissAlphaY.resize(numberOfInteraction);
    twissBetaY.resize(numberOfInteraction);
    twissAlphaZ.resize(numberOfInteraction);
    twissBetaZ.resize(numberOfInteraction);
    
    twissDispX.resize(numberOfInteraction);				
    twissDispPX.resize(numberOfInteraction);;          
    twissDispY.resize(numberOfInteraction);;				
    twissDispPY.resize(numberOfInteraction);;             

//    interActionLocations.resize(numberOfInteraction);
    vacuumPressure.resize(numberOfInteraction);
    temperature.resize(numberOfInteraction);
    ionLineDensity.resize(numberOfInteraction);
    interactionLength.resize(numberOfInteraction);
    xPhaseAdv.resize(numberOfInteraction);
    yPhaseAdv.resize(numberOfInteraction);
    zPhaseAdv.resize(numberOfInteraction);
    
    macroIonNumber.resize(numberOfInteraction);
    ionNumber.resize(numberOfInteraction);
    macroIonCharge.resize(numberOfInteraction);
    
    ionPositionX.resize(numberOfInteraction);
    ionPositionY.resize(numberOfInteraction);
    ionVelocityX.resize(numberOfInteraction);
    ionVelocityY.resize(numberOfInteraction);
    
    ionAccumuFx.resize(numberOfInteraction);
    ionAccumuFy.resize(numberOfInteraction);
    
    ionAccumuNumber.resize(numberOfInteraction);
    ionAccumuPositionX.resize(numberOfInteraction);
    ionAccumuPositionY.resize(numberOfInteraction);
    ionAccumuVelocityX.resize(numberOfInteraction);
    ionAccumuVelocityY.resize(numberOfInteraction);
    
    ionAccumuAverX.resize(numberOfInteraction);
    ionAccumuAverY.resize(numberOfInteraction);
    ionAccumuRMSX.resize(numberOfInteraction);
    ionAccumuRMSY.resize(numberOfInteraction); 
    
    
    
    
    xTransferMatrix.resize(numberOfInteraction);
    yTransferMatrix.resize(numberOfInteraction);
    zTransferMatrix.resize(numberOfInteraction);
    
    transferMatrix.resize(numberOfInteraction);
    
  

    for(int i=0;i<numberOfInteraction;i++)
    {
        xTransferMatrix[i].resize(4);
        yTransferMatrix[i].resize(4);
        zTransferMatrix[i].resize(4);
    }

}



void LatticeInterActionPoint::InitialLattice(ReadInputSettings &inputParameter)
{


    ifstream fin("InterPointParameter.dat");
    if (! fin.is_open())
    {
        cerr<< "Error opening file InterPointParameter.dat"; exit (1);
    }
    
    
    string str;
    vector<string> strVec;
    
    getline(fin,str);
    bool flag;
    
    int index = 0;
    
    while (!fin.eof())
    {
        getline(fin,str);

        StringSplit(str, strVec);
        
        string stringTest;
        for(int i=0;i<str.size();i++)
        {
            stringTest.push_back(' ');
        }
        
        if(stringTest == str )  continue;
        
        for(int i=0;i<numberOfInteraction;i++)
        {
            twissAlphaX[i] = stod(strVec[0]);
            twissAlphaY[i] = stod(strVec[1]);
            twissAlphaZ[i] = stod(strVec[2]);
            twissBetaX[i]  = stod(strVec[3]);
            twissBetaY[i]  = stod(strVec[4]); 
            twissBetaZ[i]  = stod(strVec[5]);
            
            twissDispX[i]  = stod(strVec[6]);
            twissDispPX[i] = stod(strVec[7]);
            twissDispY[i]  = stod(strVec[8]);
            twissDispPY[i] = stod(strVec[9]);
            
            temperature[i] = stod(strVec[10]) ;
            vacuumPressure[i] =  stod(strVec[11])* 1.0E-9 * 133.3224;;
            
        }

        index++;
    }



    if((index)!=inputParameter.numberofIonBeamInterPoint)
    {
        cerr<<"data of numberofIonBeamInterPoint in input.dat and InterPointParameter.dat files does not match";
        exit(0);
    }


    fin.close();

    double betaX1;
    double betaX2;
    double alphaX1;
    double alphaX2;
    
    double betaY1;
    double betaY2;
    double alphaY1;
    double alphaY2;
  
	double betaZ1;
    double betaZ2;
    double alphaZ1;
    double alphaZ2;
  
    
    
    double phaseAdvanceX;
    double phaseAdvanceY;
    double phaseAdvanceZ;
    

    
    double workQx = inputParameter.workQx;
    double workQy = inputParameter.workQy;
    double workQz = inputParameter.workQz;

    
    

    for(int i=0;i<numberOfInteraction;i++)
    {


        xPhaseAdv[i]         = 2*PI*workQx/numberOfInteraction; // [rad]
        yPhaseAdv[i]         = 2*PI*workQy/numberOfInteraction; // [rad]    
        zPhaseAdv[i]         = 2*PI*workQz/numberOfInteraction; // [rad]   
        interactionLength[i] = circRing/numberOfInteraction;    // [m]


        alphaX1 = twissAlphaX[i];
        betaX1  = twissBetaX[i];
        alphaY1 = twissAlphaY[i];
        betaY1  = twissBetaY[i];
        alphaZ1 = twissAlphaZ[i];
        betaZ1	= twissBetaZ[i];
        
        
        phaseAdvanceX = xPhaseAdv[i];
        phaseAdvanceY = yPhaseAdv[i];
        phaseAdvanceZ = zPhaseAdv[i];
        
        if(i<numberOfInteraction-1)
        {
            alphaX2 = twissAlphaX[i+1];
            betaX2  = twissBetaX[i+1];
            alphaY2 = twissAlphaY[i+1];
            betaY2  = twissBetaY[i+1];
            alphaZ2 = twissAlphaZ[i+1];
            betaZ2  = twissBetaZ[i+1];
        }
        else
        {
            alphaX2 = twissAlphaX[0];
            betaX2  = twissBetaX[0];
            alphaY2 = twissAlphaY[0];
            betaY2  = twissBetaY[0];
            alphaZ2 = twissAlphaZ[0];
            betaZ2  = twissBetaZ[0];
        }
               
        xTransferMatrix[i][0]  = sqrt(betaX2 / betaX1) * (cos(phaseAdvanceX) + alphaX1 * sin(phaseAdvanceX));
        xTransferMatrix[i][1]  = sqrt(betaX2 * betaX1) * sin(phaseAdvanceX);
        xTransferMatrix[i][2]  = -(1 + alphaX1 * alphaX2)/sqrt(betaX1 * betaX2) * sin(phaseAdvanceX) 
                                 +(    alphaX1 - alphaX2)/sqrt(betaX1 * betaX2) * cos(phaseAdvanceX);
        xTransferMatrix[i][3]  = sqrt(betaX1 / betaX2) * (cos(phaseAdvanceX) - alphaX2 * sin(phaseAdvanceX));

        yTransferMatrix[i][0]  = sqrt(betaY2 / betaY1) * (cos(phaseAdvanceY) + alphaY1 * sin(phaseAdvanceY));
        yTransferMatrix[i][1]  = sqrt(betaY2 * betaY1) * sin(phaseAdvanceY);
        yTransferMatrix[i][2]  = -(1 + alphaY1 * alphaY2)/sqrt(betaY1 * betaY2) * sin(phaseAdvanceY) 
                                 +(    alphaY1 - alphaY2)/sqrt(betaY1 * betaY2) * cos(phaseAdvanceY);
        yTransferMatrix[i][3]  = sqrt(betaY1 / betaY2) * (cos(phaseAdvanceY) - alphaY2 * sin(phaseAdvanceY));
		
	
        zTransferMatrix[i][0]  = sqrt(betaZ2 / betaZ1) * (cos(phaseAdvanceZ) + alphaZ1 * sin(phaseAdvanceZ));
        zTransferMatrix[i][1]  = sqrt(betaZ2 * betaZ1) * sin(phaseAdvanceZ);
        zTransferMatrix[i][2]  = -(1 + alphaZ1 * alphaZ2)/sqrt(betaZ1 * betaZ2) * sin(phaseAdvanceZ) 
                                 +(    alphaZ1 - alphaZ2)/sqrt(betaZ1 * betaZ2) * cos(phaseAdvanceZ);
        zTransferMatrix[i][3]  = sqrt(betaZ1 / betaZ2) * (cos(phaseAdvanceZ) - alphaZ2 * sin(phaseAdvanceZ));	


    }




    for(int k=0; k<numberOfInteraction;k++)
    {
        macroIonNumber[k]=inputParameter.macroIonNumberGeneratedPerIP;          // at each interaction point, one macroIon is generated 
    }



//  at kth interaction point, macroIonNumber[k]  ions are generated
    for(int k=0;k<numberOfInteraction;k++)
    {
        ionPositionX[k].resize(macroIonNumber[k]);
        ionPositionY[k].resize(macroIonNumber[k]);
        ionVelocityX[k].resize(macroIonNumber[k]);
        ionVelocityY[k].resize(macroIonNumber[k]);
        
    }





}


void LatticeInterActionPoint::IonGenerator(double rmsRx, double rmsRy, double xAver,double yAver, int k)
{


    macroIonCharge[k] = ionNumber[k]/macroIonNumber[k];



//    srand(time(0)+randomIndex);

//    std::random_device myseed0{};
//    std::random_device myseed1{};

//    double ix   = double(rand())/RAND_MAX;
//    double iy   = double(rand())/RAND_MAX;
//    
//    mt19937 engine1(myseed0()) ;
//    mt19937 engine2(myseed1()) ;

//    cout<<myseed0()<<endl;
//    cout<<myseed1()<<endl;
//    getchar();

//    std::normal_distribution<double> nx(xAver, rmsRx);
//    std::normal_distribution<double> ny(yAver, rmsRy);



    std::random_device rd{};
    std::mt19937 gen{rd()};
 

    std::normal_distribution<> dx{xAver,rmsRx};
    std::normal_distribution<> dy{yAver,rmsRy};
    
    double tempx;
    double tempy;


    int i=0;
    while(i<macroIonNumber[k])
    {


//        tempx = Gaussrand(rmsRx,xAver,i);
//        tempy = Gaussrand(rmsRy,yAver,i);

        tempx = dx(gen);
        tempy = dy(gen);



        if( pow( (tempx-xAver)/rmsRx, 2)  + pow( (tempy-yAver)/rmsRy, 2) > 9.E0  ) // beam is 3 sigma truncted.
        {
            continue;
        }
        


//        ionPositionX[k][i]=tempx;
//        ionPositionY[k][i]=tempy;
        
        ionPositionX[k][i]=rmsRx;
        ionPositionY[k][i]=rmsRy;
        
        
        ionVelocityX[k][i]=0.E0;
        ionVelocityY[k][i]=0.E0;

        i++;
    }
    

    

}

void LatticeInterActionPoint::IonsUpdate(int k)
{
    vector<int> lossIonIndex;

//****************2019-04-22--original method, ionAccumuNumber[k] increaes montotnously during calculation
//    for(int i=0;i<ionAccumuNumber[k];i++)
//    {
//        if(pow(ionAccumuPositionX[k][i],2)+pow(ionAccumuPositionY[k][i],2)>pow(PipeAperatureR,2))
//        {
//            lossIonIndex.push_back(i);
//        }
//    }


//    for(int i=0;i<macroIonNumber[k];i++)
//    {
//        if(i<lossIonIndex.size())
//        {
//            ionAccumuPositionX[k][lossIonIndex[i]] = ionPositionX[k][i];
//            ionAccumuPositionY[k][lossIonIndex[i]] = ionPositionY[k][i];
//            ionAccumuVelocityX[k][lossIonIndex[i]] = ionVelocityX[k][i];
//            ionAccumuVelocityY[k][lossIonIndex[i]] = ionVelocityY[k][i];
//            ionAccumuFx[k][lossIonIndex[i]]        = 0.E0;
//            ionAccumuFy[k][lossIonIndex[i]]        = 0.E0;
//        }
//        else
//        {
//            ionAccumuPositionX[k].push_back(ionPositionX[k][i]);
//            ionAccumuPositionY[k].push_back(ionPositionY[k][i]);
//            ionAccumuVelocityX[k].push_back(ionVelocityX[k][i]);
//            ionAccumuVelocityY[k].push_back(ionVelocityY[k][i]);
//            ionAccumuFx[k].push_back(0.E0);
//            ionAccumuFy[k].push_back(0.E0);
//        }

//    }
//    
//    ionAccumuNumber[k]  =  ionAccumuPositionX[k].size();
//
//
//************************************** end **********************************************




// ********modified 2019-04-22, the ionAccumuNumber[k] dynamically changed with calculation*******

    for(int i=0;i<macroIonNumber[k];i++)
    {
            ionAccumuPositionX[k].push_back(ionPositionX[k][i]);
            ionAccumuPositionY[k].push_back(ionPositionY[k][i]);
            ionAccumuVelocityX[k].push_back(ionVelocityX[k][i]);
            ionAccumuVelocityY[k].push_back(ionVelocityY[k][i]);
            ionAccumuFx[k].push_back(0.E0);
            ionAccumuFy[k].push_back(0.E0);
    }


    int i=0;
    
    while(i<ionAccumuPositionX[k].size())
    {
        if(pow(ionAccumuPositionX[k][i],2)+pow(ionAccumuPositionY[k][i],2)>pow(pipeAperatureR,2))
        {
            ionAccumuPositionX[k].erase(ionAccumuPositionX[k].begin()+i);
            ionAccumuPositionY[k].erase(ionAccumuPositionY[k].begin()+i);
            ionAccumuVelocityX[k].erase(ionAccumuVelocityX[k].begin()+i);
            ionAccumuVelocityY[k].erase(ionAccumuVelocityY[k].begin()+i);
            ionAccumuFx[k].erase(ionAccumuFx[k].begin()+i);
            ionAccumuFy[k].erase(ionAccumuFy[k].begin()+i);
            i--;
        }
        
        
        i++;
    }
    
    ionAccumuNumber[k]  =  ionAccumuPositionX[k].size(); 
    

//*******************end ****************************************************************

    if(ionAccumuNumber[k]>ionMaxNumberOneInterPoint)
    {
        cerr<<"Ions accumulated are larger than setting limit "<<ionMaxNumberOneInterPoint<<endl;
        exit(0);
    }


}



void LatticeInterActionPoint::IonRMSCal(int k)
{


    double x2Aver=0.E0;
    double y2Aver=0.E0;
    
    ionAccumuAverX[k] =0.E0;
    ionAccumuAverY[k] =0.E0;
    ionAccumuRMSX[k]  =0.E0;
    ionAccumuRMSY[k]  =0.E0;
    
    ionAccumuNumber[k]=ionAccumuPositionX[k].size();

    ionAccumuAverX[k]   =   accumulate(begin(ionAccumuPositionX[k]), end(ionAccumuPositionX[k]), 0.0); 
    ionAccumuAverY[k]   =   accumulate(begin(ionAccumuPositionY[k]), end(ionAccumuPositionY[k]), 0.0); 



    ionAccumuAverX[k]   = ionAccumuAverX[k]   /ionAccumuNumber[k];
    ionAccumuAverY[k]   = ionAccumuAverY[k]   /ionAccumuNumber[k];



    for(int i=0;i<ionAccumuNumber[k];i++)
    {
        x2Aver  +=  pow(ionAccumuPositionX[k][i]-ionAccumuAverX[k] ,2);
        y2Aver  +=  pow(ionAccumuPositionY[k][i]-ionAccumuAverY[k] ,2);
        
    }

    x2Aver  = x2Aver  /ionAccumuNumber[k];
    y2Aver  = y2Aver  /ionAccumuNumber[k];
    
    ionAccumuRMSX[k] = sqrt(x2Aver);
    ionAccumuRMSY[k] = sqrt(y2Aver);

}


void LatticeInterActionPoint::IonTransferDueToBunch(int bunchGap,int k)
{

    
    for(int i=0;i<ionAccumuNumber[k];i++)
    {

        // kick----drift model to update ion velocity and position, time step is the bunchGap
        

        ionAccumuVelocityX[k][i] +=  ionAccumuFx[k][i]  ; 
        ionAccumuVelocityY[k][i] +=  ionAccumuFy[k][i]  ;
        

        ionAccumuPositionX[k][i] +=  ionAccumuVelocityX[k][i] * circRing/harmonics*bunchGap/CLight;
        ionAccumuPositionY[k][i] +=  ionAccumuVelocityY[k][i] * circRing/harmonics*bunchGap/CLight;
        

        
    }


}




