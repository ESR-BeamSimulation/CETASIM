#pragma once

#include"ReadInputSettings.h"

#include<string>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include"Global.h"
#include <algorithm>
using namespace std;
using std::vector;

ReadInputSettings::ReadInputSettings()
{

}

int ReadInputSettings::ParamRead()
{


    ifstream fin("input.dat");
  
 
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"------------ the basic paprameter of the ring- --------"<<endl;
 
 
    vector<string> strVec;
    string         str;
    
    if (! fin.is_open())
    {
        cerr<< "Error opening file input.dat"; exit (1);
    }
    while (!fin.eof())
    {
        getline(fin,str);
        
        if(str.length()==0 || str.find("=")==-1 )  continue;
        StringVecSplit(str, strVec);
        
        if(strVec[0]=="circring")
        {
          circRing = stod(strVec[1]);
          cout<<"the circular of the ring is                    :"<< stod(strVec[1])<<"[m]"<<endl;
        }
        
        if(strVec[0]=="workqx")
        {
          workQx = stod(strVec[1]);
          cout<<"the working point Qx                           :"<< stod(strVec[1])<<endl;
        }

        if(strVec[0]=="workqy")
        {
          workQy = stod(strVec[1]);
          cout<<"the working point Qy                           :"<< stod(strVec[1])<<endl;
        }
        
        if(strVec[0]=="workqz")
        {
          workQz = stod(strVec[1]);
          cout<<"the working point Qz                           :"<< stod(strVec[1])<<endl;
        }
                
        if(strVec[0]=="rfbasefrequency")
        {
          rfBaseFrequency = stod(strVec[1]);
          cout<<"the baisc RF frequency is                      :"<< stod(strVec[1])<<"[Hz]"<<endl;
        }        
 
    
        if(strVec[0]=="corsssectionei")
        {
          corssSectionEI = stod(strVec[1]);
          cout<<"the corssseciton of ion-e interaction          :"<< stod(strVec[1])<<endl;
        }

        if(strVec[0]=="pipeaperaturex")
        {
          pipeAperatureX = stod(strVec[1]);
          cout<<"half of pipe aperture in  x direction (tangl boundary)     :"<< stod(strVec[1])<<"[m]"<<endl;
        }
        
        if(strVec[0]=="pipeaperaturey")
        {
          workQz = stod(strVec[1]);
          cout<<"half of pipe aperture in  y direction  (tangl boundary)    :"<< stod(strVec[1])<<"[m]"<<endl;
        }
                
        if(strVec[0]=="pipeaperaturer")
        {
          pipeAperatureR = stod(strVec[1]);
          cout<<"the rauius of the pipe  (round boundary)                   :"<< stod(strVec[1])<<"[m]"<<endl;
        }   
        if(strVec[0]=="ionmaxnumber")
        {
          ionMaxNumber = stod(strVec[1]);
          cout<<"Max number of ions can be stored                           :"<< stod(strVec[1])<<endl;
        }   
        
        if(strVec[0]=="numberofionbeaminterpoint")
        {
          numberofIonBeamInterPoint = stod(strVec[1]);
          cout<<"the number of ion beam interaction point is                 :"<< stod(strVec[1])<<endl;
        } 
        if(strVec[0]=="electronbeamenergy")
        {
          electronBeamEnergy = stod(strVec[1]);
          cout<<"the energy of the electron beam is                          :"<< stod(strVec[1])<<endl;
        }
        
        if(strVec[0]=="totbunchnumber")
        {
          totBunchNumber = stod(strVec[1]);
          cout<<"the total  electron number of beam bunches is               :"<< stod(strVec[1])<<endl;
        }
        
        if(strVec[0]=="trainnumber")
        {
          trainNumber = stod(strVec[1]);
          cout<<"the number of beam trains in the ring  is                   :"<< stod(strVec[1])<<endl;
        }
        if(strVec[0]=="printinterval")
        {
          printInterval = stod(strVec[1]);
          cout<<"The inteval of truns for data printing                      :"<< stod(strVec[1])<<endl;
        }
        if(strVec[0]=="synchraddamptimex")
        {
          synchRadDampTime.push_back( stod(strVec[1]));
          cout<<"synchraddamptimex in x direction (turns)                     :"<<stod(strVec[1])  <<endl;
        }
        if(strVec[0]=="synchraddamptimey")
        {
          synchRadDampTime.push_back( stod(strVec[1]));
          cout<<"synchraddamptimex in y direction (turns)                     :" <<  stod(strVec[1]) <<endl;
        }
        if(strVec[0]=="synchraddamptimez")
        {
          synchRadDampTime.push_back( stod(strVec[1]));
          cout<<"synchraddamptimez in z direction (turns)                     :"<< stod(strVec[1]) <<endl;
        }
        if(strVec[0]=="macroelenumperbunch")
        {
          macroEleNumPerBunch = stod(strVec[1]);
          cout<<"macro Electron beam Number  Per Bunch                        :"<< stod(strVec[1]) <<endl;
        }
        if(strVec[0]=="bunchbinnumberz")
        {
          bunchBinNumberZ = stod(strVec[1]);
          cout<<"Longitudinal impedance binsize for                           :"<< stod(strVec[1]) <<endl;
        }
        if(strVec[0]=="current")
        {
          current = stod(strVec[1]);
          cout<<"beam current                                                 :"<< stod(strVec[1])<<"[A]" <<endl;
        }
        if(strVec[0]=="emittancex")
        {
          emittanceX = stod(strVec[1]);
          cout<<"beam emittance in x direction                                :"<< stod(strVec[1])<<"[m.rad]"  <<endl;
        }
        if(strVec[0]=="kappa")
        {
          kappa = stod(strVec[1]);
          cout<<"beam emittance coupling factor                                :"<< stod(strVec[1]) <<endl;
        }
        
        if(strVec[0]=="rmsbunchlength")
        {
          rmsBunchLength = stod(strVec[1]);
          cout<<"longitudinal rms bunch length                                 :"<< stod(strVec[1])<<"[m]" <<endl;
        }
        
        if(strVec[0]=="rmsenergyspread")
        {
          rmsEnergySpread = stod(strVec[1]);
          cout<<"longitudinal rms bunch energy spread                                 :"<< stod(strVec[1]) <<endl;
        }
        
        if(strVec[0]=="distributiontype")
        {
          distributionType = stod(strVec[1]);
          cout<<"beam distribution type in transverse                                 :"<< stod(strVec[1]) <<endl;
        }
        
        if(strVec[0]=="nturns")
        {
          nTurns = stod(strVec[1]);
          cout<<"the total calculation turns is                                  :"<< stod(strVec[1]) <<endl;
        }
        
        if(strVec[0]=="calsettings")
        {
          calSettings = stod(strVec[1]);
          cout<<"calculation settins is                                          :"<< stod(strVec[1]) <<endl;
        }
        if(strVec[0]=="synraddampingflag")
        {
          synRadDampingFlag = stod(strVec[1]);
          cout<<"to the    SynRadDamping ?                                      :"<< stod(strVec[1]) <<endl;
        }
        if(strVec[0]=="intevalofturnsiondataprint")
        {
          intevalofTurnsIonDataPrint = stod(strVec[1]);
          cout<<"to the    SynRadDamping ?                                      :"<< stod(strVec[1]) <<endl;
        }
        if(strVec[0]=="firbunchbybunchfeedbackflag")
        {
          fIRBunchByBunchFeedbackFlag = stod(strVec[1]);
          cout<<"WITH FIR bunch by bunch feedback ?                              :"<< stod(strVec[1]) <<endl;
        }
        if(strVec[0]=="macroionnumbergeneratedperip")
        {   
            macroIonNumberGeneratedPerIP = stod(strVec[1]);
            cout<<"macro Ion generated per interaction                             :"<< stod(strVec[1]) <<endl;
        }
        if(strVec[0]=="initialdisdx")
        {   
            initialDisDx = stod(strVec[1]);
            cout<<"Max initial bunch displacement   error   Displacement x          :"<< stod(strVec[1]) <<endl;
        }
       if(strVec[0]=="initialdisdy")
        {   
            initialDisDy = stod(strVec[1]);
            cout<<"Max initial bunch displacement   error   Displacement y          :"<< stod(strVec[1]) <<endl;
        }
    }
    
    emittanceY = emittanceX * kappa;
    cout<<"beam emittance in y direction                                         :"<< emittanceY<<"[m.rad]"   <<endl;
    
    cout<<"-------------------------------------------------------"<<endl;

    fin.close();
    
    
    double rGamma;
    double rBeta;

    rGamma = electronBeamEnergy/ElectronMassEV;
    rBeta  = sqrt(1.E0-1.E0/rGamma);
    
    t0     = circRing/rBeta/CLight; 
    omegas = 2*PI/t0;
    harmonics = round(2*PI*rfBaseFrequency/(omegas));

    return 0;
}

