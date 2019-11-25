

#include <vector>
#include <complex>
#include "Train.h"
#include <iostream>
#include<fstream>

using namespace std;
using std::vector;
using std::complex;



Train::Train()
{

}


Train::~Train()
{
    
}

void Train::Initial(ReadInputSettings &inputParameter)
{

    trainNumber = inputParameter.trainNumber;


    bunchNumberPerTrain.resize(trainNumber);
    trainGaps.resize(trainNumber);
    trainStart.resize(trainNumber);
    

    int harmonics = inputParameter.harmonics;

    

//  need to re-write this part to have a flexibility to cover various problems
    
    ifstream fin("input.dat");
    
    vector<string> strVec;
    string         str;
    
    while (!fin.eof())
    {
        getline(fin,str);
        if(str.length()==0 )  continue;

        StringSplit(str, strVec);
        
        transform(strVec[0].begin(), strVec[0].end(), strVec[0].begin(), ::tolower);
        

        if(strVec[0]=="bunchnumberpertrain")
        {
            cout<<"the number of bunches in trains                            ";
            for(int i=1;i<strVec.size();i++)
            {
                bunchNumberPerTrain[i-1] = stod(strVec[i]);
                cout<<bunchNumberPerTrain[i-1]<<"       ";
            }
            cout<<endl;
        }
        if(strVec[0]=="traingaps")
        {
            cout<<"the gaps betweens bunch trains                             ";
            for(int i=1;i<strVec.size();i++)
            {
                trainGaps[i-1]           = stod(strVec[i]);
                cout<<trainGaps[i-1]<<"     ";
            }
            cout<<endl;
        }
    }

     cout<<"--------------------------------------------------------------------"  <<endl;

    int bunchNumberCheck = 0;
       int gapsNumberCheck = 0;
    for(int i=0;i<trainNumber;i++)
    {
         bunchNumberCheck +=  bunchNumberPerTrain[i];
         gapsNumberCheck  +=  trainGaps[i];
    }

    totBunchNum = inputParameter.totBunchNumber;
    
//    cout<<bunchNumberCheck<<endl;
//    cout<<gapsNumberCheck<<endl;
//    getchar();
    
    if(bunchNumberCheck!=totBunchNum || gapsNumberCheck != harmonics -bunchNumberCheck )
    {
        cerr<<"mistakes in settings for bunch filling pattern "<<endl;
        exit(0);
    }

    for(int i=0;i<trainNumber;i++)
    {
        trainStart[i] = 0;
        
        if(i>0)
        {
            for(int j=0;j<i;j++)
            {
                trainStart[i] = trainStart[i] +  bunchNumberPerTrain[j] + trainGaps[j];
            }
        }
    }

}
