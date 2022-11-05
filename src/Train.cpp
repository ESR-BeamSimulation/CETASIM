//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             


#include <vector>
#include <complex>
#include "Train.h"
#include <iostream>
#include<fstream>
#include<iomanip>
	
#include <numeric>
#include <cmath>

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

    int harmonics = inputParameter.ringParBasic->harmonics;
    trainNumber   = inputParameter.ringFillPatt->trainNumber; 
    totBunchNum   = inputParameter.ringFillPatt->totBunchNumber;  

    bunchNumberPerTrain.resize(trainNumber);
    trainGaps.resize(trainNumber);
    bunchGaps.resize(trainNumber);
    trainStart.resize(trainNumber);
    

    for (int i=0; i<trainNumber; i++)
    {
        bunchNumberPerTrain[i] = inputParameter.ringFillPatt->bunchNumberPerTrain[i];
        trainGaps[i]           = inputParameter.ringFillPatt->trainGaps[i];
        bunchGaps[i]           = inputParameter.ringFillPatt->bunchGaps[i];            
    }
    
 
    int bunchNumberCheck = 0;
    int harmonicCheck=0;

    cout<<"---------------------- bunch filling pattern--------------------------"<<endl;

    for(int i=0;i<trainNumber;i++)
    {
        bunchNumberCheck +=  bunchNumberPerTrain[i]; 
        harmonicCheck    +=  bunchNumberPerTrain[i] * (bunchGaps[i]+1) + trainGaps[i];
        
        cout<<"train "<<i<<": "
            <<setw(5)<<bunchNumberPerTrain[i] <<" bunches   "
            <<setw(5)<<(bunchGaps[i]+1) * bunchNumberPerTrain[i] <<" buckets   "
            <<setw(5)<<trainGaps[i]<<" buckets-train-gaps"
            <<endl;        
    }
    cout<<"-----------------------------------------------------------------------"<<endl;
    
        
    if(bunchNumberCheck!=totBunchNum || harmonicCheck != harmonics)
    {
        cerr<<"mistakes in settings for bunch filling pattern, sum of gaps + totBunchNumber has to be harmonics"<<endl;
        exit(0);
    }


    for(int i=0;i<trainNumber;i++)
    {
        trainStart[i] = 0;        
        if(i>0)
        {
            for(int j=0;j<i;j++)
            {
                trainStart[i] = trainStart[i] +  bunchNumberPerTrain[j] * (bunchGaps[i]+1)  + trainGaps[j];
            }
        }
    }
}





void Train::InitialDesy(ReadInputSettings &inputParameter)
{
	trainNumber = inputParameter.ringFillPatt->trainNumber;

    bunchNumberPerTrain.resize(trainNumber);
    trainGaps.resize(trainNumber);
    trainStart.resize(trainNumber);
    
    int harmonics = inputParameter.ringParBasic->harmonics;

	totBunchNum = inputParameter.ringFillPatt->totBunchNumber;
	
	
	for(int trainIndex=0;trainIndex<trainNumber;trainIndex++)
	{
		bunchNumberPerTrain[trainIndex]=1;
		trainGaps[trainIndex]=1;
		if((trainIndex+1)%20==0)
		{
			trainGaps[trainIndex]=9;
		}
	}
	
	


    int bunchNumberCheck = 0;
    int gapsNumberCheck = 0;
    for(int i=0;i<trainNumber;i++)
    {
         bunchNumberCheck +=  bunchNumberPerTrain[i];
         gapsNumberCheck  +=  trainGaps[i];
    }


   
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


void Train::InitialUSSR310(ReadInputSettings &inputParameter)
{
	
	trainNumber = inputParameter.ringFillPatt->trainNumber;

    bunchNumberPerTrain.resize(trainNumber);
    trainGaps.resize(trainNumber);
    trainStart.resize(trainNumber);
    
    int harmonics = inputParameter.ringParBasic->harmonics;

	totBunchNum = inputParameter.ringFillPatt->totBunchNumber;
	
	// 310*4 = 1240  
	for(int trainIndex=0;trainIndex<trainNumber;trainIndex++)
	{
		bunchNumberPerTrain[trainIndex]=1;			
		trainGaps[trainIndex]=3;		
	}
	
	//trainGaps[trainNumber-1]=343;


    int bunchNumberCheck = 0;
    int gapsNumberCheck = 0;
    for(int i=0;i<trainNumber;i++)
    {
         bunchNumberCheck +=  bunchNumberPerTrain[i];
         gapsNumberCheck  +=  trainGaps[i];
    }



   
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











