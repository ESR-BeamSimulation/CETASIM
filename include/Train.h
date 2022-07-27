//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#ifndef TRAIN_H
#define TRAIN_H


#include "Global.h"
#include <vector>
#include <complex>
#include "ReadInputSettings.h"


using namespace std;
using std::vector;
using std::complex;


class Train
{


public:
    Train();
    ~Train();
    
    int trainNumber;
    int totBunchNum;
    vector<int> trainGaps;     //  number of buckets between adjacent trains 
    vector<int> bunchGaps;     //  number of buckets between adjacent bunch 
    
    vector<int> trainStart;    //  start index of each train 
    
    
    vector<double> bunchNumberPerTrain;
    void Initial(ReadInputSettings &inputParameter);
	void InitialDesy(ReadInputSettings &inputParameter);
	void InitialUSSR310(ReadInputSettings &inputParameter);

private:

};




#endif
