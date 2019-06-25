

#include <vector>
#include <complex>
#include "Train.h"
#include <iostream>

using namespace std;
using std::vector;
using std::complex;



Train::Train()
{
    trainNumber = 1;
    bunchNumberPerTrain.resize(trainNumber);
    trainGaps.resize(trainNumber);
    trainStart.resize(trainNumber);
}


Train::~Train()
{
    
}

void Train::Initial()
{

//  need to re-write this part to have a flexibility to cover various problems
    for(int i=0;i<trainNumber;i++)
    {
        bunchNumberPerTrain[i] = 680;
        trainGaps[i]           = Harmonics/trainNumber - bunchNumberPerTrain[i];
    }



    for(int i=0;i<trainNumber;i++)
    {
        totBunchNum = totBunchNum + bunchNumberPerTrain[i]; // has to be less than harmonics
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
//        cout<<trainStart[i]<<endl;
    }


}
