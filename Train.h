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
    int totBunchNum=0;
    vector<int> trainGaps;     //  number of buckets between adjacent trains 
    vector<int> trainStart;    //  start index of each train 
    
    
    vector<double> bunchNumberPerTrain;
    void Initial(ReadInputSettings &inputParameter);

private:

};




#endif
