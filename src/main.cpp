//*************************************************************************
//Copyright (c) 2020 IHEP                                                  
//Copyright (c) 2021 DESY                                                  
//This program is free software; you can redistribute it and/or modify     
//it under the terms of the GNU General Public License                     
//Author: chao li, li.chao@desy.de                                         
//*************************************************************************
#pragma once                                                             
#include "stdio.h"
#include "string.h"
//#include "mpi.h"
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <omp.h>
#include "Global.h"
//#include "SDDS.h"    // refer to Makefile, only ensure compliable
#include"ReadInputSettings.h"
#include "LatticeInterActionPoint.h"
#include "Train.h"
#include "CavityResonator.h"
#include "SPBeam.h"
#include "MPBeam.h"


using namespace std;
using std::vector;
using std::complex;

int main(int argc,char *argv[])
{

//    MPI_Init(&argc,&argv);
//    MPI_Status status;

//    MPI_Comm_size(MPI_COMM_WORLD,&numProcess);
//    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

    /*
    int num=omp_get_max_threads();
    #pragma omp parallel private(num)
    {
        int num = omp_get_thread_num();
        printf("%d\n",num);
    }
    getchar();
    */
        
    ReadInputSettings inputParameter;
    inputParameter.ParamRead(argc, argv);
 
    LatticeInterActionPoint latticeInterActionPoint;
    latticeInterActionPoint.Initial(inputParameter);

    
    
    Train train;
    train.Initial(inputParameter);

    CavityResonator cavityResonator;
    cavityResonator.Initial(inputParameter);

    

    if(inputParameter.ringRun->calSetting==1 && inputParameter.ringBunchPara->macroEleNumPerBunch==1)   // bunch is rigid represneted by only single particle...
    {        
        SPBeam spbeam;
        spbeam.Initial(train,latticeInterActionPoint,inputParameter);
        spbeam.InitialcavityResonator(inputParameter,cavityResonator); 
        spbeam.Run(train,latticeInterActionPoint,inputParameter,cavityResonator);   
    
    }
    else if ((inputParameter.ringRun->calSetting==2 && inputParameter.ringBunchPara->macroEleNumPerBunch!=1))                      
    {
        MPBeam mpbeam;
        mpbeam.Initial(train,latticeInterActionPoint,inputParameter);
        mpbeam.InitialcavityResonator(inputParameter,cavityResonator); 
        mpbeam.Run(train,latticeInterActionPoint,inputParameter,cavityResonator);      
    }
    else
    {
        cerr<<"wrong settings about SP and MP tracking...  "<<endl;
        exit(0);
    }
    


//    MPI_Finalize();
    return 0;
}
