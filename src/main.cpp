#include "stdio.h"
#include "string.h"
//#include "mpi.h"

#include <fstream>
#include <stdlib.h> 
#include <iostream>
#include <vector>
#include "Global.h"
#include"ReadInputSettings.h"


#include "LatticeInterActionPoint.h"
#include "Train.h"

#include "Beam.h"
/*
#include "LongImpSingalBunch.h"
*/

using namespace std;
using std::vector;
using std::complex;

int main(int argc,char *argv[])
{

//    MPI_Init(&argc,&argv);
//    MPI_Status status;

//    MPI_Comm_size(MPI_COMM_WORLD,&numProcess);
//    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

    ReadInputSettings inputParameter;    
    inputParameter.ParamRead();
    	
    
    LatticeInterActionPoint latticeInterActionPoint;     	
    latticeInterActionPoint.Initial(inputParameter);
    
    latticeInterActionPoint.InitialLattice(inputParameter);
       
    Train train;	
	train.Initial(inputParameter);    
	//train.InitialDesy(inputParameter);
	//train.InitialUSSR310(inputParameter);
    
    Beam beam;        
    beam.Initial(train,latticeInterActionPoint,inputParameter);

    beam.Run(train,latticeInterActionPoint,inputParameter);
       
//    MPI_Finalize();
    return 0;
}











