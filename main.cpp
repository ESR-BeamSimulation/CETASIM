#include "stdio.h"
#include "string.h"
#include "mpi.h"
#include "LongImpSingalBunch.h"
#include <fstream>
#include <stdlib.h> 
#include <iostream>
#include <vector>

#include "Global.h"
#include "Bunch.h"
#include "Train.h"
#include "LatticeInterActionPoint.h"
#include "Beam.h"



using namespace std;
using std::vector;
using std::complex;

int main(int argc,char *argv[])
{

    MPI_Init(&argc,&argv);
    MPI_Status status;
    
    MPI_Comm_size(MPI_COMM_WORLD,&numProcess);
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

    
    LatticeInterActionPoint latticeInterActionPoint;
    latticeInterActionPoint.Initial();

    Train train;
    train.Initial();

	

    Beam beam;
    beam.Initial(train,latticeInterActionPoint);


    beam.Run(train,latticeInterActionPoint);
    

    MPI_Finalize();
    return 0;
}











