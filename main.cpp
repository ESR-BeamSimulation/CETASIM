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
#include"ReadInputSettings.h"


using namespace std;
using std::vector;
using std::complex;

int main(int argc,char *argv[])
{

    MPI_Init(&argc,&argv);
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD,&numProcess);
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);


    ReadInputSettings inputParameter;
    inputParameter.ParamRead();


    LatticeInterActionPoint latticeInterActionPoint;
    latticeInterActionPoint.Initial(inputParameter);
    latticeInterActionPoint.InitialLattice(inputParameter);

    Train train;
    train.Initial(inputParameter);


    Beam beam;
    beam.Initial(train,latticeInterActionPoint,inputParameter);


    beam.Run(train,latticeInterActionPoint,inputParameter);

    MPI_Finalize();
    return 0;
}











