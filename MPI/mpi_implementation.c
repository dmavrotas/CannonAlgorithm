#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "common.h"

#define TABLE_SIZE 400

void GetAndValidateParameters(int argc, char* argv[]) {

}

int main(int argc, char* argv[]) {
    Matrix* matrixa = NULL;
    Matrix* matrixb = NULL;
    Matrix* matrixc = NULL;

    GetAndValidateParameters(argc, argv);

    matrixa = CreateMatrix(TABLE_SIZE, TABLE_SIZE, CreateRandomDataSet());
    matrixb = CreateMatrix(TABLE_SIZE, TABLE_SIZE, CreateRandomDataSet());
    matrixc = CreateMatrix(TABLE_SIZE, TABLE_SIZE, CreateRandomDataSet());


    MPI_Init(&argc, &argv);
    
    int processors;

    MPI_Comm_size(MPI_COMM_WORLD, &processors); 
    
    MPI_Finalize();
    
    return 0;
}