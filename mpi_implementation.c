#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include "common.h"

#define TABLE_SIZE 400


double** CreateBlockTables(int sizex, int sizey) {
    double** table = NULL;
    int i = 0;
    
    if(sizex == 0 || sizey == 0) return NULL;

    table = malloc(sizex * sizeof(double*));
    if(table == NULL) return NULL;

    for(i = 0; i < sizex; i++) {
        table[i] = malloc(sizey * sizeof(double));

        if(table[i] == NULL) return NULL;
    }

    return table;
}

void ImplementSlaveMPI(Matrix* matrix, Matrix* matrixa, Matrix* matrixb) {
    
}

void ImplementMasterMPI() {

}

int main(int argc, char* argv[]) {
    int processors = 0;
    int dimension = 0;
    int maxRowBlock = 0;
    int maxColumnBlock = 0;
    int rowBlock = 0;
    int columnBlock = 0;
    int i = 0;
    int j = 0;
    int blocks = 0;
    int rank = 0;
    int current_block = 0;

    double** blocksa = NULL;
    double** blocksb = NULL;
    double** blocksc = NULL;

    double* matrixa = NULL;
    double* matrixb = NULL;
    double* matrixc = NULL;

    GetAndValidateParameters(argc, argv);

    matrixa = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE);
    matrixb = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE);
    matrixc = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processors); 

    if(processors == 0 ) return -1;

    dimension = TABLE_SIZE*TABLE_SIZE/processors;
    blocks = TABLE_SIZE/dimension;

    blocksa = CreateBlockTables(processors, dimension*dimension);
    blocksb = CreateBlockTables(processors, dimension*dimension);
    blocksc = CreateBlockTables(processors, dimension*dimension);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    maxRowBlock = dimension;
    maxColumnBlock = dimension;

    columnBlock = rank % blocks;
    rowBlock = (rank - columnBlock) / blocks;

    for(i = rowBlock * dimension; i < rowBlock * dimension + dimension; i++) {
        for(j = columnBlock * dimension; j < columnBlock * dimension + dimension; j++) {
            
        }
    }
    
    MPI_Finalize();
    
    return 0;
}