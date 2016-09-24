#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
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

void RunMPI(int size, double* blocksa, double* blocksb, double* blocksc, double* matrix, MPI_Comm comm) {
    int i = 0;
    int local = 0;
    int np = 0;
    int dimensions[2];
    int periods[2];
    int rank = 0;
    int rank2D = 0;
    int coords[2];
    int uprank = 0;
    int downrank = 0;
    int rightrank = 0;
    int leftrank = 0;
    int coords2[2];
    int shiftright = 0;
    int shiftleft = 0;
    int indentation = 0;
    int size_n = size;
    int indentation_column = 0;

    printf("Inside RUNMPI \n");

    MPI_Status status;
    MPI_Comm commloc;

    MPI_Comm_size(comm, &np);
    MPI_Comm_rank(comm, &rank);

    dimensions[0] = dimensions[1] = sqrt(np);
    periods[0] = periods[1] = 1;

    MPI_Cart_create(comm, 2, dimensions, periods, 1, &commloc);

    MPI_Comm_rank(commloc, &rank2D);
    MPI_Cart_coords(commloc, rank2D, 2, coords);

    MPI_Cart_shift(commloc, 1, -1, &rightrank, &leftrank);
    MPI_Cart_shift(commloc, 0, -1, &downrank, &uprank);

    local = size/dimensions[0];

    MPI_Cart_shift(commloc, 1, -coords[0], &shiftright, &shiftleft);
	MPI_Sendrecv_replace(blocksa, local*local, MPI_INT, shiftleft, 1, shiftright, 1, commloc, &status);
	MPI_Cart_shift(commloc, 0, -coords[1], &shiftright, &shiftleft);
	MPI_Sendrecv_replace(blocksb, local*local, MPI_INT, shiftleft, 1, shiftright, 1, commloc, &status);

    for (i=0; i<dimensions[0]; i++)
	{
		MatrixMultiply(local, blocksa, blocksb, blocksc);
		
		MPI_Sendrecv_replace(blocksa, local*local, MPI_INT, leftrank, 1, rightrank, 1, commloc, &status);
		
		MPI_Sendrecv_replace(blocksb, local*local, MPI_INT, uprank, 1, downrank, 1, commloc, &status);
	}

    MPI_Cart_shift(commloc, 1, +coords[0], &shiftright, &shiftleft);
	MPI_Sendrecv_replace(blocksa, local*local, MPI_INT,shiftleft, 1, shiftright, 1, commloc, &status);
	MPI_Cart_shift(commloc, 0, +coords[1], &shiftright, &shiftleft);
	MPI_Sendrecv_replace(blocksb, local*local, MPI_INT, shiftleft, 1, shiftright, 1, commloc, &status);
	MPI_Comm_free(&commloc);

    for(indentation = 0; indentation < local; indentation++) {
        for(indentation_column = 0; indentation_column < local; indentation_column++) {
            matrix[(coords[0]*local+indentation) + (coords[1]*local+indentation_column)] =
                    blocksc[indentation*local+indentation_column];
        }
    }

    if(rank != 0) {
        MPI_Reduce(matrix, matrix, size_n*size_n, MPI_INT, MPI_SUM, 0, comm);
    }
    else {
        MPI_Reduce(MPI_IN_PLACE, matrix, size_n*size_n, MPI_INT, MPI_SUM, 0, comm);
    }
}

int ImplementMasterMPI(int rank, double* matrixa, double* matrixb, double* matrixc) {
    if(rank == 0) {
        double* matrixd = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE);
        matrixd = MultiplyMatrixesDouble(matrixa, matrixb, TABLE_SIZE);

        if(matrixd == NULL) {
            return -2; 
        }

        int equity = CheckIfEquals(matrixd, matrixc, TABLE_SIZE, TABLE_SIZE);

        if(equity == 1) return 1;
    }

    return -1;
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
    double start = 0;
    double end = 0;

    double** blocksa = NULL;
    double** blocksb = NULL;
    double** blocksc = NULL;

    double* matrixa = NULL;
    double* matrixb = NULL;
    double* matrixc = NULL;

    matrixa = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE);
    matrixb = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE);
    matrixc = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE);
    printf("Start \n");
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processors); 
    printf("Initialization successful \n");
    if(processors == 0 ) return -1;

    dimension = TABLE_SIZE/sqrt(processors);
    printf("Dimension : %d \n", dimension);
    blocks = TABLE_SIZE/dimension;

    blocksa = CreateBlockTables(processors, dimension*dimension);
    blocksb = CreateBlockTables(processors, dimension*dimension);
    blocksc = CreateBlockTables(processors, dimension*dimension);

    printf("Blocks created \n");
    printf("Blocks : %d \n", blocks);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    maxRowBlock = dimension;
    maxColumnBlock = dimension;

    columnBlock = rank % blocks;
    rowBlock = (rank - columnBlock) / blocks;

    printf("Ranks given \n");

    for(i = rowBlock * dimension; i < rowBlock * dimension + dimension; i++) {
        for(j = columnBlock * dimension; j < columnBlock * dimension + dimension; j++) {
            blocksa[rank][current_block] = matrixa[i * TABLE_SIZE + j];
            blocksb[rank][current_block] = matrixb[i * TABLE_SIZE + j];
            blocksc[rank][current_block] = 0;
            current_block++;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    start = MPI_Wtime();
    
    RunMPI(TABLE_SIZE, blocksa[rank], blocksb[rank], blocksc[rank], &matrixc[0], MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    return ImplementMasterMPI(rank, matrixa, matrixb, matrixc);

    MPI_Finalize();
    
    return 0;
}