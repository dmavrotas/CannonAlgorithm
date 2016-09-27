#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "common.h"

#define TABLE_SIZE 400
#define BLOCK_SIZE 10

double** CreateBlockTables(int processors, int process_block) {
    double*** table = NULL;
    int i = 0;
    int j = 0;
    
    if(sizex == 0 || sizey == 0) return NULL;

    table = malloc(processors * sizeof(double**));
    if(table == NULL) return NULL;

    for(i = 0; i < processors; i++) {
        table[i] = malloc(sizeof(double*) * process_block);
        
        if(table[i] == NULL) return NULL;

        for(j = 0; j < process_block; j++) {
            table[i][j] = malloc((BLOCK_SIZE*BLOCK_SIZE) * sizeof(double));
            
            if(table[i][j] == NULL) return NULL;
        }
    }

    return table;
}

double* Copy_Block(double* matrix, int size) {
    int i = 0;

    double* dest_matrix = malloc((size * size) * sizeof(double));

    if(dest_matrix == NULL) return NULL;

    for(i = 0; i < size * size; i++) {
        dest_matrix[i] = matrix[i];
    }

    return dest_matrix;
}

double** Copy_Blocks(double** matrix, int sizex, int sizey) {
    int i = 0;
    int j = 0;
    double** dest_matrix = NULL;

    dest_matrix = malloc(sizex * sizeof(double*));

    if(dest_matrix == NULL) return NULL;

    for(i = 0; i < sizey*sizey; i++) {
        dest_matrix[i] = malloc((sizey*sizey) * sizeof(double));

        if(dest_matrix[i] == NULL) return NULL;
    }

    for(i = 0; i < sizex; i++) {
        for(j = 0; j < sizey*sizey; j++) {
            dest_matrix[i][j] = matrix[i][j];
        }
    }
}

double** Create_Block(int process_block) {
    int i = 0;
    double** matrix = NULL;

    matrix = malloc(process_block * sizeof(double*));
    if(matrix == NULL) return NULL;

    for(i = 0; i < process_block; i++) {
        matrix[i] = malloc((BLOCK_SIZE * BLOCK_SIZE) * sizeof(double));

        if(matrix[i] == NULL) return NULL;
    }

    return matrix;
}

void Shift_Left(double** matrix, int sizex, int sizey) {
    int i = 0;
    int j = 0;
    int k = 0;
    int step = 0;

    int row = 0;
    int col = 0;
    int position = 0;

    step = (int)sqrt(sizex);
    int out[step];

    for(i = 0, col = 0; i < sizex; i += step, col++) {
        for(k = 0; k < sizey * sizey; k++) {
            for(j = i, row = 0; j < (i + step); j++, row++) {
                position = j + 1;
                if(position >= (step + i)) {
                    position = col * step;
                }

                out[r] = matrix[position][k];
            }
            for(j = i, row = 0; j < i + step; j++, row++) {
                matrix[j][k] = out[row];
            }
        }
    }
}

void Shift_Up(double** matrix, int sizex, int sizey) {
    int i = 0;
    int j = 0;
    int k = 0;
    int step = 0;

    int row = 0;
    int col = 0;ß

    step = (int)sqrt(sizex);
    int out[step];

    for(k = 0; k < sizey * sizey; k++) {
        for(i = 0; i < sizex; i++) {
            out[i] = matrix[(i+step) % sizex][k];
        }
        for(i = 0; i < sizex; i++) {
            matrix[j][k] = out[i];
        }
    }
}
 
void RunMPI(int size, double** blocksa, double** blocksb, double** blocksc, 
            double* matrix, int process_block, MPI_Comm comm) {
    int i = 0;
    int j = 0;
    int limit = 0;
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

    double** matrixa_buffer = NULL;
    double** matrixb_buffer = NULL;
    double** matrix_copy = NULL;

    double** matrixa = NULL;
    double** matrixb = NULL;

    double* dest_matrixa = NULL;
    double* dest_matrixb = NULL;

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

    matrixa_buffer = Create_Block(process_block);
    matrixb_buffer = Create_Block(process_block);
    matrix_copy  = Creat_Block(process_block);

    matrix_copy = Copy_Blocks(&blocksb[0], process_block, BLOCK_SIZE);

    limit = (int)sqrt(process_block);

    matrixa = malloc((BLOCK_SIZE * BLOCK_SIZE) * sizeof(double));
    if(matrixa == NULL) return;

    matrixb = malloc((BLOCK_SIZE * BLOCK_SIZE) * sizeof(double));
    if(matrixb == NULL) return;

    dest_matrixa = malloc((BLOCK_SIZE * BLOCK_SIZE) * sizeof(double));
    if(dest_matrixa == NULL) return;

    dest_matrixb = malloc((BLOCK_SIZE * BLOCK_SIZE) * sizeof(double));
    if(dest_matrixb == NULL) return;

    matrixa_buffer = Copy_Blocks(&blocksa[0], process_block, BLOCK_SIZE);

    for(i = 0; i < limit; i++) {
        matrixb_buffer = Copy_Blocks(&matrix_copy[0], process_block, BLOCK_SIZE);
        int rowdisplay = [coords[0] + i * sqrt(np));

        if((coords[1] - rowdisplay) < 0) {
            Shift_Left(&matrixa_buffer[0], process_blocks, BLOCK_SIZE);
        }

        for(j = 0; j < limit; j++) {
            int columndisplay = (coords[1] + j * sqrt(np));

            if((coords[0] - columndisplay) > 0) {
                Shift_Up(&matrixb_buffer[0], process_block, BLOCK_SIZE);
            }

            matrixa = Copy_Block(&matrixa_buffer[i * limit * j][0], BLOCK_SIZE);
            matrixb = Copy_Block(&matrixb_buffer[i * limit * j][0], BLOCK_SIZE);

            MPI_Cart_shift(commloc, 1, -rowdisplay, &shiftright, &shiftleft);
            MPI_Sendrecv(matrixa, BLOCK_SIZE * BLOCK_SIZE, MPI_INT, shiftleft, 1, 
                        dest_matrixa, BLOCK_SIZE*BLOCK_SIZE, MPI_INT, shiftright, 
                        commloc, &status);
        }
    }

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
    int dimension_block = 0;
    int process_block = 0;
    int row_pq = 0;
    int column_pq = 0;
    int row_mn = 0;
    int column_mn = 0;
    int rowIndex = 0;
    int columnIndex = 0;
    int rowOffset = (TABLE_SIZE / dimension_block);
    int rowBegin = 0;
    int rowEnd = 0;
    int columnBegin = 0;
    int columnEnd = 0;
    int row = 0;
    int column = 0;
    int step = 0;
    double start = 0;
    double end = 0;

    double*** blocksa = NULL;
    double*** blocksb = NULL;
    double*** blocksc = NULL;

    double* matrixa = NULL;
    double* matrixb = NULL;
    double* matrixc = NULL;

    matrixa = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE, 1);
    matrixb = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE, 1);
    matrixc = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE, 0);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processors); 
    if(processors == 0 ) return -1;

    dimension = sqrt(processors);
    dimension_block = (TABLE_SIZE / BLOCK_SIZE);
    process_block = (dimension_block / dimension) * (dimension_block / dimension);

    blocksa = CreateBlockTables(processors, process_block);
    blocksb = CreateBlockTables(processors, process_block);
    blocksc = CreateBlockTables(processors, process_block);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    row_pq = (rank / dimension);
    column_pq = (rank % dimension);

    for(row_mn = row_pq; column_mn < dimension_block; row_mn += dimension) {
        rowBegin = row_mn * rowOffset;
        rowEnd = rowBegin + rowOffset;

        for(column_mn = column_pq; column_mn < dimension_block; column_mn += dimension, step++) {
            columnBegin = column_mn * rowOffset;
            columnEnd = columnBegin + rowOffset;

            for(rowIndex = rowBegin, row = 0; rowIndex < rowEnd; rowIndex++, row++) {
                for(columnIndex = columnBegin, step = 0; columnIndex < columnEnd; row++, step++) {
                    blocksa[rank][step][row * BLOCK_SIZE + step] = matrixa[rowIndex * TABLE_SIZE + columnIndex];
                    blocksb[rank][step][row * BLOCK_SIZE + step] = matrixb[rowIndex * TABLE_SIZE + columnIndex];
                    blocksb[rank][step][row * BLOCK_SIZE + step] = 0;
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    start = MPI_Wtime();

    end = MPI_Wtime();

    if(rank == 0) {
        double* matrixd = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE, 0);
    }

    MPI_Finalize();

    return 0;
}