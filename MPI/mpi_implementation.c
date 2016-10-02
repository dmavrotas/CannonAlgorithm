#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "common.h"

#define TABLE_SIZE 400
#define BLOCK_SIZE 10

double*** CreateBlockTables(int processors, int process_block) {
    double*** table = NULL;
    int i = 0;
    int j = 0;

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


void Copy_Block(double* matrix, double* dest_matrix, int size) {
    int i = 0;

    if(matrix == NULL) return;

    if(dest_matrix == NULL) return;

    for(i = 0; i < size * size; i++) {
        dest_matrix[i] = matrix[i];
    }
}

double** Copy_Blocks(double** matrix, int sizex, int sizey) {
    int i = 0;
    int j = 0;
    double** dest_matrix = NULL;

    dest_matrix = malloc(sizex * sizeof(double*));

    if(dest_matrix == NULL) return NULL;

    for(i = 0; i < sizex; i++) {
        dest_matrix[i] = malloc((sizey*sizey) * sizeof(double));

        if(dest_matrix[i] == NULL) return NULL;
    }

    for(i = 0; i < sizex; i++) {
        for(j = 0; j < sizey*sizey; j++) {
            dest_matrix[i][j] = matrix[i][j];
        }
    }

    return dest_matrix;
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

                out[row] = matrix[position][k];
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
    int col = 0;

    step = (int)sqrt(sizex);
    int out[sizex];

    for(k = 0; k < sizey * sizey; k++) {
        for(i = 0; i < sizex; i++) {
            out[i] = matrix[(i+step) % sizex][k];
        }
        for(i = 0; i < sizex; i++) {
            matrix[j][k] = out[i];
        }
    }
}

void Update(int size, double* matrixa, double* matrixb, double* matrixc) {
    int i = 0;
    int j = 0;
    int k = 0;

    for(i = 0; i < size; i++) {
        for(j = 0; j < size; j++) {
            for(k = 0; k < size; k++) {
                matrixc[i * size + j] += matrixa[i * size + k] * matrixb[k * size + j];
            }
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

    double* matrixa = NULL;
    double* matrixb = NULL;

    double* dest_matrixa = NULL;
    double* dest_matrixb = NULL;

    MPI_Status status;
    MPI_Comm commloc;

    MPI_Comm_size(comm, &np);
    MPI_Comm_rank(comm, &rank);

    printf("Rank achieved  \n");

    dimensions[0] = dimensions[1] = sqrt(np);
    periods[0] = periods[1] = 1;

    MPI_Cart_create(comm, 2, dimensions, periods, 1, &commloc);

    MPI_Comm_rank(commloc, &rank2D);
    MPI_Cart_coords(commloc, rank2D, 2, coords);

    MPI_Cart_shift(commloc, 1, -1, &rightrank, &leftrank);
    MPI_Cart_shift(commloc, 0, -1, &downrank, &uprank);

    printf("MPI World initialized \n");

    matrixa_buffer = Create_Block(process_block);
    matrixb_buffer = Create_Block(process_block);
    matrix_copy  = Create_Block(process_block);

    matrix_copy = Copy_Blocks(&blocksb[0], process_block, BLOCK_SIZE);

    printf("Step 1 \n");

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

    printf("Step 2 \n");

    for(i = 0; i < limit; i++) {
        printf("i = %d \n", i);
        matrixb_buffer = Copy_Blocks(&matrix_copy[0], process_block, BLOCK_SIZE);
        int rowdisplay = (coords[0] + i * sqrt(np));

        if((coords[1] - rowdisplay) < 0) {
            Shift_Left(&matrixa_buffer[0], process_block, BLOCK_SIZE);
        }

        printf("First Shift successful \n");

        for(j = 0; j < limit; j++) {
            int columndisplay = (coords[1] + j * sqrt(np));

            if((coords[0] - columndisplay) < 0) {
                Shift_Up(&matrixb_buffer[0], process_block, BLOCK_SIZE);
            }

            printf("i = %d - j = %d \n", i, j);

            printf("Copy Block 1 \n");
            Copy_Block(&matrixa_buffer[i * limit + j][0], &matrixa[0], BLOCK_SIZE);
            Copy_Block(&matrixb_buffer[i * limit + j][0], &matrixb[0], BLOCK_SIZE);
            printf("Copy Block 2 \n");

            MPI_Cart_shift(commloc, 1, -rowdisplay, &shiftright, &shiftleft);
            MPI_Sendrecv(matrixa, BLOCK_SIZE * BLOCK_SIZE, MPI_INT, shiftleft, 1,
                        dest_matrixa, BLOCK_SIZE*BLOCK_SIZE, MPI_INT, shiftright, 1,
                        commloc, &status);

            MPI_Cart_shift(commloc, 0, -columndisplay, &shiftleft, &shiftright);
            MPI_Sendrecv(matrixb, BLOCK_SIZE * BLOCK_SIZE, MPI_INT, shiftleft, 1,
                        dest_matrixb, BLOCK_SIZE * BLOCK_SIZE, MPI_INT, shiftright,
                        1, commloc, &status);
            printf("Copy Block 3 \n");
            Copy_Block(&dest_matrixa[0], &blocksa[i * limit + j][0], BLOCK_SIZE);
            Copy_Block(&dest_matrixb[0], &blocksb[i * limit + j][0], BLOCK_SIZE);
            printf("Copy Block 4 \n");
            printf("Main MPI successful , i = %d - j = %d \n", i, j);
        }
    }

    printf("Next step \n");

    int sqrtIndex = (int)sqrt(np);
    int dimension = size / BLOCK_SIZE;
    int maximum = dimension / sqrtIndex;

    int ia = 0;
    int ja = 0;
    int i_prime = 0;
    int j_prime = 0;

    for(ia = 0; i < sqrtIndex; ia++) {
        for(i = 0; i < maximum; i++) {
            for(j = 0; j < maximum; j++) {
                for(ja = 0; ja < maximum; ja++) {
                    i_prime = (j % dimension + ja) % (dimension / sqrtIndex);
                    j_prime = (i % dimension + ja) % (dimension / sqrtIndex);
                    Update(BLOCK_SIZE, &blocksa[i * maximum + j_prime][0],
                            &blocksb[i_prime * maximum + j][0],
                            &blocksc[i * maximum + j][0]);
                }
            }
        }

        if(coords[1] == 0) {
            Shift_Left(&blocksa[0], process_block, BLOCK_SIZE);
        }

        for(i = 0; i < process_block; i++) {
            MPI_Sendrecv_replace(blocksa[i], BLOCK_SIZE * BLOCK_SIZE, MPI_INT,
                                leftrank, 1, rightrank, 1, commloc, &status);
        }

        if(coords[0] == 0) {
            Shift_Up(&blocksb[0], process_block, BLOCK_SIZE);
        }

        for(i = 0; i < process_block; i++) {
            MPI_Sendrecv_replace(blocksb[i], BLOCK_SIZE * BLOCK_SIZE, MPI_INT,
                                uprank, 1, downrank, 1, commloc, &status);
        }
    }

    MPI_Comm_free(&commloc);

    int process_number = 0;
    int blockIndex = 0;
    int position = 0;

    process_number = rank;

    int column_process = process_number % (int)sqrt(np);
    int row_process = (process_number - column_process) / sqrt(np);

    for(blockIndex = 0; blockIndex < process_block; blockIndex++) {
        int column_block = blockIndex % (int)sqrt(process_block);
        int row_block = (blockIndex - column_block) / sqrt(process_block);

        for(position = 0; position < (BLOCK_SIZE * BLOCK_SIZE); position++) {
            int column_position = position % BLOCK_SIZE;
            int row_position = (position - column_position) / BLOCK_SIZE;

            int big_row = row_block * sqrt(np) + row_process;
            int big_column = column_block * sqrt(np) + column_process;
            matrix[big_row * size * BLOCK_SIZE + row_position * size + big_column * BLOCK_SIZE + column_position] =
                blocksc[blockIndex][position];
        }
    }

    if(rank != 0) {
        MPI_Reduce(matrix, matrix, np * process_block * BLOCK_SIZE * BLOCK_SIZE, MPI_INT, MPI_SUM, 0, comm);
    }
    else {
        MPI_Reduce(MPI_IN_PLACE, matrix, np * process_block * BLOCK_SIZE * BLOCK_SIZE, MPI_INT, MPI_SUM, 0, comm);
    }
}

int ImplementMasterMPI(int rank, double* matrixa, double* matrixb, double* matrixc) {
    if(rank == 0) {
        double* matrixd = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE, 0);
        matrixd = MultiplyMatrixesDouble(matrixa, matrixb, TABLE_SIZE);

        if(matrixd == NULL) {
            return -2;
        }

        int equity = CheckIfEquals(matrixd, matrixc, TABLE_SIZE, TABLE_SIZE);

        if(equity == 1) return 0;
        return 0;
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
    int rowBegin = 0;
    int rowEnd = 0;
    int columnBegin = 0;
    int columnEnd = 0;
    int row = 0;
    int column = 0;
    int step = 0;
    double start = 0;
    double end = 0;
    int rank = 0;

    double*** blocksa = NULL;
    double*** blocksb = NULL;
    double*** blocksc = NULL;

    double* matrixa = NULL;
    double* matrixb = NULL;
    double* matrixc = NULL;
    printf("Starting... \n");
    matrixa = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE, 1);
    matrixb = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE, 1);
    matrixc = CreateHorizontalMatrix(TABLE_SIZE, TABLE_SIZE, 0);
    printf("Tables created \n");
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processors);
    if(processors == 0 ) return -1;
    printf("Initialization \n");
    dimension = sqrt(processors);
    dimension_block = (TABLE_SIZE / BLOCK_SIZE);
    process_block = (dimension_block / dimension) * (dimension_block / dimension);

    int rowOffset = (TABLE_SIZE / dimension_block);

    blocksa = CreateBlockTables(processors, process_block);
    blocksb = CreateBlockTables(processors, process_block);
    blocksc = CreateBlockTables(processors, process_block);

    printf("Blocks created!!! \n");

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
                for(columnIndex = columnBegin, column = 0; columnIndex < columnEnd; columnIndex++, column++) {
                    blocksa[rank][step][row * BLOCK_SIZE + column] = matrixa[rowIndex * TABLE_SIZE + columnIndex];
                    blocksb[rank][step][row * BLOCK_SIZE + column] = matrixb[rowIndex * TABLE_SIZE + columnIndex];
                    blocksb[rank][step][row * BLOCK_SIZE + column] = 0;
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    start = MPI_Wtime();

    printf("GO MPI \n");

    RunMPI(TABLE_SIZE, &blocksa[rank][0], &blocksb[rank][0], &blocksc[rank][0], 
            &matrixc[0], process_block, MPI_COMM_WORLD);

    end = MPI_Wtime();

    int result = 0;

    if(rank == 0) {
        result = ImplementMasterMPI(rank, matrixa, matrixb, matrixc);
    }

    MPI_Finalize();

    return result;
}