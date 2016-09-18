#include <stdio.h>
#include <stdlib.h>

/* This function creates a matrix and returns it */
Matrix* CreateMatrix(int sizex, int sizey, double** dataz) {
    int i = 0;
    int j = 0;

    if(sizex == 0 || sizey == 0) return NULL;

    if(dataz == NULL) return NULL;

    Matrix* matrix = NULL;

    matrix = malloc(sizeof(Matrix));
    if(matrix == NULL) return NULL;

    matrix->sizex = sizex;
    matrix->sizey = sizey;
    matrix->items = malloc(sizex * sizeof(MatrixItem*));
    if(matrix->items == NULL) return NULL;

    for(i = 0; i < sizex; i++) {
        matrix->items[i] = malloc(sizey * sizeof(MatrixItem));
        if(matrix->items[i] == NULL) return NULL;
        for(j = 0; j < sizey; j++) {
            if(matrix->items[i][j] != NULL) {
                matrix->items[i][j]->row = i;
                matrix->items[i][j]->column = j;
                matrix->items[i][j]->data = dataz[i][j];
            } 
        }
    }

    return matrix;
}

/* This function inserts a MatrixItem inside a Matrix and returns the matrix */
Matrix* InsertMatrixItem(Matrix* matrix, MatrixItem* item) {
    int i = 0;
    int j = 0;
    
    if(matrix == NULL) return NULL;

    if(item == NULL) return NULL;

    for(i = 0; i < matrix->sizex; i++) {
        for(j = 0; j < matrix->sizey; j++) {
            matrix->items[i][j] = item;
        }
    }

    return matrix;
}

/* This function just prints a matrix */
void PrintMatrix(Matrix* matrix) {
    int i = 0;
    int j = 0;

    if(matrix == NULL) return;

    if(matrix->sizex == 0 || matrix->sizey == 0) return;

    for(i = 0; i < matrix->sizex; i++) {
        for(j = 0; j < matrix->sizey; j++) {
            if((matrix->items[i][j])->data[0] != NULL) {
                printf(" %2f ", (matrix->items[i][j])->data[0]);
            }
        }
        printf("\n");
    }
}

/* This function gets the data matrix from the Matrix */
void ExtractMatrixInformation(Matrix* matrix, double** data) {
    int i = 0;
    int j = 0;

    for(i = 0; i < matrix->sizex; i++) {
        for(j = 0; j < matrix->sizey; j++) {
            matrix->items[i][j]->data = data[i][j];
        }
    }
}

/* This function shifts the matrix left */
void LeftShiftMatrix(Matrix* matrix, int blockSize, int initialBlock) {
    int i = 0;
    int j = 0;
    int k = 0;
    int s = 0;
    int step = blockSize;
    Matrix* middle = NULL;

    middle = CreateMatrix(1, matrix->sizey, matrix);
    for(k = 0, s = 0; k < matrix->sizey; k+= blockSize, s++) {
        for(i = k; i < (k + blockSize); i++) {
            if(initialBlock > 0) {
                step = s * blockSize;
            }
            for(j = 0; j < matrix->sizey; j++) {
                middle->items[0][j]->data = matrix->items[i][(j + step)
                                            % matrix->sizey];
            }
            for(j = 0; j < matrix->sizey; j++) {
                matrix->items[i][j]->data = middle->items[0][j]->data;
            }
        }
    }
}

/* This function shifts the matrix up */
void UpShiftMatrix(Matrix* matrix, int blockSize, int initialBlock) {
    int i = 0;
    int j = 0;
    int k = 0;
    int s = 0;
    int step = blockSize;
    Matrix* middle = NULL;

    middle = CreateMatrix(1, matrix->sizex, matrix);
    for(k = 0, s = 0; k < matrix->sizex; k+= blockSize, s++) {
        for(i = k; i < (k + blockSize); i++) {
            if(initialBlock > 0) {
                step = s * blockSize;
            }
            for(j = 0; j < matrix->sizex; j++) {
                middle->items[0][j]->data = matrix->items[(j + step)
                                            % matrix->sizex][i];
            }
            for(j = 0; j < matrix->sizex; j++) {
                matrix->items[i][j]->data = middle->items[0][j]->data;
            }
        }
    }
}

/* This function shifts the matrix right */
void RightShiftMatrix(Matrix* matrix, int blockSize, int initialBlock) {
    int i = 0;
    int j = 0;
    int k = 0;
    int s = 0;
    int step = blockSize;
    Matrix* middle = NULL;

    middle = CreateMatrix(1, matrix->sizey, matrix);
    for(k = 0, s = 0; k < matrix->sizey; k+= blockSize, s++) {
        for(i = k; i < (k + blockSize); i++) {
            if(initialBlock > 0) {
                step = s * blockSize;
            }
            for(j = 0; j < matrix->sizey; j++) {
                middle->items[matrix->sizey - 1][j]->data = matrix->items[(j + step)
                                            % matrix->sizey][i];
            }
            for(j = 0; j < matrix->sizey; j++) {
                matrix->items[matrix->sizey - 1][j]->data = middle->items[0][j]->data;
            }
        }
    }
}
 
 /* This function shifts the matrix down */
void DownShiftMatrix(Matrix* matrix, int blockSize, int initialBlock) {
    int i = 0;
    int j = 0;
    int k = 0;
    int s = 0;
    int step = blockSize;
    Matrix* middle = NULL;

    middle = CreateMatrix(1, matrix->sizex, matrix);
    for(k = 0, s = 0; k < matrix->sizex; k+= blockSize, s++) {
        for(i = k; i < (k + blockSize); i++) {
            if(initialBlock > 0) {
                step = s * blockSize;
            }
            for(j = 0; j < matrix->sizex; j++) {
                middle->items[matrix->sizex - 1][j]->data = matrix->items[(j + step)
                                            % matrix->sizex][i];
            }
            for(j = 0; j < matrix->sizex; j++) {
                matrix->items[matrix->sizex - 1][j]->data = middle->items[0][j]->data;
            }
        }
    }
}

/* This function calculates the product of the table elements */
void MultiplyMatrix(Matrix* matrix, Matrix* matrixa, Matrix* matrixb) {
    int i = 0;
    int j = 0;
    int k = 0;

    for(i = 0; i < matrixa->row; i++) {
        for(j = 0; j < matrixb->column; j++) {
            for(k = 0; k < matrixa->column) {
                (matrix->items[i][j])->data += (matrixa->items[i][k])->data * 
                                                (matrixb->items[k][j])->data;
            }
        }
    }
}

/* This function creates an array that can be matched with a matrix */
double* CreateArrayAsMatrix(int i, int j) {
    double* arr = malloc(i * j sizeof(double));
    if(arr == NULL) return NULL;
    
    return arr;
}

/* This function randomizes an array */
double* RandomizeArray(double* arr, int i, int j) {
    int k = 0;

    if(arr == NULL) return NULL;

    for(k = 0; k < i*j; k++) {
        arr[k] = rand() % 10 + 1;
    }
}

/* Simple AsEquals function for arrays */
int ArrayAsEquals(double* arr1, double* arr2, int i, int j) {
    int k = 0;

    if(arr1 == NULL) return NULL;
    if(arr2 == NULL) return NULL;

    for(k = 0; k < i*j; k++) {
        if(arr1[i] != arr2[i]) {
            return 0;
        }
    }

    return 1;
}
