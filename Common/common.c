#include <stdio.h>
#include <stdlib.h>

/* This function creates a matrix and returns it */
Matrix* CreateMatrix(int sizex, int sizey, double* dataz) {
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
                matrix->items[i][j]->data = dataz;
            } 
        }
    }
}

/* This function inserts a MatrixItem inside a Matrix and returns the matrix */
Matrix* InsertMatrixItem(Matrix** matrix, MatrixItem* item) {

}

/* This function just prints a matrix */
void PrintMatrix(Matrix** matrix) {
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

