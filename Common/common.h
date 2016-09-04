/* Represantation of the Matrix inside the program, along with its functions */

struct MatrixItem {
    int row;
    int column;
    double* data;
};

struct Matrix {
    int sizex;
    int sizey;
    MatrixItem** items;
};

/* Matrix functions */

Matrix* CreateMatrix(int sizex, int sizey, double** dataz);
Matrix* InsertMatrixItem(Matrix* matrix, MatrixItem* item)
void PrintMatrix(Matrix** matrix);

