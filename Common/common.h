/* Represantation of the Matrix inside the program, along with its functions */

struct MatrixItem {
    int row;
    int column;
    double data;
};

struct Matrix {
    int sizex;
    int sizey;
    MatrixItem** items;
};

/* Matrix functions */

Matrix* CreateMatrix(int sizex, int sizey, double** dataz);
Matrix* InsertMatrixItem(Matrix* matrix, MatrixItem* item);
void PrintMatrix(Matrix* matrix);
void ExtractMatrixInformation(Matrix* matrix, double** data);
void LeftShiftMatrix(Matrix* matrix, int blockSize, int initialBlock);
void UpShiftMatrix(Matrix* matrix, int blockSize, int initialBlock);
void RightShiftMatrix(Matrix* matrix, int blockSize, int initialBlock);
void DownShiftMatrix(Matrix* matrix, int blockSize, int initialBlock);
void MultiplyMatrix(Matrix* matrix, Matrix* matrixa, Matrix* matrixb);
double* CreateArrayAsMatrix(int i, int j);
double* RandomizeArray(double* arr, int i, int j);
int ArrayAsEquals(double* arr1, double* arr2, int i, int j);
double** CreateRandomDataSet(int sizex, int sizey);
double* MultiplyMatrixesDouble(double* matrixa, double* matrixb, int size);
double* CreateHorizontalMatrix(int sizex, int sizey);
Matrix* MultiplyMatrixes(Matrix* matrixa, Matrix* matrixb);
double** CreateEmptyDataSet(int sizex, int sizey);
int CheckIfEquals(double* matrixa, double* matrixb, int sizex, int sizey);