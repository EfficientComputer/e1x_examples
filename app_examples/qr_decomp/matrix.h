typedef struct _matrix
{
    int height;
    int width;
    int *data;
} matrix;

matrix *makeMatrix(int width, int height);
void freeMatrix(matrix *m);
matrix *scaleMatrix(matrix *m, int value);
matrix *copyMatrix(matrix *m);
void printMatrix(matrix *m);
