//homework2.h
//Daisuke Tanaka

#include "netpbm.h"

void edgeDetection(char *inputFilename, char *sobelFilename, char *cannyFilename);
Image sobel(Image img);
Image canny(Image img);
Matrix convolve(Matrix m1, Matrix m2);
double filterSum(Matrix *input, Matrix *filter,int x, int y);
int isEightNeighbor(Matrix *m1, int x, int y, int value);
void printMatrix(Matrix matrix, int startx, int starty, int width, int height);
