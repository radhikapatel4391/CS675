// netpbm_hough.c
// Test and demo program for the Hough transform using the netpbm.c library.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "netpbm.h"

#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define PI 3.141592653589793


// Build a Hough parameter map for matrix mxSpatial for detecting straight lines.
// Rows in this map represent the normal alpha and columns represent the distance d from the origin.
// Increasing the size of the map in each dimension improves the resolution of the corresponding parameter.
Matrix houghTransformCircle(Matrix mxSpatial, int rmin, int rmax)
{
	int m, n, angle, dist,x,y,mapHeight=360,mapWidth=360;
	double alpha, maxD = sqrt((double) (SQR(mxSpatial.height) + SQR(mxSpatial.width)));
	Matrix mxParam = createMatrix(mapHeight, mapWidth);
	Matrix sincos = createMatrix(mapHeight, 2);

	// Generate lookup table for sin and cos values to speed up subsequent computation.
	for (angle = 0; angle < mapHeight; angle++)
	{
		alpha = -0.5*PI + 1.5*PI*(double) angle/(double) mapHeight;
		sincos.map[angle][0] = (double) mapWidth/maxD*sin(alpha);
		sincos.map[angle][1] = (double) mapWidth/maxD*cos(alpha);
	}

	for (m = 0; m < mxSpatial.height; m++)
		for (n = 0; n < mxSpatial.width; n++)
			for (angle = 0; angle < mapHeight; angle++)
			{

                    for(int r=rmin;r<rmax+1;r++)
                    {
                        x= (int) ((m-(r*sincos.map[angle][1]))+0.5);
                        y= (int) ((n-(r*sincos.map[angle][0]))+0.5);
                        if(x>=1 && y>=1)
                            mxParam.map[x][y] += mxSpatial.map[m][n];
                    }
			}

	deleteMatrix(sincos);
	return mxParam;
}



// Test whether entry (m, n) in matrix mx is a local maximum, i.e., is not exceeded by any of its
// maximally 8 neighbors. Return 1 if true, 0 otherwise.
int isLocalMaximum(Matrix mx, int m, int n)
{
	double strength = mx.map[m][n];
	int i, j;
	int iMin = (m == 0)? 0:(m - 1);
	int iMax = (m == mx.height -1)? m:(m + 1);
	int jMin = (n == 0)? 0:(n - 1);
	int jMax = (n == mx.width -1)? n:(n + 1);

	for (i = iMin; i <= iMax; i++)
		for (j = jMin; j <= jMax; j++)
			if (mx.map[i][j] > strength)
				return 0;
	return 1;
}


// Insert a new entry, consisting of vPos, hPos, and strength, into the list of maxima mx.
void insertMaxEntry(Matrix mx, int vPos, int hPos, double strength)
{
	int m, n = mx.width - 1;

	while (n > 0 && mx.map[2][n - 1] < strength)
	{
		for (m = 0; m < 3; m++)
			mx.map[m][n] = mx.map[m][n - 1];
		n--;
	}
	mx.map[0][n] = (double) vPos;
	mx.map[1][n] = (double) hPos;
	mx.map[2][n] = strength;
}


// Delete entry number i from the list of maxima mx.
void deleteMaxEntry(Matrix mx, int i)
{
	int m, n;

	for (n = i; n < mx.width - 1; n++)
		for (m = 0; m < 3; m++)
			mx.map[m][n] = mx.map[m][n + 1];

	mx.map[2][mx.width - 1] = -1.0;
}


// Find the <number> highest maxima in a Hough parameter map that are separated by a Euclidean distance
// of at least <minSeparation> in the map. The result is a three-row matrix with each column representing
// the row, the column, and the strength of one maximum, in descending order of strength.
Matrix findHoughMaxima(Matrix mx, int number, double minSeparation)
{
	int j, m, n, k, l, r, index, do_not_insert;
	double minSepSquare = SQR(minSeparation);
	double strength;
	Matrix maxima = createMatrix(3, number);

	for (j = 0; j < number; j++)
		maxima.map[2][j] = -1.0;

	for (m = 0; m < mx.height; m++)
		for (n = 0; n < mx.width; n++)
		{
			strength = mx.map[m][n];
			index = 0;
			do_not_insert = 0;
			if (strength > 0.0 && strength > maxima.map[2][number - 1] && isLocalMaximum(mx, m, n))
			{
				while (index < number && !do_not_insert)
				{
					if (maxima.map[2][index] > 0.0 && SQR(maxima.map[0][index] - (double) m) + SQR(maxima.map[1][index] - (double) n) < minSepSquare)
						if (strength > maxima.map[2][index])
							deleteMaxEntry(maxima, index);
						else
							do_not_insert = 1;
					index++;
				}
				if (!do_not_insert)
					insertMaxEntry(maxima, m, n, strength);
				index = 0;
			}
		}
	return maxima;
}

/**for Edge detection*/
/** Return sum of filter at point x,y */
double filterSum(Matrix *input, Matrix *filter, int x, int y) {
	double sum = 0.0;
	int anchorx = filter->width % 2 == 0 ? filter->width / 2 - 1 : filter->width / 2;
	int anchory = filter->height % 2 == 0 ? filter->height / 2 - 1 : filter->height / 2;
	for (int row = 0; row < filter->height; row++) {
		for (int column = 0; column < filter->width; column++)
			sum += filter->map[row][column] * input->map[y + row - anchory][x + column - anchorx];
	}
	//Subtract anchor point
	sum -= filter->map[anchory][anchorx] * input->map[y][x];
	return sum;
}
Matrix convolve(Matrix m1, Matrix m2) {
	Matrix output = createMatrix(m1.height, m1.width);
	int anchorx = m2.width % 2 == 0 ? m2.width / 2 - 1 : m2.width / 2;
	int anchory = m2.height % 2 == 0 ? m2.height / 2 - 1 : m2.height / 2;
	for (int y = anchory; y < m1.height - anchory; y++) {
		for (int x = anchorx; x < m1.width - anchorx; x++)
			output.map[y][x] = filterSum(&m1,&m2,x,y);
	}
	return output;
}
Image sobel(Image img) {
	Matrix i = image2Matrix(img);
	Matrix j = image2Matrix(img);
	Matrix m = image2Matrix(img);
	double sobeli[3][3] = {
		{1,2,1},
		{0,0,0},
		{-1,-2,-1}
	};
	double sobelj[3][3] = {
		{1,0,-1},
		{2,0,-2},
		{1,0,-1}
	};
	i = convolve(m, createMatrixFromArray(sobeli,3,3));
	j = convolve(m, createMatrixFromArray(sobelj,3,3));
    //createMatrixFromArray
	for (int y = 0; y < img.height; y++) {
		for (int x = 0; x < img.width; x++)
			m.map[y][x] = sqrt(pow(i.map[y][x],2) + pow(j.map[y][x],2));
	}

	return matrix2Image(m, 1, 1.0);
}
/** Convolve input matrix m1 with filter m2 */

// Read image "desk.ppm" and write Hough transform related output images.
void main()
{
	int i, m, n, m1, n1, m2, n2,k=0,BigCircle=6,smallCircle=4,r1=32,r2=26,max = 0, index,x,y,r; //30,25
	int r1min=0,r1max=0,r2min=10,r2max=10;
    int temp[250][4],large[10][3];

	double gaussFilter[3][3] = {{1.0, 2.0, 1.0}, {2.0, 4.0, 2.0}, {1.0, 2.0, 1.0}};
	Matrix gauss = createMatrixFromArray(&gaussFilter[0][0], 3, 3);
	Image inputImage = readImage("coin.pgm");
	Matrix inputMatrix = image2Matrix(inputImage);
	Image edgeImage, houghImage;
	Matrix edgeMatrix = createMatrix(inputImage.height, inputImage.width);
	Matrix houghMatrix, maxMatrix;
	double maxLength, alpha, dist;

	// Add code for generating edge matrix here!!!
    edgeImage = sobel(inputImage);
    edgeMatrix = image2Matrix(edgeImage);

	writeImage(edgeImage, "coin_edges.pgm");

/**Detection of big coin */
	houghMatrix = houghTransformCircle(edgeMatrix,r1,r1);//Big circle in image //for all radius projection on one image pass rmin,rmax

	houghImage = matrix2Image(houghMatrix, 1, 1.0);
	writeImage(houghImage, "coin_Bigcoin_hough.pgm");

	maxMatrix = findHoughMaxima(houghMatrix, BigCircle, 50.0);
	for (i = 0; i < BigCircle; i++)
		ellipse(houghImage, maxMatrix.map[0][i], maxMatrix.map[1][i], 15, 15, 2, 10, 7, 0, 0, 255, 0);

	writeImage(houghImage, "coin_Bigcoin_hough_max.ppm");

	for (i = 0; i < BigCircle; i++)
	{
        ellipse(inputImage, maxMatrix.map[0][i], maxMatrix.map[1][i], r1, r1, 2, 10, 7, 0, 255, 0, 0);
	}

		writeImage(inputImage, "coin_Bigcoin_hough_circle.ppm");
/** detection of small coin */
     inputImage = readImage("coin.pgm");
    houghMatrix = houghTransformCircle(edgeMatrix,r2,r2);//small circle in image //for all radius projection on one image pass rmin,rmax

	houghImage = matrix2Image(houghMatrix, 1, 1.0);
	writeImage(houghImage, "coin_Smallcoin_hough.pgm");

	maxMatrix = findHoughMaxima(houghMatrix, smallCircle, 50.0);
	for (i = 0; i < smallCircle; i++)
		ellipse(houghImage, maxMatrix.map[0][i], maxMatrix.map[1][i], 15, 15, 2, 10, 7, 0, 0, 255, 0);

	writeImage(houghImage, "coin_Smallcoin_hough_max.ppm");

	maxLength = sqrt((double) (SQR(inputImage.height) + SQR(inputImage.width)));

	for (i = 0; i < smallCircle; i++)
	{
        ellipse(inputImage, maxMatrix.map[0][i], maxMatrix.map[1][i], r2, r2, 2, 10, 7, 255, 0, 0, 0);
	}

		writeImage(inputImage, "coin_Smallcoin_hough_circle.ppm");

/** all Circle detection */
     inputImage = readImage("coin.pgm");
    houghMatrix = houghTransformCircle(edgeMatrix,22,35);//Big circle and small circle togather in one image so pass range and get data

	houghImage = matrix2Image(houghMatrix, 1, 1.0);
	writeImage(houghImage, "coin_hough.pgm");

	maxMatrix = findHoughMaxima(houghMatrix, 10, 50.0);
	for (i = 0; i < 10; i++)
		ellipse(houghImage, maxMatrix.map[0][i], maxMatrix.map[1][i], 15, 15, 2, 10, 7, 0, 0, 255, 0);

	writeImage(houghImage, "coin_hough_max.ppm");

	maxLength = sqrt((double) (SQR(inputImage.height) + SQR(inputImage.width)));

	for (i = 0; i < 10; i++)
	{
        ellipse(inputImage, maxMatrix.map[0][i], maxMatrix.map[1][i], 30, 30, 2, 10, 7, 255, 0, 0, 0);
	}

		writeImage(inputImage, "coin_hough_circle.ppm");


/**....................*/





	//For detailed work.........

    printf("Radius range 0 to 50\n");
	printf("Maxima: radius: x: y:\n");

	//for different radius value let's check maxima and then among them let's find maxima
    for(int j=0;j<25;j++)
    {
        r1min = r1min + 2;
        r1max = r1max + 2;

        houghMatrix = houghTransformCircle(edgeMatrix,r1min,r1max);
        houghImage = matrix2Image(houghMatrix, 1, 1.0);
        maxMatrix = findHoughMaxima(houghMatrix, 10, 50.0);

        for (i = 0; i < 10; i++)
        {
            printf("%0.0lf : %d : %0.0lf :  %0.0lf\n",maxMatrix.map[2][i],r1max,maxMatrix.map[0][i],maxMatrix.map[1][i]);
            temp[k][0]= maxMatrix.map[2][i];
            temp[k][1]= r1max;
            temp[k][2]=maxMatrix.map[0][i];
            temp[k][3]=maxMatrix.map[1][i];
            k++;
        }
        printf("------------------------------\n");

    }
    //find maxima among all possible radius value top 10
    for (int j = 0; j < 10; j++) {
        max = temp[0][0];
        x = temp[0][2];
        y = temp[0][3];
        r = temp[0][1];
        index = 0;
        for (int i =1; i < 250; i++) {
            if (max < temp[i][0]) {
                max = temp[i][0];
                x = temp[i][2];
                y = temp[i][3];
                r = temp[i][1];
                index = i;
            }
        }
        large[j][0]= r;
        large[j][1]= x;
        large[j][2]= y;
        temp[index][0]=0;

    }
    //Print top 10 maxima and respective radius, x,y
     for (int j = 0; j < 10; j++) {
       // printf("%d : %d : %d :  \n",large[j][0],large[j][1],large[j][2]);

     }

    deleteMatrix(edgeMatrix);
	deleteMatrix(houghMatrix);
	deleteImage(inputImage);
	deleteImage(edgeImage);
	deleteImage(houghImage);
}
