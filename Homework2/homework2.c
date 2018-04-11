/**
	HW 2 CS 675 Daisuke Tanaka
	Type make to compile the program.
	./homework1 to run the assignment
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "homework2.h"

void main()
{
    printf("Radhika");
	edgeDetection("desk.ppm", "lennasobel.pgm", "lennacanny.pbm");
	//edgeDetection("sample.ppm", "samplesobel.pgm", "samplecanny.pbm");
	//edgeDetection("valve.ppm", "valvesobel.pgm", "valvecanny.pbm");
}

/** Output sobel and canny output */
void edgeDetection(char *inputFilename, char *sobelFilename, char *cannyFilename) {
	Image input = readImage(inputFilename);
	printf("RD");
	writeImage(sobel(input), sobelFilename);
	writeImage(canny(input), cannyFilename);
	deleteImage(input);
}

/** Apply Sobel filter to image */
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

	for (int y = 0; y < img.height; y++) {
		for (int x = 0; x < img.width; x++)
			m.map[y][x] = sqrt(pow(i.map[y][x],2) + pow(j.map[y][x],2));
	}
	printMatrix(m,50,50,50,50);
	return matrix2Image(m, 1, 1.0);
}

/** Apply Canny Edge Detection to image */
Image canny(Image img) {
	Matrix output = image2Matrix(img);
	Matrix magnitude = createMatrix(img.height, img.width);
	Matrix orientation = createMatrix(img.height, img.width);
	Matrix p = createMatrix(img.height, img.width);
	Matrix q = createMatrix(img.height, img.width);
	Matrix e = createMatrix(img.height, img.width);
	Matrix sectors = createMatrix(img.height, img.width);

	//Hysteresis thresholds
	double lowthreshold = 50;
	double highthreshold = 90;

	double gaussian[5][5] = {
		{(double)1/273,(double)4/273,(double)7/273,(double)4/273,(double)1/273},
		{(double)4/273,(double)16/273,(double)26/273,(double)16/273,(double)4/273},
		{(double)7/273,(double)26/273,(double)41/273,(double)26/273,(double)7/273},
		{(double)4/273,(double)16/273,(double)26/273,(double)16/273,(double)4/273},
		{(double)1/273,(double)4/273,(double)7/273,(double)4/273,(double)1/273}
	};
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

	//Smooth Image
	output = convolve(output, createMatrixFromArray(gaussian,5,5));

	//Sobel FilteOutput
	p = convolve(output, createMatrixFromArray(sobeli,3,3));
	q = convolve(output, createMatrixFromArray(sobelj,3,3));

	//compute magnitude, orientation and sectors of gradient
	for (int y = 0; y < img.height; y++) {
		for (int x = 0; x < img.width; x++) {
			double degree = atan2(q.map[y][x], p.map[y][x]) * 180.0 / M_PI;
			while(degree < 0.0)
				degree += 360;
			degree = fmod(degree,180.0);
			//printf("degree: %f mod: %f \n", degree, mod);
			magnitude.map[y][x] = sqrt(pow(p.map[y][x],2) + pow(q.map[y][x],2));
			orientation.map[y][x] = degree;
			//printf("p: %f q: %f magnitude: %f degree: %f\n", p.map[y][x], q.map[y][x], magnitude.map[y][x], degree);
			if ((degree >=0.0 && degree <= 22.5) || (degree >= 157.5 && degree <= 180))
				sectors.map[y][x] = 0;
			else if (degree >=22.5 && degree <= 67.5)
				sectors.map[y][x] = 1;
			else if (degree >=67.5 && degree <= 112.5)
				sectors.map[y][x] = 2;
			else if (degree >=112.5 && degree <= 157.5)
				sectors.map[y][x] = 3;
		}
	}

	//Non-Maxima supression
	for (int y = 1; y < img.height - 1; y++) {
		for (int x = 1; x < img.width - 1; x++) {
			int sector = sectors.map[y][x];
			e.map[y][x] = magnitude.map[y][x];
			if (sector == 0) {
				if (magnitude.map[y-1][x] > magnitude.map[y][x])
					 e.map[y][x] = 0;
				if (magnitude.map[y+1][x] > magnitude.map[y][x])
					 e.map[y][x] = 0;
			}
			if (sector == 1) {
				if (magnitude.map[y-1][x-1] > magnitude.map[y][x])
					 e.map[y][x] = 0;
				if (magnitude.map[y+1][x+1] > magnitude.map[y][x])
					 e.map[y][x] = 0;
			}
			if (sector == 2) {
				if (magnitude.map[y][x-1] > magnitude.map[y][x])
					 e.map[y][x] = 0;
				if (magnitude.map[y][x+1] > magnitude.map[y][x])
					 e.map[y][x] = 0;
			}
			if (sector == 3) {
				if (magnitude.map[y+1][x-1] > magnitude.map[y][x])
					 e.map[y][x] = 0;
				if (magnitude.map[y-1][x+1] > magnitude.map[y][x])
					 e.map[y][x] = 0;
			}
		}
	}

	//Hysteresis thresholds
	for (int y = 1; y < img.height - 1; y++) {
		for (int x = 1; x < img.width - 1; x++) {
			if (e.map[y][x] > highthreshold)
				e.map[y][x] = 255;
			else if (e.map[y][x] < lowthreshold)
				e.map[y][x] = 0;
			else
				e.map[y][x] = -1;
		}
	}

	//Check Edge Candidates
	for (int y = 1; y < img.height - 1; y++) {
		for (int x = 1; x < img.width - 1; x++) {
			if (e.map[y][x] == -1) {
				if (isEightNeighbor(&e, x, y, -1))
					output.map[y][x] = 255;
				else
					output.map[y][x] = 0;
			}
			else
				output.map[y][x] = e.map[y][x];
		}
	}
	return matrix2Image(output,0,1.0);
}

/** Convolve input matrix m1 with filter m2 */
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

/** Return 1 if 8-Neightbor at x,y contains value */
int isEightNeighbor(Matrix *m1, int x, int y, int value) {
	if (m1->map[y-1][x-1] == value)
		return 1;
	if (m1->map[y][x-1] == value)
		return 1;
	if (m1->map[y+1][x-1] == value)
		return 1;
	if (m1->map[y-1][x] == value)
		return 1;
	if (m1->map[y+1][x] == value)
		return 1;
	if (m1->map[y-1][x+1] == value)
		return 1;
	if (m1->map[y][x+1] == value)
		return 1;
	if (m1->map[y+1][x+1] == value)
		return 1;
	return 0;
}


/** Tests */
int testConvolve() {
	double a[4][4] = {
		{0,1,2,3},
		{5,5,6,7},
		{8,9,0,1},
		{2,3,4,5}
	};

	double b[5][5] = {
		{0,1,2,3,4},
		{5,5,6,7,8},
		{8,9,0,1,2},
		{2,3,4,5,6}
	};

	double filtereven[2][2] = {
		{-1,1},
		{0,1}
	};

	double filterodd[3][3] = {
		{1,0,-1},
		{5,2,0},
		{-2,3,4}
	};

	Matrix amatrix = createMatrixFromArray(a,4,4);
	Matrix filterevenmatrix = createMatrixFromArray(filtereven,2,2);
	Matrix filteroddmatrix = createMatrixFromArray(filterodd,3,3);

	if (filterSum(&amatrix,&filterevenmatrix,2,2) == 6.0)
		printf("filterSum test passed!\n");
	else
		printf("filterSum test failed!\n");

	if (filterSum(&amatrix,&filteroddmatrix,2,2) == 69.0)
		printf("filterSum test passed!\n");
	else
		printf("filterSum test failed!\n");

	if (filterSum(&amatrix,&filteroddmatrix,1,1) == 34.0)
		printf("filterSum test passed!\n");
	else
		printf("filterSum test failed!\n");

	// Matrix convolved = convolve(amatrix, filterevenmatrix);
	// printMatrix(amatrix, 0,0,4,4);
	// printMatrix(filterevenmatrix, 0,0,2,2);
	// printMatrix(convolved, 0,0,4,4);

	//printMatrix(amatrix, 0,0,4,4);
	//printMatrix(filterevenmatrix, 0,0,2,2);
	//printf("%f",filterSum(&amatrix,&filteroddmatrix,1,1));
}

/** Debug Function */
void printMatrix(Matrix matrix, int startx, int starty, int width, int height) {
	for (int y = starty; y < starty + height; y++) {
		for (int x = startx; x < startx + width; x++)
			printf("%.*f ",0, matrix.map[y][x]);
		printf("\n");
	}
}
