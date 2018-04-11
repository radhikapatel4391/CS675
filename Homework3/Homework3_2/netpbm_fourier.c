// netpbm_fourier.c
// Useful functions for computing discrete Fourier transforms with the netpbm.c library.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "netpbm.h"

#define MIN(X,Y) ((X)<(Y)?(X):(Y))
#define MAX(X,Y) ((X)>(Y)?(X):(Y))


// Return a new Matrix that is the transpose of the input Matrix.
Matrix transposeMatrix(Matrix inputMatrix)
{
	int i, j;
	Matrix outputMatrix = createMatrix(inputMatrix.width, inputMatrix.height);

	for (i = 0; i < inputMatrix.height; i++)
		for (j = 0; j < inputMatrix.width; j++)
			outputMatrix.map[j][i] = inputMatrix.map[i][j];
	return outputMatrix;
}


// Create a 2D lookup table for a given function; used to accelerate the Fourier transform.
Matrix createLookupTable(int size, double (*func)(double), int offset)
{
	Matrix table = createMatrix(size, size);
	int i, j;

	for (i = 0; i < size; i++)
		for (j = 0; j < size; j++)
			table.map[i][(j + offset)%size] = (*func)(2.0*PI*(double) i*(double) j/(double) size);

	return table;
}


// Apply the complex discrete 1D Fourier transform to each row of the matrix,
// and then transpose the result. Performing this procedure twice results
// in a 2D Fourier transform.
void fourierRow(Matrix inReal, Matrix inImag, Matrix *outReal, Matrix *outImag)
{
	int i, j, k, m = inReal.height, n = inImag.width;
	double sumReal, sumImag;
	Matrix newReal = createMatrix(m, n);
	Matrix newImag = createMatrix(m, n);
	Matrix tableCos = createLookupTable(n, cos, 0);
	Matrix tableSin = createLookupTable(n, sin, 0);

	for (k = 0; k < m; k++)
		for (i = 0; i < n; i++)
		{
			sumReal = sumImag = 0.0;
			for (j = 0; j < n; j++)
			{
				sumReal +=  inReal.map[k][j]*tableCos.map[i][j] + inImag.map[k][j]*tableSin.map[i][j];
				sumImag += -inReal.map[k][j]*tableSin.map[i][j] + inImag.map[k][j]*tableCos.map[i][j];
			}
			newReal.map[k][(i + n/2)%n] = sumReal;
			newImag.map[k][(i + n/2)%n] = sumImag;
		}

	*outReal = transposeMatrix(newReal);
	*outImag = transposeMatrix(newImag);
	deleteMatrix(tableCos);
	deleteMatrix(tableSin);
	deleteMatrix(newReal);
	deleteMatrix(newImag);
}


// Fourier transform the Matrix mxSpace and output magnitude and phase matrices.
void fourier(Matrix mxSpace, Matrix *mxMagni, Matrix *mxPhase)
{
	Matrix mxReal1, mxImag1, mxReal2, mxImag2;
	Matrix mxZero = createMatrix(mxSpace.height, mxSpace.width);
	Matrix mxMag = createMatrix(mxSpace.height, mxSpace.width);
	Matrix mxPha = createMatrix(mxSpace.height, mxSpace.width);
	int i, j;

	fourierRow(mxSpace, mxZero, &mxReal1, &mxImag1);
	fourierRow(mxReal1, mxImag1, &mxReal2, &mxImag2);

	for (i = 0; i < mxSpace.height; i++)
		for (j = 0; j < mxSpace.width; j++)
		{
			mxMag.map[i][j] = sqrt(SQR(mxReal2.map[i][j]) + SQR(mxImag2.map[i][j]));
			mxPha.map[i][j] = atan2(mxImag2.map[i][j], mxReal2.map[i][j]);
		}

	*mxMagni = mxMag;
	*mxPhase = mxPha;

	deleteMatrix(mxReal1);
	deleteMatrix(mxImag1);
	deleteMatrix(mxReal2);
	deleteMatrix(mxImag2);
	deleteMatrix(mxZero);
}


// Apply the complex discrete 1D inverse Fourier transform to each row of the matrix,
// and then transpose the result. Performing this procedure twice in a row results
// in a 2D inverse Fourier transform.
void invFourierRow(Matrix inReal, Matrix inImag, Matrix *outReal, Matrix *outImag)
{
	int i, j, k, m = inReal.height, n = inImag.width;
	double sumReal, sumImag;
	Matrix newReal = createMatrix(m, n);
	Matrix newImag = createMatrix(m, n);
	Matrix tableCos = createLookupTable(n, cos, n/2);
	Matrix tableSin = createLookupTable(n, sin, n/2);

	for (k = 0; k < m; k++)
		for (i = 0; i < n; i++)
		{
			sumReal = sumImag = 0.0;
			for (j = 0; j < n; j++)
			{
				sumReal += inReal.map[k][j]*tableCos.map[i][j] - inImag.map[k][j]*tableSin.map[i][j];
				sumImag += inReal.map[k][j]*tableSin.map[i][j] + inImag.map[k][j]*tableCos.map[i][j];
			}
			newReal.map[k][i] = sumReal/(double) n;
			newImag.map[k][i] = sumImag/(double) n;
		}

	*outReal = transposeMatrix(newReal);
	*outImag = transposeMatrix(newImag);
	deleteMatrix(tableCos);
	deleteMatrix(tableSin);
	deleteMatrix(newReal);
	deleteMatrix(newImag);
}


// Perform an inverse Fourier transform on given magnitude and phase matrices
// and output a matrix containing the data in the spatial domain.
Matrix invFourier(Matrix mxMagni, Matrix mxPhase)
{
	Matrix mxReal1, mxImag1, mxReal2, mxImag2;
	Matrix mxCos = createMatrix(mxMagni.height, mxMagni.width);
	Matrix mxSin = createMatrix(mxMagni.height, mxMagni.width);
	int i, j;

	for (i = 0; i < mxMagni.height; i++)
		for (j = 0; j < mxMagni.width; j++)
		{
			mxCos.map[i][j] = mxMagni.map[i][j]*cos(mxPhase.map[i][j]);
			mxSin.map[i][j] = mxMagni.map[i][j]*sin(mxPhase.map[i][j]);
		}

	invFourierRow(mxCos, mxSin, &mxReal1, &mxImag1);
	invFourierRow(mxReal1, mxImag1, &mxReal2, &mxImag2);

	deleteMatrix(mxReal1);
	deleteMatrix(mxImag1);
	deleteMatrix(mxImag2);
	deleteMatrix(mxCos);
	deleteMatrix(mxSin);
	return mxReal2;
}


// Perform a gamma transform on all entries of Matrix mx so that the new
// maximum value is given by <scale>.
Matrix gammaTransform(Matrix mx, double scale, double gamma)
{
	int m, n;
	double maxVal = -DBL_MAX;
	Matrix result = createMatrix(mx.height, mx.width);

	for (m = 0; m < mx.height; m++)
		for (n = 0; n < mx.width; n++)
			maxVal = MAX(maxVal, mx.map[m][n]);

	if (maxVal < 1e-10)
		maxVal = 1.0;

	for (m = 0; m < mx.height; m++)
		for (n = 0; n < mx.width; n++)
			result.map[m][n] = scale*pow(mx.map[m][n]/maxVal, gamma);

	return result;
}


// Perform a Fourier transform on a given image file and write image files
// containing the magnitude and phase data. Note that the magnitude data are
// gamma transformed with gamma = 0.2 so that their extreme dynamic range
// is reduced and even small values are still visible. Phase data are scaled
// linearly so that values 0 and 255 correspond to -PI and PI.
void fourierDemo(char *inputName, char *outputMagniName, char *outputPhaseName)
{
	Image inputImage = readImage(inputName);
	Matrix inputMatrix = image2Matrix(inputImage);
	Matrix mxMagni, mxPhase;
	Image imMagni, imPhase = createImage(inputImage.height, inputImage.width);
	int i, j;

	fourier(inputMatrix, &mxMagni, &mxPhase);
	imMagni = matrix2Image(mxMagni, 1, 0.2);
	writeImage(imMagni, outputMagniName);

	for (i = 0; i < mxPhase.height; i++)
		for (j = 0; j < mxPhase.width; j++)
			mxPhase.map[i][j] = mxPhase.map[i][j]*127.0/PI + 128.0;

	imPhase = matrix2Image(mxPhase, 0, 0.0);
	writeImage(imPhase, outputPhaseName);

	deleteMatrix(inputMatrix);
	deleteMatrix(mxMagni);
	deleteMatrix(mxPhase);
	deleteImage(inputImage);
	deleteImage(imMagni);
	deleteImage(imPhase);
}


// Perform an inverse Fourier transform on given image files for magnitude and
// phase data and write the resulting spatial image into a file. Note the way
// in which magnitude and phase data are stored as described above.
void invFourierDemo(char *inputMagniName, char *inputPhaseName, char *outputName)
{
	Image magniInImage = readImage(inputMagniName);
	Matrix magniTemp = image2Matrix(magniInImage);
	Matrix magniInMatrix = gammaTransform(magniTemp, 1000.0, 5.0);
	Image phaseInImage = readImage(inputPhaseName);
	Matrix phaseTemp = image2Matrix(phaseInImage);
	Matrix phaseInMatrix = createMatrix(phaseTemp.height, phaseTemp.width);
	Matrix mxOutput;
	Image imOutput;
	int i, j;

	for (i = 0; i < phaseTemp.height; i++)
		for (j = 0; j < phaseTemp.width; j++)
			phaseInMatrix.map[i][j] = (phaseTemp.map[i][j] - 128.0)/127.0*PI;

	mxOutput = invFourier(magniInMatrix, phaseInMatrix);
	imOutput = matrix2Image(mxOutput, 1, 1.0);
	writeImage(imOutput, outputName);

	deleteMatrix(magniTemp);
	deleteMatrix(magniInMatrix);
	deleteMatrix(phaseTemp);
	deleteMatrix(phaseInMatrix);
	deleteMatrix(mxOutput);
	deleteImage(magniInImage);
	deleteImage(phaseInImage);
	deleteImage(imOutput);
}


void main()
{
	//fourierDemo("umb_noisy.pgm", "umb_magni.pgm", "umb_phase.pgm");
	invFourierDemo("umb_magni2.pgm", "umb_phase.pgm", "umb_restored.pgm");
}
