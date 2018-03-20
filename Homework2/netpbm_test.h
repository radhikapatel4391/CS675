// netpbm_test.c
// Test and demo program for netpbm. Reads a sample image and creates several output images.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "netpbm.h"


void main()
{
	printf("Welcome");
	Image inputImage, rotatedImage, invertedImage;
	int x, y;

	inputImage = readImage("sample.ppm");
	rotatedImage = createImage(inputImage.height, inputImage.width);
	invertedImage = createImage(inputImage.height, inputImage.width);

	for (x = 0; x < inputImage.width; x++)
		for (y = 0; y < inputImage.height; y++)
		{
			rotatedImage.map[y][x] = inputImage.map[inputImage.height - y - 1][inputImage.width - x - 1];
			invertedImage.map[y][x].r = 255 - inputImage.map[y][x].r;
			invertedImage.map[y][x].g = 255 - inputImage.map[y][x].g;
			invertedImage.map[y][x].b = 255 - inputImage.map[y][x].b;
			// Let's just ignore 'i' here; it's irrelevant if we want to save image as PPM.
		}
	writeImage(rotatedImage, "rotated.ppm");
	writeImage(invertedImage, "inverted.ppm");
	writeImage(inputImage, "gray.pgm");
	writeImage(inputImage, "black-white.pbm");
	deleteImage(inputImage);
	deleteImage(rotatedImage);
	deleteImage(invertedImage);
}

