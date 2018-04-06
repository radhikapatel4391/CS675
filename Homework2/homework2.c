
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "netpbm.h"

Image sobel(Image inputImg);
Image canny(Image inputImg);
void edgeDetection(char *inputFilename, char *sobelFilename, char
                   *cannyFilename);

Matrix convolve(Matrix imageM, Matrix filterM);
Matrix linearScaling(Matrix input,int a,int b,double min,double max);
Matrix gussianSmoothing(Matrix iImage);      //with 5*5 guassian filter with 1.4.
void createGussianFilter(double gKernel[5][5],double sigma);

Matrix sobelfilterXderivative(Matrix forXDiv);
Matrix sobelfilterYderivative(Matrix forYDiv);
Matrix sobelMagnitude(Matrix xMag,Matrix yMag);

Matrix arctan(Matrix xMag, Matrix yMag);
double angle(double x, double y);
Matrix sectorAllocation(Matrix arcTan);

Matrix nonMaximalSupression(Matrix mag, Matrix sectorAll);

Matrix hysteresisThresholding(Matrix nonMax,int lowTh,int hightTh);

void printMatrix(Matrix m,int sx,int sy,int ex, int ey);


void main()
{

    edgeDetection("1.pgm","sobel_output_1.pgm","canny_output_1.pgm");

}
void edgeDetection(char *inputFilename, char *sobelFilename, char
                   *cannyFilename)
{

    Image img;
    img = readImage(inputFilename);
    writeImage(sobel(img),sobelFilename);
    printf("Sobel Done \n");
    writeImage(canny(img),cannyFilename);
    printf("Done \n");
}
Image sobel(Image inputImg)
{
    Matrix inputMatrix,outputMatrix;

    inputMatrix = image2Matrix(inputImg);

    outputMatrix = sobelMagnitude(sobelfilterXderivative(inputMatrix),sobelfilterYderivative(inputMatrix));

    return matrix2Image(outputMatrix,1,1);
}
Image canny(Image inputImg)
{

    printf("Canny: \n\n");
    Matrix oM,xD,yD;

//Step 0 image to Matrix
    oM = image2Matrix(inputImg);
    printf("input\n");
    printMatrix(oM,75,75,100,100);
//Step 1 smoothing
    printf("Step 1 Gaussian Smoothing with sigma 1.4 and window size 5: \n\n");
    oM = gussianSmoothing(oM);
   printf("gussian output \n");
   printMatrix(oM,75,75,100,100);
    writeImage(matrix2Image(oM,1,1), "GussianBlurOutPut_1.pgm");

//Step 2 claculate magnitute and arc and do non maxima suppression
    printf("Step 2 Magnitude Calculation and Alpha: \n\n");
    xD = sobelfilterXderivative(oM);
    yD = sobelfilterYderivative(oM);
   printf("X derivative output \n");
   printMatrix(xD,75,75,100,100);
    printf("y derivative output \n");
   printMatrix(yD,75,75,100,100);
//Step 3
    printf("Step 3 non Maximal Supression: \n\n");
//    printf("arc\n");
//    printMatrix(arctan(xD,yD),75,75,100,100);
//    printf("Sector allocationop \n");
//    printMatrix(sectorAllocation(arctan(xD,yD)),75,75,100,100);
//    printf("magnitude allocationop \n");
//    printMatrix(sobelMagnitude(xD,yD),75,75,100,100);
    oM = nonMaximalSupression(sobelMagnitude(xD,yD), sectorAllocation(arctan(xD,yD)));
    writeImage(matrix2Image(oM,1,1), "MaximaSuppression_atCanny_output_1.pgm");
//    printf("maximal op \n");
//    printMatrix(oM,75,75,100,100);
//step 4 thresholding

    printf("Step 4 hysteresis Thresholding: with threshold 40,100 \n\n");
    oM = hysteresisThresholding(oM,40,60);

    return matrix2Image(oM,1,1);
}
Matrix convolve(Matrix imageM, Matrix filterM)
{

    int imageWidth,imageHeight,fx,fy;
    Matrix outPutMatrix;
    imageWidth = imageM.width;
    imageHeight = imageM.height;

    outPutMatrix = createMatrix(imageHeight,imageWidth);

    fx = floor(filterM.height/2);
    fy = floor(filterM.width/2);


    for (int x = fx; x <(imageHeight-fx); x++)
    {
        for (int y = fy; y <(imageWidth-fy); y++)
        {

            for (int i = -fx; i <= fx; i++)
            {
                for (int j = -fy; j <= fy; j++)
                {

                    outPutMatrix.map[x][y] += filterM.map[i + fx][j + fy] * imageM.map[x + i][y + j];


                }
            }

        }
    }
    return outPutMatrix;


}
Matrix linearScaling(Matrix input,int a,int b,double min,double max)
{
    int h,w;
    h = input.height;
    w = input.width;

    for(int x=0; x<h; x++)
    {
        for(int y = 0; y<w; y++)
        {
            input.map[x][y] = ((b-a)*((input.map[x][y])- min)/(max-min)) + a;

        }

    }

    return input;
}
Matrix sobelfilterXderivative(Matrix m)
{

    Matrix fx;
    int xa[9] = {-1,0,1,-2,0,2,-1,0,1};//array for sobel filter x direction
    int k = 0;
    fx = createMatrix(3,3);//x filter

    //create x direction and y direction filter matrix from array
    for(int x=0; x<3; x++)
    {
        for(int y = 0; y<3; y++)
        {
            fx.map[x][y] = xa[k];
            k++;
        }
    }

    return convolve(m,fx);

}
Matrix sobelfilterYderivative(Matrix m)
{

    Matrix fy;
    int k=0;
    int ya[9] = {-1,-2,-1,0,0,0,1,2,1};//array for sobel filter y direction

    fy = createMatrix(3,3);//y filter

    //create x direction and y direction filter matrix from array
    for(int x=0; x<3; x++)
    {
        for(int y = 0; y<3; y++)
        {

            fy.map[x][y] = ya[k];
            k++;
        }
    }


    return convolve(m,fy);
}
Matrix sobelMagnitude(Matrix xM,Matrix yM)
{
    int h,w;
    h = xM.height;
    w = xM.width;
    Matrix sop = createMatrix(h,w);
    double min=0, max=0;


    for(int x=0; x<h; x++)
    {
        for(int y = 0; y<w; y++)
        {
            double temp = sqrt(pow(xM.map[x][y],2)+ pow(yM.map[x][y],2));
            sop.map [x][y] = temp;
            min = temp < min ? temp: min;
            max = temp < max ? max: temp;
        }

    }

    return linearScaling(sop,0,255,min,max);
    //return sop;
}


Matrix arctan(Matrix xM, Matrix yM)
{
    int h,w;
    h = xM.height;
    w = xM.width;

    for(int x=0; x<h; x++)
    {
        for(int y = 0; y<w; y++)
        {
            xM.map [x][y] = angle(xM.map[x][y],yM.map[x][y]);

        }

    }
    return xM;
}
double angle(double x, double y)
{
    double ang;
    ang = (atan2(y,x))*(180/3.14159265);
    if(ang<0)
    {
        ang +=360;
    }
    return round(ang);
}
Matrix sectorAllocation(Matrix arcTan)
{
    printf("arc....\n");
    printMatrix(arcTan,75,75,100,100);
    int h,w;
    h = arcTan.height;
    w = arcTan.width;
    Matrix sAllOP = createMatrix(h,w);
    for(int x=0; x<h; x++)
    {
        for(int y = 0; y<w; y++)
        {
            double temp = arcTan.map[x][y];
            if((temp > 157.5 && temp < 202.5) || (temp < 22.5 || temp >337.5))
            {
                sAllOP.map[x][y] = 0;

            }
            else if((temp > 22.5 && temp < 67.5) || (temp > 202.5 && temp < 247.5))
            {
                sAllOP.map[x][y] = 1;

            }
            else if((temp > 67.5 && temp < 112.2) || (temp > 247.5 && temp < 292.5))
            {
                sAllOP.map[x][y] = 2;

            }
            else
            {
                sAllOP.map[x][y] = 3;

            }
        }

    }
   printf("sector....\n");
    printMatrix(sAllOP,75,75,100,100);
    return sAllOP;
}


Matrix nonMaximalSupression(Matrix magnitudeMatrix, Matrix sectorAllocateMatrix)
{

    printf("Magnitude....\n");
    printMatrix(magnitudeMatrix,75,75,100,100);

    printf("sector allo....\n");
    printMatrix(sectorAllocateMatrix,75,75,100,100);
    int h,w;
    h = magnitudeMatrix.height;
    w = magnitudeMatrix.width;

    Matrix oM;
    oM = createMatrix(h,w);

    for(int x=1; x<(h-1); x++)
    {
        for(int y = 1; y<(w-1); y++)
        {
            if(sectorAllocateMatrix.map[x][y] == 0)
            {
                if((magnitudeMatrix.map[x][y] < magnitudeMatrix.map[x-1][y]) || (magnitudeMatrix.map[x][y] < magnitudeMatrix.map[x+1][y]))
                {
                    oM.map[x][y] = 0;
                }
                else
                {
                    oM.map[x][y] = magnitudeMatrix.map[x][y];

                }
            }

            else if(sectorAllocateMatrix.map[x][y] == 1)
            {

                if((magnitudeMatrix.map[x][y] < magnitudeMatrix.map[x-1][y-1]) || (magnitudeMatrix.map[x][y] < magnitudeMatrix.map[x+1][y+1]))
                {

                    oM.map[x][y] = 0;
                }
                else
                {

                    oM.map[x][y] = magnitudeMatrix.map[x][y];

                }
            }
            else if(sectorAllocateMatrix.map[x][y] == 2)
            {

                if((magnitudeMatrix.map[x][y] < magnitudeMatrix.map[x][y-1]) || (magnitudeMatrix.map[x][y] < magnitudeMatrix.map[x][y+1]))
                {

                    oM.map[x][y] = 0;
                }
                else
                {

                    oM.map[x][y] = magnitudeMatrix.map[x][y];
                     //  printf("%lf : %lf \n",sectorAllocateMatrix.map[x][y], oM.map[x][y]);

                }

            }
            else
            {

                if((magnitudeMatrix.map[x][y] < magnitudeMatrix.map[x+1][y-1]) || (magnitudeMatrix.map[x][y] < magnitudeMatrix.map[x-1][y+1]))
                {
                    oM.map[x][y] = 0;
                }
                else
                {
                    oM.map[x][y] = magnitudeMatrix.map[x][y];

                }

            }

        }

    }

     printf("Non maximal suppression allo....\n");
    printMatrix(oM,75,75,100,100);
    return oM;
}
Matrix hysteresisThresholding(Matrix m,int lowTh,int hightTh)
{

    Matrix oM;
    int h,w;
    h = m.height;
    w = m.width;
    int count = 0,tcount = (h-2)*(w-2);

    oM = createMatrix(h,w);
    for(int x=0; x<h; x++)
    {
        for(int y = 0; y<w; y++)
        {
            if(m.map[x][y] > hightTh)
            {
                m.map[x][y] = 255;
                oM.map[x][y] = 255;
            }
            else if(m.map[x][y] < lowTh)
            {
                m.map[x][y] = 0;
                oM.map[x][y] = 0;
            }
            else
            {
                oM.map[x][y] = m.map[x][y];
            }
        }
    }
    //second stage of thresholding
    bool flag = true;
    while(flag==true)
    {
        count = 0;
        for(int x=1; x<(h-1); x++)
        {
            for(int y = 1; y<(w-1); y++)
            {
                if(m.map[x][y] > 0 && m.map[x][y] < 255)
                {
                    if((m.map[x-1][y] == 255)||(m.map[x+1][y] == 255)|| (m.map[x][y-1] == 255)|| (m.map[x][y+1] == 255)|| (m.map[x-1][y-1] == 255)|| (m.map[x+1][y+1] == 255)|| (m.map[x-1][y+1] == 255)|| (m.map[x+1][y-1] == 255))
                    {

                        m.map[x][y] = 255;
                    }
                    else
                    {
                        m.map[x][y] = 0;
                    }
                }
                else
                {
                    count ++;

                }
            }
        }

        if (count == tcount )
        {
            flag = false;
        }
    }
    return m;
}

Matrix gussianSmoothing(Matrix m)
{

    Matrix fy;
    int k=0;
    int ya[25] = {2,4,5,4,2,4,9,12,9,4,5,12,15,12,5,4,9,12,9,4,2,4,5,4,2};

    fy = createMatrix(5,5);

    //create x direction and y direction filter matrix from array
    for(int x=0; x<5; x++)
    {
        for(int y = 0; y<5; y++)
        {

            fy.map[x][y] = (ya[k]);
            k++;
        }
    }


    m =  convolve(m,fy);
    //for scaling.....................
    int h,w;
    double min = 0,max=0;
    h = m.height;
    w = m.width;
    for(int x=0; x<h; x++)
    {
        for(int y = 0; y<w; y++)
        {
            double temp = m.map [x][y];
            min = temp < min ? temp: min;
            max = temp < max ? max: temp;
        }

    }

    m = linearScaling(m,0,255,min,max);


    return m;
}
void createGussianFilter(double gKernel[5][5],double sigma)
{
    // set standard deviation to 1.0

    double r, s = 2.0 * sigma * sigma;

    // sum is for normalization
    double sum = 0.0;

    // generate 5x5 kernel
    for (int x = -2; x <= 2; x++)
    {
        for(int y = -2; y <= 2; y++)
        {
            r = (x*x + y*y);
            gKernel[x + 2][y + 2] = (exp(-(r)/s))/(3.14 * s);
            sum += gKernel[x + 2][y + 2];
            // printf("I am sum  %.1lf ",sum);
        }
    }
    // normalize the Kernel
//    for(int i = 0; i < 5; ++i)
//        for(int j = 0; j < 5; ++j){
//            gKernel[i][j] /= sum;
//            gKernel[i][j] = round(gKernel[i][j]);
//        }

}

void printMatrix(Matrix m,int sx,int sy,int ex, int ey)
{
    printf("MAtrix: \n\n");
    for(int x = sx; x<ex; x++)
    {
        for (int y = sy ; y<ey; y++)
        {
            printf(" %0.0f ",m.map[x][y]);
        }
        printf("\n");
    }

}



