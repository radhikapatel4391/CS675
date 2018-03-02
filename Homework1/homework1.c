//Remove salt and paper noise and find connected component for high resolution text image.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "netpbm.h"

int random_number(int min_num, int max_num);
Image expandingORshrinking (Image inputImage,int eORs);
void imageCleaning(Image inputImage);
Image labelingToConnectedComponnent(Image inputImage,int threshold);

void main()
{
    Image inputImage, outputImage, cleanImage, temp;
    Matrix inputMatrix, outputMatrix;
    printf("Just keep patience...\n");
    printf("It will take time ...\n");
    inputImage = readImage("text_image.pbm");
    imageCleaning(inputImage); //Q1 answer.

    cleanImage = readImage("text_cleaned.pbm");
    outputImage = labelingToConnectedComponnent(cleanImage,150);//Q-3 answer. threshold set to 150.
    writeImage(outputImage, "text_colored.ppm");
    //ouputImage = matrix2Image(outputMatrix,1,1);
    //writeImage(ouputImage, "matrixtoImage.pgm");
}


//function for expanding(255) and shrinking(0) based on 4NN
    Image expandingORshrinking (Image inputImage,int eORs)
    {
        //printf("For 0 do shrinking for 255 do expanding \n");
        int x,y;
        int temp1 = eORs;
        int temp2 = 255 - eORs;
        Image outputImage;
        Matrix inputMatrix, outputMatrix;
        inputMatrix = image2Matrix(inputImage);
        outputMatrix = image2Matrix(inputImage);
        //Look 4-NN neighbor value and based on that change the current value if any neighbor has different value then current.
        //Modified value doesn't count for the next discussion.
        for (int x = 1; x < (inputImage.width - 1); x++)
            for (int y = 1; y < (inputImage.height - 1); y++)
            {
                if(inputMatrix.map[y][x] == temp1)
                {
                    if(inputMatrix.map[y-1][x] == temp2 || inputMatrix.map[y+1][x] == temp2 || inputMatrix.map[y][x-1] == temp2 || inputMatrix.map[y][x+1] == temp2 )
                    {
                        outputMatrix.map[y][x] = temp2;
                    }

                }
            }
        //set border value to background value.
        for (int x = 0; x<1 ; x++)
            for (int y = 0; y < inputImage.height; y++)
            {
                outputMatrix.map[y][x] = 255;
            }
        for (x = 0; x <inputImage.width; x++)
            for (y = 0; y <1; y++)
            {
                outputMatrix.map[y][x] = 255;
            }
        for (x = inputImage.width-1; x <inputImage.width; x++)
            for (y = 0; y < inputImage.height; y++)
            {
                outputMatrix.map[y][x] = 255;
            }
        for (x = 0; x <inputImage.width; x++)
            for (y = inputImage.height-1; y <inputImage.height; y++)
            {
                outputMatrix.map[y][x] = 255;
            }
        outputImage = matrix2Image(outputMatrix,1,1);
        return outputImage;
    } //end of function

//function for cleaning Image Sequences of shrinking nd expanding
    void imageCleaning(Image inputImage)
    {
        Image ouputImage;
        ouputImage = expandingORshrinking(inputImage, 255);
       // writeImage(ouputImage, "f1.pbm");
        ouputImage = expandingORshrinking(ouputImage, 0);
       // writeImage(ouputImage, "f2.pbm");
        ouputImage = expandingORshrinking(ouputImage, 0);
       // writeImage(ouputImage, "f3.pbm");
        ouputImage = expandingORshrinking(ouputImage, 255);
       // writeImage(ouputImage, "f4.pbm");
        ouputImage = expandingORshrinking(ouputImage, 255);
       // writeImage(ouputImage, "f5.pbm");
        ouputImage = expandingORshrinking(ouputImage, 0);
        writeImage(ouputImage, "text_cleaned.pbm");
        //Other approach with less operation work more or less same.
        ouputImage = expandingORshrinking(inputImage, 0);
       // writeImage(ouputImage, "s1.pbm");
        ouputImage = expandingORshrinking(ouputImage, 255);
       // writeImage(ouputImage, "s2.pbm");
        ouputImage = expandingORshrinking(ouputImage, 255);
       // writeImage(ouputImage, "s3.pbm");
        ouputImage = expandingORshrinking(ouputImage, 0);
        writeImage(ouputImage, "secondText_clean.pbm");
    }

//function for labelingToComponent for 4NN connected components.

    Image labelingToConnectedComponnent(Image inputImage,int threshold)
    {

        writeImage(inputImage,"temp.ppm");
        Image temp = readImage("temp.ppm");
        Matrix inputMatrix = createMatrix(inputImage.height,inputImage.width);
        int label = 1;
        int totalText = 0;
        int x,y;
//traversing on image.
        for (x = 1; x < (inputMatrix.width - 1); x++)
        {
            for (y = 1; y < (inputMatrix.height - 1); y++)
            {
                //is it text pixcel?
                if(inputImage.map[y][x].i == 0)
                {
                    // are there both neighbor have label?
                    if(inputMatrix.map[y-1][x]!= 0 && inputMatrix.map[y][x-1]!=0)
                    {
                        // are there both neighbor contain same label?
                        if(inputMatrix.map[y-1][x] == inputMatrix.map[y][x-1])
                        {
                            //assign any one.
                            inputMatrix.map[y][x] = inputMatrix.map[y-1][x];
                        }
                        else if(inputMatrix.map[y-1][x] < inputMatrix.map[y][x-1])
                        {
                            // upper label has small value so we will use that for current value and will change large value to low value  every where also.
                            //instead of Equivalence table when i found both same i simply change it.
                            inputMatrix.map[y][x] = inputMatrix.map[y-1][x];
                            int mLabel = inputMatrix.map[y][x-1];
                            int nLabel = inputMatrix.map[y-1][x];
                            for(int xi =0; xi<(x+1); xi++)
                                for(int yi = 0; yi<(y+1); yi++)
                                {
                                    if(inputMatrix.map[yi][xi] == mLabel)
                                    {
                                        inputMatrix.map[yi][xi] = nLabel;
                                    }
                                }
                        }
                        else
                        {
                            // second is small so same as above.
                            inputMatrix.map[y][x] = inputMatrix.map[y][x-1];
                            int mLabel = inputMatrix.map[y-1][x];
                            int nLabel = inputMatrix.map[y][x-1];
                            for(int xi =0; xi<(x+1); xi++)
                            {
                                for(int yi = 0; yi<(y+1); yi++)
                                {
                                    if(inputMatrix.map[yi][xi] == mLabel)
                                    {

                                        inputMatrix.map[yi][xi] = nLabel;
                                    }
                                }
                            }

                        }
                    }
                    else if(inputMatrix.map[y-1][x]!= 0 || inputMatrix.map[y][x-1]!=0)
                    {
                        //only one has label and we will use it for current.
                        if((inputMatrix.map[y-1][x]!= 0))
                        {
                            inputMatrix.map[y][x] = inputMatrix.map[y-1][x];
                        }
                        else
                        {
                            inputMatrix.map[y][x] = inputMatrix.map[y][x-1];
                        }
                    }
                    else
                    {
                        //no one has label so let's assign new one.
                        inputMatrix.map[y][x] = label;
                        label++;
                    }
                }
            }
        }//finish the labeling task
 //let's start remove elements which are smaller in size.
 //to store each label's size took array size.
        int size[label];
//for each label value let's traverse Matrix.
        for(int l = 0; l <label; l++)
        {
            size[l] = 0;
            for (x = 1; x < (inputMatrix.width - 1); x++)
            {
                for (y = 1; y < (inputMatrix.height - 1); y++)
                {
                    if(inputMatrix.map[y][x] == (l+1))
                    {
                        size[l]++;
                    }
                }
            }
    //if size for specific label is smaller then threshold then simply discard(assign 0 value as label) it's label.
            if(size[l] < threshold)
            {
                for (x = 1; x < (inputMatrix.width - 1); x++)
                {
                    for (y = 1; y < (inputMatrix.height - 1); y++)
                    {
                        if(inputMatrix.map[y][x] == (l+1))
                        {
                            temp.map[y][x].r = 0;
                            temp.map[y][x].g = 0;
                            temp.map[y][x].b = 0;
                            inputMatrix.map[y][x] = 0;
                        }
                    }
                }
            }
    //if number of pixels for label is higher then make it's label to specific number.
            else
            {
                totalText++;
                int t = random_number(1,254);
                double ra = t;
                double rg = 255-t;
                double rb = totalText;

                for (x = 1; x < (inputMatrix.width - 1); x++)
                {
                    for (y = 1; y < (inputMatrix.height - 1); y++)
                    {
                        if(inputMatrix.map[y][x] == (l+1))
                        {
                           //assign green to all component which have higher size then threshold.
                            temp.map[y][x].r = ra;
                            temp.map[y][x].g = rg;
                            temp.map[y][x].b = rb;
                           //temp.map[y][x].i = totalText*1.4;
                            inputMatrix.map[y][x] = totalText*1.3;

                        }
                    }
                }
            }
        }

        printf("Total texts(Connected components) in image except , and . are  %d \n",totalText);
        return temp;
    }


int random_number(int min_num, int max_num)
{
    int result = 0, low_num = 0, hi_num = 0;

    if (min_num < max_num)
    {
        low_num = min_num;
        hi_num = max_num + 1; // include max_num in output
    } else {
        low_num = max_num + 1; // include max_num in output
        hi_num = min_num;
    }

    srand(time(NULL));
    result = (rand() % (hi_num - low_num)) + low_num;
    return result;
}



