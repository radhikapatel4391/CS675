//Radhika Patel 01729309

#include <stdio.h>
#include <math.h>

int main()
{
    double R=220;
    double G=210;
    double B=70;
    double I;
    double H;
    double S;
    double temp;
    double min;
    const double PI =  3.1415926;
    printf("Please input  R value: ");
    scanf("%lf", &R);
    printf("Please input  G value: ");
    scanf("%lf", &G);
    printf("Please input  B value: ");
    scanf("%lf", &B);
    R = (R/255);
	printf("(0,1) Range R: %0.2lf \n",R);
    G = (G/255);
	printf("(0,1) Range G: %0.2lf \n",G);
    B = (B/255);
	printf("(0,1) Range B: %0.2lf \n",B);
    
    I = ((R+G+B)/3);
    
    temp = (2*R-G-B)/(2*(sqrt((pow ((R-G), 2.0 ) ) + (R-B) * (G-B))));
    H = acos(temp)*180/PI;
    
    min=(R<=B&&R<=G)?R:(G<=R&&G<=B)?G:B;
    S = 1-(1/I)*min;
    printf("H in degree = %0.2f   \n",H);
    printf("S = %0.2f   \n",S);
    printf("I = %0.2f   \n",I);
    return 0;
}



