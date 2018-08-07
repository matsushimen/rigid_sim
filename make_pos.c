#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI M_PI
#define R 1.0
#define RR 30.0
#define large 50
#define small 7
int i,j,J = 0;
double x[large+small+large+small], y[large+small+large+small], z[large+small+large+small],r[large+small+large+small], scale[large+small+large+small];

double dist(double x,double y,double z);
double radius(int i);
int main(void)
{
	printf("%d\n",large+small+large+small);
    printf("%d\n", large);
    printf("%d\n", small);
    printf("%d\n", 0);//step
    printf("%d\n", 0);//save
    for(i = 0; i < large; i++)//0-10の11こ
    {
        x[i] = RR*cos(i*PI/(large-1)+PI/2.0);
        y[i] = RR*sin(i*PI/(large-1)+PI/2.0);
        z[i] = 0.0;
        r[i] = R;
        scale[i] = 1.0;
        
    }
    for(i = 0; i < small; i++)//11-131の121こ
    {
        x[i+large] = x[large-1]  + 2*(R) + (R) * 2 * i ;
        y[i+large] = y[large-1];
        z[i+large] = 0.0;
        r[i+large] = R;
        scale[i+large] = 0.0;
    }
    for(i = 0; i < large; i++)//132-142の11こ
    {
        x[i+large+small] = RR*cos(i*PI/(large-1)+3*PI/2.0) + x[large+small-1] + R + R;
        y[i+large+small] = RR*sin(i*PI/(large-1)+3*PI/2.0);
        z[i+large+small] = 0.0;
        r[i+large+small] = R;
        scale[i+large+small] = 1.0;
    }
    for(i = 0; i < small; i++)//143-263の121こ
    {
        x[i+large+small+large] = x[large+small+large-1] - R - R - R * 2 * i;
        y[i+large+small+large] = y[large+small+large-1];
        z[i+large+small+large] = 0.0;
        r[i+large+small+large] = R;
        scale[i+large+small+large] = 0.0;
    }
    
    for(i = 0; i < large+small+large+small; i++)
    {
        printf("%d %lf %lf %lf %lf %lf %lf\n", i, x[i], y[i], z[i],r[i], 1.0, scale[i]);//番号　x座標　y座標　z座標　半径　質量　サイズ
    }
}

double dist(double x,double y,double z)
{
	return(sqrt(x*x+y*y+z*z));
}


