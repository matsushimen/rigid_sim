//
//  myfunc.h
//  
//
//  Created by 松島　佑樹 on 2016/06/16.
//
//

#ifndef _myfunc_h
#define _myfunc_h


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mol_def.h"

double norm(double x, double y, double z)//2点間の距離
{
    return(sqrt(x*x+y*y+z*z));
}

double dist(int i, int j)//2粒子の中心の距離
{
    return(norm(mlcl[i].p[X] - mlcl[j].p[X],mlcl[i].p[Y] - mlcl[j].p[Y],mlcl[i].p[Z] - mlcl[j].p[Z]));
}

double EVEx(int i, int j, double sc)//x方向の排除体積効果
{
    double f=0.0,d,R;
    R = mlcl[i].r + mlcl[j].r;//バネの自然長
    d =dist(i,j);
    if((d>0.0)&& (d < (mlcl[i].r + mlcl[j].r)))
    {
        f = -sc*(mlcl[i].p[X] - mlcl[j].p[X])*(d - R)/d;
    }
    return(f);
}

double EVEy(int i, int j, double sc)//y方向の排除体積効果
{
    double f=0.0,d,R;
    R = mlcl[i].r + mlcl[j].r;//バネの自然長
    d =dist(i,j);
    if((d>0.0)&& (d < (mlcl[i].r + mlcl[j].r)))
    {
        f = -sc*(mlcl[i].p[Y] - mlcl[j].p[Y])*(d - R)/d;
    }
    return(f);
}

double EVEz(int i, int j, double sc)//z方向の排除体積効果
{
    double f=0.0,d,R;
    R = mlcl[i].r + mlcl[j].r;//バネの自然長
    d =dist(i,j);
    if((d>0.0)&& (d < (mlcl[i].r + mlcl[j].r)))
    {
        f = -sc*(mlcl[i].p[Z] - mlcl[j].p[Z])*(d - R)/d;
    }
    return(f);
}


double LSx(int i, int j, double sc)//線形バネで接続
{
    double R, d;
    R = mlcl[i].r + mlcl[j].r;//バネの自然長
    R = R * 0.9;
    d =dist(i,j);

    return(-sc*(mlcl[i].p[X] - mlcl[j].p[X])*(d - R)/d);
}
double LSy(int i, int j, double sc)
{
    double R, d;
    R = mlcl[i].r + mlcl[j].r;//バネの自然長
    R = R * 0.9;
    d = dist(i,j);

    return(-sc*(mlcl[i].p[Y] - mlcl[j].p[Y])*(d - R)/d);
}
double LSz(int i, int j, double sc)
{
    double R, d;
    R = mlcl[i].r + mlcl[j].r;//バネの自然長
    R = R * 0.9;
    d = dist(i,j);

    return(-sc*(mlcl[i].p[Z] - mlcl[j].p[Z])*(d - R)/d);
}

double LS2x(int i, int j, double sc)//２個隣と線形バネで接続
{
    double R,d;
    R =  4.0;
    d = dist(i,j);
    return(-sc*(mlcl[i].p[X] - mlcl[j].p[X])*(d - R)/d);
}
double LS2y(int i, int j, double sc)
{
    double R,d;
    R =  4.0;
    d = dist(i,j);
    return(-sc*(mlcl[i].p[Y] - mlcl[j].p[Y])*(d - R)/d);
}
double LS2z(int i, int j, double sc)
{
    double R,d;
    R =  4.0;
    d = dist(i,j);
    return(-sc*(mlcl[i].p[Z] - mlcl[j].p[Z])*(d - R)/d);
}

#endif
