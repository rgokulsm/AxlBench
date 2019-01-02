/*
 * kinematics.cpp
 * 
 *  Created on: Sep. 10 2013
 *			Author: Amir Yazdanbakhsh <yazdanbakhsh@wisc.edu>
 */
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include <softinj.hh>

#include <cmath>
#include "kinematics.hpp"
#include</research/sgokul/gem5-stable/gem5-stable//util/m5/m5op.h>


extern int v[92];
extern int t[92];

extern int tv(int v, int t);

extern double fadd_volt_approx[7][9], fmul_volt_approx[7][9], dadd_volt_approx[7][9], dmul_volt_approx[7][9];
extern double gok_ADD(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);
extern double gok_SUB(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);
extern double gok_MUL(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);


float l1 = 0.5 ;
float l2 = 0.5 ;

void forwardk2j(float theta1, float theta2, float* x, float* y) {
        *x = l1 * cos(theta1) + l2 * cos(theta1 + theta2) ;
        *y = l1 * sin(theta1) + l2 * sin(theta1 + theta2) ;

}

void inversek2j(float x, float y, float* theta1, float* theta2) {

	double dataIn[2];
	dataIn[0] = x;
	dataIn[1] = y;

	double dataOut[2];

        float temp1 = gok_MUL(x,x,tv(v[0],t[0]),20,22);
        float temp2 = gok_MUL(y,y,tv(v[1],t[1]),24,26);
        float temp3 = gok_MUL(l1,l1,tv(v[2],t[2]),28,30);
        float temp4 = gok_MUL(l2,l2,tv(v[3],t[3]),32,34);
        float temp5 = gok_MUL(l1,l2,tv(v[4],t[4]),36,38);
        float temp6 = gok_ADD(temp1 , temp2,tv(v[5],t[5]),40,42);
        float temp7 = gok_ADD(temp3 , temp4,tv(v[6],t[6]),44,46);
        float temp8 = gok_SUB(temp6 , temp7,tv(v[7],t[7]),48,50);
        float temp9 = gok_MUL(2.0,temp5,tv(v[8],t[8]),52,54);

        //*theta2 = (float)acos(((x * x) +(y * y) - (l1 * l1) - (l2 * l2))/(2 * (l1 * l2)));
        *theta2 = (float)(temp8/temp9);//acos

        temp1 = (*theta2);//cos
        temp2 = (*theta2);//sin
        temp3 = gok_MUL(l2,temp1,tv(v[11],t[11]),56,58);
        temp4 = gok_MUL(l2,temp2,tv(v[12],t[12]),60,62);
        temp5 = gok_ADD(l1 , temp3,tv(v[13],t[13]),64,66);
        temp6 = gok_MUL(y,temp5,tv(v[14],t[14]),68,70);
        temp7 = gok_MUL(x,temp4,tv(v[15],t[15]),72,74);
        temp8 = gok_SUB(temp5 , temp7,tv(v[16],t[16]),76,78);
        temp9 = gok_MUL(x,x,tv(v[17],t[17]),80,82);
        float temp10 = gok_MUL(y,y,tv(v[18],t[18]),84,86);
        float temp11 = gok_ADD(temp9 , temp10,tv(v[19],t[19]),88,90);

        //*theta1 = (float)asin((( y * ((l1 + l2 * cos(*theta2)))) - (x * (l2 * sin(*theta2))))/((x * x) + (y * y)));
        *theta1 = (float)(temp8/temp11);//asin


	dataOut[0] = (*theta1);
	dataOut[1] = (*theta2);



	*theta1 = dataOut[0];
	*theta2 = dataOut[1];
}
