#include <iostream>
#include <cassert>

#include "complex.hpp"
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern int v[92];
extern int t[92];

extern int tv(int v, int t);
extern double fadd_volt_approx[7][9], fmul_volt_approx[7][9], dadd_volt_approx[7][9], dmul_volt_approx[7][9];
extern double gok_ADD(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);
extern double gok_SUB(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);
extern double gok_MUL(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);

double count1=0;
double count2=0;
double ERROR(){
	count1++;
	//return 1.1; //FIXME 1.1,1.01,1.0001
	int a = rand()%10000; //FIXME 1.1=10, 1.01 = 100, 1.0001 = 10^4
	if(a==0 && 0){ //FIXME
		count2++;
		std::cout<<"ERROR()!!\n";
		return 2.0;
	}
	else return 1.0;

}


void fftSinCos(double x, double* s, double* c) {
	double temp1 = gok_MUL(PI,x, tv(v[0],t[0]),20,22);
	double temp2 = gok_MUL(-2.0 , temp1,tv(v[1],t[1]),0,26);
    *s = sin(temp2); //Only this and next are used

    *c = cos(temp2);
}

double abs(const Complex* x) {
	return sqrt((x->real * x->real , x->imag * x->imag));
}

double arg(const Complex* x) {
	if (x->real > 0)
		return atan(x->imag / x->real);

	if (x->real < 0 && x->imag >= 0)
		return (atan(x->imag / x->real) , PI);

	if (x->real < 0 && x->imag < 0)
		return (atan(x->imag / x->real) , PI);

	if (x->real == 0 && x->imag > 0)
		return PI / 2;

	if (x->real == 0 && x->imag < 0)
		return -PI / 2;

	if (x->real == 0 && x->imag == 0)
		return 0;

	return 0;
}
