#ifndef __COMPLEX_HPP__
#define __COMPLEX_HPP__

#define PI 3.1415926535897931

typedef struct {
   double real;
   double imag;
} Complex;

void fftSinCos(double x, double* s, double* c);
double abs(const Complex* x);
double arg(const Complex* x);

#endif
