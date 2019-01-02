#include</research/sgokul/gem5-stable/gem5-stable//util/m5/m5op.h>
#include "fourier.hpp"
#include <cmath>
#include <fstream>
#include <map>


extern int v[92];
extern int t[92];

extern int tv(int v, int t);

extern double fadd_volt_approx[7][9], fmul_volt_approx[7][9], dadd_volt_approx[7][9], dmul_volt_approx[7][9];
extern double gok_ADD(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);
extern double gok_SUB(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);
extern double gok_MUL(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);

void calcFftIndices(int K, int* indices)
{
	int i, j ;
	int N ;

	N = (int)log2f(K) ;

	indices[0] = 0 ;
	indices[1 << 0] = 1 << (N - (0 + 1)) ;
	for (i = 1; i < N; ++i)
	{
		for(j = (1 << i) ; j < (1 << (i + 1)); ++j)
		{
			indices[j] = indices[j - (1 << i)] + (1 << (N - (i + 1))) ;
		}
	}
}

void radix2DitCooleyTykeyFft(int K, int* indices, Complex* x, Complex* f)
{

	calcFftIndices(K, indices) ;

	int step ;
	double arg ;
	int eI ;
	int oI ;

	double fftSin;
	double fftCos;

	Complex bunt;
	int i ;
	int N ;
	int j ;
	int k ;

	double dataIn[1];
	double dataOut[2];

	for(i = 0, N = 1 << (i + 1); N <= K ; i++, N = 1 << (i + 1))
	{
		for(j = 0 ; j < K ; j += N)
		{
			step = N >> 1 ;
			for(k = 0; k < step ; k++)
			{
				arg = (double)k / N ;
				eI = j + k ; 
				oI = j + step + k ;

				dataIn[0] = arg;

//#pragma parrot(input, "fft", [1]dataIn)
//m5_start_approx();
//m5_end_approx(); //FIXME: Switch between these 2 - needed for serialization
				fftSinCos(arg, &fftSin, &fftCos);
				//fftSin = arg;
				//double temp0 = gok_MUL(arg,arg,tv(v[2],t[2]));	
				//fftCos = gok_SUB(1.0,temp0,tv(v[3],t[3]));

				dataOut[0] = fftSin;
				dataOut[1] = fftCos;

//#pragma parrot(output, "fft", [2]<0.0; 2.0>dataOut)
//m5_end_approx();
				fftSin = dataOut[0];
				fftCos = dataOut[1];


				// Non-approximate
				bunt =  x[indices[eI]] ;
				x[indices[eI]].real = gok_ADD(bunt.real , gok_SUB(gok_MUL(x[indices[oI]].real , fftCos,tv(v[4],t[4]),28,30) , gok_MUL(x[indices[oI]].imag , fftSin,tv(v[5],t[5]),32,34), tv(v[6],t[6]),36,38), tv(v[7],t[7]),40,42);
                		x[indices[eI]].imag = gok_ADD(bunt.imag , gok_ADD(gok_MUL(x[indices[oI]].imag , fftCos,tv(v[8],t[8]),44,46) , gok_MUL(x[indices[oI]].real , fftSin,tv(v[9],t[9]),48,50), tv(v[10], t[10]),52,54), tv(v[11], t[11]),56,58);

                		x[indices[oI]].real = gok_SUB(bunt.real , gok_SUB(gok_MUL(x[indices[oI]].real , fftCos,tv(v[12],t[12]),60,62) , gok_MUL(x[indices[oI]].imag , fftSin,tv(v[13],t[13]),64,66), tv(v[14], t[14]),68,70), tv(v[15], t[15]),72,74);
                		x[indices[oI]].imag = gok_SUB(bunt.imag , gok_ADD(gok_MUL(x[indices[oI]].imag , fftCos,tv(v[16],t[16]),76,78) , gok_MUL(x[indices[oI]].real , fftSin,tv(v[17],t[17]),80,82), tv(v[18], t[18]),84,86), tv(v[19], t[19]),88,90);
			}
		}
	}

	for (int i = 0 ; i < K ; i++)
	{
		f[i] = x[indices[i]] ;
	}
}
