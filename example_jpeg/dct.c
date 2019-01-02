/*
 * dct.c
 * 
 * Created on: Sep 9, 2013
 * 			Author: Amir Yazdanbakhsh <a.yazdanbakhsh@gatech.edu>
 */

#include "datatype.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stdio.h"
#include "jpegconfig.h"
float count1=0;
float count2=0;
float ERROR(){
	count1++;
	//return 1.1; //FIXME 1.1,1.01,1.0001
	int a = rand()%10000; //FIXME 1.1=10, 1.01 = 100, 1.0001 = 10^4
	if(a==0 && 0){
		count2++;
		std::cout<<"ERROR()!!\n";
		return 2.0;
	}
	else return 1.0;

}

/* Level shifting to get 8 bit SIGNED values for the data  */
void levelShift (INT16* const data)
{
	INT16 i;

	for (i=63; i>=0; i--)
		data [i] -= 128;
}

float C(int u) {
	if (u == 0)
		return (1.0 / sqrt(2.0));
	else
		return 1.0;
}

/* DCT for One block(8x8) */
void dct (INT16 *data)
{

	UINT16 i;
	INT32 x0, x1, x2, x3, x4, x5, x6, x7, x8;

/*	All values are shifted left by 10
	and rounded off to nearest integer */

	static const UINT16 c1=1420;	/* cos PI/16 * root(2)	*/
	static const UINT16 c2=1338;	/* cos PI/8 * root(2)	*/
	static const UINT16 c3=1204;	/* cos 3PI/16 * root(2)	*/
	static const UINT16 c5=805;		/* cos 5PI/16 * root(2)	*/
	static const UINT16 c6=554;		/* cos 3PI/8 * root(2)	*/
	static const UINT16 c7=283;		/* cos 7PI/16 * root(2)	*/

	static const UINT16 s1=3;
	static const UINT16 s2=10;
	static const UINT16 s3=13;

	for (i=8; i>0; i--)
	{
		x8 = ERROR() * (data [0] + data [7]);
		x0 = ERROR() * (data [0] - data [7]);

		x7 = ERROR() * (data [1] + data [6]);
		x1 = ERROR() * (data [1] - data [6]);

		x6 = ERROR() * (data [2] + data [5]);
		x2 = ERROR() * (data [2] - data [5]);

		x5 = ERROR() * (data [3] + data [4]);
		x3 = ERROR() * (data [3] - data [4]);

		x4 = ERROR() * (x8 + x5);
		x8 -= (ERROR() * x5);

		x5 = ERROR() * (x7 + x6);
		x7 -= (ERROR() * x6);

		data [0] = (INT16) (ERROR() * (x4 + x5));
		data [4] = (INT16) (ERROR() * (x4 - x5));

		data [2] = (INT16) (((1*x8*c2 + 1*x7*c6)) >> s2);
		data [6] = (INT16) (((1*x8*c6 - 1*x7*c2)) >> s2);

		data [7] = (INT16) (((1*1*x0*c7 - 1*1*x1*c5 + 1*1*x2*c3 - 1*1*x3*c1)) >> s2);
		data [5] = (INT16) (((1*1*x0*c5 - 1*1*x1*c1 + 1*1*x2*c7 + 1*1*x3*c3)) >> s2);
		data [3] = (INT16) (((1*1*x0*c3 - 1*1*x1*c7 - 1*1*x2*c1 - 1*1*x3*c5)) >> s2);
		data [1] = (INT16) (((1*1*x0*c1 + 1*1*x1*c3 + 1*1*x2*c5 + 1*1*x3*c7)) >> s2);

		data += 8;
	}

	data -= 64;

	for (i=8; i>0; i--)
	{
		x8 = ERROR()*(data [0] + data [56]);
		x0 = ERROR()*(data [0] - data [56]);

		x7 = ERROR()*(data [8] + data [48]);
		x1 = ERROR()*(data [8] - data [48]);

		x6 = ERROR()*(data [16] + data [40]);
		x2 = ERROR()*(data [16] - data [40]);

		x5 = ERROR()*(data [24] + data [32]);
		x3 = ERROR()*(data [24] - data [32]);

		x4 = ERROR()*(x8 + x5);
		x8 -= (ERROR()*x5);

		x5 = ERROR()*(x7 + x6);
		x7 -= (ERROR()*x6);

		data [0] = (INT16) (((1*x4 + 1*x5)) >> s1);
		data [32] = (INT16) (((1*x4 - 1*x5)) >> s1);

		data [16] = (INT16) (((1*1*x8*c2 + 1*1*x7*c6)) >> s3);
		data [48] = (INT16) (((1*1*x8*c6 - 1*1*x7*c2)) >> s3);

		data [56] = (INT16) (((1*1*x0*c7 - 1*1*x1*c5 + 1*1*x2*c3 - 1*1*x3*c1)) >> s3);
		data [40] = (INT16) (((1*1*x0*c5 - 1*1*x1*c1 + 1*1*x2*c7 + 1*1*x3*c3)) >> s3);
		data [24] = (INT16) (((1*1*x0*c3 - 1*1*x1*c7 - 1*1*x2*c1 - 1*1*x3*c5)) >> s3);
		data [8] = (INT16) (((1*1*x0*c1 + 1*1*x1*c3 + 1*1*x2*c5 + 1*1*x3*c7)) >> s3);

		data++;
	}
}
