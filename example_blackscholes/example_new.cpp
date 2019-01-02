// Copyright (c) 2007 Intel Corp.

// Black-Scholes
// Analytical method for calculating European Options
//
// 
// Reference Source: Options, Futures, and Other Derivatives, 3rd Edition, Prentice 
// Hall, John C. Hull,

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
//#include</research/sgokul/gem5-stable/gem5-stable//util/m5/m5op.h>
#include <time.h>
//double max_otype, min_otype ;
//double max_sptprice, min_sptprice;
//double max_strike, min_strike;
//double max_rate, min_rate ;
//double max_volatility, min_volatility;
//double max_otime, min_otime ;
//double max_out_price, min_out_price;
#include <unistd.h>

#include <cassert>
//#include <softinj.hh>


const int variables = 180; //Should become ~200
int v[variables];
int t[variables];
int pc[variables]={
68006,//0
68238,//1
68294,//2
68350,//3
68374,//4
68422,//5
68456,//6
68724,//7
68772,//8
68836,//9
69010,//10
69100,//11
0,//12
0,//13
0,//14
69036,//15
69060,//16
0,//17
67946,//18
67978,//19
68046,//20
68070,//21
68094,//22
68118,//23
68146,//24
68178,//25
68210,//26
68266,//27
68322,//28
68398,//29
68676,//30
68700,//31
68748,//32
68796,//33
68962,//34
68986,//35

68006-8,68006-8,0,0,//36
68238-8,68238-8,68238-4,68238-4,//40
68294-8,68294-8,68294-4,68294-4,//44
68350-8,68350-8,68350-4,68350-4,//48
68374-8,68374-8,68374-4,68374-4,//52
0,0,68422-4,68422-4,//56
0,0,68456-4,68456-4,//60
68724-8,68724-8,68724-4,68724-4,//64
68772-8,68772-8,68772-4,68772-4,//68
68836-8,68836-8,68836-4,68836-4,//72
69010-8,69010-8,69010-4,69010-4,//76
0,0,0,0,//80
0,0,0,0,//84
0,0,0,0,//88
0,0,0,0,//92
0,0,69036-4,69036-4,//96
0,0,69060-4,69060-4,//100
0,0,0,0,//104
0,0,67946-4,67946-4,//108
0,0,67978-4,67978-4,//112
68046-8,68046-8,68046-4,68046-4,//116
68070-8,68070-8,68070-4,68070-4,//120
68094-8,68094-8,68094-4,68094-4,//124
68118-8,68118-8,68118-4,68118-4,//128
0,0,68146-4,68146-4,//132
0,0,68178-4,68178-4,//136
0,0,68210-4,68210-4,//140
0,0,68266-4,68266-4,//144
0,0,68322-4,68322-4,//148
68398-8,68398-8,68398-4,68398-4,//152
68676-8,68676-8,68676-4,68676-4,//156
68700-8,68700-8,0,0,//160
68748-8,68748-8,68748-4,68748-4,//164
68796-8,68796-8,68796-4,68796-4,//168
68962-8,68962-8,68962-4,68962-4,//172
68986-8,68986-8,68986-4,68986-4//176

};//LVAshouldappeartwice,startsfrom67998

int tv(int v, int t){ //Voltage based on threshold
    int rand_val = rand()%100;
    int V4Ven;
    if(rand_val < t){
	V4Ven = v+1;
    }
    else {
    	V4Ven = v;
    }
    if(V4Ven > 12) V4Ven = 12;
return V4Ven;
}




double approx_table[200][4];//Not more than 100 approx loads


void clear_approx_table()
{
	for(int i =0; i<200; i++)
		for(int j=0; j <4; j++)
			approx_table[i][j]=-123;

}


double read_approx_table(int val)
{
//	printf("read_approx_table: val = %d, [0] = %lf, [1] = %lf, [2] = %lf, [3] = %lf, out=%lf\n",val,approx_table[val][0],approx_table[val][1],approx_table[val][2],approx_table[val][3],((approx_table[val][0]+approx_table[val][1]+approx_table[val][2]+approx_table[val][3])/4.0));
	return ((approx_table[val][0]+approx_table[val][1]+approx_table[val][2]+approx_table[val][3])/4.0); //Average
}

void write_approx_table(double write_val, int val)
{

//	printf("write_approx_table: val = %d, (BEFORE) [0] = %lf, [1] = %lf, [2] = %lf, [3] = %lf ",val,approx_table[val][0],approx_table[val][1],approx_table[val][2],approx_table[val][3]);
	for(int i = 0; i<3;i++) approx_table[val][i] = approx_table[val][i+1];
	approx_table[val][3] = write_val;

//	printf("(AFTER) [0] = %lf, [1] = %lf, [2] = %lf, [3] = %lf \n",approx_table[val][0],approx_table[val][1],approx_table[val][2],approx_table[val][3]);
}


double fadd_volt_approx[7][9] =
{{0.000698,	0.055853,	0.577385,	5.370308,	5.673313,	16.134663,	24.162198,	43.227072,	4.798509},
{0.000000,	0.077753,	0.693472,	5.749510,	5.984870,	16.494816,	23.990614,	42.521014,	4.487952},
{0.014666,	0.958181,	2.963936,	11.768444,	8.127776,	17.738218,	21.709920,	34.171858,	2.547001},
{22.048257,	28.328576,	14.072892,	18.456390,	4.987415,	2.395281,	6.067999,	3.608740,	0.034449},
{92.849817,	3.320782,	0.172818,	0.131370,	0.037936,	0.077979,	2.320402,	1.088194,	0.000703},
{99.182326,	0.385853,	0.151834,	0.112134,	0.002786,	0.056415,	0.067559,	0.041093,	0.000000},
{100.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000}};


double fmul_volt_approx[7][9] =
{{0.001396,	0.218527,	1.184791,	7.575123,	6.661221,	17.406725,	23.461238,	39.770442,	3.720537},
{0.003502,	0.614318,	2.088120,	10.306108,	7.878257,	18.122723,	22.313673,	35.837069,	2.836229},
{0.004889,	0.574769,	2.086069,	9.795514,	7.628433,	17.847166,	22.649943,	36.360589,	3.052630},
{0.004218,	0.518146,	2.082425,	10.053572,	7.929667,	17.917856,	22.411029,	36.133804,	2.949282},
{0.054094,	2.309865,	4.123052,	15.440546,	10.010116,	18.378458,	19.130148,	28.890872,	1.662850},
{71.375838,	11.593698,	2.901559,	3.384920,	2.583265,	3.701124,	2.259399,	2.189054,	0.011144},
{100.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000}};

double dadd_volt_approx[7][9] =
{{0.000000,	0.043286,	0.537589,	5.269074,	5.782227,	16.062751,	24.127988,	43.317834,	4.859249},
{0.000000,	0.075651,	0.680163,	5.813253,	5.944943,	16.417064,	23.993416,	42.623984,	4.451527},
{0.000000,	0.915580,	2.959745,	11.779618,	7.981116,	17.710981,	21.776266,	34.279409,	2.597285},
{2.213895,	48.395647,	14.195222,	18.391007,	4.905862,	2.331304,	5.959729,	3.576400,	0.030934},
{84.834839,	11.346297,	0.160173,	0.155958,	0.037233,	0.065334,	2.294409,	1.105054,	0.000703},
{97.849253,	1.713354,	0.153923,	0.131636,	0.002786,	0.041789,	0.070345,	0.036914,	0.000000},
{100.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000}};


double dmul_volt_approx[7][9] =
{{0.000000,	0.182920,	1.178508,	7.862768,	6.791080,	17.308283,	23.408875,	39.486986,	3.780580},
{0.000000,	0.568086,	2.146960,	10.240964,	7.939899,	18.202578,	22.229616,	35.813953,	2.857943},
{0.000000,	0.533564,	2.107719,	9.817163,	7.655669,	17.838785,	22.517250,	36.462553,	3.067296},
{0.000000,	0.577905,	2.048679,	10.099973,	7.707504,	18.087290,	22.325961,	36.190048,	2.962640},
{0.000000,	2.355528,	4.272688,	15.615472,	9.901227,	18.396723,	18.828067,	28.985008,	1.645287},
{39.002493,	44.012314,	2.887629,	3.301341,	2.664754,	3.605009,	2.303278,	2.214128,	0.009054},
{100.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000,	0.000000}};


void printBits(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;

    for (i=size-1;i>=0;i--)
    {
        for (j=7;j>=0;j--)
        {
            byte = (b[i] >> j) & 1;
            printf("%u", byte);
        }
    }
    puts("");
}



int sig_bit_lo(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;
    int count = 63;
    int LSM = -1;

    for (i=size-1;i>=0;i--)
    {
        for (j=7;j>=0;j--)
        {
            if(count < 52) byte = (b[i] >> j) & 1;
            else byte = 0;
            count--;
            //printf("%u", byte);
            if(byte) LSM = count; //if byte is high
        }
    }
    //puts("");
    if (LSM < 0) LSM = 0;
    if (LSM > 64) LSM = 64;
    return LSM;
}



int sig_bit_hi(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;
    int count = 0;
    int MSM = -1;

    for (i=0;i<=size-1;i++)
    {
        for (j=0;j<=7;j++)
        {
            if(count < 52) byte = (b[i] >> j) & 1;
            else byte = 0;
            count++;
            //printf("%u", byte);
            if(byte) MSM = count; //if byte is high
        }
    }
    //puts("");
    if (MSM < 0) MSM = 0;
    if (MSM > 64) MSM = 64;
    return MSM;
}


double gok_MUL(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2)
{

	//Stuff for LVA
	int use1, refresh1, use2, refresh2;

	if(dop1_val >=0 ){
	use1 = tv(v[dop1_val],t[dop1_val]);
	refresh1 = tv(v[dop1_val+1],t[dop1_val+1]);
	}
	if(dop2_val>=0){
	use2 = tv(v[dop2_val],t[dop2_val]);
	refresh2 = tv(v[dop2_val+1],t[dop2_val+1]);
	}

	int rand_use1 = 6+rand()%6;
	int rand_use2 = 6+rand()%6;
	int rand_ref1 = 6+rand()%6;
	int rand_ref2 = 6+rand()%6;

	double dop1_approx, dop2_approx;

	if(approx_table[dop1_val][0] == -123){
		dop1_approx = dop1;
		write_approx_table(dop1,dop1_val);
	}
	else{
		if(dop1_val >=0 && rand_use1 >= use1) dop1_approx = read_approx_table(dop1_val);
		else dop1_approx = dop1;

		if(dop1_val >=0 && rand_ref1 < refresh1) write_approx_table(dop1,dop1_val);
	}



	if(approx_table[dop2_val][0] == -123){
		dop2_approx = dop2;
		write_approx_table(dop2,dop2_val);
	}
	else{
		if(dop2_val >=0 && rand_use2 >= use2) dop2_approx = read_approx_table(dop2_val);
		else dop2_approx = dop2;

		if(dop2_val >=0 && rand_ref2 < refresh2) write_approx_table(dop2,dop2_val);
	}

//	printf("gok_MUL: op1 = %lf, op2 = %lf, volt = %d, dop1_val=%d, dop2_val=%d, dop1_approx=%lf, dop2_approx=%lf\n",dop1,dop2,volt,dop1_val,dop2_val,dop1_approx,dop2_approx);

	//LVA: Stuff above

	double op1 = dop1_approx;
	double op2 = dop2_approx;

	int randi = rand()%100;
	int randj = rand()%100;
	int ap = 0;
	double approx = 1.0;
	for(int i =0; i < 9; i++){
		double sum = 0;
		for(int j =0; j<=i;j++){
			sum+=fmul_volt_approx[volt-6][j];
		}
		if(randi < (int)sum) 
		{
			ap=i;
			break;
		}
	}

	if(ap == 0) approx = 1.0;
	else if (ap == 1){
		approx = ((double)randj/100.0)/0.9999  + ((double)(100.0-randj)/100.0)*0.9999;
	}
	else if (ap == 2){
		approx = ((double)randj/100.0)/0.999  + ((double)(100.0-randj)/100.0)*0.999;
	}
	else if (ap == 3){
		approx = ((double)randj/100.0)/0.99  + ((double)(100.0-randj)/100.0)*0.99;
	}
	else if (ap == 4){
		approx = ((double)randj/100.0)/0.98  + ((double)(100.0-randj)/100.0)*0.98;
	}
	else if (ap == 5){
		approx = ((double)randj/100.0)/0.95  + ((double)(100.0-randj)/100.0)*0.95;
	}
	else if (ap == 6){
		approx = ((double)randj/100.0)/0.9  + ((double)(100.0-randj)/100.0)*0.9;
	}
	else if (ap == 7){
		approx = ((double)randj/100.0)/0.75  + ((double)(100.0-randj)/100.0)*0.75;
	}
	else if (ap == 8){
		approx = ((double)randj/100.0)/0.5  + ((double)(100.0-randj)/100.0)*0.5;
	}
//	printf("volt=%d, randi=%d, ap = %d, approx_amt=%lf, pd=%lf, approx_pd=%lf\n",volt,randi,ap,approx,(op1*op2),approx*(op1*op2));
//	printf("ap = %d\n",ap);
	return approx*(op1*op2);

//	float fop1 = (float)op_1;
//	float fop2 = (float)op_2;
//	double op1 = (double)fop1;
//	double op2 = (double)fop2;
//
//	//double op1=op_1;
//	//double op2=op_2;
//	
//	float fsum = fop1*fop2;
//	double sum = (double)fsum;
//	
//	printf("gok_MUL: op1 = %lf, op2 = %lf, volt = %d, hi=%d, lo=%d, out=%lf\n",op1,op2,volt,sig_bit_hi(sizeof(sum),&sum),sig_bit_lo(sizeof(sum),&sum),gok_MUL(op1, op2, tv, sig_bit_lo(sizeof(sum),&sum), sig_bit_hi(sizeof(sum),&sum)));
//	return gok_MUL(op1, op2, tv, sig_bit_lo(sizeof(sum),&sum), sig_bit_hi(sizeof(sum),&sum));
//	//return gok_MUL(op1, op2, tv, 64,0);
}


double gok_ADD(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2)
{

	//Stuff for LVA
	int use1, refresh1, use2, refresh2;

	if(dop1_val >=0 ){
	use1 = tv(v[dop1_val],t[dop1_val]);
	refresh1 = tv(v[dop1_val+1],t[dop1_val+1]);
	}
	if(dop2_val>=0){
	use2 = tv(v[dop2_val],t[dop2_val]);
	refresh2 = tv(v[dop2_val+1],t[dop2_val+1]);
	}

	int rand_use1 = 6+rand()%6;
	int rand_use2 = 6+rand()%6;
	int rand_ref1 = 6+rand()%6;
	int rand_ref2 = 6+rand()%6;

	double dop1_approx, dop2_approx;


	if(approx_table[dop1_val][0] == -123){
		dop1_approx = dop1;
		write_approx_table(dop1,dop1_val);
	}
	else{
		if(dop1_val >=0 && rand_use1 >= use1) dop1_approx = read_approx_table(dop1_val);
		else dop1_approx = dop1;

		if(dop1_val >=0 && rand_ref1 < refresh1) write_approx_table(dop1,dop1_val);
	}



	if(approx_table[dop2_val][0] == -123){
		dop2_approx = dop2;
		write_approx_table(dop2,dop2_val);
	}
	else{
		if(dop2_val >=0 && rand_use2 >= use2) dop2_approx = read_approx_table(dop2_val);
		else dop2_approx = dop2;

		if(dop2_val >=0 && rand_ref2 < refresh2) write_approx_table(dop2,dop2_val);
	}
	

	
//	printf("gok_ADD: op1 = %lf, op2 = %lf, volt = %d, dop1_val=%d, dop2_val=%d, dop1_approx=%lf, dop2_approx=%lf\n",dop1,dop2,volt,dop1_val,dop2_val,dop1_approx,dop2_approx);

	//LVA: Stuff above

	double op1 = dop1_approx;
	double op2 = dop2_approx;

	int randi = rand()%100;
	int randj = rand()%100;
	int ap = 0;
	double approx = 1.0;
	for(int i =0; i < 9; i++){
		double sum = 0;
		for(int j =0; j<=i;j++){
			sum+=fadd_volt_approx[volt-6][j];
		}
		if(randi < (int)sum) 
		{
			ap=i;
			break;
		}
	}

	if(ap == 0) approx = 1.0;
	else if (ap == 1){
		approx = ((double)randj/100.0)/0.9999  + ((double)(100.0-randj)/100.0)*0.9999;
	}
	else if (ap == 2){
		approx = ((double)randj/100.0)/0.999  + ((double)(100.0-randj)/100.0)*0.999;
	}
	else if (ap == 3){
		approx = ((double)randj/100.0)/0.99  + ((double)(100.0-randj)/100.0)*0.99;
	}
	else if (ap == 4){
		approx = ((double)randj/100.0)/0.98  + ((double)(100.0-randj)/100.0)*0.98;
	}
	else if (ap == 5){
		approx = ((double)randj/100.0)/0.95  + ((double)(100.0-randj)/100.0)*0.95;
	}
	else if (ap == 6){
		approx = ((double)randj/100.0)/0.9  + ((double)(100.0-randj)/100.0)*0.9;
	}
	else if (ap == 7){
		approx = ((double)randj/100.0)/0.75  + ((double)(100.0-randj)/100.0)*0.75;
	}
	else if (ap == 8){
		approx = ((double)randj/100.0)/0.5  + ((double)(100.0-randj)/100.0)*0.5;
	}
//	printf("volt=%d, randi=%d, ap = %d, approx_amt=%lf, pd=%lf, approx_pd=%lf\n",volt,randi,ap,approx,(op1*op2),approx*(op1*op2));
//	printf("ap = %d\n",ap);
	return approx*(op1+op2);

//	float fop1 = (float)op_1;
//	float fop2 = (float)op_2;
//	double op1 = (double)fop1;
//	double op2 = (double)fop2;
//
//	//double op1=op_1;
//	//double op2=op_2;
//	
//	float fsum = fop1*fop2;
//	double sum = (double)fsum;
//	
//	printf("gok_MUL: op1 = %lf, op2 = %lf, volt = %d, hi=%d, lo=%d, out=%lf\n",op1,op2,volt,sig_bit_hi(sizeof(sum),&sum),sig_bit_lo(sizeof(sum),&sum),gok_MUL(op1, op2, tv, sig_bit_lo(sizeof(sum),&sum), sig_bit_hi(sizeof(sum),&sum)));
//	return gok_MUL(op1, op2, tv, sig_bit_lo(sizeof(sum),&sum), sig_bit_hi(sizeof(sum),&sum));
//	//return gok_MUL(op1, op2, tv, 64,0);
}




double gok_SUB(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2)
{

	//Stuff for LVA
	int use1, refresh1, use2, refresh2;

	if(dop1_val >=0 ){
	use1 = tv(v[dop1_val],t[dop1_val]);
	refresh1 = tv(v[dop1_val+1],t[dop1_val+1]);
	}
	if(dop2_val>=0){
	use2 = tv(v[dop2_val],t[dop2_val]);
	refresh2 = tv(v[dop2_val+1],t[dop2_val+1]);
	}

	int rand_use1 = 6+rand()%6;
	int rand_use2 = 6+rand()%6;
	int rand_ref1 = 6+rand()%6;
	int rand_ref2 = 6+rand()%6;

	double dop1_approx, dop2_approx;

	if(approx_table[dop1_val][0] == -123){
		dop1_approx = dop1;
		write_approx_table(dop1,dop1_val);
	}
	else{
		if(dop1_val >=0 && rand_use1 >= use1) dop1_approx = read_approx_table(dop1_val);
		else dop1_approx = dop1;

		if(dop1_val >=0 && rand_ref1 < refresh1) write_approx_table(dop1,dop1_val);
	}



	if(approx_table[dop2_val][0] == -123){
		dop2_approx = dop2;
		write_approx_table(dop2,dop2_val);
	}
	else{
		if(dop2_val >=0 && rand_use2 >= use2) dop2_approx = read_approx_table(dop2_val);
		else dop2_approx = dop2;

		if(dop2_val >=0 && rand_ref2 < refresh2) write_approx_table(dop2,dop2_val);
	}
	


//	printf("gok_SUB: op1 = %lf, op2 = %lf, volt = %d, dop1_val=%d, dop2_val=%d, dop1_approx=%lf, dop2_approx=%lf\n",dop1,dop2,volt,dop1_val,dop2_val,dop1_approx,dop2_approx);
	//LVA: Stuff above

	double op1 = dop1_approx;
	double op2 = dop2_approx;

	int randi = rand()%100;
	int randj = rand()%100;
	int ap = 0;
	double approx = 1.0;
	for(int i =0; i < 9; i++){
		double sum = 0;
		for(int j =0; j<=i;j++){
			sum+=fadd_volt_approx[volt-6][j];
		}
		if(randi < (int)sum) 
		{
			ap=i;
			break;
		}
	}

	if(ap == 0) approx = 1.0;
	else if (ap == 1){
		approx = ((double)randj/100.0)/0.9999  + ((double)(100.0-randj)/100.0)*0.9999;
	}
	else if (ap == 2){
		approx = ((double)randj/100.0)/0.999  + ((double)(100.0-randj)/100.0)*0.999;
	}
	else if (ap == 3){
		approx = ((double)randj/100.0)/0.99  + ((double)(100.0-randj)/100.0)*0.99;
	}
	else if (ap == 4){
		approx = ((double)randj/100.0)/0.98  + ((double)(100.0-randj)/100.0)*0.98;
	}
	else if (ap == 5){
		approx = ((double)randj/100.0)/0.95  + ((double)(100.0-randj)/100.0)*0.95;
	}
	else if (ap == 6){
		approx = ((double)randj/100.0)/0.9  + ((double)(100.0-randj)/100.0)*0.9;
	}
	else if (ap == 7){
		approx = ((double)randj/100.0)/0.75  + ((double)(100.0-randj)/100.0)*0.75;
	}
	else if (ap == 8){
		approx = ((double)randj/100.0)/0.5  + ((double)(100.0-randj)/100.0)*0.5;
	}
//	printf("volt=%d, randi=%d, ap = %d, approx_amt=%lf, pd=%lf, approx_pd=%lf\n",volt,randi,ap,approx,(op1*op2),approx*(op1*op2));
//	printf("ap = %d\n",ap);
	return approx*(op1-op2);

//	float fop1 = (float)op_1;
//	float fop2 = (float)op_2;
//	double op1 = (double)fop1;
//	double op2 = (double)fop2;
//
//	//double op1=op_1;
//	//double op2=op_2;
//	
//	float fsum = fop1*fop2;
//	double sum = (double)fsum;
//	
//	printf("gok_MUL: op1 = %lf, op2 = %lf, volt = %d, hi=%d, lo=%d, out=%lf\n",op1,op2,volt,sig_bit_hi(sizeof(sum),&sum),sig_bit_lo(sizeof(sum),&sum),gok_MUL(op1, op2, tv, sig_bit_lo(sizeof(sum),&sum), sig_bit_hi(sizeof(sum),&sum)));
//	return gok_MUL(op1, op2, tv, sig_bit_lo(sizeof(sum),&sum), sig_bit_hi(sizeof(sum),&sum));
//	//return gok_MUL(op1, op2, tv, 64,0);
}








#define DIVIDE 120.0


//Precision to use for calculations
#define fptype float

#define NUM_RUNS 1

typedef struct OptionData_ {
        fptype s;          // spot price
        fptype strike;     // strike price
        fptype r;          // risk-free interest rate
        fptype divq;       // dividend rate
        fptype v;          // volatility
        fptype t;          // time to maturity or option expiration in years 
                           //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)  
        char OptionType;   // Option type.  "P"=PUT, "C"=CALL
        fptype divs;       // dividend vals (not used in this test)
        fptype DGrefval;   // DerivaGem Reference Value
} OptionData;

OptionData *data;
fptype *prices;
int numOptions;
fptype *precisedata;

int    * otype;
fptype * sptprice;
fptype * strike;
fptype * rate;
fptype * volatility;
fptype * otime;
int numError = 0;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244
#define inv_sqrt_2xPI 0.39894228040143270286

fptype CNDF ( fptype InputX ) 
{
    int sign;

    fptype OutputX;
    fptype xInput;
    fptype xNPrimeofX;
    fptype expValues;
    fptype xK2;
    fptype xK2_2, xK2_3;
    fptype xK2_4, xK2_5;
    fptype xLocal, xLocal_1;
    fptype xLocal_2, xLocal_3;

    // Check for negative value of InputX
    if (InputX < 0.0) {
        InputX = -InputX;
        sign = 1;
    } else 
        sign = 0;

    xInput = InputX;
 
    // Compute NPrimeX term common to both four & six decimal accuracy calcs
    expValues = exp(-0.5f * InputX * InputX);
    xNPrimeofX = expValues;
    xNPrimeofX = gok_MUL(xNPrimeofX, inv_sqrt_2xPI, tv(v[18],t[18]),-2,110);

    xK2 =  gok_MUL(0.2316419 , xInput, tv(v[19],t[19]),-2,114);
    xK2 = gok_ADD((fptype)1.0 , xK2, tv(v[0],t[0]),-2,36);
    xK2 = 1.0 / xK2;
    xK2_2 = gok_MUL(xK2 , xK2, tv(v[20],t[20]),116,118);
    xK2_3 = gok_MUL(xK2_2 , xK2, tv(v[21],t[21]),120,122);
    xK2_4 = gok_MUL(xK2_3 , xK2, tv(v[22],t[22]),124,126);
    xK2_5 = gok_MUL(xK2_4 , xK2, tv(v[23],t[23]),128,130);
    
    xLocal_1 = gok_MUL(xK2 , 0.319381530, tv(v[24],t[24]),-2,134);
    xLocal_2 = gok_MUL(xK2_2 , (-0.356563782), tv(v[25],t[25]),-2,138);
    xLocal_3 = gok_MUL(xK2_3 , 1.781477937, tv(v[26],t[26]),-2,142);
    xLocal_2 = gok_ADD(xLocal_2 , xLocal_3, tv(v[1],t[1]),40,42);
    xLocal_3 = gok_MUL(xK2_4 , (-1.821255978), tv(v[27],t[27]),-2,146);
    xLocal_2 = gok_ADD(xLocal_2 , xLocal_3, tv(v[2],t[2]),44,46);
    xLocal_3 = gok_MUL(xK2_5 , 1.330274429, tv(v[28],t[28]),-2,150);
    xLocal_2 = gok_ADD(xLocal_2 , xLocal_3, tv(v[3],t[3]),48,50);

    xLocal_1 = gok_ADD(xLocal_2 , xLocal_1, tv(v[4],t[4]),52,54);
    xLocal   = gok_MUL(xLocal_1 , xNPrimeofX, tv(v[29],t[29]),152,154);

    //printf("# xLocal: %10.10f\n", xLocal);



    xLocal   = gok_SUB((fptype)1.0 , xLocal, tv(v[5],t[5]),-2,58);

    OutputX  = xLocal;

    //printf("# Output: %10.10f\n", OutputX);
    
    if (sign) {
        OutputX = gok_SUB((fptype)1.0 , OutputX, tv(v[6],t[6]),-2,62);
    }
    
    return OutputX;
} 

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
fptype BlkSchlsEqEuroNoDiv( fptype sptprice,
                            fptype strike, fptype rate, fptype volatility,
                            fptype time, int otype, float timet, fptype* N1, fptype* N2)
{
    fptype OptionPrice;

    // local private working variables for the calculation
    //fptype xStockPrice;
    //fptype xStrikePrice;
    fptype xRiskFreeRate;
    fptype xVolatility;
    fptype xTime;
    fptype xSqrtTime;

    fptype logValues;
    fptype xLogTerm;
    fptype xD1; 
    fptype xD2;
    fptype xPowerTerm;
    fptype xDen;
    fptype d1;
    fptype d2;
    fptype FutureValueX;
    fptype NofXd1;
    fptype NofXd2;
    fptype NegNofXd1;
    fptype NegNofXd2;  
    
    //xStockPrice = sptprice;
    //xStrikePrice = strike;
    xRiskFreeRate = rate;
    xVolatility = volatility;
    xTime = time;


    xSqrtTime = sqrt(xTime);

    logValues = log( (sptprice / strike) );
        
    xLogTerm = logValues;
        
    
    xPowerTerm = gok_MUL(xVolatility,xVolatility, tv(v[30],t[30]),156,158);
    xPowerTerm = gok_MUL(xPowerTerm , 0.5, tv(v[31],t[31]),160,-2);
        
    xD1 = gok_ADD(xRiskFreeRate, xPowerTerm, tv(v[7],t[7]),64,66);
    //if(gok_ADD(xRiskFreeRate, xPowerTerm) != (xRiskFreeRate + xPowerTerm)) printf("Hola!\n");
    xD1 = gok_MUL(xD1,xTime, tv(v[32],t[32]),164,166);
    xD1 = gok_ADD(xD1 , xLogTerm, tv(v[8],t[8]),68,70);

    

    xDen = gok_MUL(xVolatility , xSqrtTime, tv(v[33],t[33]),168,170);
    xD1 = (xD1 / xDen);
    xD2 = gok_SUB(xD1,xDen, tv(v[9],t[9]),72,74);

    d1 = xD1;
    d2 = xD2;
    
    NofXd1 = CNDF( d1 );

    //if(NofXd1 > 1.0)
        //std::cerr << "Greater than one!" << std::endl ; //GOKUL
    //printf("# d1: %10.10f\n", NofXd1);

    NofXd2 = CNDF( d2 );
    //if(NofXd2 > 1.0)
         //std::cerr << "Greater than one!" << std::endl ; //GOKUL
    //printf("# d2: %10.10f\n", NofXd2);

    *N1 = NofXd1 ;
    *N2 = NofXd2 ;

    FutureValueX = (strike * ( exp( -(rate)*(time) ) ));  

    if (otype == 0) {            
        OptionPrice = gok_SUB( gok_MUL(sptprice , NofXd1, tv(v[34],t[34]),172,174) , gok_MUL(FutureValueX , NofXd2, tv(v[35],t[35]),176,178), tv(v[10],t[10]),76,78);
        
    } else { 
        NegNofXd1 = gok_SUB((fptype)1.0, NofXd1, tv(v[15],t[15]),-2,98);//1.0
        NegNofXd2 = gok_SUB((fptype)1.0 ,NofXd2, tv(v[16],t[16]),-2,102);//1.0
        OptionPrice =  gok_SUB((FutureValueX * NegNofXd2) , (sptprice * NegNofXd1), tv(v[11],t[11]));
    }
    
    return OptionPrice;
}


double normalize(double in, double min, double max, double min_new, double max_new)
{
    return gok_ADD(((gok_SUB(in , min, tv(v[17],t[17])) / gok_SUB(max , min, tv(v[13],t[13]))) * gok_SUB(max_new , min_new, tv(v[12],t[12]))) , min_new, tv(v[14],t[14])) ;
}

int bs_thread(void *tid_ptr) {
    int i, j;

    int tid = *(int *)tid_ptr;
    int start = tid * (numOptions);
    int end = start + (numOptions);
    fptype price_orig;

    for (j=0; j<NUM_RUNS; j++) {
        for (i=start; i<end; i++) {
            /* Calling main function to calculate option value based on 
             * Black & Scholes's equation.
             */
            fptype price;
            fptype N1, N2;

            double dataIn[6];
            double dataOut[1];

            dataIn[0]   = sptprice[i];
            dataIn[1]   = strike[i];
            dataIn[2]   = rate[i];
            dataIn[3]   = volatility[i];
            dataIn[4]   = otime[i];
            dataIn[5]   = otype[i];

//#pragma parrot(input, "blackscholes", [6]dataIn)
//m5_start_approx();
//m5_end_approx(); //FIXME: Switch between these 2 - needed for serialization

                price_orig = BlkSchlsEqEuroNoDiv( sptprice[i], strike[i],
                                         rate[i], volatility[i], otime[i], 
                                         otype[i], 0, &N1, &N2);
                dataOut[0] = price_orig;

//#pragma parrot(output, "blackscholes", [1]<0.1; 0.9>dataOut)
//m5_end_approx();
                price_orig = dataOut[0];
                prices[i] = price_orig;
        }
    }
    return 0;
}


char* appendCharToCharArray(char* array, char a)
{
    size_t len = strlen(array);

    char* ret = new char[len+2];

    strcpy(ret, array);    
    ret[len] = a;
    ret[len+1] = '\0';

    return ret;
}


 
char* itoa(int num, char* str, int base)
{
    int i = 0;
    bool isNegative = false;

    /* Handle 0 explicitely, otherwise empty string is printed for 0 */
    if (num == 0)
    {
        str[i++] = '0';
        str[i] = '\0';
        return str;
    }

    // In standard itoa(), negative numbers are handled only with
    // base 10. Otherwise numbers are considered unsigned.
    if (num < 0 && base == 10)
    {
        isNegative = true;
        num = -num;
    }

    // Process individual digits
    while (num != 0)
    {
        int rem = num % base;
        str[i++] = (rem > 9)? (rem-10) + 'a' : rem + '0';
        num = num/base;
    }

    // If number is negative, append '-'
    if (isNegative)
        str[i++] = '-';

    str[i] = '\0'; // Append string terminator

    // Reverse the string
    //std::reverse(str);

    return str;
}





int main (int argc, char **argv)
{
//
//
//
//
//
//
//
//




    FILE *file;
    FILE *file_2, *file_ap, *file_gem5;
    int i;
    int loopnum;
    fptype * buffer;
    int * buffer2;
    int rv, rv_2;


	fflush(NULL);


    char *inputFile = argv[1];
    char *outputFile_base = argv[2];
    char *inputFile_2 = argv[3];
    char* outputFile;
    char approxFile[33];
    char buffy[33];
//	printf("Gokul 1 \n");
    //Read input data from file
    file = fopen(inputFile, "r");
    file_2 =  fopen(inputFile_2, "r");
    
    if(file == NULL) {
      printf("ERROR(): Unable to open file `%s'.\n", inputFile);
      exit(1);
    }
    rv = fscanf(file, "%i", &numOptions);
    if(rv != 1) {
      printf("ERROR(): Unable to read from file `%s'.\n", inputFile);
      fclose(file);
      exit(1);
    }

    if(file_2 == NULL) {
      printf("ERROR(): Unable to open file `%s'.\n", inputFile_2);
      exit(1);
    }

//	printf("Gokul 2 \n");

    // alloc spaces for the option data
    data = (OptionData*)malloc(numOptions*sizeof(OptionData));
    prices = (fptype*)malloc(numOptions*sizeof(fptype));

    precisedata = (fptype*)malloc(numOptions*sizeof(fptype));
    for ( loopnum = 0; loopnum < numOptions; ++ loopnum )
    {
        //rv = fscanf(file, "%lf %lf %lf %lf %lf %lf %c %lf %lf", &data[loopnum].s, &data[loopnum].strike, &data[loopnum].r, &data[loopnum].divq, &data[loopnum].v, &data[loopnum].t, &data[loopnum].OptionType, &data[loopnum].divs, &data[loopnum].DGrefval);
        rv = fscanf(file, "%f %f %f %f %f %f %c %f %f", &data[loopnum].s, &data[loopnum].strike, &data[loopnum].r, &data[loopnum].divq, &data[loopnum].v, &data[loopnum].t, &data[loopnum].OptionType, &data[loopnum].divs, &data[loopnum].DGrefval);
        if(rv != 9) {
          printf("ERROR(): Unable to read from file `%s'.\n", inputFile);
          fclose(file);
          exit(1);
        }
    }
    rv = fclose(file);
    if(rv != 0) {
      printf("ERROR(): Unable to close file `%s'.\n", inputFile);
      exit(1);
    }

//	printf("Gokul 2.5 \n");

    for ( int k = 0; k < numOptions; ++k )
    {


//	printf("Gokul 2.6 \n");
        //rv_2 = fscanf(file_2, "%lf", &precisedata[k]);
        rv_2 = fscanf(file_2, "%f", &precisedata[k]);

//	printf("Gokul 2.7 \n");
        if(rv_2 != 1) {
          printf("ERROR(): Unable to read from file `%s'.\n", inputFile_2);
          fclose(file_2);
          exit(1);
        }
	
//	printf("Gokul 2.8 \n");
    }
    rv_2 = fclose(file_2);
    if(rv_2 != 0) {
      printf("ERROR(): Unable to close file `%s'.\n", inputFile_2);
      exit(1);
    }



//	printf("Gokul 3 \n");

#define PAD 256
#define LINESIZE 64



//Type 6: Gradient Descent w/ gem5 calls
int itn = 0, update_itn=0;
float ipc[variables];
float ipc_old;
float ipc_older;
float loads[variables];
float loads_old;
float vf[variables];
float grad[variables];
int step[variables];
float cost = 0, cost_old=0;
int man = 0; //index being manipulated
float error_threshold = 5.0;
float error[variables]; // This variable has to be set at the end of each iteration - its % error
float error_old;
int g=0;
float delta_count;

for(int i=0; i < variables; i++){
        v[i]= 9+rand()%4; //9 to 12
        vf[i] = 1.0*(float)v[i];
	t[i]=0;
}


for(int i=0; i< variables; i++)
{
	step[i] = 1;
	ipc[i] = 0;
	ipc_old = 0;
        ipc_older = 0;
	loads[i]=0;
	loads_old=0;
	error[i] = 0;
	error_old = 0;
	delta_count = 0;
}
int var_update = 0;
float delta = 0;
while((update_itn < 50) || (update_itn < 100 && delta_count < 5)){ //GOKUL-1 There should also be a convergence condition here //FIXME
g++;
itn++;
int man = 0;
printf("Gok-1\n");
//Step 2: if time for update, officially update the vars
if(var_update){
	if (ipc_older > 0){
 		delta = (fabs(ipc_old - ipc_older) / ipc_older);
		if(update_itn >= 10 && delta < 0.01) delta_count++;
		else delta_count=0;
	}
	//First scale the gradients - this is equivalent of learning rate
	float grad_max = 0;
	for(int k=0; k < variables; k++){
		if(fabs(grad[k]) > grad_max) grad_max = fabs(grad[k]);	
	}
	if(grad_max != 0){
		//Always keep grad_max as 1 (this may or may not be great). note that grad-max is positive thanks to above
		float scale = 1/grad_max; //0.5 or 1 or 2?
		for(int k=0; k < variables; k++){
			grad[k] *= scale;	
		}
	}

	//Next, do the updates
	for(int k=0; k < variables; k++){
		vf[k] += (step[k]*grad[k]*1.0);
		if(vf[k] > 12) vf[k] = 12;
		if(vf[k] < 6) vf[k] = 6;
		v[k] = (int) vf[k];
		//Now set t[k] - t[k] indicates the % of time we will use v[k]+1
		t[k] = (int)(100.0*(vf[k] - v[k]*1.0));
		printf("%f - %d - %d,",vf[k], v[k], t[k]);
	}
	update_itn++;
	printf("Update_itn #%d, delta is %f\n",update_itn, delta);
	var_update = 0;
}

printf("Gok-2\n");
//Put below stuff into a loop inorder to be able to parallelize
//TODO
for(man=-1; man < variables; man++)
{

//Below should be in both loops
//Step 3: Perform step change //TODO this only works if the vars are independent??
if(man > -1){
	//Note that we are not perturbing the t[k]
	step[man] = -1;
	v[man] += step[man];
	if(v[man] < 6) v[man] = 6;
	if(v[man] > 12) v[man] = 12;
}

printf("Gok-3\n");


//Step 4: Check Perf for new setting - this requires first writing to approx.txt, then calling gem5, then reading from gem5
//4.1: Write to approx.txt

    printf("We are here\n");

    //approxFile = strcat("approx-",itoa(man),".txt");
    if(man > -1) snprintf(approxFile, sizeof(approxFile), "approx-%d.txt", man);
    else snprintf(approxFile, sizeof(approxFile), "approx-%d.txt", man+200);
    file_ap = fopen(approxFile, "w");
    if(file_ap == NULL) {
      printf("ERROR(): Unable to open file approx.txt.\n");
      exit(1);
    }
    printf("\n");
    for(i=0; i<variables; i++) {
      printf("%d,",v[i]);
      rv = fprintf(file_ap, "%d %d %d\n", pc[i], v[i], t[i]); //pc[i] is obtained from separate pass
      if(rv < 0) {
        printf("ERROR(): Unable to write to file approx.txt.\n");
        fclose(file_ap);
        exit(1);
      }
    }
	fclose(file_ap);
	printf("We are here-2\n");



//4.2 call gem5
   char job[12];
   snprintf(job, sizeof(job), "job_%d.sh", man);
   char outdir[35];
   snprintf(outdir, sizeof(outdir), "m5out_test_%d", man);
   //strcat(outdir, job);
   printf("We are here-3\n");
   char front[1000] = "echo /research/sgokul/gem5-stable/gem5-stable//build/ARM/gem5.opt --outdir=";
   strcat(front,outdir);

   printf("We are here-3.1\n");
   char middle[110];
   if(man > -1) snprintf(middle, sizeof(middle)," /research/sgokul/gem5-stable/gem5-stable//configs/example/se.py --approx-num=%d",man);
   else snprintf(middle, sizeof(middle)," /research/sgokul/gem5-stable/gem5-stable//configs/example/se.py --approx-num=%d",man+200);
   strcat(front, middle);
   //strcat(front,approxFile);
   
   printf("We are here-2\n");
   char back[500] = " --cpu-type=arm_detailed  --restore-with-cpu=arm_detailed --checkpoint-restore=1 --checkpoint-dir=/research/sgokul/bHive/example_4_gem5/o3_ckpt_1k --mem-size=1GB --caches --l2cache  --cpu-clock=2GHz --sys-clock=2GHz --cmd=/research/sgokul/bHive/example_4_gem5/example.out --options=\"'blackscholesTrain_1K.data out'\" > "; 
   strcat(back, job);
   strcat(front,back);

   for(i=0; i<variables; i++){
		printf("%d (%d), ",v[i],t[i]);
   }
   printf("Phase1  \n");

   system( front); //Writes to job.sh



   printf("We are here-4\n");


//Below should be in both loops
//Step 6: Reset previous man variable to old value
if(man > -1){
	v[man] -= step[man];
	//Ideally below 2 will never trigger
	if(v[man] < 6) v[man] = 6;
	if(v[man] > 12) v[man] = 12;
}
}
//Launch condor_script

system ("rm -rf m5out*");
system("condor_submit condor_script_1.scr");
usleep(60000000*2);
system("condor_submit condor_script_2.scr");
usleep(60000000*2);
system("condor_submit condor_script_3.scr");
usleep(60000000*2);
system("condor_submit condor_script_4.scr");
printf("Launched Runs!\n");

for(man = -1; man < variables; man++){
//4.3 read from gem5/stats.txt


//Below should be in both loops
//Step 3: Perform step change //TODO this only works if the vars are independent??

if(man > -1){
	//Note we ain't no peturbing t[k]
	step[man] = -1;
	v[man] += step[man];
	if(v[man] < 6) v[man] = 6;
	if(v[man] > 12) v[man] = 12;
}

printf("Gok-4\n");

if(man > -1){
	ipc[man] = 0;
	loads[man]=0;
	char outdir[35];
	float ticks=0;
	snprintf(outdir, sizeof(outdir), "m5out_test_%d", man);
	while(!ticks){
		if(man==0) usleep(60000000); //One minute
		else usleep(1000000); //One second
		ticks=0;
		char front[220] = "grep sim_ticks ";
		char back[110] = "/stats.txt | awk '{print $2}'";
		strcat(front, outdir);
		strcat(front, back);
		file_gem5 = popen(front,"r"); 
		fscanf(file_gem5,"%f",&ticks);
		ipc[man] = 1000000000/ticks;
		fclose(file_gem5);

		char front2[220] = "grep commit.loads ";
		char back2[110] = "/stats.txt | awk '{print $2}'";
		strcat(front2, outdir);
		strcat(front2, back2);
		file_gem5 = popen(front2,"r"); 
		fscanf(file_gem5,"%f",&loads[man]);
		fclose(file_gem5);

	}

	printf("Final ipc is %f\n",ipc[man]);
}
else{
        ipc_older = ipc_old;
	ipc_old = 0;
	loads_old=0;
	char outdir[35];
	float ticks=0;
	snprintf(outdir, sizeof(outdir), "m5out_test_%d", man);
	while(!ticks){
		usleep(60000000); //One minute
		ticks=0;
		char front[220] = "grep sim_ticks ";
		char back[110] = "/stats.txt | awk '{print $2}'";
		strcat(front, outdir);
		strcat(front, back);
		file_gem5 = popen(front,"r"); 
		fscanf(file_gem5,"%f",&ticks);
		ipc_old = 1000000000/ticks;
		fclose(file_gem5);
		
		char front2[220] = "grep commit.loads ";
		char back2[110] = "/stats.txt | awk '{print $2}'";
		strcat(front2, outdir);
		strcat(front2, back2);
		file_gem5 = popen(front2,"r"); 
		fscanf(file_gem5,"%f",&loads_old);
		fclose(file_gem5);
	}

	printf("Final ipc_old is %f\n",ipc_old);


}


//Step 5: Check accuracy by running the below portion of the code 

    int g1 = 0;
    printf("{%d} V(T): ",man);
    for(int n=0; n < variables; n++){
	printf("[%d] %d (%d,%d),",n,pc[n],v[n],t[n]);
    }
    printf("\n");

  //  std::cout<<"Here 1\n";

    buffer = (fptype *) malloc(5 * numOptions * sizeof(fptype) + PAD);
    sptprice = (fptype *) (((unsigned long long)buffer + PAD) & ~(LINESIZE - 1));
    strike = sptprice + numOptions;
    rate = strike + numOptions;
    volatility = rate + numOptions;
    otime = volatility + numOptions;

    buffer2 = (int *) malloc(numOptions * sizeof(fptype) + PAD);
    otype = (int *) (((unsigned long long)buffer2 + PAD) & ~(LINESIZE - 1));

    //std::cout<<"Here 2\n";
    for (i=0; i<numOptions; i++) {
        otype[i]      = (data[i].OptionType == 'P') ? 1 : 0;
        sptprice[i]   = data[i].s / DIVIDE;
        strike[i]     = data[i].strike / DIVIDE;
        rate[i]       = data[i].r;
        volatility[i] = data[i].v;    
        otime[i]      = data[i].t;
    }
//std::cout<<"Here 3\n";
    //serial version
    int tid=0;
    std::cout<<"Bada bing!!!!\n";
    //m5_checkpoint(0,0);
    std::cout<<"I am the one who knocks!!!!\n";
    bs_thread(&tid);
    std::cout<<"Done Done-ah Done!!!!\n";
    //m5_exit(0);


    //Write prices to output file
    outputFile = strcat(itoa(g,buffy,10),outputFile_base);
    //g++;
    file = fopen(outputFile, "w");
    if(file == NULL) {
      printf("ERROR(): Unable to open file `%s'.\n", outputFile);
      exit(1);
    }
    //rv = fprintf(file, "%i\n", numOptions);
    if(rv < 0) {
      printf("ERROR(): Unable to write to file `%s'.\n", outputFile);
      fclose(file);
      exit(1);
    }

//	printf("Gokul 5 \n");
    for(i=0; i<numOptions; i++) {
      rv = fprintf(file, "%.18f\n", prices[i]);
      if(rv < 0) {
        printf("ERROR(): Unable to write to file `%s'.\n", outputFile);
        fclose(file);
        exit(1);
      }
    }
    rv = fclose(file);
    if(rv != 0) {
      printf("ERROR(): Unable to close file `%s'.\n", outputFile);
      exit(1);
    }

//	printf("Gokul 6 \n");
    std::cout<<"Calculating error:\n";
    float error_inter= 0;
    for(int l=0; l <numOptions; l++){
	if(precisedata[l]) 
		if (fabs((prices[l]  - precisedata[l]) / precisedata[l]) < 1) error_inter += fabs((prices[l]  - precisedata[l]) / precisedata[l]); //TODO might need to change error metric
		else error_inter += 1;

 	//printf("Inter Error %lf \n",error_inter);
   }
if(man > -1){
   error[man] = 100*error_inter/(1.0*numOptions);
 printf(" Error is %lf \n", error[man]);
}
else {
error_old = 100*error_inter/(1.0*numOptions);
printf(" Error-old is %lf \n", error_old);
}




//Below should be in both loops
//Step 6: Reset previous man variable to old value

if(man > -1){
	v[man] -= step[man];
	//Ideally below 2 will never trigger
	if(v[man] < 6) v[man] = 6;
	if(v[man] > 12) v[man] = 12;


	printf("Gok-5\n");



	//Step 7: Caluclate gradient of man
	//Note: We say ipc but actually use 1000000000/execution-ticks
	//We should never go backwards unless error threshold is breached
//	if((error_old > error_threshold) && (error[man] > error_threshold)) grad[man] = ((error_threshold - error[man])/error_threshold); //Does this make sense
//	else if (error[man] > error_threshold) grad[man] = 0;
//	else if(error[man] != error_old && ipc_old && error_old){ 
//		//ipc - ipc_old is a bit finicky because small variations seem to appear randomly
//		//if(ipc[man] > ipc_old) grad[man] = ((ipc[man] - ipc_old)/ipc_old)/fabs((error[man] - error_old)/(100-error_old));
//		if(ipc[man] > ipc_old) grad[man] = (ipc[man]/error[man]);
//		else  grad[man] = 0;
//	}
//	else if(ipc[man] > ipc_old) grad[man] = 1;
//	else  grad[man] = 0;

	//Clean gradient ipc%*accuracy%
	if(error_old > error_threshold) grad[man] = (error_threshold - error[man])/error_threshold;	
	else if(error[man] > error_threshold) grad[man] = 0;
	else if(ipc[man] > ipc_old) grad[man] = ((ipc[man] - ipc_old)/ipc_old);
	//else if(loads[man] < loads_old) grad[man] = ((loads_old - loads[man])/loads_old);
//	else if (ipc[man] > ipc_old && error[man] <= error_old) grad[man] = 1;
	else grad[man]=0;



	printf("grad: %f, error: %f. error_old: %f, ipc: %f, ipc_old: %f, loads: %f, loads_old: %f\n",grad[man], error[man], error_old, ipc[man], ipc_old, loads[man], loads_old);
}
//For next itn:
//ipc_old[man] = ipc[man];
//error_old[man] = error[man];


printf("Gok-6\n");

//Step 1: Iterate over each variable and perturb it one by one
if(man==variables - 1){ // variables - 1 would have been done in the previous iteration
	//Time for update v update based on gradient, after which var_update goes back to 0
	var_update = 1;
}

}

}//Gokul

//
//
///
//
    free(data);
    free(prices);
    free(outputFile);

    return 0;
}

