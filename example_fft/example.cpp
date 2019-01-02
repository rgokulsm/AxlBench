#include <cstdio>
#include <iostream>
#include "fourier.hpp"
#include <fstream>
#include <time.h>
//#include</research/sgokul/gem5-stable/gem5-stable//util/m5/m5op.h>
#include <cmath>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>
//#include</research/sgokul/gem5-stable/gem5-stable//util/m5/m5op.h>
#include <unistd.h>

const int variables = 92;
int v[variables]; //Not using variables because of extern
int t[variables];
int pc[variables] = {
67862,//0
67886,//1
0,//2
0,//3
69894,//4
69938,//5
69966,//6
70010,//7
70054,//8
70098,//9
70126,//10
70170,//11
70214,//12
70258,//13
70286,//14
70330,//15
70374,//16
70418,//17
70446,//18
70490,//19


67862-12,67862-12,67862-4,67862-4,//20
0,0,67886-8,67886-8,//24
69894-8,69894-8,69894-4,69894-4,//28
69938-8,69938-8,69938-4,69938-4,//32
69966-8,69966-8,69966-4,69966-4,//36
70010-12,70010-12,70010-4,70010-4,//40
70054-8,70054-8,70054-4,70054-4,//44
70098-8,70098-8,70098-4,70098-4,//48
70126-8,70126-8,70126-4,70126-4,//52
70170-12,70170-12,70170-4,70170-4,//56
70214-8,70214-8,70214-4,70214-4,//60
70258-8,70258-8,70258-4,70258-4,//64
70286-8,70286-8,70286-4,70286-4,//68
70330-12,70330-12,70330-4,70330-4,//72
70374-8,70374-8,70374-4,70374-4,//76
70418-8,70418-8,70418-4,70418-4,//80
70446-8,70446-8,70446-4,70446-4,//84
70490-12,70490-12,70490-4,70490-4//88

};

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




static int* indices;
static Complex* x;
static Complex* f;
static Complex* precise;


int main(int argc, char* argv[])
{




	int i ;

	int n 						= atoi(argv[1]);
	char* outputFilename_base 	= argv[2];
	char *inputFile_2 = argv[3];
    	char* outputFilename;
    	char buffy[33];
        char approxFile[33];
        FILE *file;
    	FILE *file_2, *file_ap, *file_gem5;
	int rv, rv_2;


    precise = (Complex*)malloc(n * sizeof (Complex));
    file_2 =  fopen(inputFile_2, "r");
    if(file_2 == NULL) {
      printf("ERROR(): Unable to open file `%s'.\n", inputFile_2);
      exit(1);
    }

    for ( int k = 0; k < n; ++k ) //"n" is reassigned as "K" later
    {


//	printf("Gokul 2.6 \n");
        rv_2 = fscanf(file_2, "%lf %lf", &precise[k].real, &precise[k].imag);

//	printf("Gokul 2.7 \n");
        if(rv_2 != 2) {
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



	// create the arrays
	x 		= (Complex*)malloc(n * sizeof (Complex));
	f 		= (Complex*)malloc(n * sizeof (Complex));
	indices = (int*)malloc(n * sizeof (int));

	if(x == NULL || f == NULL || indices == NULL)
	{
		std::cout << "cannot allocate memory for the triangle coordinates!" << std::endl;
		return -1 ;
	}


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
float error_threshold = 25.0;
float error[variables]; // This variable has to be set at the end of each iteration - its % error
float error_old;
int g=0;
float delta_count;


for(int i=0; i < variables; i++){
        v[i]=12;
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
   char back[500] = " --cpu-type=arm_detailed  --restore-with-cpu=arm_detailed --checkpoint-restore=1 --checkpoint-dir=/research/sgokul/bHive/example_fft/example_4_gem5/o3_ckpt_slite --mem-size=1GB --caches --l2cache  --cpu-clock=2GHz --sys-clock=2GHz --cmd=/research/sgokul/bHive/example_fft/example_4_gem5/example.out --options=\"'256 out'\" > "; 
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
system("rm -rf m5out*");
system("condor_submit condor_script_1.scr");
usleep(60000000*2);
system("condor_submit condor_script_2.scr");
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
        ipc[man] = 0;//0
	loads[man]=0;
        char outdir[35];
	float ticks=0;
        snprintf(outdir, sizeof(outdir), "m5out_test_%d", man);
        while(!ticks){
                if(man==0) usleep(60000000); //One minute
                else usleep(1000000); //One second

                char front[220] = "grep sim_ticks ";
                char back[110] = "/stats.txt | awk '{print $2}'";
                strcat(front, outdir);
                strcat(front, back);
                file_gem5 = popen(front,"r");
                fscanf(file_gem5,"%f",&ticks);
		ipc[man] = 1000000000/ticks;	
                fclose(file_gem5);
                //printf("ipc is %f\n",ipc[man]);
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
        ipc_old = 0;//0
	loads_old = 0;
        char outdir[35];
	float ticks=0;
        snprintf(outdir, sizeof(outdir), "m5out_test_%d", man);
        while(!ticks){
                usleep(60000000); //One minute

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


	int K = n;

        for(i = 0;i < K ; i++)
        {
                x[i].real = i;
                x[i].imag = 0 ;
        }



    std::cout<<"g2 - "<<g<<"\n";
    outputFilename = strcat(itoa(g,buffy,10),outputFilename_base);
    //g++;
    std::cout<<"g3 - "<<g<<"\n";
	// prepare the output file for writting the theta values
	std::ofstream outputFileHandler;
	outputFileHandler.open(outputFilename);
	outputFileHandler.precision(5);

	//m5_checkpoint(0,0);
	std::cout<<"I am the one who knocks!!!!\n";	
	radix2DitCooleyTykeyFft(K, indices, x, f) ;
	std::cout<<"Done Done-ah Done!!!!\n";
	//m5_exit(0);

//	for(int xn = 0; xn < 20; xn++){
//	outputFileHandler<<v[xn]<<",";
//	}
//	outputFileHandler<<std::endl;
	
	for(i = 0;i < K ; i++)
	{
		outputFileHandler << f[i].real << " " << f[i].imag << std::endl;
	}


	outputFileHandler.close();



//	printf("Gokul 6 \n");
    std::cout<<"Calculating error:\n";
    float error_inter= 0;
    for(int l=0; l <K; l++){
		if(precise[l].real && precise[l].imag){
			double diffReal        = f[l].real - precise[l].real;
			double diffImag        = f[l].imag - precise[l].imag;
			//if(diffReal > 1.0) diffReal = 1.0;
			//if(diffImag > 1.0) diffImag = 1.0;
			
			//double approx  = pow(f[l].real*f[l].real + f[l].imag*f[l].imag, 0.5);
			//double correct = pow(precise[l].real*precise[l].real+ precise[l].imag*precise[l].imag, 0.5);
			double numerator = sqrt(diffReal*diffReal + diffImag*diffImag);
			//double denominator = fabs(correct);
			double denominator =	sqrt(precise[l].real*precise[l].real + precise[l].imag*precise[l].imag);
		
			if (fabs(numerator/denominator) < 1) error_inter += fabs(numerator/denominator); //TODO might need to change error metric
			else error_inter += 1;
		}

 	//	printf("f.real - %lf, p.real - %lf, f.imag - %lf, p.imag - %lf, error.inter - %lf\n", f[l].real,precise[l].real,f[l].imag,precise[l].imag,error_inter);
 	//printf("Inter Error %lf \n",error_inter);
   }
if(man > -1){
   error[man] = 100*error_inter/(K*1.0);
 printf(" Error is %lf \n", error[man]);
}
else {
error_old = 100*error_inter/(K*1.0);
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
	else if(loads[man] < loads_old) grad[man] = ((loads_old - loads[man])/loads_old);
//	else if (ipc[man] > ipc_old && error[man] <= error_old) grad[man] = 1;
	else grad[man]=0;



	printf("grad: %f, error: %f. error_old: %f, ipc: %f, ipc_old: %f, loads: %f, loads_old: %f\n",grad[man], error[man], error_old, ipc[man], ipc_old, loads[man], loads_old);


        //We should never go backwards unless error threshold is breached
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

	return 0 ;
}
