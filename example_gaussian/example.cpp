
#include <iostream>
#include <cassert>
#include <math.h>

#include <softinj.hh>

#define ITER_LIMIT 1000000

/* Note: the type of the 'sum' and 'i' variables can be either
   char, short, int, long, float, or double.
   However, when type is float or double the bitwise logic operations
   can not be used.
   Moreover, mixed type wrappers are not support. 
*/

double fRand()
{
    return  (double)rand() * ((double)rand() / RAND_MAX);
//    return fMin + f * (fMax - fMin);
}



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


int
original()
{
  double sum = 0;
  double i;

  for (i=0; i<ITER_LIMIT; i=i+0.001) {
    sum = sum + i;
  }

//  sum = sum - i;
//  sum = sum ^ i;
//  sum = sum | i;
//  sum = sum & i;


  return sum;
}

int
modified()
{
  int level_array[15];
  float bit_array[52];
  int count_array[52];
  for (int j = 0; j < 15; j++){
	level_array[j] = 0;
  }
  for (int j = 0; j < 52; j++){
	bit_array[j] = 0;
	count_array[j]=0;
  }
  printf("Starting...\n");
  double a, b, sum_approx, sum_accurate;
  int level, start_bit, end_bit;
  start_bit = 51;
  end_bit = 51;
  for (int i=0; i<ITER_LIMIT; i++) {

//  softinj::initialize(SOFTINJ_STATUS,
//		      SOFTINJ_VOLTAGE,
//		      SOFTINJ_START_BIT, SOFTINJ_END_BIT,
//		      SOFTINJ_RANDOM_SEED,
//		      SOFTINJ_DEBUG_LEVEL);
	int k =rand();
  	a = fRand();
	b = fRand();


	if(k%100<100) sum_accurate = a+b; //Only does +
//	else if (k%100<34) sum_accurate = a-b;
//	else if (k%100<51) sum_accurate = a^b;
//	else if (k%100<68) sum_accurate = a|b;
//	else if (k%100<84) sum_accurate = a&b;
//	else sum_accurate = a;
	printf("accurate: %lf, ",sum_accurate);

	printBits(sizeof(sum_accurate), &sum_accurate);

	//end_bit = sig_bit_hi(sizeof(sum_accurate), &sum_accurate);
	//start_bit = sig_bit_lo(sizeof(sum_accurate), &sum_accurate);
	end_bit = 51;
	start_bit = 0;
	//end_bit = rand()%50 + 2;
	//start_bit = rand()%(end_bit-1) + 1;


//	start_bit--;
//	if(start_bit ==0){
//		//end_bit--;
//		//start_bit = end_bit-1;
//		break;
//	}
	//if(end_bit ==0) break;
	printf(" hi: %d: lo: %d,",end_bit, start_bit);

	for (int v = 12; v >=6; v--){
		if(k%100<100) sum_approx = softinj::ADD(a, b, v, start_bit, end_bit); //Only does +
//		else if (k%100<34) sum_approx = softinj::SUB(a, b, v);
//		else if (k%100<51) sum_approx = softinj::XOR(a, b, v);
//		else if (k%100<68) sum_approx = softinj::OR(a, b, v);
//		else if (k%100<84) sum_approx = softinj::AND(a, b, v);
//		else sum_approx = softinj::MOV(a, v);
//		printf("app[%d]: %lf,",v,sum_approx);

//		printBits(sizeof(sum_approx), &sum_approx);
		if(sum_approx >= 0.9999*sum_accurate && sum_approx <= 1.0001*sum_accurate) level = v;
		else break;
	}
	sum_approx = softinj::ADD(a, b, level, start_bit, end_bit);

	printf(",L = %d\n",level);
	level_array[level]++;
	//bit_array[start_bit]+=level;
	
	count_array[end_bit-start_bit]++;
	bit_array[end_bit-start_bit] = (((bit_array[end_bit-start_bit]*(count_array[end_bit-start_bit]-1))+level)/count_array[end_bit-start_bit]);
	
  }
  
 
  for (int j = 0; j < 15; j++){
  	printf("%d %d\n",j,level_array[j]);
  }
  printf("\n");
  for (int j = 51; j <= 51; j++){
  	printf("%d %f\n",j,bit_array[j]);
  }
//  sum = softinj::SUB(sum, i);
//  sum = softinj::XOR(sum, i);
//  sum = softinj::OR(sum, i);
//  sum = softinj::AND(sum, i);


  return 0;
}


//int
//modified_2ops()
//{
//  int level_array[15];
//  for (int j = 0; j < 15; j++){
//	level_array[j] = 0;
//  }
//  int a, b, sum_approx, sum_accurate, level;
//  int a2, b2, sum_approx2, sum_accurate2, level2;
//  for (int i=0; i<ITER_LIMIT; i++) {
//	int k =rand();
//  	a = rand();
//	b = rand();
//	if(k%100<17) sum_accurate = a+b;
//	else if (k%100<34) sum_accurate = a-b;
//	else if (k%100<51) sum_accurate = a^b;
//	else if (k%100<68) sum_accurate = a|b;
//	else if (k%100<84) sum_accurate = a&b;
//	else sum_accurate = a;
//
//	int k2 = rand();
//  	a2 = rand();
//	b2 = rand();
//	if(k2%100<17) sum_accurate2 = a2+b2;
//	else if (k2%100<34) sum_accurate2 = a2-b2;
//	else if (k2%100<51) sum_accurate2 = a2^b2;
//	else if (k2%100<68) sum_accurate2 = a2|b2;
//	else if (k2%100<84) sum_accurate2 = a2&b2;
//	else sum_accurate2 = a2;
//
//	for (int v = 12; v >=6; v--){
//		if(k%100<17) sum_approx = softinj::ADD(a, b, v);
//		else if (k%100<34) sum_approx = softinj::SUB(a, b, v);
//		else if (k%100<51) sum_approx = softinj::XOR(a, b, v);
//		else if (k%100<68) sum_approx = softinj::OR(a, b, v);
//		else if (k%100<84) sum_approx = softinj::AND(a, b, v);
//		else sum_approx = softinj::MOV(a, v);
//
//
//		if(sum_accurate == sum_approx) level = v;
//		else break;
//	}
//
//	for (int v = 12; v >=6; v--){
//		if(k2%100<17) sum_approx2 = softinj::ADD(a2, b2, v);
//		else if (k2%100<34) sum_approx2 = softinj::SUB(a2, b2, v);
//		else if (k2%100<51) sum_approx2 = softinj::XOR(a2, b2, v);
//		else if (k2%100<68) sum_approx2 = softinj::OR(a2, b2, v);
//		else if (k2%100<84) sum_approx2 = softinj::AND(a2, b2, v);
//		else sum_approx2 = softinj::MOV(a2, v);
//
//
//		if(sum_accurate2 == sum_approx2) level2 = v;
//		else break;
//	}
//
//
//	//printf("%d\n",level);
//	level_array[(level + level2)/2 + 1]++;
//	
//  }
//  
//  for (int j = 0; j < 15; j++){
//  	printf("%d\n",level_array[j]);
//  }
//  printf("\n");
////  sum = softinj::SUB(sum, i);
////  sum = softinj::XOR(sum, i);
////  sum = softinj::OR(sum, i);
////  sum = softinj::AND(sum, i);
//
//
//  return 0;
//}
//
//int
//modified_4ops()
//{
//  int level_array[15];
//  for (int j = 0; j < 15; j++){
//	level_array[j] = 0;
//  }
//  int a, b, sum_approx, sum_accurate, level;
//  int a2, b2, sum_approx2, sum_accurate2, level2;
//  int a3, b3, sum_approx3, sum_accurate3, level3;
//  int a4, b4, sum_approx4, sum_accurate4, level4;
//  for (int i=0; i<ITER_LIMIT; i++) {
//
////  softinj::initialize(SOFTINJ_STATUS,
////		      SOFTINJ_VOLTAGE,
////		      SOFTINJ_START_BIT, SOFTINJ_END_BIT,
////		      SOFTINJ_RANDOM_SEED,
////		      SOFTINJ_DEBUG_LEVEL);
//
//	int k=rand();
//  	a = rand();
//	b = rand();
//	if(k%100<17) sum_accurate = a+b;
//	else if (k%100<34) sum_accurate = a-b;
//	else if (k%100<51) sum_accurate = a^b;
//	else if (k%100<68) sum_accurate = a|b;
//	else if (k%100<84) sum_accurate = a&b;
//	else sum_accurate = a;
//
//	int k2=rand();
//  	a2 = rand();
//	b2 = rand();
//	if(k2%100<17) sum_accurate2 = a2+b2;
//	else if (k2%100<34) sum_accurate2 = a2-b2;
//	else if (k2%100<51) sum_accurate2 = a2^b2;
//	else if (k2%100<68) sum_accurate2 = a2|b2;
//	else if (k2%100<84) sum_accurate2 = a2&b2;
//	else sum_accurate2 = a2;
//
//	int k3=rand();
//  	a3 = rand();
//	b3 = rand();
//	if(k3%100<17) sum_accurate3 = a3+b3;
//	else if (k3%100<34) sum_accurate3 = a3-b3;
//	else if (k3%100<51) sum_accurate3 = a3^b3;
//	else if (k3%100<68) sum_accurate3 = a3|b3;
//	else if (k3%100<84) sum_accurate3 = a3&b3;
//	else sum_accurate3 = a3;
//
//	int k4=rand();
//  	a4 = rand();
//	b4 = rand();
//	if(k4%100<17) sum_accurate4 = a4+b4;
//	else if (k4%100<34) sum_accurate4 = a4-b4;
//	else if (k4%100<51) sum_accurate4 = a4^b4;
//	else if (k4%100<68) sum_accurate4 = a4|b4;
//	else if (k4%100<84) sum_accurate4 = a4&b4;
//	else sum_accurate4 = a4;
//
//	
//	for (int v = 12; v >=6; v--){
//  		sum_approx = softinj::ADD(a, b, v);
//		if(k%100<17) sum_approx = softinj::ADD(a, b, v);
//		else if (k%100<34) sum_approx = softinj::SUB(a, b, v);
//		else if (k%100<51) sum_approx = softinj::XOR(a, b, v);
//		else if (k%100<68) sum_approx = softinj::OR(a, b, v);
//		else if (k%100<84) sum_approx = softinj::AND(a, b, v);
//		else sum_approx = softinj::MOV(a, v);
//		
//		if(sum_accurate == sum_approx) level = v;
//		else break;
//	}
//
//	for (int v = 12; v >=6; v--){
//		if(k2%100<17) sum_approx2 = softinj::ADD(a2, b2, v);
//		else if (k2%100<34) sum_approx2 = softinj::SUB(a2, b2, v);
//		else if (k2%100<51) sum_approx2 = softinj::XOR(a2, b2, v);
//		else if (k2%100<68) sum_approx2 = softinj::OR(a2, b2, v);
//		else if (k2%100<84) sum_approx2 = softinj::AND(a2, b2, v);
//		else sum_approx2 = softinj::MOV(a2, v);
//
//
//		if(sum_accurate2 == sum_approx2) level2 = v;
//		else break;
//	}
//
//	for (int v = 12; v >=6; v--){
//		if(k3%100<17) sum_approx3 = softinj::ADD(a3, b3, v);
//		else if (k3%100<34) sum_approx3 = softinj::SUB(a3, b3, v);
//		else if (k3%100<51) sum_approx3 = softinj::XOR(a3, b3, v);
//		else if (k3%100<68) sum_approx3 = softinj::OR(a3, b3, v);
//		else if (k3%100<84) sum_approx3 = softinj::AND(a3, b3, v);
//		else sum_approx3 = softinj::MOV(a3, v);
//
//
//		if(sum_accurate3 == sum_approx3) level3 = v;
//		else break;
//	}
//
//	for (int v = 12; v >=6; v--){
//		if(k4%100<17) sum_approx4 = softinj::ADD(a4, b4, v);
//		else if (k4%100<34) sum_approx4 = softinj::SUB(a4, b4, v);
//		else if (k4%100<51) sum_approx4 = softinj::XOR(a4, b4, v);
//		else if (k4%100<68) sum_approx4 = softinj::OR(a4, b4, v);
//		else if (k4%100<84) sum_approx4 = softinj::AND(a4, b4, v);
//		else sum_approx4 = softinj::MOV(a4, v);
//
//
//		if(sum_accurate4 == sum_approx4) level4 = v;
//		else break;
//	}
//
//
//	//printf("%d\n",level);
//	level_array[(level + level2 + level3 + level4)/4 + 1]++;
//	
//  }
//  
//  for (int j = 0; j < 15; j++){
//  	printf("%d\n",level_array[j]);
//  }
//  printf("\n");
////  sum = softinj::SUB(sum, i);
////  sum = softinj::XOR(sum, i);
////  sum = softinj::OR(sum, i);
////  sum = softinj::AND(sum, i);
//
//
//  return 0;
//}
//

int
main(int argc, char **argv)
{

  // The initialization function of the library needs to be
  // before any use of the wrapping functions
  softinj::initialize(SOFTINJ_STATUS,
		      SOFTINJ_VOLTAGE,
		      SOFTINJ_START_BIT, SOFTINJ_END_BIT,
		      SOFTINJ_RANDOM_SEED,
		      SOFTINJ_DEBUG_LEVEL);

  
 // int o = original();
  int m = modified();
    
  //printf("%d\n",m);

  // The call on the finalize function of the library is not
  // essential for the use of the library. However, it may be
  // useful for recording the statistics collected during the
  // execution of the application
  // finalize() prints in stderr
  softinj::finalize();
}

