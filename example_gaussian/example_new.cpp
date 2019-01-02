
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
  int level_array[15][15];
 // float bit_array[52];
  //int count_array[52];

  for (int j = 0; j < 15; j++){
	for(int i=0; i<15;i++)	level_array[i][j] = 0;
  }

 // for (int j = 0; j < 52; j++){
	//bit_array[j] = 0;
	//count_array[j]=0;
  //}
  printf("Starting...\n");
  double a, b, sum_approx, sum_accurate;
  int level, start_bit, end_bit;
  for (int i=0; i<ITER_LIMIT; i++) {

	int k =rand();
  	a = fRand();
	b = fRand();
	int v = 12 - rand()%7;


	if(k%100<100) sum_accurate = a*b; //Only does +
//	printf("accurate: %lf, ",sum_accurate);

	//printBits(sizeof(sum_accurate), &sum_accurate);

	//end_bit = sig_bit_hi(sizeof(sum_accurate), &sum_accurate);
	//start_bit = sig_bit_lo(sizeof(sum_accurate), &sum_accurate);
	end_bit = 51;
	start_bit = 29;//fp
	//end_bit = rand()%50 + 2;
	//start_bit = rand()%(end_bit-1) + 1;


	//printf(" hi: %d: lo: %d,",end_bit, start_bit);

	//for (int v = 12; v >=6; v--){
		if(k%100<100) sum_approx = softinj::MUL(a, b, v, start_bit, end_bit); //Only does +
//		printf("app[%d]: %lf,",v,sum_approx);

//		printBits(sizeof(sum_approx), &sum_approx);

		if(sum_approx == sum_accurate) level = 0;
		else if(sum_approx >= 0.9999*sum_accurate && sum_approx <= sum_accurate/0.9999) level = 1;
		else if(sum_approx >= 0.999*sum_accurate && sum_approx <= sum_accurate/0.999) level = 2;
		else if(sum_approx >= 0.99*sum_accurate && sum_approx <= sum_accurate/0.99) level = 3;
		else if(sum_approx >= 0.98*sum_accurate && sum_approx <= sum_accurate/0.98) level = 4;
		else if(sum_approx >= 0.95*sum_accurate && sum_approx <= sum_accurate/0.95) level = 5;
		else if(sum_approx >= 0.9*sum_accurate && sum_approx <= sum_accurate/0.9) level = 6;
		else if(sum_approx >= 0.75*sum_accurate && sum_approx <= sum_accurate/0.75) level = 7;
		else if(sum_approx >= 0.5*sum_accurate && sum_approx <= sum_accurate/0.5) level = 8;
//		else if(sum_approx >= 0.9*sum_accurate && sum_approx <= sum_accurate/0.9999) level = 9;
		else level = 9;
	//}
	//sum_approx = softinj::ADD(a, b, level, start_bit, end_bit);

	//printf(",L = %d\n",level);
	//level_array[level]++;
	//bit_array[start_bit]+=level;
	
	//count_array[v]++;
	//level_array[v] = (((level_array[v]*(count_array[v]-1))+level)/count_array[v]);
	level_array[v][level]++;
	
  }
  
  printf("\taccurate,\t0.01,\t0.1,\t1,\t2,\t5,\t10,\t25,\t50,\tmore,\n" ); 
  for (int j = 6; j < 13; j++){
	printf("[v=%d]:\t",j);
	int total = 0;
	for(int i = 0; i < 10; i++) total+=level_array[j][i];
	for(int i = 0; i < 10; i++) printf("%f,\t",(100.0*(float)level_array[j][i]/(float)total));
	printf("\n");
  }
  printf("\n");
//  for (int j = 51; j <= 51; j++){
//  	printf("%d %f\n",j,bit_array[j]);
//  }
//  sum = softinj::SUB(sum, i);
//  sum = softinj::XOR(sum, i);
//  sum = softinj::OR(sum, i);
//  sum = softinj::AND(sum, i);


  return 0;
}



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

