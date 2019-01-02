
#include <softinj.hh>

#include "pp_scan.h"
//#include "../error_func.c"

//int gok_itn = 0, gok_acc = 0, gok_temp = 0, gok_error = 0;
//int LENGTH = 24;


int voltage = 11;
int threshold = 100;

void print(int *a, int size)
{
	int i;

	for (i = 0; i < size; i++)
		printf("%u\t", a[i]);
}

//typedef enum { false, true } bool;
void local_scan(TYPE bucket[BUCKETSIZE])
{
  int radixID, i;
  loop1_outter:for (radixID = 0; radixID < SCAN_RADIX; ++radixID)
    loop1_inner:for (i = 1; i < SCAN_BLOCK; ++i)
    {
      //gok_temp = calc_length(bucket[radixID * SCAN_BLOCK + i],bucket[radixID * SCAN_BLOCK + i - 1]);
      //gok_acc += gok_temp;
      //if(gok_temp > LENGTH) {
	//printf("Err Length = %d\n",gok_temp);
	//gok_error++;
      //}
      //gok_itn +=1;

      bucket[radixID * SCAN_BLOCK + i] += bucket[radixID * SCAN_BLOCK + i - 1];	
    }
}

void local_scan_m(TYPE bucket[BUCKETSIZE])
{
  int radixID, i;
  loop1_outter:for (radixID = 0; radixID < SCAN_RADIX; ++radixID)
    loop1_inner:for (i = 1; i < SCAN_BLOCK; ++i)
    {
      //gok_temp = calc_length(bucket[radixID * SCAN_BLOCK + i],bucket[radixID * SCAN_BLOCK + i - 1]);
      //gok_acc += gok_temp;
      //if(gok_temp > LENGTH) {
	//printf("Err Length = %d\n",gok_temp);
	//gok_error++;
      //}
      //gok_itn +=1;
				int v;
				if(rand() % 100 < threshold) v = voltage;
				else v = 12;

      bucket[radixID * SCAN_BLOCK + i] = softinj::ADD(bucket[radixID * SCAN_BLOCK + i], bucket[radixID * SCAN_BLOCK + i - 1],v);	
    }
}



void sum_scan(TYPE sum[SCAN_RADIX], TYPE bucket[BUCKETSIZE])
{
  int radixID;
  sum[0] = 0;
  loop2:for (radixID = 1; radixID < SCAN_RADIX; ++radixID){
      //gok_temp = calc_length(sum[radixID -1], bucket[radixID * SCAN_BLOCK - 1]);
      //gok_acc += gok_temp;
      //if(gok_temp > LENGTH) {
	//printf("Err Length = %d\n",gok_temp);
	//gok_error++;
      //}
      //gok_itn +=1;

     sum[radixID] = sum[radixID -1] + bucket[radixID * SCAN_BLOCK - 1];
  }
}

void sum_scan_m(TYPE sum[SCAN_RADIX], TYPE bucket[BUCKETSIZE])
{
  int radixID;
  sum[0] = 0;
  loop2:for (radixID = 1; radixID < SCAN_RADIX; ++radixID){
      //gok_temp = calc_length(sum[radixID -1], bucket[radixID * SCAN_BLOCK - 1]);
      //gok_acc += gok_temp;
      //if(gok_temp > LENGTH) {
	//printf("Err Length = %d\n",gok_temp);
	//gok_error++;
      //}
      //gok_itn +=1;

				int v;
				if(rand() % 100 < threshold) v = voltage;
				else v = 12;


     sum[radixID] = softinj::ADD(sum[radixID -1] , bucket[radixID * SCAN_BLOCK - 1], v);
  }
}



void last_step_scan(TYPE bucket[BUCKETSIZE], TYPE bucket2[BUCKETSIZE], TYPE sum[SCAN_RADIX])
{
  int radixID, i;
  loop3_outter:for (radixID = 0; radixID < SCAN_RADIX; ++radixID)
    loop3_inner:for (i = 0; i < SCAN_BLOCK; ++i)
    {

      //gok_temp = calc_length(bucket[radixID * SCAN_BLOCK + i ],sum[radixID]);
      //gok_acc += gok_temp;
      //if(gok_temp > LENGTH) {
	//printf("Err Length = %d\n",gok_temp);
	//gok_error++;
      //}
      //gok_itn +=1;

      bucket2[radixID * SCAN_BLOCK + i] =
        bucket[radixID * SCAN_BLOCK + i ] + sum[radixID];
    }

}

void last_step_scan_m(TYPE bucket[BUCKETSIZE], TYPE bucket2[BUCKETSIZE], TYPE sum[SCAN_RADIX])
{
  int radixID, i;
  loop3_outter:for (radixID = 0; radixID < SCAN_RADIX; ++radixID)
    loop3_inner:for (i = 0; i < SCAN_BLOCK; ++i)
    {

      //gok_temp = calc_length(bucket[radixID * SCAN_BLOCK + i ],sum[radixID]);
      //gok_acc += gok_temp;
      //if(gok_temp > LENGTH) {
	//printf("Err Length = %d\n",gok_temp);
	//gok_error++;
      //}
      //gok_itn +=1;
				int v;
				if(rand() % 100 < threshold) v = voltage;
				else v = 12;


      bucket2[radixID * SCAN_BLOCK + i] = softinj::ADD(bucket[radixID * SCAN_BLOCK + i ], sum[radixID],v);
    }

}




void pp_scan(TYPE bucket[BUCKETSIZE], TYPE bucket2[BUCKETSIZE],
  TYPE sum[SCAN_RADIX])
{
  local_scan(bucket);
  sum_scan(sum, bucket);
  last_step_scan(bucket, bucket2, sum);
}
void pp_scan_m(TYPE bucket[BUCKETSIZE], TYPE bucket2[BUCKETSIZE],
  TYPE sum[SCAN_RADIX])
{
  local_scan_m(bucket);
  sum_scan_m(sum, bucket);
  last_step_scan_m(bucket, bucket2, sum);
}


int main()
{
	int i;
  TYPE  *bucket, *bucket_o;
  TYPE  *bucket2, *bucket2_o;
  bucket = (TYPE*) malloc(sizeof(TYPE) * N);
  bucket2 = (TYPE*) malloc(sizeof(TYPE) * N);
  bucket_o = (TYPE*) malloc(sizeof(TYPE) * N);
  bucket2_o = (TYPE*) malloc(sizeof(TYPE) * N);
  
  TYPE sum[SCAN_RADIX];
  TYPE sum_o[SCAN_RADIX];
  TYPE max, min;
	srand(8650341L);
  max = 2248;
  min = 0;

  softinj::initialize(SOFTINJ_STATUS,
		      SOFTINJ_VOLTAGE,
		      SOFTINJ_START_BIT, SOFTINJ_END_BIT,
		      SOFTINJ_RANDOM_SEED,
		      SOFTINJ_DEBUG_LEVEL);


  for(i=0; i<N; i++){
    bucket[i] = (TYPE)(((double) rand() / (RAND_MAX)) * (max-min) + min);
    bucket2[i] = (TYPE)(((double) rand() / (RAND_MAX)) * (max-min) + min);
    bucket_o[i] = bucket[i];
    bucket2_o[i] = bucket2[i];
  }

	//print(&bucket[0], 1);
	//printf("\n");

	pp_scan_m(bucket, bucket2, sum);
	
	pp_scan(bucket_o, bucket2_o, sum_o);
	
        //printf("Avg Length = %f, Errors = %d, Itns = %d\n",(float)gok_acc/gok_itn,gok_error,gok_itn); 
	//print(&bucket[0], 2);
	//printf("\n");
	for (i = 0; i < N; i++){ 

		printf("%lf %lf\n", bucket2[i], bucket2_o[i]);
	}

  softinj::finalize();
	return 0;
}
