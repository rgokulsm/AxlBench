
#include <iostream>
#include <cassert>
#include <cmath>
#include <softinj.hh>

#define ITER_LIMIT 1000000
#define VOLTAGE 7
#define SHIFT 24
#define VAL 16

/* Note: the type of the 'sum' and 'i' variables can be either
   char, short, int, long, float, or double.
   However, when type is float or double the bitwise logic operations
   can not be used.
   Moreover, mixed type wrappers are not support. 
*/


//assumes little endian
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
  float sum1 = 0.0, sum2=0, sum3=0, sum4=0, sum5=0, sum6=0;
  float err1 = 0, err2=0;
  float errx;
  int a,b;
  int i;

  for (i=0; i<ITER_LIMIT; i=i+1) {
//     if(i<200) softinj::gokul_reinitializeFaultModels(12);
//     else if (i < 500) softinj::gokul_reinitializeFaultModels(8);
//     else softinj::gokul_reinitializeFaultModels(6);
  a = rand()%VAL;
 // int ay = a>>SHIFT;
 // int az = ay<<SHIFT;
 // short ax = ay;
 // int aw = ay;
// printf("a: %d, az: %d, ay: %d, ax: %d\n",a,az, ay, ax);
//       printBits(sizeof(a), &a);
//	printBits(sizeof(az), &az);
//        printBits(sizeof(ay), &ay);
//        printBits(sizeof(ax), &ax);
//printf("\n");
  b = rand()%VAL;
 // int by = b>>SHIFT;
 // int bz = by<<SHIFT;
 // short bx = by;
 // int bw = by;
    //if(sum + i !=  softinj::ADD(sum, i))    printf("%f: %f - %.64f - %.64f \n",i,sum,sum+i,softinj::ADD(sum, i));
  
    sum1 = (a+b);

  sum2 = softinj::ADD(a, b, VOLTAGE);
 // sumx = softinj::ADD(ax, bx, VOLTAGE);
 // sumw = softinj::ADD(aw, bw, VOLTAGE);
 // sumy = (int) sumx;
 // sumz = sumy << SHIFT;
 // sumwz = sumw << SHIFT;
  //sum3 = softinj::ADD(a, b, 10);
  //sum4 = softinj::ADD(a, b, 9);
  //sum5 = softinj::ADD(a, b, 8);
  //sum6 = softinj::ADD(a, b, 7);

if(fabs((sum2*1.0 - sum1*1.0)/sum1*1.0)<1) errx = fabs((sum2*1.0 - sum1*1.0)/sum1*1.0);
else errx = 1;

err2 += errx;
//err3 += abs((sumz - sum1)*100.0/sum1);
//err4 += abs((sumwz - sum1)*100.0/sum1);
//printf("%f\n",err2);

//if(err2 == 0) err1 +=1;
//else if(err2 < 1)err3 += 1;
//else if(err2 < 5)err4 += 1;
//else if(err2 < 10)err5 += 1;
//else if(err2 < 25)err6 += 1;
//else err7+=1;

}

//  sum = softinj::SUB(sum, i);
//  sum = softinj::XOR(sum, i);
//  sum = softinj::OR(sum, i);
//  sum = softinj::AND(sum, i);

printf("%f-%f \n",err2,100*err2/ITER_LIMIT);
//printf("%f \n",err3/ITER_LIMIT);
//printf("%f \n",err4/ITER_LIMIT);
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

  
  int o = original();
  int m = modified();
    
  printf("original: %d, modified: %d\n", o, m);

  // The call on the finalize function of the library is not
  // essential for the use of the library. However, it may be
  // useful for recording the statistics collected during the
  // execution of the application
  // finalize() prints in stderr
  softinj::finalize();
}

