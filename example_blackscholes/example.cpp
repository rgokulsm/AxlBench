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


#include <cassert>
#include <softinj.hh>

const int variables = 18;
int v[variables];
int t[variables];

int tv(int v, int t){ //Voltage based on threshold
    int rand_val = rand()%100;
    int V4Ven;
    if(rand_val < t){
	V4Ven = v;
    }
    else {
    	V4Ven = 12;
    }
return V4Ven;
}




#define DIVIDE 120.0


//Precision to use for calculations
#define fptype double

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
    xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

    xK2 = 0.2316419 * xInput;
    xK2 = softinj::ADD((fptype)1.0 , xK2, tv(v[0],t[0]));
    xK2 = 1.0 / xK2;
    xK2_2 = xK2 * xK2;
    xK2_3 = xK2_2 * xK2;
    xK2_4 = xK2_3 * xK2;
    xK2_5 = xK2_4 * xK2;
    
    xLocal_1 = xK2 * 0.319381530;
    xLocal_2 = xK2_2 * (-0.356563782);
    xLocal_3 = xK2_3 * 1.781477937;
    xLocal_2 = softinj::ADD(xLocal_2 , xLocal_3, tv(v[1],t[1]));
    xLocal_3 = xK2_4 * (-1.821255978);
    xLocal_2 = softinj::ADD(xLocal_2 , xLocal_3, tv(v[2],t[2]));
    xLocal_3 = xK2_5 * 1.330274429;
    xLocal_2 = softinj::ADD(xLocal_2 , xLocal_3, tv(v[3],t[3]));

    xLocal_1 = softinj::ADD(xLocal_2 , xLocal_1, tv(v[4],t[4]));
    xLocal   = xLocal_1 * xNPrimeofX;

    //printf("# xLocal: %10.10f\n", xLocal);



    xLocal   = softinj::SUB((fptype)1.0 , xLocal, tv(v[5],t[5]));

    OutputX  = xLocal;

    //printf("# Output: %10.10f\n", OutputX);
    
    if (sign) {
        OutputX = softinj::SUB((fptype)1.0 , OutputX, tv(v[6],t[6]));
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
        
    
    xPowerTerm = (xVolatility*xVolatility);
    xPowerTerm = (xPowerTerm * 0.5);
        
    xD1 = softinj::ADD(xRiskFreeRate, xPowerTerm, tv(v[7],t[7]));
    //if(softinj::ADD(xRiskFreeRate, xPowerTerm) != (xRiskFreeRate + xPowerTerm)) printf("Hola!\n");
    xD1 = (xD1*xTime);
    xD1 = softinj::ADD(xD1 , xLogTerm, tv(v[8],t[8]));

    

    xDen = (xVolatility * xSqrtTime);
    xD1 = (xD1 / xDen);
    xD2 = softinj::SUB(xD1,xDen, tv(v[9],t[9]));

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
        OptionPrice = softinj::SUB((sptprice * NofXd1) , (FutureValueX * NofXd2), tv(v[10],t[10]));
        
    } else { 
        NegNofXd1 = softinj::SUB((fptype)1.0, NofXd1, tv(v[15],t[15]));//1.0
        NegNofXd2 = softinj::SUB((fptype)1.0 ,NofXd2, tv(v[16],t[16]));//1.0
        OptionPrice =  softinj::SUB((FutureValueX * NegNofXd2) , (sptprice * NegNofXd1), tv(v[11],t[11]));
    }
    
    return OptionPrice;
}


double normalize(double in, double min, double max, double min_new, double max_new)
{
    return softinj::ADD(((softinj::SUB(in , min, tv(v[17],t[17])) / softinj::SUB(max , min, tv(v[13],t[13]))) * softinj::SUB(max_new , min_new, tv(v[12],t[12]))) , min_new, tv(v[14],t[14])) ;
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

 softinj::initialize(SOFTINJ_STATUS,
                      SOFTINJ_VOLTAGE,
                      SOFTINJ_START_BIT, SOFTINJ_END_BIT,
                      SOFTINJ_RANDOM_SEED,
                      SOFTINJ_DEBUG_LEVEL);



    FILE *file;
    FILE *file_2;
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
    char buffy[33];
//	printf("Gokul 1 \n");
    //Read input data from file
    file = fopen(inputFile, "r");
    file_2 =  fopen(inputFile_2, "r");
    
    if(file == NULL) {
      printf("ERROR()1: Unable to open file `%s'.\n", inputFile);
      exit(1);
    }
    rv = fscanf(file, "%i", &numOptions);
    if(rv != 1) {
      printf("ERROR(): Unable to read from file `%s'.\n", inputFile);
      fclose(file);
      exit(1);
    }

    if(file_2 == NULL) {
      printf("ERROR()2: Unable to open file `%s'.\n", inputFile_2);
      exit(1);
    }

//	printf("Gokul 2 \n");

    // alloc spaces for the option data
    data = (OptionData*)malloc(numOptions*sizeof(OptionData));
    prices = (fptype*)malloc(numOptions*sizeof(fptype));

    precisedata = (fptype*)malloc(numOptions*sizeof(fptype));
    for ( loopnum = 0; loopnum < numOptions; ++ loopnum )
    {
        rv = fscanf(file, "%lf %lf %lf %lf %lf %lf %c %lf %lf", &data[loopnum].s, &data[loopnum].strike, &data[loopnum].r, &data[loopnum].divq, &data[loopnum].v, &data[loopnum].t, &data[loopnum].OptionType, &data[loopnum].divs, &data[loopnum].DGrefval);
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
        rv_2 = fscanf(file_2, "%lf", &precisedata[k]);

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

int man = 0; //index being manipulated
float error_threshold = 10; //TODO
float error=0; // This variable has to be set at the end of each iteration - its % error
int g=0;
int upper = 0;
for(int i=0; i < variables; i++){
        v[i]=12;
	t[i]=100;
}


//Type 0: Just for accurate. No dynamic %.
while(man < variables  ){//GOKUL - Need to improve this condition
g++;
	printf("Count:%d\n",g);
	if(error > error_threshold) // Dial back on current variable and move to next variable
	{
		v[man]++;
		if(v[man] > 12) v[man] = 12;
		//man++;
		upper = 1;
		//Might cause a few duplicates
	}
	else if(v[man] == 6 || upper ==1 ) //6 is min?
	{
		man++;
		upper = 0;
	}
	else
	{
		v[man]--;
	}





////Type 1: Decreasing assignment. No dynamic %.
//while(error < error_threshold || man < variables  ){//GOKUL - Need to improve this condition
//g++;
//	printf("Count:%d\n",g);
//	if(error > error_threshold) // Dial back on current variable and move to next variable
//	{
//		v[man]++;
//		if(v[man] > 12) break;
//		man++;
//		//Might cause a few duplicates
//	}
//	else if(v[man] == 6) //6 is min?
//	{
//		man++;
//	}
//	else
//	{
//		v[man]--;
//	}
//
//



////Type 2: Decereasing as a whole. No dynamic %.
//int loop_break = 0;
//while(!loop_break && (error < error_threshold || man < variables)  ){//GOKUL - Need to improve this condition
//g++;
//        printf("Count:%d\n",g);
//	printf("Man=%d\n",man);
//        if(error < error_threshold) // Try reducing all vars
//        {
//                for(int p=man; p < variables; p++)
//		{
//			v[p]--;
//			if(v[p] < 6){
//				loop_break = 1;
//		 		break;
//			}
//		}
//        }
//        else
//        {	//Dial it back up by one if error becomes too big
//		for(int p=man; p < variables; p++)
//                {
//                        v[p]++;
//                        if(v[p] > 12){ 
//				loop_break = 1;
//				break;
//			}
//                }
//		man++; // Now in next iterartion we try only on subset
//		if(man >= variables) break;
//                
//        }
//
//

////Type 3: Solo Test.
//for (int z=0; z<1; z++){
//g++;
//        printf("Count:%d\n",g);
//
////v[0] = 10;
////v[1] = 10;
////v[2] = 9;
////v[3] = 10;
////v[4] = 11;
////v[5] = 11;
////v[6] = 12;
////v[7] = 8;
////v[8] = 9;
////v[9] = 11;
////v[10] = 8;
////v[11] = 12;
////v[12] = 12;
////v[13] = 12;
////v[14] = 12;
////v[15] = 12;
////v[16] = 12;
////v[17] = 12;
////
////
//
//v[0] = 11;
//v[1] = 11;
//v[2] = 10;
//v[3] = 10;
//v[4] = 10;
//v[5] = 10;
//v[6] = 10;
//v[7] = 10;
//v[8] = 10;
//v[9] = 10;
//v[10] = 10;
//v[11] = 10;
//v[12] = 10;
//v[13] = 10;
//v[14] = 10;
//v[15] = 10;
//v[16] = 10;
//v[17] = 10;
//
//
////v[0] = 12;
////v[1] = 12;
////v[2] = 12;
////v[3] = 12;
////v[4] = 12;
////v[5] = 12;
////v[6] = 12;
////v[7] = 12;
////v[8] = 12;
////v[9] = 12;
////v[10] = 12;
////v[11] = 12;
////v[12] = 12;
////v[13] = 12;
////v[14] = 12;
////v[15] = 12;
////v[16] = 12;
////v[17] = 12;
////
//




////Type 4: Decreasing as a whole. w/ dynamic %. Where, first try reducing dynamic %, then do whole
//int loop_break = 0;
//while(!loop_break && (error < error_threshold || man < variables)  ){//GOKUL - Need to improve this condition
//g++;
//        printf("Count:%d\n",g);
//        if(error < error_threshold) // Try reducing all vars
//	{
//		int flag=0;
//		for(int p=man; p < variables; p++)
//		{
//			if(t[p] < 100) t[p] = t[p]+25;
//			else flag = 1;
//		}
//
//		if (flag){
//			for(int p=man; p < variables; p++)
//			{
//				v[p]--;
//				t[p] = 100;
//				if(v[p] < 6){ loop_break = 1; break;}
//			}
//		}
//	
//	}
//        else //error > threshold
//        {
//		
//		int flag=0;
//		if(1){ // In this algo design we vary % first - this is not necessarily better, and alternate should also be explored
//			for(int p=man; p < variables; p++)
//			{
//				if(t[p] > 25) t[p] = t[p]-25;
//				else flag = 1;
//			}
//		}
//
//		if(flag){
//
//			for(int p=man; p < variables; p++)
//			{
//				v[p]++;
//				t[p] = 25; //Might loste some opportunity but prevents issues
//
//				if(v[p] > 12){ loop_break = 1; break;}
//			}
//		}
//
//		man++; // Now in next iterartion we try only on subset
//		if(man >= variables){
//			printf("She's the man\n");
//			break;
//		}
//
//
//        }
//
//

////Type 5: Decreasing assignment. w/ dynamic %., dynamic first
//while(error < error_threshold || man < variables  ){//GOKUL - Need to improve this condition
//g++;
//	printf("Count:%d\n",g);
//	if(error > error_threshold) // Dial back on current variable and move to next variable
//	{
//		if(t[man] > 25) {
//			t[man] -= 25;
//		}
//		else {
//			v[man]++;
//			//t[man] = 100; //potentially losing out a bit without doing this, but otherwise might cause issues
//			if(v[man] > 12) break;
//		}
//
//		man++;
//		//Might cause a few duplicates
//	}
//	else if(v[man] == 6) //6 is min?
//	{
//		man++;
//	}
//	else // error < threshold
//	{
//		if(t[man] < 100) t[man] += 25;
//		else {
//			v[man]--;
//			t[man] = 25; // Start from min
//		}
//	}
//
//





    int g1 = 0;
printf("V(T): ");
for(int n=0; n < 18; n++){
printf("%d (%d),",v[n],t[n]);
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
    file = fopen(outputFile, "w");
    if(file == NULL) {
      printf("ERROR()3: Unable to open file `%s'.\n", outputFile);
      exit(1);
    }
    //rv = fprintf(file, "%i\n", numOptions);
    if(rv < 0) {
      printf("ERROR(): Unable to write to file `%s'.\n", outputFile);
      fclose(file);
      exit(1);
    }

//	printf("Gokul 5 \n");
    //rv = fprintf(file, "%d-%d-%d-%d-%d-%d-%d-%d-%d-%d-%d-%d-%d-%d-%d-%d-%d-%d\n", v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18);
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
		if (abs((prices[l]  - precisedata[l]) / precisedata[l]) < 1) error_inter += abs((prices[l]  - precisedata[l]) / precisedata[l]); //TODO might need to change error metric
		else error_inter += 1;

 	//printf("Inter Error %lf \n",error_inter);
   }
   error = 100*error_inter/numOptions;
 printf(" Error is %lf \n", error);


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

