const int variables = VARIABLES;
int v[VARIABLES];
int t[VARIABLES];


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





void application_approx()
{

}


void application_accurate()
{

}


int
main(int argc, char **argv)
{

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
    else snprintf(approxFile, sizeof(approxFile), "approx-%d.txt", man+100);
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
   char job[10];
   snprintf(job, sizeof(job), "job_%d.sh", man);
   char outdir[33];
   snprintf(outdir, sizeof(outdir), "m5out_test_%d", man);
   //strcat(outdir, job);
   printf("We are here-3\n");
   char front[1000] = "echo /research/sgokul/gem5-stable/gem5-stable//build/ARM/gem5.opt --outdir=";
   strcat(front,outdir);

   printf("We are here-3.1\n");
   char middle[100];
   if(man > -1) snprintf(middle, sizeof(middle)," /research/sgokul/gem5-stable/gem5-stable//configs/example/se.py --approx-num=%d",man);
   else snprintf(middle, sizeof(middle)," /research/sgokul/gem5-stable/gem5-stable//configs/example/se.py --approx-num=%d",man+100);
   strcat(front, middle);
   //strcat(front,approxFile);

   printf("We are here-2\n");
   char back[500] = " --cpu-type=arm_detailed  --restore-with-cpu=arm_detailed --checkpoint-restore=1 --checkpoint-dir=/research/sgokul/bHive/example_gemm/example_4_gem5/o3_ckpt --mem-size=16GB --caches --l2cache  --cpu-clock=2GHz --sys-clock=2GHz --cmd=/research/sgokul/bHive/example_gemm/example_4_gem5/example.out --options=\"'out'\" > ";
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
system("condor_submit condor_script.scr");
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
        char outdir[33];
        snprintf(outdir, sizeof(outdir), "m5out_test_%d", man);
        while(!ipc[man]){
                if(man==0) usleep(60000000); //One minute
                else usleep(1000000); //One second

                char front[200] = "grep ipc_total ";
                char back[100] = "/stats.txt | awk '{print $2}'";
                strcat(front, outdir);
                strcat(front, back);
                file_gem5 = popen(front,"r");
                fscanf(file_gem5,"%f",&ipc[man]);
                fclose(file_gem5);
                //printf("ipc is %f\n",ipc[man]);
        }

        printf("Final ipc is %f\n",ipc[man]);
}
else{
        ipc_older = ipc_old;
        ipc_old = 0;//0
        char outdir[33];
        snprintf(outdir, sizeof(outdir), "m5out_test_%d", man);
        while(!ipc_old){
                usleep(60000000); //One minute

                char front[200] = "grep ipc_total ";
                char back[100] = "/stats.txt | awk '{print $2}'";
                strcat(front, outdir);
                strcat(front, back);
                file_gem5 = popen(front,"r");
                fscanf(file_gem5,"%f",&ipc_old);
                fclose(file_gem5);
        }

        printf("Final ipc_old is %f\n",ipc_old);


}

//Step 5: Check accuracy by running the below portion of the code 

    int g1 = 0;
    printf("V(T): ");
    for(int n=0; n < variables; n++){
        printf("%d (%d),",v[n],t[n]);
    }
    printf("\n");

  //  std::cout<<"Here 1\n";


  for(i=0; i<N; i++){
    z[i]=0;
    z_o[i]=0;
  }


  softinj::initialize(SOFTINJ_STATUS,
		      SOFTINJ_VOLTAGE,
		      SOFTINJ_START_BIT, SOFTINJ_END_BIT,
		      SOFTINJ_RANDOM_SEED,
		      SOFTINJ_DEBUG_LEVEL);


	outputFile = strcat(itoa(g,buffy,10),outputFile_base);
        file=fopen(outputFile, "w");
	int r;
	bb_gemm(x, y, z);
	bb_gemm_o(x,y, z_o);
	float error_inter = 0;
	for (int i = 0; i < N; i++){ 

		fprintf(file,"%lf %lf\n", z[i], z_o[i]);

		if(z_o[i]){
			if (fabs(z[i] - z_o[i])/fabs(z_o[i]) < 1) error_inter += fabs(z[i] - z_o[i])/fabs(z_o[i]); //TODO might need to change error metric
			else error_inter += 1;
		}
	}

        //printf("Inter Error %lf \n",error_inter);



if(man > -1){
   error[man] = 100*error_inter/(N*1.0);
 printf("Inter: %lf, N: %d",error_inter, N);
 printf(" Error is %lf \n", error[man]);
}
else {
error_old = 100*error_inter/(N*1.0);
 printf("Inter: %lf, N: %d",error_inter, N);
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
        //We should never go backwards unless error threshold is breached
        if(error[man] > 2*error_threshold) grad[man] = ((error_threshold - error[man])/error_threshold); //Does this make sense
        else if (error[man] > error_threshold) grad[man] = 0;
        else if(error[man] != error_old && ipc_old && error_old && error[man]){
                //ipc - ipc_old is a bit finicky because small variations seem to appear randomly
                //if(ipc[man] > ipc_old) grad[man] = ((ipc[man] - ipc_old)/ipc_old)/fabs((error[man] - error_old)/(100-error_old));
                if(ipc[man] >= ipc_old) grad[man] = (ipc[man]/error[man]);
                else  grad[man] = 0;
        }
        else if(ipc[man] > ipc_old) grad[man] = 1;
        else  grad[man] = 0;

        printf("grad: %f, error: %f. error_old: %f, ipc: %f. ipc_old: %f\n",grad[man], error[man], error_old, ipc[man], ipc_old);
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


	
		//fprintf(file,"%d\n", z[i]);


//}

//printf("%f\n",err);
 
  // The call on the finalize function of the library is not
  // essential for the use of the library. However, it may be
  // useful for recording the statistics collected during the
  // execution of the application
  // finalize() prints in stderr
  softinj::finalize();
  free(x);
  free(y);
  free(z);
  free(z_o);
  return 0;
}

