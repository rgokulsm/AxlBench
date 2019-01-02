/*
 * distance.c
 * 
 * Created on: Sep 9, 2013
 * 			Author: Amir Yazdanbakhsh <a.yazdanbakhsh@gatech.edu>
 */

//#include</research/sgokul/gem5-stable/gem5-stable//util/m5/m5op.h>
#include "distance.h"
#include <math.h>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <softinj.hh>

extern int v[49];
extern int t[49];

extern int tv(int v, int t);


extern double fadd_volt_approx[7][9], fmul_volt_approx[7][9], dadd_volt_approx[7][9], dmul_volt_approx[7][9];
extern double gok_ADD(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);
extern double gok_SUB(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);
extern double gok_MUL(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);


int count = 0;
#define MAX_COUNT 1200000


float euclideanDistance_p(RgbPixel* p, Centroid* c1) {
	float r;

	r = 0;
	double r_tmp;
	
	double dataIn[6];
	dataIn[0] = p->r;
	dataIn[1] = p->g;
	dataIn[2] = p->b;
	dataIn[3] = c1->r;
	dataIn[4] = c1->g;
	dataIn[5] = c1->b;

//#pragma parrot(input, "kmeans", [6]dataIn)
//m5_start_approx();
//m5_end_approx(); //FIXME: Switch between these 2 - needed for serialization

	r += ((p->r - c1->r) * (p->r - c1->r));
	r += ((p->g - c1->g) * (p->g - c1->g));
	r += ((p->b - c1->b) * (p->b - c1->b));

	r_tmp = sqrt(r);

//m5_end_approx();
//#pragma parrot(output, "kmeans", <0.0; 1.0>r_tmp)

	return r_tmp;
}








float euclideanDistance(RgbPixel* p, Centroid* c1) {
	float r;

	r = 0;
	double r_tmp;
	
	double dataIn[6];
	dataIn[0] = p->r;
	dataIn[1] = p->g;
	dataIn[2] = p->b;
	dataIn[3] = c1->r;
	dataIn[4] = c1->g;
	dataIn[5] = c1->b;

//#pragma parrot(input, "kmeans", [6]dataIn)
//m5_start_approx();
//m5_end_approx(); //FIXME: Switch between these 2 - needed for serialization
	float t1, t2;

	t1 = gok_SUB(p->r , c1->r,tv(v[4],t[4]),13,15);
	t2 = gok_MUL(t1 , t1,tv(v[5],t[5]),17,19);
	r = gok_ADD(r,t2,tv(v[6],t[6]),21,23);
	t1 = gok_SUB(p->g , c1->g,tv(v[7],t[7]),25,27);
	t2 = gok_MUL(t1,t1,tv(v[8],t[8]),29,31);
	r = gok_ADD(r,t2,tv(v[9],t[9]),33,35);
	t1 = gok_SUB(p->b , c1->b,tv(v[10],t[10]),37,39);
	t2 = gok_MUL(t1,t1,tv(v[11],t[11]),41,43);
	r = gok_ADD(r,t2,tv(v[12],t[12]),45,47);

	//r += ((p->r - c1->r) * (p->r - c1->r));
	//r += ((p->g - c1->g) * (p->g - c1->g));
	//r += ((p->b - c1->b) * (p->b - c1->b));

	r_tmp = sqrt(r);

//m5_end_approx();
//#pragma parrot(output, "kmeans", <0.0; 1.0>r_tmp)

	return r_tmp;
}

int pickCluster(RgbPixel* p, Centroid* c1) {
	float d1;

	d1 = euclideanDistance(p, c1);

	if (p->distance <= d1)
		return 0;

	return 1;
}



void assignCluster_p(RgbPixel* p, Clusters* clusters) {
	int i = 0;

	float d;
	d = euclideanDistance_p(p, &clusters->centroids[i]);
	p->distance = d;

	for(i = 1; i < clusters->k; ++i) {
		d = euclideanDistance_p(p, &clusters->centroids[i]);
		if (d < p->distance) {
			p->cluster = i;
			p->distance = d;
		}
	}
}



void assignCluster(RgbPixel* p, Clusters* clusters) {
	int i = 0;

	float d;
	d = euclideanDistance(p, &clusters->centroids[i]);
	p->distance = d;

	for(i = 1; i < clusters->k; ++i) {
		d = euclideanDistance(p, &clusters->centroids[i]);
		if (d < p->distance) {
			p->cluster = i;
			p->distance = d;
		}
	}
}

