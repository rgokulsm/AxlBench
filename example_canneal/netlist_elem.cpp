// netlist_elem.cpp
//
// Created by Daniel Schwartz-Narbonne on 14/04/07.
//
// Copyright 2007 Princeton University
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.


#include <assert.h>
#include <math.h>

#include "annealer_types.h"
#include "location_t.h"
#include "netlist_elem.h"


extern int v[180];
extern int t[180];

extern int tv(int v, int t);


extern double fadd_volt_approx[7][9], fmul_volt_approx[7][9], dadd_volt_approx[7][9], dmul_volt_approx[7][9];
extern double gok_ADD(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);
extern double gok_SUB(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);
extern double gok_MUL(double dop1, double dop2, int volt, int dop1_val=-2, int dop2_val=-2);



netlist_elem::netlist_elem()
:present_loc(NULL)//start with the present_loc as nothing at all.  Filled in later by the netlist
{
}

//*****************************************************************************************
// Calculates the routing cost using the manhatten distance
// I make sure to get the pointer in one operation, and then use it
// SYNC: Do i need to make this an atomic operation?  i.e. are there misaligned memoery issues that can cause this to fail
//       even if I have atomic writes?
//*****************************************************************************************
routing_cost_t netlist_elem::routing_cost_given_loc(location_t loc)
{
	routing_cost_t fanin_cost = 0;
	routing_cost_t fanout_cost = 0;
	
	for (int i = 0; i< fanin.size(); ++i){
		location_t* fanin_loc = fanin[i]->present_loc.Get();

		//asm volatile("orn r0, r0, #0");
		routing_cost_t temp_x = gok_SUB(loc.x , fanin_loc->x,tv(v[0],t[0]),26,28);
		//asm volatile("orn r0, r1, #0");
		temp_x = fabs(temp_x);


		//asm volatile("orn r0, r0, #1");
		routing_cost_t temp_y = gok_SUB(loc.y , fanin_loc->y,tv(v[1],t[1]),30,32);
		//asm volatile("orn r0, r1, #1");
		temp_y = fabs(temp_y);
		
		//asm volatile("orn r0, r0, #2");
		fanin_cost = gok_ADD(fanin_cost , temp_x,tv(v[2],t[2]),34,36);
		//asm volatile("orn r0, r1, #2");
		//asm volatile("orn r0, r0, #3");
		fanin_cost = gok_ADD(fanin_cost , temp_y,tv(v[3],t[3]),38,40);
		//asm volatile("orn r0, r1, #3");
	}

	for (int i = 0; i< fanout.size(); ++i){
		location_t* fanout_loc = fanout[i]->present_loc.Get();


		//asm volatile("orn r0, r0, #4");
		routing_cost_t temp_x = gok_SUB(loc.x , fanout_loc->x,tv(v[4],t[4]),42,44);
		//asm volatile("orn r0, r1, #4");
		temp_x = fabs(temp_x);


		//asm volatile("orn r0, r0, #5");
		routing_cost_t temp_y = gok_SUB(loc.y, fanout_loc->y,tv(v[5],t[5]),46,48);
		//asm volatile("orn r0, r1, #5");
		temp_y = fabs(temp_y);
		
		//asm volatile("orn r0, r0, #6");
		fanout_cost = gok_ADD(fanout_cost , temp_x,tv(v[6],t[6]),50,52);
		//asm volatile("orn r0, r1, #6");
		//asm volatile("orn r0, r0, #7");
		fanout_cost = gok_ADD(fanout_cost , temp_y,tv(v[7],t[7]),54,56);
		//asm volatile("orn r0, r1, #7");

	}
	
	//asm volatile("orn r0, r0, #8");
	routing_cost_t total_cost = gok_ADD(fanin_cost , fanout_cost,tv(v[8],t[8]),58,60);
	//asm volatile("orn r0, r1, #8");

	return total_cost;
}

//*****************************************************************************************
//  Get the cost change of swapping from our present location to a new location
//*****************************************************************************************
routing_cost_t netlist_elem::swap_cost(location_t* old_loc, location_t* new_loc)
{
	routing_cost_t no_swap = 0;
	routing_cost_t yes_swap = 0;
	
	for (int i = 0; i< fanin.size(); ++i){
		location_t* fanin_loc = fanin[i]->present_loc.Get();

		//asm volatile("orn r0, r0, #9");
		routing_cost_t temp_x = gok_SUB(old_loc->x , fanin_loc->x,tv(v[9],t[9]),62,64);
		//asm volatile("orn r0, r1, #9");
	
		temp_x = fabs(temp_x);

		//asm volatile("orn r0, r0, #10");
		routing_cost_t temp_y = gok_SUB(old_loc->y , fanin_loc->y,tv(v[10],t[10]),66,68);
		//asm volatile("orn r0, r1, #10");
        
	        temp_y = fabs(temp_y);
		

		//asm volatile("orn r0, r0, #11");
		no_swap = gok_ADD(no_swap , temp_x,tv(v[11],t[11]),70,72);
		//asm volatile("orn r0, r1, #11");
	
		//asm volatile("orn r0, r0, #12");
		no_swap = gok_ADD(no_swap , temp_y,tv(v[12],t[12]),74,76);
		//asm volatile("orn r0, r1, #12");
	

		//asm volatile("orn r0, r0, #13");
		temp_x = gok_SUB(new_loc->x , fanin_loc->x,tv(v[13],t[13]),78,80);
		//asm volatile("orn r0, r1, #13");
	
		temp_x = fabs(temp_x);

		//asm volatile("orn r0, r0, #14");
		temp_y = gok_SUB(new_loc->y , fanin_loc->y,tv(v[14],t[14]),82,84);
		//asm volatile("orn r0, r1, #14");
        
	        temp_y = fabs(temp_y);
		

		//asm volatile("orn r0, r0, #15");
		yes_swap = gok_ADD(yes_swap , temp_x,tv(v[15],t[15]),86,88);
		//asm volatile("orn r0, r1, #15");
	
		//asm volatile("orn r0, r0, #16");
		yes_swap = gok_ADD(yes_swap , temp_y,tv(v[16],t[16]),90,92);
		//asm volatile("orn r0, r1, #16");

	
	}
	
	for (int i = 0; i< fanout.size(); ++i){
		location_t* fanout_loc = fanout[i]->present_loc.Get();



		//asm volatile("orn r0, r0, #17");
		routing_cost_t temp_x = gok_SUB(old_loc->x , fanout_loc->x,tv(v[17],t[17]),94,96);
		//asm volatile("orn r0, r1, #17");
	
		temp_x = fabs(temp_x);

		//asm volatile("orn r0, r0, #18");
		routing_cost_t temp_y = gok_SUB(old_loc->y , fanout_loc->y,tv(v[18],t[18]),98,100);
		//asm volatile("orn r0, r1, #18");
        
	        temp_y = fabs(temp_y);
		

		//asm volatile("orn r0, r0, #19");
		no_swap = gok_ADD(no_swap,temp_x,tv(v[19],t[19]),102,104);
		//asm volatile("orn r0, r1, #19");
	
		//asm volatile("orn r0, r0, #20");
		no_swap = gok_ADD(no_swap,temp_y,tv(v[20],t[20]),106,108);
		//asm volatile("orn r0, r1, #20");
	

		//asm volatile("orn r0, r0, #21");
		temp_x = gok_SUB(new_loc->x , fanout_loc->x,tv(v[21],t[21]),110,112);
		//asm volatile("orn r0, r1, #21");
	
		temp_x = fabs(temp_x);

		//asm volatile("orn r0, r0, #22");
		temp_y = gok_SUB(new_loc->y , fanout_loc->y,tv(v[22],t[22]),114,116);
		//asm volatile("orn r0, r1, #22");
        
	        temp_y = fabs(temp_y);
		

		//asm volatile("orn r0, r0, #23");
		yes_swap = gok_ADD(yes_swap,temp_x,tv(v[23],t[23]),118,120);
		//asm volatile("orn r0, r1, #23");
	
		//asm volatile("orn r0, r0, #24");
		yes_swap = gok_ADD(yes_swap,temp_y,tv(v[24],t[24]),122,124);
		//asm volatile("orn r0, r1, #24");


	}
	
	//asm volatile("orn r0, r0, #25");
	routing_cost_t return_val = gok_SUB(yes_swap , no_swap,tv(v[25],t[25]),126,128);	
	//asm volatile("orn r0, r1, #25");

	return return_val;
}

