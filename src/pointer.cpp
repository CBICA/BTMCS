///////////////////////////////////////////////////////////////////////////////////////
// pointer.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"

// Version 1.1 of STSolver by F. Abraham and G. Biros 
// August 2005


void rowcolindex(int i, int j, int k, MatStencil *row, MatStencil *col,int jacobian_flag)
{
	if (nsd == 2) {
		row[0].i = i;row[0].j=j;row[0].c=0;
		row[1].i = i;row[1].j=j;row[1].c=1;
		row[2].i = i+1;row[2].j=j;row[2].c=0;
		row[3].i = i+1;row[3].j=j;row[3].c=1;
		row[4].i = i+1;row[4].j=j+1;row[4].c=0;
		row[5].i = i+1;row[5].j=j+1;row[5].c=1;
		row[6].i = i;row[6].j=j+1;row[6].c = 0;
		row[7].i = i;row[7].j=j+1;row[7].c=1;

		if (jacobian_flag ==1){
			col[0].i = i;col[0].j=j;col[0].c=0;
			col[1].i = i;col[1].j=j;col[1].c=1;
			col[2].i = i+1;col[2].j=j;col[2].c=0;
			col[3].i = i+1;col[3].j=j;col[3].c=1;
			col[4].i = i+1;col[4].j=j+1;col[4].c=0;
			col[5].i = i+1;col[5].j=j+1;col[5].c=1;
			col[6].i = i;col[6].j=j+1;col[6].c = 0;
			col[7].i = i;col[7].j=j+1;col[7].c=1;
		}
	} else if (nsd == 3) {
		row[0].i = i;row[0].j=j;row[0].k=k;row[0].c=0; 
		row[1].i = i;row[1].j=j;row[1].k=k;row[1].c=1; 
		row[2].i = i;row[2].j=j;row[2].k=k;row[2].c=2; 

		row[3].i = i+1;row[3].j=j;row[3].k=k;row[3].c=0;
		row[4].i = i+1;row[4].j=j;row[4].k=k;row[4].c=1;
		row[5].i = i+1;row[5].j=j;row[5].k=k;row[5].c=2;

		row[6].i = i+1;row[6].j=j+1;row[6].k=k;row[6].c=0;
		row[7].i = i+1;row[7].j=j+1;row[7].k=k;row[7].c=1;
		row[8].i = i+1;row[8].j=j+1;row[8].k=k;row[8].c=2;

		row[9].i = i;row[9].j=j+1;row[9].k=k;row[9].c = 0; 
		row[10].i = i;row[10].j=j+1;row[10].k=k;row[10].c=1; 
		row[11].i = i;row[11].j=j+1;row[11].k=k;row[11].c=2; 


		row[12].i = i;row[12].j=j;row[12].k=k+1;row[12].c=0; 
		row[13].i = i;row[13].j=j;row[13].k=k+1;row[13].c=1; 
		row[14].i = i;row[14].j=j;row[14].k=k+1;row[14].c=2; 

		row[15].i = i+1;row[15].j=j;row[15].k=k+1;row[15].c=0; 
		row[16].i = i+1;row[16].j=j;row[16].k=k+1;row[16].c=1; 
		row[17].i = i+1;row[17].j=j;row[17].k=k+1;row[17].c=2; 

		row[18].i = i+1;row[18].j=j+1;row[18].k=k+1; row[18].c=0;
		row[19].i = i+1;row[19].j=j+1;row[19].k=k+1; row[19].c=1;
		row[20].i = i+1;row[20].j=j+1;row[20].k=k+1; row[20].c=2;

		row[21].i = i;row[21].j=j+1;row[21].k=k+1; row[21].c=0;
		row[22].i = i;row[22].j=j+1;row[22].k=k+1; row[22].c=1;
		row[23].i = i;row[23].j=j+1;row[23].k=k+1; row[23].c=2;

		if (jacobian_flag == 1) {
			col[0].i = i;col[0].j=j;col[0].k=k;col[0].c=0; 
			col[1].i = i;col[1].j=j;col[1].k=k;col[1].c=1; 
			col[2].i = i;col[2].j=j;col[2].k=k;col[2].c=2; 

			col[3].i = i+1;col[3].j=j;col[3].k=k;col[3].c=0;
			col[4].i = i+1;col[4].j=j;col[4].k=k;col[4].c=1;
			col[5].i = i+1;col[5].j=j;col[5].k=k;col[5].c=2;

			col[6].i = i+1;col[6].j=j+1;col[6].k=k;col[6].c=0;
			col[7].i = i+1;col[7].j=j+1;col[7].k=k;col[7].c=1;
			col[8].i = i+1;col[8].j=j+1;col[8].k=k;col[8].c=2;

			col[9].i = i;col[9].j=j+1;col[9].k=k;col[9].c = 0; 
			col[10].i = i;col[10].j=j+1;col[10].k=k;col[10].c=1; 
			col[11].i = i;col[11].j=j+1;col[11].k=k;col[11].c=2; 


			col[12].i = i;col[12].j=j;col[12].k=k+1;col[12].c=0; 
			col[13].i = i;col[13].j=j;col[13].k=k+1;col[13].c=1; 
			col[14].i = i;col[14].j=j;col[14].k=k+1;col[14].c=2; 

			col[15].i = i+1;col[15].j=j;col[15].k=k+1;col[15].c=0; 
			col[16].i = i+1;col[16].j=j;col[16].k=k+1;col[16].c=1; 
			col[17].i = i+1;col[17].j=j;col[17].k=k+1;col[17].c=2; 

			col[18].i = i+1;col[18].j=j+1;col[18].k=k+1; col[18].c=0;
			col[19].i = i+1;col[19].j=j+1;col[19].k=k+1; col[19].c=1;
			col[20].i = i+1;col[20].j=j+1;col[20].k=k+1; col[20].c=2;

			col[21].i = i;col[21].j=j+1;col[21].k=k+1; col[21].c=0;
			col[22].i = i;col[22].j=j+1;col[22].k=k+1; col[22].c=1;
			col[23].i = i;col[23].j=j+1;col[23].k=k+1; col[23].c=2;
		}
	}
}
