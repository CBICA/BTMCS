///////////////////////////////////////////////////////////////////////////////////////
// OriginalLandmarks.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"

//landmarks (manually-placed, typically voxels)

//InputLandmarkFileUndef/Def: txt files containing manually-placed landmarks (original and target)


void  OriginalLandmarks(char InputLandmarkFileUndef[800],char InputLandmarkFileDef[800], double *undeflmark, double *deflmarkobj)
{
	int ii,jj;
	double tmp1,tmp2,res_x,res_y,res_z;  
	FILE  *fp1,*fp2;

	res_x=gres_x;
	res_y=gres_y;
	res_z=gres_z;

	fp1=fopen(InputLandmarkFileUndef,"r"); //undeformed landmarks coordinates, written in the input file as x,y,z
	fp2=fopen(InputLandmarkFileDef,"r"); //objective (deformed) landmarks coordinates, written in the input file as x,y,z

	for(jj=0;jj<nblmark;jj++) {

		//printf("lmark %d\n", jj);

		for(ii=0; ii<3;ii++) {

			fscanf(fp1,"%lf",&tmp1);
			fscanf(fp2,"%lf",&tmp2);

			//printf("undeflmark %g\n",tmp1);
			//printf("deflmarkobj %g\n",tmp2);


			if (ii==0){
				undeflmark[jj*3+ii]=tmp1*res_x;
				deflmarkobj[jj*3+ii]=tmp2*res_x;
			}else if (ii==1){
				undeflmark[jj*3+ii]=tmp1*res_y;
				deflmarkobj[jj*3+ii]=tmp2*res_y;
			}else if (ii==2){
				undeflmark[jj*3+ii]=tmp1*res_z;
				deflmarkobj[jj*3+ii]=tmp2*res_z;
			}
		}//end for ii=0,2
	}//end for jj

	fclose(fp1);
	fclose(fp2);
} //end function OriginalLandmarks


