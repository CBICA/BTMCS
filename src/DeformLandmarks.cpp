///////////////////////////////////////////////////////////////////////////////////////
// DeformLandmarks.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "function.h"


#undef __FUNCT__
#define __FUNCT__ "DeformLandmarks"
void DeformLandmarks(int nblmarks, double *undeflmark, double *vel, double dt, int xm, int ym, int zm, double *deflmark)
{
	int i,j,k,ii,jj;
	int knl1,knl2,knl3,knl4,knl5,knl6,knl7,knl8;
	int ic,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ic8;
	//FILE *fp;
	double wx,wy,wz; //weights for trilinear interpolation
	double *velxlmark, *velylmark, *velzlmark;

	velxlmark = (double*)malloc(nblmarks * sizeof(double));
	velylmark = (double*)malloc(nblmarks * sizeof(double));
	velzlmark = (double*)malloc(nblmarks * sizeof(double));

	nnc=xm*ym*zm;
	ne=(xm-1)*(ym-1)*(zm-1);
	Hx = Lx/(xm-1); Hy = Ly/(ym-1); Hz = Lz/(zm-1);


	for(jj=0;jj<nblmarks;jj++) {

		//printf("lmark %d\n",jj);

		for(ii=0; ii<3;ii++) {

			if (ii==0){

				i=(int)(undeflmark[jj*3+ii]/Hx);//left most corner x-index on the regular grid
				wx=(undeflmark[jj*3+ii]-i*Hx)/Hx;

				//printf("undeflmark wx %g\n",undeflmark[jj*3+ii],wx);

			}else if (ii==1){

				j=(int)(undeflmark[jj*3+ii]/Hy);
				wy=(undeflmark[jj*3+ii]-j*Hy)/Hy;

				//printf("wy %g\n",wy);

			}else if (ii==2){

				k=(int)(undeflmark[jj*3+ii]/Hz);
				wz=(undeflmark[jj*3+ii]-k*Hz)/Hz;

				//printf("wz %g\n",wz);

			}

		}//end for ii=0,2


		ic=i+j*ym+k*xm*ym;
		ic1=ic;
		ic2=ic1+1;
		ic3=ic+xm;
		ic4=ic3+1;

		ic5=ic+xm*ym;
		ic6=ic5+1;
		ic7=ic5+xm;
		ic8=ic7+1;  

		knl1 = ic1 * nsd;
		knl2 = ic2 * nsd;
		knl3 = ic3 * nsd;
		knl4 = ic4 * nsd;
		knl5 = ic5 * nsd;
		knl6 = ic6 * nsd;
		knl7 = ic7 * nsd;
		knl8 = ic8 * nsd;


		//printf("knl1-8 %d %d %d %d %d %d %d %d\n",knl1,knl2,knl3,knl4,knl5,knl6,knl7,knl8);

		wx=1.0-wx;
		wy=1.0-wy;
		wz=1.0-wz;       		

		velxlmark[jj]=wz*(wy*(wx*vel[knl1]+(1.0-wx)*vel[knl2])+(1.0-wy)*(wx*vel[knl3]+(1.0-wx)*vel[knl4]))+(1.0-wz)*(wy*(wx*vel[knl5]+(1.0-wx)*vel[knl6])+(1.0-wy)*(wx*vel[knl7]+(1.0-wx)*vel[knl8]));

		velylmark[jj]=wz*(wy*(wx*vel[knl1+1]+(1.0-wx)*vel[knl2+1])+(1.0-wy)*(wx*vel[knl3+1]+(1.0-wx)*vel[knl4+1]))+(1.0-wz)*(wy*(wx*vel[knl5+1]+(1.0-wx)*vel[knl6+1])+(1.0-wy)*(wx*vel[knl7+1]+(1.0-wx)*vel[knl8+1]));

		velzlmark[jj]=wz*(wy*(wx*vel[knl1+2]+(1.0-wx)*vel[knl2+2])+(1.0-wy)*(wx*vel[knl3+2]+(1.0-wx)*vel[knl4+2]))+(1.0-wz)*(wy*(wx*vel[knl5+2]+(1.0-wx)*vel[knl6+2])+(1.0-wy)*(wx*vel[knl7+2]+(1.0-wx)*vel[knl8+2]));			   
	} //end for jj=0,number of landmarks


	//fp=fopen("DefLmark.txt","w");

	for(jj=0;jj<nblmarks;jj++) {

		for(ii=0; ii<3;ii++) {

			//printf("Inisde DeformLandmarks: jj undeflmark[jj*3+ii] velx vely velz %d %g %g %g %g\n",jj,undeflmark[jj*3+ii],velxlmark[jj],velylmark[jj],velzlmark[jj]);

			if (ii==0){
				deflmark[jj*3+ii]= undeflmark[jj*3+ii]+dt*velxlmark[jj];  
				//fprintf(fp,"%lf ",deflmark[jj*3+ii]);			
			}else if (ii==1){
				deflmark[jj*3+ii]= undeflmark[jj*3+ii]+dt*velylmark[jj]; 
				//fprintf(fp,"%lf ",deflmark[jj*3+ii]); 
			}else if (ii==2){
				deflmark[jj*3+ii]= undeflmark[jj*3+ii]+dt*velzlmark[jj];  
				//fprintf(fp,"%lf \n",deflmark[jj*3+ii]);
			}
		}
	}

	//fclose(fp);			

	free(velxlmark);
	free(velylmark);
	free(velzlmark);
} //end function DefLmark
