///////////////////////////////////////////////////////////////////////////////////////
// CumulateDisplacement.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#undef __FUNCT__
#define __FUNCT__ "TrackParticlesLagrangian"
void TrackParticlesLagrangian(double *vel, double timespan, int xm, int ym, int zm, double Hx, double Hy, double Hz, double *deflmark)
{
	int nnc,ne,nsd;
	int i,j,k,ii,jj,knl;
	int knl1,knl2,knl3,knl4,knl5,knl6,knl7,knl8;
	int counti,countj,countk,ic,ic1,ic2,ic3,ic4,ic5,ic6,ic7,ic8;
	double wx,wy,wz; //weights for trilinear interpolation
	double *velxlmark, *velylmark, *velzlmark;

	nnc=xm*ym*zm;
	nsd=3;
	ne=(xm-1)*(ym-1)*(zm-1);

	velxlmark = (double*)malloc(nnc * sizeof(double));
	velylmark = (double*)malloc(nnc * sizeof(double));
	velzlmark = (double*)malloc(nnc * sizeof(double));

	for(jj=0;jj<nnc;jj++) {
		for(ii=0; ii<3;ii++) {
			if (ii==0){
				//left most corner x-index on the regular grid
				i=(int)(deflmark[jj*3+ii]/Hx);
				wx=(deflmark[jj*3+ii]-i*Hx)/Hx;
			} else if (ii==1){
				j=(int)(deflmark[jj*3+ii]/Hy);
				wy=(deflmark[jj*3+ii]-j*Hy)/Hy;
			} else if (ii==2){
				k=(int)(deflmark[jj*3+ii]/Hz);
				wz=(deflmark[jj*3+ii]-k*Hz)/Hz;
			}
		} //end for ii={0,1,2}

		if ((i<xm-1) && (j<ym-1) && (k<zm-1)) {
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

			wx=1.0-wx;
			wy=1.0-wy;
			wz=1.0-wz;                

			velxlmark[jj]=wz*(wy*(wx*vel[knl1]+(1.0-wx)*vel[knl2])+(1.0-wy)*(wx*vel[knl3]+(1.0-wx)*vel[knl4]))+(1.0-wz)*(wy*(wx*vel[knl5]+(1.0-wx)*vel[knl6])+(1.0-wy)*(wx*vel[knl7]+(1.0-wx)*vel[knl8]));
			velylmark[jj]=wz*(wy*(wx*vel[knl1+1]+(1.0-wx)*vel[knl2+1])+(1.0-wy)*(wx*vel[knl3+1]+(1.0-wx)*vel[knl4+1]))+(1.0-wz)*(wy*(wx*vel[knl5+1]+(1.0-wx)*vel[knl6+1])+(1.0-wy)*(wx*vel[knl7+1]+(1.0-wx)*vel[knl8+1]));
			velzlmark[jj]=wz*(wy*(wx*vel[knl1+2]+(1.0-wx)*vel[knl2+2])+(1.0-wy)*(wx*vel[knl3+2]+(1.0-wx)*vel[knl4+2]))+(1.0-wz)*(wy*(wx*vel[knl5+2]+(1.0-wx)*vel[knl6+2])+(1.0-wy)*(wx*vel[knl7+2]+(1.0-wx)*vel[knl8+2]));               
		} else {
			velxlmark[jj]=0;
			velylmark[jj]=0;
			velzlmark[jj]=0;               
		} //end if i,j,k      
	} //end for jj=0,number of landmarks

	counti = 1; countj = 1; countk = 1;
	for (jj=0;jj<nnc;jj++){
		knl = jj * nsd;
		i = counti-1; j = countj-1; k = countk-1;

		deflmark[knl  ] = deflmark[knl  ]+timespan*velxlmark[jj];
		deflmark[knl+1] = deflmark[knl+1]+timespan*velylmark[jj];
		deflmark[knl+2] = deflmark[knl+2]+timespan*velzlmark[jj];

		counti++;
		if (counti==xm+1){ countj++; counti=1; }
		if (countj==ym+1){ counti=1; countj=1; countk++;}
	}

	free(velxlmark);
	free(velylmark);
	free(velzlmark);      
} //end function TrackParticlesLagrangian

#undef __FUNCT__
#define __FUNCT__ "CumulateDisplacement"
void CumulateDisplacement(double *vel, double timespan, int xm, int ym, int zm, double Hx, double Hy, double Hz, double *dispc)
{
	int nnc,ne,nsd;
	double xcoord,ycoord,zcoord;
	double  *deflmark;  
	int jj, counti,countj,countk,i,j,k,knl;

	nnc=xm*ym*zm;
	nsd=3;
	ne=(xm-1)*(ym-1)*(zm-1);

	deflmark=(double*)malloc(nnc*nsd*sizeof(double));

	counti = 1; countj = 1; countk = 1;
	for (jj=0;jj<nnc;jj++) {
		knl = jj * nsd;
		i = counti-1; j = countj-1; k = countk-1;

		//previous (deformed) configuration of the regular Cartesian grid used to solve the FEM problem
		deflmark[knl  ] = i * Hx + dispc[knl  ];
		deflmark[knl+1] = j * Hy + dispc[knl+1];
		deflmark[knl+2] = k * Hz + dispc[knl+2];

		counti++;
		if (counti==xm+1){ countj++; counti=1; }
		if (countj==ym+1){ counti=1; countj=1; countk++;}
	}

	// update the (current) cumulative displacement dispc
	TrackParticlesLagrangian(vel, timespan, xm, ym, zm, Hx, Hy, Hz, deflmark);

	counti = 1; countj = 1; countk = 1;
	for (jj=0;jj<nnc;jj++) {
		knl = jj * nsd;
		i = counti-1; j = countj-1; k = countk-1;

		xcoord = i * Hx;
		ycoord = j * Hy;
		zcoord = k * Hz;

		//always compute displacement with respect to the initial undeformed Cartesian grid
		dispc[knl  ]=deflmark[knl  ]-xcoord;
		dispc[knl+1]=deflmark[knl+1]-ycoord;
		dispc[knl+2]=deflmark[knl+2]-zcoord;

		counti++;
		if (counti==xm+1){ countj++; counti=1; }
		if (countj==ym+1){ counti=1; countj=1; countk++;}
	}

	free(deflmark);
} //end function CumulateDisplacement
