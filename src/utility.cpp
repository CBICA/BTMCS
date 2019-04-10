///////////////////////////////////////////////////////////////////////////////////////
// utility.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include "function.h"


#ifndef rint
double rint(double x)
{
	//middle value point test
	if (ceil(x+0.5) == floor(x+0.5))
	{
		int a = (int)ceil(x);
		if (a%2 == 0) { return ceil(x); }
		else { return floor(x); }
	}
	else return floor(x+0.5);
}
#endif


void dirichlet_identifier(int *dirichlet_flag, MatStencil *row, int mx, int my, int mz)
{
	int inl, knl;

	for (inl = 0; inl < nen; inl++) {
		knl = inl *ndf;
		dirichlet_flag[knl] = 1;
		dirichlet_flag[knl+1] = 1;
		if (nsd == 3) dirichlet_flag[knl+2] = 1;
		if (boundaryid[0].idf == 0) {if (row[knl].j == 0) {dirichlet_flag[knl]   = boundaryid[0].idf;}}
		if (boundaryid[0].jdf == 0) {if (row[knl].j == 0) dirichlet_flag[knl+1] = boundaryid[0].jdf;}

		if (boundaryid[2].idf == 0) {if (row[knl].j == my-1) dirichlet_flag[knl]   = boundaryid[2].idf;}
		if (boundaryid[2].jdf == 0) {if (row[knl].j == my-1) dirichlet_flag[knl+1] = boundaryid[2].jdf;}

		if (boundaryid[3].idf == 0) {if (row[knl].i == 0) dirichlet_flag[knl]   = boundaryid[3].idf;}
		if (boundaryid[3].jdf == 0) {if (row[knl].i == 0) dirichlet_flag[knl+1] = boundaryid[3].jdf;}

		if (boundaryid[1].idf == 0) {if (row[knl].i == mx-1) dirichlet_flag[knl]    = boundaryid[1].idf;}
		if (boundaryid[1].jdf == 0) {if (row[knl].i == mx-1) dirichlet_flag[knl+1]  = boundaryid[1].jdf;}

		if (nsd== 3){
			if (boundaryid[4].idf == 0) {if (row[knl].k == 0)    dirichlet_flag[knl]   = boundaryid[4].idf;}
			if (boundaryid[4].jdf == 0) {if (row[knl].k == 0)    dirichlet_flag[knl+1] = boundaryid[4].jdf;}
			if (boundaryid[5].idf == 0) {if (row[knl].k == mz-1) dirichlet_flag[knl]   = boundaryid[5].idf;}
			if (boundaryid[5].jdf == 0) {if (row[knl].k == mz-1) dirichlet_flag[knl+1] = boundaryid[5].jdf;}
			if (boundaryid[0].kdf == 0) {if (row[knl].j == 0)    dirichlet_flag[knl+2] = boundaryid[0].kdf;}
			if (boundaryid[1].kdf == 0) {if (row[knl].i == mx-1)    dirichlet_flag[knl+2] = boundaryid[1].kdf;}
			if (boundaryid[2].kdf == 0) {if (row[knl].j == my-1) dirichlet_flag[knl+2] = boundaryid[2].kdf;}
			if (boundaryid[3].kdf == 0) {if (row[knl].i == 0) dirichlet_flag[knl+2] = boundaryid[3].kdf;}
			if (boundaryid[4].kdf == 0) {if (row[knl].k == 0)    dirichlet_flag[knl+2] = boundaryid[4].kdf;}
			if (boundaryid[5].kdf == 0) {if (row[knl].k == mz-1) dirichlet_flag[knl+2] = boundaryid[5].kdf;}
		}
	}
}

void hacks()
{
	int inc, ie, i, inl, j, ie1;
	PetscMalloc(nnc*sizeof(PetscReal),&hack1n);
	PetscMalloc(ne*nen*sizeof(PetscReal),&hack1e);

	for (inc=0; inc<nnc;inc++){ hack1n[inc] =0.0; }
	for (ie=0; ie<ne*nen;ie++){ hack1e[ie] =1.0; }
	for (ie=0; ie<ne;ie++){
		ie1 = ie * nen; 
		for (inl=0;inl<nen;inl++){
			inc = ien[ie1 + inl];
			hack1n[inc] +=	hack1e[ie1+inl];
		}		
	}
	for (ie=0; ie<ne*nen;ie++){ hack1e[ie] =0.0; } 
	for (i=0;i<ne;i++){
		ie =  i * nen;
		for (j=0;j<nen;j++){
			inc = ien[ie+j];
			hack1e[ie+j] = hack1n[inc];
		}
	}
}

void createien(int xm, int ym)
{
	int counti, countj, countk, i,j,k, ie, inc;
	int rowi, rowj,rowk;

	counti=1; countj=1; countk=1;

	PetscMalloc(ne*nen*sizeof(int),&ien);
	for (ie=0;ie<ne;ie++){
		i = counti-1; j = countj-1; k = countk-1;
		inc = ie*nen;

		rowi = i; rowj =j;rowk =k;
		ien[inc] = node_index(rowi,rowj,rowk, xm, ym);
		rowi = i+1; rowj =j;rowk =k;
		ien[inc+1] = node_index(rowi,rowj,rowk,xm,ym);
		rowi = i+1; rowj =j+1;rowk =k;
		ien[inc+2] = node_index(rowi,rowj,rowk,xm,ym);
		rowi = i; rowj =j+1;rowk =k;
		ien[inc+3] = node_index(rowi,rowj,rowk,xm,ym);

		if (nsd ==3){
			rowi = i; rowj =j;rowk =k+1;
			ien[inc+4] = node_index(rowi,rowj,rowk,xm,ym);
			rowi = i+1; rowj =j;rowk =k+1;
			ien[inc+5] = node_index(rowi,rowj,rowk,xm,ym);
			rowi = i+1; rowj =j+1;rowk =k+1;
			ien[inc+6] = node_index(rowi,rowj,rowk,xm,ym);
			rowi = i; rowj =j+1;rowk =k+1;
			ien[inc+7] = node_index(rowi,rowj,rowk,xm,ym);
		}
		counti++;
		if (counti==xm){ countj++; counti=1; }
		if (countj==ym){counti=1; countj=1; countk++; }
	}
}

#undef __FUNCT__
#define __FUNCT__ "node_index"
int node_index(int rowi, int rowj, int rowk, int xd, int yd)
{
	int index1;

	index1 = rowj * xd+ rowi;
	if (ndf ==3) index1 += rowk*(xd)*(yd);
	return index1;
}

void gatherreal1(double *xr, double *yr)
{
	int i,j, inl,ie;

	for (i=0;i<ne;i++) {
		ie =  i * nen;
		for (j=0;j<nen;j++){
			inl = ien[ie+j];
			yr[ie+j] = xr[inl];
		}
	}
}

void gatherreal(double *xr, double *yr)
{
	int i,j,inl,inl1, ie1,ie;

	for (i=0;i<ne;i++) {
		ie1 = i * nen * nsd;
		ie  = i * nen;
		for (j=0;j<nen;j++){
			inl  = ien[ie+j];
			inl1 = inl * nsd;
			yr[ie1]   = xr[inl1];
			yr[ie1+1] = xr[inl1+1];
			if (nsd==3) yr[ie1+2] = xr[inl1+2];
			ie1 += nsd;
		}
	}
}

void CartesianDirichlet(double *xr, int xm, int ym, int zm)
{
	int inc, counti, countj, countk, i, j, k, ind;

	if (titlecase == 1)
	{
		counti=1; countj=1; countk=1;
		for (inc =0;inc<nnc;inc++)
		{
			ind = inc * nsd + 2;
			i = counti-1; j = countj-1; k = countk-1;
			if (k == zm-1){
				if (tchoice ==1)      xr[ind] = -0.4;
				else if (tchoice ==2) xr[ind] = -1.0;
				else if (tchoice ==3) xr[ind] = -1.46;
			}
			counti++;
			if (counti==xm+1){ countj++; counti=1; }
			if (countj==ym+1){ counti=1; countj=1; countk++;}
		}
	} else {
		counti=1; countj=1; countk=1;
		for (inc =0;inc<nnc;inc++){
			i = counti-1; j = countj-1; k = countk-1;
			ind = inc * nsd;
			xr[ind] = 0.0;
			xr[ind+1] = 0.0;
			if ( nsd ==3)xr[ind+2] = 0.0;
			counti++;
			if (counti==xm+1){ countj++; counti=1; }
			if (countj==ym+1){ counti=1; countj=1; countk++;}
		}
	}
}

void inhomogeneous_material(double *lambda, double *mu, int ne, int xm, int ym, double Hx, double Hy, double Hz)
{
	int ie, ie1, ix, jy, kz, i, j, k, counti, countj, countk, mxfm1;
	double xmult, xc, yc, zc, lamnode, munode;
	int notfinelevel = 1;

	xmult = (mxf-1) * (myf-1);
	mxfm1 = mxf -1 ;

	notfinelevel = 1;
	if (fabs(Hx-Hxf)<eps5 && fabs(Hy-Hyf)<eps5 && fabs(Hz-Hzf) <eps5) notfinelevel = 0;

	if (material_projection_required == 1 && notfinelevel ==1)
	{
		counti=1; countj=1; countk=1;

		for (ie=0;ie < ne;ie++)
		{
			ie1 = ie * nen * nsd;
			i = counti-1; j = countj-1; k = countk-1;

			lamnode = 0.0; munode = 0.0;
			xc = i*Hx; yc = j*Hy; zc = k*Hz;
			ix =(int)(xc/Hxf + 0.5); jy =(int)(yc/Hyf + 0.5); kz = (int)(zc/Hzf + 0.5); 
			ie1 = (int)(kz * xmult + jy * mxfm1 + ix); 
			munode = munode + muf[ie1];
			lamnode = lamnode + lambdaf[ie1];

			//if(ie==30221) printf("ie1 lambdaf[ie1] %d %g\n",ie1,lambdaf[ie1]);

			xc = (i+1) * Hx; yc = j * Hy; zc = k * Hz;
			ix = (int)(xc/Hxf + 0.5)-1; jy = (int)(yc/Hyf + 0.5); kz =(int)(zc/Hzf + 0.5);
			ie1 = (int)(kz * xmult + jy * mxfm1 + ix);
			munode = munode + muf[ie1];
			lamnode = lamnode + lambdaf[ie1];

			//if(ie==30221) printf("ie1 lambdaf[ie1] %d %g\n",ie1,lambdaf[ie1]);

			xc = (i+1) * Hx; yc = (j+1) * Hy; zc = k * Hz;
			ix = (int)(xc/Hxf + 0.5)-1; jy = (int)(yc/Hyf + 0.5)-1; kz = (int)(zc/Hzf + 0.5);
			ie1 = (int)(kz * xmult + jy * mxfm1 + ix);
			munode = munode + muf[ie1];
			lamnode = lamnode + lambdaf[ie1];

			//if(ie==30221) printf("ie1 lambdaf[ie1] %d %g\n",ie1,lambdaf[ie1]);

			xc = i * Hx; yc = (j+1) * Hy; zc = k * Hz;
			ix = (int)(xc/Hxf + 0.5); jy = (int)(yc/Hyf + 0.5)-1; kz = (int)(zc/Hzf + 0.5);
			ie1 = (int)(kz * xmult + jy * mxfm1 + ix);
			munode = munode + muf[ie1];
			lamnode = lamnode + lambdaf[ie1];

			//if(ie==30221) printf("ie1 lambdaf[ie1] %d %g\n",ie1,lambdaf[ie1]);

			xc = i * Hx; yc = j * Hy; zc = (k+1) * Hz;
			ix = (int)(xc/Hxf + 0.5); jy = (int)(yc/Hyf + 0.5); kz = (int)(zc/Hzf + 0.5)-1;
			ie1 = (int)(kz * xmult + jy * mxfm1 + ix);
			munode = munode + muf[ie1];
			lamnode = lamnode + lambdaf[ie1];

			// if(ie==30221) printf("ie1 lambdaf[ie1] %d %g\n",ie1,lambdaf[ie1]);

			xc = (i+1) * Hx; yc = j * Hy; zc = (k+1) * Hz;
			ix = (int)(xc/Hxf + 0.5)-1; jy = (int)(yc/Hyf + 0.5); kz = (int)(zc/Hzf + 0.5)-1;
			ie1 = (int)(kz * xmult + jy * mxfm1 + ix);
			munode = munode + muf[ie1];
			lamnode = lamnode + lambdaf[ie1];

			//if(ie==30221) printf("ie1 lambdaf[ie1] %d %g\n",ie1,lambdaf[ie1]);

			xc = (i+1) * Hx; yc = (j+1) * Hy; zc = (k+1) * Hz;
			ix = (int)(xc/Hxf + 0.5)-1;jy = (int)(yc/Hyf + 0.5)-1; kz = (int)(zc/Hzf + 0.5)-1;
			ie1 = (int)(kz * xmult + jy * mxfm1 + ix);
			munode = munode + muf[ie1];
			lamnode = lamnode + lambdaf[ie1];

			// if(ie==30221) printf("ie1 lambdaf[ie1] %d %g\n",ie1,lambdaf[ie1]);

			xc = i * Hx; yc = (j+1) * Hy; zc = (k+1) * Hz;
			ix = (int)(xc/Hxf + 0.5); jy = (int)(yc/Hyf + 0.5)-1; kz = (int)(zc/Hzf + 0.5)-1;
			ie1 = (int)(kz * xmult + jy * mxfm1 + ix);
			munode = munode + muf[ie1];
			lamnode = lamnode + lambdaf[ie1];

			//if(ie==30221) printf("ie1 lambdaf[ie1] %d %g\n",ie1,lambdaf[ie1]);

			mu[ie] = munode/nen;
			lambda[ie] = lamnode/nen;

			//if(lambdaf[ie1]<651724) printf("utility ie1 lambdaf %d %g\n", ie1, lambdaf[ie1]);
			//if(lambda[ie]<651724) printf("utility ie lambda %d %g\n", ie, lambda[ie]);

			//if(ie==30221) printf("ie lambda[ie] %d %g\n",ie,lambda[ie]);

			counti++;
			if (counti==xm){ countj++; counti=1; }
			if (countj==ym){ counti=1; countj=1; countk++; }
		}
	} else {
		for (ie=0;ie < ne;ie++){
			mu[ie] = muf[ie];
			lambda[ie] = lambdaf[ie];

			//if(lambdaf[ie]<651724) printf("utility ie lambda %d %g\n", ie,lambda[ie]);
		}
	}
}
