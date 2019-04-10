///////////////////////////////////////////////////////////////////////////////////////
// penalized_neumann_mat.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include "function.h"


// Version 1.1 of STSolver by F. Abraham and G. Biros 
// August 2005


extern void penalized_neumann_boundary(); // ReadFiles.c


#undef __FUNCT__
#define __FUNCT__ "Penalized_Neumann_Contribution_mat"
int Penalized_Neumann_Contribution_mat(stsDMMG dmmg,Mat jac)
{
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	DA da = (DA)dmmg->dm;
#else
	DM da = dmmg->dm;
#endif
	MatStencil   *row, *col;
	PetscScalar h;
	int count, counti, countj, countk, ie;
	int ierr,i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,ie2,xm1,ym1;
	int index1, ie1,inl,JNL,knl,iec;
	int jacobian_flag, iflag, iflag2;

	PetscInt *ielist, *ixlist, *iylist, *izlist;
	PetscReal **sh, *v;
	PetscReal invpenalizedneumanneps;

	PetscReal xq2, yq2, zq2, gx, gy, gz, shxi, shyi, shzi, Nainl;
	PetscReal shxj, shyj, NaJNL, shzj;
	PetscReal xt1, yt1, zt1;
	int ix, iy, iz;
	PetscReal Hxd, **coordinl, domd;
	double shxi1, shyi1,factor;
	int ineup, listsize;

#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	ierr = DAGetInfo((DA)dmmg->dm,0,&mx,&my,&mz,0,0,0,0,0,0,0); CHKERRQ(ierr);
#else
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	ierr = DMDAGetInfo(dmmg->dm,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
#endif
	PetscFunctionBegin;
	invpenalizedneumanneps = 1.0/penalizedneumanneps;

	xm1 = xm -1;
	ym1 = ym -1;
	ne = (my-1) * (mx-1);
	if (nsd == 3) ne = ne * (mz-1);
	nnc = (xm) * (ym) ;
	if (nsd == 3) nnc = nnc * (zm);
	Hx = Lx/ (PetscReal)(mx-1); Hy = Ly/ (PetscReal)(my-1); Hz = Lz/ (PetscReal)(mz-1);
	mxg = mx; myg = my; mzg = mz;
	if (nsd == 2) Hz = 0.0;

	if (nsd == 2) listsize = 9;
	if (nsd == 3) listsize = 27;

	ierr = PetscMalloc(nen*sizeof(PetscReal),&sh); CHKERRQ(ierr);
	ierr = PetscMalloc(blksize*blksize*sizeof(PetscReal),&v); CHKERRQ(ierr);
	ierr = PetscMalloc(blksize*sizeof(MatStencil),&row); CHKERRQ(ierr);
	ierr = PetscMalloc(blksize*sizeof(MatStencil),&col); CHKERRQ(ierr);

	ierr = PetscMalloc(nen*sizeof(PetscReal),&coordinl); CHKERRQ(ierr);
	for (j=0;j<nen;j++) {
		ierr = PetscMalloc(nsd*sizeof(PetscReal),&coordinl[j]); CHKERRQ(ierr);
		ierr = PetscMalloc(nen*sizeof(PetscReal),&sh[j]); CHKERRQ(ierr);
	}
	analnode(sh);

	if (nsd == 2) {
		ierr = PetscMalloc(9*sizeof(int),&ielist);
		ierr = PetscMalloc(9*sizeof(int),&ixlist);
		ierr = PetscMalloc(9*sizeof(int),&iylist);
		ierr = PetscMalloc(9*sizeof(int),&izlist);
	} else {
		ierr = PetscMalloc(27*sizeof(int),&ielist);
		ierr = PetscMalloc(27*sizeof(int),&ixlist);
		ierr = PetscMalloc(27*sizeof(int),&iylist);
		ierr = PetscMalloc(27*sizeof(int),&izlist);
	}

	h = 0.0;
	counti=1; countj=1; countk=1;
	i = xs; j = ys; k = zs;
	jacobian_flag =1;

	h = 0.0;
	i = xs; j = ys; k = zs;
	counti=1; countj=1; countk=1;
	Hxd = Hx * Hx + Hy * Hy + Hz * Hz;
	Hxd = sqrt(Hxd);
	if (nsd == 2) {
		domd = Hxd;
	} else domd = Hxd * Hxd;
	penalized_neumann_boundary();
	for (ineup=0;ineup< npenneupoints;ineup++)
	{
		xq2 = penneupts[ineup][0];
		yq2 = penneupts[ineup][1];
		if (nsd == 3) {
			zq2 = penneupts[ineup][2];
			iz = (int)(zq2/Hz);
		}
		if (normal_based_dirichlet_bc == 1){
			gx  = normal_vector[ineup][0];
			gy  = normal_vector[ineup][1];
			if (nsd == 3) gz  = normal_vector[ineup][2];
		}

		ix = (int)(xq2/Hx);
		iy = (int)(yq2/Hy);
		if (nsd == 2) {
			i = ix-1;j =iy-1;k = 0;
			ixlist[0] = i; iylist[0] = j; 
			ie = i + j * xm1 + k * xm1 * ym1;
			ielist[0] = ie;
			ixlist[1] = i+1; iylist[1] = j; izlist[1] = k;
			ielist[1] = ie+1;
			ixlist[2] = i+2; iylist[2] = j; izlist[2] = k;
			ielist[2] = ie+2;
			ie = ie + xm1;
			ixlist[3] = i; iylist[3] = j+1; izlist[3] = k;
			ielist[3] = ie;
			ixlist[4] = i+1; iylist[4] = j+1; izlist[4] = k;
			ielist[4] = ie+1;
			ixlist[5] = i+2; iylist[5] = j+1; izlist[5] = k;
			ielist[5] = ie+2;
			ie = ie + xm1;
			ixlist[6] = i; iylist[6] = j+2; izlist[6] = k;
			ielist[6] = ie;
			ixlist[7] = i+1; iylist[7] = j+2; izlist[7] = k;
			ielist[7] = ie+1;
			ixlist[8] = i+2; iylist[8] = j+2; izlist[8] = k;
			ielist[8] = ie+2;
		}

		if (nsd == 3) {
			i = ix-1;j =iy-1;k = iz-1;
			ixlist[0] = i; iylist[0] = j; izlist[0] = k;
			ie = i + j * xm1 + k * xm1 * ym1;
			ielist[0] = ie;
			ixlist[1] = i+1; iylist[1] = j; izlist[1] = k;
			ielist[1] = ie+1;
			ixlist[2] = i+2; iylist[2] = j; izlist[2] = k;
			ielist[2] = ie+2;
			ie = ie + xm1;
			ixlist[3] = i; iylist[3] = j+1; izlist[3] = k;
			ielist[3] = ie;
			ixlist[4] = i+1; iylist[4] = j+1; izlist[4] = k;
			ielist[4] = ie+1;
			ixlist[5] = i+2; iylist[5] = j+1; izlist[5] = k;
			ielist[5] = ie+2;
			ie = ie + xm1;
			ixlist[6] = i; iylist[6] = j+2; izlist[6] = k;
			ielist[6] = ie;
			ixlist[7] = i; iylist[7] = j+2; izlist[7] = k;
			ielist[7] = ie+1;
			ixlist[8] = i+1; iylist[8] = j+2; izlist[8] = k;
			ielist[8] = ie+2;

			ie = i + j * xm1 + (k+1) * xm1 * ym1;
			ixlist[9] = i; iylist[9] = j; izlist[9] = k+1;
			ielist[9] = ie;
			ixlist[10] = i+1; iylist[10] = j; izlist[10] = k+1;
			ielist[10] = ie+1;
			ixlist[11] = i+2; iylist[11] = j; izlist[11] = k+1;
			ielist[11] = ie+2;
			ie = ie + xm1;
			ixlist[12] = i; iylist[12] = j+1; izlist[12] = k+1;
			ielist[12] = ie;
			ixlist[13] = i+1; iylist[13] = j+1; izlist[13] = k+1;
			ielist[13] = ie+1;
			ixlist[14] = i+2; iylist[14] = j+1; izlist[14] = k+1;
			ielist[14] = ie+2;
			ie = ie + xm1;
			ixlist[15] = i; iylist[15] = j+2; izlist[15] = k+1;
			ielist[15] = ie;
			ixlist[16] = i+1; iylist[16] = j+2; izlist[16] = k+1;
			ielist[16] = ie+1;
			ixlist[17] = i+2; iylist[17] = j+2; izlist[17] = k+1;
			ielist[17] = ie+2;

			ie = i + j * xm1 + (k+2) * xm1 * ym1;
			ixlist[18] = i; iylist[18] = j; izlist[18] = k+2;
			ielist[18] = ie;
			ixlist[19] = i+1; iylist[10] = j; izlist[19] = k+2;
			ielist[19] = ie+1;
			ixlist[20] = i+2; iylist[20] = j; izlist[20] = k+2;
			ielist[20] = ie+2;
			ie = ie + xm1;
			ixlist[21] = i; iylist[21] = j+1; izlist[21] = k+2;
			ielist[21] = ie;
			ixlist[22] = i+1; iylist[22] = j+1; izlist[22] = k+2;
			ielist[22] = ie+1;
			ixlist[23] = i+2; iylist[23] = j+1; izlist[23] = k+2;
			ielist[23] = ie+2;
			ie = ie + xm1;
			ixlist[24] = i; iylist[24] = j+2; izlist[24] = k+2;
			ielist[24] = ie;
			ixlist[25] = i+1; iylist[25] = j+2; izlist[25] = k+2;
			ielist[25] = ie+1;
			ixlist[26] = i+2; iylist[26] = j+2; izlist[26] = k+2;
			ielist[26] = ie+2;
		}

		/* Find the list of Elements which belong to or touch this point */

		/* Form list of neighbouring elements 8 for 2D, 36 for 3D*/
		/* For each of these elements make sure that this point is within ot touches it */

		for (iec=0;iec<listsize;iec++)
		{
			ie = ielist[iec];

			ie1 = ie * nen;
			ie2 = ie1* ndf;
			i = ixlist[iec]; j = iylist[iec]; k = izlist[iec];

			for (count=0;count<blksize*blksize;count++)
				v[count]=0;
			iflag = 0;
			if (nsd == 3 && ixlist[iec] > 0 && ixlist[iec] < mx && iylist[iec] > 0 && iylist[iec]< my && izlist[iec] > 0 && izlist[iec] < mz){
				iflag = 1;
			}
			if (nsd == 2 && ixlist[iec] > 0 && ixlist[iec] < mx && iylist[iec] > 0 && iylist[iec]< my ){
				iflag = 1;
			}
			if (iflag == 1) {
				rowcolindex(i,j,k,row,col,jacobian_flag);

				coordinl[0][0] = i * Hx;	
				coordinl[0][1] = j * Hy;	

				coordinl[1][0] = (i+1) * Hx;	
				coordinl[1][1] = j * Hy;	

				coordinl[2][0] = (i+1) * Hx;	
				coordinl[2][1] = (j+1) * Hy;	

				coordinl[3][0] = i * Hx;	
				coordinl[3][1] = (j+1) * Hy;	

				if ( nsd == 3){
					coordinl[4][0] = i * Hx;	
					coordinl[4][1] = j * Hy;	

					coordinl[5][0] = (i+1) * Hx;	
					coordinl[5][1] = j * Hy;	

					coordinl[6][0] = (i+1) * Hx;	
					coordinl[6][1] = (j+1) * Hy;	

					coordinl[7][0] = i * Hx;	
					coordinl[7][1] = (j+1) * Hy;	

					coordinl[0][2] = k * Hz;	
					coordinl[1][2] = k * Hz;	
					coordinl[2][2] = k * Hz;	
					coordinl[3][2] = k * Hz;	

					coordinl[4][2] = (k+1) * Hz;	
					coordinl[5][2] = (k+1) * Hz;	
					coordinl[6][2] = (k+1) * Hz;	
					coordinl[7][2] = (k+1) * Hz;	
				}
				xt1 = xq2 -i*Hx -Hx/2;yt1=yq2 -j*Hy-Hy/2; zt1 = zq2-k*Hz-Hz/2;
				xt1 = xt1/Hx * 2; yt1 = yt1/Hy * 2; zt1 = zt1/Hz * 2;
				iflag2 = 0;
				if (nsd == 3 &&( xt1 > 1 || xt1 < -1 ||yt1 >1 || yt1 < -1 ||zt1>1||zt1<-1)) iflag2 = 1;
				if (nsd == 2 &&( xt1 > 1 || xt1 < -1 ||yt1 >1 || yt1 < -1 )) iflag2 = 1;

				if (iflag2 == 0) {
					for (inl=0;inl < nen;inl++) {
						knl = inl * ndf;
						shxi = sh[inl][xsd];
						shyi = sh[inl][ysd];
						if (nsd == 3) {
							shzi = sh[inl][zsd];
						}
						if (nsd == 3) {
							Nainl = 0.125 * (1 + shxi * xt1) * (1 + shyi * yt1) * (1 + shzi * zt1);
						}else{
							Nainl = 0.25 * (1 + shxi * xt1) * (1 + shyi * yt1);
						}
						for (JNL=0;JNL<nen;JNL++) {
							index1 = (inl * blksize + JNL)* ndf;
							shxj = sh[JNL][xsd];
							shyj = sh[JNL][ysd];
							if (nsd == 3){
								//shzi = sh[JNL][zsd];
								shzj = sh[JNL][zsd];
								NaJNL = 0.125 * (1 + shxj * xt1) * (1 + shyj * yt1) * (1 + shzj * zt1);
							} else {
								NaJNL = 0.25 * (1 + shxj * xt1) * (1 + shyj * yt1);
							}

							shxi1 = 0.5 * (1+shxi * xt1);
							shyi1 = 0.5 * (1+shyi * yt1);

							factor = 1.0;
							v[index1]     = factor * Nainl * NaJNL* Hpenalizedneumann * invpenalizedneumanneps;
							if (normal_based_dirichlet_bc == 1) v[index1] = v[index1] * gx;
							index1 = index1 + blksize;
							v[index1+1]   = factor * Nainl * NaJNL* Hpenalizedneumann * invpenalizedneumanneps; 
							if (normal_based_dirichlet_bc == 1) v[index1+1] = v[index1+1] * gy;
							index1 = index1 + blksize;
							if (nsd == 3) {
								v[index1+2]   = factor * Nainl * NaJNL * Hpenalizedneumann * invpenalizedneumanneps;
								if (normal_based_dirichlet_bc == 1) v[index1+2] = v[index1+2] * gz;
							}

						}
					}
				}
				ierr = MatSetValuesStencil(jac,blksize,row,blksize,col,v,ADD_VALUES);CHKERRQ(ierr);
			}
		}
	}
	
	ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatSetBlockSize(jac, nsd); CHKERRQ(ierr);

	ierr = PetscFree(sh); CHKERRQ(ierr);
	//ierr = PetscFree(idx); CHKERRQ(ierr);

	PetscFree(v);
	PetscFree(row);
	PetscFree(col);
	PetscFree(ixlist);
	PetscFree(iylist);
	PetscFree(izlist);
	PetscFree(ielist);

	PetscFunctionReturn(0);
}

