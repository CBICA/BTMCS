///////////////////////////////////////////////////////////////////////////////////////
// matrixfree.c
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


static stsData arr[100]; //static int jaccnt=0;
static PetscErrorCode ComputeMaterialProperties(stsData *data);

PetscScalar *lambdamf=0,*mumf=0,*vbacmf=0;

//#define MATVEC_PROF

//#define node_index3(rowi,rowj,rowk,xd,yd) ( (rowj)*(xd) + (rowi) + (rowk)*(xd)*(yd))


#undef __FUNCT__
#define __FUNCT__ "matrixFreeMatVec"
int matrixFreeMatVec(stsData *data, Vec Vecin, Vec Vecout)
{
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	DA da = (DA)(data->dmmg->dm);
#else
	DM da = data->dmmg->dm;
#endif
	int ie1,ierr,i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
	int counti, countj, countk, ie,inl,jnl,knl;
	int index1,index2, jacobian_flag, *inc1,inc;
	int rowi, rowj, rowk, index3, index4;
	int iflagx, iflagy, iflagz;

	double shyshy, shxshx, shzshz, tmp1;
	double shdshd, sh0sh0;
	double shxshy,shxshz, shyshx,shyshz, shzshy,shzshx;
	double h;

	PetscReal lamnode, munode, lampmunode;
	int Nfebs;

	PetscReal *vinarray,*voutarray;
	PetscReal *vloc,*vloc1;

#ifdef MATVEC_PROF 
	PetscLogEventBegin(MatVecFlop,0,0,0,0);
#endif

#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0);CHKERRQ(ierr);  
	ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
#else
	ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);  
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
#endif

	inhomogenity_constant = 1;

	nnc = (xm) * (ym) ;
	if (nsd == 3) nnc = nnc * (zm);
	ne = (my-1) * (mx-1);
	if (nsd ==3) ne = ne * (mz-1);

	Hx = Lx/ (PetscReal)(mx-1); Hy = Ly/ (PetscReal)(my-1); Hz = Lz/ (PetscReal)(mz-1);

	//printf(" xm ym Hx Hy %d %d %g %g\n",my,mx,Hx, Hy);

	VecGetSize(data->lambdavec,&Nfebs);
	h = 0.0;
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	ierr = VecSet(&h, Vecout);CHKERRQ(ierr);
#else
	ierr = VecSet(Vecout, h); CHKERRQ(ierr);
#endif

	ierr = VecGetArray1d(Vecin,nnc*nsd,0, &vinarray);
	ierr = VecGetArray1d(Vecout,nnc*nsd,0,&voutarray);

	ierr = VecGetArray1d(data->muvec,ne,0,&mumf);
	ierr = VecGetArray1d(data->lambdavec,ne,0,&lambdamf);
	ierr = VecGetArray1d(data->vbac,nnc*nsd,0,&vbacmf);
	//ierr = VecGetArray1d(data->vjacprec,nnc*nsd,0,&vjacprec);

	//ierr=VecView(lambdamf,PETSC_VIEWER_STDOUT_WORLD);
	PetscMalloc(nen*nsd*sizeof(PetscReal),&vloc);
	PetscMalloc(nen*nsd*sizeof(PetscReal),&vloc1);
	PetscMalloc(nen*sizeof(int),&inc1);

	counti=1; countj=1; countk=1;
	for (inc =0;inc<nnc;inc++){
		i = counti-1; j = countj-1; k = countk-1;
		knl = inc * nsd;
		iflagx = 0;iflagy = 0;iflagz = 0;

		if (counti==1 && boundaryid[3].idf == 0) iflagx=1;
		if (counti==1 && boundaryid[3].jdf == 0) iflagy=1;
		if (nsd ==3 && counti==1 && boundaryid[3].kdf == 0) iflagz=1;

		if (counti==mx && boundaryid[1].idf == 0) iflagx=1;
		if (counti==mx && boundaryid[1].jdf == 0) iflagy=1;
		if (nsd ==3 && counti==mx && boundaryid[1].kdf == 0) iflagz=1;

		if (countj==1 && boundaryid[0].idf == 0) iflagx=1;
		if (countj==1 && boundaryid[0].jdf == 0) iflagy=1;
		if (nsd ==3 && countj==1 && boundaryid[0].kdf == 0) iflagz=1;

		if (countj==my && boundaryid[2].idf == 0) iflagx=1;
		if (countj==my && boundaryid[2].jdf == 0) iflagy=1;
		if (nsd == 3 && countj==my && boundaryid[2].kdf == 0) iflagz=1;

		if ( nsd ==3){
			if (countk==1 && boundaryid[4].idf == 0) iflagx=1;
			if (countk==1 && boundaryid[4].jdf == 0) iflagy=1;
			if (countk==1 && boundaryid[4].kdf == 0) iflagz=1;

			if (countk==mz && boundaryid[5].idf == 0) iflagx=1;
			if (countk==mz && boundaryid[5].jdf == 0) iflagy=1;
			if (countk==mz && boundaryid[5].kdf == 0) iflagz=1;
		}

		if (iflagx == 0) vbacmf[knl]   = vinarray[knl];
		else vbacmf[knl]   = 0.0;

		if (iflagy == 0) vbacmf[knl+1]   = vinarray[knl+1];
		else vbacmf[knl+1] = 0.0;

		if (iflagz == 0) vbacmf[knl+2]   = vinarray[knl+2];
		else vbacmf[knl+2] = 0.0;

		counti++;
		if (counti==mx+1){ countj++; counti=1; }
		if (countj==my+1){ counti=1; countj=1; countk++; }
	}

	jacobian_flag = 1;
	counti=1; countj=1; countk=1;
	i = xs; j = ys; k = zs;
	for (ie=0;ie<ne;ie++){
		ie1 = ie * nen;
		i = counti-1; j = countj-1; k = countk-1;

		rowi = i; rowj =j;rowk =k;
		inc1[0] = node_index3(rowi,rowj,rowk, xm, ym);
		rowi = i+1; rowj =j;rowk =k;
		inc1[1] = node_index3(rowi,rowj,rowk,xm,ym);
		rowi = i+1; rowj =j+1;rowk =k;
		inc1[2] = node_index3(rowi,rowj,rowk,xm,ym);
		rowi = i; rowj =j+1;rowk =k;
		inc1[3] = node_index3(rowi,rowj,rowk,xm,ym);

		rowi = i; rowj =j;rowk =k+1;
		inc1[4] = node_index3(rowi,rowj,rowk,xm,ym);
		rowi = i+1; rowj =j;rowk =k+1;
		inc1[5] = node_index3(rowi,rowj,rowk,xm,ym);
		rowi = i+1; rowj =j+1;rowk =k+1;
		inc1[6] = node_index3(rowi,rowj,rowk,xm,ym);
		rowi = i; rowj =j+1;rowk =k+1;
		inc1[7] = node_index3(rowi,rowj,rowk,xm,ym);

		for (inl=0;inl<nen;inl++){
			inc = inc1[inl];
			inc = inc *nsd;
			knl = inl * nsd;
			vloc[knl]   = vbacmf[inc];
			vloc[knl+1] = vbacmf[inc+1];
			vloc[knl+2] = vbacmf[inc+2];
			vloc1[knl]   = 0.0;
			vloc1[knl+1] = 0.0;
			vloc1[knl+2] = 0.0;
		}

		lamnode     = lambdamf[ie];
		munode      = mumf[ie];

		//	printf("matrixfree lamnode %g\n", lamnode);
		//	printf("matrixfree munode %g\n", munode);

		lampmunode  = lamnode + munode;
#ifdef MATVEC_PROF 
		PetscLogFlops(1);
#endif
		sh0sh0 = 0.0;
		for (inl=0;inl<nen;inl++){
			index2 = inl * nsd;
			index3 = index2 + 1;
			index4 = index2 + 2;
			for (jnl=0;jnl<nen;jnl++){
				index1 = jnl * nsd;

				shxshx = data->vshxshx[inl][jnl];
				shxshy = data->vshxshy[inl][jnl];
				shxshz = data->vshxshz[inl][jnl];
				shyshy = data->vshyshy[inl][jnl];
				shyshz = data->vshyshz[inl][jnl];
				shyshx = data->vshyshx[inl][jnl];
				shzshx = data->vshzshx[inl][jnl];
				shzshy = data->vshzshy[inl][jnl];
				shzshz = data->vshzshz[inl][jnl];
				shdshd = data->vshdshd[inl][jnl];

				tmp1   = (shdshd*munode + shxshx * lampmunode);
				tmp1   += sh0sh0;

				vloc1[index2]   += tmp1 * vloc[index1];
				tmp1 = (shyshx * munode + shxshy * lamnode);
				vloc1[index2]   += tmp1 * vloc[index1+1];
				tmp1 = (shzshx * munode + shxshz * lamnode);
				vloc1[index2]   += tmp1 * vloc[index1+2];

				tmp1   = (shxshy * munode + shyshx * lamnode);
				vloc1[index3]   += tmp1 * vloc[index1];
				tmp1 = (shdshd *munode + shyshy * lampmunode );
				tmp1 += sh0sh0;
				vloc1[index3]   += tmp1 * vloc[index1+1];
				tmp1 = (shzshy * munode + shyshz * lamnode);
				vloc1[index3]   += tmp1 * vloc[index1+2];


				tmp1 = (shxshz * munode + shzshx * lamnode);
				vloc1[index4]   += tmp1 * vloc[index1];
				tmp1 = (shyshz * munode + shzshy * lamnode);
				vloc1[index4]   += tmp1 * vloc[index1+1];
				tmp1 = (shdshd *munode + shzshz * lampmunode);
				tmp1 += sh0sh0;
				vloc1[index4]   += tmp1 * vloc[index1+2];

#ifdef MATVEC_PROF 
				PetscLogFlops(48);
#endif
			}
		}

		for (inl=0;inl<nen;inl++){
			inc = inc1[inl];
			inc = inc * nsd;
			knl = inl * nsd;
			voutarray[inc]   +=  vloc1[knl];
			voutarray[inc+1] +=  vloc1[knl+1];
			voutarray[inc+2] +=  vloc1[knl+2];

#ifdef MATVEC_PROF 
			PetscLogFlops(3);
#endif

		}
		counti++;
		if (counti==xm){ countj++; counti=1; } 
		if (countj==ym){ counti=1; countj=1; countk++; }
	}

	counti=1; countj=1; countk=1;
	for (inc =0;inc<nnc;inc++){
		i = counti-1; j = countj-1; k = countk-1;
		knl = inc * nsd;
		iflagx = 0; iflagy = 0; iflagz = 0;
		if (counti==1 && boundaryid[3].idf == 0) iflagx=1;
		if (counti==1 && boundaryid[3].jdf == 0) iflagy=1;
		if (nsd ==3 && counti==1 && boundaryid[3].kdf == 0) iflagz=1;

		if (counti==mx && boundaryid[1].idf == 0) iflagx=1;
		if (counti==mx && boundaryid[1].jdf == 0) iflagy=1;
		if (nsd ==3 && counti==mx && boundaryid[1].kdf == 0) iflagz=1;

		if (countj==1 && boundaryid[0].idf == 0) iflagx=1;
		if (countj==1 && boundaryid[0].jdf == 0) iflagy=1;
		if (nsd ==3 && countj==1 && boundaryid[0].kdf == 0) iflagz=1;

		if (countj==my && boundaryid[2].idf == 0) iflagx=1;
		if (countj==my && boundaryid[2].jdf == 0) iflagy=1;
		if (nsd ==3 && countj==my && boundaryid[2].kdf == 0) iflagz=1;

		if ( nsd ==3){
			if (countk==1 && boundaryid[4].idf == 0) iflagx=1;
			if (countk==1 && boundaryid[4].jdf == 0) iflagy=1;
			if (countk==1 && boundaryid[4].kdf == 0) iflagz=1;

			if (countk==mz && boundaryid[5].idf == 0) iflagx=1;
			if (countk==mz && boundaryid[5].jdf == 0) iflagy=1;
			if (countk==mz && boundaryid[5].kdf == 0) iflagz=1;
		}

		if (iflagx == 1){
			if (mgnlevels>1){
				voutarray[knl]   = vinarray[knl] * youngs_global;
#ifdef MATVEC_PROF 
				PetscLogFlops(1);
#endif
			}else{
				voutarray[knl]   = vinarray[knl];
			}
		}
		if (iflagy == 1){
			if (mgnlevels>1){
				voutarray[knl+1]   = vinarray[knl+1] * youngs_global;
#ifdef MATVEC_PROF 
				PetscLogFlops(1);
#endif
			}else{
				voutarray[knl+1]   = vinarray[knl+1];
			}
		}
		if (iflagz == 1){
			if (mgnlevels>1){
				voutarray[knl+2]   = vinarray[knl+2] * youngs_global;
#ifdef MATVEC_PROF 
				PetscLogFlops(1);
#endif
			}else{
				voutarray[knl+2]   = vinarray[knl+2];
			}
		}

		counti++;
		if (counti==mx+1){ countj++; counti=1; }
		if (countj==my+1){ counti=1; countj=1; countk++; }
	}

	//   if (jacobiprecflag == 1){
	//	for (inc =0;inc<nnc;inc++){
	//	    knl = inc * nsd;
	//	    voutarray[knl]     = voutarray[knl]   * vjacprec[knl];
	//	    voutarray[knl+1]   = voutarray[knl+1] * vjacprec[knl+1];
	//	    voutarray[knl+2]   = voutarray[knl+2] * vjacprec[knl+2];
	//	}
	//   }

	VecRestoreArray1d(Vecin,nnc*ndf,0,&vinarray);
	VecRestoreArray1d(Vecout,nnc*ndf,0,&voutarray);
	VecRestoreArray1d(data->muvec,ne,0,&mumf);
	VecRestoreArray1d(data->lambdavec,ne,0,&lambdamf);
	VecRestoreArray1d(data->vbac,nnc*nsd,0,&vbacmf);
	//   VecRestoreArray1d(data->vjacprec,nnc*nsd,0,&vjacprec);

	PetscFree(vloc);
	PetscFree(vloc1);
	PetscFree(inc1);

#ifdef MATVEC_PROF
	PetscLogEventEnd(MatVecFlop,0,0,0,0);
#endif

	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "CreateJacobianMatrixFree"
PetscErrorCode CreateJacobianMatrixFree(stsDMMG dmmg, Mat *jac)
{
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	DA da = (DA)dmmg->dm;
#else
	DM da = dmmg->dm;
#endif
	PetscErrorCode ierr;
	//PetscViewer myviewer;
	PetscInt  N;
	int mx, my, mz;

	//printf("Matrixfree: jaccnt _____________________________________ %d\n",jaccnt);

	assert(jaccnt<MAX_DMMG_LEVELS);
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0); CHKERRQ(ierr);
#else
	ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
#endif
	N = mx*my*mz*3;

	if(dmmg->nlevels<mgnlevels)

	{
		// create matrix-free at all levels except for the coarsest one
		shellData[jaccnt].dmmg = dmmg;
		ierr = MatCreateShell(PETSC_COMM_WORLD,N,N,PETSC_DECIDE,PETSC_DECIDE,(  void*)(&shellData[jaccnt]), jac);CHKERRQ(ierr);
		ierr = MatShellSetOperation(*jac ,MATOP_MULT, (void(*)(void)) JacobianMatVecMatrixFree);CHKERRQ(ierr);
		ierr=MatShellSetOperation(*jac ,MATOP_GET_DIAGONAL, (void(*)(void)) JacobianMatGetDiagonal);CHKERRQ(ierr);
		jaccnt++;

		//printf("Inside CreateJacobianMatrixFree Level = %d\n ",dmmg->nlevels);

	}

	else

	{// create matrix-based jacobian at the coarsest level

#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
		ierr = DAGetMatrix(da, MATAIJ,jac);CHKERRQ(ierr);
#elif (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 2))
		ierr = DMGetMatrix(da, MATAIJ,jac);CHKERRQ(ierr);
#elif (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 4))
		ierr = DMCreateMatrix(da, MATAIJ,jac);CHKERRQ(ierr);
#else
		ierr = DMCreateMatrix(da, jac);CHKERRQ(ierr);
#endif

		//printf("Inside CreateJacobianMatrixBased Level = %d\n ",dmmg->nlevels);


	}


	//  if (dmmg->nlevels==2) {
	//  MatComputeExplicitOperator(*jac,&B);
	//  iC(PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Mf11\0",&myviewer));
	//  iC(PetscViewerSetFormat(myviewer,PETSC_VIEWER_ASCII_MATLAB));
	//  iC(MatView(B,myviewer));
	//  }
	//  else {
	//  iC(PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Mf22\0",&myviewer));
	//  iC(PetscViewerSetFormat(myviewer,PETSC_VIEWER_ASCII_MATLAB));
	//  iC(MatView(B,myviewer));
	//  }
	//  
	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "ComputeJacobianMatrixFree"
PetscErrorCode ComputeJacobianMatrixFree(stsDMMG dmmg,Mat J)
{
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	DA da = (DA)dmmg->dm;
#else
	DM da = dmmg->dm;
#endif
	PetscInt mx,my,mz;
	stsData *data;

	PetscErrorCode ierr;
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0); CHKERRQ(ierr);
#else
	ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
#endif
	nnc = (mx) * (my) ;
	if (nsd == 3) nnc = nnc * (mz);

	if(dmmg->nlevels<mgnlevels)
	{
		// compute matrix-free jacobian at all levels except for the coarsest one
		ierr = MatShellGetContext( J, (void **)&data);CHKERRQ(ierr);
		ierr = ComputeMaterialProperties(data);
		//ierr = jacobiprec(data);
		//printf("Inside ComputeJacobianMatrixFree Level = %d\n ",dmmg->nlevels);
	}
	else
	{// compute matrix-based jacobian at the coarsest level
		ierr = ComputeJacobian(dmmg,J);CHKERRQ(ierr);
		//printf("Inside CreateJacobianMatrixBased Level = %d\n ",dmmg->nlevels);
	}

	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "ComputeJacobianMatDiagonal"
int ComputeJacobianMatDiagonal(stsData *data,Vec Vecdiag)
{
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	DA da = (DA)(data->dmmg->dm);
#else
	DM da = data->dmmg->dm;
#endif
	int    ie1,ierr,i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
	int    counti, countj, countk, ie,inl,knl;
	int index2, jacobian_flag, *inc1,inc;
	int rowi, rowj, rowk, index3, index4;
	int iflagx, iflagy, iflagz;

	double shyshy, shxshx, shzshz, tmp1;
	double shdshd,sh0sh0;
	double shxshy,shxshz, shyshx,shyshz, shzshy,shzshx;

	double h;
	int Nfebs;

	PetscReal *diagarray;
	//PetscReal *lambdamf, *mumf;
	PetscReal lamnode, munode, lampmunode;                                                                                                                           
	PetscReal *vloc1;

#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0);CHKERRQ(ierr);  
	ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
#else
	ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);  
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
#endif

	inhomogenity_constant = 1;

	nnc = (xm) * (ym) ;
	if (nsd == 3) nnc = nnc * (zm);
	ne = (my-1) * (mx-1);
	if (nsd ==3) ne = ne * (mz-1);

	Hx = Lx/ (PetscReal)(mx-1); Hy = Ly/ (PetscReal)(my-1); Hz = Lz/ (PetscReal)(mz-1);

	VecGetSize(data->lambdavec,&Nfebs);
	h = 0.0;
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	ierr = VecSet(&h, Vecdiag);CHKERRQ(ierr);
#else
	ierr = VecSet(Vecdiag, h); CHKERRQ(ierr); 
#endif


	//printf("matrixfree-CJMD %d %d\n",ne,Nfebs);

	ierr = VecGetArray1d(Vecdiag,nnc*nsd,0,&diagarray);

	ierr = VecGetArray1d(data->muvec,ne,0,&mumf);
	ierr = VecGetArray1d(data->lambdavec,ne,0,&lambdamf);

	PetscMalloc(nen*nsd*sizeof(PetscReal),&vloc1);
	PetscMalloc(nen*sizeof(int),&inc1);

	//  for (inc =0;inc<nnc;inc++){
	//	knl = inc * nsd;
	//	diagarray[knl]=0.0;
	//	diagarray[knl+1]=0.0;
	//	diagarray[knl+2]=0.0;
	//  }

	jacobian_flag = 1;
	counti=1; countj=1; countk=1;
	i = xs; j = ys; k = zs;
	/* THE DIAGONAL is being computed here*/
	for (ie=0;ie<ne;ie++){
		ie1 = ie * nen;
		i = counti-1; j = countj-1; k = countk-1;

		rowi = i; rowj =j;rowk =k;
		inc1[0] = node_index3(rowi,rowj,rowk, xm, ym);
		rowi = i+1; rowj =j;rowk =k;
		inc1[1] = node_index3(rowi,rowj,rowk,xm,ym);
		rowi = i+1; rowj =j+1;rowk =k;
		inc1[2] = node_index3(rowi,rowj,rowk,xm,ym);
		rowi = i; rowj =j+1;rowk =k;
		inc1[3] = node_index3(rowi,rowj,rowk,xm,ym);

		rowi = i; rowj =j;rowk =k+1;
		inc1[4] = node_index3(rowi,rowj,rowk,xm,ym);
		rowi = i+1; rowj =j;rowk =k+1;
		inc1[5] = node_index3(rowi,rowj,rowk,xm,ym);
		rowi = i+1; rowj =j+1;rowk =k+1;
		inc1[6] = node_index3(rowi,rowj,rowk,xm,ym);
		rowi = i; rowj =j+1;rowk =k+1;
		inc1[7] = node_index3(rowi,rowj,rowk,xm,ym);

		lamnode     = lambdamf[ie];

		//printf("matfree.c ie lamnode %d %g\n", ie, lamnode);

		munode      = mumf[ie];
		lampmunode  = lamnode + munode;
		sh0sh0 = 0.0;
		for (inl=0;inl<nen;inl++){
			index2 = inl * nsd;
			index3 = index2 + 1;
			index4 = index2 + 2;

			shxshx = data->vshxshx[inl][inl];
			shxshy = data->vshxshy[inl][inl];
			shxshz = data->vshxshz[inl][inl];
			shyshy = data->vshyshy[inl][inl];
			shyshz = data->vshyshz[inl][inl];
			shyshx = data->vshyshx[inl][inl];
			shzshx = data->vshzshx[inl][inl];
			shzshy = data->vshzshy[inl][inl];
			shzshz = data->vshzshz[inl][inl];
			shdshd = data->vshdshd[inl][inl];

			tmp1   = (shdshd*munode + shxshx * lampmunode);
			vloc1[index2]   = tmp1;

			tmp1   = (shdshd*munode + shyshy * lampmunode);
			vloc1[index3]   = tmp1; 

			tmp1   = (shdshd*munode + shzshz * lampmunode);
			vloc1[index4]   = tmp1; 

		}
		for (inl=0;inl<nen;inl++){
			inc = inc1[inl];
			inc = inc * nsd;
			knl = inl * nsd;
			diagarray[inc]   +=  vloc1[knl];
			diagarray[inc+1] +=  vloc1[knl+1];
			diagarray[inc+2] +=  vloc1[knl+2];
		}
		counti++;
		if (counti==xm){ countj++; counti=1; } 
		if (countj==ym){ counti=1; countj=1; countk++; }
	}
	/* THE DIAGONAL is being computed here*/

	counti=1; countj=1; countk=1;
	for (inc =0;inc<nnc;inc++){
		i = counti-1; j = countj-1; k = countk-1;
		knl = inc * nsd;
		iflagx = 0;iflagy = 0;iflagz = 0;

		if (counti==1 && boundaryid[3].idf == 0) iflagx=1;
		if (counti==1 && boundaryid[3].jdf == 0) iflagy=1;
		if (nsd ==3 && counti==1 && boundaryid[3].kdf == 0) iflagz=1;

		if (counti==mx && boundaryid[1].idf == 0) iflagx=1;
		if (counti==mx && boundaryid[1].jdf == 0) iflagy=1;
		if (nsd ==3 && counti==mx && boundaryid[1].kdf == 0) iflagz=1;

		if (countj==1 && boundaryid[0].idf == 0) iflagx=1;
		if (countj==1 && boundaryid[0].jdf == 0) iflagy=1;
		if (nsd ==3 && countj==1 && boundaryid[0].kdf == 0) iflagz=1;

		if (countj==my && boundaryid[2].idf == 0) iflagx=1;
		if (countj==my && boundaryid[2].jdf == 0) iflagy=1;
		if (nsd == 3 && countj==my && boundaryid[2].kdf == 0) iflagz=1;

		if ( nsd ==3){
			if (countk==1 && boundaryid[4].idf == 0) iflagx=1;
			if (countk==1 && boundaryid[4].jdf == 0) iflagy=1;
			if (countk==1 && boundaryid[4].kdf == 0) iflagz=1;

			if (countk==mz && boundaryid[5].idf == 0) iflagx=1;
			if (countk==mz && boundaryid[5].jdf == 0) iflagy=1;
			if (countk==mz && boundaryid[5].kdf == 0) iflagz=1;
		}

		if (iflagx != 0){
			if (mgnlevels>1) {diagarray[knl]  = 1.0*youngs_global;}
			else {diagarray[knl]  = 1.0;}
		}

		if (iflagy != 0){
			if (mgnlevels>1) {diagarray[knl+1]  = 1.0*youngs_global;}
			else {diagarray[knl+1]  = 1.0;}
		}

		if (iflagz != 0){
			if (mgnlevels>1) {diagarray[knl+2]  = 1.0*youngs_global;}
			else {diagarray[knl+2]  = 1.0;}
		}

		counti++;
		if (counti==mx+1){ countj++; counti=1; }
		if (countj==my+1){ counti=1; countj=1; countk++; }
	}

	VecRestoreArray1d(Vecdiag,nnc*nsd,0,&diagarray);
	VecRestoreArray1d(data->muvec,ne,0,&mumf);
	VecRestoreArray1d(data->lambdavec,ne,0,&lambdamf);

	PetscFree(inc1);
	PetscFree(vloc1);

	//ierr=VecView(Vecdiag,PETSC_VIEWER_STDOUT_WORLD);

	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "JacobianMatGetDiagonal"
PetscErrorCode JacobianMatGetDiagonal(Mat J, Vec diag)
{
	PetscErrorCode ierr;
	stsData *data;
	//int diagsize;
	//PetscViewer myviewer;

	ierr = MatShellGetContext( J, (void **)&data);CHKERRQ(ierr);
	//ierr = ComputeMaterialProperties(data);
	ierr = ComputeJacobianMatDiagonal(data, diag); CHKERRQ(ierr);

	//  if (data->dmmg->nlevels==1) {
	//  iC(PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Diag1\0",&myviewer));
	//  iC(PetscViewerSetFormat(myviewer,PETSC_VIEWER_ASCII_MATLAB));
	//  iC(VecView(diag,myviewer));
	//  }
	//  else {
	//  iC(PetscViewerASCIIOpen(PETSC_COMM_WORLD,"Diag2\0",&myviewer));
	//  iC(PetscViewerSetFormat(myviewer,PETSC_VIEWER_ASCII_MATLAB));
	//  iC(VecView(diag,myviewer));
	//  }

	//ierr=VecView(diag,PETSC_VIEWER_STDOUT_WORLD);

	//ierr=VecGetSize(diag,&diagsize);
	//printf("JACMGD %d\n", diagsize);

	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "ComputeMaterialProperties"
PetscErrorCode ComputeMaterialProperties(stsData *data)
{
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	DA da = (DA) ( data->dmmg->dm);
#else
	DM da = data->dmmg->dm;
#endif
	PetscErrorCode ierr;
	int mx,my,mz,xm,ym,zm,xs,ys,zs;
	//PetscScalar *lambdamf,*mumf, *vbacmf;

#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0);CHKERRQ(ierr);  
	ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
#else
	ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);  
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
#endif

	nnc = (mx) * (my) ;
	if (nsd == 3) nnc = nnc * (mz);
	ne = (my-1) * (mx-1);
	if (nsd ==3) ne = ne * (mz-1);

	Hx = Lx/ (PetscReal)(mx-1); 
	Hy = Ly/ (PetscReal)(my-1); 
	Hz = Lz/ (PetscReal)(mz-1);

	mxg = mx; myg  = my; mzg = mz;

	PetscMalloc(ne*sizeof(PetscScalar),&lambdamf);
	PetscMalloc(ne*sizeof(PetscScalar),&mumf);
	PetscMalloc(nnc*nsd*sizeof(PetscScalar),&vbacmf);
	createien(xm, ym);
	inhomogeneous_material(lambdamf, mumf, ne,xm,ym,Hx,Hy,Hz);

	//printf("matrixfree ne %d\n", ne);


#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 2))
	VecCreateSeqWithArray(PETSC_COMM_SELF,ne,PETSC_NULL,&data->lambdavec);
	VecCreateSeqWithArray(PETSC_COMM_SELF,ne,PETSC_NULL,&data->muvec);
	VecCreateSeqWithArray(PETSC_COMM_SELF,nnc*nsd,PETSC_NULL,&data->vbac);
#else
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, ne, PETSC_NULL, &data->lambdavec); CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, ne, PETSC_NULL, &data->muvec); CHKERRQ(ierr);
	ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, nnc*nsd, PETSC_NULL, &data->vbac); CHKERRQ(ierr);
#endif
	ierr=VecPlaceArray(data->lambdavec, lambdamf); CHKERRQ(ierr);
	ierr=VecPlaceArray(data->muvec,  mumf); CHKERRQ(ierr);
	ierr=VecPlaceArray(data->vbac, vbacmf); CHKERRQ(ierr);
	denseelemmat(data, mx,my,mz);

	PetscFree(ien);

	return 0;
}

#undef __FUNCT__
#define __FUNCT__ "JacobianMatVec"
PetscErrorCode JacobianMatVecMatrixFree(Mat J, Vec in, Vec out)
{
	PetscErrorCode ierr;
	stsData *ptr;

	ierr = MatShellGetContext( J, (void **)&ptr);CHKERRQ(ierr);
	ierr = matrixFreeMatVec(ptr, in, out); CHKERRQ(ierr);
	return 0;
}


void denseelemmat(stsData *data, int mx, int my, int mz)
{
	int i,j;
	int count, inl,jnl,knl,lnl;
	int iq;

	double facx, facy, facz;
	double facHxsq, facHysq, facHzsq, facHxHy, facHxHz,facHyHz;
	double **xq,**sq,***sh,*det,*wq,**sh1;
	double shyshy, shxshx, shzshz, sh0i,sh0j;
	double shyi,shxi,shzi,shdshd, shyj,shxj,shzj;
	double shxshy,shxshz, shyshx,shyshz, shzshy,shzshx,eff0;
	int quadrature_based;
	quadrature_based = 0;

	/* QUADARATURE */
	Hx = Lx/ (PetscReal)(mx-1); 
	Hy = Ly/ (PetscReal)(my-1); 
	Hz = Lz/ (PetscReal)(mz-1);

	facHxsq = Hy*Hz/ (16*Hx);
	facHysq = Hx*Hz/ (16*Hy);
	facHzsq = Hx*Hy/ (16*Hz);

	facHxHy = Hz/ 16;
	facHxHz = Hy/ 16;
	facHyHz = Hx/ 16;

	for (inl=0;inl<nen;inl++){
		for (jnl=0;jnl<nen;jnl++){
			data->vshxshx[inl][jnl] = 0.0;
			data->vshxshy[inl][jnl] = 0.0;
			data->vshxshz[inl][jnl] = 0.0;
			data->vshyshx[inl][jnl] = 0.0;
			data->vshyshz[inl][jnl] = 0.0;
			data->vshyshy[inl][jnl] = 0.0;
			data->vshzshx[inl][jnl] = 0.0;
			data->vshzshy[inl][jnl] = 0.0;
			data->vshzshz[inl][jnl] = 0.0;
			data->vshdshd[inl][jnl] = 0.0;
		}
	}
	if (quadrature_based == 1){
		PetscMalloc(nquad*sizeof(PetscReal),&det);
		PetscMalloc(nquad*sizeof(PetscReal),&wq);
		PetscMalloc(nquad*sizeof(PetscReal),&xq);
		PetscMalloc(nquad*sizeof(PetscReal),&sh);
		PetscMalloc(nen*sizeof(PetscReal),&sq);

		for (j=0;j<nen;j++) PetscMalloc(nquad*sizeof(PetscReal),&sq[j]);
		for (j=0;j<nquad;j++){
			PetscMalloc(nsd*sizeof(PetscReal),&xq[j]);
			PetscMalloc(nen*sizeof(PetscReal),&sh[j]);
			for (i=0;i<nen;i++) PetscMalloc(nsd*sizeof(PetscReal),&sh[j][i]);
		}
		quad(xq, sq, sh,det,wq,Hx, Hy, Hz);

		for (iq=0;iq<nquad;iq++){
			count = 0;
			eff0 = wq[iq] * det[iq];
			for (inl=0;inl < nen;inl++){
				knl = inl * ndf;
				sh0i = sq[inl][iq];
				shxi = sh[iq][inl][xsd] * eff0;
				shyi = sh[iq][inl][ysd] * eff0;
				shzi = sh[iq][inl][zsd] * eff0;
				for (jnl=0;jnl<nen;jnl++){
					lnl = jnl * ndf;
					sh0j = sq[jnl][iq];
					shxj = sh[iq][jnl][xsd];
					shyj = sh[iq][jnl][ysd];

					shxshx = shxi * shxj;
					shyshy = shyi * shyj;
					shxshy = shxi * shyj;
					shyshx = shyi * shxj;
					shdshd = shxshx + shyshy;

					shzj = sh[iq][jnl][zsd];
					shxshz = shxi * shzj;
					shzshz = shzi * shzj;
					shyshz = shyi * shzj;
					shzshx = shzi * shxj;
					shzshy = shzi * shyj;
					shdshd +=shzshz;

					data->vshxshx[inl][jnl] += shxshx;
					data->vshxshy[inl][jnl] += shxshy;
					data->vshxshz[inl][jnl] += shxshz;
					data->vshyshx[inl][jnl] += shyshx;
					data->vshyshz[inl][jnl] += shyshz;
					data->vshyshy[inl][jnl] += shyshy;
					data->vshzshx[inl][jnl] += shzshx;
					data->vshzshy[inl][jnl] += shzshy;
					data->vshzshz[inl][jnl] += shzshz;
					data->vshdshd[inl][jnl] += shdshd;
				}
			}
		}
		PetscFree(det);
		PetscFree(wq);
		PetscFree(xq);
		PetscFree(sh);
		PetscFree(sq);
	}else{
		PetscMalloc(nen*sizeof(PetscReal),&sh1);
		for (j=0;j<nen;j++){
			PetscMalloc(nsd*sizeof(PetscReal),&sh1[j]);
		}
		analnode(sh1);
		for (inl=0;inl < nen;inl++){
			knl = inl * ndf;
			shxi = sh1[inl][xsd];
			shyi = sh1[inl][ysd];
			shzi = sh1[inl][zsd];

			for (jnl=0;jnl<nen;jnl++){
				lnl = jnl * ndf;
				shxj = sh1[jnl][xsd];
				shyj = sh1[jnl][ysd];

				facx = (1 + shxi * shxj / 3);
				facy = (1 + shyi * shyj / 3);

				shxshx = facHxsq * facy * shxi * shxj;
				shyshy = facHysq * facx * shyi * shyj;
				shxshy = facHxHy * shxi * shyj;
				shyshx = facHxHy * shyi * shxj;

				shzj   = sh1[jnl][zsd];
				facz   = (1 + shzi * shzj/3);

				shxshy = shxshy * facz;
				shyshx = shyshx * facz;
				shxshx = shxshx * facz;
				shyshy = shyshy * facz;

				shxshz = facHxHz * shxi * shzj * facy;
				shzshx = facHxHz * shzi * shxj * facy;
				shyshz = facHyHz * shyi * shzj * facx;
				shzshy = facHyHz * shzi * shyj * facx;
				shzshz = facHzsq * facx * facy * shzi * shzj;

				shdshd = shxshx + shyshy + shzshz;

				data->vshxshx[inl][jnl] = shxshx;
				data->vshxshy[inl][jnl] = shxshy;
				data->vshxshz[inl][jnl] = shxshz;
				data->vshyshx[inl][jnl] = shyshx;
				data->vshyshz[inl][jnl] = shyshz;
				data->vshyshy[inl][jnl] = shyshy;
				data->vshzshx[inl][jnl] = shzshx;
				data->vshzshy[inl][jnl] = shzshy;
				data->vshzshz[inl][jnl] = shzshz;
				data->vshdshd[inl][jnl] = shdshd;
			}
		}
		PetscFree(sh1);
	}
}
