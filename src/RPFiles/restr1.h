///////////////////////////////////////////////////////////////////////////////////////
// restr1.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef __RESTR1_H
#define __RESTR1_H

//Restriction is simple injection, interpolation is by bilinear

#undef __FUNCT__
#define __FUNCT__ "AddRestriction1MatVec"
PetscErrorCode AddRestriction1MatVec(Mat R, Vec v1, Vec v2, Vec v3)
{
	PetscScalar one = 1.0,zero = 0.0;
	Vec tmp;
	PetscFunctionBegin;
	iC(VecDuplicate(v3,&tmp));
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	iC(VecSet(&zero,tmp));
#else
	iC(VecSet(tmp,zero));
#endif
	iC(Restriction1MatVec(R, v1, tmp));
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	//VecAXPY(PetscScalar *a,Vec x, Vec y);
	iC(VecAXPY(&one,v2,tmp));//tmp = v2 + tmp = v2 + R*v1
#else
	//VecAXPY(Vec y,PetscScalar a,Vec x);
	iC(VecAXPY(tmp, one, v2));//tmp = v2 + tmp = v2 + R*v1
#endif
	iC(VecCopy(tmp,v3));//v3 =  tmp
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	iC(VecDestroy(tmp));
#else
	iC(VecDestroy(&tmp));
#endif
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Restriction1MatVec"
PetscErrorCode Restriction1MatVec(Mat R, Vec in, Vec out)
{
	IMFreeData *data;
	PetscScalar  *u;//,*ep;
	PetscInt     mx,my,mz,fmx,fmy,fmz,ci,cj,ck,CNref,FNref,count;
	PetscScalar  hx,hy,hz;
	PetscScalar  zero =0.0;
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	DA da;
#else
	DM da;
#endif
	stsDMMG dmmg;
	PetscFunctionBegin;
	iC(MatShellGetContext(R, (void **)&data));
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	dmmg = data->dmmg; da = (DA)dmmg->dm; 
	iC(DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0));//mx,my,mz = number of nodes in each direction
#else
	dmmg = data->dmmg; da = dmmg->dm; 
	iC(DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0));//mx,my,mz = number of nodes in each direction
#endif
	fmx = (2*(mx-1)) + 1;fmy = (2*(my-1)) + 1;fmz = (2*(mz-1)) + 1;
	hx=Lx/(((PetscScalar)mx)-1.0);hy=Ly/(((PetscScalar)my)-1.0);
	hz=Lz/(((PetscScalar)mz)-1.0);
	//iC(VecGetArray(data->epsilon,&ep));
	iC(VecGetArray(in,&u));//fine grid vector
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	iC(VecSet(&zero,out)); //Initialize out
#else
	iC(VecSet(out,zero)); //Initialize out
#endif
	//Loop through All Nodes on coarse grid.
	// loc(ci,cj,ck) gives the coarse grid node number.
	for (ck=0; ck<=mz-1; ck++){
		for (cj=0; cj<=my-1; cj++){
			for(ci=0; ci<=mx-1; ci++){

				//CH: try to account for 3 degrees of freedom/node
				for (count=0;count<=2;count++){

					CNref = nsd*(loc(ci,cj,ck)-1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
					//injection from shared nodes: i-->2*fi
					FNref =  nsd*(floc(2*ci,2*cj,2*ck)-1)+count;//floc(i,j,k):=0-based global node number of this node on the fine grid
					iC(VecSetValue(out,CNref,u[FNref],INSERT_VALUES));

				}//end count (degreess of freedom/node)

			}//end ci
		}//end cj
	}//end ck
	//iC(VecRestoreArray(data->epsilon,&ep));
	iC(VecRestoreArray(in,&u));//fine grid vector
	iC(VecAssemblyBegin(out));//coarse grid vector
	iC(VecAssemblyEnd(out));
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CreateInterpolation1"
PetscErrorCode CreateInterpolation1(stsDMMG dmmg, Mat *interp)
{
	//RS: Note the input is dmmg instead of da, because assignMAterial requires dmmg.
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	DA da = (DA)(dmmg->dm);
#else
	DM da = dmmg->dm;
#endif
	PetscInt  M,N,m,n,MD,ND,PD,md,nd,pd;    
	PetscFunctionBegin;
	//printf("Inside CreateInterpolation 1. Level = %d , cnt = %d\n",dmmg->nlevels,iCnt);
	//RS: This is coarse grid info
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	iC(DAGetInfo(da,PETSC_NULL,&MD,&ND,&PD,&md,&nd,&pd,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL));
#else
	iC(DMDAGetInfo(da,PETSC_NULL,&MD,&ND,&PD,&md,&nd,&pd,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL));
#endif
	n= N= nsd*MD*ND*PD;//cols = coarse
	m= M = nsd*(2*(MD-1)+1)*(2*(ND-1)+1)*(2*(PD-1)+1);//rows = fine
	assert(iCnt <100);
	idata[iCnt].dmmg = dmmg;
	iC(MatCreateShell(PETSC_COMM_WORLD, m ,n,M,N,(void*)(&idata[iCnt]),interp));
	iC(MatShellSetOperation(*interp,MATOP_MULT, (void(*)(void))Interpolation1MatVec));
	iC(MatShellSetOperation(*interp,MATOP_MULT_ADD, (void(*)(void))AddInterpolation1MatVec));
	iC(MatShellSetOperation(*interp,MATOP_MULT_TRANSPOSE,(void(*)(void))Restriction1MatVec));//Injection
	iC(MatShellSetOperation(*interp,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))AddRestriction1MatVec));//Injection
	iCnt ++;
	//printf("Finished CreateInterpolation 1. Level = %d, cnt = %d\n",dmmg->nlevels,iCnt);
	//printf("M=%d, N= %d\n",M,N);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeInterpolation1"
PetscErrorCode ComputeInterpolation1(stsDMMG dmmg,Mat I)
{
	PetscFunctionBegin;
	//printf("Inside ComputeInterpolation. Level = %d\n",dmmg->nlevels);
	/*levels*/
	//printf("Finished ComputeInterpolation. Level = %d\n",dmmg->nlevels);
	PetscFunctionReturn(0);
}
#endif

