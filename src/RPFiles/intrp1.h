///////////////////////////////////////////////////////////////////////////////////////
// intrp1.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef __INTRP1_H
#define __INTRP1_H

#undef __FUNCT__
#define __FUNCT__ "AddInterpolationMatVec"
PetscErrorCode AddInterpolation1MatVec(Mat I, Vec v1, Vec v2, Vec v3)
{
	PetscScalar one = 1.0,zero = 0.0;
	Vec tmp;
	PetscFunctionBegin;
	iC(VecDuplicate(v3,&tmp));
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	iC(VecSet(&zero,tmp));
#else
	iC(VecSet(tmp, zero));
#endif
	iC(Interpolation1MatVec(I,v1,tmp));//tmp = I*v1
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	//VecAXPY(PetscScalar *a,Vec x, Vec y);
	iC(VecAXPY(&one,v2,tmp));//tmp = v2 + tmp = v2 + I*v1
#else
	//VecAXPY(Vec y,PetscScalar a,Vec x);
	iC(VecAXPY(tmp, one, v2)); //tmp = v2 + tmp = v2 + I*v1
#endif
	iC(VecCopy(tmp,v3));//v3 = tmp
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	iC(VecDestroy(tmp));
#else
	iC(VecDestroy(&tmp));
#endif

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InterpolationMatVec"
PetscErrorCode Interpolation1MatVec(Mat I, Vec in, Vec out)
{
	IMFreeData *data;
	PetscScalar  *u;//,*ep;
	PetscInt     mx,my,mz,fmx,fmy,fmz,ci,cj,ck,CNref,count;
	PetscScalar  hx,hy,hz;
	PetscScalar  zero =0.0;
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	DA da;
#else
	DM da;
#endif
	stsDMMG dmmg;
	PetscFunctionBegin;
	iC(MatShellGetContext(I, (void **)&data));
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
	iC(VecGetArray(in,&u));//coarse grid vector
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	iC(VecSet(&zero,out)); //Initialize out
#else
	iC(VecSet(out,zero)); //Initialize out
#endif
	//Loop through All Nodes on coarse grid.
	// loc(ci,cj,ck) gives the coarse grid node number.
	//cinternal Nodes: ci=1:mx-2,cj=1:my-2,ck=1:mz-2
	for (ck=1; ck<=mz-2; ck++){
		for (cj=1; cj<=my-2; cj++){
			for(ci=1; ci<=mx-2; ci++){

				//CH: try to account for 3 degrees of freedom/node
				for (count=0;count<=2;count++){

					PetscInt FNref[27]; PetscScalar Vals[27];
					CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
					//floc(i,j,k):=0-based global node number of this node on the fine grid
					//injection for old nodes: i-->2*fi
					FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
					//The contribution is 0.5 for (i+di/2,j,k)
					FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
					FNref[2] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
					FNref[3] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
					FNref[4] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[4] = 0.5*u[CNref];
					FNref[5] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[5] = 0.5*u[CNref];
					FNref[6] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[6] = 0.5*u[CNref];
					//The contribution is 0.25 for (i+di/2,j+dj/2,k)
					FNref[7] =  nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Vals[7] = 0.25*u[CNref];
					FNref[8] =  nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Vals[8] = 0.25*u[CNref];
					FNref[9] =  nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Vals[9] = 0.25*u[CNref];
					FNref[10] = nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Vals[10] = 0.25*u[CNref];
					FNref[11] = nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Vals[11] = 0.25*u[CNref];
					FNref[12] = nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Vals[12] = 0.25*u[CNref];
					FNref[13] = nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Vals[13] = 0.25*u[CNref];
					FNref[14] = nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Vals[14] = 0.25*u[CNref];
					FNref[15] = nsd*( floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[15] = 0.25*u[CNref];
					FNref[16] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[16] = 0.25*u[CNref];
					FNref[17] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[17] = 0.25*u[CNref];
					FNref[18] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[18] = 0.25*u[CNref];
					//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
					FNref[19] =  nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[19] = 0.125*u[CNref];
					FNref[20] =  nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[20] = 0.125*u[CNref];
					FNref[21] =  nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[21] = 0.125*u[CNref];
					FNref[22] =  nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[22] = 0.125*u[CNref];
					FNref[23] =  nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[23] = 0.125*u[CNref];
					FNref[24] =  nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[24] = 0.125*u[CNref];
					FNref[25] =  nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[25] = 0.125*u[CNref];
					FNref[26] =  nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[26] = 0.125*u[CNref];
					iC(VecSetValues(out,27,FNref,Vals,ADD_VALUES));

				}//end count (degreess of freedom/node)

			}//end for ci
		}//end for cj
	}//end for ck

	//Surface Nodes: z=0
	ck=0;
	for (cj=1; cj<=my-2; cj++){
		for(ci=1; ci<=mx-2; ci++){
			PetscInt FNref[18]; PetscScalar Vals[18];

			//CH: try to account for 3 degrees of freedom/node
			for (count=0;count<=2;count++){

				CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
				//floc(i,j,k):=0-based global node number of this node on the fine grid
				//injection for old nodes: i-->2*fi
				FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
				//The contribution is 0.5 for (i+di/2,j,k)
				FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
				FNref[2] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
				FNref[3] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
				FNref[4] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[4] = 0.5*u[CNref];
				FNref[5] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[5] = 0.5*u[CNref];
				//The contribution is 0.25 for (i+di/2,j+dj/2,k)
				FNref[6] =  nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Vals[6] = 0.25*u[CNref];
				FNref[7] =  nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Vals[7] = 0.25*u[CNref];
				FNref[8] =  nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Vals[8] = 0.25*u[CNref];
				FNref[9] =  nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Vals[9] = 0.25*u[CNref];
				FNref[10] = nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Vals[10] = 0.25*u[CNref];
				FNref[11] = nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Vals[11] = 0.25*u[CNref];
				FNref[12] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[12] = 0.25*u[CNref];
				FNref[13] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[13] = 0.25*u[CNref];
				//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
				FNref[14] =  nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[14] = 0.125*u[CNref];
				FNref[15] =  nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[15] = 0.125*u[CNref];
				FNref[16] =  nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[16] = 0.125*u[CNref];
				FNref[17] =  nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[17] = 0.125*u[CNref];
				iC(VecSetValues(out,18,FNref,Vals,ADD_VALUES));

			}//end count (degreess of freedom/node)

		}//end for ci
	}//end for cj

	//Surface Nodes: z=mz-1
	ck=mz-1;
	for (cj=1; cj<=my-2; cj++){
		for(ci=1; ci<=mx-2; ci++){
			PetscInt FNref[18]; PetscScalar Vals[18];

			//CH: try to account for 3 degrees of freedom/node
			for (count=0;count<=2;count++){

				CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
				//floc(i,j,k):=0-based global node number of this node on the fine grid
				//injection for old nodes: i-->2*fi
				FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
				//The contribution is 0.5 for (i+di/2,j,k)
				FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
				FNref[2] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
				FNref[3] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
				FNref[4] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[4] = 0.5*u[CNref];
				FNref[5] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[5] = 0.5*u[CNref];
				//The contribution is 0.25 for (i+di/2,j+dj/2,k)
				FNref[6] =  nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Vals[6] = 0.25*u[CNref];
				FNref[7] =  nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Vals[7] = 0.25*u[CNref];
				FNref[8] =  nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Vals[8] = 0.25*u[CNref];
				FNref[9] =  nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Vals[9] = 0.25*u[CNref];
				FNref[10] = nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Vals[10] = 0.25*u[CNref];
				FNref[11] = nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Vals[11] = 0.25*u[CNref];
				FNref[12] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[12] = 0.25*u[CNref];
				FNref[13] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[13] = 0.25*u[CNref];
				//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
				FNref[14] = nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[14] = 0.125*u[CNref];
				FNref[15] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[15] = 0.125*u[CNref];
				FNref[16] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[16] = 0.125*u[CNref];
				FNref[17] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[17] = 0.125*u[CNref];
				iC(VecSetValues(out,18,FNref,Vals,ADD_VALUES));

			}//end count (degreess of freedom/node)

		}//end for ci
	}//end for cj

	//Surface Nodes: y=0
	cj=0;
	for (ck=1; ck<=my-2; ck++){
		for(ci=1; ci<=mx-2; ci++){
			PetscInt FNref[18]; PetscScalar Vals[18];

			//CH: try to account for 3 degrees of freedom/node
			for (count=0;count<=2;count++){

				CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
				//floc(i,j,k):=0-based global node number of this node on the fine grid
				//injection for old nodes: i-->2*fi
				FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
				//The contribution is 0.5 for (i+di/2,j,k)
				FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
				FNref[2] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
				FNref[3] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
				FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[4] = 0.5*u[CNref];
				FNref[5] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[5] = 0.5*u[CNref];
				//The contribution is 0.25 for (i+di/2,j+dj/2,k)
				FNref[6] =  nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Vals[6] = 0.25*u[CNref];
				FNref[7] =  nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Vals[7] = 0.25*u[CNref];
				FNref[8] =  nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Vals[8] = 0.25*u[CNref];
				FNref[9] =  nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Vals[9] = 0.25*u[CNref];
				FNref[10] = nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Vals[10] = 0.25*u[CNref];
				FNref[11] = nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Vals[11] = 0.25*u[CNref];
				FNref[12] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[12] = 0.25*u[CNref];
				FNref[13] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[13] = 0.25*u[CNref];
				//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
				FNref[14] = nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[14] = 0.125*u[CNref];
				FNref[15] = nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[15] = 0.125*u[CNref];
				FNref[16] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[16] = 0.125*u[CNref];
				FNref[17] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[17] = 0.125*u[CNref];
				iC(VecSetValues(out,18,FNref,Vals,ADD_VALUES));

			}//end count (degreess of freedom/node)

		}//end for ci
	}//end for ck

	//Surface Nodes: y=my-1
	cj=my-1;
	for (ck=1; ck<=my-2; ck++){
		for(ci=1; ci<=mx-2; ci++){
			PetscInt FNref[18]; PetscScalar Vals[18];

			//CH: try to account for 3 degrees of freedom/node
			for (count=0;count<=2;count++){

				CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
				//floc(i,j,k):=0-based global node number of this node on the fine grid
				//injection for old nodes: i-->2*fi
				FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
				//The contribution is 0.5 for (i+di/2,j,k)
				FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
				FNref[2] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
				FNref[3] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
				FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[4] = 0.5*u[CNref];
				FNref[5] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[5] = 0.5*u[CNref];
				//The contribution is 0.25 for (i+di/2,j+dj/2,k)
				FNref[6] = nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Vals[6] = 0.25*u[CNref];
				FNref[7] = nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Vals[7] = 0.25*u[CNref];
				FNref[8] = nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Vals[8] = 0.25*u[CNref];
				FNref[9] = nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Vals[9] = 0.25*u[CNref];
				FNref[10] = nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Vals[10] = 0.25*u[CNref];
				FNref[11] = nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Vals[11] = 0.25*u[CNref];
				FNref[12] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[12] = 0.25*u[CNref];
				FNref[13] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[13] = 0.25*u[CNref];
				//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
				FNref[14] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[14] = 0.125*u[CNref];
				FNref[15] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[15] = 0.125*u[CNref];
				FNref[16] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[16] = 0.125*u[CNref];
				FNref[17] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[17] = 0.125*u[CNref];
				iC(VecSetValues(out,18,FNref,Vals,ADD_VALUES));

			}//end count (degreess of freedom/node)

		}//end for ci
	}//end for ck


	//Surface Nodes: x=0
	ci=0;
	for (ck=1; ck<=my-2; ck++){
		for(cj=1; cj<=mx-2; cj++){
			PetscInt FNref[18]; PetscScalar Vals[18];

			//CH: try to account for 3 degrees of freedom/node
			for (count=0;count<=2;count++){

				CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
				//floc(i,j,k):=0-based global node number of this node on the fine grid
				//injection for old nodes: i-->2*fi
				FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
				//The contribution is 0.5 for (i+di/2,j,k)
				FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
				FNref[2] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
				FNref[3] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
				FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[4] = 0.5*u[CNref];
				FNref[5] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[5] = 0.5*u[CNref];
				//The contribution is 0.25 for (i+di/2,j+dj/2,k)
				FNref[6] =  nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Vals[6] = 0.25*u[CNref];
				FNref[7] =  nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Vals[7] = 0.25*u[CNref];
				FNref[8] =  nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Vals[8] = 0.25*u[CNref];
				FNref[9] =  nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Vals[9] = 0.25*u[CNref];
				FNref[10] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[10] = 0.25*u[CNref];
				FNref[11] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[11] = 0.25*u[CNref];
				FNref[12] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[12] = 0.25*u[CNref];
				FNref[13] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[13] = 0.25*u[CNref];
				//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
				FNref[14] = nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[14] = 0.125*u[CNref];
				FNref[15] = nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[15] = 0.125*u[CNref];
				FNref[16] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[16] = 0.125*u[CNref];
				FNref[17] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[17] = 0.125*u[CNref];
				iC(VecSetValues(out,18,FNref,Vals,ADD_VALUES));

			}//end count (degreess of freedom/node)


		}//end for cj
	}//end for ck

	//Surface Nodes: x=mx-1
	ci=mx-1;
	for (ck=1; ck<=my-2; ck++){
		for(cj=1; cj<=mx-2; cj++){

			PetscInt FNref[18]; PetscScalar Vals[18];

			//CH: try to account for 3 degrees of freedom/node
			for (count=0;count<=2;count++){

				CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
				//floc(i,j,k):=0-based global node number of this node on the fine grid
				//injection for old nodes: i-->2*fi
				FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
				//The contribution is 0.5 for (i+di/2,j,k)
				FNref[1] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
				FNref[2] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
				FNref[3] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
				FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[4] = 0.5*u[CNref];
				FNref[5] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[5] = 0.5*u[CNref];
				//The contribution is 0.25 for (i+di/2,j+dj/2,k)
				FNref[6] = nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Vals[6] = 0.25*u[CNref];
				FNref[7] = nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Vals[7] = 0.25*u[CNref];
				FNref[8] = nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Vals[8] = 0.25*u[CNref];
				FNref[9] = nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Vals[9] = 0.25*u[CNref];
				FNref[10] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[10] = 0.25*u[CNref];
				FNref[11] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[11] = 0.25*u[CNref];
				FNref[12] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[12] = 0.25*u[CNref];
				FNref[13] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[13] = 0.25*u[CNref];
				//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
				FNref[14] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[14] = 0.125*u[CNref];
				FNref[15] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[15] = 0.125*u[CNref];
				FNref[16] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[16] = 0.125*u[CNref];
				FNref[17] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[17] = 0.125*u[CNref];
				iC(VecSetValues(out,18,FNref,Vals,ADD_VALUES));

			}//end count (degreess of freedom/node)

		}//end for cj
	}//end for ck


	//Edge Nodes: x=0, y=0
	ci=0;cj=0;
	for (ck=1; ck<=my-2; ck++){
		PetscInt FNref[12]; PetscScalar Vals[12];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[3] = 0.5*u[CNref];
			FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[4] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[5] =  nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] =  nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Vals[6] = 0.25*u[CNref];
			FNref[7] =  nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Vals[7] = 0.25*u[CNref];
			FNref[8] =  nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[8] = 0.25*u[CNref];
			FNref[9] =  nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[9] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[10] =  nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[10] = 0.125*u[CNref];
			FNref[11] =  nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[11] = 0.125*u[CNref];
			iC(VecSetValues(out,12,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}//end for ck

	//Edge Nodes: x=0, y=my-1
	//n= 4,8
	ci=0;cj=my-1;
	for (ck=1; ck<=my-2; ck++){
		PetscInt FNref[12]; PetscScalar Vals[12];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[3] = 0.5*u[CNref];
			FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[4] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[5] =  nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] =  nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Vals[6] = 0.25*u[CNref];
			FNref[7] =  nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Vals[7] = 0.25*u[CNref];
			FNref[8] =  nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[8] = 0.25*u[CNref];
			FNref[9] =  nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[9] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[10] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[10] = 0.125*u[CNref];
			FNref[11] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[11] = 0.125*u[CNref];
			iC(VecSetValues(out,12,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}//end for ck

	//Edge Nodes: x=0, z=0
	//n= 1,4
	ci=0;ck=0;
	for (cj=1; cj<=my-2; cj++){
		PetscInt FNref[12]; PetscScalar Vals[12];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
			FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[4] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[5] = nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Vals[6] = 0.25*u[CNref];
			FNref[7] = nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Vals[7] = 0.25*u[CNref];
			FNref[8] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[8] = 0.25*u[CNref];
			FNref[9] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[9] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[10] = nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[10] = 0.125*u[CNref];
			FNref[11] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[11] = 0.125*u[CNref];
			iC(VecSetValues(out,12,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}//end for cj

	//Edge Nodes: x=0, z=mz-1
	//n= 5,8
	ci=0;ck=mz-1;
	for (cj=1; cj<=my-2; cj++){
		PetscInt FNref[12]; PetscScalar Vals[12];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
			FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[4] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[5] = nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Vals[6] = 0.25*u[CNref];
			FNref[7] = nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Vals[7] = 0.25*u[CNref];
			FNref[8] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[8] = 0.25*u[CNref];
			FNref[9] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[9] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[10] = nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[10] = 0.125*u[CNref];
			FNref[11] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[11] = 0.125*u[CNref];
			iC(VecSetValues(out,12,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}//end for cj

	//Edge Nodes: x=mx-1,y=0
	//n= 2,6
	ci=mx-1;cj=0;
	for (ck=1; ck<=mz-2; ck++){
		PetscInt FNref[12]; PetscScalar Vals[12];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[3] = 0.5*u[CNref];
			FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[4] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[5] = nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Vals[6] = 0.25*u[CNref];
			FNref[7] = nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Vals[7] = 0.25*u[CNref];
			FNref[8] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[8] = 0.25*u[CNref];
			FNref[9] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[9] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[10] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[10] = 0.125*u[CNref];
			FNref[11] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[11] = 0.125*u[CNref];
			iC(VecSetValues(out,12,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}//end for ck

	//Edge Nodes: x=mx-1,y=my-1
	//n= 3,7
	ci=mx-1;cj=my-1;
	for (ck=1; ck<=mz-2; ck++){
		PetscInt FNref[12]; PetscScalar Vals[12];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[3] = 0.5*u[CNref];
			FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[4] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[5] = nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Vals[6] = 0.25*u[CNref];
			FNref[7] = nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Vals[7] = 0.25*u[CNref];
			FNref[8] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[8] = 0.25*u[CNref];
			FNref[9] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[9] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[10] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[10] = 0.125*u[CNref];
			FNref[11] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[11] = 0.125*u[CNref];
			iC(VecSetValues(out,12,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}//end for ck

	//Edge Nodes: x=mx-1,z=0
	//n= 2,3
	ci=mx-1;ck=0;
	for (cj=1; cj<=my-2; cj++){
		PetscInt FNref[12]; PetscScalar Vals[12];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
			FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[4] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[5] = nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Vals[6] = 0.25*u[CNref];
			FNref[7] = nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Vals[7] = 0.25*u[CNref];
			FNref[8] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[8] = 0.25*u[CNref];
			FNref[9] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[9] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[10] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[10] = 0.125*u[CNref];
			FNref[11] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[11] = 0.125*u[CNref];
			iC(VecSetValues(out,12,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}//end for cj

	//Edge Nodes: x=mx-1,z=mz-1
	//n= 6,7
	ci=mx-1;ck=mz-1;
	for (cj=1; cj<=my-2; cj++){
		PetscInt FNref[12]; PetscScalar Vals[12];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
			FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[4] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[5] = nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Vals[6] = 0.25*u[CNref];
			FNref[7] = nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Vals[7] = 0.25*u[CNref];
			FNref[8] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[8] = 0.25*u[CNref];
			FNref[9] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[9] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[10] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[10] = 0.125*u[CNref];
			FNref[11] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[11] = 0.125*u[CNref];
			iC(VecSetValues(out,12,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}//end for cj

	//Edge Nodes: y=0, z= 0
	//n=1,2
	cj=0;ck=0;
	for (ci=1; ci<=mx-2; ci++){
		PetscInt FNref[12]; PetscScalar Vals[12];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
			FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[4] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[5] = nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Vals[6] = 0.25*u[CNref];
			FNref[7] = nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Vals[7] = 0.25*u[CNref];
			FNref[8] = nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Vals[8] = 0.25*u[CNref];
			FNref[9] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[9] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[10] = nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[10] = 0.125*u[CNref];
			FNref[11] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[11] = 0.125*u[CNref];
			iC(VecSetValues(out,12,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}//end for ci

	//Edge Nodes: y=0, z= mz-1
	//n=5,6
	cj=0;ck=mz-1;
	for (ci=1; ci<=mx-2; ci++){
		PetscInt FNref[12]; PetscScalar Vals[12];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
			FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[4] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[5] = nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Vals[6] = 0.25*u[CNref];
			FNref[7] = nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Vals[7] = 0.25*u[CNref];
			FNref[8] = nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Vals[8] = 0.25*u[CNref];
			FNref[9] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[9] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[10] = nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[10] = 0.125*u[CNref];
			FNref[11] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[11] = 0.125*u[CNref];
			iC(VecSetValues(out,12,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}//end for ci

	//Edge Nodes: y=my-1, z=0
	cj=my-1;ck=0;
	for (ci=1; ci<=mx-2; ci++){
		PetscInt FNref[12]; PetscScalar Vals[12];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
			FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[4] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[5] = nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Vals[6] = 0.25*u[CNref];
			FNref[7] = nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Vals[7] = 0.25*u[CNref];
			FNref[8] = nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Vals[8] = 0.25*u[CNref];
			FNref[9] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[9] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[10] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[10] = 0.125*u[CNref];
			FNref[11] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[11] = 0.125*u[CNref];
			iC(VecSetValues(out,12,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}//end for ci

	//Edge Nodes: y=my-1, z=mz-1
	cj=my-1;ck=mz-1;
	for (ci=1; ci<=mx-2; ci++){
		PetscInt FNref[12]; PetscScalar Vals[12];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[3] = 0.5*u[CNref];
			FNref[4] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[4] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[5] = nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Vals[6] = 0.25*u[CNref];
			FNref[7] = nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Vals[7] = 0.25*u[CNref];
			FNref[8] = nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Vals[8] = 0.25*u[CNref];
			FNref[9] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[9] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[10] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[10] = 0.125*u[CNref];
			FNref[11] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[11] = 0.125*u[CNref];
			iC(VecSetValues(out,12,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}//end for ci

	//Corner x=0,y=0,z=0
	//n=1
	ci=0;cj=0;ck=0;
	{
		PetscInt FNref[8]; PetscScalar Vals[8];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[3] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[4] = nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Vals[4] = 0.25*u[CNref];
			FNref[5] = nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[6] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[7] = nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[7] = 0.125*u[CNref];
			iC(VecSetValues(out,8,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}

	//Corner x=mx-1,y=0,z=0
	//n=2
	ci=mx-1;cj=0;ck=0;
	{
		PetscInt FNref[8]; PetscScalar Vals[8];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[3] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[4] = nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Vals[4] = 0.25*u[CNref];
			FNref[5] = nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[6] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[7] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1)+count; Vals[7] = 0.125*u[CNref];
			iC(VecSetValues(out,8,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}

	//Corner x=mx-1,y=my-1,z=0
	//n=3
	ci=mx-1;cj=my-1;ck=0;
	{
		PetscInt FNref[8]; PetscScalar Vals[8];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[3] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[4] = nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Vals[4] = 0.25*u[CNref];
			FNref[5] = nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[6] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[7] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[7] = 0.125*u[CNref];
			iC(VecSetValues(out,8,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}

	//Corner x=0,y=my-1,z=0
	//n=4
	ci=0;cj=my-1;ck=0;
	{
		PetscInt FNref[8]; PetscScalar Vals[8];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Vals[3] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[4] = nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Vals[4] = 0.25*u[CNref];
			FNref[5] = nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[6] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[7] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1)+count; Vals[7] = 0.125*u[CNref];
			iC(VecSetValues(out,8,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}

	//Corner x=0,y=0,z=mz-1
	//n=5
	ci=0;cj=0;ck=mz-1;
	{
		PetscInt FNref[8]; PetscScalar Vals[8];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[3] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[4] = nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Vals[4] = 0.25*u[CNref];
			FNref[5] = nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[6] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[7] = nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[7] = 0.125*u[CNref];
			iC(VecSetValues(out,8,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}

	//Corner x=mx-1,y=0,z=mz-1
	//n=6
	ci=mx-1;cj=0;ck=mz-1;
	{
		PetscInt FNref[8]; PetscScalar Vals[8];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[3] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[4] = nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Vals[4] = 0.25*u[CNref];
			FNref[5] = nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[6] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[7] = nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1)+count; Vals[7] = 0.125*u[CNref];
			iC(VecSetValues(out,8,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}

	//Corner x=mx-1,y=my-1,z=mz-1
	//n=7
	ci=mx-1;cj=my-1;ck=mz-1;
	{
		PetscInt FNref[8]; PetscScalar Vals[8];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[3] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[4] = nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Vals[4] = 0.25*u[CNref];
			FNref[5] = nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[6] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[7] = nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[7] = 0.125*u[CNref];
			iC(VecSetValues(out,8,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}

	//Corner x=0,y=my-1,z=mz-1
	//n=8
	ci=0;cj=my-1;ck=mz-1;
	{
		PetscInt FNref[8]; PetscScalar Vals[8];

		//CH: try to account for 3 degrees of freedom/node
		for (count=0;count<=2;count++){

			CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
			//floc(i,j,k):=0-based global node number of this node on the fine grid
			//injection for old nodes: i-->2*fi
			FNref[0] = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Vals[0] = u[CNref];
			//The contribution is 0.5 for (i+di/2,j,k)
			FNref[1] = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Vals[1] = 0.5*u[CNref];
			FNref[2] = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Vals[2] = 0.5*u[CNref];
			FNref[3] = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Vals[3] = 0.5*u[CNref];
			//The contribution is 0.25 for (i+di/2,j+dj/2,k)
			FNref[4] = nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Vals[4] = 0.25*u[CNref];
			FNref[5] = nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Vals[5] = 0.25*u[CNref];
			FNref[6] = nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[6] = 0.25*u[CNref];
			//The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
			FNref[7] = nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1)+count; Vals[7] = 0.125*u[CNref];
			iC(VecSetValues(out,8,FNref,Vals,ADD_VALUES));

		}//end count (degreess of freedom/node)

	}
	//iC(VecRestoreArray(data->epsilon,&ep));
	iC(VecRestoreArray(in,&u));//coarse grid vector
	iC(VecAssemblyBegin(out));//fine grid vector
	iC(VecAssemblyEnd(out));
	PetscFunctionReturn(0);
}

#endif
