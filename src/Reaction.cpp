///////////////////////////////////////////////////////////////////////////////////////
// Reaction.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include "function.h"


// Reaction term in the tumor equation - forward problem.
// Can be solved analytically over time dt.
// Input: c0 ('initial' condition at time t), paramtum (here only rho=tumor growth rate, here assumed constant), time step dt
// Output: c (updated quantity at time t+dt)
// This is always performed *node-wise* here

#undef __FUNCT__
#define __FUNCT__ "ReactionForward"
void ReactionForward(Vec c0, double rho, double dt, Vec c)
{
	//stencil (box) created in Initialize
	int countx, county, countz;
	int  ic,i,j,k;
	//int size;
	Vec c0dup;
	PetscScalar *cold, *phi;
	double zero, one, cnew, fac;

	zero=0.0;
	one=1.0;
	//fac=exp(rho*dt)-one;

	// (fine) grid information
	countx=mxf; //number of *nodes*
	county=myf;
	countz=mzf;

	//VecGetSize(c0,&size); //get size of the block-vector c0
	//printf("size %d\n",size);

	//reaction -  node-wise

	VecDuplicate(c0,&c0dup);
	VecCopy(c0,c0dup);

	VecGetArray(c0dup,&cold);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	VecSet(&zero,c); //Initialize output vector c
#else
	VecSet(c, zero); //Initialize output vector c
#endif

	VecGetArray(gphi,&phi);

	//loop over the original grid *node-wise*
	for (k=0;k<countz;k++) {
		for (j=0;j<county;j++) {
			for (i=0;i<countx;i++) { 
				ic=i+j*countx+k*countx*county; //current node   	

				//update c (closed-form solution for the quadratic term reaction) 

				//iff global phase-field variable gphi for tracking ventricles!!
				fac=exp(rho*phi[ic]*dt)-one;

				cnew=one-(one-cold[ic])/(cold[ic]*fac+one);

				VecSetValue(c,ic,cnew,INSERT_VALUES);
			} //end loop i 	  
		}//end loop j       
	}//end loop k	

	VecAssemblyBegin(c);
	VecAssemblyEnd(c); 

	VecRestoreArray(c0dup,&cold);
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	VecDestroy(c0dup);
#else
	VecDestroy(&c0dup);
#endif

	VecRestoreArray(gphi,&phi);
} //end function ReactionForward
