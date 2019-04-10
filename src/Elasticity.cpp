///////////////////////////////////////////////////////////////////////////////////////
// Elasticity.c
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
#include "string.h"
#include "common.h"
#include "function.h"
#include "RPinclude.h"
#include "PCShellFiles/pcShellHeader.h"
#include "PCShellFiles/pcShell1.h"
#include "PCShellFiles/pcShell2.h"
#include "PCShellFiles/PCShellMain.h"

#include "PCvar.h"


#define DESTROY_DMMG_HACK

//#define SAVE_MATPROP
//#define SAVE_DEFJAC
#define SAVE_DEFFIELD
#define OPT_CASE_2	//if serial scans/manual landmarks available

extern PetscScalar *lambdamf,*mumf,*vbacmf;


// Solves (static) linear elasticity at each time step
// Input: Vec lambda, Vec mu (elastic mat. properties on the *fine* level, defined *element-wise*)
// Input: Vec force (*nodal* force)
// Output: Vec disp (*nodal* displacement x,y,z)

#undef __FUNCT__
#define __FUNCT__ "SolveElasticity"
int SolveElasticity(Vec lambda, Vec mu, Vec force, Vec disp)
{
	int ic,i,j,k,countx,county,countz,knl;
	PetscReal zero;
	PetscReal *f, *soldisp;

	/*// Initialize dmmg context *if* necessary
	if ( nsd ==3){
	ierr = stsDMMGCreate(PETSC_COMM_WORLD,mgnlevels,PETSC_NULL,&dmmg);CHKERRQ(ierr);
	ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,ndimx,ndimy,ndimz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,3,1,0,0,0,&da); */
	iCnt=0; //CH, March, 23 2007
	//RS: This must be called before calling stsDMMGSetDM
	iC(stsDMMGSetInterpolationMatrixFree(dmmg,CreateInterpolationMatrixFree,ComputeInterpolationMatrixFree));
	/*}else{
	ierr = stsDMMGCreate(PETSC_COMM_WORLD,mgnlevels,PETSC_NULL,&dmmg);CHKERRQ(ierr);
	ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,ndimx,ndimy,PETSC_DECIDE,PETSC_DECIDE,2,1,0,0,&da); 
	}
	ierr = stsDMMGSetDM(dmmg,(DM)da);
	ierr = DADestroy(da);CHKERRQ(ierr);*/


	/***************************************************************************************************/
	/***** Fine mesh statistics ************************************************************************/
	/***************************************************************************************************/
	/* myf = (ndimy-1) * pow(2,mgnlevels-1) + 1;
	mxf = (ndimx-1) * pow(2,mgnlevels-1) + 1;
	mzf = (ndimz-1) * pow(2,mgnlevels-1) + 1;
	xmf = mxf; ymf = myf;zmf = mzf;

	Hxf = Lx/ (PetscReal)(mxf-1);
	Hyf = Ly/ (PetscReal)(myf-1);
	Hzf = Lz/ (PetscReal)(mzf-1);

	ne_fine = (myf-1) * (mxf-1);
	if (nsd ==3) ne_fine = ne_fine * (mzf-1);

	nnc_fine = xmf * ymf ;
	if (nsd == 3) nnc_fine = nnc_fine * zmf;

	for (ilevel=0;ilevel<mgnlevels;ilevel++){
	mglev 		= ilevel + 1;
	levelmx[ilevel] =   (ndimx-1) * pow(2,mglev-1) + 1;
	levelmy[ilevel] =   (ndimy-1) * pow(2,mglev-1) + 1;
	levelmz[ilevel] =   (ndimz-1) * pow(2,mglev-1) + 1;
	}*/
	/***************************************************************************************************/
	/***************************************************************************************************/

	countx=mxf; //number of *nodes*
	county=myf;
	countz=mzf;

	zero=0.0;

	VecGetArray(lambda,&lambdaf); //lambdaf is *global* variable here: elast. mat. prop. lambda at the finest level, element-wise!
	VecGetArray(mu,&muf); //muf is *global* variable here: elast. mat. prop. mu at the finest level, element-wise!
	VecGetArray(force,&f); 

	ierr=PetscMalloc(nnc_fine*nsd*sizeof(PetscReal),&forcenode); //if done here, then remove from ComputeRHS function!!!

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	VecSet(&zero,disp); //Initialize output vector 
#else
	VecSet(disp, zero); // Initialize output vector 
#endif

	/* For matrixfree linear system solve*/
	if(pcShellOpt) {
		switch(pcShellOpt) {
		case 1: {
			iC(PetscStrallocpy("PCShell1:Identity (Same as None)\0",&(pcShellName)));
			break;
				}
		case 2: {
			iC(PetscStrallocpy("PCShell2:InverseDiagonal\0",&(pcShellName)));
			break;
				}
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
		default: SETERRQ(PETSC_ERR_USER,"Check the option -pcShell!");
#else
		default: SETERRQ(MPI_COMM_SELF, PETSC_ERR_USER, "Check the option -pcShell!");
#endif
		}//end switch

		//RS: This must be called before calling stsDMMGSetKSP
		iC(stsDMMGSetAllLevelsPCShell(dmmg,CreatePCShellContext,ApplyPCShell,SetUpPCShell,DestroyPCShell,pcShellName,&pcShell));
	}

	// include given nodal force (fine level) in the RHS  - here through global variable "forcenode"
	for (k=0;k<countz;k++) {
		for (j=0;j<county;j++) {
			for (i=0;i<countx;i++) { 

				ic=i+j*countx+k*countx*county; //current node 
				knl=ic*nsd;

				forcenode[knl]=f[knl]; //x-component
				forcenode[knl+1]=f[knl+1]; //y-component
				forcenode[knl+2]=f[knl+2]; //z-component

				//if (ic==111634){printf("forcex forcey forcez %g %g %g\n",forcenode[knl],forcenode[knl+1],forcenode[knl+2]);}

			} //end loop i 	  
		}//end loop j       
	}//end loop k			

	ierr = stsDMMGSetKSP(dmmg,CreateJacobianMatrixFree,ComputeRHS,ComputeJacobianMatrixFree);CHKERRQ(ierr);	
	ierr = stsDMMGSolve(dmmg);CHKERRQ(ierr);

	VecGetArray1d(dmmg[mgnlevels-1]->x,nnc_fine*nsd,0,&soldisp); //elasticity solution, matrix-free (multigrid)

	/*// DORIN ADDITION START
	iC(FreeCosminaIntermediates(dmmg));
	// DORIN ADDITION END*/

	// destroy the KSP context after each step
	iC(stsDMMGDestroyKSPAndPCShell(dmmg));

	if(lambdamf!=0){PetscFree(lambdamf);}
	if(mumf!=0){PetscFree(mumf);}
	if(vbacmf!=0){PetscFree(vbacmf);}

	//construct output Vec disp
	for (k=0;k<countz;k++) {
		for (j=0;j<county;j++) {
			for (i=0;i<countx;i++) { 
				ic=i+j*countx+k*countx*county; //current node 
				knl=ic*nsd;
				VecSetValue(disp,knl,soldisp[knl],INSERT_VALUES); //x-component
				VecSetValue(disp,knl+1,soldisp[knl+1],INSERT_VALUES); //y-component
				VecSetValue(disp,knl+2,soldisp[knl+2],INSERT_VALUES); //z-component

				//printf("ic soldispx soldispy soldislz %d %g %g %g\n",ic,soldisp[knl],soldisp[knl+1],soldisp[knl+2]);

			} //end loop i 	  
		}//end loop j       
	}//end loop k	

	VecAssemblyBegin(disp);
	VecAssemblyEnd(disp); 

	VecRestoreArray(lambda,&lambdaf);   
	VecRestoreArray(mu,&muf); 
	VecRestoreArray(force,&f);	

	return 0;
}//end function SolveElasticity
