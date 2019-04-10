///////////////////////////////////////////////////////////////////////////////////////
// InterpMF.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef __INTERPMF_H
#define __INTERPMF_H


/* ********************************************************************** */
/* Matrix free routines*/
/* ********************************************************************** */
#undef __FUNCT__
#define __FUNCT__ "CreateInterpolationMatrixFree"
PetscErrorCode CreateInterpolationMatrixFree(stsDMMG dmmg, Mat *I) {
	PetscInt       resType =2;//interpolation option 
	PetscFunctionBegin;
	iC(PetscOptionsGetInt(0,"-restype",&resType,0));    

	switch(resType) {
	case 1: iC(CreateInterpolation1(dmmg,I));break;
	case 2: iC(CreateInterpolation2(dmmg,I));break;
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	default: SETERRQ(PETSC_ERR_USER,"Check the option -restype!");	
#else
	default: SETERRQ(MPI_COMM_SELF, PETSC_ERR_USER, "Check the option -restype!");
#endif
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeInterpolationMatrixFree"
PetscErrorCode ComputeInterpolationMatrixFree(stsDMMG dmmg,Mat I) {
	PetscInt       resType =2;//interpolation option 
	PetscFunctionBegin;
	iC(PetscOptionsGetInt(0,"-restype",&resType,0));

	switch(resType) {
	case 1: {iC(ComputeInterpolation1(dmmg,I)); break;}
	case 2: {iC(ComputeInterpolation2(dmmg,I)); break;}
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	default: SETERRQ(PETSC_ERR_USER,"Check the option -restype!");
#else
	default: SETERRQ(MPI_COMM_SELF, PETSC_ERR_USER, "Check the option -restype!");
#endif
	}
	PetscFunctionReturn(0);
}

#endif
