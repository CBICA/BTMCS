///////////////////////////////////////////////////////////////////////////////////////
// pcShell1.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef __PC_SHELL_1_H
#define __PC_SHELL_1_H

//typedef struct{} PC1ShellCtx;

#undef __FUNCT__
#define __FUNCT__ "CreatePC1Context"
PetscErrorCode CreatePC1Context(PC1ShellCtx** pc1Sh) {
	PC1ShellCtx *tmpPc1;
	PetscFunctionBegin;
	//printf("Inside CreatePC1Context.\n");
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 4))
	iC(PetscNew(PC1ShellCtx,&tmpPc1));
#else
	iC(PetscNew(&tmpPc1));
#endif
	//printf("Finished CreatePC1Context.\n");
	*pc1Sh = tmpPc1;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetUpPC1"
PetscErrorCode SetUpPC1(PC1ShellCtx* pc1Sh,Mat B) {
	PetscFunctionBegin;
	//printf("Inside SetUpPC1.\n");
	//printf("Finished SetUpPC1.\n");
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ApplyPC1"
PetscErrorCode ApplyPC1(PC1ShellCtx* pc1Sh,Vec in,Vec out) {
	PetscFunctionBegin;
	//printf("Inside ApplyPC1.\n");
	iC(VecCopy(in,out));
	//printf("Finished ApplyPC1.\n");
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DestroyPC1"
PetscErrorCode DestroyPC1(PC1ShellCtx* pc1Sh) {
	PetscFunctionBegin;
	//printf("Inside DestroyPC1.\n");
	PetscFree(pc1Sh);
	//printf("Finished DestroyPC1.\n");
	PetscFunctionReturn(0);
}

#endif
