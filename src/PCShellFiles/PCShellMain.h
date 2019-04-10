///////////////////////////////////////////////////////////////////////////////////////
// PCShellMain.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef __PC_SHELL_MAIN_H
#define __PC_SHELL_MAIN_H

// {void *pcShellStruct; } PCShellCtx;

#undef __FUNCT__
#define __FUNCT__ "CreatePCShellContext"
PetscErrorCode CreatePCShellContext(void** pcSh) {
	PetscInt pcShellOpt = 1;//pcShell option
	PCShellCtx *tmpPcSh;
	PetscFunctionBegin;
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 4))
	iC(PetscNew(PCShellCtx,&tmpPcSh));
#else
	iC(PetscNew(&tmpPcSh));
#endif
	iC(PetscOptionsGetInt(0,"-pcShell",&pcShellOpt,0));
	switch(pcShellOpt) {
	case 1: {
		iC(CreatePC1Context((PC1ShellCtx**)(&(tmpPcSh->pcShellStruct))));
		break;
			}
	case 2: {
		iC(CreatePC2Context((PC2ShellCtx**)(&(tmpPcSh->pcShellStruct))));
		break;
			}
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	default: SETERRQ(PETSC_ERR_USER,"Check the option -pcShell!");
#else
	default: SETERRQ(MPI_COMM_SELF, PETSC_ERR_USER, "Check the option -pcShell!");
#endif
	}
	*pcSh = tmpPcSh;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetUpPCShell"
PetscErrorCode SetUpPCShell(void* pcSh,Mat B) {
	PetscInt pcShellOpt = 1;//pcShell option
	PetscFunctionBegin;
	iC(PetscOptionsGetInt(0,"-pcShell",&pcShellOpt,0));
	switch(pcShellOpt) {
	case 1: {
		iC(SetUpPC1((PC1ShellCtx*)(((PCShellCtx*)(pcSh))->pcShellStruct),B));
		break;
			}
	case 2: {
		iC(SetUpPC2((PC2ShellCtx*)(((PCShellCtx*)(pcSh))->pcShellStruct),B));
		break;
			}
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	default: SETERRQ(PETSC_ERR_USER,"Check the option -pcShell!");
#else
	default: SETERRQ(MPI_COMM_SELF, PETSC_ERR_USER, "Check the option -pcShell!");
#endif
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ApplyPCShell"
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 0))
PetscErrorCode ApplyPCShell(void* pcSh,Vec in,Vec out) {
	PetscInt pcShellOpt = 1;//pcShell option
	PetscFunctionBegin;
#else
PetscErrorCode ApplyPCShell(PC pc, Vec in, Vec out) {
	PetscInt pcShellOpt = 1;//pcShell option
	void* pcSh;
	PetscFunctionBegin;
	PCShellGetContext(pc, &pcSh);
#endif
	iC(PetscOptionsGetInt(0,"-pcShell",&pcShellOpt,0));
	switch(pcShellOpt) {
	case 1: {
		iC(ApplyPC1((PC1ShellCtx*)(((PCShellCtx*)(pcSh))->pcShellStruct),in,out));
		break;
			}
	case 2: {
		iC(ApplyPC2((PC2ShellCtx*)(((PCShellCtx*)(pcSh))->pcShellStruct),in,out));
		break;
			}
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	default: SETERRQ(PETSC_ERR_USER,"Check the option -pcShell!");
#else
	default: SETERRQ(MPI_COMM_SELF, PETSC_ERR_USER, "Check the option -pcShell!");
#endif
	}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DestroyPCShell"
PetscErrorCode DestroyPCShell(void* pcSh) {
	PetscInt pcShellOpt = 1;//pcShell option
	PetscFunctionBegin;
	iC(PetscOptionsGetInt(0,"-pcShell",&pcShellOpt,0));
	switch(pcShellOpt) {
	case 1: {
		iC(DestroyPC1((PC1ShellCtx*)(((PCShellCtx*)(pcSh))->pcShellStruct)));
		break;
			}
	case 2: {
		iC(DestroyPC2((PC2ShellCtx*)(((PCShellCtx*)(pcSh))->pcShellStruct)));
		break;
			}
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	default: SETERRQ(PETSC_ERR_USER,"Check the option -pcShell!");
#else
	default: SETERRQ(MPI_COMM_SELF, PETSC_ERR_USER, "Check the option -pcShell!");
#endif
	}
	PetscFunctionReturn(0);
}

#endif
