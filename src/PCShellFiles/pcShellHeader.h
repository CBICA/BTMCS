///////////////////////////////////////////////////////////////////////////////////////
// pcShellHeader.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef __PC_SHELL_HEADER_H
#define __PC_SHELL_HEADER_H

typedef struct {
	void *pcShellStruct; 
} PCShellCtx;

//typedef struct {} PC1ShellCtx;
typedef struct {char nothing[1];} PC1ShellCtx;

typedef struct{ Vec diag;} PC2ShellCtx;

extern PetscErrorCode CreatePC1Context(PC1ShellCtx**);
extern PetscErrorCode SetUpPC1(PC1ShellCtx*,Mat);
extern PetscErrorCode ApplyPC1(PC1ShellCtx*,Vec,Vec);
extern PetscErrorCode DestroyPC1(PC1ShellCtx*);

extern PetscErrorCode CreatePC2Context(PC2ShellCtx**);
extern PetscErrorCode SetUpPC2(PC2ShellCtx*,Mat );
extern PetscErrorCode ApplyPC2(PC2ShellCtx*,Vec,Vec);
extern PetscErrorCode DestroyPC2(PC2ShellCtx*);

extern PetscErrorCode CreatePCShellContext(void**);
extern PetscErrorCode SetUpPCShell(void*,Mat);
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 0))
extern PetscErrorCode ApplyPCShell(void*,Vec,Vec);
#else
extern PetscErrorCode ApplyPCShell(PC,Vec,Vec);
#endif
extern PetscErrorCode DestroyPCShell(void*);

#endif
