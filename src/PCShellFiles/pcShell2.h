///////////////////////////////////////////////////////////////////////////////////////
// pcShell2.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef __PC_SHELL_2_H
#define __PC_SHELL_2_H

//typedef struct{ Vec diag;} PC2ShellCtx;

#undef __FUNCT__
#define __FUNCT__ "CreatePC2Context"
PetscErrorCode CreatePC2Context(PC2ShellCtx** pc2Sh){
	PC2ShellCtx *tmpPc2;
	PetscFunctionBegin;
	//printf("Inside CreatePC2Context.\n");
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 4))
	iC(PetscNew(PC2ShellCtx,&tmpPc2));
#else
	iC(PetscNew(&tmpPc2));
#endif
	tmpPc2->diag =0;
	*pc2Sh = tmpPc2;
	//printf("Finished CreatePC2Context.\n");
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetUpPC2"
PetscErrorCode SetUpPC2(PC2ShellCtx* pc2Sh,Mat B ){
	Vec diag;
	int i,sizediag1;
	//int sizediag2;
	PetscReal *diagarray;
	PetscFunctionBegin;
	//printf("Inside SetUpPC2.\n");
	iC(MatGetVecs(B,&diag,0));
	iC(MatGetDiagonal(B,diag));
	//MatGetSize(B,&sizediag1,&sizediag2);
	VecGetSize(diag,&sizediag1);
	VecGetArray1d(diag,sizediag1,0, &diagarray);
	for (i=0;i<sizediag1;i++){if (fabs(diagarray[i])<=0.000001) {printf("Found zero on main diagonal sizediag=  i=  diagarray[i]=  %d %d %g\n",sizediag1,i,diagarray[i]);};}
	//if (sizediag1==3*9*9*9) {VecView(diag,0);}
	//printf("PC SIZEDIAG %d %d\n",sizediag1,sizediag2);
	iC(VecReciprocal(diag));
	pc2Sh->diag = diag;
	//printf("Finished SetUpPC2.\n");
	//if (sizediag1==3*17*17*17) {VecView(pc2Sh->diag,0);}
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ApplyPC2"
PetscErrorCode ApplyPC2(PC2ShellCtx* pc2Sh,Vec in,Vec out){
	int size;
	PetscFunctionBegin;
	//printf("Inside ApplyPC2.\n");
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)  
	//VecPointwiseMult(Vec x,Vec y, Vec w);
	iC(VecPointwiseMult(in,pc2Sh->diag,out));
#else
	//VecPointwiseMult(Vec w,Vec x,Vec y);
	iC(VecPointwiseMult(out,in,pc2Sh->diag));
#endif
	VecGetSize(out,&size);
	//VecView(pc2Sh->diag,0);
	//printf("PC SIZE %d\n",size);
	//printf("Finished ApplyPC2.\n");
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DestroyPC2"
PetscErrorCode DestroyPC2(PC2ShellCtx* pc2Sh){
	PetscFunctionBegin;
	//printf("Inside DestroyPC2.\n");
	if(pc2Sh->diag) {
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
		iC(VecDestroy(pc2Sh->diag));
#else
		iC(VecDestroy(&pc2Sh->diag));
#endif
	}
	iC(PetscFree(pc2Sh));
	//printf("Finished DestroyPC2.\n");
	PetscFunctionReturn(0);
}

#endif
