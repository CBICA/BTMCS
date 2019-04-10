///////////////////////////////////////////////////////////////////////////////////////
// rpHeader.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef __RP_HEADER_H
#define __RP_HEADER_H


typedef struct{stsDMMG dmmg;} IMFreeData;
static IMFreeData idata[100]; //static PetscInt iCnt =0;


static PetscInt iCnt; //CH, March 23, 2007

extern PetscErrorCode Interpolation1MatVec(Mat, Vec, Vec);
extern PetscErrorCode AddInterpolation1MatVec(Mat, Vec, Vec, Vec);

extern PetscErrorCode Restriction1MatVec(Mat, Vec, Vec);
extern PetscErrorCode AddRestriction1MatVec(Mat, Vec, Vec, Vec);
extern PetscErrorCode CreateInterpolation1(stsDMMG ,Mat *);
extern PetscErrorCode ComputeInterpolation1(stsDMMG ,Mat );

extern PetscErrorCode Restriction2MatVec(Mat, Vec, Vec);
extern PetscErrorCode AddRestriction2MatVec(Mat, Vec, Vec, Vec);
extern PetscErrorCode CreateInterpolation2(stsDMMG ,Mat *);
extern PetscErrorCode ComputeInterpolation2(stsDMMG ,Mat );

extern PetscErrorCode CreateInterpolationMatrixFree(stsDMMG ,Mat *);
extern PetscErrorCode ComputeInterpolationMatrixFree(stsDMMG ,Mat );

#endif
