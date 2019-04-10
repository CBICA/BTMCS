///////////////////////////////////////////////////////////////////////////////////////
// function.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _FUNCTIONS
#define _FUNCTIONS

//#include "ConservLaw.h"

extern int ComputeJacobian(stsDMMG,Mat);
extern int AnComputeJacobian(stsDMMG,Mat);
extern int ComputeRHS(stsDMMG,Vec);
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
extern int RedComputeRHS(DMMG,Vec);
#else
extern int RedComputeRHS(stsDMMG,Vec);
#endif
void nodalforcevector(double*, int, int, int, double, double, double);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
extern int LagComputeJacobian(DMMG,Mat);
extern int LagComputeRHS(DMMG,Vec);
#else
extern int LagComputeJacobian(stsDMMG,Mat);
extern int LagComputeRHS(stsDMMG,Vec);
#endif
extern int lagrangian_init();

/*******************************************************/
/*  Embedded boundary conditions if any *****/
extern int NeumannContibution(stsDMMG,Vec);
void neumann_boundary();
extern int DirichletContibution(stsDMMG,Vec);
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
extern int Dirichlet_set(DMMG);
#else
extern int Dirichlet_set(stsDMMG);
#endif
/*******************************************************/

/*******************************************************/
/* Matrix Free Routines */
extern int ComputeMatFreeJacobian(Vec Vecin, Vec Vecout);
extern int CreateJacobian(stsDMMG,Mat *);
extern int CreateJacobianMatrixFree(stsDMMG dmmg,Mat *jac);
extern int ComputeJacobianMatrixFree(stsDMMG dmmg,Mat J);
extern int JacobianMatVecMatrixFree(Mat,Vec,Vec);
extern int JacobianMatGetDiagonal(Mat, Vec);
extern int matrixFreeMatVec(stsData *, Vec Vecin, Vec Vecout);
extern int ComputeJacobianMatDiagonal(stsData *,Vec Vecdiag);
/*******************************************************/

/*******************************************************/
/*********** Cartesian Dirichlet boundary identifiers **/
extern void dirichlet_boundary();
extern void dirichlet_identifier(int *dirichlet_flag, MatStencil *row,int mx, int my, int mz);
/*******************************************************/

/*******************************************************/
/***** Element level integration ***********************/
extern void quad(double **xq, double **sq, double ***sh, double *det,double *wq, double , double , double );
extern void quadlag(double **xq, double **sq, double ***sh, double *det,double *wq, int ie);
extern void analnode(double **sh);
/*******************************************************/

extern PetscErrorCode FormFunctionNLinear(SNES snes,Vec X,Vec F,void* ptr);
extern int matrixFreeMatVecNLinear(stsData *data, Vec Vecin, Vec Vecout);
extern int FormInitialGuess(SNES snes,Vec X,void *ptr);


extern void initialize();
extern void readinput(int argc,char **argv);

extern void denseelemmat(stsData*, int,int,int);
extern void rowcolindex(int i, int j, int k, MatStencil *row, MatStencil *col,int jacobian_flag);
extern void inhomogeneous_material(double*, double*, int, int , int, double, double, double );
extern void material_props_finelevel(double*, double*);
extern void createien(int xm, int ym);
extern void gatherreal(double*,  double*);
extern void gatherreal1(double*, double* );
extern void hacks();
extern int node_index(int rowi, int rowj, int rowk,int xd,int yd);
/*********************************************************************/
/* Test case specific *********************************************/
extern void testcases();
extern void truthcube_data();
extern void truthcube_compare();
extern void abaqus_data();
extern void abaqus_compare();
/*********************************************************************/

extern int jacobiprec(stsData *data);

/**********************************************************************/
/*Level Set Stuff*/
extern void marknarrowband(int*, double*, int, int, int);
extern void reinitlevset(double*,int, int, int);

/*Tumor stuff*/
extern double ComputeMass(double *disp, int xm, int ym, int zm);
extern double ComputeTumVolume(int xm, int ym, int zm);
extern void  OriginalMatProp(char InputImageFile[800], 
	double stiffwm, double stiffgm, double stiffvent, double stiffcsf,double diffwm, double diffgm, double diffvent,
	double compresbrain, double compresvent, double *lambda, double *mu, double *dif);
extern void  OriginalLandmarks(char InputLandmarkFileUndef[800],char InputLandmarkFileDef[800], double *undeflmark, double *deflmarkobj);
extern void DeformLandmarks(int nblmarks, double *undeflmark, double *vel, double dt, int xm, int ym, int zm, double *deflmark);


/****************************************************************/
/*Matrix allocation/freeing*/
extern double** my_malloc_matrix(int rows, int cols);
extern void my_free_matrix(double ** arr, int rows, int cols);


/*Diffusion 3D*/
extern void InterpolateVecNodalToElement(Vec vnod, Vec velem);
extern void InterpolateElementToNodal(Vec velem, Vec vnod);
extern void InterpolateNodalToElement(Vec vnod, Vec velem);
extern void RetrieveVecs(Vec m, Vec *msep);
extern void ReactionForward(Vec c0, double rho, double dt, Vec c);
//extern void AdvectionGrad(struct ConservationLaw *cl, Vec m0, Vec v, Vec s, double dt, double cfln,  Vec m);
//extern void AdvectionDiv(struct ConservationLaw *cl, Vec m0, Vec v, Vec s, double dt, double cfln,  Vec m);
//extern void DiffusionExplicit(struct ConservationLaw *cl, Vec c0, Vec D, Vec s, double dt, double cfln,  Vec c);
extern void ComputeForceForward(Vec c, double *paramforce, Vec force);
extern int SolveElasticity(Vec lambda, Vec mu, Vec force, Vec disp);
extern void ComputeVelocity(Vec dispold, Vec dispnew, double dt, Vec v);
extern double ComputeObjectiveFuncDensitySpace(Vec c, Vec cobj);
extern void TrackParticlesLagrangian(double *vel, double timespan, int xm, int ym, int zm, double Hx, double Hy, double Hz, double *deflmark);

extern int InitializeForwardSolver(int argc,char **argv);
extern void InitializeModel(char InputImageFile[800], Vec m, Vec c, Vec disp, Vec v, Vec cobj, Vec phi);
extern void OutputAuxiliary(Vec c, double *cumdisparray, Vec v, Vec lambda, Vec mu, Vec D, Vec phi, int counter);

extern void CumulateDisplacement(double *vel, double timespan, int xm, int ym, int zm, double Hx, double Hy, double Hz, double *dispc);

#define node_index3(rowi,rowj,rowk,xd,yd) ( (rowj)*(xd) + (rowi) + (rowk)*(xd)*(yd))

#endif 
