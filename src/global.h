///////////////////////////////////////////////////////////////////////////////////////
// global.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "petsc.h"
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
#include "petscda.h"
#else
#include "petscdmda.h"
#define PetscTruth PetscBool
#define EXTERN extern
#endif
#include "petscksp.h"
#include "petsclog.h"
#include "stsdamg.h"
#include "petscsnes.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "assert.h"


#ifndef _GLOBALVARS
#define _GLOBALVARS

#define MAX_DMMG_LEVELS 100
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif
#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef rint
//#define rint(a) (int)(a)
double rint(double x);
#endif

static int _internal_ierr = 0;
#define iC(fun) {_internal_ierr = fun; CHKERRQ(_internal_ierr);}


//*** global flag for tumor initialization


extern int jaccnt;

extern int TRUE, FALSE;
extern PetscReal eps5, eps16;
extern PetscReal one, mone, RTOD, pi, poisson_ratio;

extern int *ien;
extern int ndim, ndimx, ndimy, ndimz;
// 20120523 djk
extern int imgx, imgy, imgz;
extern double imgdx, imgdy, imgdz;
//
extern PetscReal Hx, Hy, Hz;
extern PetscReal Lx, Ly, Lz;
extern int xsd, ysd, zsd;

extern double T;
extern int ntimesteps,nstore;
extern PetscReal cfln;

extern int nblmark;


// dump them here temporarily - though it's not a good idea
extern double gstiffwm,gstiffgm,gstiffvent,gstiffcsf;
extern double gcompresbrain, gcompresvent;
extern double gdiffwm,gdiffgm,gdiffvent;

extern double grho, gp1, gp2, gse;
extern double gcinit, gcs;
extern double gxc, gyc, gzc, gsigsq;

extern double gres_x,gres_y,gres_z;


//these inputs have been added by Ali.
extern double gstop_mass;
extern char gfileTumorInput[800];
extern PetscTruth  tumorFileFlag;
extern PetscTruth  stopMassFlag;
//
extern char gfileInput[800];
extern char gfileInputLmarkUndef[800];
extern char gfileInputLmarkDef[800];

//extern Vec cobj;
extern Vec gphi;

extern int MatVecFlop;


/***********************************************************************/
/***********************************************************************/
/* Tytpe of equations used */
/***********************************************************************/
/***********************************************************************/
extern int lagrangian;
extern int dynamic; 	      /* dynamic elatsicity equation active ? */
extern int nonlinear; 	      /* using nonlinear/linear equations? */
extern int steady;  	      /* steady or unsteady ? */
extern int needlesimulation; /* Are force terms due to the needle */
extern int meshmoving;       /* Is the mesh moving or a fixed grid*/
extern int matrixfree;       /* Is linear solver used matrix free */
extern int ntimestep;        /* Number of timesteps */
extern int itimestep; 	      /* current timestep */
extern PetscReal dt; 	      /* what is the size of the time step */
extern PetscReal soltime;    /* What is the current solution time */
/* Factors for the dynamic elasticity equations */
extern PetscReal velocityalpha; 
extern PetscReal positiongamma; 
extern PetscReal density; 
extern PetscReal damping_factor;
/***********************************************************************/
/***********************************************************************/




/***********************************************************************/
/***********************************************************************/
/* Needle specific variables */
/***********************************************************************/
/***********************************************************************/
extern int num_needle_nodes; /* number of needle nodes */
extern int needle_rigid ;    /* is needle flexible oir rigid ? */
extern PetscReal **coord_needle_nodes;  /* Nodes of the needle */
extern PetscReal *needle_velocity_direction; /* Meedle velocity direction */
extern PetscReal **force_needle_nodes; /* force on needle nodes */
/***********************************************************************/
/***********************************************************************/



/***********************************************************************/
/***********************************************************************/
/* Testcase-specific variables */
/***********************************************************************/
/***********************************************************************/
extern int titlecase; /* what is the specific problem that we are solving;
					  code user defined factes depending on title */
extern int tchoice;
extern int nnabaqus;
extern PetscReal **truthcubed1, **truthcubed2; 
extern PetscReal **truthcubed3, **truthcubed4;
extern PetscReal *truthcubesol, *truthcubelem;
extern PetscReal **abaqusnode, **abaqusnode1;
extern PetscReal *abaqussol, *abaquselem;
/***********************************************************************/
/***********************************************************************/

/***********************************************************************/
/********************* Variables for elements and nodes ****************/
/*************** these are partition specific **************************/
/***********************************************************************/
extern int ne;  /* number of elements in the local partition */
extern int nen; /* number of nodes per element */
extern int nsd; /* number of spatial dimensions */
extern int ndf; /* number of degrees of freedom; for compressible
				elasticity will be eual to nsd */
extern int nquad; /* number of quadrature points for integration scheme
				  if used */
extern int blksize; /* nen * ndf */
extern int nnc; /* number of nodes within partition */
extern int rednquad; /* number of quadratures for reduced integration */
extern int nboundary; /* number of oundaries for the Cartesian grid */
extern int mgnlevels; /* Number of multigrid levels */
extern int mxg, myg, mzg; 
extern int nnc_fine; /* number of nodes within partiton for finest level */
extern int mxf, myf, mzf, xmf, ymf, zmf; /* Fine mesh details*/
extern int ne_fine;
/***********************************************************************/

extern PetscReal Hxf, Hyf, Hzf;


extern PetscReal *diriel,*dirivec;
extern PetscReal *forcenode;
extern PetscReal *hack1n, *hack1e;
extern Vec resid_corr, DiriVec;
extern FILE *filep;

/***********************************************************************/
/*********** Material Properties ***************************************/
/***********************************************************************/
extern int material_projection_required;
extern int inhomogenity_constant;
extern PetscReal *nodemu, *nodelam;
extern PetscReal *lambdaf, *muf, *diff;
extern PetscReal youngs_global,  youngs_min;
extern PetscReal youngs_global1, youngs_global2;
extern PetscReal youngs_global3, youngs_global4;
/***********************************************************************/


/***********************************************************************/
/**********************************DIRICHLET B.C. **********************/
/***********************************************************************/
extern int ndirpoints;
extern double **dirpts, **dirbc;
/***********************************************************************/

/***********************************************************************/
/********************************* NEUMANN B.C.*************************/
extern int neumann_bc_flag_on;
extern int nneupoints;
extern PetscReal **neupts; 
extern PetscReal **neubc;
extern PetscReal neumanneps, Hneumann;
/***********************************************************************/

/********************************* PENALIZED NEUMANN B.C.*************/
extern int penalized_neumann_bc_flag_on;
extern int normal_based_dirichlet_bc;
extern int npenneupoints;
extern PetscReal **penneupts; 
extern PetscReal **penneubc;
extern PetscReal **normal_vector;
extern PetscReal penalizedneumanneps, Hpenalizedneumann;
/***********************************************************************/

/**********************************************************************/
/************** Variables to denote the Cartesian boundaries ***********/
typedef struct {
	int idf, jdf, kdf;
} boundid;
extern boundid *boundaryid;
/**********************************************************************/
/**********************************************************************/

/*********************************************************************/
/****************  Matrix-free variables ******************************/
typedef struct{
	Vec vbac; 
	Vec lambdavec; 
	Vec muvec; 
	Vec vjacprec;
	stsDMMG dmmg;
	PetscReal vshxshx[8][8], vshxshy[8][8], vshxshz[8][8]; 
	PetscReal vshyshx[8][8], vshyshz[8][8], vshyshy[8][8];
	PetscReal vshzshy[8][8], vshzshz[8][8], vshzshx[8][8]; 
	PetscReal vshdshd[8][8];
	Mat J;
} stsData;
extern stsData shellData[MAX_DMMG_LEVELS], data1[MAX_DMMG_LEVELS];
extern int globallevelcnt, jacobiprecflag;
extern int levelmx[MAX_DMMG_LEVELS], levelmy[MAX_DMMG_LEVELS], levelmz[MAX_DMMG_LEVELS];
/*********************************************************************/
/*********************************************************************/

/**************************** Nonlinear material *********************/
typedef struct {
	PetscScalar u,v,w;
} Field;
extern int material_type;
/*********************************************************************/


#endif  /*_GLOBALVARS*/


