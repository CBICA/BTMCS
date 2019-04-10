///////////////////////////////////////////////////////////////////////////////////////
// global.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"


int jaccnt;

int TRUE, FALSE;
PetscReal eps5, eps16;
PetscReal one, mone, RTOD, pi, poisson_ratio;

int *ien;
int ndim, ndimx, ndimy, ndimz;
// 20120523 djk
int imgx, imgy, imgz;
double imgdx, imgdy, imgdz;
//
PetscReal Hx, Hy, Hz;
PetscReal Lx, Ly, Lz;
int xsd, ysd, zsd;

double T;
int ntimesteps,nstore;
PetscReal cfln;

int nblmark;


// dump them here temporarily - though it's not a good idea
double gstiffwm,gstiffgm,gstiffvent,gstiffcsf;
double gcompresbrain, gcompresvent;
double gdiffwm,gdiffgm,gdiffvent;

double grho, gp1, gp2, gse;
double gcinit, gcs;
double gxc, gyc, gzc, gsigsq;

double gres_x,gres_y,gres_z;


//these inputs have been added by Ali.
double gstop_mass;
char gfileTumorInput[800];
PetscTruth tumorFileFlag;
PetscTruth stopMassFlag;
//
char gfileInput[800];
char gfileInputLmarkUndef[800];
char gfileInputLmarkDef[800];

//Vec cobj;
Vec gphi;

int MatVecFlop;


/***********************************************************************/
/***********************************************************************/
/* Tytpe of equations used */
/***********************************************************************/
/***********************************************************************/
int lagrangian;
int dynamic; 	      /* dynamic elatsicity equation active ? */
int nonlinear; 	      /* using nonlinear/linear equations? */
int steady;  	      /* steady or unsteady ? */
int needlesimulation; /* Are force terms due to the needle */
int meshmoving;       /* Is the mesh moving or a fixed grid*/
int matrixfree;       /* Is linear solver used matrix free */
int ntimestep;        /* Number of timesteps */
int itimestep; 	      /* current timestep */
PetscReal dt; 	      /* what is the size of the time step */
PetscReal soltime;    /* What is the current solution time */
/* Factors for the dynamic elasticity equations */
PetscReal velocityalpha; 
PetscReal positiongamma; 
PetscReal density; 
PetscReal damping_factor;
/***********************************************************************/
/***********************************************************************/




/***********************************************************************/
/***********************************************************************/
/* Needle specific variables */
/***********************************************************************/
/***********************************************************************/
int num_needle_nodes; /* number of needle nodes */
int needle_rigid ;    /* is needle flexible oir rigid ? */
PetscReal **coord_needle_nodes;  /* Nodes of the needle */
PetscReal *needle_velocity_direction; /* Meedle velocity direction */
PetscReal **force_needle_nodes; /* force on needle nodes */
/***********************************************************************/
/***********************************************************************/



/***********************************************************************/
/***********************************************************************/
/* Testcase-specific variables */
/***********************************************************************/
/***********************************************************************/
int titlecase; /* what is the specific problem that we are solving;
			   code user defined factes depending on title */
int tchoice;
int nnabaqus;
PetscReal **truthcubed1, **truthcubed2; 
PetscReal **truthcubed3, **truthcubed4;
PetscReal *truthcubesol, *truthcubelem;
PetscReal **abaqusnode, **abaqusnode1;
PetscReal *abaqussol, *abaquselem;
/***********************************************************************/
/***********************************************************************/

/***********************************************************************/
/********************* Variables for elements and nodes ****************/
/*************** these are partition specific **************************/
/***********************************************************************/
int ne;  /* number of elements in the local partition */
int nen; /* number of nodes per element */
int nsd; /* number of spatial dimensions */
int ndf; /* number of degrees of freedom; for compressible
		 elasticity will be eual to nsd */
int nquad; /* number of quadrature points for integration scheme
		   if used */
int blksize; /* nen * ndf */
int nnc; /* number of nodes within partition */
int rednquad; /* number of quadratures for reduced integration */
int nboundary; /* number of oundaries for the Cartesian grid */
int mgnlevels; /* Number of multigrid levels */
int mxg, myg, mzg; 
int nnc_fine; /* number of nodes within partiton for finest level */
int mxf, myf, mzf, xmf, ymf, zmf; /* Fine mesh details*/
int ne_fine;
/***********************************************************************/

PetscReal Hxf, Hyf, Hzf;


PetscReal *diriel,*dirivec;
PetscReal *forcenode;
PetscReal *hack1n, *hack1e;
Vec resid_corr, DiriVec;
FILE *filep;

/***********************************************************************/
/*********** Material Properties ***************************************/
/***********************************************************************/
int material_projection_required;
int inhomogenity_constant;
PetscReal *nodemu, *nodelam;
PetscReal *lambdaf, *muf, *diff;
PetscReal youngs_global,  youngs_min;
PetscReal youngs_global1, youngs_global2;
PetscReal youngs_global3, youngs_global4;
/***********************************************************************/


/***********************************************************************/
/**********************************DIRICHLET B.C. **********************/
/***********************************************************************/
int ndirpoints;
double **dirpts, **dirbc;
/***********************************************************************/

/***********************************************************************/
/********************************* NEUMANN B.C.*************************/
int neumann_bc_flag_on;
int nneupoints;
PetscReal **neupts; 
PetscReal **neubc;
PetscReal neumanneps, Hneumann;
/***********************************************************************/

/********************************* PENALIZED NEUMANN B.C.*************/
int penalized_neumann_bc_flag_on;
int normal_based_dirichlet_bc;
int npenneupoints;
PetscReal **penneupts; 
PetscReal **penneubc;
PetscReal **normal_vector;
PetscReal penalizedneumanneps, Hpenalizedneumann;
/***********************************************************************/

/**********************************************************************/
/************** Variables to denote the Cartesian boundaries ***********/
boundid *boundaryid;
/**********************************************************************/
/**********************************************************************/

/*********************************************************************/
/****************  Matrix-free variables ******************************/

stsData shellData[MAX_DMMG_LEVELS], data1[MAX_DMMG_LEVELS];
int globallevelcnt, jacobiprecflag;
int levelmx[MAX_DMMG_LEVELS], levelmy[MAX_DMMG_LEVELS], levelmz[MAX_DMMG_LEVELS];
/*********************************************************************/
/*********************************************************************/

/**************************** Nonlinear material *********************/
int material_type;
/*********************************************************************/
