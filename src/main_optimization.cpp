///////////////////////////////////////////////////////////////////////////////////////
// main_optimization.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "function.h"
#include "ConservLaw.h"


//#define OPT_CASE_1 //if serial scans available, for tumor density-based optimization
#define OPT_CASE_2 //if serial scans/manual landmarks available, for landmark-based optimization


extern stsDMMG *dmmg;
extern int ierr;

static char help[] = "Solves Optimization Problem: 3D coupled diffusion/elasticity framework.\n\n";

void version()
{
	printf("==========================================================================\n");
	printf("LmarkObjective\n");
#ifdef SW_VER
	printf("  Version %s\n", SW_VER);
#endif
#ifdef SW_REV
	printf("  Revision %s\n", SW_REV);
#endif
	printf("Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.\n");
	printf("See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.\n");
	printf("==========================================================================\n");
}
void usage()
{
	printf("\n");
	printf("Usage:\n\n");
	printf("  LmarkObjective  -gfileInputLmarkUndef <string>  -gfileInputLmarkDef <string>\n");
	printf("                  [input]  [output]  [tag]  [tag_F]\n");
	printf("\n");
	printf("[input], [output], [tag], [tag_F] are appended by APPSPACK or HOPSPACK.");
	printf("\n");
	printf("For the complete list of parameters, check the software manual.\n");
	printf("\n");
	printf("Input parameters could be specified in ForwardSolver.in file or the command\n");
	printf("line arguments of the software. The ForwardSolver.in file must be in the\n");
	printf("current directory that executes the software.\n");
	printf("\n\n");
}


// n is the length of the optimization vector parameter x
// here x=(c0,D0,rho,p1) 
#undef __FUNCT__
#define __FUNCT__ "ForwardSolver_feval"
double ForwardSolver_feval(int n, double *x)
{
	struct ConservationLaw cl;
	Vec m,v,c,cobj;
	Vec sm,sc,sphi;
	Vec *msep,lambda, mu, D, force;
	Vec disp,dispold;
	PetscScalar paramtum,paramforce[3];
	double zero;
	int nblock,countxn,countyn,countzn,countn,countxe,countye,countze,counte;
	int ntimestep;

	double f=0.0; //initialize output value of the objective function f

	//FILE *fp1,*fp2;
#ifdef OPT_CASE_2
	double *undeflmark, *deflmark, *deflmarkobj;
	double *vel;
	int i, jj;
#endif

	// setup transport options
	InitConservationLaw(&cl);

	//setup initial conditions (t=0) for the model field variables
	nblock=3;
	countxn=mxf;
	countyn=myf;
	countzn=mzf;
	countn=countxn*countyn*countzn;
	countxe=countxn-1;
	countye=countyn-1;
	countze=countzn-1;
	counte=countxe*countye*countze;

	zero=0.0;

	//create all the work vectors required by the program
	VecCreate(PETSC_COMM_WORLD,&m);
	VecSetSizes(m,PETSC_DECIDE,nblock*counte);
	VecSetFromOptions(m);

	VecCreate(PETSC_COMM_WORLD,&c);
	VecSetSizes(c,PETSC_DECIDE,countn);
	VecSetFromOptions(c);

	VecCreate(PETSC_COMM_WORLD,&disp);
	VecSetSizes(disp,PETSC_DECIDE,countn*nsd);
	VecSetFromOptions(disp);

	VecCreate(PETSC_COMM_WORLD,&lambda);
	VecSetSizes(lambda,PETSC_DECIDE,counte);
	VecSetFromOptions(lambda);

	VecDuplicate(disp,&v);
	VecDuplicate(c,&cobj);
	VecDuplicate(c,&gphi);
	VecDuplicate(lambda,&mu);
	VecDuplicate(lambda,&D);
	VecDuplicate(disp,&force);
	VecDuplicate(disp,&dispold);
	VecDuplicate(m,&sm);
	VecDuplicate(c,&sc);
	VecDuplicate(c,&sphi);
	VecDuplicateVecs(lambda,nblock,&msep);

	gcinit=x[0];
	gdiffwm=x[1];
	gdiffgm=gdiffwm/5.0;
	grho=x[2];
	gp1=x[3];
	//gp2=x[4]; //if parameter p2 is also being considered for optimization

	paramtum=grho;
	paramforce[0]=gp1;
	paramforce[1]=gp2;
	paramforce[2]=gse;

	//printf("rho p1 p2 s %g %g %g %g\n",paramtum,paramforce[0],paramforce[1],paramforce[2]);

	InitializeModel(gfileInput,m,c,disp,v,cobj,gphi);

#ifdef OPT_CASE_2
	PetscMalloc(nblmark*nsd*sizeof(double),&undeflmark);
	PetscMalloc(nblmark*nsd*sizeof(double),&deflmark);
	PetscMalloc(nblmark*nsd*sizeof(double),&deflmarkobj);

	OriginalLandmarks(gfileInputLmarkUndef,gfileInputLmarkDef,undeflmark,deflmarkobj); //target landmarks for optimization
#endif

	// start time loop for Tumor solver
	for (ntimestep=1;ntimestep<=ntimesteps;ntimestep++){
		//printf("current time step %d\n",ntimestep);

		// advect material properties (elastic+diffusivity) (and eventually particle path psi if landmark-based approach)

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
		VecSet(&zero,sm);
#else
		VecSet(sm, zero);
#endif
		(*(cl.setGradAdvectionOn))(&cl);
		(*(cl.transport))(&cl, m, v, sm, dt, cfln, m);  

		// from the collective block vector m, retrieve the individual components (lambda, mu, D (and psi *iff* landmark-based optimization))
		RetrieveVecs(m,msep);

		VecCopy(msep[0],lambda);
		//printf("done lambda\n");
		VecCopy(msep[1],mu);
		//printf("done mu\n");
		VecCopy(msep[2],D);
		//printf("done D\n");	

		// advect phase-field variable phi for tracking the ventricles
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
		VecSet(&zero,sphi);
#else
		VecSet(sphi, zero);
#endif
		(*(cl.setGradAdvectionOn))(&cl);
		(*(cl.transport))(&cl, gphi, v, sphi, dt, cfln, gphi);  

		// solve tumor equation: fractional time steps advect->diffuse->react
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
		VecSet(&zero,sc);
#else
		VecSet(sc, zero);
#endif

		(*(cl.setDivAdvectionOn))(&cl);
		(*(cl.transport))(&cl, c, v, sc, dt, cfln, c); 

		//printf("done TumDiv\n");	 	

		(*(cl.setDiffusionOn))(&cl);
		(*(cl.transport))(&cl, c, D, sc, dt, cfln, c);

		//printf("done TumDiff\n");

		ReactionForward(c,paramtum,dt,c);

		//printf("done TumReact\n");

		// once c updated, compute here the term "force" (Vec) to be used in SolveElasticity
		ComputeForceForward(c, paramforce, force);

		//printf("done ComputeForce\n");

		//store previous displacement field, for computing velocity
		VecCopy(disp,dispold);

		// solve elasticity
		jaccnt=0;
		SolveElasticity(lambda, mu, force, disp);

		//printf("done SolveElasticity\n");

		// update velocity field v (node-wise)
		ComputeVelocity(dispold, disp, dt, v);

		//printf("done ComputeVelocity\n");

		/*#ifdef OPT_CASE_1
		spaceintobj=ComputeObjectiveFuncDensitySpace(c,cobj); //cobj is set in InitializeModel!
		printf("done ComputeObjFunc\n");
		#endif*/

#ifdef OPT_CASE_2
		// compute the model-estimated displacement of the manually placed landmarks
		VecGetArray(v,&vel);

		DeformLandmarks(nblmark, undeflmark, vel, dt, mxf, myf, mzf, deflmark);

		for(i=0;i<nblmark*nsd;i++){undeflmark[i]=deflmark[i];}

		VecRestoreArray(v,&vel);
#endif
	}//end time loop

#ifdef OPT_CASE_2 
	//update the objective functional - Euclidean distance minimization, landmark-based

	//L1-norm (Euclidean!)
	for(jj=0;jj<nblmark;jj++) {
		f=f+sqrt((deflmark[jj*3]-deflmarkobj[jj*3])*(deflmark[jj*3]-deflmarkobj[jj*3])+\
			(deflmark[jj*3+1]-deflmarkobj[jj*3+1])*(deflmark[jj*3+1]-deflmarkobj[jj*3+1])+\
			(deflmark[jj*3+2]-deflmarkobj[jj*3+2])*(deflmark[jj*3+2]-deflmarkobj[jj*3+2]));
	}

	//re-scale f in *mm* here - it's in *m* originally, consistent with everything else!!
	f=1000.0*f;

	//printf("Main: objective function f %g\n",f);

	/* //L2-norm (Euclidean!)
	for(jj=0;jj<nblmark;jj++) {

	f=f+((deflmark[jj*3]-deflmarkobj[jj*3])*(deflmark[jj*3]-deflmarkobj[jj*3])+\
	(deflmark[jj*3+1]-deflmarkobj[jj*3+1])*(deflmark[jj*3+1]-deflmarkobj[jj*3+1])+\
	(deflmark[jj*3+2]-deflmarkobj[jj*3+2])*(deflmark[jj*3+2]-deflmarkobj[jj*3+2]));

	}

	//re-scale f in *mm* here - it's in *m* originally, consistent with everything else!!
	f=1.0/2.0*1000.0*1000.0*f; //here mm^2!!

	printf("Main: objective function f %g\n",f); */
#endif


#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	VecDestroyVecs(msep,nblock);
	VecDestroy(lambda);
	VecDestroy(mu);
	VecDestroy(D);
	VecDestroy(force);
	VecDestroy(disp);
	VecDestroy(dispold);

	VecDestroy(m);
	VecDestroy(v);
	VecDestroy(c);
	VecDestroy(cobj);
	VecDestroy(sm);
	VecDestroy(sc);
#else
	VecDestroyVecs(nblock, &msep);
	VecDestroy(&lambda);
	VecDestroy(&mu);
	VecDestroy(&D);
	VecDestroy(&force);
	VecDestroy(&disp);
	VecDestroy(&dispold);

	VecDestroy(&m);
	VecDestroy(&v);
	VecDestroy(&c);
	VecDestroy(&cobj);
	VecDestroy(&sm);
	VecDestroy(&sc);
#endif

#ifdef OPT_CASE_2
	PetscFree(undeflmark);
	PetscFree(deflmark);
	PetscFree(deflmarkobj);
#endif

	return(f);
}//end function ForwardSolver_feval

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	//double *undeflmark, *deflmark;
	//double *testdeflmark;
	int n;			// number of dimensions
	double* x;		// n-dimensional vector
	double y;		// solution of f(x)
#ifdef OPT_CASE_2
	FILE* fp;		// file pointer
	int i;
#endif
	PetscReal timebegin, timeend;

	// parse command line
	{
		if (argc == 1) {
			printf("use option -h or --help for help\n");
			exit(EXIT_FAILURE);
		}
		if (argc == 2) {
			if        (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
				version();
				usage();
				exit(EXIT_FAILURE);
			} else if (strcmp(argv[1], "-u") == 0 || strcmp(argv[1], "--usage") == 0) {
				version();
				usage();
				exit(EXIT_FAILURE);
			} else if (strcmp(argv[1], "-V") == 0 || strcmp(argv[1], "--version") == 0) {
				version();
				exit(EXIT_FAILURE);
			} else {
				printf("error: wrong arguments\n");
				printf("use option -h or --help for help\n");
				exit(EXIT_FAILURE);
			}
		}
	}

#ifdef OPT_CASE_2 
	//Check input arguments 
	//if (argc < 2) {
	//	fprintf (stderr, "usage: %s <input file> <output file>\n", argv[0]);
	//	return -1;
	//}

	//Open input file 

	// if the 2 landmark files (original and target) need to be specified in the .apps file
	// (important particularly for the cross-validation tests)
	if ((fp = fopen(argv[5], "r")) == NULL) {
		fprintf(stderr, "%s - Error opening input file %s.\n", argv[4], argv[5]);
		return -1;
	}

	//Read size of x 

	if ((fscanf(fp, "%d", &n)) != 1) {
		fprintf(stderr, "%s - Error reading n.\n", argv[4]);
		return -1;
	}

	//Allocate memory for x 

	if ((x = (double*)malloc(n * sizeof(double))) == NULL) {
		fprintf(stderr, "%s - Error allocating space for x.\n", argv[4]);
		return -1;
	}

	//Read x 

	for (i = 0; i < n; i ++)
		if ((fscanf(fp, "%le", &x[i])) != 1) {
			fprintf(stderr, "%s - Error reading x[%d].\n", argv[4], i);
			return -1;
		}

	//Close input file 

	fclose(fp);

	/*//else
	if ((fp = fopen(argv[1], "r")) == NULL) {
	fprintf(stderr, "%s - Error opening input file %s.\n", argv[0], argv[1]);
	return -1;
	}

	//Read size of x 

	if ((fscanf(fp, "%d", &n)) != 1) {
	fprintf(stderr, "%s - Error reading n.\n", argv[0]);
	return -1;
	}

	//Allocate memory for x 

	if ((x = malloc(n * sizeof(double))) == NULL) {
	fprintf(stderr, "%s - Error allocating space for x.\n", argv[0]);
	return -1;
	}

	//Read x 

	for (i = 0; i < n; i ++)
	if ((fscanf(fp, "%le", &x[i])) != 1) {
	fprintf(stderr, "%s - Error reading x[%d].\n", argv[0], i);
	return -1;
	}

	//Close input file 

	fclose(fp);*/
#endif

	printf("==========================================================================\n");
	printf("LmarkObjective\n");
#ifdef SW_VER
	printf("  Version %s\n", SW_VER);
#endif
#ifdef SW_REV
	printf("  Revision %s\n", SW_REV);
#endif
	printf("Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.\n");
	printf("See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.\n");
	printf("==========================================================================\n");

	/* Evaluate function at x */

	// initialize PETSC
	PetscInitialize(&argc,&argv,"ForwardSolver.in",help);

	//initialize context (DA, dmmg, etc.)
	InitializeForwardSolver(argc,argv);

	//printf("ndimx ndimy ndimz %d %d %d\n", ndimx,ndimy,ndimz);

	//printf("Entering Forward Solver\n");
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 3))
	ierr = PetscGetTime(&timebegin);CHKERRQ(ierr);
#else
	ierr = PetscTime(&timebegin); CHKERRQ(ierr);
#endif

	y=ForwardSolver_feval(n,x);

#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 3))
	ierr = PetscGetTime(&timeend);CHKERRQ(ierr);
#else
	ierr = PetscTime(&timeend); CHKERRQ(ierr);
#endif
	timeend = timeend-timebegin;
	//printf("Exiting Forward Solver \n");
	//ierr = PetscPrintf(PETSC_COMM_WORLD,"Time of execution %g\n",timeend);CHKERRQ(ierr); 

	ierr = stsDMMGDestroy(dmmg);CHKERRQ(ierr);
	ierr = PetscFinalize();CHKERRQ(ierr);
	//printf("Finished with PETSC\n"); 

#ifdef OPT_CASE_2 
	//Open output file 

	// if the 2 landmark files (original and target) need to be specified in the .apps file
	// (important particularly for the cross-validation tests)
	if ((fp = fopen(argv[6], "w")) == NULL) {
		fprintf(stderr, "%s - Error opening output file.\n", argv[4]);
		return -1;
	}

	/*//else
	if ((fp = fopen(argv[2], "w")) == NULL) {
	fprintf(stderr, "%s - Error opening output file.\n", argv[0]);
	return -1;
	}*/

	//Write function value to output file 

	fprintf(fp, "%e\n", y);

	//Close output file 

	fclose(fp);


	//Re-save the input in file (APPSPACK can overwrite it!)

	//this is if 2 lmark. file specified in the .apps file; else must be replaced with   fp=fopen(argv[1],"w");
	fp=fopen(argv[5],"w"); 

	for (i = 0; i < n; i ++) {fprintf(fp,"%e\n",x[i]);}

	fclose(fp);

	//Release memory 
	free(x);
#endif

	return 0;
}

