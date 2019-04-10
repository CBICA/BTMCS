///////////////////////////////////////////////////////////////////////////////////////
// main_forward_solver.c
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


//#define MODEL_OPT_CASE_1 //run forward model once to generate synthetic target tumor density
//#define MODEL_OPT_CASE_2 //run forward model once to generate synthetic target landmarks  


static char help[] = "Solves Forward Problem: 3D coupled diffusion/elasticity framework.\n\n";

void version()
{
	printf("==========================================================================\n");
	printf("ForwardSolverDiffusion\n");
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
	printf("Parameters:\n\n");
	printf("  -gfileInput <string>  Input segmented brain image (with or without tumor).\n");
	printf("  -T <double>           Represents the length of the time interval over\n");
	printf("                        which the tumor grows (e.g., number of days).\n");
	printf("  -gdiffgm <double>     Tumor cell diffusivity in the grey matter.\n");
	printf("  -gdiffwm <double>     Tumor cell diffusivity in the white matter.\n");
	printf("  -grho <double>        Tumor growth rate.\n");
	printf("  -gp1 <double>         Mass-effect parameters.\n");
	printf("  -gp2 <double>         Mass-effect parameters.\n");
	printf("  -gxc <double>         Initial tumor seed x-coordinate.\n");
	printf("  -gyc <double>         Initial tumor seed y-coordinate.\n");
	printf("  -gzc <double>         Initial tumor seed z-coordinate.\n");
	printf("\n");
	printf("For the complete list of parameters, check the software manual.\n");
	printf("\n");
	printf("Input parameters could be specified in ForwardSolver.in file or the command\n");
	printf("line arguments of the software. The ForwardSolver.in file must be in the\n");
	printf("current directory that executes the software.\n");
	printf("\n\n");
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	struct ConservationLaw cl;
	Vec m,v,c,cobj;
	Vec sm,sc,sphi;
	Vec *msep,lambda, mu, D, force;
	Vec disp,dispold,dispc;
	PetscScalar paramtum,paramforce[3];
	double *cumdisp,*vel, zero;
	int nblock,countxn,countyn,countzn,countn,countxe,countye,countze,counte;
	int n;
	PetscReal c_norm;
	PetscInt c_size;
	PetscInt ierr;

#ifdef MODEL_OPT_CASE_2 
	double f=0.0; //initialize output value of the objective function f

	FILE *fp1,*fp2;
	double *undeflmark, *deflmark, *deflmarkobj;
	double tmp1,tmp2;
	int jj,ii,i;
#endif

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

	printf("==========================================================================\n");
	printf("ForwardSolverDiffusion\n");
#ifdef SW_VER
	printf("  Version %s\n", SW_VER);
#endif
#ifdef SW_REV
	printf("  Revision %s\n", SW_REV);
#endif
	printf("Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.\n");
	printf("See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.\n");
	printf("==========================================================================\n");

	// initialize PETSC
	ierr = PetscInitialize(&argc, &argv, "ForwardSolver.in", help); CHKERRQ(ierr);

	//initialize context (DA, dmmg, etc.)
	InitializeForwardSolver(argc, argv);

	printf("ndimx ndimy ndimz dt %d %d %d %g\n", ndimx,ndimy,ndimz,dt);

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
	VecDuplicate(disp,&dispc);
	VecDuplicate(m,&sm);
	VecDuplicate(c,&sc);
	VecDuplicate(c,&sphi);
	VecDuplicateVecs(lambda,nblock,&msep);

	//printf("I am calling InitializeModel\n");
	InitializeModel(gfileInput,m,c,disp,v,cobj,gphi);

	printf("-Hxf %g\n",Hxf);
	printf("-Hyf %g\n",Hyf);
	printf("-Hzf %g\n",Hzf);
	printf("%s %d\n", gfileInput, ne_fine);

	paramtum=grho;
	paramforce[0]=gp1;
	paramforce[1]=gp2;
	paramforce[2]=gse;

	printf("rho p1 p2 s %g %g %g %g\n",paramtum,paramforce[0],paramforce[1],paramforce[2]);

	VecCopy(disp,dispc);
	VecGetArray(dispc,&cumdisp);

	/*// write original (initial) material properties in files: this is optional
	RetrieveVecs(m,msep);
	VecCopy(msep[0],lambda);
	VecCopy(msep[1],mu);
	VecCopy(msep[2],D);   
	OutputAuxiliary(c,cumdisp,v,lambda,mu,D,gphi,0); // this writes output into files*/


#ifdef MODEL_OPT_CASE_2 
	PetscMalloc(nblmark*nsd*sizeof(double),&undeflmark);
	PetscMalloc(nblmark*nsd*sizeof(double),&deflmark);
	PetscMalloc(nblmark*nsd*sizeof(double),&deflmarkobj);

	fp1=fopen(gfileInputLmarkUndef,"r"); //undeformed landmarks coordinates, *physical* , written in the input file as x,y,z
	for(jj=0;jj<nblmark;jj++) {
		//printf("lmark %d\n", jj);
		for(ii=0; ii<3;ii++) {
			fscanf(fp1,"%lf",&tmp1);

			if (ii==0){
				undeflmark[jj*3+ii]=tmp1*gres_x;
			}else if (ii==1){
				undeflmark[jj*3+ii]=tmp1*gres_y;
			}else if (ii==2){
				undeflmark[jj*3+ii]=tmp1*gres_z;
			}
			//printf("jj undeflmark %d %g\n",jj,undeflmark[jj*3+ii]);


		}//end for ii=0,2
	}//end for jj
	fclose(fp1);
#endif

	// this variable is intented to conditionally stop the loop
	c_norm = 0.0;
	c_size = 1;
	VecGetSize(c, &c_size);

	// start time loop for Tumor solver
	//
	// check the termination condition based on the input options
	// if stop_mass is specified, it will be given higher priority
	for (n = 1; n <= ntimesteps || (stopMassFlag && c_norm < gstop_mass); n++) {
		printf("time step %d\n", n);

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
		printf("done lambda\n");
		VecCopy(msep[1],mu);
		printf("done mu\n");
		VecCopy(msep[2],D);
		printf("done D\n"); 

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
		printf("done TumDiv\n");        
		(*(cl.setDiffusionOn))(&cl);
		(*(cl.transport))(&cl, c, D, sc, dt, cfln, c);
		printf("done TumDiff\n");
		ReactionForward(c,paramtum,dt,c);
		printf("done TumReact\n");

		//computing the total mass of tumor
		VecNorm(c,NORM_1,&c_norm);  
		c_norm = c_norm/((double)c_size);

		// once c updated, compute here the term "force" (Vec) to be used in SolveElasticity
		ComputeForceForward(c, paramforce, force);
		printf("done ComputeForce\n");

		//store previous displacement field, for computing velocity
		VecCopy(disp,dispold);

		// solve elasticity
		jaccnt=0;
		SolveElasticity(lambda, mu, force, disp);
		printf("done SolveElasticity\n");

		// update velocity field v (node-wise)
		ComputeVelocity(dispold, disp, dt, v);
		printf("done ComputeVelocity\n");

		// update the cumulative displacement field    
		VecGetArray(v,&vel);
		CumulateDisplacement(vel,dt,xmf,ymf,zmf,Hxf,Hyf,Hzf,cumdisp);
		VecRestoreArray(v,&vel);
		printf("done CumulateDisplacement\n");

		// store/save output every other "nstore" number of time steps
		if (fmod((double)n, (double)nstore) == 0) {
			OutputAuxiliary(c,cumdisp,v,lambda,mu,D,gphi, nstore > 0 ? n / nstore : n);
			printf("done Output\n");
		}
		printf("mass %f\n", c_norm);
	}

#ifdef MODEL_OPT_CASE_2  
	fp2=fopen("DefLmarkModel.txt","w"); //objective (deformed) landmarks coordinates, *physical* , written in the input file as x,y,z
	for(jj=0;jj<nblmark;jj++) {
		for(ii=0; ii<3;ii++) {
			if (ii==0){  
				fprintf(fp2,"%lf ",deflmark[jj*3+ii]);          
			}else if (ii==1){
				fprintf(fp2,"%lf ",deflmark[jj*3+ii]); 
			}else if (ii==2){  
				fprintf(fp2,"%lf \n",deflmark[jj*3+ii]);
			}           

		}
	}
	fclose(fp2);

	//if desired, compute the corresponding objective functional - Euclidean distance minimization, landmark-based

	//L1-norm (Euclidean!)
	for(jj=0;jj<nblmark;jj++) {
		f=f+sqrt((deflmark[jj*3]-deflmarkobj[jj*3])*(deflmark[jj*3]-deflmarkobj[jj*3])+\
			(deflmark[jj*3+1]-deflmarkobj[jj*3+1])*(deflmark[jj*3+1]-deflmarkobj[jj*3+1])+\
			(deflmark[jj*3+2]-deflmarkobj[jj*3+2])*(deflmark[jj*3+2]-deflmarkobj[jj*3+2]));
	}

	//re-scale f in *mm* here - it's in *m* originally, consistent with everything else!!
	f=1000.0*f;
	printf("Main: objective function f %g\n",f);

	/* //L2-norm (Euclidean!)
	for(jj=0;jj<nblmark;jj++) {
	f=f+((deflmark[jj*3]-deflmarkobj[jj*3])*(deflmark[jj*3]-deflmarkobj[jj*3])+\
	(deflmark[jj*3+1]-deflmarkobj[jj*3+1])*(deflmark[jj*3+1]-deflmarkobj[jj*3+1])+\
	(deflmark[jj*3+2]-deflmarkobj[jj*3+2])*(deflmark[jj*3+2]-deflmarkobj[jj*3+2]));
	}

	//re-scale f in *mm* here - it's in *m* originally, consistent with everything else!!
	f=1.0/2.0*1000.0*1000.0*f; //here mm^2!!

	printf("Main: objective function f %g\n",f); */

	PetscFree(undeflmark);
	PetscFree(deflmark);
	PetscFree(deflmarkobj);
#endif

	VecRestoreArray(dispc,&cumdisp);

#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	VecDestroyVecs(msep,nblock);
	VecDestroy(lambda);
	VecDestroy(mu);
	VecDestroy(D);
	VecDestroy(force);
	VecDestroy(disp);
	VecDestroy(dispold);
	VecDestroy(dispc);

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
	VecDestroy(&dispc);

	VecDestroy(&m);
	VecDestroy(&v);
	VecDestroy(&c);
	VecDestroy(&cobj);
	VecDestroy(&sm);
	VecDestroy(&sc);
#endif

	PetscFinalize();

	return 0;
}
