///////////////////////////////////////////////////////////////////////////////////////
// InitializeAll.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h> // exit()

#include "global.h"
#include "common.h"
#include "PCvar.h"

stsDMMG *dmmg;
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
DA da;
#else
DM da;
#endif
int ierr, ilevel, mglev;

PCShellCtx pcShell;
char *pcShellName = 0;
PetscInt pcShellOpt = 0;

extern void initialize(); // initialize.c
extern void readinput(int argc,char **argv); // initialize.c
extern void OriginalMatProp(char InputImageFile[800],
    double stiffwm, double stiffgm, double stiffvent, double stiffcsf,double diffwm, double diffgm, double diffvent,
    double compresbrain, double compresvent, double *lambda, double *mu, double *dif); // OriginalMatProp.c
extern void InterpolateElementToNodal(Vec velem, Vec vnod); // ComputeQuant.c
extern void InitializeTumorFromFile(char InputTumorImageFile[800], Vec c);

#undef __FUNCT__
#define __FUNCT__ "InitializeForwardSolver"
int InitializeForwardSolver(int argc,char **argv)
{
	PetscFunctionBegin;
	ierr=PetscOptionsGetInt(0,"-pcShell",&pcShellOpt,0);

	initialize();
	readinput(argc, argv);

	if (nsd == 3) {
		ierr = stsDMMGCreate(PETSC_COMM_WORLD,mgnlevels,PETSC_NULL,&dmmg);CHKERRQ(ierr);
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
		ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,ndimx,ndimy,ndimz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,3,1,0,0,0,&da); 
#elif (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 4))
		ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,ndimx,ndimy,ndimz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,3,1,0,0,0,&da); 
#else
		ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,ndimx,ndimy,ndimz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,3,1,0,0,0,&da); 
#endif
		//RS: This must be called before calling stsDMMGSetDM
		iC(stsDMMGSetInterpolationMatrixFree(dmmg,CreateInterpolationMatrixFree,ComputeInterpolationMatrixFree));
	} else {
		ierr = stsDMMGCreate(PETSC_COMM_WORLD,mgnlevels,PETSC_NULL,&dmmg);CHKERRQ(ierr);
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
		ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,ndimx,ndimy,PETSC_DECIDE,PETSC_DECIDE,2,1,0,0,&da); 
#elif (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 4))
		ierr = DMDACreate2d(PETSC_COMM_WORLD,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,DMDA_STENCIL_BOX,ndimx,ndimy,PETSC_DECIDE,PETSC_DECIDE,2,1,0,0,&da); 
#else
		ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,ndimx,ndimy,PETSC_DECIDE,PETSC_DECIDE,2,1,0,0,&da); 
#endif
	}
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	ierr = stsDMMGSetDM(dmmg,(DM)da);
	ierr = DADestroy(da);CHKERRQ(ierr);
#else
	ierr = stsDMMGSetDM(dmmg, da);
	ierr = DMDestroy(&da); CHKERRQ(ierr);
#endif

	/***************************************************************************************************/
	/***** Fine mesh statistics ************************************************************************/
	/***************************************************************************************************/
	mxf = (int)((ndimx-1) * pow((double)2, mgnlevels-1) + 1);
	myf = (int)((ndimy-1) * pow((double)2, mgnlevels-1) + 1);
	mzf = (int)((ndimz-1) * pow((double)2, mgnlevels-1) + 1);
	xmf = mxf;
	ymf = myf;
	zmf = mzf;

	Hxf = Lx/ (PetscReal)(mxf-1);
	Hyf = Ly/ (PetscReal)(myf-1);
	Hzf = Lz/ (PetscReal)(mzf-1);

	ne_fine = (mxf-1) * (myf-1);
	if (nsd == 3) ne_fine = ne_fine * (mzf-1);

	nnc_fine = xmf * ymf;
	if (nsd == 3) nnc_fine = nnc_fine * zmf;

	for (ilevel=0;ilevel<mgnlevels;ilevel++){
		//mglev 		= mgnlevels - ilevel;
		mglev 		= ilevel + 1;
		levelmx[ilevel] = (int)((ndimx-1) * pow((double)2, mglev-1) + 1);
		levelmy[ilevel] = (int)((ndimy-1) * pow((double)2, mglev-1) + 1);
		levelmz[ilevel] = (int)((ndimz-1) * pow((double)2, mglev-1) + 1);
	}

	/***************************************************************************************************/
	/***************************************************************************************************/

	dt=T/ntimesteps; //time-step; T and N shall be given as input parameters (global variables)

	PetscFunctionReturn(0);
}//end function InitializeForwardSolver

// Input: original image, segmented;
// Output: initial material properties (global arrays lambdaf, muf and diff) assigned element-wise
// - also stored in the Vec m for consistency; initial tumor density (Vec c), node-wise

#undef __FUNCT__
#define __FUNCT__ "InitializeModel"
void InitializeModel(char InputImageFile[800], Vec m, Vec c, Vec disp, Vec v, Vec cobj, Vec phi)
{
	int i,j,k,l,countxn,countyn,countzn,countn,ic,countxe,countye,countze,counte,ie,iel;
	int knl,nblock;
	PetscReal valc, valcobj, valm, valdispx, valdispy, valdispz, valvx, valvy, valvz, valphi;
	double zero, one, Hx, Hy, Hz, xcoord, ycoord, zcoord;
	Vec phielem;
	int size;

	//double *phielemarray;
	//float tmpf;
	//FILE *fp7;

	printf("ne_fine nnc_fine %d %d \n", ne_fine,nnc_fine);

	zero=0.0;
	one=1.0;

	//initial material properties: lambda, mu, D (on the *fine* mesh); 

	PetscMalloc(ne_fine*sizeof(double),&lambdaf); //globals
	PetscMalloc(ne_fine*sizeof(double),&muf);
	PetscMalloc(ne_fine*sizeof(double),&diff);

	OriginalMatProp(InputImageFile,gstiffwm,gstiffgm,gstiffvent,gstiffcsf,gdiffwm,gdiffgm,gdiffvent,gcompresbrain,gcompresvent, lambdaf, muf, diff);

	countxn=mxf;
	countyn=myf;
	countzn=mzf;

	countn=countxn*countyn*countzn;

	Hx=Hxf;
	Hy=Hyf;
	Hz=Hzf;

	countxe=countxn-1;
	countye=countyn-1;
	countze=countzn-1;

	//this line seems so suspicious to me
	//counte=countxe*countye*countxe;
	//Changed by Ali:
	counte=countxe*countye*countze;

	nblock=3;

	//VecCreate(PETSC_COMM_WORLD, &m);
	//printf("%d, %d\n", nblock, counte); 
	//VecSetSizes(m,PETSC_DECIDE,nblock*counte);
	//VecSetFromOptions(m);

	VecGetSize(m,&size);
	//printf("size %d\n",size);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	VecSet(&zero,m);
#else
	VecSet(m, zero);
#endif

	for (k=0;k<countze;k++) {
		for (j=0;j<countye;j++) {
			for (i=0;i<countxe;i++) { 

				ie=i+j*countxe+k*countxe*countye;

				//printf("ie %d\n", ie);

				for (l=0;l<nblock;l++) {

					iel=ie+l*counte;

					if (l==0) valm=lambdaf[ie];
					else if (l==1) valm=muf[ie];
					else if (l==2) valm=diff[ie];


					//if (diff[ie]>1.3e-9) printf("valm %g\n", valm);

					VecSetValue(m,iel,valm,INSERT_VALUES);

				} //end l loop
			}
		}
	}

	VecAssemblyBegin(m);
	VecAssemblyEnd(m); 

	//initial conditions on the tumor density c, elastic displacement u (disp) and velocity (v);

	//VecCreate(PETSC_COMM_WORLD,c);
	//VecSetSizes(*c,PETSC_DECIDE,countn);
	//VecSetFromOptions(*c);

	//VecCreate(PETSC_COMM_WORLD,disp);
	//VecSetSizes(*disp,PETSC_DECIDE,countn*nsd);
	//VecSetFromOptions(*disp);

	//VecCreate(PETSC_COMM_WORLD,v);
	//VecSetSizes(*v,PETSC_DECIDE,countn*nsd);
	//VecSetFromOptions(*v);

	//VecCreate(PETSC_COMM_WORLD,cobj);
	//VecSetSizes(*cobj,PETSC_DECIDE,countn);
	//VecSetFromOptions(*cobj);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	VecSet(&zero,c);
	VecSet(&zero,disp);
	VecSet(&zero,v);

	VecSet(&zero,cobj); // Vec cobj is *global* variable here
#else
	VecSet(c, zero);
	VecSet(disp, zero);
	VecSet(v, zero);

	VecSet(cobj, zero); // Vec cobj is *global* variable here
#endif

	//Changed by Ali:
	//To see whether we should initialize using a file or not:
	//see global.h
	if(!tumorFileFlag)
	{
		printf("Initilizing tumor from a Gaussian profile.\n");   
		for (k=0;k<countzn;k++) 
			for (j=0;j<countyn;j++) 
				for (i=0;i<countxn;i++) { 
					//set initial tumor density - try a Gaussian distribution with three parameters here
					xcoord=i*Hx;
					ycoord=j*Hy;
					zcoord=k*Hz;
					ic=i+j*countxn+k*countxn*countyn;
					valc=gcinit*exp(-((xcoord-gxc)*(xcoord-gxc)+(ycoord-gyc)*(ycoord-gyc)+(zcoord-gzc)*(zcoord-gzc))/(2*gsigsq));
					VecSetValue(c,ic,valc,INSERT_VALUES);
				}
	}
	else
	{
		printf("Initilizing tumor from tumor file:%s\n",gfileTumorInput);   
		InitializeTumorFromFile(gfileTumorInput, c); 
	} 

	for (k=0;k<countzn;k++) {
		for (j=0;j<countyn;j++) {
			for (i=0;i<countxn;i++) { 

				ic=i+j*countxn+k*countxn*countyn;
				knl=nsd*ic;

				//printf("ic %d\n", ic);

				//set initial tumor density - try a Gaussian distribution with three parameters here
				xcoord=i*Hx;
				ycoord=j*Hy;
				zcoord=k*Hz;

				//		valc=gcinit*exp(-((xcoord-gxc)*(xcoord-gxc)+(ycoord-gyc)*(ycoord-gyc)+(zcoord-gzc)*(zcoord-gzc))/(2*gsigsq));

				//printf("ic valc %d %g\n",ic,valc);


				valcobj=0.0;


				valdispx=0.0;
				valdispy=0.0;
				valdispz=0.0;
				valvx=0.0;
				valvy=0.0;
				valvz=0.0;

				//		VecSetValue(c,ic,valc,INSERT_VALUES);

				VecSetValue(cobj,ic,valcobj,INSERT_VALUES);

				VecSetValue(disp,knl,valdispx,INSERT_VALUES);

				VecSetValue(disp,knl+1,valdispy,INSERT_VALUES);

				VecSetValue(disp,knl+2,valdispz,INSERT_VALUES);

				VecSetValue(v,knl,valvx,INSERT_VALUES);

				VecSetValue(v,knl+1,valvy,INSERT_VALUES);

				VecSetValue(v,knl+2,valvz,INSERT_VALUES);

			}
		}
	}

	VecAssemblyBegin(c);
	VecAssemblyEnd(c); 

	VecAssemblyBegin(cobj);
	VecAssemblyEnd(cobj); 

	VecAssemblyBegin(disp);
	VecAssemblyEnd(disp); 

	VecAssemblyBegin(v);
	VecAssemblyEnd(v); 

	//initial condition on the phase-field variable phi introduced to track the ventricles
	//initially, based on material properties assignment, it will be defined element-wise;
	//however, since we need it *node-wise* inside the diffusion equation (advective+reaction terms)
	//will immediately interpolate to nodal from the very beginning
	//phi is set equal to 1 in the brain and zero in the ventricles

	VecCreate(PETSC_COMM_WORLD,&phielem);
	VecSetSizes(phielem,PETSC_DECIDE,counte);
	VecSetFromOptions(phielem);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	VecSet(&one,phielem);
#else
	VecSet(phielem, one);
#endif

	for (k=0;k<countze;k++) {
		for (j=0;j<countye;j++) {
			for (i=0;i<countxe;i++) { 

				ie=i+j*countxe+k*countxe*countye;

				if (diff[ie]==gdiffvent) {

					valphi=zero;  //printf("ie %d\n",ie);
					VecSetValue(phielem,ie,valphi,INSERT_VALUES);

				}
			}
		}
	}

	VecAssemblyBegin(phielem);
	VecAssemblyEnd(phielem); 

	/*fp7=fopen("InitialPhase-field", "wb");  
	VecGetArray(phielem,&phielemarray);
	for (k=0;k<countze;k++) {
	for (j=0;j<countye;j++) {
	for (i=0;i<countxe;i++) { 

	ie=i+j*countxe+k*countxe*countye;

	tmpf=phielemarray[ie];
	//printf("ie tmpf %d %g\n",ie,tmpf);
	fwrite(&tmpf, sizeof(float), 1, fp7);

	}
	}
	}
	fclose(fp7);
	VecRestoreArray(phielem,&phielemarray);*/

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	VecSet(&one,phi);
#else
	VecSet(phi, one);
#endif
	InterpolateElementToNodal(phielem, phi);

#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	VecDestroy(phielem);
#else
	VecDestroy(&phielem);
#endif
} //end function InitializeAll

#undef __FUNCT__
#define __FUNCT__ "InitilizeTumorFromFile"
/*Input: initial tumor distribution;
Output: Vec c, node-wise*/
void InitializeTumorFromFile(char InputTumorImageFile[800], Vec c)
{
	int ie;
	Vec celem;
	float valc;
	FILE  *fb;
	fb = fopen(InputTumorImageFile,"rb");//original segmented image
	if (fb == NULL) {
		perror("Could not open tumor initialization file\n");
		exit(0);
	}

	VecCreate(PETSC_COMM_WORLD,&celem);
	VecSetSizes(celem,PETSC_DECIDE,ne_fine);
	VecSetFromOptions(celem);


	for (ie=0;ie<ne_fine; ie++){

		fread(&valc,sizeof(float),1,fb);
		VecSetValue(celem,ie,valc,INSERT_VALUES);
	}
	//    printf("ie: %d, ne_fine: %d\n",ie,ne_fine);

	//these must be called right after VecSetValue  
	VecAssemblyBegin(celem);
	VecAssemblyEnd(celem); 
	//interpolating to node-wise Vec
	InterpolateElementToNodal(celem, c);
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	VecDestroy(celem);
#else
	VecDestroy(&celem);
#endif
}
