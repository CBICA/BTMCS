///////////////////////////////////////////////////////////////////////////////////////
// initialize.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include "function.h"


void readinput(int argc,char **argv)
{
	int i;

	PetscTruth  flg,PreLoad = PETSC_FALSE;

	//these inputs are made by Ali.
	PetscOptionsGetReal(0,"-gstop_mass",&gstop_mass,&stopMassFlag);
	PetscOptionsGetString(PETSC_NULL,"-gfileTumorInput",gfileTumorInput,800,&tumorFileFlag);
	//

	PetscOptionsGetReal(0,"-T",&T,0);
	PetscOptionsGetInt(0,"-ntimesteps",&ntimesteps,0);
	PetscOptionsGetInt(0,"-nstore",&nstore,0);

	PetscOptionsGetReal(0,"-cfln",&cfln,0);

	PetscOptionsGetReal(0,"-gstiffwm",&gstiffwm,0);
	PetscOptionsGetReal(0,"-gstiffgm",&gstiffgm,0);
	PetscOptionsGetReal(0,"-gstiffvent",&gstiffvent,0);
	PetscOptionsGetReal(0,"-gstiffcsf",&gstiffcsf,0);

	PetscOptionsGetReal(0,"-gdiffwm",&gdiffwm,0);
	PetscOptionsGetReal(0,"-gdiffgm",&gdiffgm,0);
	PetscOptionsGetReal(0,"-gdiffvent",&gdiffvent,0);

	PetscOptionsGetReal(0,"-gcompresbrain",&gcompresbrain,0);
	PetscOptionsGetReal(0,"-gcompresvent",&gcompresvent,0);

	PetscOptionsGetString(PETSC_NULL,"-gfileInput",gfileInput,800,&flg);
	PetscOptionsGetString(PETSC_NULL,"-gfileInputLmarkUndef",gfileInputLmarkUndef,800,&flg);
	PetscOptionsGetString(PETSC_NULL,"-gfileInputLmarkDef",gfileInputLmarkDef,800,&flg);

	PetscOptionsGetReal(0,"-gres_x",&gres_x,0);
	PetscOptionsGetReal(0,"-gres_y",&gres_y,0);
	PetscOptionsGetReal(0,"-gres_z",&gres_z,0);
	PetscOptionsGetInt(0,"-nblmark",&nblmark,0);

	PetscOptionsGetReal(0,"-grho",&grho,0);
	PetscOptionsGetReal(0,"-gp1",&gp1,0);
	PetscOptionsGetReal(0,"-gp2",&gp2,0);
	PetscOptionsGetReal(0,"-gse",&gse,0);

	PetscOptionsGetReal(0,"-gcinit",&gcinit,0);
	//PetscOptionsGetReal(0,"-gcs",&gcs,0);
	PetscOptionsGetReal(0,"-gxc",&gxc,0);
	PetscOptionsGetReal(0,"-gyc",&gyc,0);
	PetscOptionsGetReal(0,"-gzc",&gzc,0);
	PetscOptionsGetReal(0,"-gsigsq",&gsigsq,0);

	PetscOptionsGetInt(0,"-ndim",&ndim,0);
	PetscOptionsGetInt(0,"-matrixfree",&matrixfree,0);
	PetscOptionsGetInt(0,"-jacobiprec",&jacobiprecflag,0);
	PetscOptionsGetInt(0,"-titlecase",&titlecase,0);
	ndimx = ndim; ndimy = ndim; ndimz = ndim;
	PetscOptionsGetReal(0,"-youngs_min",&youngs_min,0);
	PetscOptionsGetReal(0,"-youngs_global",&youngs_global,0);
	PetscOptionsGetReal(0,"-density",&density,0);
	PetscOptionsGetReal(0,"-Lx",&Lx,0);
	PetscOptionsGetReal(0,"-Ly",&Ly,0);
	PetscOptionsGetReal(0,"-Lz",&Lz,0);
	PetscOptionsGetReal(0,"-poisson",&poisson_ratio,0);
	PetscOptionsGetInt(0,"-material_type",&material_type,0);
	PetscOptionsGetInt(0,"-inhomogenity_constant",&inhomogenity_constant,0);
	PetscOptionsGetInt(0,"-neumann_bc_flag_on",&neumann_bc_flag_on,0);
	PetscOptionsGetInt(0,"-material_projection_required",&material_projection_required,0);
	PetscOptionsGetReal(0,"-neumanneps",&neumanneps,0);
	PetscOptionsGetInt(0,"-penalized_neumann_bc_flag_on",&penalized_neumann_bc_flag_on,0);
	PetscOptionsGetReal(0,"-penalizedneumanneps",&penalizedneumanneps,0);

	PetscOptionsGetInt(0,"-tchoice",&tchoice,0);
	PetscOptionsGetInt(0,"-ndimx",&ndimx,0);
	PetscOptionsGetInt(0,"-ndimy",&ndimy,0);
	PetscOptionsGetInt(0,"-ndimz",&ndimz,0);
	// 20120523 djk
	imgx = imgy = imgz = 0;
	imgdx = imgdy = imgdz = 0;  
	PetscOptionsGetInt(0,"-imgx",&imgx,0);
	PetscOptionsGetInt(0,"-imgy",&imgy,0);
	PetscOptionsGetInt(0,"-imgz",&imgz,0);
	PetscOptionsGetReal(0,"-imgdx",&imgdx,0);
	PetscOptionsGetReal(0,"-imgdy",&imgdy,0);
	PetscOptionsGetReal(0,"-imgdz",&imgdz,0);
	//
	PetscOptionsGetInt(0,"-nsd",&nsd,0);
	if (nsd==2) nboundary = 4;
	if (nsd==3) nboundary = 6;
	PetscOptionsGetInt(0,"-needlesimulation",&needlesimulation,0);
	PetscOptionsGetReal(0,"-dt",&dt,0);
	PetscOptionsGetReal(0,"-positiongamma",&positiongamma,0);
	PetscOptionsGetReal(0,"-velocityalpha",&velocityalpha,0);
	PetscOptionsGetReal(0,"-damping_factor",&damping_factor,0);
	PetscMalloc(nboundary*sizeof(boundid),&boundaryid);
	for (i=0;i<nboundary;i++){
		boundaryid[i].idf = 0;
		boundaryid[i].jdf = 0;
		if (nsd ==3) boundaryid[i].kdf = 0;
	}
	PetscOptionsGetInt(0,"-lagrangian",&lagrangian,0);
	PetscOptionsGetInt(0,"-meshmoving",&meshmoving,0);
	PetscOptionsGetInt(0,"-steady",&steady,0);
	PetscOptionsGetInt(0,"-dynamic",&dynamic,0);
	PetscOptionsGetInt(0,"-ntimestep",&ntimestep,0);
	PetscOptionsGetInt(0,"-mgnlevels",&mgnlevels,0);
	// some of these will be invalid
	PetscOptionsGetInt(0,"-boundary1x",&boundaryid[0].idf,0);
	PetscOptionsGetInt(0,"-boundary1y",&boundaryid[0].jdf,0);
	PetscOptionsGetInt(0,"-boundary2x",&boundaryid[1].idf,0);
	PetscOptionsGetInt(0,"-boundary2y",&boundaryid[1].jdf,0);
	PetscOptionsGetInt(0,"-boundary3x",&boundaryid[2].idf,0);
	PetscOptionsGetInt(0,"-boundary3y",&boundaryid[2].jdf,0);
	PetscOptionsGetInt(0,"-boundary4x",&boundaryid[3].idf,0);
	PetscOptionsGetInt(0,"-boundary4y",&boundaryid[3].jdf,0);
	if (nsd ==3){
		PetscOptionsGetInt(0,"-boundary1z",&boundaryid[0].kdf,0);
		PetscOptionsGetInt(0,"-boundary2z",&boundaryid[1].kdf,0);
		PetscOptionsGetInt(0,"-boundary3z",&boundaryid[2].kdf,0);
		PetscOptionsGetInt(0,"-boundary4z",&boundaryid[3].kdf,0);
		PetscOptionsGetInt(0,"-boundary5x",&boundaryid[4].idf,0);
		PetscOptionsGetInt(0,"-boundary5y",&boundaryid[4].jdf,0);
		PetscOptionsGetInt(0,"-boundary5z",&boundaryid[4].kdf,0);
		PetscOptionsGetInt(0,"-boundary6x",&boundaryid[5].idf,0);
		PetscOptionsGetInt(0,"-boundary6y",&boundaryid[5].jdf,0);
		PetscOptionsGetInt(0,"-boundary6z",&boundaryid[5].kdf,0);
	}
	if (steady==0) soltime = 0.0;

	if (matrixfree == 1){
		/*PetscMalloc(nen*sizeof(PetscReal),&vshxshx);
		PetscMalloc(nen*sizeof(PetscReal),&vshxshy);
		PetscMalloc(nen*sizeof(PetscReal),&vshxshz);
		PetscMalloc(nen*sizeof(PetscReal),&vshyshx);
		PetscMalloc(nen*sizeof(PetscReal),&vshyshy);
		PetscMalloc(nen*sizeof(PetscReal),&vshyshz);
		PetscMalloc(nen*sizeof(PetscReal),&vshzshx);
		PetscMalloc(nen*sizeof(PetscReal),&vshzshy);
		PetscMalloc(nen*sizeof(PetscReal),&vshzshz);
		PetscMalloc(nen*sizeof(PetscReal),&vshdshd);
		for (i=0;i<nen;i++){
		PetscMalloc(nen*sizeof(PetscReal),&vshxshx[i]);
		PetscMalloc(nen*sizeof(PetscReal),&vshxshy[i]);
		PetscMalloc(nen*sizeof(PetscReal),&vshxshz[i]);
		PetscMalloc(nen*sizeof(PetscReal),&vshyshx[i]);
		PetscMalloc(nen*sizeof(PetscReal),&vshyshy[i]);
		PetscMalloc(nen*sizeof(PetscReal),&vshyshz[i]);
		PetscMalloc(nen*sizeof(PetscReal),&vshzshx[i]);
		PetscMalloc(nen*sizeof(PetscReal),&vshzshy[i]);
		PetscMalloc(nen*sizeof(PetscReal),&vshzshz[i]);
		PetscMalloc(nen*sizeof(PetscReal),&vshdshd[i]);
		}*/
	}

	if (nsd ==2){ndf =2;nen=4;nquad=4;}
	else if (nsd ==3){ndf =3;nen=8;nquad=8;}
	blksize = nen*ndf;
}

void initialize()
{
	pi               = 3.141592653589793238462643383276;
	RTOD             = 180.0/pi;
	rednquad         = 1;
	matrixfree       = 0;
	jacobiprecflag   = 0;
	tchoice   	   = 1;
	damping_factor   = 1.0;
	density          = 1000.0;
	needlesimulation = 0;
	lagrangian       = 0;
	meshmoving       = 0;
	eps5             = 0.00001;
	mone             = -1.0;
	one              = 1.0;
	TRUE             = 1;
	FALSE            = 0;
	xsd = 0;ysd =1;zsd =2;
	steady           = 1;
	soltime          = 1.0;
	dt               = 0.0;
	ntimestep        = 1;
	dynamic          = 0;
	positiongamma    = 1.0;
	velocityalpha    = 0.5;
	inhomogenity_constant = 1;
	Lx = 1.0;
	Ly = 1.0;
	Lz = 1.0;
	material_projection_required = 1;

	imgx = imgy = imgz = 0;
	imgdx = imgdy = imgdz = 0;  
}
