///////////////////////////////////////////////////////////////////////////////////////
// ReadFiles.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include "function.h"

// Version 1.1 of STSolver by F. Abraham and G. Biros 
// August 2005


void material_props_finelevel(double *lambdaf, double *muf)
{
	int ie;
	FILE *fp, *fp1;

	fp  = fopen("materialmu.txt","r");
	fp1 = fopen("materiallambda.txt","r");

	for (ie=0;ie<ne_fine;ie++){fscanf(fp,"%lf",&muf[ie])    ; }
	for (ie=0;ie<ne_fine;ie++){fscanf(fp1,"%lf",&lambdaf[ie]);}

	printf("ne_fine %d\n",ne_fine);

	fclose(fp);
	fclose(fp1);
}

void levsetinit(double cx,double cy, double cz, double radius,double *lsnew)
{
	int counti,countj,countk,i,j,k,inc;
	double radius_force,xcoord,ycoord,zcoord;


	//  //Brain_Paper1 case
	//  cx=0.095;
	//  cy=0.07;
	//  cz=0.186/2.0;
	//  radius_force=0.0075;

	//  //Brain_Paper2 case
	//  cx=0.145;
	//  cy=0.18;
	//  cz=0.186/2.0;
	//  radius_force=0.006;

	//Dog DC1 case
	// cx=76*Hxf;
	// cy=72*Hyf;
	// cz=49*Hzf;
	// radius_force=0.003+max(Hzf,max(Hxf,Hyf));

	//Dog DC2 case
	//cx=72*Hxf;
	//cy=86*Hyf;
	//cz=59*Hzf;
	radius_force=radius+max(Hzf,max(Hxf,Hyf));

	//Dog DC3 case
	// cx=62*Hxf;
	// cy=47*Hyf;
	// cz=32*Hzf;
	// radius_force=0.003+max(Hzf,max(Hxf,Hyf));

	//Dog DC3 case, 64^3 image
	// cx=58.28/1000.0;
	// cy=44.18/1000.0;
	// cz=62.0/1000.0;
	// radius_force=0.003+max(Hzf,max(Hxf,Hyf));

	// //Human HC case
	//   cx=107.16/1000.0;
	//   cy=91.18/1000.0;
	//   cz=157.5/1000.0;
	//   //radius_force=0.0055+max(Hzf,max(Hxf,Hyf));
	//   radius_force=0.0055+(0.0165-0.0055)/2.0; //case edema

	counti=1;countj=1;countk=1;

	for (inc=0;inc<nnc_fine;inc++) {
		i = counti-1; j = countj-1; k = countk-1;  

		xcoord=i*Hxf;
		ycoord=j*Hyf;
		zcoord=k*Hzf;

		//start with a sphere of radius radius_force and center cx,cy,cz and initialize the level set accordingly
		//inc=i+j*mxf+k*mxf*myf;
		lsnew[inc]=sqrt((xcoord-cx)*(xcoord-cx)+(ycoord-cy)*(ycoord-cy)+(zcoord-cz)*(zcoord-cz))-radius_force;

		counti++;
		if (counti==mxf+1){ countj++; counti=1; }  
		if (countj==myf+1){ counti=1; countj=1; countk++; }  
	}
}

void penalized_neumann_boundary()
{
	int ineup,j;
	FILE *fp;

	normal_based_dirichlet_bc = 0;
	fp  = fopen("penalizedneumanndata.txt","r");
	fscanf(fp,"%d",&normal_based_dirichlet_bc);
	if (normal_based_dirichlet_bc == 0){
		fscanf(fp,"%d",&npenneupoints);
		fscanf(fp,"%lf",&Hpenalizedneumann);
		PetscMalloc(npenneupoints*sizeof(PetscReal),&penneupts);
		for (j=0;j<npenneupoints;j++) PetscMalloc(nsd*sizeof(PetscReal),&penneupts[j]);
		PetscMalloc(npenneupoints*sizeof(PetscReal),&penneubc);
		for (j=0;j<npenneupoints;j++) PetscMalloc(nsd*sizeof(PetscReal),&penneubc[j]);

		for (ineup =0;ineup<npenneupoints;ineup++){
			fscanf(fp,"%lf",&penneupts[ineup][0]);
			fscanf(fp,"%lf",&penneupts[ineup][1]);
			if ( nsd ==3 ) fscanf(fp,"%lf",&penneupts[ineup][2]);
		}
		for (ineup =0;ineup<npenneupoints;ineup++){
			fscanf(fp,"%lf",&penneubc[ineup][0]);
			fscanf(fp,"%lf",&penneubc[ineup][1]);
			if ( nsd ==3 ) fscanf(fp,"%lf",&penneubc[ineup][2]);
		}

		fclose(fp);
	} else {
		fscanf(fp,"%d",&npenneupoints);
		fscanf(fp,"%lf",&Hpenalizedneumann);
		PetscMalloc(npenneupoints*sizeof(PetscReal),&penneupts);
		for (j=0;j<npenneupoints;j++) PetscMalloc(nsd*sizeof(PetscReal),&penneupts[j]);
		PetscMalloc(npenneupoints*sizeof(PetscReal),&normal_vector);
		for (j=0;j<npenneupoints;j++) PetscMalloc(nsd*sizeof(PetscReal),&normal_vector[j]);
		PetscMalloc(npenneupoints*sizeof(PetscReal),&penneubc);
		for (j=0;j<npenneupoints;j++) PetscMalloc(1*sizeof(PetscReal),&penneubc[j]);

		for (ineup =0;ineup<npenneupoints;ineup++){
			fscanf(fp,"%lf",&penneupts[ineup][0]);
			fscanf(fp,"%lf",&penneupts[ineup][1]);
			if ( nsd ==3 ) fscanf(fp,"%lf",&penneupts[ineup][2]);
		}
		for (ineup =0;ineup<npenneupoints;ineup++){
			fscanf(fp,"%lf",&penneubc[ineup][0]);
		}
		for (ineup =0;ineup<npenneupoints;ineup++){
			fscanf(fp,"%lf",&normal_vector[ineup][0]);
			fscanf(fp,"%lf",&normal_vector[ineup][1]);
			if ( nsd ==3 ) fscanf(fp,"%lf",&normal_vector[ineup][2]);
		}
	}
}

