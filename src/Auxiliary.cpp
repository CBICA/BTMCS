///////////////////////////////////////////////////////////////////////////////////////
// Auxiliary.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include <string.h> // strcpy()
#include "global.h"


extern void InterpolateNodalToElement(Vec vnod, Vec velem); // ComputeQuant.c


#define SAVE_C
#define SAVE_DF
//#define SAVE_VF
//#define SAVe_LAMDA
//#define SAVE_MU
//#define SAVE_D


// Writes the ouput of the forward problem into separate files for c, disp, v, lambda, mu, D respectively;
// Output: binary files, here *float* format (can be double, if necessary, but it's 2x more expensive!)

#undef _FUNCT_
#define _FUNCT_ "OutputAuxiliary"
void OutputAuxiliary(Vec c, double *cumdisparray, Vec v, Vec lambda, Vec mu, Vec D, Vec phi, int counter)
{
	//stencil (box) created in Initialize
	int countx, county, countz, counte;
	int i,j,k,ic,ie,knl;
	PetscScalar *carray, *varray, *lambdarray, *muarray, *Darray, *carrayelem, *phiarray;
	Vec celem;
	//double res_x,res_y,res_z; //original image physical resolution in each dir. - eventually if deformed landmarks stored in a file

	//PetscScalar *cumdisparray;
	//Vec cumdispdup;

#ifdef SAVE_C
	float tmpfc;
	FILE *fp1, *fp1h;
#endif
#ifdef SAVE_DF
	float tmpfdisp;
	FILE *fp2, *fp2h;
#endif
#ifdef SAVE_VF
	float tmpfv;
	FILE *fp3, *fp3h;
#endif
#ifdef SAVe_LAMDA
	float tmpflam;
	FILE *fp4, *fp4h;
#endif
#ifdef SAVE_MU
	float tmpfmu;
	FILE *fp5, *fp5h;
#endif
#ifdef SAVE_D
	float tmpfdiff;
	FILE *fp6, *fp6h;
#endif
	//FILE *fp7;

	char hdrName[128];
	char imgName[128];
	char file_prefix_tum[128];
	char file_prefix_def[128];
	char file_prefix_vel[128];
	char file_prefix_lam[128];
	char file_prefix_mu[128];
	char file_prefix_diff[128];
	char file_prefix_gphi[128];

	strcpy(file_prefix_tum, "TumorDensity");
	strcpy(file_prefix_def, "DeformationField");
	strcpy(file_prefix_vel, "VelocityField");
	strcpy(file_prefix_lam, "Lambda");
	strcpy(file_prefix_mu, "Mu");
	strcpy(file_prefix_diff, "Diffusivity");
	strcpy(file_prefix_gphi, "Phase-field");

	// (fine) grid information (global, from Initialize)
	countx=mxf; //number of *nodes*
	county=myf;
	countz=mzf;

	counte=(countx-1)*(county-1)*(countz-1); //number of *elements* (fine grid)

	VecGetArray(c,&carray);    
	VecGetArray(v,&varray);  
	VecGetArray(lambda,&lambdarray);  
	VecGetArray(mu,&muarray);
	VecGetArray(D,&Darray);

	VecGetArray(phi,&phiarray);

#ifdef SAVE_C
	// tumor density header file
	sprintf(hdrName, "%s.%.3d.mhd", file_prefix_tum, counter);
	fp1h=fopen(hdrName, "w");
	sprintf(imgName, "%s.%.3d.fld", file_prefix_tum, counter);	                 
	fp1=fopen(imgName, "wb");
	fprintf(fp1h, "ObjectType = Image\nNDims = %d\nBinaryData = True\nBinaryDataByteOrderMSB = False\n", nsd);
	fprintf(fp1h, "TransformMatrix = %d %d %d %d %d %d %d %d %d\n",1,0,0,0,-1,0,0,0,1);
	fprintf(fp1h, "AnatomicalOrientation = RPI\n");
	if (imgdx == 0) {
		fprintf(fp1h, "ElementSpacing = %f %f %f\nDimSize = %d %d %d\n", 240.0/(countx-1), 240.0/(county-1), 192.0/(countz-1), countx-1, county-1, countz-1);
	} else {
		fprintf(fp1h, "ElementSpacing = %f %f %f\nDimSize = %d %d %d\n", imgdx, imgdy, imgdz, countx-1, county-1, countz-1);
	}
	fprintf(fp1h, "ElementNumberOfChannels = 1\nElementType = MET_FLOAT\nElementDataFile = %s\n", imgName);
	fclose(fp1h);
#endif

#ifdef SAVE_DF  
	// displacement header file
	sprintf(hdrName, "%s.%.3d.mhd", file_prefix_def, counter);
	fp2h=fopen(hdrName, "w");
	sprintf(imgName, "%s.%.3d.fld", file_prefix_def, counter);
	fp2=fopen(imgName, "wb"); 
	fprintf(fp2h, "ObjectType = Image\nNDims = %d\nBinaryData = True\nBinaryDataByteOrderMSB = False\n", nsd);
	if (imgdx == 0) {
		fprintf(fp2h, "ElementSpacing = %f %f %f\nDimSize = %d %d %d\n", 240.0/(countx-1), 240.0/(county-1), 192.0/(countz-1), countx, county, countz);
	} else {
		fprintf(fp2h, "ElementSpacing = %f %f %f\nDimSize = %d %d %d\n", imgdx, imgdy, imgdz, countx, county, countz);
	}
	fprintf(fp2h, "ElementNumberOfChannels = 3\nElementType = MET_FLOAT\nElementDataFile = %s\n", imgName);
	fclose(fp2h);
#endif

#ifdef SAVE_VF  
	// velocity header file
	sprintf(hdrName, "%s.%.3d.mhd", file_prefix_vel, counter);
	fp3h=fopen(hdrName, "w");
	sprintf(imgName, "%s.%.3d.fld", file_prefix_vel, counter);
	fp3=fopen(imgName, "wb");
	fprintf(fp3h, "ObjectType = Image\nNDims = %d\nBinaryData = True\nBinaryDataByteOrderMSB = False\n", nsd);
	if (imgdx == 0) {
		fprintf(fp3h, "ElementSpacing = %f %f %f\nDimSize = %d %d %d\n", 3.75, 3.75, 3.0, countx, county, countz);
	} else {
		fprintf(fp3h, "ElementSpacing = %f %f %f\nDimSize = %d %d %d\n", imgdx, imgdy, imgdz, countx, county, countz);
	}
	fprintf(fp3h, "ElementNumberOfChannels = 3\nElementType = MET_FLOAT\nElementDataFile = %s\n", imgName);
	fclose(fp3h);
#endif

#ifdef SAVE_LAMBDA  
	// lambda header file
	sprintf(hdrName, "%s.%.3d.mhd", file_prefix_lam, counter);
	fp4h=fopen(hdrName, "w");
	sprintf(imgName, "%s.%.3d.fld", file_prefix_lam, counter);	                 
	fp4=fopen(imgName, "wb"); 
	fprintf(fp4h, "ObjectType = Image\nNDims = %d\nBinaryData = True\nBinaryDataByteOrderMSB = False\n", nsd);
	if (imgdx == 0) {
		fprintf(fp4h, "ElementSpacing = %f %f %f\nDimSize = %d %d %d\n", 3.75, 3.75, 3.0, countx-1, county-1, countz-1);
	} else {
		fprintf(fp4h, "ElementSpacing = %f %f %f\nDimSize = %d %d %d\n", imgdx, imgdy, imgdz, countx-1, county-1, countz-1);
	}
	fprintf(fp4h, "ElementNumberOfChannels = 1\nElementType = MET_FLOAT\nElementDataFile = %s\n", imgName);
	fclose(fp4h);
#endif

#ifdef SAVE_MU
	// mu header file
	sprintf(hdrName, "%s.%.3d.mhd", file_prefix_mu, counter);
	fp5h=fopen(hdrName, "w");
	sprintf(imgName, "%s.%.3d.fld", file_prefix_mu, counter);	                 
	fp5=fopen(imgName, "wb"); 
	fprintf(fp5h, "ObjectType = Image\nNDims = %d\nBinaryData = True\nBinaryDataByteOrderMSB = False\n", nsd);
	if (imgdx == 0) {
		fprintf(fp5h, "ElementSpacing = %f %f %f\nDimSize = %d %d %d\n", 1.0, 1.0, 1.0, countx-1, county-1, countz-1);
	} else {
		fprintf(fp5h, "ElementSpacing = %f %f %f\nDimSize = %d %d %d\n", imgdx, imgdy, imgdz, countx-1, county-1, countz-1);
	}
	fprintf(fp5h, "ElementNumberOfChannels = 1\nElementType = MET_FLOAT\nElementDataFile = %s\n", imgName);
	fclose(fp5h);
#endif

#ifdef SAVE_D
	// D header file
	sprintf(hdrName, "%s.%.3d.mhd", file_prefix_diff, counter);
	fp6h=fopen(hdrName, "w");
	sprintf(imgName, "%s.%.3d.fld", file_prefix_diff, counter);	                 
	fp6=fopen(imgName, "wb"); 
	fprintf(fp6h, "ObjectType = Image\nNDims = %d\nBinaryData = True\nBinaryDataByteOrderMSB = False\n", nsd);
	if (imgdx == 0) {
		fprintf(fp6h, "ElementSpacing = %f %f %f\nDimSize = %d %d %d\n", 1.0, 1.0, 1.0, countx-1, county-1, countz-1);
	} else {
		fprintf(fp6h, "ElementSpacing = %f %f %f\nDimSize = %d %d %d\n", imgdx, imgdy, imgdz, countx-1, county-1, countz-1);
	}
	fprintf(fp6h, "ElementNumberOfChannels = 1\nElementType = MET_FLOAT\nElementDataFile = %s\n", imgName);
	fclose(fp6h);
#endif

	//sprintf(imgName, "%s.%.3d.fld", file_prefix_gphi, counter);
	//fp7=fopen(imgName, "wb"); 

	// save the *node-wise* quantities
	for (k=0;k<countz;k++) {
		for (j=0;j<county;j++) {
			for (i=0;i<countx;i++) { 

				ic=i+j*countx+k*countx*county; //current node   
				knl=nsd*ic;

				//save the tumor density *node-wise*; it's normalized (between 0 and ~1)

				//tmpfc=carray[ic];
				//printf("ic carray %d %f\n",ic,carray[ic]);
				//fwrite(&tmpfc, sizeof(float), 1, fp1);

				//tmpfc=phiarray[ic];
				//fwrite(&tmpfc, sizeof(float), 1, fp7);

#ifdef SAVE_DF
				// save displacement field, (x,y,z) order, here in *mm*!!

				tmpfdisp = (float)(cumdisparray[knl]*1000.0); // x-component
				fwrite(&tmpfdisp, sizeof(float), 1, fp2);
				tmpfdisp = (float)(cumdisparray[knl+1]*1000.0); // y-component
				fwrite(&tmpfdisp, sizeof(float), 1, fp2);
				tmpfdisp = (float)(cumdisparray[knl+2]*1000.0); // z-component
				fwrite(&tmpfdisp, sizeof(float), 1, fp2);
#endif

#ifdef SAVE_VF
				// save velocity field, (x,y,z) order, here in *mm/<time>*!!

				tmpfv=varray[knl]*1000.0; // x-component
				fwrite(&tmpfv, sizeof(float), 1, fp3);
				tmpfv=varray[knl+1]*1000.0; // y-component
				fwrite(&tmpfv, sizeof(float), 1, fp3);
				tmpfv=varray[knl+2]*1000.0; // z-component
				fwrite(&tmpfv, sizeof(float), 1, fp3);
#endif		
			}
		}
	}

	// save the *element-wise* quantities

	VecCreate(PETSC_COMM_WORLD,&celem);
	VecSetSizes(celem,PETSC_DECIDE,counte);
	VecSetFromOptions(celem);

	InterpolateNodalToElement(c,celem);
	VecGetArray(celem,&carrayelem); 

	for (k=0;k<countz-1;k++) {
		for (j=0;j<county-1;j++) {
			for (i=0;i<countx-1;i++) { 

				ie=i+j*(countx-1)+k*(countx-1)*(county-1); //current element 

#ifdef SAVE_C
				//save the tumor density *element-wise*; it's normalized (between 0 and ~1)

				tmpfc = (float)carrayelem[ie];
				//printf("ic carray %d %f\n",ic,carray[ic]);
				fwrite(&tmpfc, sizeof(float), 1, fp1);
#endif

#ifdef SAVE_LAMBDA
				//save lambda

				tmpflam=lambdarray[ie];
				fwrite(&tmpflam, sizeof(float), 1, fp4);
#endif

#ifdef SAVE_MU
				//save mu

				tmpfmu=muarray[ie];
				fwrite(&tmpfmu, sizeof(float), 1, fp5);
#endif

#ifdef SAVE_D
				//save D

				tmpfdiff=Darray[ie];
				fwrite(&tmpfdiff, sizeof(float), 1, fp6);
#endif
			}
		}
	}


#ifdef SAVE_C
	fclose(fp1);
#endif
#ifdef SAVE_DF
	fclose(fp2);
#endif
#ifdef SAVE_VF
	fclose(fp3);
#endif
#ifdef SAVE_LAMBDA
	fclose(fp4);
#endif
#ifdef SAVE_MU
	fclose(fp5);
#endif
#ifdef SAVE_D
	fclose(fp6);
#endif
	//fclose(fp7);

	VecRestoreArray(c,&carray);
	VecRestoreArray(v,&varray);
	VecRestoreArray(lambda,&lambdarray);
	VecRestoreArray(mu,&muarray);
	VecRestoreArray(D,&Darray);

	VecRestoreArray(phi,&phiarray);
	VecRestoreArray(celem,&carrayelem);
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	VecDestroy(celem);
#else
	VecDestroy(&celem);
#endif
} // end function OutputAuxiliary
