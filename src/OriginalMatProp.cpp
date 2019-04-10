///////////////////////////////////////////////////////////////////////////////////////
// OriginalMatProp.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h> // exit()
#include "global.h"


//this function assignes material properties based on a given 3D segmented image

//InputImageFile: 3D segmented image, properly resampled to the size of the desired regular grid for FE solution of displacement, here
//given by the global variables mxf,myf,mzf

//Material properties (linear elastic): stiffness E and compressibility nu - white matter, gray matter, CSF and ventricles

//Initial tumor parameters: location (tumcx,tumcy,tumcz) and size (tumradius)

void  OriginalMatProp(char InputImageFile[800], 
	double stiffwm, double stiffgm, double stiffvent, double stiffcsf,double diffwm, double diffgm, double diffvent,
	double compresbrain, double compresvent, double *lambda, double *mu, double *dif)
{

	int counti, countj, countk, ie;
	double stiffskull, diffskull, youngsfmodulus, pr, diffcoef;
	unsigned char grayscale;
	FILE  *fb;

	neumanneps=0.01;
	stiffskull=100.0*stiffgm;
	diffskull=0.01*diffgm;

	//printf("-Hxf %g\n",Hxf);
	//printf("-Hyf %g\n",Hyf);
	//printf("-Hzf %g\n",Hzf);

	fb = fopen(InputImageFile,"rb");//original segmented image
	if (fb == NULL) {
		perror("Could not open file\n");
		exit(0);
	}

	//printf("%s %d\n", InputImageFile, ne_fine);

	counti = 1; countj = 1; countk = 1; 

	for (ie=0;ie<ne_fine; ie++)
	{
		//printf("%d\r", ie);

		fread(&grayscale,sizeof(char),1,fb);
		//fscanf(fb,"%1c",&grayscale);
		//if(grayscale)
		//printf("ie grayscale %d %c\n",ie,grayscale);

		if (grayscale == 250){
			youngsfmodulus = stiffwm; /* white matter */
			pr=compresbrain;
			diffcoef=diffwm;		
		}else if (grayscale == 150){
			youngsfmodulus = stiffgm; /* gray matter */
			pr=compresbrain;
			diffcoef=diffgm;
		}else if (grayscale == 50){
			youngsfmodulus = stiffvent; /* venticular CSF*/
			pr=compresvent;
			diffcoef=diffvent;
		}else if (grayscale == 10){
			youngsfmodulus = stiffcsf ; /* CSF*/
			pr=compresbrain;
			//diffcoef=diffwm;
			diffcoef=diffvent;
		}else if (grayscale == 0){
			youngsfmodulus = stiffskull; /* skull*/
			pr=compresbrain;
			diffcoef=diffskull;
		}else if (grayscale == 200){
			youngsfmodulus = stiffwm; /* initialize tumor with white matter *if* initial tumor given with label 200 and spherical mask used */
			pr=compresbrain;
			diffcoef=diffwm;		
		}

		mu[ie]     = youngsfmodulus/(2.0*(1+pr));
		lambda[ie] = pr * youngsfmodulus/((1+pr) *(1-2.0*pr));
		diff[ie]   = diffcoef;

		counti++;
		if (counti==xmf){ countj++; counti=1; }
		if (countj==ymf){ counti=1; countj=1; countk++; }
	}

	fclose(fb);

	//For ParaView visualization of the initial material properties
	// fp = fopen("materialmu.fld", "wb");
	// printf("negp is %d\n", negp);
	// fwrite(mu, sizeof(double), negp, fp);
	// fclose(fp); 
	// fp = fopen("materiallambda.fld", "wb");
	// fwrite(lambda, sizeof(double), negp, fp);
	// fclose(fp); 
}
