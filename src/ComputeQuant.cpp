///////////////////////////////////////////////////////////////////////////////////////
// ComputeQuant.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"


// Computes the force term in the elasticity equation - forward problem
// Input: Vec c (tumor cell normalized density),paramforce=(p1,p2,s);
// Output: *nodal* force

#undef _FUNCT_
#define _FUNCT_ "ComputeForceForward"
void ComputeForceForward(Vec c, double *paramforce, Vec force)
{
	//stencil (box) created in Initialize

	int countx, county, countz;
	int ic,icxb,icxf,icyl,icyr,iczd,iczu,knl;
	int i,j,k;
	PetscScalar *carray;
	double qx,qy,qz,vforce,valforcex,valforcey,valforcez;
	double zero;

	zero=0.0;

	// (fine) grid information
	qx=1.0/(2.0*Hxf); //grid spacing
	qy=1.0/(2.0*Hyf);
	qz=1.0/(2.0*Hzf);
	countx=mxf; //number of *nodes*
	county=myf;
	countz=mzf;

	//force node-wise

	VecGetArray(c,&carray);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	VecSet(&zero,force); //Initialize output vector force
#else
	VecSet(force, zero); //Initialize output vector force
#endif

	for (k=1;k<countz-1;k++) {
		for (j=1;j<county-1;j++) {
			for (i=1;i<countx-1;i++) { 

				ic=i+j*countx+k*countx*county; //current node   
				knl=nsd*ic;		

				//printf("ic carray p1 p2 s %d %g %g %g %g\n",ic,carray[ic],paramforce[0],paramforce[1],paramforce[2]);	  

				// the 6 grid neighbors in x,y,and z for the spatial stencil         	
				icxb=ic-1;
				icxf=ic+1;

				icyl=ic-countx;
				icyr=ic+county;

				iczd=ic-countx*county;
				iczu=ic+countx*county;


				if (paramforce[1]>zero)
				{vforce=-paramforce[0]*exp(-paramforce[1]/(pow(carray[ic],paramforce[2])))*exp(-paramforce[1]/(pow(2.0-carray[ic],paramforce[2])));}
				else
				{vforce=-paramforce[0];}


				valforcex=vforce*qx*(carray[icxf]-carray[icxb]);	
				valforcey=vforce*qy*(carray[icyr]-carray[icyl]);
				valforcez=vforce*qz*(carray[iczu]-carray[iczd]);	

				VecSetValue(force,knl,valforcex,INSERT_VALUES);
				VecSetValue(force,knl+1,valforcey,INSERT_VALUES);
				VecSetValue(force,knl+2,valforcez,INSERT_VALUES);

				//if ((fabs(valforcex)>100.0) || (fabs(valforcey)>100.0) || (fabs(valforcez)>100.0))
				//{printf("ic carray[ic] vforce fx fy fz %d %g %g %g %g %g\n",ic,carray[ic],vforce,valforcex,valforcey,valforcez);}

				//if (ic==243134)
				//{printf("k j i c cf cb cl cr cd cu %d %d %d %g %g %g %g %g %g %g\n",
				//k,j,i,carray[ic],carray[icxf],carray[icxb],carray[icyl],carray[icyr],carray[iczd],carray[iczu]);}
			}
		}
	}

	VecAssemblyBegin(force);
	VecAssemblyEnd(force); 

	VecRestoreArray(c,&carray);

}// end function ComputeForceForward


// Computes the velocity *node-wise*, given the nodal displacement
// Input: Vec dispold, dispnew (nodal displacement at times t and t+dt), time step dt;
// Output: *nodal* velocity Vec v

#undef _FUNCT_
#define _FUNCT_ "ComputeVelocity"
void ComputeVelocity(Vec dispold, Vec dispnew, double dt, Vec v){

	//stencil (box) created in Initialize

	int ic,knl;
	int countx, county,countz, i,j,k;
	PetscScalar *darrayold, *darraynew;
	double qt,valvx,valvy,valvz,zero;

	zero=0.0;

	// (fine) grid information
	qt=1.0/dt;
	countx=mxf; //number of *nodes*
	county=myf;
	countz=mzf;

	//displacement node-wise

	VecGetArray(dispold,&darrayold);
	VecGetArray(dispnew,&darraynew);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	VecSet(&zero,v); //Initialize output vector v
#else
	VecSet(v, zero); //Initialize output vector v
#endif

	for (k=1;k<countz-1;k++) {
		for (j=1;j<county-1;j++) {
			for (i=1;i<countx-1;i++) { 

				ic=i+j*countx+k*countx*county; //current node   
				knl=nsd*ic;	

				valvx=qt*(darraynew[knl]-darrayold[knl]);
				valvy=qt*(darraynew[knl+1]-darrayold[knl+1]);
				valvz=qt*(darraynew[knl+2]-darrayold[knl+2]);

				VecSetValue(v,knl,valvx,INSERT_VALUES);	
				VecSetValue(v,knl+1,valvy,INSERT_VALUES);	
				VecSetValue(v,knl+2,valvz,INSERT_VALUES);		


			}
		}
	}

	VecAssemblyBegin(v);
	VecAssemblyEnd(v); 

	VecRestoreArray(dispold,&darrayold);
	VecRestoreArray(dispnew,&darraynew);

} //end function ComputeVelocity


// Computes the 3D spatial integral in the objective functional based on tumor density (L2 - least squares in *space* only!)
// Input: Vec c (tumor density yielded by the model; nodal); Vec cobj (the target tumor density to be matched);
// Output: *spatial* integral value spaceintobj

#undef _FUNCT_
#define _FUNCT_ "ComputeObjectiveFuncDensitySpace"
double ComputeObjectiveFuncDensitySpace(Vec c, Vec cobj)
{
	//stencil (box) created in Initialize
	int ic;
	int countx,county,countz, i,j,k;
	PetscScalar *carray, *cobjarray;
	double qx,qy,qz,zero;
	double spaceintobj, spaceintobjx, spaceintobjxy, trapfacx, trapfacy, trapfacz;

	zero=0.0;

	// (fine) grid information
	qx=Hxf; //grid spacing
	qy=Hyf;
	qz=Hzf;
	countx=mxf; //number of *nodes*
	county=myf;
	countz=mzf;

	VecGetArray(c,&carray);
	VecGetArray(cobj,&cobjarray);

	spaceintobj=zero;

	for (k=0;k<countz;k++) {

		spaceintobjxy=zero;

		for (j=0;j<county;j++) {

			spaceintobjx=zero;

			for (i=0;i<countx;i++) { 

				ic=i+j*countx+k*countx*county; //current node  

				if (i==0 || i==countx-1) trapfacx=1.0;
				else trapfacx=2.0;

				spaceintobjx=spaceintobjx+trapfacx*(carray[ic]-cobjarray[ic])*(carray[ic]-cobjarray[ic]);


			}; //end i loop

			if (j==0 || j==county-1) trapfacy=1.0;
			else trapfacy=2.0;

			spaceintobjxy=spaceintobjxy+trapfacy*spaceintobjx;

		}; // end j loop

		if (k==0 || k==countz-1) trapfacz=1.0;
		else trapfacz=2.0;

		spaceintobj=spaceintobj+trapfacz*spaceintobjxy;

	} // end k loop

	spaceintobj=1.0/2.0*spaceintobj*Hxf/2.0*Hyf/2.0*Hzf/2.0;

	VecRestoreArray(c,&carray);
	VecRestoreArray(cobj,&cobjarray);

	return(spaceintobj);
} //end function ComputeObjectiveFuncDensitySpace


// Given a  *scalar* quantity (e.g. phase-field variable phi) *element-wise*, InterpolateVecElementToNodal returns the quantity *node-wise*
// at interior nodes; boundary nodes (here the box's 6 faces) must be accounted for separately, based on additional information;
// here, they will be set to zero

#undef __FUNCT__
#define __FUNCT__ "InterpolateElementToNodal"
void InterpolateElementToNodal(Vec velem, Vec vnod)
{
	//stencil (box) created in Initialize
	int countx, county, countz;
	int ie, ic, ie1, ie2, ie3, ie4, ie5, ie6, ie7, ie8, i,j,k;
	PetscScalar  *vecelem;
	double zero,vecn;
	int size;

	zero=0.0;

	VecGetSize(velem,&size);
	printf("velem size %d\n",size);

	VecGetArray(velem,&vecelem); // node-wise 
	//VecGetSize(vecelem,&size);
	//printf("velem size %d\n",size);

	// (fine) grid information  
	countx=mxf-1; //number of *elements* (fine grid)
	county=myf-1;
	countz=mzf-1;

	//printf("countx:%d  county:%d  countz:%d\n",countx, county, countz);

	//VecSet(&zero,vnod); //Initialize output vector, node-wise, *iff* desired; here - it will be initialized upon calling the function

	//VecGetSize(vnod,&size);
	//printf("vnod size %d\n",size);

	//loop over grid, *interior* nodes only     
	for (k=1;k<countz;k++) {
		for (j=1;j<county;j++) {
			for (i=1;i<countx;i++) { 

				ic=i+j*(countx+1)+k*(countx+1)*(county+1); //current interior node
				ie=(i-1)+(j-1)*countx+(k-1)*countx*county; //left-most corner of the current element on the staggered grid with center ic

				//printf("ie %d\n",ie);

				ie1=ie;
				ie2=ie1+1;
				ie3=ie+countx;
				ie4=ie3+1;

				ie5=ie+countx*county;
				ie6=ie5+1;
				ie7=ie5+countx;
				ie8=ie7+1;

				//if( isnan(vecelem[ie1]) || isnan(vecelem[ie2]) || isnan(vecelem[ie3]) || isnan(vecelem[ie4])
				//    || isnan(vecelem[ie5]) || isnan(vecelem[ie6]) || isnan(vecelem[ie7]) || isnan(vecelem[ie8]))
				//  printf("nan\n");
				//if( size <= ie1 || size <= ie2 || size <= ie3 || size <= ie4
				//    || size <= ie5 || size <= ie6 || size <= ie7 || size <= ie8 )
				//		  printf("ooo\n");
				vecn=1.0/8.0*(vecelem[ie1]+vecelem[ie2]+vecelem[ie3]+vecelem[ie4]+vecelem[ie5]+vecelem[ie6]+vecelem[ie7]+vecelem[ie8]);   
				//printf("%f %f %f %f %f %f %f %f\n",vecelem[ie1],vecelem[ie2],vecelem[ie3],vecelem[ie4],vecelem[ie5],vecelem[ie6],vecelem[ie7],vecelem[ie8]);

				VecSetValue(vnod,ic,vecn,INSERT_VALUES); 			   

				//if ((fabs(vecxe)>1e-6) || (fabs(vecye)>1e-6) || (fabs(vecze)>1e-6))
				//{printf("ie vecxe vecye vecze %d %g %g %g\n",ie,vecxe,vecye,vecze);}	
			} //end i loop
		} //end j loop
	} //end k loop

	printf("_________________________From inside of Interpolate___________\n");                

	VecAssemblyBegin(vnod);
	VecAssemblyEnd(vnod); 

	VecRestoreArray(velem,&vecelem);

	//printf("size %d\n",size);
	//printf("out of interpolate\n");
}//end function InterpolateElementToNodal


// Given a *scalar* quantity (e.g. density) *node-wise*, InterpolateNodalToElement returns the quantity *element-wise
#undef __FUNCT__
#define __FUNCT__ "InterpolateNodalToElement"
void InterpolateNodalToElement(Vec vnod, Vec velem)
{
	//stencil (box) created in Initialize
	int countx, county, countz;
	int ie, ic, ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, i,j,k;
	//Vec vnoddup;
	PetscScalar  *vecnod;
	double zero,vece;
	int counte,size;

	zero=0.0;

	VecGetSize(vnod,&size);
	//printf("size %d\n",size);

	//VecDuplicate(vnod, &vnoddup);
	//VecCopy(vnod,vnoddup);

	//VecGetSize(vnoddup,&size);
	//printf("size %d\n",size);

	VecGetArray(vnod,&vecnod); // node-wise 

	// (fine) grid information  
	countx=mxf; //number of *nodes* (fine grid)
	county=myf;
	countz=mzf;

	counte=(countx-1)*(county-1)*(countz-1); //number of *elements* (fine grid)

	//VecCreate(PETSC_COMM_WORLD,&velem);
	//VecSetSizes(velem,PETSC_DECIDE,nsd*counte);
	//VecSetFromOptions(velem);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	VecSet(&zero,velem); //Initialize output vector, element-wise
#else
	VecSet(velem, zero); //Initialize output vector, element-wise
#endif

	VecGetSize(velem,&size);
	//printf("size %d\n",size);

	//loop over staggered grid (element center); diffusivity defined *element-wise*       
	for (k=0;k<countz-1;k++) {
		for (j=0;j<county-1;j++) {
			for (i=0;i<countx-1;i++) { 

				ie=i+j*(countx-1)+k*(countx-1)*(county-1); //current element
				ic=i+j*countx+k*countx*county; //left-most corner of the current element with center ie

				//printf("ie %d\n",ie);

				ic1=ic;
				ic2=ic1+1;
				ic3=ic+countx;
				ic4=ic3+1;

				ic5=ic+countx*county;
				ic6=ic5+1;
				ic7=ic5+countx;
				ic8=ic7+1;

				vece=1.0/8.0*(vecnod[ic1]+vecnod[ic2]+vecnod[ic3]+vecnod[ic4]+vecnod[ic5]+vecnod[ic6]+vecnod[ic7]+vecnod[ic8]);   

				VecSetValue(velem,ie,vece,INSERT_VALUES); 			   

				//if ((fabs(vecxe)>1e-6) || (fabs(vecye)>1e-6) || (fabs(vecze)>1e-6))
				//{printf("ie vecxe vecye vecze %d %g %g %g\n",ie,vecxe,vecye,vecze);}	
			} //end i loop
		} //end j loop
	} //end k loop

	VecAssemblyBegin(velem);
	VecAssemblyEnd(velem); 

	VecRestoreArray(vnod,&vecnod);

	//VecDestroy(vnoddup);

	//printf("size %d\n",size);
	//printf("out of interpolate\n");
}//end function InterpolateNodalToElement
