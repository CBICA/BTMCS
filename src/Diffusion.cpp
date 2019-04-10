///////////////////////////////////////////////////////////////////////////////////////
// Diffusion.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include "function.h"

#include "ConservLaw.h"


// Linear diffusion - variable diffusivity - over time interval dt; fully explicit scheme.
// Input: c0 ('initial' condition at time t), D (diffusivity at time t) (include Vec s <-> source)
// for consistency with the ConservationLaw prototype;
// Output: c (updated quantity at time t+dt)
// This is always performed *node-wise* here

#undef __FUNCT__
#define __FUNCT__ "DiffusionExplicit"
void DiffusionExplicit(struct ConservationLaw *cl, Vec c0, Vec D, Vec s, double dt, double cfln,  Vec c)
{
	//stencil (box) created in Initialize
	int countx, county, countz;
	int ic, ie, ie1, ie2, ie3, ie4, ie5, ie6, ie7, ie8, icxb,icxf,icyl,icyr,iczd,iczu;
	int count, nsteps, i,j,k;
	int size;
	Vec c0dup;
	PetscScalar *cold, *cnew, *Diff, *source;
	double qx,qy,qz,cflthrsh,ratio,zero;
	double Dxb, Dxf, Dyl, Dyr, Dzd, Dzu;

	zero=0.0;

	// (fine) grid information
	qx=1.0/(Hxf*Hxf); //grid spacing
	qy=1.0/(Hyf*Hyf);
	qz=1.0/(Hzf*Hzf);
	countx=mxf; //number of *nodes*
	county=myf;
	countz=mzf;

	VecGetSize(c0,&size); //get size of the block-vector c0
	//printf("size %d\n",size);

	//diffusion node-wise

	VecDuplicate(c0,&c0dup);
	VecCopy(c0,c0dup);

	VecGetArray(c0dup,&cold);
	VecGetArray(D,&Diff); 
	VecGetArray(s,&source); 
	PetscMalloc(size*sizeof(double),&cnew);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	VecSet(&zero,c); //Initialize output vector c
#else
	VecSet(c, zero); //Initialize output vector c
#endif

	//compute the "optimal" number of sub-steps in the interval (t,t+dt) to obey the stability condition
	//loop over staggered grid (element center); diffusivity defined *element-wise* 
	nsteps=1;

	for (k=0;k<countz-1;k++) {
		for (j=0;j<county-1;j++) {
			for (i=0;i<countx-1;i++) { 

				//set cfl with respect to diffusivity 

				ie=i+j*(countx-1)+k*(countx-1)*(county-1); //current element   			  

				cflthrsh=dt*fabs(Diff[ie])*(qx+qy+qz);

				//printf("ie Diff[ie] cflthrsh %d %g %g\n",ie,Diff[ie],cflthrsh);

				if (cflthrsh>cfln) {nsteps=(int)max(nsteps,rint(cflthrsh/cfln+0.5));}

			} //end loop i 	  
		}//end loop j       
	}//end loop k	

	//printf("nsteps %d\n",nsteps);
	ratio=dt/nsteps;

	for (count=1;count<=nsteps;count++){
		//loop over the original grid: diffusion *node-wise*; leave the box (fictitious domain) boundaries outside loop here
		for (k=1;k<countz-1;k++) {
			for (j=1;j<county-1;j++) {
				for (i=1;i<countx-1;i++) { 

					ic=i+j*countx+k*countx*county; //current node   			  

					// the 6 grid neighbors in x,y,and z for the spatial stencil         	
					icxb=ic-1;
					icxf=ic+1;

					icyl=ic-countx;
					icyr=ic+county;

					iczd=ic-countx*county;
					iczu=ic+countx*county;

					//neighboring elements needed to interpolate diffusivities 
					// note: can be re-written more efficient in terms of ic
					ie1=(i-1)+(j-1)*(countx-1)+(k-1)*(countx-1)*(county-1);
					ie2=i+(j-1)*(countx-1)+(k-1)*(countx-1)*(county-1);
					ie3=(i-1)+j*(countx-1)+(k-1)*(countx-1)*(county-1);
					ie4=i+j*(countx-1)+(k-1)*(countx-1)*(county-1);

					ie5=(i-1)+(j-1)*(countx-1)+k*(countx-1)*(county-1);
					ie6=i+(j-1)*(countx-1)+k*(countx-1)*(county-1);
					ie7=(i-1)+j*(countx-1)+k*(countx-1)*(county-1);
					ie8=i+j*(countx-1)+k*(countx-1)*(county-1);


					// interpolated diffusivities
					Dxb=1.0/4.0*(Diff[ie1]+Diff[ie3]+Diff[ie5]+Diff[ie7]);	
					Dxf=1.0/4.0*(Diff[ie2]+Diff[ie4]+Diff[ie6]+Diff[ie8]);
					Dyl=1.0/4.0*(Diff[ie1]+Diff[ie2]+Diff[ie5]+Diff[ie6]);
					Dyr=1.0/4.0*(Diff[ie3]+Diff[ie4]+Diff[ie7]+Diff[ie8]);
					Dzd=1.0/4.0*(Diff[ie1]+Diff[ie2]+Diff[ie3]+Diff[ie4]);
					Dzu=1.0/4.0*(Diff[ie5]+Diff[ie6]+Diff[ie7]+Diff[ie8]);


					//update c (fully explicit scheme)        		
					cnew[ic]=cold[ic]+ratio*\
						(qx*(Dxf*(cold[icxf]-cold[ic])-Dxb*(cold[ic]-cold[icxb]))+\
						qy*(Dyr*(cold[icyr]-cold[ic])-Dyl*(cold[ic]-cold[icyl]))+\
						qz*(Dzu*(cold[iczu]-cold[ic])-Dzd*(cold[ic]-cold[iczd]))+\
						source[ic]); 




				} //end loop i 	  
			}//end loop j       
		}//end loop k	


		//homogeneous Neumann BC for c on the box boundary

		k=0;
		for (j=0;j<county;j++) {
			for (i=0;i<countx;i++) {
				ic=i+j*countx+k*countx*county;
				cnew[ic]=cnew[ic+countx*county];	        
			} 
		}	

		k=countz-1;
		for (j=0;j<county;j++) {
			for (i=0;i<countx;i++) {
				ic=i+j*countx+k*countx*county;
				cnew[ic]=cnew[ic-countx*county];	        
			} 
		}

		j=0;
		for (k=0;k<countz;k++) {
			for (i=0;i<countx;i++) {
				ic=i+j*countx+k*countx*county;
				cnew[ic]=cnew[ic+countx];	        
			} 
		}

		j=county-1;
		for (k=0;k<countz;k++) {
			for (i=0;i<countx;i++) {
				ic=i+j*countx+k*countx*county;
				cnew[ic]=cnew[ic-countx];	        
			} 
		}

		i=0;
		for (k=0;k<countz;k++) {
			for (j=0;j<county;j++) {
				ic=i+j*countx+k*countx*county;
				cnew[ic]=cnew[ic+1];	        
			} 
		}

		i=countx-1;
		for (k=0;k<countz;k++) {
			for (j=0;j<county;j++) {
				ic=i+j*countx+k*countx*county;
				cnew[ic]=cnew[ic-1];	        
			} 
		}				  

		for (k=0;k<countz;k++) {
			for (j=0;j<county;j++) {
				for (i=0;i<countx;i++) { 
					ic=i+j*countx+k*countx*county;
					cold[ic]=cnew[ic];
				}
			}
		}
	} //end count loop

	//construct the output vector Vec c
	for (k=0;k<countz;k++) {
		for (j=0;j<county;j++) {
			for (i=0;i<countx;i++) { 
				ic=i+j*countx+k*countx*county;
				VecSetValue(c,ic,cnew[ic],INSERT_VALUES);
			}
		}
	}

	VecAssemblyBegin(c);
	VecAssemblyEnd(c); 

	VecRestoreArray(c0dup,&cold);
	VecRestoreArray(D,&Diff);
	VecRestoreArray(s,&source);
	PetscFree(cnew);
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	VecDestroy(c0dup);
#else
	VecDestroy(&c0dup);
#endif
} //end function AdvectionDiv
