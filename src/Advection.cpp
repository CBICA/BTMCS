///////////////////////////////////////////////////////////////////////////////////////
// Advection.c
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


// Given a *vectorial* quantity (e.g. velocity) *node-wise*, InterpolateVecNodalToElement returns the quantity *element-wise

#undef __FUNCT__
#define __FUNCT__ "InterpolateVecNodalToElement"
void InterpolateVecNodalToElement(Vec vnod, Vec velem)
{
	//stencil (box) created in Initialize
	int countx, county, countz;
	int ie, ic, ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8, knl, knl1, knl2, knl3, knl4, knl5, knl6, knl7, knl8, i,j,k;
	//Vec vnoddup;
	PetscScalar  *vecnod;
	double zero,vecxe,vecye,vecze;
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

				knl1 = ic1 * nsd;
				knl2 = ic2 * nsd;
				knl3 = ic3 * nsd;
				knl4 = ic4 * nsd;
				knl5 = ic5 * nsd;
				knl6 = ic6 * nsd;
				knl7 = ic7 * nsd;
				knl8 = ic8 * nsd;

				knl=ie*nsd;

				vecxe=1.0/8.0*(vecnod[knl1]+vecnod[knl2]+vecnod[knl3]+vecnod[knl4]+vecnod[knl5]+vecnod[knl6]+vecnod[knl7]+vecnod[knl8]);   
				vecye=1.0/8.0*(vecnod[knl1+1]+vecnod[knl2+1]+vecnod[knl3+1]+vecnod[knl4+1]+vecnod[knl5+1]+vecnod[knl6+1]+vecnod[knl7+1]+vecnod[knl8+1]);
				vecze=1.0/8.0*(vecnod[knl1+2]+vecnod[knl2+2]+vecnod[knl3+2]+vecnod[knl4+2]+vecnod[knl5+2]+vecnod[knl6+2]+vecnod[knl7+2]+vecnod[knl8+2]);

				VecSetValue(velem,knl,vecxe,INSERT_VALUES); 			   

				VecSetValue(velem,knl+1,vecye,INSERT_VALUES); 

				VecSetValue(velem,knl+2,vecze,INSERT_VALUES); 

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
}//end function Interpolate


// Linear advection - "color equation" - over time interval dt, with CFL number cfln;
// Input: m0 ('initial' condition at time t),v (and eventually s <-> right hand side, if source term present (e.g. adjoint eqns.)
// Output: m (updated quantity at time t+dt)
// Note: in the implementation below, we assume the velocity v given at the *nodes*
// consistent with the elasticity solution u; interpolating for the velocity
// per element is done in a separate function whenever necessary

#undef __FUNCT__
#define __FUNCT__ "AdvectionGrad"
void AdvectionGrad(struct ConservationLaw *cl, Vec m0, Vec v, Vec s, double dt, double cfln,  Vec m)
{
	//stencil (box) created in Initialize
	int countx, county, countz;
	int is, count, countspace, nsteps,i,j,k,l, knl;
	int isl, isxbl, isxfl, isyll, isyrl, iszdl, iszul, size, nblock;
	Vec m0dup;
	Vec velem;
	//Vec velemdup, vdup, sdup;
	PetscScalar *mold, *mnew, *vel, *source;
	double qx,qy,qz,cflthrsh,velx,vely,velz,zero,ratio;

	zero=0.0;

	// (fine) grid information
	qx=1.0/Hxf; //grid spacing
	qy=1.0/Hyf;
	qz=1.0/Hzf;

	VecGetSize(m0,&size);
	//VecView(m0,PETSC_VIEWER_STDOUT_WORLD);

	if (fmod((double)size, (double)ne_fine) == 0) {
		nblock=size/ne_fine; //if advection element-wise
		countspace=ne_fine;
		VecCreate(PETSC_COMM_WORLD,&velem);
		VecSetSizes(velem,PETSC_DECIDE,nsd*countspace);
		VecSetFromOptions(velem);
		InterpolateVecNodalToElement(v,velem);
		//VecGetSize(velem,&size);
		//printf("size %d\n",size);
		//VecDuplicate(velem,&velemdup);
		//VecCopy(velem,velemdup);
		VecGetArray(velem,&vel); //advection velocity: element-wise, after Interpolation
		countx=mxf-1;
		county=myf-1;
		countz=mzf-1;	
	}
	else if (fmod((double)size, (double)nnc_fine) == 0) {
		nblock=size/nnc_fine;  //if advection node-wise
		countspace=nnc_fine;
		//VecDuplicate(v,&vdup);
		//VecCopy(v,vdup);
		VecGetArray(v,&vel);
		countx=mxf;
		county=myf;
		countz=mzf;
	}

	//printf("before duplicate\n");
	VecDuplicate(m0,&m0dup);
	//printf("after duplicate\n");
	//VecDuplicate(s,&sdup);  

	VecCopy(m0,m0dup);
	//VecCopy(s,sdup);

	VecGetArray(m0dup,&mold); 
	VecGetArray(s,&source); 
	PetscMalloc(size*sizeof(double),&mnew);

	//Initialize output vector m
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	VecSet(&zero,m); 
#else
	VecSet(m, zero);
#endif

	//compute the "optimal" number of sub-steps in the interval (t,t+dt) to obey the CFL condition: nsteps

	nsteps=1;

	if (cfln > 0) {
		for (k=0;k<countz;k++) {
			for (j=0;j<county;j++) {
				for (i=0;i<countx;i++) {

					is=i+j*countx+k*countx*county;   //current spatial location (element or node) 	
					knl = is * nsd;

					//set cfl with respect to velocity 			   
					cflthrsh=dt*(qx*fabs(vel[knl])+qy*fabs(vel[knl+1])+qz*fabs(vel[knl+2]));

					if (cflthrsh>cfln) {nsteps=(int)max(nsteps,rint(cflthrsh/cfln+0.5));}

				} //end loop i 	  
			}//end loop j       
		}//end loop k
	}

	ratio=dt/nsteps;

	//printf("nsteps nblock dt ratio qx qy qz %d %d %g %g %g %g %g\n",nsteps,nblock,dt,ratio,qx,qy,qz);

	for (count=1;count<=nsteps;count++){

		//loop over the spatial grid; leave the box (fictitious domain) boundaries outside loop here

		for (k=1;k<countz-1;k++) {
			for (j=1;j<county-1;j++) {
				for (i=1;i<countx-1;i++) {

					is=i+j*countx+k*countx*county;   //current spatial location (element or node) 	
					knl =is*nsd;

					velx=vel[knl];	   		   			   
					vely=vel[knl+1];
					velz=vel[knl+2];

					for (l=0;l<nblock;l++) {

						isl=is+l*countspace;

						isxbl=isl-1;
						isxfl=isl+1;

						isyll=isl-countx;
						isyrl=isl+countx;

						iszdl=isl-countx*county;
						iszul=isl+countx*county;

						mnew[isl]=mold[isl]-ratio*\
							(qx*(max(velx,zero)*(mold[isl]-mold[isxbl])+min(velx,zero)*(mold[isxfl]-mold[isl]))+\
							qy*(max(vely,zero)*(mold[isl]-mold[isyll])+min(vely,zero)*(mold[isyrl]-mold[isl]))+\
							qz*(max(velz,zero)*(mold[isl]-mold[iszdl])+min(velz,zero)*(mold[iszul]-mold[isl]))+\
							source[isl]); 

						//if ((fabs(velx)>0.001) || (fabs(vely)>0.001) || (fabs(velz)>0.001))
						//{printf("is l velx vely velz mold[isl] mnew[isl]	%d %d %g %g %g %g %g\n",is,l,velx,vely,velz,mold[isl],mnew[isl]);}	


					} //end loop l


				} //end loop i 	  
			}//end loop j       
		}//end loop k	   

		// artificial (numerical) BC on the box boundary for the advected quantity; try fixed Dirichlet here
		// if problems, implement/use reflective boundary condition instead

		for (k=1;k<countz-1;k++) {
			for (j=1;j<county-1;j++) {
				for (i=1;i<countx-1;i++) {

					is=i+j*countx+k*countx*county;

					for (l=0;l<nblock;l++) {

						isl=is+l*countspace;
						mold[isl]=mnew[isl];

					} 
				}  	  
			}      
		}
	} //end loop count	

	//printf("Out of main loop AdvectionGrad nblock %d\n",nblock);

	// construct output vector

	for (k=0;k<countz;k++) {
		for (j=0;j<county;j++) {
			for (i=0;i<countx;i++) {
				is=i+j*countx+k*countx*county;
				for (l=0;l<nblock;l++) {
					isl=is+l*countspace;
					//printf("is l mold[isl] %d %d %g\n",is,l,mold[isl]);
					VecSetValue(m,isl,mold[isl],INSERT_VALUES);
				}
			} 	  
		}      
	}

	VecAssemblyBegin(m);
	VecAssemblyEnd(m); 

	VecRestoreArray(m0dup,&mold);
	if (fmod((double)size, (double)ne_fine) == 0) {
		VecRestoreArray(velem,&vel); 
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
		VecDestroy(velem);
#else
		VecDestroy(&velem);
#endif
	}      
	else if (fmod((double)size, (double)nnc_fine) == 0) {
		VecRestoreArray(v,&vel);
	}    
	VecRestoreArray(s,&source);
	PetscFree(mnew);

#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	VecDestroy(m0dup);
	//VecDestroy(sdup);
#else
	VecDestroy(&m0dup);
	//VecDestroy(&sdup);
#endif
} //end function AdvectionGrad


// Linear advection - conservative mass equation - over time interval dt, with CFL number cfln;
// Input: m0 ('initial' condition at time t),v (and eventually s <-> right hand side, if source term present (e.g. adjoint eqns.)
// Output: m (updated quantity at time t+dt)
// Note: in the implementation below, we assume the velocity v given at the *nodes*
// consistent with the elasticity solution u

#undef __FUNCT__
#define __FUNCT__ "AdvectionDiv"
void AdvectionDiv(struct ConservationLaw *cl, Vec m0, Vec v, Vec s, double dt, double cfln,  Vec m)
{
	//stencil (box) created in Initialize
	int countx, county, countz; 
	int count, countspace, nsteps,i,j,k,l, knl;
	int is, isxb, isxf, isyl, isyr, iszd, iszu;
	int isl, isxbl, isxfl, isyll, isyrl, iszdl, iszul, size, nblock;
	PetscScalar *mold, *mnew, *vel, *source, *phi;
	Vec m0dup,velem ;
	//Vec velemdup, vdup, sdup;
	double qx,qy,qz,cflthrsh,velx,vely,velz,velxb,velyl,velzd,velxf,velyr,velzu,zero,ratio;

	zero=0.0;
	qx=1.0/Hxf;
	qy=1.0/Hyf;
	qz=1.0/Hzf;

	VecGetSize(m0,&size); //get size of the block-vector m0
	//printf("size %d\n",size);

	if (fmod((double)size, (double)ne_fine) == 0) {
		nblock=size/ne_fine; //if advection element-wise
		countspace=ne_fine;
		VecCreate(PETSC_COMM_WORLD,&velem);
		VecSetSizes(velem,PETSC_DECIDE,nsd*countspace);
		VecSetFromOptions(velem);
		InterpolateVecNodalToElement(v,velem);
		//VecDuplicate(velem,&velemdup);
		//VecCopy(velem,velemdup);
		VecGetArray(velem,&vel); //advection velocity: element-wise, after Interpolation
		countx=mxf-1;
		county=myf-1;	
		countz=mzf-1;
	}
	else if (fmod((double)size, (double)nnc_fine) == 0) {
		nblock=size/nnc_fine;  //if advection node-wise
		countspace=nnc_fine;
		//VecDuplicate(v,&vdup);
		//VecCopy(v,vdup);
		VecGetArray(v,&vel);
		countx=mxf;
		county=myf;
		countz=mzf;
	}

	//incorporate the global phase-field variable gphi here - iff tracking ventricles necessary !!
	VecGetArray(gphi,&phi);

	for (k=0;k<countz;k++) {
		for (j=0;j<county;j++) {
			for (i=0;i<countx;i++) {
				is=i+j*countx+k*countx*county;   //current spatial location (element or node) 	
				knl = is * nsd;

				vel[knl]=vel[knl]*phi[is];
				vel[knl+1]=vel[knl+1]*phi[is];
				vel[knl+2]=vel[knl+2]*phi[is];
			} //end loop i 	  
		}//end loop j       
	}//end loop k	

	VecDuplicate(m0,&m0dup);
	//VecDuplicate(s,&sdup);
	VecCopy(m0,m0dup);
	//VecCopy(s,sdup);

	VecGetArray(m0dup,&mold); 
	VecGetArray(s,&source); 
	PetscMalloc(size*sizeof(double),&mnew);

	//VecCreate(PETSC_COMM_WORLD,m);
	//VecSetSizes(*m,PETSC_DECIDE,size);
	//VecSetFromOptions(*m);

#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	VecSet(&zero,m); //Initialize output vector m
#else
	VecSet(m, zero); //Initialize output vector m
#endif

	//compute the "optimal" number of sub-steps in the interval (t,t+dt) to obey the CFL condition: nsteps

	nsteps=1;
	for (k=0;k<countz;k++) {
		for (j=0;j<county;j++) {
			for (i=0;i<countx;i++) {
				is=i+j*countx+k*countx*county;   //current spatial location (element or node) 	
				knl = is * nsd;

				//set cfl with respect to velocity 			   
				cflthrsh=dt*(qx*fabs(vel[knl])+qy*fabs(vel[knl+1])+qz*fabs(vel[knl+2]));

				if (cflthrsh>cfln) {nsteps=(int)max(nsteps,rint(cflthrsh/cfln+0.5));}

			} //end loop i 	  
		}//end loop j       
	}//end loop k	

	//printf("nsteps %d\n",nsteps);

	ratio=dt/nsteps;

	for (count=1;count<=nsteps;count++){
		//loop over the spatial grid; leave the box (fictitious domain) boundaries outside loop here
		for (k=1;k<countz-1;k++) {
			for (j=1;j<county-1;j++) {
				for (i=1;i<countx-1;i++) {
					is=i+j*countx+k*countx*county;   //current spatial location (element or node) 	
					knl =is*nsd;

					velx=vel[knl];	   		   			   
					vely=vel[knl+1];
					velz=vel[knl+2];

					isxb=is-1;
					isxf=is+1;

					isyl=is-countx;
					isyr=is+county;

					iszd=is-countx*county;
					iszu=is+countx*county;

					velxb=1.0/2.0*(vel[isxb*nsd]+velx);
					velxf=1.0/2.0*(vel[isxf*nsd]+velx);
					velyl=1.0/2.0*(vel[isyl*nsd]+vely);
					velyr=1.0/2.0*(vel[isyr*nsd]+vely);
					velzd=1.0/2.0*(vel[iszd*nsd]+velz);
					velzu=1.0/2.0*(vel[iszu*nsd]+velz);

					for (l=0;l<nblock;l++) {
						isl=is+l*countspace;

						isxbl=isl-1;
						isxfl=isl+1;

						isyll=isl-countx;
						isyrl=isl+countx;

						iszdl=isl-countx*county;
						iszul=isl+countx*county;

						mnew[isl]=mold[isl]-ratio*\
							(qx*((max(velxf,zero)*mold[isl]+min(velxf,zero)*mold[isxfl])-(max(velxb,zero)*mold[isxbl]+min(velxb,zero)*mold[isl]))+\
							qy*((max(velyr,zero)*mold[isl]+min(velyr,zero)*mold[isyrl])-(max(velyl,zero)*mold[isyll]+min(velyl,zero)*mold[isl]))+\
							qz*((max(velzu,zero)*mold[isl]+min(velzu,zero)*mold[iszul])-(max(velzd,zero)*mold[iszdl]+min(velzd,zero)*mold[isl]))+\
							source[isl]); 
					} //end loop l
				} //end loop i 	  
			}//end loop j       
		}//end loop k	   

		// artificial (numerical) BC on the box boundary for the advected quantity; try fixed Dirichlet here
		for (k=1;k<countz-1;k++) {
			for (j=1;j<county-1;j++) {
				for (i=1;i<countx-1;i++) {
					is=i+j*countx+k*countx*county;

					for (l=0;l<nblock;l++) {
						isl=is+l*countspace;
						mold[isl]=mnew[isl];
					}
				}  	  
			}      
		}
	} //end loop count	

	// construct output vector
	for (k=0;k<countz;k++) {
		for (j=0;j<county;j++) {
			for (i=0;i<countx;i++) {
				is=i+j*countx+k*countx*county;
				for (l=0;l<nblock;l++) {
					isl=is+l*countspace;
					VecSetValue(m,isl,mold[isl],INSERT_VALUES);
				}
			} 	  
		}      
	}

	VecAssemblyBegin(m);
	VecAssemblyEnd(m); 

	VecRestoreArray(m0dup,&mold);
	if (fmod((double)size, (double)ne_fine) == 0) {
		VecRestoreArray(velem,&vel); 
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
		VecDestroy(velem);
#else
		VecDestroy(&velem);
#endif
	}
	else if (fmod((double)size, (double)nnc_fine) == 0) {VecRestoreArray(v,&vel);}
	//VecRestoreArray(sdup,&source);
	PetscFree(mnew);

	VecRestoreArray(gphi,&phi);

#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	VecDestroy(m0dup); 
	//VecDestroy(sdup);
#else
	VecDestroy(&m0dup); 
	//VecDestroy(&sdup);
#endif
} //end function AdvectionDiv


// Input: the "block" vector m; nvecs=number of compounding Vecs (e.g. collective material properties m=(lambda,mu,D); nvecs=3)
// Output: array of vectors *msep<->nvecs*Vec

#undef __FUNCT__
#define __FUNCT__ "RetrieveVecs"
void RetrieveVecs(Vec m, Vec *msep)
{
	//stencil (box) created in Initialize
	int countx, county, countz; 
	int countspace,i,j,k,l;
	int is, isl, size, nblock;
	PetscScalar *marray;
	//Vec  mdup;
	double zero, msepval;

	zero=0.0;

	VecGetSize(m,&size); //get size of the block-vector m

	if (fmod((double)size, (double)ne_fine) == 0) {
		nblock=size/ne_fine; //if advection element-wise
		countspace=ne_fine;
		countx=mxf-1;
		county=myf-1;	
		countz=mzf-1;
	}
	else if (fmod((double)size, (double)nnc_fine) == 0) {
		nblock=size/nnc_fine;  //if advection node-wise
		countspace=nnc_fine;
		countx=mxf;
		county=myf;
		countz=mzf;
	}

	//printf("sizem in RetrieveVecs %d %d\n",size,nblock);

	//VecDuplicate(m,&mdup);
	//VecCopy(m,mdup);

	VecGetArray(m,&marray); 

	//Initialize output array of vectors msep
	for (l=0;l<nblock;l++) {	  		
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
		VecSet(&zero,msep[l]); 	
#else
		VecSet(msep[l], zero);
#endif
		//VecGetSize(msep[l],&size);
		//printf("l size %d %d\n",l,size);		
	} 

	//loop over the spatial grid
	for (k=0;k<countz;k++) {
		for (j=0;j<county;j++) {
			for (i=0;i<countx;i++) {
				is=i+j*countx+k*countx*county;   //current spatial location (element or node) 	

				//printf("is %d\n",is);

				//if (is>0) break;

				for (l=0;l<nblock;l++) {
					isl=is+l*countspace;

					msepval=marray[isl]; 

					VecSetValue(msep[l],is,msepval,INSERT_VALUES);
				} //end loop l
			} //end loop i 	  
		}//end loop j       
	}//end loop k

	for (l=0;l<nblock;l++) {   
		VecAssemblyBegin(msep[l]);
		VecAssemblyEnd(msep[l]);
	} 

	VecRestoreArray(m,&marray);	
} //end function RetrieveVecs 
