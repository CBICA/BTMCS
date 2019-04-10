/*!
 * \file AnJacobian.c
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See COPYING file or https://www.rad.upenn.edu/sbia/software/license.html.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include "global.h"
#include "function.h"

/* Version 1.1 of STSolver by F. Abraham and G. Biros 
August 2005*/

#undef __FUNCT__
#define __FUNCT__ "LagComputeJacobian"
int AnComputeJacobian(stsDMMG dmmg,Mat jac)
{
  DA           da = (DA)dmmg->dm;

  int    ie1,ierr,i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
  int count, counti, countj, countk, ie,inl,jnl,knl,lnl,isd;
  int iq, index1,index2, jacobian_flag, *dirichlet_flag;

  PetscReal *v,HxHydHz,HyHzdHx,HxHzdHy;
  PetscReal shyshy, shxshx, shzshz, sh0i,sh0j,tmp1;
  PetscReal shyi,shxi,shzi,shdshd,sh0sh0, shyj,shxj,shzj;
  PetscReal shxshy,shxshz, shyshx,shyshz, shzshy,shzshx,eff0;
  PetscReal facHxsq, facHysq, facHzsq, facHxHy, facHxHz,facHyHz;
  PetscReal *sq,**sh,detinv, *lambda,*mu, lamnode, munode, lampmunode;
  PetscReal stiffness_mat_fac, facx, facy, facz;
                                                                                                                             
  MatStencil   *row, *col;

  ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0);CHKERRQ(ierr);  
  ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  ne = (my-1) * (mx-1);
  if (nsd ==3) ne = ne * (mz-1);
  nnc = (xm) * (ym) ;
  Hx = Lx/ (PetscReal)(mx-1); 
  Hy = Ly/ (PetscReal)(my-1); 

  if (nsd == 3) {
	nnc = nnc * (zm);
  	Hz = Lz/ (PetscReal)(mz-1);
  }
  if ( nsd ==2) Hz = Hx;
 
  facHxsq = Hy*Hz/ (16*Hx);
  facHysq = Hx*Hz/ (16*Hy);
  facHzsq = Hy*Hx/ (16*Hz);

  facHxHy = Hz/ 16;
  facHxHz = Hy/ 16;
  facHyHz = Hx/ 16;



  PetscMalloc(nen*sizeof(PetscReal),&sh);
  PetscMalloc(ne*sizeof(PetscReal),&lambda);
  PetscMalloc(ne*sizeof(PetscReal),&mu);
  ierr=PetscMalloc(nen*ndf*sizeof(int),&dirichlet_flag);CHKERRQ(ierr);
  for (j=0;j<nen;j++){
          PetscMalloc(nsd*sizeof(PetscReal),&sh[j]);
  }

  PetscMalloc(blksize*blksize*sizeof(PetscReal),&v);
  PetscMalloc(blksize*sizeof(MatStencil),&row);
  PetscMalloc(blksize*sizeof(MatStencil),&col);

  analnode(sh);
  createien(xm, ym);
  hacks();
  if (inhomogenity_constant == 1) inhomogeneous_material(lambda, mu, ne,xm,ym, Hx, Hy, Hz);



   
  jacobian_flag = 1;
  counti=1; countj=1; countk=1;
  i = xs; j = ys; k = zs;
  for (ie=0;ie<ne;ie++){
	ie1 = ie *nen;
	i = counti-1; j = countj-1; k = countk-1;
        rowcolindex(i,j,k,row,col,jacobian_flag);
	for (count=0;count<blksize*blksize;count++)v[count]=0;
	dirichlet_identifier(dirichlet_flag, row,mx,my,mz);
	count = 0;
        lamnode = lambda[ie];
        munode  = mu[ie];
        
        //printf("An Jacobian lamnode %g\n", lamnode);

       	for (inl=0;inl < nen;inl++){

		 knl = inl * ndf;
                 shxi = sh[inl][xsd];
                 shyi = sh[inl][ysd];
                 if (ndf ==3) shzi = sh[inl][zsd];

         	 for (jnl=0;jnl<nen;jnl++){
			  lnl = jnl * ndf;
			  shxj = sh[jnl][xsd];
			  shyj = sh[jnl][ysd];

			  facx = (1 + shxi * shxj / 3);
			  facy = (1 + shyi * shyj / 3);

			  shxshx = facHxsq * facy * shxi * shxj;
			  shyshy = facHysq * facx * shyi * shyj;
			  shxshy = facHxHy * shxi * shyj;
			  shyshx = facHxHy * shyi * shxj;


			  if (ndf ==3){

			       shzj   = sh[jnl][zsd];
			       facz   = (1 + shzi * shzj/3); 

			       shxshy = shxshy * facz;
			       shyshx = shyshx * facz;
			       shxshx = shxshx * facz;
			       shyshy = shyshy * facz;

			       shxshz = facHxHz * shxi * shzj * facy;
			       shzshx = facHxHz * shzi * shxj * facy;
			       shyshz = facHyHz * shyi * shzj * facx;
			       shzshy = facHyHz * shzi * shyj * facx;
			       shzshz = facHzsq * facx * facy * shzi * shzj;

			  } else shzshz = 0.0;

			  shdshd = shxshx + shyshy + shzshz;
			  shdshd = munode * shdshd;

			  index1 = (inl * blksize + jnl)* ndf;

			  if ( dynamic ==1) stiffness_mat_fac = 0.5 * positiongamma * dt * dt;
			  else stiffness_mat_fac = 1.0;
		          sh0sh0 = 0.0;
	
			  if (dirichlet_flag[knl] == 1){
			       if (dirichlet_flag[lnl] == 1){
				    v[index1]   += (shdshd + shxshx * (lamnode + munode)) * stiffness_mat_fac;
				    v[index1]   +=  sh0sh0;
			       }
			       if (dirichlet_flag[lnl+1] == 1){
				    v[index1+1] += (shyshx * munode + shxshy * lamnode)   * stiffness_mat_fac;
			       }
			       if (dirichlet_flag[lnl+2] == 1 && ndf ==3 ){
				    v[index1+2] += (shzshx * munode + shxshz * lamnode)   * stiffness_mat_fac;
			       }
			  } else if (inl==jnl){
                                v[index1] = 1.0/hack1e[ie1+inl];
                                if (mgnlevels >1) v[index1] = youngs_global/hack1e[ie1+inl];
                                v[index1+1] = 0.0;
                                if ( nsd==3) v[index1+2] = 0.0;
                          }

			    
			  index2 = index1 + blksize;

			  if (dirichlet_flag[knl+1] == 1){
			       if (dirichlet_flag[lnl] == 1){
				    v[index2]   += (shxshy * munode + shyshx * lamnode)   * stiffness_mat_fac;;
			       }
			       if (dirichlet_flag[lnl+1] == 1){
				    v[index2+1] += (shdshd + shyshy * (lamnode + munode)) * stiffness_mat_fac;
				    v[index2+1] += sh0sh0;
			       }
			       if (dirichlet_flag[lnl+2] == 1 && ndf ==3){
				    v[index2+2] += (shzshy * munode + shyshz * lamnode)   * stiffness_mat_fac;
				}
			  }else if (inl==jnl){
                                v[index2]   = 0.0;
                                v[index2+1] = 1.0/hack1e[ie1+inl];
                                if (mgnlevels >1) v[index2+1] = youngs_global/hack1e[ie1+inl];
                                if ( nsd==3) v[index2+2] = 0.0;
			  }

			  if (nsd ==3){
				index2 = index2 + blksize;
			    	if (dirichlet_flag[knl+2] == 1){
			       	   if (dirichlet_flag[lnl] == 1){
					    v[index2]   += (shxshz * munode + shzshx * lamnode)   * stiffness_mat_fac;;
				    }
			       	    if (dirichlet_flag[lnl+1] == 1){
					    v[index2+1] += (shyshz * munode + shzshy * lamnode)   * stiffness_mat_fac;;
				    }
			       	    if (dirichlet_flag[lnl+2] == 1){
					    v[index2+2] += (shdshd + shzshz * (lamnode + munode)) * stiffness_mat_fac;;
				    	    v[index2+2] += sh0sh0;
				    }
				 } else if (inl==jnl){
					v[index2]   = 0.0;
					v[index2+1] = 0.0;
                                        v[index2+2] = 1.0/hack1e[ie1+inl];
                                        if (mgnlevels >1) v[index2+2] = youngs_global/hack1e[ie1+inl];

				 }
			  }			   		
		 }
	}
	ierr = MatSetValuesStencil(jac,blksize,row,blksize,col,v,ADD_VALUES);CHKERRQ(ierr);
	counti++;
	if (counti==xm){ countj++; counti=1; } 
	if (countj==ym){ counti=1; countj=1; countk++; }
  }
		
  ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatSetBlockSize(jac, nsd);CHKERRQ(ierr);

  ierr = PetscFree(sh);CHKERRQ(ierr);
  ierr = PetscFree(lambda);CHKERRQ(ierr);
  ierr = PetscFree(mu);CHKERRQ(ierr);
  ierr = PetscFree(dirichlet_flag);CHKERRQ(ierr);


  return 0;
}
