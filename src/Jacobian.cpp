///////////////////////////////////////////////////////////////////////////////////////
// Jacobian.c
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


extern int Penalized_Neumann_Contribution_mat(stsDMMG, Mat);


#undef __FUNCT__
#define __FUNCT__ "ComputeJacobian"
PetscInt ComputeJacobian(stsDMMG dmmg,Mat jac)
{
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	DA da = (DA)dmmg->dm;
#else
	DM da = dmmg->dm;
#endif
	PetscReal *v;

	PetscReal shyshy, shxshx, shzshz, sh0i,sh0j;
	PetscReal shxshy, shxshz, shyshx, shyshz, shzshy,shzshx,eff0;
	PetscReal shyi, shxi, shzi, shdshd, sh0sh0, shyj,shxj,shzj;
	PetscReal **xq,**sq, ***sh,*det,*wq;
	PetscReal *lambda,*mu, lamnode, munode;
	PetscReal stiffness_mat_fac;

	PetscInt    iq, index1,index2, jacobian_flag, *dirichlet_flag;
	PetscInt    ie1,ierr,i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
	PetscInt    count, counti, countj, countk, ie,inl,JNL,knl,lnl;
	MatStencil   *row, *col;

#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0);CHKERRQ(ierr);  
	ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
#else
	ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);  
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
#endif

	ne = (my-1) * (mx-1);
	if (nsd ==3) ne = ne * (mz-1);
	nnc = (xm) * (ym) ;
	if (nsd == 3) nnc = nnc * (zm);
	Hx = Lx/ (PetscReal)(mx-1); Hy = Ly/ (PetscReal)(my-1); 
	if (nsd==3)Hz = Lz/ (PetscReal)(mz-1);

	shzshz = 0.0; shxshz = 0.0; shzshx = 0.0;shyshz = 0.0; shzshy = 0.0; shzi = 0.0;

	PetscMalloc(nquad*sizeof(PetscReal),&det);
	PetscMalloc(nquad*sizeof(PetscReal),&wq);

	PetscMalloc(nquad*sizeof(PetscReal),&xq);
	PetscMalloc(nquad*sizeof(PetscReal),&sh);
	PetscMalloc(nen*sizeof(PetscReal),&sq);

	ierr=PetscMalloc(ne*sizeof(PetscReal),&lambda);CHKERRQ(ierr);
	ierr=PetscMalloc(ne*sizeof(PetscReal),&mu);CHKERRQ(ierr);
	ierr=PetscMalloc(nen*ndf*sizeof(PetscInt),&dirichlet_flag);CHKERRQ(ierr);

	for (j=0;j<nen;j++) PetscMalloc(nquad*sizeof(PetscReal),&sq[j]);

	for (j=0;j<nquad;j++){
		PetscMalloc(nsd*sizeof(PetscReal),&xq[j]);
		PetscMalloc(nen*sizeof(PetscReal),&sh[j]);
		for (i=0;i<nen;i++) PetscMalloc(nsd*sizeof(PetscReal),&sh[j][i]);
	}

	PetscMalloc(blksize*blksize*sizeof(PetscReal),&v);
	PetscMalloc(blksize*sizeof(MatStencil),&row);
	PetscMalloc(blksize*sizeof(MatStencil),&col);

	createien(xm, ym);
	hacks();
	quad(xq, sq, sh,det,wq,Hx, Hy, Hz);
	/* Piecewise constant interpolation */
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
		lamnode = lambda[ie];
		munode  = mu[ie];

		for (iq=0;iq<nquad;iq++){
			count = 0;
			eff0 = wq[iq] * det[iq];
			for (inl=0;inl < nen;inl++){
				knl = inl * ndf;
				sh0i = sq[inl][iq];
				shxi = sh[iq][inl][xsd] * eff0;
				shyi = sh[iq][inl][ysd] * eff0;
				if (ndf ==3) shzi = sh[iq][inl][zsd] * eff0;
				for (JNL=0;JNL<nen;JNL++){
					lnl = JNL * ndf;
					sh0j = sq[JNL][iq];
					shxj = sh[iq][JNL][xsd];
					shyj = sh[iq][JNL][ysd];

					shxshx = shxi * shxj;
					shyshy = shyi * shyj;
					shxshy = shxi * shyj;
					shyshx = shyi * shxj;
					shdshd = shxshx + shyshy;
					if (dynamic == 1){
						sh0sh0 = sh0i * sh0j;
						sh0sh0 = sh0sh0 * (density + damping_factor * velocityalpha * dt);
					}else{
						sh0sh0 = 0.0;
					}

					if (ndf ==3){
						shzj = sh[iq][JNL][zsd];
						shxshz = shxi * shzj;
						shzshz = shzi * shzj;
						shyshz = shyi * shzj;
						shzshx = shzi * shxj;
						shzshy = shzi * shyj;
						shdshd +=shzshz;
					}
					shdshd =  munode * shdshd;

					index1 = (inl * blksize + JNL)* ndf;

					if ( dynamic ==1) stiffness_mat_fac = 0.5 * positiongamma * dt * dt;
					else stiffness_mat_fac = 1.0;

					if (dirichlet_flag[knl] == 1){
						if (dirichlet_flag[lnl] == 1){
							v[index1]   += (shdshd + shxshx * (lamnode + munode)) * stiffness_mat_fac;
							v[index1]   += sh0sh0;
						}
						if (dirichlet_flag[lnl+1] == 1){
							v[index1+1] += (shyshx * munode + shxshy * lamnode) * stiffness_mat_fac;
						}
						if (dirichlet_flag[lnl+2] == 1 && ndf ==3 ){
							v[index1+2] += (shzshx * munode + shxshz * lamnode) * stiffness_mat_fac;
						}
					} else if (inl==JNL){
						v[index1] = 1.0/hack1e[ie1+inl];
						if (mgnlevels >1) v[index1] = youngs_global/hack1e[ie1+inl];
						v[index1+1] = 0.0;
						if ( nsd==3) v[index1+2] = 0.0;
					}

					index2 = index1 + blksize;
					if (dirichlet_flag[knl+1] == 1){
						if (dirichlet_flag[lnl] == 1){
							v[index2]   += (shxshy * munode + shyshx * lamnode) * stiffness_mat_fac;;
						}
						if (dirichlet_flag[lnl+1] == 1){
							v[index2+1] += (shdshd + shyshy * (lamnode + munode)) * stiffness_mat_fac;
							v[index2+1] += sh0sh0;
						}
						if (dirichlet_flag[lnl+2] == 1 && ndf ==3){
							v[index2+2] += (shzshy * munode + shyshz * lamnode) * stiffness_mat_fac;
						}
					} else if (inl==JNL){
						v[index2]   = 0.0;
						v[index2+1] = 1.0/hack1e[ie1+inl];
						if (mgnlevels >1) v[index2+1] = youngs_global/hack1e[ie1+inl];
						if ( nsd==3) v[index2+2] = 0.0;
					}

					if (nsd ==3){
						index2 = index2 + blksize;
						if (dirichlet_flag[knl+2] == 1){
							if (dirichlet_flag[lnl] == 1){
								v[index2]   += (shxshz * munode + shzshx * lamnode) * stiffness_mat_fac;;
							}
							if (dirichlet_flag[lnl+1] == 1){
								v[index2+1] += (shyshz * munode + shzshy * lamnode) *  stiffness_mat_fac;;
							}
							if (dirichlet_flag[lnl+2] == 1){
								v[index2+2] += (shdshd + shzshz * (lamnode + munode)) * stiffness_mat_fac;;
								v[index2+2] += sh0sh0;
							}
						} else if (inl==JNL){
							v[index2]   = 0.0;
							v[index2+1] = 0.0;
							v[index2+2] = 1.0/hack1e[ie1+inl];
							if (mgnlevels >1) v[index2+2] = youngs_global/hack1e[ie1+inl];
						}
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

	if (penalized_neumann_bc_flag_on == 1) Penalized_Neumann_Contribution_mat(dmmg,jac);

	ierr = PetscFree(det); CHKERRQ(ierr);
	ierr = PetscFree(wq); CHKERRQ(ierr);
	ierr = PetscFree(xq); CHKERRQ(ierr);
	ierr = PetscFree(sh); CHKERRQ(ierr);
	ierr = PetscFree(sq); CHKERRQ(ierr);
	ierr = PetscFree(lambda); CHKERRQ(ierr);
	ierr = PetscFree(mu); CHKERRQ(ierr);
	ierr = PetscFree(dirichlet_flag); CHKERRQ(ierr);
	ierr = PetscFree(v); CHKERRQ(ierr);
	ierr = PetscFree(row); CHKERRQ(ierr);
	ierr = PetscFree(col); CHKERRQ(ierr);

	//ierr=PetscFree(dirichlet_flag);CHKERRQ(ierr);

	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "CreateJacobian"
PetscErrorCode CreateJacobian(stsDMMG dmmg,Mat *jac)
{
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	DA da = (DA)dmmg->dm;
#else
	DM da = dmmg->dm;
#endif
	int mx, my, mz;
	int ierr;

#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	ierr = DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DAGetMatrix(da,MATAIJ,jac);CHKERRQ(ierr);
#elif (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 2))
	ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
	ierr = DMGetMatrix(da, MATAIJ, jac); CHKERRQ(ierr);
#elif (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 4))
	ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
	ierr = DMCreateMatrix(da, MATAIJ, jac); CHKERRQ(ierr);
#else
	ierr = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
	ierr = DMCreateMatrix(da, jac); CHKERRQ(ierr);
#endif
	return 0;
}
