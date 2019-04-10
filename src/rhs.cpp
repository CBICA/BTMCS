///////////////////////////////////////////////////////////////////////////////////////
// rhs.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include "function.h"
#include "petscviewer.h"


// Version 1.1 of STSolver by F. Abraham and G. Biros 
// August 2005

extern void penalized_neumann_boundary(); // ReadFile.c
extern int Penalized_Neumann_Contribution(stsDMMG, Vec);


#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
int ComputeRHS(stsDMMG dmmg,Vec b)
{
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	DA da = (DA)dmmg->dm;
#else
	DM da = dmmg->dm;
#endif
	MatStencil   *row, *col;
	PetscScalar h;
	int count, counti, countj, countk, ie;
	int ierr,i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs,ie2;
	int index1, index2, iq,ie1,inl,lnl,JNL,knl;
	int *idx,*globindex,jacobian_flag,*dirichlet_flag;

	Vec b1, b2;

	double norm;
	double *v, **xq,**sq,***sh,*det,*wq;
	double shyshy, shxshx, shzshz, sh0i,sh0j;
	double shyi,shxi,shzi,shdshd,sh0sh0, shyj,shxj,shzj;
	double shxshy,shxshz, shyshx,shyshz, shzshy,shzshx,eff0;
	double du,dv,dw;
	double yc1;
	double nodalforcex,nodalforcey,nodalforcez;
	double *lambda,*mu, lamnode, munode, lampmunode,*barray;
	//double *vjacprec;
	double stiffness_mat_fac;

	//test Feby's force
	double *forceelem;

	stsData *ptr;
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
	ierr = DAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DAGetInfo((DA)dmmg->dm,0,&mx,&my,&mz,0,0,0,0,0,0,0);CHKERRQ(ierr);
#else
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); CHKERRQ(ierr);
	ierr = DMDAGetInfo(dmmg->dm,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0); CHKERRQ(ierr);
#endif

	PetscFunctionBegin;

	ne = (my-1) * (mx-1);
	if (nsd ==3) ne = ne * (mz-1);
	nnc = (xm) * (ym) ;
	if (nsd == 3) nnc = nnc * (zm);
	Hx = Lx/ (PetscReal)(mx-1); 
	Hy = Ly/ (PetscReal)(my-1); 
	Hz = Lz/ (PetscReal)(mz-1);
	mxg = mx; myg = my; mzg = mz;

	//test against Feby's force
	ierr = PetscMalloc(ne*nen*nsd*sizeof(PetscReal),&forceelem); CHKERRQ(ierr);

	ierr = PetscMalloc(nen*sizeof(int),&globindex); CHKERRQ(ierr);

	//ierr = PetscMalloc(nnc*nsd*sizeof(PetscReal),&forcenode); CHKERRQ(ierr);

	ierr = PetscMalloc(ne*sizeof(PetscReal),&lambda); CHKERRQ(ierr);
	ierr = PetscMalloc(ne*sizeof(PetscReal),&mu); CHKERRQ(ierr);

	ierr = PetscMalloc(nquad*sizeof(PetscReal),&det); CHKERRQ(ierr);
	ierr = PetscMalloc(nquad*sizeof(PetscReal),&wq); CHKERRQ(ierr);
	ierr = PetscMalloc(nquad*sizeof(PetscReal),&xq); CHKERRQ(ierr);
	ierr = PetscMalloc(nquad*sizeof(PetscReal),&sh); CHKERRQ(ierr);
	ierr = PetscMalloc(nen*sizeof(PetscReal),&sq); CHKERRQ(ierr);
	for (j=0;j<nen;j++) {
		ierr = PetscMalloc(nquad*sizeof(PetscReal),&sq[j]); CHKERRQ(ierr);
	}
	for (j=0;j<nquad;j++) {
		ierr = PetscMalloc(nsd*sizeof(PetscReal),&xq[j]); CHKERRQ(ierr);
		ierr = PetscMalloc(nen*sizeof(PetscReal),&sh[j]); CHKERRQ(ierr);
		for (i=0;i<nen;i++) {
			ierr = PetscMalloc(nsd*sizeof(PetscReal),&sh[j][i]); CHKERRQ(ierr);
		}
	}

	ierr = PetscMalloc(nen*nsd*sizeof(int),&dirichlet_flag); CHKERRQ(ierr);
	ierr = PetscMalloc(blksize*blksize*sizeof(int),&idx); CHKERRQ(ierr);
	ierr = PetscMalloc(blksize*blksize*sizeof(PetscReal),&v); CHKERRQ(ierr);
	ierr = PetscMalloc(blksize*sizeof(MatStencil),&row); CHKERRQ(ierr);
	ierr = PetscMalloc(blksize*sizeof(MatStencil),&col); CHKERRQ(ierr);

	quad(xq, sq, sh,det,wq,Hx,Hy,Hz); 
	inhomogeneous_material(lambda, mu, ne, xm, ym, Hx, Hy, Hz);
	//printf("rhs ne  %d\n",ne);
	createien(xm, ym);

	if (matrixfree == 1){
		ptr = &shellData[mgnlevels-1];
		//ierr = VecGetArray1d(ptr->vjacprec,nnc*nsd,0,&vjacprec); CHKERRQ(ierr);
	}

	h = 0.0;
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
	ierr = VecSet(&h, b);CHKERRQ(ierr);
#else
	ierr = VecSet(b, h); CHKERRQ(ierr);
#endif
	counti=1; countj=1; countk=1;
	i = xs; j = ys; k = zs;
	jacobian_flag =0;

	shzshz = 0.0; shyshz = 0.0; shzshy = 0.0;shxshz = 0.0;
	shzshx = 0.0; shzi = 0.0;
	sh0i   = 0.0; yc1 = 0.0;
	jacobian_flag = 0;

	/* FORCE TERMS - CH */

	//  nodalforcevector(forcenode,nnc,mx,my,Hx,Hy,Hz);
	//  counti=1; countj=1; countk=1;
	//  for (ie=0;ie < ne;ie++){
	//  
	//         ie1 = ie * nen * nsd;     
	//         i = counti-1; j = countj-1; k = countk-1;
	//         rowcolindex(i,j,k,row,col,jacobian_flag);
	//         
	//         globindex[0]=i+j*mx+k*mx*my;
	//         globindex[1]=(i+1)+j*mx+k*mx*my;
	//         globindex[2]=i+(j+1)*mx+k*mx*my;
	//         globindex[3]=(i+1)+(j+1)*mx+k*mx*my;
	//         globindex[4]=i+j*mx+(k+1)*mx*my;
	//         globindex[5]=(i+1)+j*mx+(k+1)*mx*my;
	//         globindex[6]=i+(j+1)*mx+(k+1)*mx*my;
	//         globindex[7]=(i+1)+(j+1)*mx+(k+1)*mx*my;
	//         
	//         for (count=0;count<blksize;count++)v[count]=0;
	//         for (iq=0;iq<nquad;iq++){
	//                nodalforcex = 0.0;
	//                nodalforcey = 0.0;
	//                nodalforcez = 0.0;
	//                for (inl=0;inl<nen;inl++){
	//                                
	//                     sh0i = sq[inl][iq];
	//                  
	//		     vecglobindex=nsd*globindex[inl];
	//		     	
	//		     nodalforcex += sh0i * forcenode[vecglobindex];
	//                     nodalforcey += sh0i * forcenode[vecglobindex+1];
	//                     nodalforcez += sh0i * forcenode[vecglobindex+2];
	//		
	//                }
	//                eff0 = wq[iq] * det[iq];
	//                for (inl=0;inl<nen;inl++){
	//                      knl    = inl * ndf;
	//                      index1 = row[knl].j*mx+row[knl].i;
	//                      if (ndf ==3) index1 += row[knl].k*mx*my;
	//                      index1 = index1 * ndf;
	//                      index2 = knl;
	//                      idx[index2]   =  index1;
	//                      idx[index2+1] =  index1 + 1;
	//                      v[index2]     += sh0i * nodalforcex * eff0;
	//                      v[index2+1]   += sh0i * nodalforcey * eff0;
	//                      if (nsd == 3){
	//			 idx[index2+2] =  index1 + 2;
	//                      	 v[index2+2]   += sh0i * nodalforcez * eff0;
	//		      }
	//                }
	//             }
	//             ierr = VecSetValues(b, blksize, idx, v, ADD_VALUES);CHKERRQ(ierr);
	//             counti++;
	//             if (counti==xm){ countj++; counti=1; }
	//             if (countj==ym){ counti=1; countj=1; countk++; }
	//  }

	/* FORCE TERMS - FA */
	//test Feby's force
	//nodalforcevector(forcenode,nnc,mx,my,Hx,Hy,Hz);

	//printf("forcenodex forcenodey forcenodez %g %g %g\n",forcenode[334902],forcenode[334903],forcenode[334904]);

	gatherreal(forcenode, forceelem);
	counti=1; countj=1; countk=1;
	for (ie=0;ie < ne;ie++){
		ie1 = ie * nen * nsd;
		i = counti-1; j = countj-1; k = countk-1;
		rowcolindex(i,j,k,row,col,jacobian_flag);
		for (count=0;count<blksize;count++)v[count]=0;
		for (iq=0;iq<nquad;iq++){
			nodalforcex = 0.0;
			nodalforcey = 0.0;
			nodalforcez = 0.0;
			for (inl=0;inl<nen;inl++){
				ie2 = ie1 + inl * nsd;
				sh0i = sq[inl][iq];
				nodalforcex += sh0i * forceelem[ie2];
				nodalforcey += sh0i * forceelem[ie2+1];
				nodalforcez += sh0i * forceelem[ie2+2];
			}
			eff0 = wq[iq] * det[iq];
			for (inl=0;inl<nen;inl++){
				knl    = inl * ndf;
				index1 = row[knl].j*mx+row[knl].i;
				if (ndf ==3) index1 += row[knl].k*mx*my;
				index1 = index1 * ndf;
				index2 = knl;
				idx[index2]   =  index1;
				idx[index2+1] =  index1 + 1;
				v[index2]     += sh0i * nodalforcex * eff0;
				v[index2+1]   += sh0i * nodalforcey * eff0;
				if (nsd == 3){
					idx[index2+2] =  index1 + 2;
					v[index2+2]   += sh0i * nodalforcez * eff0;
				}
			}
		}
		ierr = VecSetValues(b, blksize, idx, v, ADD_VALUES);CHKERRQ(ierr);
		counti++;
		if (counti==xm){ countj++; counti=1; }
		if (countj==ym){ counti=1; countj=1; countk++; }
	}

	VecAssemblyBegin(b);
	VecAssemblyEnd(b);

	VecCreateSeq(PETSC_COMM_WORLD,nnc*ndf,&b1);
	ierr = VecGetArray1d(b1,nnc*nsd,0,&barray);

	/*Penalized Embedded Neumann BC */
	if (penalized_neumann_bc_flag_on == 1){
		VecCreateSeq(PETSC_COMM_WORLD,nnc*ndf,&b2);
		penalized_neumann_boundary();
		Penalized_Neumann_Contribution(dmmg,b2);
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
		//VecAXPY(PetscScalar *a,Vec x, Vec y);
		ierr = VecAXPY(&one,b2,b);CHKERRQ(ierr);
#else
		//VecAXPY(Vec y,PetscScalar a,Vec x);
		ierr = VecAXPY(b, one, b2); CHKERRQ(ierr);
#endif
		VecNorm(b2,NORM_2,&norm);
	}

	/* CARTESIAN BOUNDARY DIRICHLET BC*/

	/*CartesianDirichlet(dirivec,xm, ym,zm);
	for (inc =0;inc<nnc;inc++){
	knl = inc * nsd;
	barray[knl]     = dirivec[knl];
	barray[knl+1]   = dirivec[knl+1];
	if ( nsd ==3)barray[knl+2]   = dirivec[knl+2];
	}

	ierr = VecRestoreArray1d(b1,nnc*nsd,0,&barray);
	ierr = VecAXPY(&one,b1,b);CHKERRQ(ierr);
	*/

	/* Dirichlet CORRECTION: This is for the way matrix-based methods are implemented.
	Note that Dirichlet refers only to Cartresian boundary walls
	and not embedded Dirichlet walls */

	if (matrixfree == 0)
	{
		ierr=PetscMalloc(nnc*nsd*sizeof(PetscReal),&dirivec);
		ierr=PetscMalloc(ne*nen*ndf*sizeof(PetscReal),&diriel);

		gatherreal(dirivec, diriel);
		h = 0.0;
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
		ierr = VecSet(&h, b1);CHKERRQ(ierr);
#else
		ierr = VecSet(b1, h); CHKERRQ(ierr);
#endif

		i = xs; j = ys; k = zs;
		counti=1; countj=1; countk=1;
		for (ie=0;ie<ne;ie++){
			ie1 = ie * nen;
			ie2 = ie1* ndf;
			i = counti-1; j = countj-1; k = countk-1;
			rowcolindex(i,j,k,row,col,jacobian_flag);
			for (count=0;count<blksize;count++)v[count]=0;
			dirichlet_identifier(dirichlet_flag, row,mx,my,mz);
			for (iq=0;iq<nquad;iq++){
				count = 0;
				eff0 = wq[iq] * det[iq];
				lamnode = lambda[ie];
				munode = mu[ie];
				lampmunode = lamnode + munode;

				du = 0.0; dv = 0.0; dw = 0.0;
				for (inl=0;inl < nen;inl++){
					knl = inl * ndf;
					index1 = row[knl].j*xm+row[knl].i;
					if (ndf ==3) index1 += row[knl].k*xm*ym;
					index1 = index1 * ndf;
					index2 = knl;
					idx[index2]   =  index1;
					idx[index2+1] =  index1 + 1;
					if (nsd==3) idx[index2+2] = index1 + 2;

					sh0i = sq[inl][iq];
					shxi = sh[iq][inl][xsd] * eff0;
					shyi = sh[iq][inl][ysd] * eff0;
					if (ndf ==3) shzi = sh[iq][inl][zsd] * eff0;
					for (JNL=0;JNL<nen;JNL++){
						lnl = JNL * ndf;
						if (mgnlevels>1){
							du  = diriel[ie2+lnl]    * youngs_global ;
							dv  = diriel[ie2+lnl+1]  * youngs_global ;
							dw  = diriel[ie2+lnl+2]  * youngs_global ;
						}else{
							du  = diriel[ie2+lnl];
							dv  = diriel[ie2+lnl+1];
							dw  = diriel[ie2+lnl+2];
						}

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
						shdshd = munode * shdshd;


						if ( dynamic ==1) stiffness_mat_fac = 0.5 * positiongamma * dt * dt;
						else stiffness_mat_fac = 1.0;

						if (dirichlet_flag[knl] == 1){
							if (dirichlet_flag[lnl] == 0){
								v[index2]   -= (shdshd + shxshx * lampmunode)   * stiffness_mat_fac * du;
								v[index2]   -= sh0sh0 * du;
							}
							if (dirichlet_flag[lnl+1] == 0){
								v[index2] -= (shyshx * munode + shxshy * lamnode) * stiffness_mat_fac * dv;
							}
							if (dirichlet_flag[lnl+2] == 0 && ndf ==3){
								v[index2] -= (shzshx * munode + shxshz * lamnode) * stiffness_mat_fac *dw;
							}
						} 
						if (dirichlet_flag[knl+1] == 1){
							if (dirichlet_flag[lnl] == 0){
								v[index2+1]   -= (shxshy * munode + shyshx * lamnode) * stiffness_mat_fac * du;
							}
							if (dirichlet_flag[lnl+1] == 0){
								v[index2+1] -= (shdshd + shyshy * lampmunode) * stiffness_mat_fac * dv;
								v[index2+1] -= sh0sh0 * dv;
							}
							if (dirichlet_flag[lnl+2] == 0 && ndf ==3){
								v[index2+1] -= (shzshy * munode + shyshz * lamnode) * stiffness_mat_fac * dw;
							}
						}
						if (nsd ==3){
							if (dirichlet_flag[knl+2] == 1){
								if (dirichlet_flag[lnl] == 0){
									v[index2+2] -= (shxshz * munode + shzshx * lamnode) * stiffness_mat_fac * du;
								}
								if (dirichlet_flag[lnl+1] == 0){
									v[index2+2] -= (shyshz * munode + shzshy * lamnode) * stiffness_mat_fac * dv;
								}
								if (dirichlet_flag[lnl+2] == 0){
									v[index2+2] -= (shdshd + shzshz * lampmunode) * stiffness_mat_fac * dw;
									v[index2+2] -= sh0sh0 * dw;
								}
							}
						}			   		
					}
				}
			}
			ierr = VecSetValues(b1, blksize, idx, v, ADD_VALUES);CHKERRQ(ierr);
			counti++;
			if (counti==xm){ countj++; counti=1; }
			if (countj==ym){ counti=1; countj=1; countk++; }

		}

		VecAssemblyBegin(b1);
		VecAssemblyEnd(b1);
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
		//VecAXPY(PetscScalar *a,Vec x, Vec y);
		ierr = VecAXPY(&one,b1,b);CHKERRQ(ierr);
#else
		//VecAXPY(Vec y,PetscScalar a,Vec x);
		ierr = VecAXPY(b, one, b1);CHKERRQ(ierr);
#endif

		PetscFree(dirivec);
		PetscFree(diriel);
	}

	// if (matrixfree == 1 && jacobiprecflag == 1){
	//  ierr = VecGetArray1d(b,nnc*nsd,0,&barray);
	//  for (inc =0;inc<nnc;inc++){
	//    knl = inc * nsd;
	//    barray[knl]     = barray[knl]   * vjacprec[knl];
	//    barray[knl+1]   = barray[knl+1] * vjacprec[knl+1];
	//    barray[knl+2]   = barray[knl+2] * vjacprec[knl+2];
	//   }
	//   VecRestoreArray1d(ptr->vjacprec,nnc*nsd,0,&vjacprec);
	//  }

	VecRestoreArray1d(b,nnc*nsd,0,&barray);

	//  ierr=VecView(b,PETSC_VIEWER_STDOUT_WORLD);

	//  ierr = VecGetArray1d(b,nnc*nsd,0,&barray);
	//  xval = 0.0; yval = 0.0; zval = 0.0;
	//  fp=fopen("rhs.txt","w");
	//  for (inc =0;inc<nnc;inc++){
	//    knl = inc * nsd;
	//    fprintf(fp,"%g\n",sqrt(barray[knl]*barray[knl]+barray[knl+1]*barray[knl+1]+barray[knl+2]*barray[knl+2]));
	//    xval += barray[knl];;
	//    yval += barray[knl+1];;
	//    zval += barray[knl+2];;
	//  }
	//  fclose(fp);
	//   VecRestoreArray1d(b,nnc*nsd,0,&barray);
	//  printf("Xerr Yerr Zerr %g %g %g\n",xval,yval,zval);
	// 
	//ierr = VecGetArray1d(b,nnc*nsd,0,&barray);
	// xval = 1000000000000000.0; xval1 = -10000000000000.0 ;
	// for (inc =0;inc<nnc;inc++){
	//   knl = inc * nsd;
	//   if ( xval > barray[knl]){
	//        incout = inc;
	//        xval = barray[knl];
	//   }
	//   if ( xval1 < barray[knl]){
	//        incout1 = inc;
	//        xval1 = barray[knl];
	//   }
	// }
	//  printf("incout %d\n",incout);
	//  printf("incout1 %d\n", incout1);
	//  VecRestoreArray1d(b,nnc*nsd,0,&barray);

	PetscFree(det);
	PetscFree(wq);
	PetscFree(xq);
	PetscFree(sh);
	PetscFree(sq);
	PetscFree(globindex);
	PetscFree(lambda);
	PetscFree(mu);
	PetscFree(dirichlet_flag);
	PetscFree(idx);
	PetscFree(v);
	PetscFree(row);
	PetscFree(col);
	//PetscFree(forcenode);
	PetscFree(ien);

	//test Feby's force
	PetscFree(forceelem);

	PetscFunctionReturn(0);
}
