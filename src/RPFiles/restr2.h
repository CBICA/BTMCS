///////////////////////////////////////////////////////////////////////////////////////
// restr2.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef __RESTR2_H
#define __RESTR2_H

//interpolation is bilinear, restriciton =P'
#undef __FUNCT__
#define __FUNCT__ "CreateInterpolation2"
PetscErrorCode CreateInterpolation2(stsDMMG dmmg, Mat *interp)
{
    //RS: Note the input is dmmg instead of da, because assignMAterial requires dmmg.
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
    DA da = (DA)(dmmg->dm);
#else
    DM da = dmmg->dm;
#endif
    PetscInt  M,N,m,n,MD,ND,PD,md,nd,pd;    
    PetscFunctionBegin;
    printf("Inside CreateInterpolation 2. Level = %d , cnt = %d\n",dmmg->nlevels,iCnt);
   //RS: This is coarse grid info
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
    iC(DAGetInfo(da,PETSC_NULL,&MD,&ND,&PD,&md,&nd,&pd,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL));
#else
    iC(DMDAGetInfo(da,PETSC_NULL,&MD,&ND,&PD,&md,&nd,&pd,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL));
#endif
    n= N= nsd*MD*ND*PD;//cols = coarse
    m= M = nsd*(2*(MD-1)+1)*(2*(ND-1)+1)*(2*(PD-1)+1);//rows = fine
    printf("restr2.h iCnt _________________________________________ %d\n",iCnt);
    assert(iCnt <100);
    idata[iCnt].dmmg = dmmg;
    iC(MatCreateShell(PETSC_COMM_WORLD, m, n, PETSC_DETERMINE,PETSC_DETERMINE,(void*)(&idata[iCnt]),interp));
    iC(MatShellSetOperation(*interp,MATOP_MULT, (void(*)(void))Interpolation1MatVec));
    iC(MatShellSetOperation(*interp,MATOP_MULT_ADD, (void(*)(void))AddInterpolation1MatVec));
    iC(MatShellSetOperation(*interp,MATOP_MULT_TRANSPOSE,(void(*)(void))Restriction2MatVec));//Injection
    iC(MatShellSetOperation(*interp,MATOP_MULT_TRANSPOSE_ADD,(void(*)(void))AddRestriction2MatVec));//Injection
    iCnt ++;
    printf("Finished CreateInterpolation 2. Level = %d, cnt = %d\n",dmmg->nlevels,iCnt);
    printf("M=%d, N= %d\n",M,N);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeInterpolation2"
PetscErrorCode ComputeInterpolation2(stsDMMG dmmg,Mat I)
{
    PetscFunctionBegin;
    //printf("Inside ComputeInterpolation 2. Level = %d\n",dmmg->nlevels);
  /*levels*/
    //printf("Finished ComputeInterpolation 2. Level = %d\n",dmmg->nlevels);
    PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "AddRestriction2MatVec"
PetscErrorCode AddRestriction2MatVec(Mat R, Vec v1, Vec v2, Vec v3)
{
    PetscScalar one = 1.0,zero = 0.0;
    Vec tmp;
    PetscFunctionBegin;
    iC(VecDuplicate(v3,&tmp));
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
    iC(VecSet(&zero,tmp));
#else
    iC(VecSet(tmp,zero));
#endif
    iC(Restriction2MatVec(R,v1,tmp));
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
    //VecAXPY(PetscScalar *a,Vec x, Vec y);
    iC(VecAXPY(&one,v2,tmp));//v2 = v2 + tmp = v2 + I*v1
#else
    //VecAXPY(Vec y,PetscScalar a,Vec x);
    iC(VecAXPY(tmp, one, v2));//v2 = v2 + tmp = v2 + I*v1
#endif
    iC(VecCopy(tmp,v3));
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
    iC(VecDestroy(tmp));
#else
    iC(VecDestroy(&tmp));
#endif
    PetscFunctionReturn(0);
}

//This is simply I', so even for boundaries it is not simple injection, they will be scaled
#undef __FUNCT__
#define __FUNCT__ "Restriction2MatVec"
PetscErrorCode Restriction2MatVec(Mat R, Vec in, Vec out)
{
    IMFreeData *data;
    PetscScalar *u;
	//PetscScalar *ep;
    //PetscInt size;
    PetscInt mx,my,mz,fmx,fmy,fmz,ci,cj,ck,CNref,FNref,count;
    PetscScalar Val,hx,hy,hz;
    PetscScalar zero =0.0;
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
    DA da;
#else
    DM da;
#endif
    stsDMMG dmmg;
    PetscFunctionBegin;
    iC(MatShellGetContext(R, (void **)&data));
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
    dmmg = data->dmmg; da = (DA)dmmg->dm; 
    iC(DAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0));//mx,my,mz = number of nodes in each direction
#else
    dmmg = data->dmmg; da = dmmg->dm; 
    iC(DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0));//mx,my,mz = number of nodes in each direction
#endif
    fmx = (2*(mx-1)) + 1;fmy = (2*(my-1)) + 1;fmz = (2*(mz-1)) + 1;
    hx=Lx/(((PetscScalar)mx)-1.0);hy=Ly/(((PetscScalar)my)-1.0);
    hz=Lz/(((PetscScalar)mz)-1.0);
    //iC(VecGetArray(data->epsilon,&ep));
    

//    iC(VecGetSize(in,&size));   
//    printf("Rest2MatVec size %d %x\n",size, *((void**)(in->data))); // *((PetscScalar **)(in)->data));
//    printf("Rest2MatVec size %d %d\n",size, in->petscnative);
    
    iC(VecGetArray(in,&u));//fine grid vector
     
    
#if (PETSC_VERSION_MAJOR == 2) && (PETSC_VERSION_MINOR <= 2)
    iC(VecSet(&zero,out)); //Initialize out
#else
    iC(VecSet(out,zero)); //Initialize out
#endif
    //Loop through All Nodes on coarse grid.
    // loc(ci,cj,ck) gives the coarse grid node number.    
   //cinternal Nodes: ci=1:mx-2,cj=1:my-2,ck=1:mz-2
    for (ck=1; ck<=mz-2; ck++){
        for (cj=1; cj<=my-2; cj++){
            for(ci=1; ci<=mx-2; ci++){
            
            	//CH: try to account for 3 degrees of freedom/node
        	for (count=0;count<=2;count++){
            
                CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
                //floc(i,j,k):=0-based global node number of this node on the fine grid
                //injection for old nodes: i-->2*fi
                FNref = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Val = u[FNref];
                //The contribution is 0.5 for (i+di/2,j,k)
                FNref = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Val += 0.5*u[FNref];
                FNref = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Val += 0.5*u[FNref];
                FNref = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Val += 0.5*u[FNref];
                FNref = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Val += 0.5*u[FNref];
                FNref = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Val += 0.5*u[FNref];
                FNref = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Val += 0.5*u[FNref];
                //The contribution is 0.25 for (i+di/2,j+dj/2,k)
                FNref =  nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Val += 0.25*u[FNref];
                FNref =  nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Val += 0.25*u[FNref];
                FNref =  nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Val += 0.25*u[FNref];
                FNref =  nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Val += 0.25*u[FNref];
                FNref =  nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Val += 0.25*u[FNref];
                FNref =  nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Val += 0.25*u[FNref];
                FNref =  nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Val += 0.25*u[FNref];
                FNref =  nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Val += 0.25*u[FNref];
                FNref =  nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Val += 0.25*u[FNref];
                FNref =  nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Val += 0.25*u[FNref];
                FNref =  nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Val += 0.25*u[FNref];
                FNref =  nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Val += 0.25*u[FNref];
                //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
                FNref =  nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1)+count; Val += 0.125*u[FNref];
                FNref =  nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1)+count; Val += 0.125*u[FNref];
                FNref =  nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1)+count; Val += 0.125*u[FNref];
                FNref =  nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1)+count; Val += 0.125*u[FNref];
                FNref =  nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1)+count; Val += 0.125*u[FNref];
                FNref =  nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1)+count; Val += 0.125*u[FNref];
                FNref =  nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1)+count; Val += 0.125*u[FNref];
                FNref =  nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1)+count; Val += 0.125*u[FNref];
                iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
                
                }//end count (degreess of freedom/node)
                
            }//end for ci
        }//end for cj
    }//end for ck
    
    //Surface Nodes: z=0
    ck=0;
    for (cj=1; cj<=my-2; cj++){
        for(ci=1; ci<=mx-2; ci++){
        
            //CH: try to account for 3 degrees of freedom/node
            for (count=0;count<=2;count++){
        
            CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
            //floc(i,j,k):=0-based global node number of this node on the fine grid
            //injection for old nodes: i-->2*fi
            FNref = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Val = u[FNref];
            //The contribution is 0.5 for (i+di/2,j,k)
            FNref = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Val += 0.5*u[FNref];
            FNref = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Val += 0.5*u[FNref];
            FNref = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Val += 0.5*u[FNref];
            FNref = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Val += 0.5*u[FNref];
            FNref = nsd*(floc(2*ci,2*cj,(2*ck)+1) - 1)+count; Val += 0.5*u[FNref];
            //The contribution is 0.25 for (i+di/2,j+dj/2,k)
            FNref =  nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc((2*ci)+1,2*cj,(2*ck)+1) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc((2*ci)-1,2*cj,(2*ck)+1) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc(2*ci,(2*cj)+1,(2*ck)+1) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc(2*ci,(2*cj)-1,(2*ck)+1) - 1)+count; Val += 0.25*u[FNref];
            //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
            FNref =  nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1)+count; Val += 0.125*u[FNref];
            FNref =  nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1)+count; Val += 0.125*u[FNref];
            FNref =  nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1)+count; Val += 0.125*u[FNref];
            FNref =  nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1)+count; Val += 0.125*u[FNref];
            iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
            
            }//end count (degreess of freedom/node)
            
        }//end for ci
    }//end for cj
    
    //Surface Nodes: z=mz-1
    ck=mz-1;
    for (cj=1; cj<=my-2; cj++){
        for(ci=1; ci<=mx-2; ci++){
        
            //CH: try to account for 3 degrees of freedom/node
            for (count=0;count<=2;count++){
        
            CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
            //floc(i,j,k):=0-based global node number of this node on the fine grid
            //injection for old nodes: i-->2*fi
            FNref = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Val = u[FNref];
            //The contribution is 0.5 for (i+di/2,j,k)
            FNref = nsd*(floc((2*ci)+1,2*cj,2*ck) - 1)+count; Val += 0.5*u[FNref];
            FNref = nsd*(floc((2*ci)-1,2*cj,2*ck) - 1)+count; Val += 0.5*u[FNref];
            FNref = nsd*(floc(2*ci,(2*cj)+1,2*ck) - 1)+count; Val += 0.5*u[FNref];
            FNref = nsd*(floc(2*ci,(2*cj)-1,2*ck) - 1)+count; Val += 0.5*u[FNref];
            FNref = nsd*(floc(2*ci,2*cj,(2*ck)-1) - 1)+count; Val += 0.5*u[FNref];
            //The contribution is 0.25 for (i+di/2,j+dj/2,k)
            FNref =  nsd*(floc((2*ci)+1,(2*cj)+1,2*ck) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc((2*ci)+1,(2*cj)-1,2*ck) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc((2*ci)+1,2*cj,(2*ck)-1) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc((2*ci)-1,(2*cj)+1,2*ck) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc((2*ci)-1,(2*cj)-1,2*ck) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc((2*ci)-1,2*cj,(2*ck)-1) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc(2*ci,(2*cj)+1,(2*ck)-1) - 1)+count; Val += 0.25*u[FNref];
            FNref =  nsd*(floc(2*ci,(2*cj)-1,(2*ck)-1) - 1)+count; Val += 0.25*u[FNref];
            //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
            FNref =  nsd*(floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1)+count; Val += 0.125*u[FNref];
            FNref =  nsd*(floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1)+count; Val += 0.125*u[FNref];
            FNref =  nsd*(floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1)+count; Val += 0.125*u[FNref];
            FNref =  nsd*(floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1)+count; Val += 0.125*u[FNref];
            iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
            
            }//end count (degreess of freedom/node)
            
        }//end for ci
    }//end for cj
    
    //Surface Nodes: y=0
    cj=0;
    for (ck=1; ck<=my-2; ck++){
        for(ci=1; ci<=mx-2; ci++){
        
            //CH: try to account for 3 degrees of freedom/node
            for (count=0;count<=2;count++){
        
            CNref = nsd*(loc(ci,cj,ck) - 1)+count;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
            //floc(i,j,k):=0-based global node number of this node on the fine grid
            //injection for old nodes: i-->2*fi
            FNref = nsd*(floc(2*ci,2*cj,2*ck) - 1)+count; Val = u[FNref];
            //The contribution is 0.5 for (i+di/2,j,k)
            FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            //The contribution is 0.25 for (i+di/2,j+dj/2,k)
            FNref =  floc((2*ci)+1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)+1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)+1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)-1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)-1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)-1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc(2*ci,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc(2*ci,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
            FNref =  floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            FNref =  floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            FNref =  floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            FNref =  floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
            
            }//end count (degreess of freedom/node)
            
        }//end for ci
    }//end for ck
    
    //Surface Nodes: y=my-1
    cj=my-1;
    for (ck=1; ck<=my-2; ck++){
        for(ci=1; ci<=mx-2; ci++){
        
            //CH: try to account for 3 degrees of freedom/node
            for (count=0;count<=2;count++){
        
            CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
            //floc(i,j,k):=0-based global node number of this node on the fine grid
            //injection for old nodes: i-->2*fi
            CNref=nsd*CNref+count;
            
            FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
            //The contribution is 0.5 for (i+di/2,j,k)
            FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            //The contribution is 0.25 for (i+di/2,j+dj/2,k)
            FNref =  floc((2*ci)+1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)+1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)+1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)-1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)-1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)-1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc(2*ci,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc(2*ci,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
            FNref =  floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            FNref =  floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            FNref =  floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            FNref =  floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
            
            }//end count (degreess of freedom/node)
            
        }//end for ci
    }//end for ck
    
    
    //Surface Nodes: x=0
    ci=0;
    for (ck=1; ck<=my-2; ck++){
        for(cj=1; cj<=mx-2; cj++){
        
        
            //CH: try to account for 3 degrees of freedom/node
            for (count=0;count<=2;count++){
        
            CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
            //floc(i,j,k):=0-based global node number of this node on the fine grid
            //injection for old nodes: i-->2*fi
            CNref=nsd*CNref+count;
            
            FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
            //The contribution is 0.5 for (i+di/2,j,k)
            FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            //The contribution is 0.25 for (i+di/2,j+dj/2,k)
            FNref =  floc((2*ci)+1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)+1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)+1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)+1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc(2*ci,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc(2*ci,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc(2*ci,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc(2*ci,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
            FNref =  floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            FNref =  floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            FNref =  floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            FNref =  floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
            
            }//end count (degreess of freedom/node)
            
        }//end for cj
    }//end for ck
    
    //Surface Nodes: x=mx-1
    ci=mx-1;
    for (ck=1; ck<=my-2; ck++){
        for(cj=1; cj<=mx-2; cj++){
        
            //CH: try to account for 3 degrees of freedom/node
            for (count=0;count<=2;count++){
        
            CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
            //floc(i,j,k):=0-based global node number of this node on the fine grid
            //injection for old nodes: i-->2*fi
            CNref=nsd*CNref+count;
            
            FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
            //The contribution is 0.5 for (i+di/2,j,k)
            FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
            //The contribution is 0.25 for (i+di/2,j+dj/2,k)
            FNref =  floc((2*ci)-1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)-1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)-1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc((2*ci)-1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc(2*ci,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc(2*ci,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc(2*ci,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            FNref =  floc(2*ci,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
            //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
            FNref =  floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            FNref =  floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            FNref =  floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            FNref =  floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
            iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
            
            }//end count (degreess of freedom/node)
            
        }//end for cj
    }//end for ck
    
    
    //Edge Nodes: x=0, y=0
    ci=0;cj=0;
    for (ck=1; ck<=my-2; ck++){
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)+1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        FNref =  floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }//end for ck
    
    //Edge Nodes: x=0, y=my-1
    //n= 4,8
    ci=0;cj=my-1;
    for (ck=1; ck<=my-2; ck++){
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)+1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        FNref =  floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }//end for ck
    
    //Edge Nodes: x=0, z=0
    //n= 1,4
    ci=0;ck=0;
    for (cj=1; cj<=my-2; cj++){
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)+1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        FNref =  floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }//end for cj
    
    //Edge Nodes: x=0, z=mz-1
    //n= 5,8
    ci=0;ck=mz-1;
    for (cj=1; cj<=my-2; cj++){
    
         //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)+1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        FNref =  floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
        
    }//end for cj
    
    //Edge Nodes: x=mx-1,y=0
    //n= 2,6
    ci=mx-1;cj=0;
    for (ck=1; ck<=mz-2; ck++){
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)-1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }//end for ck
    
    //Edge Nodes: x=mx-1,y=my-1
    //n= 3,7
    ci=mx-1;cj=my-1;
    for (ck=1; ck<=mz-2; ck++){
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
        
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)-1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }//end for ck
    
    //Edge Nodes: x=mx-1,z=0
    //n= 2,3
    ci=mx-1;ck=0;
    for (cj=1; cj<=my-2; cj++){
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)-1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }//end for cj
    
    //Edge Nodes: x=mx-1,z=mz-1
    //n= 6,7
    ci=mx-1;ck=mz-1;
    for (cj=1; cj<=my-2; cj++){
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)-1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }//end for cj
    
    //Edge Nodes: y=0, z= 0
    //n=1,2
    cj=0;ck=0;
    for (ci=1; ci<=mx-2; ci++){
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)+1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }//end for ci
    
    //Edge Nodes: y=0, z= mz-1
    //n=5,6
    cj=0;ck=mz-1;
    for (ci=1; ci<=mx-2; ci++){
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)+1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }//end for ci
    
    //Edge Nodes: y=my-1, z=0
    cj=my-1;ck=0;
    for (ci=1; ci<=mx-2; ci++){
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)+1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }//end for ci
    
    //Edge Nodes: y=my-1, z=mz-1
    cj=my-1;ck=mz-1;
    for (ci=1; ci<=mx-2; ci++){
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)+1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        FNref =  floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }//end for ci
    
    //Corner x=0,y=0,z=0
    //n=1
    ci=0;cj=0;ck=0;
    {
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)+1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)+1,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }
    
    //Corner x=mx-1,y=0,z=0
    //n=2
    ci=mx-1;cj=0;ck=0;
    {
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)-1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)-1,(2*cj)+1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
        
    }
    
    //Corner x=mx-1,y=my-1,z=0
    //n=3
    ci=mx-1;cj=my-1;ck=0;
    {
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)-1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)-1,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }
    
    //Corner x=0,y=my-1,z=0
    //n=4
    ci=0;cj=my-1;ck=0;
    {
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)+1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)+1,(2*cj)-1,(2*ck)+1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }
    
    //Corner x=0,y=0,z=mz-1
    //n=5
    ci=0;cj=0;ck=mz-1;
    {
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
         
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)+1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)+1,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
    }
    
    //Corner x=mx-1,y=0,z=mz-1
    //n=6
    ci=mx-1;cj=0;ck=mz-1;
    {
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)-1,(2*cj)+1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)-1,(2*cj)+1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }
    
    //Corner x=mx-1,y=my-1,z=mz-1
    //n=7
    ci=mx-1;cj=my-1;ck=mz-1;
    {
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)-1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)-1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)-1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)-1,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }
    
    //Corner x=0,y=my-1,z=mz-1
    //n=8
    ci=0;cj=my-1;ck=mz-1;
    {
    
        //CH: try to account for 3 degrees of freedom/node
        for (count=0;count<=2;count++){
    
        CNref = loc(ci,cj,ck) - 1;//loc(i,j,k):= 0-based global node number of this node on the coarse grid
        //floc(i,j,k):=0-based global node number of this node on the fine grid
        //injection for old nodes: i-->2*fi
        CNref=nsd*CNref+count;
        
        FNref = floc(2*ci,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val = u[FNref];
        //The contribution is 0.5 for (i+di/2,j,k)
        FNref = floc((2*ci)+1,2*cj,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        FNref = floc(2*ci,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.5*u[FNref];
        //The contribution is 0.25 for (i+di/2,j+dj/2,k)
        FNref =  floc((2*ci)+1,(2*cj)-1,2*ck) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc((2*ci)+1,2*cj,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        FNref =  floc(2*ci,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.25*u[FNref];
        //The contribution is 0.125 for (i+di/2,j+dj/2,k+dk/2)
        FNref =  floc((2*ci)+1,(2*cj)-1,(2*ck)-1) - 1; FNref=nsd*FNref+count; Val += 0.125*u[FNref];
        iC(VecSetValue(out,CNref,Val,INSERT_VALUES));
        
        }//end count (degreess of freedom/node)
        
    }
    
    //iC(VecRestoreArray(data->epsilon,&ep));
    iC(VecRestoreArray(in,&u));//fine grid vector
    iC(VecAssemblyBegin(out));//coarse grid vector
    iC(VecAssemblyEnd(out));
    PetscFunctionReturn(0);
}

#endif

