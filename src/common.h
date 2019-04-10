///////////////////////////////////////////////////////////////////////////////////////
// common.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef __COMMON_H_
#define __COMMON_H_

/*
static int _i_n_element[8] = {0,1,1,0,0,1,1,0};
/static int _j_n_element[8] = {0,0,1,1,0,0,1,1};
static int _k_n_element[8] = {0,0,0,0,1,1,1,1};
*/

#define iC(fun) {_internal_ierr = fun; CHKERRQ(_internal_ierr);}
#define loc(i,j,k) (((i>=0) && (i<mx) && (j>=0) && (j<my) && (k>=0) && (k<mz)) ? ((i) + (j)*(mx) + (k)*(mx)*(my) + 1) : (-1) )
//mx,my,mz = number of nodes in each direction 
#define floc(fi,fj,fk) (((fi>=0) && (fi<fmx) && (fj>=0) && (fj<fmy) && (fk>=0) && (fk<fmz))? ((fi)+(fj)*(fmx)+(fk)*(fmx)*(fmy)+1) : (-1) )
//fmx,fmy,fmz = number of nodes in each direction 
#define iinvloc(n) fmod(n-1,mx)
#define jinvloc(n,i) fmod((n-1-i)/mx,my)
#define kinvloc(n,i,j) ((((n-1-i)/mx) - j)/my)
//The extra 1 is for 1 based indexing of node numbers. i,j,k are all 0-based.
#define liinvloc(n)_i_n_element[n-1] 
#define ljinvloc(n) _j_n_element[n-1]
#define lkinvloc(n) _k_n_element[n-1]
//l stands for local, so 1<=n<=8 (single element)
#define sliinvloc(n) fmod(n-1,3)
#define sljinvloc(n,i) fmod((n-1-i)/3,3)
#define slkinvloc(n,i,j) ((((n-1-i)/3) - j)/3)
//sl stands for local stencil, so mx=my=mz =3 (8 elements surroudnig an internal node)
#define Ni(N) iinvloc(N)
#define Nj(N) jinvloc(N,Ni(N))
#define Nk(N) kinvloc(N,Ni(N),Nj(N))
//N is global, n and m are local so the same mx,my,mz must not be used to compute i,j,k
#define nil(n) liinvloc(n)
#define njl(n) ljinvloc(n)
#define nkl(n) lkinvloc(n)

#define nis(n) sliinvloc(n)
#define njs(n) sljinvloc(n,nis(n))
#define nks(n) slkinvloc(n,nis(n),njs(n))

#define mil(m) liinvloc(m)
#define mjl(m) ljinvloc(m)
#define mkl(m) lkinvloc(m)

#define mis(m) sliinvloc(m)
#define mjs(m) sljinvloc(m,mis(m))
#define mks(m) slkinvloc(m,mis(m),mjs(m))

#define Mil(n,m,N) Ni(N) + mil(m) - nil(n)
#define Mjl(n,m,N) Nj(N) + mjl(m) - njl(n)
#define Mkl(n,m,N) Nk(N) + mkl(m) - nkl(n)

#define Mis(n,m,N) Ni(N) + mis(m) - nis(n)
#define Mjs(n,m,N) Nj(N) + mjs(m) - njs(n)
#define Mks(n,m,N) Nk(N) + mks(m) - nks(n)

#define Mref1(n,m,N) loc(Mil(n,m,N),Mjl(n,m,N),Mkl(n,m,N))
#define Mref2(n,m,N) loc(Mis(n,m,N),Mjs(n,m,N),Mks(n,m,N))

#endif //__COMMON_H_
