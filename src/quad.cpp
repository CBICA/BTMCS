///////////////////////////////////////////////////////////////////////////////////////
// quad.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include "function.h"

void quad(double **xq, double **sq, double ***sh, double *det , double *wq,double Hxg, double Hyg, double Hzg)
{ 
	int iq;
	double detinv;
	double xr11, xr12, xr13, xr21, xr22, xr23, xr31, xr32, xr33;
	double x11, x12, x13, x21, x22, x23, x31, x32, x33;
	double x14, x15, x16, x24, x25, x26, x34, x35, x36;
	double x17, x18,  x27, x28, x37, x38;
	double sq11, sq12, sq13, sq14, sq15, sq16, sq17, sq18;
	double sq21, sq22, sq23, sq24, sq25, sq26, sq27, sq28;
	double sq31, sq32, sq33, sq34, sq35, sq36, sq37, sq38;
	double cf11, cf12, cf13, cf21, cf22, cf23, cf31, cf32, cf33;

	if (nsd == 2) {
		xq[0][0] = - 0.577350269189626;
		xq[0][1] = - 0.577350269189626;
		xq[1][0] = + 0.577350269189626;
		xq[1][1] = - 0.577350269189626;
		xq[2][0] = - 0.577350269189626;
		xq[2][1] = + 0.577350269189626;
		xq[3][0] = + 0.577350269189626;
		xq[3][1] = + 0.577350269189626;
		wq[0] = 1.0;
		wq[1] = 1.0;
		wq[2] = 1.0;
		wq[3] = 1.0;
		x11 = 0; x12 = Hxg; x13 = Hxg; x14 = 0.0;
		x21 = 0; x22 = 0; x23 = Hyg; x24 = Hyg;
		for (iq=0;iq<nquad;iq++){

			sq[iq][0] = 0.25 * (1 - xq[iq][0]) * (1 - xq[iq][1]);
			sq[iq][1] = 0.25 * (1 + xq[iq][0]) * (1 - xq[iq][1]);
			sq[iq][2] = 0.25 * (1 + xq[iq][0]) * (1 + xq[iq][1]);
			sq[iq][3] = 0.25 * (1 - xq[iq][0]) * (1 + xq[iq][1]);

			sq11 = - 0.25 * (1 - xq[iq][1]);
			sq12 = + 0.25 * (1 - xq[iq][1]);
			sq13 = + 0.25 * (1 + xq[iq][1]);
			sq14 = - 0.25 * (1 + xq[iq][1]);

			sq21 = - 0.25 * (1 - xq[iq][0]);
			sq22 = - 0.25 * (1 + xq[iq][0]);
			sq23 = + 0.25 * (1 + xq[iq][0]);
			sq24 = + 0.25 * (1 - xq[iq][0]);


			xr11 = sq11 * x11 + sq12 * x12 + sq13 * x13 + sq14 * x14;
			xr12 = sq11 * x21 + sq12 * x22 + sq13 * x23 + sq14 * x24;
			xr21 = sq21 * x11 + sq22 * x12 + sq23 * x13 + sq24 * x14;
			xr22 = sq21 * x21 + sq22 * x22 + sq23 * x23 + sq24 * x24;

			cf11 = + xr22;
			cf12 = - xr12;
			detinv = xr11 * cf11 + xr21 * cf12;
			det[iq] = detinv;
			detinv = 1.0 / detinv;

			cf21 = - xr21;
			cf22 = + xr11;

			cf11 = cf11 * detinv;
			cf12 = cf12 * detinv;                         
			cf21 = cf21 * detinv;                         
			cf22 = cf22 * detinv;                         

			sh[iq][0][0] = sq11 * cf11 + sq21 * cf12;
			sh[iq][0][1] = sq11 * cf21 + sq21 * cf22;
			sh[iq][1][0] = sq12 * cf11 + sq22 * cf12;
			sh[iq][1][1] = sq12 * cf21 + sq22 * cf22;
			sh[iq][2][0] = sq13 * cf11 + sq23 * cf12;
			sh[iq][2][1] = sq13 * cf21 + sq23 * cf22;
			sh[iq][3][0] = sq14 * cf11 + sq24 * cf12;
			sh[iq][3][1] = sq14 * cf21 + sq24 * cf22;
		}

	} else if (nsd == 3) {
		xq[0][0] = - 0.577350269189626;
		xq[0][1] = - 0.577350269189626;
		xq[0][2] = - 0.577350269189626;
		xq[1][0] = + 0.577350269189626;
		xq[1][1] = - 0.577350269189626;
		xq[1][2] = - 0.577350269189626;
		xq[2][0] = - 0.577350269189626;
		xq[2][1] = + 0.577350269189626;
		xq[2][2] = - 0.577350269189626;
		xq[3][0] = + 0.577350269189626;
		xq[3][1] = + 0.577350269189626;
		xq[3][2] = - 0.577350269189626;
		xq[4][0] = - 0.577350269189626;
		xq[4][1] = - 0.577350269189626;
		xq[4][2] = + 0.577350269189626;
		xq[5][0] = + 0.577350269189626;
		xq[5][1] = - 0.577350269189626;
		xq[5][2] = + 0.577350269189626;
		xq[6][0] = - 0.577350269189626;
		xq[6][1] = + 0.577350269189626;
		xq[6][2] = + 0.577350269189626;
		xq[7][0] = + 0.577350269189626;
		xq[7][1] = + 0.577350269189626;
		xq[7][2] = + 0.577350269189626;

		wq[0] = 1.0; wq[1] = 1.0; wq[2] = 1.0; wq[3] = 1.0; wq[4] = 1.0; wq[5] = 1.0;
		wq[6] = 1.0; wq[7] = 1.0;

		x11 = 0; x12 = Hxg; x13 = Hxg; x14 = 0.0; x15 = 0.0; x16 = Hxg; x17 = Hxg; x18 = 0;
		x21 = 0; x22 = 0; x23 = Hyg; x24 = Hyg; x25 = 0; x26 = 0; x27 = Hyg; x28 = Hyg; 
		x31 = 0; x32 = 0; x33 = 0; x34 = 0; x35 = Hzg; x36 = Hzg; x37 = Hzg; x38 = Hzg;
		for (iq=0;iq<nquad;iq++){
			sq[iq][0] = 0.125 * (1 - xq[iq][0]) * (1 - xq[iq][1]) * (1 - xq[iq][2]);
			sq[iq][1] = 0.125 * (1 + xq[iq][0]) * (1 - xq[iq][1]) * (1 - xq[iq][2]);
			sq[iq][2] = 0.125 * (1 + xq[iq][0]) * (1 + xq[iq][1]) * (1 - xq[iq][2]);
			sq[iq][3] = 0.125 * (1 - xq[iq][0]) * (1 + xq[iq][1]) * (1 - xq[iq][2]);
			sq[iq][4] = 0.125 * (1 - xq[iq][0]) * (1 - xq[iq][1]) * (1 + xq[iq][2]);
			sq[iq][5] = 0.125 * (1 + xq[iq][0]) * (1 - xq[iq][1]) * (1 + xq[iq][2]);
			sq[iq][6] = 0.125 * (1 + xq[iq][0]) * (1 + xq[iq][1]) * (1 + xq[iq][2]);
			sq[iq][7] = 0.125 * (1 - xq[iq][0]) * (1 + xq[iq][1]) * (1 + xq[iq][2]);

			sq11 = - 0.125 * (1 - xq[iq][1]) * (1 - xq[iq][2]);
			sq12 = + 0.125 * (1 - xq[iq][1]) * (1 - xq[iq][2]);
			sq13 = + 0.125 * (1 + xq[iq][1]) * (1 - xq[iq][2]);
			sq14 = - 0.125 * (1 + xq[iq][1]) * (1 - xq[iq][2]);
			sq15 = - 0.125 * (1 - xq[iq][1]) * (1 + xq[iq][2]);
			sq16 = + 0.125 * (1 - xq[iq][1]) * (1 + xq[iq][2]);
			sq17 = + 0.125 * (1 + xq[iq][1]) * (1 + xq[iq][2]);
			sq18 = - 0.125 * (1 + xq[iq][1]) * (1 + xq[iq][2]);

			sq21 = - 0.125 * (1 - xq[iq][0]) * (1 - xq[iq][2]);
			sq22 = - 0.125 * (1 + xq[iq][0]) * (1 - xq[iq][2]);
			sq23 = + 0.125 * (1 + xq[iq][0]) * (1 - xq[iq][2]);
			sq24 = + 0.125 * (1 - xq[iq][0]) * (1 - xq[iq][2]);
			sq25 = - 0.125 * (1 - xq[iq][0]) * (1 + xq[iq][2]);
			sq26 = - 0.125 * (1 + xq[iq][0]) * (1 + xq[iq][2]);
			sq27 = + 0.125 * (1 + xq[iq][0]) * (1 + xq[iq][2]);
			sq28 = + 0.125 * (1 - xq[iq][0]) * (1 + xq[iq][2]);

			sq31 = - 0.125 * (1 - xq[iq][0]) * (1 - xq[iq][1]);
			sq32 = - 0.125 * (1 + xq[iq][0]) * (1 - xq[iq][1]);
			sq33 = - 0.125 * (1 + xq[iq][0]) * (1 + xq[iq][1]);
			sq34 = - 0.125 * (1 - xq[iq][0]) * (1 + xq[iq][1]);
			sq35 = + 0.125 * (1 - xq[iq][0]) * (1 - xq[iq][1]);
			sq36 = + 0.125 * (1 + xq[iq][0]) * (1 - xq[iq][1]);
			sq37 = + 0.125 * (1 + xq[iq][0]) * (1 + xq[iq][1]);
			sq38 = + 0.125 * (1 - xq[iq][0]) * (1 + xq[iq][1]);

			xr11 = sq11 * x11 + sq12 * x12 + sq13 * x13 + sq14 * x14  + sq15 * x15 + sq16 * x16 + sq17 * x17 + sq18 * x18;
			xr12 = sq11 * x21 + sq12 * x22 + sq13 * x23 + sq14 * x24  + sq15 * x25 + sq16 * x26 + sq17 * x27 + sq18 * x28;
			xr13 = sq11 * x31 + sq12 * x32 + sq13 * x33 + sq14 * x34  + sq15 * x35 + sq16 * x36 + sq17 * x37 + sq18 * x38;
			xr21 = sq21 * x11 + sq22 * x12 + sq23 * x13 + sq24 * x14  + sq25 * x15 + sq26 * x16 + sq27 * x17 + sq28 * x18;
			xr22 = sq21 * x21 + sq22 * x22 + sq23 * x23 + sq24 * x24  + sq25 * x25 + sq26 * x26 + sq27 * x27 + sq28 * x28;
			xr23 = sq21 * x31 + sq22 * x32 + sq23 * x33 + sq24 * x34  + sq25 * x35 + sq26 * x36 + sq27 * x37 + sq28 * x38;
			xr31 = sq31 * x11 + sq32 * x12 + sq33 * x13 + sq34 * x14  + sq35 * x15 + sq36 * x16 + sq37 * x17 + sq38 * x18;
			xr32 = sq31 * x21 + sq32 * x22 + sq33 * x23 + sq34 * x24  + sq35 * x25 + sq36 * x26 + sq37 * x27 + sq38 * x28;
			xr33 = sq31 * x31 + sq32 * x32 + sq33 * x33 + sq34 * x34  + sq35 * x35 + sq36 * x36 + sq37 * x37 + sq38 * x38;


			cf11 = + (xr22 * xr33 - xr32 * xr23);
			cf12 = - (xr12 * xr33 - xr32 * xr13);
			cf13 = + (xr12 * xr23 - xr22 * xr13);

			detinv = xr11 * cf11 + xr21 * cf12 + xr31 * cf13;
			det[iq] = detinv;
			detinv = 1.0 / detinv;

			cf21 = - (xr21 * xr33 - xr31 * xr23);
			cf22 = + (xr11 * xr33 - xr31 * xr13);
			cf23 = - (xr11 * xr23 - xr21 * xr13);
			cf31 = + (xr21 * xr32 - xr31 * xr22);
			cf32 = - (xr11 * xr32 - xr31 * xr12);
			cf33 = + (xr11 * xr22 - xr21 * xr12);

			cf11 = cf11 * detinv;
			cf12 = cf12 * detinv;
			cf13 = cf13 * detinv;
			cf21 = cf21 * detinv;
			cf22 = cf22 * detinv;
			cf23 = cf23 * detinv;
			cf31 = cf31 * detinv;
			cf32 = cf32 * detinv;
			cf33 = cf33 * detinv;


			sh[iq][0][0] = sq11 * cf11 + sq21 * cf12 + sq31 * cf13;
			sh[iq][0][1] = sq11 * cf21 + sq21 * cf22 + sq31 * cf23;
			sh[iq][0][2] = sq11 * cf31 + sq21 * cf32 + sq31 * cf33;
			sh[iq][1][0] = sq12 * cf11 + sq22 * cf12 + sq32 * cf13;
			sh[iq][1][1] = sq12 * cf21 + sq22 * cf22 + sq32 * cf23;
			sh[iq][1][2] = sq12 * cf31 + sq22 * cf32 + sq32 * cf33;
			sh[iq][2][0] = sq13 * cf11 + sq23 * cf12 + sq33 * cf13;
			sh[iq][2][1] = sq13 * cf21 + sq23 * cf22 + sq33 * cf23;
			sh[iq][2][2] = sq13 * cf31 + sq23 * cf32 + sq33 * cf33;
			sh[iq][3][0] = sq14 * cf11 + sq24 * cf12 + sq34 * cf13;
			sh[iq][3][1] = sq14 * cf21 + sq24 * cf22 + sq34 * cf23;
			sh[iq][3][2] = sq14 * cf31 + sq24 * cf32 + sq34 * cf33;
			sh[iq][4][0] = sq15 * cf11 + sq25 * cf12 + sq35 * cf13;
			sh[iq][4][1] = sq15 * cf21 + sq25 * cf22 + sq35 * cf23;
			sh[iq][4][2] = sq15 * cf31 + sq25 * cf32 + sq35 * cf33;
			sh[iq][5][0] = sq16 * cf11 + sq26 * cf12 + sq36 * cf13;
			sh[iq][5][1] = sq16 * cf21 + sq26 * cf22 + sq36 * cf23;
			sh[iq][5][2] = sq16 * cf31 + sq26 * cf32 + sq36 * cf33;
			sh[iq][6][0] = sq17 * cf11 + sq27 * cf12 + sq37 * cf13;
			sh[iq][6][1] = sq17 * cf21 + sq27 * cf22 + sq37 * cf23;
			sh[iq][6][2] = sq17 * cf31 + sq27 * cf32 + sq37 * cf33;
			sh[iq][7][0] = sq18 * cf11 + sq28 * cf12 + sq38 * cf13;
			sh[iq][7][1] = sq18 * cf21 + sq28 * cf22 + sq38 * cf23;
			sh[iq][7][2] = sq18 * cf31 + sq28 * cf32 + sq38 * cf33;
		}
	}
}

void quad1(double **xq, double **sq, double ***sh, double *det , double *wq,double Hxg, double Hyg, double Hzg)
{ 
	int iq;
	double detinv;
	double xr11, xr12, xr13, xr21, xr22, xr23, xr31, xr32, xr33;

	double x11, x12, x13, x21, x22, x23, x31, x32, x33;
	double x14, x15, x16, x24, x25, x26, x34, x35, x36;
	double x17, x18,  x27, x28, x37, x38;
	double sq11, sq12, sq13, sq14, sq15, sq16, sq17, sq18;
	double sq21, sq22, sq23, sq24, sq25, sq26, sq27, sq28;
	double sq31, sq32, sq33, sq34, sq35, sq36, sq37, sq38;
	double cf11, cf12, cf13, cf21, cf22, cf23, cf31, cf32, cf33;

	if (nsd == 3) {
		xq[0][0] =  0.0;
		xq[0][1] =  0.0;
		xq[0][2] =  0.0;

		wq[0] = 8.0; 

		x11 = 0; x12 = Hxg; x13 = Hxg; x14 = 0.0; x15 = 0.0; x16 = Hxg; x17 = Hxg; x18 = 0;
		x21 = 0; x22 = 0; x23 = Hyg; x24 = Hyg; x25 = 0; x26 = 0; x27 = Hyg; x28 = Hyg; 
		x31 = 0; x32 = 0; x33 = 0; x34 = 0; x35 = Hzg; x36 = Hzg; x37 = Hzg; x38 = Hzg;
		for (iq=0;iq<rednquad;iq++){
			sq[iq][0] = 0.125 * (1 - xq[iq][0]) * (1 - xq[iq][1]) * (1 - xq[iq][2]);
			sq[iq][1] = 0.125 * (1 + xq[iq][0]) * (1 - xq[iq][1]) * (1 - xq[iq][2]);
			sq[iq][2] = 0.125 * (1 + xq[iq][0]) * (1 + xq[iq][1]) * (1 - xq[iq][2]);
			sq[iq][3] = 0.125 * (1 - xq[iq][0]) * (1 + xq[iq][1]) * (1 - xq[iq][2]);
			sq[iq][4] = 0.125 * (1 - xq[iq][0]) * (1 - xq[iq][1]) * (1 + xq[iq][2]);
			sq[iq][5] = 0.125 * (1 + xq[iq][0]) * (1 - xq[iq][1]) * (1 + xq[iq][2]);
			sq[iq][6] = 0.125 * (1 + xq[iq][0]) * (1 + xq[iq][1]) * (1 + xq[iq][2]);
			sq[iq][7] = 0.125 * (1 - xq[iq][0]) * (1 + xq[iq][1]) * (1 + xq[iq][2]);

			sq11 = - 0.125 * (1 - xq[iq][1]) * (1 - xq[iq][2]);
			sq12 = + 0.125 * (1 - xq[iq][1]) * (1 - xq[iq][2]);
			sq13 = + 0.125 * (1 + xq[iq][1]) * (1 - xq[iq][2]);
			sq14 = - 0.125 * (1 + xq[iq][1]) * (1 - xq[iq][2]);
			sq15 = - 0.125 * (1 - xq[iq][1]) * (1 + xq[iq][2]);
			sq16 = + 0.125 * (1 - xq[iq][1]) * (1 + xq[iq][2]);
			sq17 = + 0.125 * (1 + xq[iq][1]) * (1 + xq[iq][2]);
			sq18 = - 0.125 * (1 + xq[iq][1]) * (1 + xq[iq][2]);

			sq21 = - 0.125 * (1 - xq[iq][0]) * (1 - xq[iq][2]);
			sq22 = - 0.125 * (1 + xq[iq][0]) * (1 - xq[iq][2]);
			sq23 = + 0.125 * (1 + xq[iq][0]) * (1 - xq[iq][2]);
			sq24 = + 0.125 * (1 - xq[iq][0]) * (1 - xq[iq][2]);
			sq25 = - 0.125 * (1 - xq[iq][0]) * (1 + xq[iq][2]);
			sq26 = - 0.125 * (1 + xq[iq][0]) * (1 + xq[iq][2]);
			sq27 = + 0.125 * (1 + xq[iq][0]) * (1 + xq[iq][2]);
			sq28 = + 0.125 * (1 - xq[iq][0]) * (1 + xq[iq][2]);

			sq31 = - 0.125 * (1 - xq[iq][0]) * (1 - xq[iq][1]);
			sq32 = - 0.125 * (1 + xq[iq][0]) * (1 - xq[iq][1]);
			sq33 = - 0.125 * (1 + xq[iq][0]) * (1 + xq[iq][1]);
			sq34 = - 0.125 * (1 - xq[iq][0]) * (1 + xq[iq][1]);
			sq35 = + 0.125 * (1 - xq[iq][0]) * (1 - xq[iq][1]);
			sq36 = + 0.125 * (1 + xq[iq][0]) * (1 - xq[iq][1]);
			sq37 = + 0.125 * (1 + xq[iq][0]) * (1 + xq[iq][1]);
			sq38 = + 0.125 * (1 - xq[iq][0]) * (1 + xq[iq][1]);

			xr11 = sq11 * x11 + sq12 * x12 + sq13 * x13 + sq14 * x14  + sq15 * x15 + sq16 * x16 + sq17 * x17 + sq18 * x18;
			xr12 = sq11 * x21 + sq12 * x22 + sq13 * x23 + sq14 * x24  + sq15 * x25 + sq16 * x26 + sq17 * x27 + sq18 * x28;
			xr13 = sq11 * x31 + sq12 * x32 + sq13 * x33 + sq14 * x34  + sq15 * x35 + sq16 * x36 + sq17 * x37 + sq18 * x38;
			xr21 = sq21 * x11 + sq22 * x12 + sq23 * x13 + sq24 * x14  + sq25 * x15 + sq26 * x16 + sq27 * x17 + sq28 * x18;
			xr22 = sq21 * x21 + sq22 * x22 + sq23 * x23 + sq24 * x24  + sq25 * x25 + sq26 * x26 + sq27 * x27 + sq28 * x28;
			xr23 = sq21 * x31 + sq22 * x32 + sq23 * x33 + sq24 * x34  + sq25 * x35 + sq26 * x36 + sq27 * x37 + sq28 * x38;
			xr31 = sq31 * x11 + sq32 * x12 + sq33 * x13 + sq34 * x14  + sq35 * x15 + sq36 * x16 + sq37 * x17 + sq38 * x18;
			xr32 = sq31 * x21 + sq32 * x22 + sq33 * x23 + sq34 * x24  + sq35 * x25 + sq36 * x26 + sq37 * x27 + sq38 * x28;
			xr33 = sq31 * x31 + sq32 * x32 + sq33 * x33 + sq34 * x34  + sq35 * x35 + sq36 * x36 + sq37 * x37 + sq38 * x38;


			cf11 = + (xr22 * xr33 - xr32 * xr23);
			cf12 = - (xr12 * xr33 - xr32 * xr13);
			cf13 = + (xr12 * xr23 - xr22 * xr13);

			detinv = xr11 * cf11 + xr21 * cf12 + xr31 * cf13;
			det[iq] = detinv;
			detinv = 1.0 / detinv;

			cf21 = - (xr21 * xr33 - xr31 * xr23);
			cf22 = + (xr11 * xr33 - xr31 * xr13);
			cf23 = - (xr11 * xr23 - xr21 * xr13);
			cf31 = + (xr21 * xr32 - xr31 * xr22);
			cf32 = - (xr11 * xr32 - xr31 * xr12);
			cf33 = + (xr11 * xr22 - xr21 * xr12);

			cf11 = cf11 * detinv;
			cf12 = cf12 * detinv;
			cf13 = cf13 * detinv;
			cf21 = cf21 * detinv;
			cf22 = cf22 * detinv;
			cf23 = cf23 * detinv;
			cf31 = cf31 * detinv;
			cf32 = cf32 * detinv;
			cf33 = cf33 * detinv;

			sh[iq][0][0] = sq11 * cf11 + sq21 * cf12 + sq31 * cf13;
			sh[iq][0][1] = sq11 * cf21 + sq21 * cf22 + sq31 * cf23;
			sh[iq][0][2] = sq11 * cf31 + sq21 * cf32 + sq31 * cf33;
			sh[iq][1][0] = sq12 * cf11 + sq22 * cf12 + sq32 * cf13;
			sh[iq][1][1] = sq12 * cf21 + sq22 * cf22 + sq32 * cf23;
			sh[iq][1][2] = sq12 * cf31 + sq22 * cf32 + sq32 * cf33;
			sh[iq][2][0] = sq13 * cf11 + sq23 * cf12 + sq33 * cf13;
			sh[iq][2][1] = sq13 * cf21 + sq23 * cf22 + sq33 * cf23;
			sh[iq][2][2] = sq13 * cf31 + sq23 * cf32 + sq33 * cf33;
			sh[iq][3][0] = sq14 * cf11 + sq24 * cf12 + sq34 * cf13;
			sh[iq][3][1] = sq14 * cf21 + sq24 * cf22 + sq34 * cf23;
			sh[iq][3][2] = sq14 * cf31 + sq24 * cf32 + sq34 * cf33;
			sh[iq][4][0] = sq15 * cf11 + sq25 * cf12 + sq35 * cf13;
			sh[iq][4][1] = sq15 * cf21 + sq25 * cf22 + sq35 * cf23;
			sh[iq][4][2] = sq15 * cf31 + sq25 * cf32 + sq35 * cf33;
			sh[iq][5][0] = sq16 * cf11 + sq26 * cf12 + sq36 * cf13;
			sh[iq][5][1] = sq16 * cf21 + sq26 * cf22 + sq36 * cf23;
			sh[iq][5][2] = sq16 * cf31 + sq26 * cf32 + sq36 * cf33;
			sh[iq][6][0] = sq17 * cf11 + sq27 * cf12 + sq37 * cf13;
			sh[iq][6][1] = sq17 * cf21 + sq27 * cf22 + sq37 * cf23;
			sh[iq][6][2] = sq17 * cf31 + sq27 * cf32 + sq37 * cf33;
			sh[iq][7][0] = sq18 * cf11 + sq28 * cf12 + sq38 * cf13;
			sh[iq][7][1] = sq18 * cf21 + sq28 * cf22 + sq38 * cf23;
			sh[iq][7][2] = sq18 * cf31 + sq28 * cf32 + sq38 * cf33;
		}

	}
}

void analnode(double **sh) 
{ 
	if (nsd == 2) {
		sh[0][0] = -1;
		sh[0][1] = -1;
		sh[1][0] = 1;
		sh[1][1] = -1;
		sh[2][0] = 1;
		sh[2][1] = 1;
		sh[3][0] = -1;
		sh[3][1] = 1;

	} else if (nsd == 3) {
		sh[0][0] = -1;
		sh[0][1] = -1;
		sh[0][2] = -1;


		sh[1][0] = 1;
		sh[1][1] = -1;
		sh[1][2] = -1;

		sh[2][0] = 1;
		sh[2][1] = 1;
		sh[2][2] = -1;

		sh[3][0] = -1;
		sh[3][1] = 1;
		sh[3][2] = -1;

		sh[4][0] = -1;
		sh[4][1] = -1;
		sh[4][2] = 1;

		sh[5][0] = 1;
		sh[5][1] = -1;
		sh[5][2] = 1;

		sh[6][0] = 1;
		sh[6][1] = 1;
		sh[6][2] = 1;

		sh[7][0] = -1;
		sh[7][1] = 1;
		sh[7][2] = 1;
	}
}
