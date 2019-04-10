///////////////////////////////////////////////////////////////////////////////////////
// ConservLaw.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef CONSERVLAW_H
#define CONSERVLAW_H

/* {{{ CONSERVATION LAW SOLVER - NO REACTION*/
/* 
	struct to implement transport (advection/diffusion) equations 
	da: the distributed array
	m0: initial conditions for transported field
	v:  velocity 
	s: right-hand side
	dt: time step
   	cfln : the CFL number for stability (e.g. 1/2)
	(the number of subsequent time sub-steps to be taken accordingly is
	being computed inside each transport function (*transport))
   	m: updated field

	update implements both gradAdvection, and DivAdvection, depending

*/

enum LawType {useDivAdvection, useGradAdvection, useDiffusion, noDiffusion};

struct ConservationLaw 
{
	//DA setupDA;

	enum LawType advt; // set advection type
	enum LawType difft; // set diffusion type

	//int (*setDA)(struct ConservationLaw *cl, DA da); // creates stencils and sets setupDA value
	int (*setGradAdvectionOn)(struct ConservationLaw *cl);
	//int (*setGradAdvectionOff)(struct ConservationLaw *cl);
	int (*setDivAdvectionOn)(struct ConservationLaw *cl);
	//int (*setDivAdvectionOff)(struct ConservationLaw *cl);
	int (*setDiffusionOn)(struct ConservationLaw *cl);

	void (*transport)(struct ConservationLaw *cl, Vec m0, Vec v, Vec s, double dt, double cfln,  Vec m);
	/* make sure da similar to setupDA */
	/* the vectors should be block vectors, so one can get the block size and work on the advection. */
};

extern void InitConservationLaw(struct ConservationLaw *cl);
extern void AdvectionGrad(struct ConservationLaw *cl, Vec m0, Vec v, Vec s, double dt, double cfln,  Vec m);
extern void AdvectionDiv(struct ConservationLaw *cl, Vec m0, Vec v, Vec s, double dt, double cfln,  Vec m);
extern void DiffusionExplicit(struct ConservationLaw *cl, Vec c0, Vec D, Vec s, double dt, double cfln,  Vec c);

#endif
