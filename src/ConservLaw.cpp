///////////////////////////////////////////////////////////////////////////////////////
// ConservLaw.c
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#include "global.h"
#include "ConservLaw.h"


int SetGradAdvectionOn(struct ConservationLaw *cl) 
{
	cl->advt = useGradAdvection;
	cl->transport = &AdvectionGrad;
	return 0;
}

int SetDivAdvectionOn(struct ConservationLaw *cl) 
{
	cl->advt = useDivAdvection;
	cl->transport = &AdvectionDiv;
	return 0;
}

int SetDiffusionOn(struct ConservationLaw *cl) 
{
	cl->difft = useDiffusion;
	cl->transport = &DiffusionExplicit;
	return 0;
}

// initialize function pointer(s) to the appropriate functions
void InitConservationLaw(struct ConservationLaw * cl) 
{
	cl->setGradAdvectionOn = &SetGradAdvectionOn;
	cl->setDivAdvectionOn = &SetDivAdvectionOn;
	cl->setDiffusionOn = &SetDiffusionOn;

	// same for the rest of set
	// bring context (DA, at df to some default status
}
