///////////////////////////////////////////////////////////////////////////////////////
// PCvar.h
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2014 University of Pennsylvania. All rights reserved.
// See http://www.cbica.upenn.edu/sbia/software/license.html or COYPING file.
//
// Contact: SBIA Group <sbia-software at uphs.upenn.edu>
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _PCVAR_H_
#define _PCVAR_H_

/*Context variables : globals*/

#include "PCShellFiles/pcShellHeader.h"
#include "RPFiles/rpHeader.h"

extern stsDMMG *dmmg;
#if (PETSC_VERSION_MAJOR <= 2) || ((PETSC_VERSION_MAJOR == 3) && (PETSC_VERSION_MINOR <= 1))
extern DA da;
#else
extern DM da;
#endif
extern int ierr, ilevel, mglev;

extern PCShellCtx pcShell;
extern char *pcShellName;
extern PetscInt pcShellOpt;

#endif
