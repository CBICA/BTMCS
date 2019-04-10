
  Section of Biomedical Image Analysis
  Department of Radiology
  University of Pennsylvania
  3600 Market Street, Suite 380
  Philadelphia, PA 19104

  Web:   http://www.cbica.upenn.edu/sbia/
  Email: sbia-software at uphs.upenn.edu

  Copyright (c) 2014 University of Pennsylvania. All rights reserved.
  See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.



INTRODUCTION
============

  Brain Tumor Modeling - Coupled Solver (BTMCS) [1] is a software package for tumor growth 
  modeling employs a DiffusionSolver approach [2,3,4]. It is a one-step further compared 
  to the previous purely mechanical, pressure-based approach employed in ElasticSolver and 
  described in [PMB2007]_. In addition to the previous elasticity-based approach to simulate 
  brain tissue deformations caused by growing tumors (mass-effect), a nonlinear reaction-
  advection-diffusion equation to describe the tumor spatio-temporal evolution is added. 
  This equation has a two-way coupling with the underlying tissue elastic deformation. In 
  this approach, the forces exerted by the tumor growth and infiltration onto the underlying 
  brain parenchyma are local ones, proportional to local tumor density gradients. 
  
  Curretnly BTMCS is used in GLioma Image SegmenTation and Registration (GLISTR) [5,6] and 
  Pre-Operative and post-Recurrence brain Tumor Registration (PORTR) [7,8].
  


PACKAGE OVERVIEW
================

  - config/             Package configuration files.
  - doc/                Software documentation such as the software manual.
  - src/                Main source code files.
  - testing/            Code and Data files for testing software.
  - AUTHORS.txt         A list of the people who contributed to this software.
  - ChangeLog.txt       A log of changes between versions.
  - CMakeLists.txt      Root CMake configuration file.
  - COPYING.txt         The copyright and license notices.
  - INSTALL.txt         Build and installation instructions.
  - README.txt          This readme file.



DOCUMENTATION
=============

  See the software manual for details on the software including a demonstration
  of how to apply the software tools provided by this package.



INSTALLATION
============

  See http://www.cbica.upenn.edu/sbia/software/tumorsimulator/installation.html or
  the installation section of the software manual.



LICENSING
=========

  See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.



REFERENCES
==========

  [1] http://www.cbica.upenn.edu/sbia/software/tumorsimulator/
  
  [2] C. Hogea, C. Davatzikos and G. Biros, "An image-driven parameter estimation
      problem for a glioma growth model with mass effects", J. Math Biol., 56(6): 
      793-825 (2008).

  [3] C. Hogea, C. Davatzikos and G. Biros, "Brain-Tumor Interaction Biophysical
      Models for Medical Image Registration", SIAM J. Scientific Computing 30(6): 
      3050-3072 (2008).

  [4] C. Hogea, C. Davatzikos and G. Biros, "Modeling Glioma Growth and Mass Effect 
      in 3D MR Images of the Brain", In Proc. MICCAI (1): 642-650 (2007).

  [5] http://www.cbica.upenn.edu/sbia/software/glistr/

  [6] D. Kwon, R.T. Shinohara, H. Akbari, C. Davatzikos, "Combining Generative Models for 
      Multifocal Glioma Segmentation and Registration", In: Proc. MICCAI (1): 763-770 (2014) 

  [7] http://www.cbica.upenn.edu/sbia/software/portr/

  [8] D. Kwon, M. Niethammer, H. Akbari, M. Bilello, C. Davatzikos, and K.M. Pohl, 
      "PORTR: Pre-Operative and Post-Recurrence Brain Tumor Registration", 
      IEEE Trans. Med. Imaging 33(3): 651-667 (2014)
