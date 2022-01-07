/*
 * $Id: qmmmrec.h,v 1.1.4.1 2008/10/07 22:49:40 lindahl Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct {
  int           nrQMatoms;      /* total nr of QM atoms              */
  rvec          *xQM;           /* shifted to center of box          */  
  rvec          *rQM;           /* for completing molecules, when using clustering */
  int           *indexQM;       /* atom i = atom indexQM[i] in mdrun */
  int           *atomicnumberQM;/* atomic numbers of QM atoms        */  
  real          *QMcharges;     /* atomic charges of QM atoms(ONIOM) */
  int           *shiftQM;
  int           QMcharge;       /* charge of the QM system           */
  int           multiplicity;   /* multipicity (no of unpaired eln)  */
  int           QMfunctional;   /* see enums.h for all functionals       */
  int           QMengine;       /* see enums.h for all engines       */
  int           QMmethod;       /* see enums.h for all methods       */
  int           QMbasis;        /* see enums.h for all bases         */
  int           nelectrons;     /* total number of elecs in QM region*/
  bool          bTS;            /* Optimize a TS, only steep, no md  */
  bool          bOPT;           /* Optimize QM subsys, only steep, no md  */
  bool          *frontatoms;    /* qm atoms on the QM side of a QM-MM bond */
  real          *pcharges;      /* Some partial charges for the system */
  int           *useQMatom;     /* can be used to remove dummy atoms in i.e. tip4p for QM */
  real          *QMcA;          /* MATRIX */
  int           *QMcB;          /* MATRIX */
  int           *QMcV;          /* VECTOR */
  int           *QMcIS;         /* MATRIX */
  int           ePBC;
  matrix        box;
  int           dI;             /* FOR RERUNS */
  int           includeQ;
  int           guessWF;        /* Guess the wavefunction with STO-3G */

  /* Gaussian specific stuff */
  int           nQMthread;
  int           nQMproc;
  int           nQMcpus;        /* no. of CPUs used for the QM calc. MPQC::PER NODE */
  int           QMmem;          /* memory for the gaussian calc.     */
  int           nQMnodes;	/* MPQC::DISTRIBUTE ON THIS AMOUNT OF NODES */
  int           accuracy;       /* convergence criterium (E(-x))     */
  bool          cpmcscf;        /* using cpmcscf(l1003)*/
  char          *gauss_dir;
  char          *gauss_exe;
  char          *devel_dir;
  real          *c6;
  real          *c12;
  /* Surface hopping stuff */
  bool          bSH;            /* surface hopping (diabatic only)   */
  real          SAon;           /* at which energy gap the SA starts */
  real          SAoff;          /* at which energy gap the SA stops  */
  int           SAsteps;        /* stepwise switchinng on the SA     */
  int           SAstep;         /* current state of SA               */
  int           CIdim;
  real          *CIvec1;
  real          *CIvec2;
  real          *CIvec1old;
  real          *CIvec2old;
  ivec          SHbasis;
  int           CASelectrons;
  int           CASorbitals;
  /* QM clustering */
  int           nQMclusterupdate; /* at which steps should we update. 1 implies every step. 0 is only the first */
  bool          bCLUSTER;         /* Do clustering */
  int           QMstep;
  int		nrQMc;            /* Number of QM clusters */
  int           *QMcluster;       /* cluster index values */
  rvec          *QMcluspos;       /* cluster center position */
  int           *QMexec;
  real          *CLR;             /* cluster radius (derived from QM atom positions) */
  real          level_shift;      // The quantum mechanical level shift ( yes: only applies for QM atoms )
  bool          bQMrestart;
  bool          bQMcheckpoint;
  bool          bQMsavestate;
  bool          bQMdispersion;    // add dispersion on Q atoms
  bool          bQMexclusion;     // exclude double counted electrostatics
  real          fQMexclusion;     // screen with this constant
  bool          bPrint_npa;       // Print normal population analysis to file npa.log
  char          mpqcappend[5];
  char          *mpihandler;
  real          Erep[11];
  int           Nrep[2];
  int  nQMcenter;
  real fQMcutoff;
  real fQMelrc;
  // S-F*exp(-K*(r-P)^Pow)
  real fQMscreenK;
  real fQMscreenP;
  real fQMscreenS;
  real fQMscreenF;
  real fQMscreenPow;
  real fMMgradMask;
  real fQMgradMask;
  bool bQMadh;                  //adhoc scaling
} t_QMrec;

typedef struct {
  int           nrMMatoms;      /* nr of MM atoms, updated every step*/
  rvec          *xMM;           /* shifted to center of box          */
  rvec          *rMM;           /* complete molecules in clusters    */
  int           *indexMM;       /* atom i = atom indexMM[I] in mdrun */
  int           *index;         /* static index */
  int           ntotatoms;      // book keep
  int           ntMM;
  real          *MMcharges;     /* MM point charges in std QMMM calc.*/
  int           *shiftMM;
  int           *MMatomtype;    /* only important for semi-emp.      */
  real          scalefactor;
  /* gaussian specific stuff */
  real          *c6;
  real          *c12;
  /* MM clustering              some not used since they are not decoupled from qm clusters*/
  int           *MMcIS;
  int           nrMMc;		/* Not used */
  int           nrQa;           /* Number of pointcharges derived from QM charges belonging to this cluster */
  int           *MMcluster;     /* Indices for the cluster */
  rvec          *MMcluspos;     /* Positions of the cluster (not used)*/
  real          *fcharges;      /* scaling constant for the charges    */
} t_MMrec;

typedef struct {
  int           QMMMscheme; /* ONIOM (multi-layer) or normal          */
  int           nrQMlayers; /* number of QM layers (total layers +1 (MM)) */
  t_QMrec       **qm;        /* atoms and run params for each QM group */
  t_MMrec       *mm;        /* there can only be one MM subsystem !   */
} t_QMMMrec;


