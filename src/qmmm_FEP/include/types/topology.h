/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

enum {
  egcTC,    egcENER,   egcACC, egcFREEZE, 
  egcUser1, egcUser2,  egcVCM, egcXTC,
  egcORFIT, egcQMMM, egcQMskip,
  egcNR 
};

typedef struct {
  char          **name;         /* Name of the molecule type  	        */
  t_atoms	atoms;		/* The atoms		       	        */
  t_ilist       ilist[F_NRE];
  t_block       cgs;            /* The charge groups                    */
  t_blocka      excls;          /* The exclusions                       */
} gmx_moltype_t;

typedef struct {
  int           type;           /* The molcule type index in mtop.moltype */
  int           nmol;           /* The number of molecules in this block  */
  int           natoms_mol;     /* The number of atoms in one molecule    */
  int           nposres_xA;     /* The number of posres coords for top A  */
  rvec          *posres_xA;     /* The posres coords for top A            */
  int           nposres_xB;     /* The number of posres coords for top B  */
  rvec          *posres_xB;     /* The posres coords for top B            */
} gmx_molblock_t;

typedef struct {
  t_grps         grps[egcNR];   /* Groups of things                     */
  int            ngrpname;      /* Number of groupnames                 */
  char           ***grpname;	/* Names of the groups		        */
  int            ngrpnr[egcNR];
  unsigned char  *grpnr[egcNR]; /* Group numbers or NULL		*/
} gmx_groups_t;

/* This macro gives the group number of group type egc for atom i.
 * This macro is useful, since the grpnr pointers are NULL
 * for group types that have all entries 0.
 */
#define ggrpnr(groups,egc,i) ((groups)->grpnr[egc] ? (groups)->grpnr[egc][i] : 0)

/* The global topology struct, based on molecule types */
typedef struct {
  char           **name;	/* Name of the topology	       	        */
  gmx_ffparams_t ffparams;
  int            nmoltype;
  gmx_moltype_t  *moltype;
  int            nmolblock;
  gmx_molblock_t *molblock;
  int            natoms;
  t_atomtypes    atomtypes;     /* Atomtype properties                  */
  t_block        mols;          /* The molecules                        */
  gmx_groups_t   groups;
  t_symtab	 symtab;        /* The symbol table			*/
} gmx_mtop_t;

/* The mdrun local topology struct, completely written out */
typedef struct {
  t_idef	idef;		/* The interaction function definition	*/
  t_atomtypes   atomtypes;      /* Atomtype properties                  */
  t_block       cgs;            /* The charge groups                    */
  t_blocka      excls;          /* The exclusions                       */
} gmx_localtop_t;

/* The old topology struct, completely written out, used in analysis tools */
typedef struct {
  char  	**name;		/* Name of the topology	       	        */
  t_idef	idef;		/* The interaction function definition	*/
  t_atoms	atoms;		/* The atoms		       	        */
  t_atomtypes   atomtypes;      /* Atomtype properties                  */
  t_block       cgs;            /* The charge groups                    */
  t_block       mols;           /* The molecules                        */
  t_blocka      excls;          /* The exclusions                       */
  t_symtab	symtab;		/* The symbol table			*/
} t_topology;
