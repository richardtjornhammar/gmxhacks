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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "assert.h"
#include "physics.h"
#include "macros.h"
#include "vec.h"
#include "force.h"
#include "invblock.h"
#include "confio.h"
#include "names.h"
#include "network.h"
#include "pbc.h"
#include "ns.h"
#include "nrnb.h"
#include "bondf.h"
#include "mshift.h"
#include "txtdump.h"
#include "copyrite.h"
#include "qmmm.h"
#include <stdio.h>
#include <string.h>
#include "gmx_fatal.h"
#include "typedefs.h"
#include <stdlib.h>
#include "cluster.h"
#include <errno.h>

/* Gaussian interface routines */

void init_gaussian(t_commrec *cr, t_QMrec *qm, t_MMrec *mm)
{
  FILE *rffile=NULL,*out=NULL;
  FILE *npalog;
  ivec
    basissets[eQMbasisNR]={{0,3,0},
			   {0,3,0},/*added for double sto-3g entry in names.c*/
			   {5,0,0},
			   {5,0,1},
			   {5,0,11},
			   {5,6,0},
			   {1,6,0},
			   {1,6,1},
			   {1,6,11},
			   {4,6,0}};
  char
    *buf;
  int
    i;
  
  /* using the ivec above to convert the basis read form the mdp file
   * in a human readable format into some numbers for the gaussian
   * route. This is necessary as we are using non standard routes to
   * do SH.
   */

  /* per layer we make a new subdir for integral file, checkpoint
   * files and such. These dirs are stored in the QMrec for
   * convenience 
   */

  if(qm->bPrint_npa){
    npalog=fopen("npa.log","w");
    fclose(npalog);
  } 

  if(!qm->nQMcpus){ /* this we do only once per layer 
		     * as we call g01 externally 
		     */

    for(i=0;i<DIM;i++)
      qm->SHbasis[i]=basissets[qm->QMbasis][i];

  /* init gradually switching on of the SA */
    qm->SAstep = 0;
  /* we read the number of cpus and environment from the environment
   * if set.  
   */
    snew(buf,MAXEL);
    buf = getenv("GAUSS_CPUS"); //GOT ERROR FROM OPENMPI WHEN ASSIGNING NCPUS
    if (buf)
      sscanf(buf,"%d",&qm->nQMcpus);
    else
      qm->nQMcpus=1;

    fprintf(stderr,"number of CPUs for gaussian = %d\n",qm->nQMcpus);
    snew(buf,MAXEL);
    buf = getenv("MEM");
    if (buf)
      sscanf(buf,"%d",&qm->QMmem);
    else
      qm->QMmem=50000000;

    fprintf(stderr,"memory for gaussian = %d\n",qm->QMmem);
    snew(buf,MAXEL);
    buf = getenv("ACC");
    if (buf)
      sscanf(buf,"%d",&qm->accuracy);
    else
      qm->accuracy=8;  

    fprintf(stderr,"accuracy in l510 = %d\n",qm->accuracy); 
    snew(buf,MAXEL);
    buf = getenv("CPMCSCF");
    if (buf)
      sscanf(buf,"%d",&qm->cpmcscf);
    else
      qm->cpmcscf=0;
    if (qm->cpmcscf)
      fprintf(stderr,"using cp-mcscf in l1003\n");
    else
      fprintf(stderr,"NOT using cp-mcscf in l1003\n"); 
    snew(buf,MAXEL);
    buf = getenv("SASTEP");
    if (buf)
      sscanf(buf,"%d",&qm->SAstep);
    else
      /* init gradually switching on of the SA */
      qm->SAstep = 0;
    /* we read the number of cpus and environment from the environment
     * if set.  
     */
    fprintf(stderr,"Level of SA at start = %d\n",qm->SAstep);
    
    /* punch the LJ C6 and C12 coefficients to be picked up by
     * gaussian and usd to compute the LJ interaction between the
     * MM and QM atoms.
     */
    if(qm->bTS||qm->bOPT){
      out = fopen("LJ.dat","w");
      for(i=0;i<qm->nrQMatoms;i++){

#ifdef GMX_DOUBLE
	fprintf(out,"%3d  %10.7lf  %10.7lf\n",
		qm->atomicnumberQM[i],qm->c6[i],qm->c12[i]);
#else
	fprintf(out,"%3d  %10.7f  %10.7f\n",
		qm->atomicnumberQM[i],qm->c6[i],qm->c12[i]);
#endif
      }
      fclose(out);
    }
    /* gaussian settings on the system */
    snew(buf,MAXEL);
    buf = getenv("GAUSS_DIR");
    fprintf(stderr,"%s",buf);

    if (buf){
      snew(qm->gauss_dir,MAXEL);
      sscanf(buf,"%s",qm->gauss_dir);
    }
    else
      gmx_fatal(FARGS,"no $GAUSS_DIR, check gaussian manual\n");
    
    snew(buf,MAXEL);    
    buf = getenv("GAUSS_EXE");
    if (buf){
      snew(qm->gauss_exe,MAXEL);
      sscanf(buf,"%s",qm->gauss_exe);
    }
    else
      gmx_fatal(FARGS,"no $GAUSS_EXE, check gaussian manual\n");
    
    snew(buf,MAXEL);
    buf = getenv("DEVEL_DIR");
    if (buf){
      snew(qm->devel_dir,MAXEL);
      sscanf(buf,"%s",qm->devel_dir);
    }
    else
      gmx_fatal(FARGS,"no $DEVEL_DIR, this is were the modified links reside.\n");
    
    
    /*  if(fr->bRF){*/
    /* reactionfield, file is needed using gaussian */
    /*    rffile=fopen("rf.dat","w");*/
    /*   fprintf(rffile,"%f %f\n",fr->epsilon_r,fr->rcoulomb/BORH2NM);*/
    /* fclose(rffile);*/
    /*  }*/
  }
  fprintf(stderr,"gaussian initialised...\n");
}

void write_gaussian_SH_input(int step,bool swap,
			     t_forcerec *fr, t_QMrec *qm, t_MMrec *mm)
{
  int
    i;
  bool
    bSA;
  FILE
    *out;
  t_QMMMrec
    *QMMMrec;
  QMMMrec = fr->qr;
  bSA = (qm->SAstep>0);

  out = fopen("input.com","w");
  /* write the route */
  fprintf(out,"%s","%scr=input\n");
  fprintf(out,"%s","%rwf=input\n");
  fprintf(out,"%s","%int=input\n");
  fprintf(out,"%s","%d2e=input\n");
/*  if(step)
 *   fprintf(out,"%s","%nosave\n");
 */
  fprintf(out,"%s","%chk=input\n");
  fprintf(out,"%s%d\n","%mem=",qm->QMmem);
  fprintf(out,"%s%3d\n","%nprocshare=",qm->nQMcpus);

  /* use the versions of
   * l301 that computes the interaction between MM and QM atoms.
   * l510 that can punch the CI coefficients
   * l701 that can do gradients on MM atoms 
   */

  /* local version */
  fprintf(out,"%s%s%s",
	  "%subst l510 ",
	  qm->devel_dir,
	  "/l510\n");
  fprintf(out,"%s%s%s",
	  "%subst l301 ",
	  qm->devel_dir,
	  "/l301\n");
  fprintf(out,"%s%s%s",
	  "%subst l701 ",
	  qm->devel_dir,
	  "/l701\n");
  
  fprintf(out,"%s%s%s",
	  "%subst l1003 ",
	  qm->devel_dir,
	  "/l1003\n");
  fprintf(out,"%s%s%s",
	  "%subst l9999 ",
	  qm->devel_dir,
	  "/l9999\n");
  /* print the nonstandard route 
   */
  fprintf(out,"%s",
	  "#P nonstd\n 1/18=10,20=1,38=1/1;\n");
  fprintf(out,"%s",
	  " 2/9=110,15=1,17=6,18=5,40=1/2;\n");
  if(mm->nrMMatoms)
    fprintf(out,
	    " 3/5=%d,6=%d,7=%d,25=1,32=1,43=1,94=-2/1,2,3;\n",
	    qm->SHbasis[0],
	    qm->SHbasis[1],
	    qm->SHbasis[2]); /*basisset stuff */
  else
    fprintf(out,
	    " 3/5=%d,6=%d,7=%d,25=1,32=1,43=0,94=-2/1,2,3;\n",
	    qm->SHbasis[0],
	    qm->SHbasis[1],
	    qm->SHbasis[2]); /*basisset stuff */
  /* development */
  if (step+1) /* fetch initial guess from check point file */
    /* hack, to always read from chk file!!!!! */
    fprintf(out,"%s%d,%s%d%s"," 4/5=1,7=6,17=",
	    qm->CASelectrons,
	    "18=",qm->CASorbitals,"/1,5;\n");
  else /* generate the first checkpoint file */
    fprintf(out,"%s%d,%s%d%s"," 4/5=0,7=6,17=",
	    qm->CASelectrons,
	    "18=",qm->CASorbitals,"/1,5;\n");
  /* the rest of the input depends on where the system is on the PES  */
  if(swap && bSA){ /* make a slide to the other surface */
    if(qm->CASorbitals>6){  /* use direct and no full diag */
      fprintf(out," 5/5=2,16=-2,17=10000000,28=2,32=2,38=6,97=100/10;\n");
    } 
    else {
      if(qm->cpmcscf){
	fprintf(out," 5/5=2,6=%d,17=31000200,28=2,32=2,38=6,97=100/10;\n",
		qm->accuracy);
	if(mm->nrMMatoms>0)
	  fprintf(out," 7/7=1,16=-2,30=1/1;\n");
	fprintf(out," 11/31=1,42=1,45=1/1;\n");
	fprintf(out," 10/6=1,10=700006,28=2,29=1,31=1,97=100/3;\n");
	fprintf(out," 7/30=1/16;\n 99/10=4/99;\n");
      }
      else{
	fprintf(out," 5/5=2,6=%d,17=11000000,28=2,32=2,38=6,97=100/10;\n",
		qm->accuracy);
	fprintf(out," 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
      }
    }
  }
  else if(bSA){ /* do a "state-averaged" CAS calculation */
    if(qm->CASorbitals>6){ /* no full diag */ 
      fprintf(out," 5/5=2,16=-2,17=10000000,28=2,32=2,38=6/10;\n");
    } 
    else {
      if(qm->cpmcscf){
	fprintf(out," 5/5=2,6=%d,17=31000200,28=2,32=2,38=6/10;\n",
		qm->accuracy);
	if(mm->nrMMatoms>0)
	  fprintf(out," 7/7=1,16=-2,30=1/1;\n");
	fprintf(out," 11/31=1,42=1,45=1/1;\n");
	fprintf(out," 10/6=1,10=700006,28=2,29=1,31=1/3;\n");
	fprintf(out," 7/30=1/16;\n 99/10=4/99;\n");
      }
      else{
      	fprintf(out," 5/5=2,6=%d,17=11000000,28=2,32=2,38=6/10;\n",
		qm->accuracy);
	fprintf(out," 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
      }
    }
  }
  else if(swap){/* do a "swapped" CAS calculation */
    if(qm->CASorbitals>6)
      fprintf(out," 5/5=2,16=-2,17=0,28=2,32=2,38=6,97=100/10;\n");
    else
      fprintf(out," 5/5=2,6=%d,17=1000000,28=2,32=2,38=6,97=100/10;\n",
	      qm->accuracy);
    fprintf(out," 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
  }
  else {/* do a "normal" CAS calculation */
    if(qm->CASorbitals>6)
      fprintf(out," 5/5=2,16=-2,17=0,28=2,32=2,38=6/10;\n");
    else
      fprintf(out," 5/5=2,6=%d,17=1000000,28=2,32=2,38=6/10;\n",
	      qm->accuracy);
    fprintf(out," 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
  }
  fprintf(out, "\ninput-file generated by gromacs\n\n");
  fprintf(out,"%2d%2d\n",qm->QMcharge,qm->multiplicity);
  for (i=0;i<qm->nrQMatoms;i++){
#ifdef GMX_DOUBLE
    fprintf(out,"%3d %10.7lf  %10.7lf  %10.7lf\n",
	    qm->atomicnumberQM[i],
	    qm->xQM[i][XX]/BORH2NM,
	    qm->xQM[i][YY]/BORH2NM,
	    qm->xQM[i][ZZ]/BORH2NM);
#else
    fprintf(out,"%3d %10.7f  %10.7f  %10.7f\n",
	    qm->atomicnumberQM[i],
	    qm->xQM[i][XX]/BORH2NM,
	    qm->xQM[i][YY]/BORH2NM,
	    qm->xQM[i][ZZ]/BORH2NM);
#endif
  }
  /* MM point charge data */
  if(QMMMrec->QMMMscheme!=eQMMMschemeoniom && mm->nrMMatoms){
    fprintf(out,"\n");
    for(i=0;i<mm->nrMMatoms;i++){
#ifdef GMX_DOUBLE
      fprintf(out,"%10.7lf  %10.7lf  %10.7lf %8.4lf\n",
	      mm->xMM[i][XX]/BORH2NM,
	      mm->xMM[i][YY]/BORH2NM,
	      mm->xMM[i][ZZ]/BORH2NM,
	      mm->MMcharges[i]);
#else
      fprintf(out,"%10.7f  %10.7f  %10.7f %8.4f\n",
	      mm->xMM[i][XX]/BORH2NM,
	      mm->xMM[i][YY]/BORH2NM,
	      mm->xMM[i][ZZ]/BORH2NM,
	      mm->MMcharges[i]);
#endif
    }
  }
  if(bSA) {/* put the SA coefficients at the end of the file */
#ifdef GMX_DOUBLE
    fprintf(out,"\n%10.8lf %10.8lf\n",
	    qm->SAstep*0.5/qm->SAsteps,
	    1-qm->SAstep*0.5/qm->SAsteps);
#else    
    fprintf(out,"\n%10.8f %10.8f\n",
	    qm->SAstep*0.5/qm->SAsteps,
	    1-qm->SAstep*0.5/qm->SAsteps);
#endif
    fprintf(stderr,"State Averaging level = %d/%d\n",qm->SAstep,qm->SAsteps);
  }
  fprintf(out,"\n");
  fclose(out);
}  /* write_gaussian_SH_input */



void write_gaussian_input(int step ,t_forcerec *fr, t_QMrec *qm, t_MMrec *mm)
{
  int     i,i0,j,qcflag=0,NRC=0,k=0,l=0,T,q=0,qq=0,cl;
  int     nrZ[MAXFILES],nrQ[MAXFILES],nrCOL=0,nrMM=0,nrQM=0,a,b;
  t_QMMMrec *QMMMrec;
  FILE    *out;
  char    wstring[MAXEL],dstring[MAXEL],tstring[MAXEL],estring[MAXEL];
  double  tchrg[MAXFILES],chrg;
  int     multip[MAXFILES],netchrg[MAXFILES];
  FILE    *clustF[MAXFILES]; 

  strcpy(tstring,"input");
  strcpy(dstring,"input");
  strcpy(estring,".com");
  NRC=qm->nrQMc;   nrQM=qm->nrQMatoms;
  nrMM=mm->nrMMatoms; nrCOL=nrMM+nrQM;

  if(NRC>1000){
    gmx_fatal(FARGS,"To many clusters. Your QM system is probably exploding");
  }
  
  QMMMrec = fr->qr;

  for(cl=0;cl<NRC;cl++){
    tchrg[cl]=0.0;multip[cl]=0.0;
  }
  for(i=0;i<nrQM;i++){
    j=qm->QMcluster[i*2]; k=qm->QMcluster[i*2+1];
    tchrg[j]+=qm->pcharges[k];
  }

  for(cl=0;cl<NRC;cl++){
    chrg=0.0;
    for(i=0;i<nrQM;i++){
      j=qm->QMcluster[i*2]; k=qm->QMcluster[i*2+1];
      if(cl==j){
	chrg+=qm->atomicnumberQM[k];
      }
    }
    /* 
     * computed for whole number of electrons 
     * NOTE: if this gets wrong then 
     * change the charge for the QM atoms in the topology 
     */
    multip[cl]=((fmod( (int)roundChrg(tchrg[cl])*(-1.0)+chrg,2.0))==0.0)?1:2;
    chrg=tchrg[cl];
    netchrg[cl]=(int)roundChrg( chrg );
  }

  for(cl=0;cl<NRC;cl++){
    sprintf(dstring,"input%d",cl);
    clustF[cl]=fopen(strcat(dstring,estring),"w");
    if(qm->QMmethod>=eQMmethodRHF)
      fprintf(clustF[cl],"%s%d\n",
	      "%chk=input",cl);
    else
      fprintf(clustF[cl],"%s",
	      "%chk=se\n");
    if(qm->nQMcpus>1)
      fprintf(clustF[cl],"%s%3d\n",
	      "%nprocshare=",qm->nQMcpus);
    fprintf(clustF[cl],"%s%d\n",
	    "%mem=",qm->QMmem);
    /* use the modified links that include the LJ contribution at the QM level */
    if(qm->bTS||qm->bOPT){
      fprintf(clustF[cl],"%s%s%s",
	      "%subst l701 ",qm->devel_dir,"/l701_LJ\n");
      fprintf(clustF[cl],"%s%s%s",
	      "%subst l301 ",qm->devel_dir,"/l301_LJ\n");
    }
    else{
      fprintf(clustF[cl],"%s%s%s",
	      "%subst l701 ",qm->devel_dir,"/l701\n");
      fprintf(clustF[cl],"%s%s%s",
	      "%subst l301 ",qm->devel_dir,"/l301\n");
    }
    fprintf(clustF[cl],"%s%s%s",
	    "%subst l9999 ",qm->devel_dir,"/l9999\n");
    if(step){
      fprintf(clustF[cl],"%s",
	      "#T ");
    }else{
      fprintf(clustF[cl],"%s",
	      "#P ");
    }
    
    if(qm->QMmethod==eQMmethodB3LYPLAN){
      fprintf(clustF[cl]," %s", 
	      "B3LYP/GEN Pseudo=Read");
    }
    else{
      fprintf(clustF[cl]," %s", 
	      eQMmethod_names[qm->QMmethod]);
      
      if(qm->QMmethod>=eQMmethodRHF){
	fprintf(clustF[cl],"/%s",
		eQMbasis_names[qm->QMbasis]);
	if(qm->QMmethod==eQMmethodCASSCF){
	  /* in case of cas, how many electrons and orbitals do we need? */
	  fprintf(clustF[cl],"(%d,%d)",
		  qm->CASelectrons,qm->CASorbitals);
	}
      }
    }
    if(QMMMrec->QMMMscheme==eQMMMschemenormal){
      if((qm->includeQ>=eIncludeYES))
        fprintf(clustF[cl]," %s","Charge ");
      if(qm->bPrint_npa)
        fprintf(clustF[cl]," %s","Pop=MK ");
    }
    if (step || qm->QMmethod==eQMmethodCASSCF){
      /* fetch guess from checkpoint file, always for CASSCF */
      fprintf(clustF[cl],"%s"," guess=read");
    }
    fprintf(clustF[cl],"\nNosymm units=bohr\n");
    
    if(qm->bTS){
      fprintf(clustF[cl],"OPT=(Redundant,TS,noeigentest,ModRedundant) Punch=(Coord,Derivatives) ");
    }
    else if (qm->bOPT){
      fprintf(clustF[cl],"OPT=(Redundant,ModRedundant) Punch=(Coord,Derivatives) ");
    }
    else{
      fprintf(clustF[cl],"FORCE Punch=(Derivatives) ");
    }
    fprintf(clustF[cl],"iop(3/33=1)\n\n");
    fprintf(clustF[cl], "input-file generated by gromacs\n\n");

    if(qm->includeQ<=eIncludeQQ){
      fprintf(clustF[cl],"%4d %4d\n",netchrg[cl],multip[cl]); //apparently not sensitive to this formating
    }
    else {
      fprintf(clustF[cl],"%2d%2d\n",qm->QMcharge,qm->multiplicity);
    }
  }
  
  for (i=0;i<qm->nrQMatoms;i++){
    j=qm->QMcluster[i*2];k=qm->QMcluster[i*2+1];
    if(qm->useQMatom[k]){
#ifdef GMX_DOUBLE
    fprintf(clustF[j],"%3d %10.7lf  %10.7lf  %10.7lf\n",
	    qm->atomicnumberQM[k],
	    qm->rQM[k][XX]/BORH2NM,
	    qm->rQM[k][YY]/BORH2NM,
	    qm->rQM[k][ZZ]/BORH2NM);
#else
    fprintf(clustF[j],"%3d %10.7f  %10.7f  %10.7f\n",
	    qm->atomicnumberQM[k],
	    qm->rQM[k][XX]/BORH2NM,
	    qm->rQM[k][YY]/BORH2NM,
	    qm->rQM[k][ZZ]/BORH2NM);
#endif
    }
  }

  for(cl=0;cl<NRC;cl++){
    /* Pseudo Potential and ECP are included here if selected (MEthod suffix LAN) */
    if(qm->QMmethod==eQMmethodB3LYPLAN){
      fprintf(clustF[cl],"\n");
      for(i=0;i<qm->nrQMatoms;i++){
	j=qm->QMcluster[i*2];k=qm->QMcluster[i*2+1];
	if(qm->atomicnumberQM[k]<21){
	  fprintf(clustF[j],"%d ",k+1); //internal coordinates for gromacs not gaussian
	}
      }
      fprintf(clustF[cl],"\n%s\n****\n",eQMbasis_names[qm->QMbasis]);
      
      for(i=0;i<qm->nrQMatoms;i++){
	j=qm->QMcluster[i*2];k=qm->QMcluster[i*2+1];
	if(qm->atomicnumberQM[k]>21){
	  fprintf(clustF[j],"%d ",k+1);
	}
      }
      fprintf(clustF[cl],"\n%s\n****\n\n","lanl2dz");    
      
      for(i=0;i<qm->nrQMatoms;i++){
	j=qm->QMcluster[i*2];k=qm->QMcluster[i*2+1];
	if(qm->atomicnumberQM[k]>21){
	  fprintf(clustF[j],"%d ",k+1);
	}
      }
      fprintf(clustF[cl],"\n%s\n","lanl2dz");    
    }    
    
    /* MM point charge data */
    if(QMMMrec->QMMMscheme!=eQMMMschemeoniom && (qm->includeQ>eIncludeNONE)){
      fprintf(clustF[cl],"\n");
      if(qm->bTS||qm->bOPT){ //  NOTE OPTIMISATION IS TOTALLY UNCHECKED !!
	/* freeze the frontier QM atoms and Link atoms. This is
	 * important only if a full QM subsystem optimization is done
	 * with a frozen MM environmeent. For dynamics, or gromacs's own
	 * optimization routines this is not important.
	 */
	for(i=0;i<qm->nrQMatoms;i++){
	  j=qm->QMcluster[i*2];k=qm->QMcluster[i*2+1];
          if(qm->useQMatom[k]){
	    if(qm->frontatoms[k]){
	      fprintf(clustF[j],"%d F\n",k+1); /* counting from 1 || this is again using internal coordinates of gromacs not gaussian */
	    }
          }
	}
	/* MM point charges include LJ parameters in case of QM optimization
	 */
	for(i=0;i<mm->nrMMatoms;i++){
#ifdef GMX_DOUBLE
	  if(mm->MMcluster[ nrCOL*cl + i ]){
	    fprintf(clustF[cl],"%10.7lf  %10.7lf  %10.7lf %8.4lf 0.0 %10.7lf %10.7lf\n",
		    mm->rMM[i][XX]/BORH2NM,
		    mm->rMM[i][YY]/BORH2NM,
		    mm->rMM[i][ZZ]/BORH2NM,
		    mm->MMcharges[i],
		    mm->c6[i],mm->c12[i]);
	  }
#else
	  if(mm->MMcluster[ nrCOL*cl + i ]){
	    fprintf(clustF[cl],"%10.7f  %10.7f  %10.7f %8.4f 0.0 %10.7f %10.7f\n",
		    mm->rMM[i][XX]/BORH2NM,
		    mm->rMM[i][YY]/BORH2NM,
		    mm->rMM[i][ZZ]/BORH2NM,
		    mm->MMcharges[i],
		    mm->c6[i],mm->c12[i]);
	  }
#endif
	} 
	fprintf(clustF[cl],"\n");
      }
      else{
        if(mm->nrMMatoms>0){
	for(i=0;i<mm->nrMMatoms;i++){
#ifdef GMX_DOUBLE
	  if(mm->MMcluster[ nrCOL*cl + i ]){
	    fprintf(clustF[cl],"%10.7lf  %10.7lf  %10.7lf %8.4lf\n",
		    mm->rMM[i][XX]/BORH2NM,
		    mm->rMM[i][YY]/BORH2NM,
		    mm->rMM[i][ZZ]/BORH2NM,
		    mm->MMcharges[i]);
	  }
#else
	  if(mm->MMcluster[ nrCOL*cl + i ]){
	    fprintf(clustF[cl],"%10.7f  %10.7f  %10.7f %8.4f\n",
		    mm->rMM[i][XX]/BORH2NM,
		    mm->rMM[i][YY]/BORH2NM,
		    mm->rMM[i][ZZ]/BORH2NM,
		    mm->MMcharges[i]);
	  }
#endif
	}
        }
	//INCLUDE OTHER QM CLUSTERS QM ATOMS AS MM ATOMS 
	for(i=0;i<nrQM;i++){
	  a=qm->QMcluster[i*2];      /* CLUSTER INDEX a:(0 -> NC-1)    */
	  b=qm->QMcluster[i*2+1];    /* QM INDEX      b:(0 -> nrQM-1)  */
	  if(mm->MMcluster[b+nrMM+nrCOL*cl]){
	    fprintf(clustF[cl],"%10.7f  %10.7f  %10.7f %8.4f\n",
		    mm->rMM[ b+nrMM+nrCOL*cl ][XX]/BORH2NM,
		    mm->rMM[ b+nrMM+nrCOL*cl ][YY]/BORH2NM,
		    mm->rMM[ b+nrMM+nrCOL*cl ][ZZ]/BORH2NM,
		    qm->pcharges[b]);
	  }
	}
      }
    }
    fprintf(clustF[cl],"\n");
    fclose(clustF[cl]);
  }
  
  /* write the route */
}  /* write_gaussian_input */

real read_gaussian_output(rvec QMgrad[],rvec MMgrad[],int step,t_QMrec *qm, t_MMrec *mm, int I)
{
  int
    i,j,k,l,m,atnum,a,b,flag,trash;
  char
    buf[300];
  real
    QMener;
  FILE *in,*g03log,*npalog;
  int  HL=0,nQMc,logi,nrHERE=0,nrMM=0,nrQM=0,nrCOL=0;
  char nameMM[100],nameQA[100],g03fnm[100];
  char line[MAXEL],hang[2*MAXEL],lp1[100],lp2[100],mulliken[MAXEL];
  FILE *out,*fmmndx,*fqandx;
  rvec GRAD;

  nQMc=qm->nrQMc;
  nrMM=mm->nrMMatoms;
  nrQM=qm->nrQMatoms;
  nrCOL=nrMM+nrQM;

  //same fort.7 for all clusters since we are calling/reading sequentially from I
  in=fopen("fort.7","r");
  /* 
   * in case of an optimization, the coordinates are printed in the
   * fort.7 file first, followed by the energy, coordinates and (if
   * required) the CI eigenvectors.
   */

//NDX FILES FOR VISUALISATION OF THE QMMM REGIME
  sprintf(nameMM,"mmatoms%d.ndx",I);
  sprintf(nameQA,"qaatoms%d.ndx",I);

  fmmndx=fopen(nameMM,"w");
  fqandx=fopen(nameQA,"w");

  if(qm->bTS||qm->bOPT){ //WARNING OPTIMISATION NOT REALLY TESTED WITH CLUSTERING!
    for(i=0;i<nrQM;i++){
      if(qm->useQMatom[i]){
      if( NULL == fgets(buf,300,in))
	{
	  gmx_fatal(FARGS,"Error reading Gaussian output - not enough atom lines?");
	}
      
#ifdef GMX_DOUBLE
      sscanf(buf,"%d %lf %lf %lf\n",
	     &atnum,
	     &qm->xQM[i][XX],
	     &qm->xQM[i][YY],
	     &qm->xQM[i][ZZ]);
#else
      sscanf(buf,"%d %f %f %f\n",
	     &atnum,
	     &qm->xQM[i][XX],
	     &qm->xQM[i][YY],
	     &qm->xQM[i][ZZ]);
#endif      
      for(j=0;j<DIM;j++){
	qm->xQM[i][j]*=BORH2NM;
      }
    }
    }
  }
  /* the next line is the energy and in the case of CAS, the energy
   * difference between the two states.
   */
  if(NULL == fgets(buf,300,in))
  {
      gmx_fatal(FARGS,"Error reading Gaussian output");
  }

#ifdef GMX_DOUBLE
  sscanf(buf,"%lf\n",&QMener);
#else
  sscanf(buf,"%f\n", &QMener);
#endif
  /* next lines contain the gradients of the QM atoms */
  for(i=0;i<nrQM;i++){
    if(qm->QMcluster[2*i]==I){
      j=qm->QMcluster[2*i+1];        //QM ATOM INDEX
      if(qm->useQMatom[j]){
      if(NULL == fgets(buf,300,in))
	{
	  gmx_fatal(FARGS,"Error reading Gaussian output");
	}
#ifdef GMX_DOUBLE
      sscanf(buf,"%lf %lf %lf\n",
	     &GRAD[XX],
	     &GRAD[YY],
	     &GRAD[ZZ]);
#else
      sscanf(buf,"%f %f %f\n",
	     &GRAD[XX],
	     &GRAD[YY],
	     &GRAD[ZZ]);
#endif     
      for(k=XX;k<=ZZ;k++)
	QMgrad[j][k]+=GRAD[k];
      }
    }
  }
  /* the next lines are the gradients of the MM atoms */
  if(qm->QMmethod>=eQMmethodRHF && qm->includeQ>eIncludeNONE){  
    l=0;
    if(nrMM>0){
    for(i=0;i<nrMM;i++){
      logi=(int)mm->MMcluster[I*nrCOL+i];
      if(logi){
        fprintf(fmmndx,"%10d",mm->indexMM[i]);
        if(++l==20){
          fprintf(fmmndx,"\n");
          l=0;
        }
	if(NULL==fgets(buf,300,in))
	  {
	    gmx_fatal(FARGS,"Error reading Gaussian output");
	  }
#ifdef GMX_DOUBLE
	sscanf(buf,"%lf %lf %lf\n",
	       &GRAD[XX],
	       &GRAD[YY],
	       &GRAD[ZZ]);
#else
	sscanf(buf,"%f %f %f\n",
	       &GRAD[XX],
	       &GRAD[YY],
	       &GRAD[ZZ]);
#endif 
	for(k=XX;k<=ZZ;k++)
	  MMgrad[i][k]+=GRAD[k];	
      }
    }
    }
    l=0;
    for(i=0;i<nrQM;i++){
      j=qm->QMcluster[2*i+1];
      if(mm->MMcluster[(j+nrMM)+I*nrCOL]){
        fprintf(fqandx,"%10d",qm->indexQM[j]);
        if(++l==20){
          fprintf(fqandx,"\n");
          l=0;
        }
	if(NULL==fgets(buf,300,in))
	  {
	    gmx_fatal(FARGS,"Error reading Gaussian output");
	  }
#ifdef GMX_DOUBLE
	sscanf(buf,"%lf %lf %lf\n",
	       &GRAD[XX],
	       &GRAD[YY],
	       &GRAD[ZZ]);
#else
	sscanf(buf,"%f %f %f\n",
	       &GRAD[XX],
	       &GRAD[YY],
	       &GRAD[ZZ]);
#endif	
	for(k=XX;k<=ZZ;k++)
	  QMgrad[j][k]+=GRAD[k]; 
	  /* MM gradients for the QM atoms as a results of interacting with other QM regions */
      }
    }
  }

  if(qm->bPrint_npa){
    if(qm->dI>0)
      sprintf(g03fnm,"output%d.log",I);
    else
      sprintf(g03fnm,"output%d_2.log",I);

    g03log=fopen(g03fnm,"r");
    npalog=fopen("npa.log","a");
    flag=1;
    while(!feof(g03log) && flag){
      if( NULL == fgets(line,MAXEL,g03log)) {               
        ;         // Error handling 
      }
      if(!strncmp(line,"            Population analysis using the SCF density.",54) && flag ){      //found Mulliken charges
        sprintf(mulliken,"%s","##");
        while(!feof(g03log)){
          if( NULL == fgets(line,MAXEL,g03log)) {               
            ;                                // Error handling 
          }
          if(line[1]=='S'){
            break;                           // end of section   
          }
          if(line[1]=='M'){
            if( NULL == fgets(line,MAXEL,g03log)) {               
              ;                                // This line contains the number 1
            } 
            for(i=0;i<qm->nrQMatoms;i++){
              a=qm->QMcluster[i*2];          // CLUSTER INDEX a:(0 -> NC-1)    
              b=qm->QMcluster[i*2+1];        //      QM INDEX b:(0 -> nrQM-1)  
              if(a==I){
                if( NULL == fgets(line,MAXEL,g03log)) {               
                  ;                        // Error handling 
                }
		fprintf(npalog,"#%s",line);
              }
            }
            break;            
          }
        }
      } 
      if(!strncmp(line,"            Electrostatic Properties Using The SCF Density",58) ){      //found MSK charges
        while(!feof(g03log)){
          if( NULL == fgets(line,MAXEL,g03log)) {               
            ;                                // Error handling 
          }
          if(line[0]=='-'){
            flag=0;
            break;                           // end of section   
          }
          if(line[5]=='1' && flag){                  // found entry    ( we have atoms belong to this cluster following )
            k=0;
            for(i=0;i<qm->nrQMatoms;i++){
              a=qm->QMcluster[i*2];          // CLUSTER INDEX a:(0 -> NC-1)    
              b=qm->QMcluster[i*2+1];        //      QM INDEX b:(0 -> nrQM-1)  
              if(a==I){
                if(!(k==0)){
                  if( NULL == fgets(line,MAXEL,g03log)) {               
                    ;                        // Error handling 
                  }
                } // We have data here
                sprintf(hang,"%d %f %f %f %s",I,qm->rQM[b][XX]/BORH2NM,qm->rQM[b][YY]/BORH2NM,qm->rQM[b][ZZ]/BORH2NM,line);
                fprintf(npalog,"%s",hang);
                k++;
              }
            }
	    flag=0;
	    break;
          }
        }
      }
    } 
    fclose(g03log);
    fclose(npalog);
  } 

  fprintf(fqandx,"\n");
  fprintf(fmmndx,"\n");
  fclose(fqandx);
  fclose(fmmndx);
  fclose(in);
  return(QMener);  
}

real read_gaussian_SH_output(rvec QMgrad[],rvec MMgrad[],int step,
			     bool swapped,t_QMrec *qm, t_MMrec *mm)
{
  int
    i;
  char
    buf[300];
  real
    QMener,DeltaE;
  FILE
    *in;
  
  in=fopen("fort.7","r");
  /* first line is the energy and in the case of CAS, the energy
   * difference between the two states.
   */
  if(NULL == fgets(buf,300,in))
  {
      gmx_fatal(FARGS,"Error reading Gaussian output");
  }

#ifdef GMX_DOUBLE
  sscanf(buf,"%lf %lf\n",&QMener,&DeltaE);
#else
  sscanf(buf,"%f %f\n",  &QMener,&DeltaE);
#endif
  
  /* switch on/off the State Averaging */
  
  if(DeltaE > qm->SAoff){
    if (qm->SAstep > 0){
      qm->SAstep--;
    }
  }
  else if (DeltaE < qm->SAon || (qm->SAstep > 0)){
    if (qm->SAstep < qm->SAsteps){
      qm->SAstep++;
    }
  }
  
  /* for debugging: */
  fprintf(stderr,"Gap = %5f,SA = %3d\n",DeltaE,(qm->SAstep>0));
  /* next lines contain the gradients of the QM atoms */
  for(i=0;i<qm->nrQMatoms;i++){
    if(NULL==fgets(buf,300,in))
    {
	gmx_fatal(FARGS,"Error reading Gaussian output");
    }

#ifdef GMX_DOUBLE
    sscanf(buf,"%lf %lf %lf\n",
	   &QMgrad[i][XX],
	   &QMgrad[i][YY],
	   &QMgrad[i][ZZ]);
#else
    sscanf(buf,"%f %f %f\n",
	   &QMgrad[i][XX],
	   &QMgrad[i][YY],
	   &QMgrad[i][ZZ]);
#endif     
  }
  /* the next lines, are the gradients of the MM atoms */
  
  for(i=0;i<mm->nrMMatoms;i++){
    if(NULL==fgets(buf,300,in))
    {
	gmx_fatal(FARGS,"Error reading Gaussian output");
    }
#ifdef GMX_DOUBLE
    sscanf(buf,"%lf %lf %lf\n",
	   &MMgrad[i][XX],
	   &MMgrad[i][YY],
	   &MMgrad[i][ZZ]);
#else
    sscanf(buf,"%f %f %f\n",
	   &MMgrad[i][XX],
	   &MMgrad[i][YY],
	   &MMgrad[i][ZZ]);
#endif	
  }
  
  /* the next line contains the two CI eigenvector elements */
  if(NULL==fgets(buf,300,in))
  {
      gmx_fatal(FARGS,"Error reading Gaussian output");
  }
  if(!step){
    sscanf(buf,"%d",&qm->CIdim);
    snew(qm->CIvec1,qm->CIdim);
    snew(qm->CIvec1old,qm->CIdim);
    snew(qm->CIvec2,qm->CIdim);
    snew(qm->CIvec2old,qm->CIdim);
  } else {
    /* before reading in the new current CI vectors, copy the current
     * CI vector into the old one.
     */
    for(i=0;i<qm->CIdim;i++){
      qm->CIvec1old[i] = qm->CIvec1[i];
      qm->CIvec2old[i] = qm->CIvec2[i];
    }
  }
  /* first vector */
  for(i=0;i<qm->CIdim;i++){
    if(NULL==fgets(buf,300,in))
    {
	gmx_fatal(FARGS,"Error reading Gaussian output");
    }
#ifdef GMX_DOUBLE
    sscanf(buf,"%lf\n",&qm->CIvec1[i]);
#else
    sscanf(buf,"%f\n", &qm->CIvec1[i]);   
#endif
  }
  /* second vector */
  for(i=0;i<qm->CIdim;i++){
    if(NULL==fgets(buf,300,in))
    {
	gmx_fatal(FARGS,"Error reading Gaussian output");
    }
#ifdef GMX_DOUBLE
    sscanf(buf,"%lf\n",&qm->CIvec2[i]);
#else
    sscanf(buf,"%f\n", &qm->CIvec2[i]);   
#endif
  }
  fclose(in);
  return(QMener);  
}

real inproduct(real *a, real *b, int n)
{
  int
    i;
  real
    dot=0.0;
  
  /* computes the inner product between two vectors (a.b), both of
   * which have length n.
   */  
  for(i=0;i<n;i++){
    dot+=a[i]*b[i];
  }
  return(dot);
}

int hop(int step, t_QMrec *qm)
{
  int
    swap = 0;
  real
    d11=0.0,d12=0.0,d21=0.0,d22=0.0;
  
  /* calculates the inproduct between the current Ci vector and the
   * previous CI vector. A diabatic hop will be made if d12 and d21
   * are much bigger than d11 and d22. In that case hop returns true,
   * otherwise it returns false.
   */  
  if(step){ /* only go on if more than one step has been done */
    d11 = inproduct(qm->CIvec1,qm->CIvec1old,qm->CIdim);
    d12 = inproduct(qm->CIvec1,qm->CIvec2old,qm->CIdim);
    d21 = inproduct(qm->CIvec2,qm->CIvec1old,qm->CIdim);
    d22 = inproduct(qm->CIvec2,qm->CIvec2old,qm->CIdim);
  }
  fprintf(stderr,"-------------------\n");
  fprintf(stderr,"d11 = %13.8f\n",d11);
  fprintf(stderr,"d12 = %13.8f\n",d12);
  fprintf(stderr,"d21 = %13.8f\n",d21);
  fprintf(stderr,"d22 = %13.8f\n",d22);
  fprintf(stderr,"-------------------\n");
  
  if((fabs(d12)>0.5)&&(fabs(d21)>0.5))
    swap = 1;
  
  return(swap);
}

void do_gaussian(int step, char *exe, int I,t_QMrec *qm)
{
  char    buf[1000],line[100],st1[20],st2[50];
  char    inname[100],outname[100],in2name[100],ut2namn[100];
  FILE    *log,*com,*newcom;
  int     nrcase=0;

  /* make the call to the gaussian binary through system()
   * The location of the binary will be picked up from the 
   * environment using getenv().
   */

  sprintf(inname,"input%d.com",I);
  sprintf(outname,"output%d.log",I);
  sprintf(in2name,"input%d_2.com",I);
  sprintf(ut2namn,"output%d_2.log",I);  
  qm->dI=(I+1);

  if(step) /* hack to prevent long inputfiles */
    sprintf(buf,"%s < %s > %s",
	    exe,
	    inname,
	    outname);
  else
    sprintf(buf,"%s < %s > %s",
	    exe,
            inname,
	    outname);
  fprintf(stderr,"Calling '%s'\n",buf);
#ifdef GMX_NO_SYSTEM
  printf("Warning-- No calls to system(3) supported on this platform.");
  gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#else
  if ( system(buf) != 0 ){
    fprintf(stderr,"FAILED::TESTING NEW CASE\n");
    log=fopen(outname,"r+");
    qm->dI=-1*(I+1);
    //What happened?
    while(!feof(log)){
      fgets(line,100,log);
      if(!strncmp(line," NtrErr Called from FileIO.",26))
	nrcase=1;
      if(!strncmp(line," The SCF is confused.",20))
	nrcase=2;
      if(!strncmp(line," Convergence failure -- run terminated.",38))
	nrcase=3;
      if(!strncmp(line," Initial guess read from the checkpoint file:",45))
	nrcase=4;
      if(!strncmp(line," RdChkP:  Unable to locate IRWF",31))
	nrcase=5;
      if(!strncmp(line," Rare condition: small coef for last iteration:",48))
	nrcase=6;
      if(!strncmp(line," Error translating basis functions",35))
        nrcase=7;
    }
    fclose(log);
    
    if(nrcase>0){
      com=fopen(inname,"r");
      newcom=fopen(in2name,"w");
      switch(nrcase){
      case 1:
	while(!feof(com)){
	  fgets(line,100,com);
	  if(line[0]=='#'){
	    sscanf(line,"%s %s",st1,st2);
          if(qm->bPrint_npa)
	    fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
          else
            fprintf(newcom,"%s %s Guess=Mix Charge\n",st1,st2);
	  }
	  else
	    fprintf(newcom,"%s",line);
	}
	break;
      case 2: 
	while(!feof(com)){
	  fgets(line,100,com);
	  if(line[0]=='#'){
	    sscanf(line,"%s %s",st1,st2);
            if(qm->bPrint_npa)
              fprintf(newcom,"%s %s Guess=Mix Pop=MK\n",st1,st2);
            else
	      fprintf(newcom,"%s %s Guess=Mix \n",st1,st2);
	  }
	  else
	    fprintf(newcom,"%s",line);
	}
	break;
      case 3:
	while(!feof(com)){
	  fgets(line,100,com);
	  if(line[0]=='#'){
	    sscanf(line,"%s %s",st1,st2);
            if(qm->bPrint_npa)
	      fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
            else
	      fprintf(newcom,"%s %s Guess=Mix Charge\n",st1,st2); //SCF=NODIIS
	  }
	  else
	    fprintf(newcom,"%s",line);
	}
	break;
      case 4:
	while(!feof(com)){
	  fgets(line,100,com);
	  if(line[0]=='#'){
	    sscanf(line,"%s %s",st1,st2);
            if(qm->bPrint_npa)
	      fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
            else
	      fprintf(newcom,"%s %s Guess=Mix Charge\n",st1,st2);
	  }
	  else
	    fprintf(newcom,"%s",line);
	}
	break;
      case 5:
	while(!feof(com)){
	  fgets(line,100,com);
	  if(line[0]=='#'){
	    sscanf(line,"%s %s",st1,st2);
            if(qm->bPrint_npa)
	      fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
            else
	      fprintf(newcom,"%s %s Guess=Mix Charge\n",st1,st2);
	  }
	  else
	    fprintf(newcom,"%s",line);
	}
	break;
      case 6:
	while(!feof(com)){
	  fgets(line,100,com);
	  if(line[0]=='#'){
	    sscanf(line,"%s %s",st1,st2);
            if(qm->bPrint_npa)
	      fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
            else
	      fprintf(newcom,"%s %s Guess=Mix Charge\n",st1,st2);
	  }
	  else
	    fprintf(newcom,"%s",line);
	}
	break;
      case 7:
	while(!feof(com)){
	  fgets(line,100,com);
	  if(line[0]=='#'){
	    sscanf(line,"%s %s",st1,st2);
            if(qm->bPrint_npa)
	      fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
            else
              fprintf(newcom,"%s %s Guess=Mix Charge\n",st1,st2);
	  }
	  else
	    fprintf(newcom,"%s",line);
	}
	break;
      default:
	while(!feof(com)){
	  fgets(line,100,com);
	  if(line[0]=='#'){
	    sscanf(line,"%s %s",st1,st2);
            if(qm->bPrint_npa)
	      fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
            else
	      fprintf(newcom,"%s %s Guess=Mix Charge\n",st1,st2);
	  }
	  else
	    fprintf(newcom,"%s",line);
	}
        break;
      }
      fclose(com);
      fclose(newcom);

      if(step) /* hack to prevent long inputfiles */
	sprintf(buf,"%s < %s > %s",
		exe,
		in2name,
		ut2namn);
      else
	sprintf(buf,"%s < %s > %s",
		exe,
		in2name,
		ut2namn);
      fprintf(stderr,"Calling '%s'\n\n",buf);
      if ( system(buf) != 0 ){
	gmx_fatal(FARGS,"Call to '%s' failed twice\n",buf);
      }
    }
    else{
      gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
    }
    //IF ONE OF THE CASES THEN WE REDO ELSE FAIL
  }
#endif
}

real call_gaussian(t_commrec *cr,  t_forcerec *fr, 
		   t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[])
{
  /* normal gaussian jobs */
  static int
    step=0;
  int
    i,j,NFILES,k,l,a,b,aa,bb;
  real
    QMener=0.0,QME0,rij;
  rvec
    *QMgrad,*MMgrad;
  char
    *exe;
  FILE *npalog;

  /* Allocate */
  snew(QMgrad,qm->nrQMatoms);
  snew(MMgrad,mm->nrMMatoms);  
  for(i=0;i<qm->nrQMatoms;i++){
    QMgrad[i][XX]=0.0; QMgrad[i][YY]=0.0; QMgrad[i][ZZ]=0.0;
  }
  for(i=0;i<mm->nrMMatoms;i++){
    MMgrad[i][XX]=0.0; MMgrad[i][YY]=0.0; MMgrad[i][ZZ]=0.0;
  }
  /* Gradient vectors */

  //BEGIN NEW QMMM CLUSTERING DEPENDENT CODE
  QMener=0.0;
  for(i=0;i<qm->nrQMatoms;i++){
    for(j=XX;j<=ZZ;j++){
      QMgrad[i][j]=0.0;
      qm->rQM[i][j]=qm->xQM[i][j];
      mm->rMM[i+mm->nrMMatoms][j]=qm->xQM[i][j];
      f[i][j]=0.0;
    }
  }
  for(i=0;i<mm->nrMMatoms;i++){
    for(j=XX;j<=ZZ;j++){
      mm->rMM[i][j]=mm->xMM[i][j];
      f[i+qm->nrQMatoms][j]=0.0;
      MMgrad[i][j]=0.0;
    }
  }
  qm->dI=0;

  if(qm->bCLUSTER){
    clustering(qm);  
    clusterMM(mm,qm);
    NFILES=qm->nrQMc; 
    //fprintf(stdout,"\nHERE_CL\n");
  }
  else{
    qm->nrQMc=1;NFILES=1;
    one_cluster(mm,qm);
  }
  if(NFILES<=0 || NFILES>MAXFILES){ /* This will only happen if the clustering has a bug */
    gmx_fatal(FARGS,"CLUSTERED GAUSSIAN FAILED");
  }
  //END NEW QMMM CLUSTERING DEPENDENT CODE
  snew(exe,1000);
  sprintf(exe,"%s/%s",qm->gauss_dir,qm->gauss_exe);


  if(qm->bPrint_npa){
    npalog=fopen("npa.log","a");
    fprintf(npalog,"@%d\n",qm->nrQMc); //every new step is marked with @ and the number of clusters in that step
    fclose(npalog);
  }
  write_gaussian_input(step,fr,qm,mm);

  for(i=0;i<NFILES;i++){
    fflush(stdout); fflush(stdin); 
    do_gaussian(step,exe,i,qm);
    QMener += read_gaussian_output(QMgrad,MMgrad,step,qm,mm,i);
  }

  // NOW ADD DISPERSION ON Q ATOMS
  if(qm->bCLUSTER)
  for(k=0;k<qm->nrQMatoms;k++){         /* ADD Q ATOMS                   */
    a=qm->QMcluster[k*2];               /* CLUSTER INDEX a:(0 -> NC-1)   */
    b=qm->QMcluster[k*2+1];             /* QM INDEX      b:(0 -> nrQM-1) */
    for(l=0;l<qm->nrQMatoms;l++){       /* DISTANCE BETWEEN POINT AND THE CLUSTER THIS QM ATOM BELONGS TO */
      aa=qm->QMcluster[l*2];
      bb=qm->QMcluster[l*2+1];
      if((qm->QMcA[b*(qm->nrQMatoms)+bb]<=(qm->fQMelrc)) && (a!=aa) ) {	// IS THIS A QM ATOM THAT BELONG TO A DIFFERENT CLUSTER ?
        for(j=XX;j<=ZZ;j++){
          QME0=qm->QMcA[b*(qm->nrQMatoms)+bb];
          QME0*=QME0; 
          rij=QME0;
          QME0*=QME0;
          rij*=QME0;
          //QMgrad[b][j]+=6.0*sqrt(qm->c6[b]*qm->c6[bb])*((qm->rQM[b][j])-(qm->rQM[bb][j]))/rij/HARTREE_BOHR2MD;
          f[b][j]+=6.0*sqrt(qm->c6[b]*qm->c6[bb])*((qm->rQM[b][j])-(qm->rQM[bb][j]))/rij;
        }
      }
    }
  }
  //END DISP.

  /*
   * put the QMMM forces in the force array and to the fshift
   */
  for(i=0;i<qm->nrQMatoms;i++){
    for(j=0;j<DIM;j++){
      f[i][j]      += HARTREE_BOHR2MD*QMgrad[i][j];
      fshift[i][j]  = HARTREE_BOHR2MD*QMgrad[i][j];
    }
  }
  for(i=0;i<mm->nrMMatoms;i++){
    for(j=0;j<DIM;j++){
      f[i+qm->nrQMatoms][j]      = HARTREE_BOHR2MD*MMgrad[i][j];      
      fshift[i+qm->nrQMatoms][j] = HARTREE_BOHR2MD*MMgrad[i][j];
    }
  }
  QMener = QMener*HARTREE2KJ*AVOGADRO;
  step++;

  /* Deallocate */
  free(exe);
  free(QMgrad);
  free(MMgrad);

  return(QMener);

} /* call_gaussian */

real call_gaussian_SH(t_commrec *cr, t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, 
		      rvec f[], rvec fshift[])
{ 
  /* a gaussian call routine intended for doing diabatic surface
   * "sliding". See the manual for the theoretical background of this
   * TSH method.  
   */
  static int
    step=0;
  int
    state,i,j;
  real
    QMener=0.0;
  static  bool
    swapped=FALSE; /* handle for identifying the current PES */
  bool
    swap=FALSE; /* the actual swap */
  rvec
    *QMgrad,*MMgrad;
  char
    *buf;
  char
    *exe;
  
  snew(exe,300);
  sprintf(exe,"%s/%s",qm->gauss_dir,qm->gauss_exe);
  /* hack to do ground state simulations */
  if(!step){
    snew(buf,20);
    buf = getenv("STATE");
    if (buf)
      sscanf(buf,"%d",&state);
    else
      state=2;
    if(state==1)
      swapped=TRUE;
  }
  /* end of hack */

  /* copy the QMMMrec pointer */
  snew(QMgrad,qm->nrQMatoms);
  snew(MMgrad,mm->nrMMatoms);
  /* at step 0 there should be no SA */
  /*  if(!step)
   * qr->bSA=FALSE;*/
  /* temporray set to step + 1, since there is a chk start */
  write_gaussian_SH_input(step,swapped,fr,qm,mm);

  do_gaussian(step,exe,0,qm);
  QMener = read_gaussian_SH_output(QMgrad,MMgrad,step,swapped,qm,mm);

  /* check for a surface hop. Only possible if we were already state
   * averaging.*/
  if(qm->SAstep>0){
    if(!swapped){
      swap    = (step && hop(step,qm));
      swapped = swap;
    } 
    else { /* already on the other surface, so check if we go back */
      swap    = (step && hop(step,qm));
      swapped = !swap; /* so swapped shoud be false again */
    }
    if (swap){/* change surface, so do another call */
      write_gaussian_SH_input(step,swapped,fr,qm,mm);
      do_gaussian(step,exe,0,qm);
      QMener = read_gaussian_SH_output(QMgrad,MMgrad,step,swapped,qm,mm);
    }
  }
  /* add the QMMM forces to the gmx force array and fshift
   */
  for(i=0;i<qm->nrQMatoms;i++){
    for(j=0;j<DIM;j++){
      f[i][j]      = HARTREE_BOHR2MD*QMgrad[i][j];
      fshift[i][j] = HARTREE_BOHR2MD*QMgrad[i][j];
    }
  }
  for(i=0;i<mm->nrMMatoms;i++){
    for(j=0;j<DIM;j++){
      f[i+qm->nrQMatoms][j]      = HARTREE_BOHR2MD*MMgrad[i][j];
      fshift[i+qm->nrQMatoms][j] = HARTREE_BOHR2MD*MMgrad[i][j];
    }
  }
  QMener = QMener*HARTREE2KJ*AVOGADRO;
  fprintf(stderr,"step %5d, SA = %5d, swap = %5d\n",
	  step,(qm->SAstep>0),swapped);
  step++;
  free(exe);
  return(QMener);

} 
/* call_gaussian_SH */    
/* end of gaussian sub routines */
