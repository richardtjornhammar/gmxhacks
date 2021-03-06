/*
 * $Id: qm_mpqc.c,v 1.4.2.1 2009/07/19 23:06:26 Tjornhammar Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.7
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
 *
 * And Hey Hey: !!NEW IMPLEMENTATIONS FOR!!
 * CLKS, HSOSKS, UKS, CLHF, HSOSHF, UHF, and MBPT2
 * 3-21G, 6-31G, 6-31G*, 6-31G**, STO-3G, STO-6G
 *
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
#include "atomNames.h"
#include <stdlib.h>
#include "cluster.h"
#include <errno.h>

#define MPQCLEN 2000
#define   sgn(a)        ((a)<0?(-1):(1))
#define   max(a,b)     ((a)>(b)?(a):(b))
double mround(double);

void init_mpqc(t_commrec *cr, t_QMrec *qm, t_MMrec *mm)
{
  char *buf,*p;
  int i,n;

  snew(buf,MPQCLEN);
  buf = getenv("MPQC_CPUS");
  if (buf){
    sscanf(buf,"%d",&qm->nQMcpus);
  }
  else{
    qm->nQMcpus=1;
  }

  buf = getenv("MPQC_NODES");
  if (buf){
    sscanf(buf,"%d",&qm->nQMnodes);
  }
  else{
    qm->nQMnodes=1;
  }

  i=0;
  buf=getenv("QM_EXE");
  if (buf){
    snew(qm->gauss_exe,MPQCLEN);
    sscanf(buf,"%s",qm->gauss_exe);
  }
  else{
    gmx_fatal(FARGS,"MUST >> export QM_EXE=[An executable with fully qualified path]\n");
  }
  fprintf(stderr,"\nMPQC initialised...\n");
}

void write_mpqc(t_forcerec *fr, t_QMrec *qm, t_MMrec *mm)
{
  // NOTE:: MPQC HAS OWN INPUT CHECKING FUNCTIONS
  char  wstring[MAXEL],dstring[MAXEL];
  char  tstring[MAXEL],estring[MAXEL];
  int   i,i0,j,qcflag=0,NRC=0,k=0,l=0,T,q=0,qq=0;
  int   nrZ[MAXFILES],nrQ[MAXFILES],nrCOL=0,nrMM=0,nrQM=0,a,b;
  double  tchrg[MAXFILES];
  FILE  *clustF[MAXFILES]; //MAX 1000 FILES CAN BE GENERATED!!!!  
  
  strcpy(tstring,"run_mpqc");
  strcpy(dstring,"run_mpqc");
  strcpy(estring,".inp");
  NRC=qm->nrQMc;   nrQM=qm->nrQMatoms;
  nrMM=mm->nrMMatoms; nrCOL=nrMM+nrQM;

  if(NRC>1000)
  {
    gmx_fatal(FARGS,"To many clusters. Your QM system is probably exploding");
  }

  for(i=0;i<NRC;i++){
    tchrg[i]=0;
    sprintf(dstring,"run_mpqc%d",i);
    clustF[i]=fopen(strcat(dstring,estring),"w+");
    fprintf(clustF[i],"%s","molecule<Molecule>:(\n unit = bohr");
    nrZ[i]=0;nrQ[i]=0;
  }

  /*
 ======>include_q::
Some of the atoms can be specified as Q and given a customizable charge. Such atoms are a point charge that do not have basis functions. If this option is true, then the Q atoms are included when computing the nuclear charge and the electric field due to the nuclear charge.
 ======>include_qq::
Some of the atoms can be specified as Q and given a customizable charge. Such atoms are a point charge that do not have basis functions. If this option is true, then the Q atoms are included when computing the nuclear repulsion energy and its derivatives.
  */

  if(nrMM >= 1 && (qm->includeQ>eIncludeNONE)) {
    if(qm->includeQ>=eIncludeQ)
      q=1;
    if(qm->includeQ==eIncludeQQ)
      qq=1;
    for(i=0;i<NRC;i++)
      fprintf(clustF[i],"\n include_q=%d\n include_qq=%d\n charge = [",q,qq);
    for(i=0;i<nrMM;i++){
      for(j=0;j<NRC;j++){
	if(mm->MMcluster[nrCOL*j + i])
	  fprintf(clustF[j],"% 4.2f ",mm->MMcharges[i]);
      }
    }
    for(i=0;i<nrQM;i++){
      a=qm->QMcluster[i*2];      /* CLUSTER INDEX a:(0 -> NC-1)   */
      b=qm->QMcluster[i*2+1];    /* QM INDEX      b:(0 -> nrQM-1) */
      for(j=0;j<NRC;j++){
	if(mm->MMcluster[nrCOL*j+b+nrMM])
	  fprintf(clustF[j],"% 4.2f ",qm->pcharges[b]);
      }      
    }
    for(i=0;i<NRC;i++)
      fprintf(clustF[i],"%s","]");
    qcflag=1;
  }

  for(i=0;i<NRC;i++){
    fprintf(clustF[i],"%s","\n{ atoms geometry } = {");
  }
  if(qcflag) {
    for(i=0;i<nrMM;i++){
      for(j=0;j<NRC;j++){
	if(mm->MMcluster[ nrCOL*j + i ]){
	  fprintf(clustF[j],"\n Q [% 6.3f % 6.3f % 6.3f]",mm->rMM[ nrCOL*j + i ][XX]/BORH2NM,mm->rMM[ nrCOL*j + i ][YY]/BORH2NM,mm->rMM[ nrCOL*j + i ][ZZ]/BORH2NM);
	  nrQ[j]++;
	}
      }
    }
    for(i=0;i<nrQM;i++){
      a=qm->QMcluster[i*2];      /* CLUSTER INDEX a:(0 -> NC-1)    */
      b=qm->QMcluster[i*2+1];    /* QM INDEX      b:(0 -> nrQM-1)  */
      for(j=0;j<NRC;j++){
	if(mm->MMcluster[b+nrMM+nrCOL*j]){
	  fprintf(clustF[j],"\n Q [% 6.3f % 6.3f % 6.3f]",mm->rMM[ b+nrMM+nrCOL*j ][XX]/BORH2NM,mm->rMM[ b+nrMM+nrCOL*j ][YY]/BORH2NM,mm->rMM[ b+nrMM+nrCOL*j ][ZZ]/BORH2NM);
	  nrQ[j]++;
	}
      }
    } 
  }
  for(i=0;i<nrQM;i++){
    j=qm->QMcluster[i*2];k=qm->QMcluster[i*2+1];
    tchrg[j]+=qm->pcharges[k];
    if(qm->atomicnumberQM[k]==0)
        fprintf(clustF[j],"\n H [% 6.3f % 6.3f % 6.3f]",qm->rQM[k][XX]/BORH2NM,qm->rQM[k][YY]/BORH2NM,qm->rQM[k][ZZ]/BORH2NM);
    else
        fprintf(clustF[j],"\n%s [% 6.3f % 6.3f % 6.3f]",getAtomName(qm->atomicnumberQM[k]),qm->rQM[k][XX]/BORH2NM,qm->rQM[k][YY]/BORH2NM,qm->rQM[k][ZZ]/BORH2NM);
    nrZ[j]++;
  }
  for(i=0;i<NRC;i++){
    fprintf(clustF[i],"%s","\n }\n )\n basis<GaussianBasisSet>: (");
    fprintf(clustF[i],"\n name = \"%s\"\n",eQMbasis_names[qm->QMbasis]);
    fprintf(clustF[i],"%s","molecule = $:molecule\n ) \n mpqc: (\n mole<");
    fprintf(clustF[i],"%s",eQMmethod_names[qm->QMmethod]);
    fprintf(clustF[i],"%s",">: (\n molecule = $:molecule \nbasis = $:basis");

    switch(qm->QMmethod) {
    case eQMmethodUKS:fprintf(clustF[i],"%s","\n functional<StdDenFunctional>: ( name = \"B3LYP\" )"); break;
    case eQMmethodCLKS:fprintf(clustF[i],"%s","\n functional<StdDenFunctional>: ( name = \"B3LYP\" )"); break;
    case eQMmethodHSOSKS:fprintf(clustF[i],"%s","\n functional<StdDenFunctional>: ( name = \"B3LYP\" )"); break;
    }
    fprintf(clustF[i], "\n total_charge = %lf \n )\n do_gradient = yes", mround(tchrg[i]));
    fprintf(clustF[i],"%s","\n)\n");
    
  }
  
  for(i=0;i<NRC;i++){
    fclose(clustF[i]);
    if(nrZ[i]==0){ 
      /* if we, during the cluster reduction phase, have generated trailing 
	 empty QM files remove them here */
	qm->nrQMc--;
    }
  }
}

real read_mpqc(rvec QMgrad[],rvec MMgrad[], t_QMrec *qm, t_MMrec *mm, int I)
{
  int  i,j,atnum,HL=0,nQMc,logi,nrHERE=0,nrMM,nrQM,nrCOL,leave=0,once=1; 
  char buf[MPQCLEN],atchar[2];
  real QMener=0.0,fei1,fei2,deltaE=0.0,tol=1;
  FILE *out,*mpqcerr,*mpqcMM,*mpqcQA;
  char txt1[20],txt2[20],txt3[20],txt4[20],txt5[20],filenmerr[50],filempqc[50],fileqa[50];
  char txt6[20],txt7[20],txt8[20];
  char filenm[100],gradtxt1[20],gradtxt2[20],enetxt1[20],enetxt2[20],enetxt3[20],enetxt4[20],hline[50],fline[MPQCLEN];
  char grdX[50],grdY[50],grdZ[50],eneQM[50];
  
  sprintf(gradtxt1,"Total"); sprintf(gradtxt2,"Gradient:");
  sprintf(enetxt1,"Value"); sprintf(enetxt2,"of"); sprintf(enetxt3,"the"); sprintf(enetxt4,"MolecularEnergy:");
  sprintf(filenm,"run_mpqc%d.out",I);
  sprintf(fileqa,"run_mpqc%d.qa",I);
  sprintf(filenmerr,"run_mpqc%d.inf",I);
  sprintf(filempqc,"run_mpqc%d.mm",I);
  
  out=fopen(filenm,"r");
  mpqcQA=fopen(fileqa,"w+");
  mpqcMM=fopen(filempqc,"w");
  mpqcerr=fopen(filenmerr,"w");
  
  nQMc=qm->nrQMc;
  nrMM=mm->nrMMatoms;
  nrQM=qm->nrQMatoms;
  nrCOL=nrMM+nrQM;
  
  fprintf(mpqcQA,"CL QM QA --->\n");
  for(i=0;i<nrQM;i++){
    fprintf(mpqcQA,"%d %d",qm->QMcluster[2*i],qm->QMcluster[2*i+1]);
    for(j=0;j<nQMc;j++){
      fprintf(mpqcQA," %d",mm->MMcluster[ ( qm->QMcluster[2*i+1]) + nrMM +j*nrCOL] );
    }
    fprintf(mpqcQA,"\n");
  }
  
  j=0;
  while(!feof(out)) {
    j++;
    if( NULL == fgets(fline,MPQCLEN,out)) {               /* Error handling */
      leave=1;
      break;
    }
    else{
      sscanf(fline," %s %s %s %s %s %s %s %s\n",txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);
      if( !strcmp(txt1,"iter")) { /* Assigns energy from each interation */
	QMener=atof(txt5);
	deltaE=atof(txt8);
      }
      if( !strcmp(txt1,gradtxt1) & !strcmp(txt2,gradtxt2) & once) {
	nrHERE=0;once=0;
        fprintf(mpqcerr,"FOUND GRADIENT SECTION!\n");
	if( nrMM>0 & qm->includeQ>eIncludeNONE){
	  for(i=0;i<nrMM;i++){
	    logi=(int)mm->MMcluster[I*nrCOL+i]; /* BOOLEAN (IS THIS MM ATOM PART OF THIS CLUSTER) */
	    if(!((logi==0)||(logi==1)))
	      gmx_fatal(FARGS,"ERROR IN LOGIC (THIS SHOULD NEVER HAPPEN)");
	    if(logi){
	      HL++; nrHERE++;
              fprintf(mpqcMM,"MM ATOM %d BELONGS TO CLUSTER %d\nMM GOT ATOM %d OUT OF %d POSSIBLE\n",i,I,nrHERE,nrMM);
	      if(NULL == fgets(fline,MPQCLEN,out))
		{
		  fprintf(mpqcerr,"MM GOT ATOM %d OUT OF %d POSSIBLE\n",i,mm->nrMMatoms);
		  gmx_fatal(FARGS,"Error reading MM MPQC output");
		} 
              fprintf(mpqcMM,"%s",fline);
	      sscanf(fline,"%s %s %s %s %s\n",atchar,hline,grdX,grdY,grdZ);
	      fprintf(mpqcMM,"MM GOT GRADIENT FOR ATOM::%s",fline);
	      
	      MMgrad[i][XX]+=atof(grdX);
	      MMgrad[i][YY]+=atof(grdY);
	      MMgrad[i][ZZ]+=atof(grdZ);
	      
	      fprintf(mpqcMM,"%f %f %f\n",MMgrad[i][XX],MMgrad[i][YY],MMgrad[i][ZZ]);
	      fprintf(mpqcMM,"%s\n\n",hline);
	    }
	  }
	  //REMOVE ALL EXTRA Q ATOMS HERE
	  nrHERE=0;
	  for(i=0;i<nrQM;i++){                
	    j=qm->QMcluster[2*i+1];
	    if(mm->MMcluster[(j+nrMM)+I*nrCOL]){ // IF WE HAVE QA ATOMS THEN REMOVE THOSE LINES
              fprintf(mpqcerr,"QA ATOM %d BELONGS TO CLUSTER %d\n",j,I);
	      if(NULL == fgets(fline,MPQCLEN,out))
		{
		  fprintf(mpqcerr,"QA nr %d :: GOT QM ATOM %d OUT OF %d POSSIBLE\n",++nrHERE,j,nrQM);
		  gmx_fatal(FARGS,"Error reading QA MPQC output");
		}
	      HL++;
	      fprintf(mpqcerr,"%s\n",fline);
	      sscanf(fline,"%s %s %s %s %s\n",atchar,hline,grdX,grdY,grdZ);
	      //WE NOW HAVE THEM BUT DONT ASSIGN THEM, WE COULD SINCE THEY ARE JUST FORCES
	    }
	  } 
	  //END EXCESS QM Q
	}
	for(i=0;i<nrQM;i++){
	  if(qm->QMcluster[2*i]==I){
	    j=qm->QMcluster[2*i+1];//QM ATOM INDEX
	    if(NULL == fgets(fline,MPQCLEN,out))
	      {
		mpqcerr=fopen(filenmerr,"w");	      
		fprintf(mpqcerr,"QM GOT ATOM %d OUT OF %d POSSIBLE\n",i,qm->nrQMatoms);
		gmx_fatal(FARGS,"Error reading QM MPQC output");
	      }
	    sscanf(fline,"%s %s %s %s %s\n",atchar,hline,grdX,grdY,grdZ);
	    if(hline[1]=='\0'){
	      hline[1]=hline[0];
	      hline[0]=' ';
	      hline[2]='\0';
	    }
	    //CHECK THAT ATOMTYPE IS CORRECT
	    if( !(( !strcmp(" H",hline) & !strcmp("LA",getAtomName(qm->atomicnumberQM[j])) ) | !strcmp(hline,getAtomName(qm->atomicnumberQM[j])) )){ 
	      fprintf(stderr,"ATOM::%s ENTRY::%s",getAtomName(qm->atomicnumberQM[j]),hline);
	      fprintf(mpqcerr,"ATOM::%s ENTRY::%s",getAtomName(qm->atomicnumberQM[j]),hline);
	      gmx_fatal(FARGS,"MISSMATCH OF QM ATOMS (THIS SHOULD NEVER HAPPEN)");
	    }
	    HL++;
	    QMgrad[j][XX]+=atof(grdX);
	    QMgrad[j][YY]+=atof(grdY);
	    QMgrad[j][ZZ]+=atof(grdZ);
	  }
	}
      }
      if(!strcmp(txt1,enetxt1) & !strcmp(txt2,enetxt2) & !strcmp(txt3,enetxt3) & !strcmp(txt4,enetxt4) ) { 
	fprintf(mpqcerr,"FOUND THE MOLECULAR ENERGY SECTION\n");
	sscanf(fline,"%s %s %s %s %s",hline,grdX,grdY,grdZ,eneQM);
	QMener=atof(eneQM);
	fprintf(mpqcerr,"--==::< QMener:: %f >::==--\n",QMener);
	break;
      }
    }
  }
  
  if(leave){
    //IF WE HAD TO LEAVE BEFORE THE MOLECULAR ENERGY THEN WE MUST HAVE GOTTEN THE GRADIENTS
    for(i=0;i<nrCOL;i++)
      if((mm->MMcluster[I*nrCOL+i]) &&  qm->includeQ>eIncludeNONE)
	HL--;
    for(i=0;i<qm->nrQMatoms;i++)
      if(qm->QMcluster[2*i]==I)
	HL--;
    fprintf(mpqcerr,"\n!!QM DIDNT PROPERLY CONVERGE!!\nEXTRACTED ENERGY AND GRADIENT::\nENERGY=%10.6f TOL=%10.6e REM=%d\n",QMener,deltaE,HL);
    fprintf(stderr ,"\n!!QM DIDNT PROPERLY CONVERGE!!\nEXTRACTED ENERGY AND GRADIENT::\nENERGY=%10.6f TOL=%10.6e REM=%d\n",QMener,deltaE,HL);
    if(deltaE>=1e-4){/* EXTREMLY OPTIMISTIC ABOUT WHAT TO INCLUDE SO BAD SIMULATIONS SHOWS IN STATISTICS */
      gmx_fatal(FARGS,"QM FAILED TO CONVERGE :: WAVEFUNCTION/ENERGY TOL=%f",deltaE);
    }
  }
  fprintf(mpqcerr,"======| QMener:: %f \n",QMener);
  
  fclose(out);
  fclose(mpqcQA);
  fclose(mpqcMM);
  fclose(mpqcerr);
  
  return(QMener); 
}

real call_mpqc(t_commrec *cr,  t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[])
{
  /* normal mpqc jobs */
  int i,j,RV=0,cl=1,NFILES=1,NFILER=1;
  real QMener=0.0;
  rvec *QMgrad,*MMgrad;
  char *exe;
  real *qmclusters,*mmclusters;
  FILE *stdver;

  /* Allocate */
  snew(QMgrad,qm->nrQMatoms);
  snew(MMgrad,mm->nrMMatoms);

  //init
  for(i=0;i<qm->nrQMatoms;i++){
    QMgrad[i][XX]=0.0; QMgrad[i][YY]=0.0; QMgrad[i][ZZ]=0.0;
    qm->rQM[i][XX]=qm->xQM[i][XX];
    qm->rQM[i][YY]=qm->xQM[i][YY];
    qm->rQM[i][ZZ]=qm->xQM[i][ZZ];
  }
  for(i=0;i<mm->nrMMatoms;i++){
    MMgrad[i][XX]=0.0; MMgrad[i][YY]=0.0; MMgrad[i][ZZ]=0.0;
    mm->rMM[i][XX]=mm->xMM[i][XX];
    mm->rMM[i][YY]=mm->xMM[i][YY];
    mm->rMM[i][ZZ]=mm->xMM[i][ZZ];
  }
  
  if(qm->bCLUSTER){
    clustering(qm);  
    clusterMM(mm,qm);
    NFILES=qm->nrQMc; 
  }
  else{
    qm->nrQMc=1;NFILES=1;
    one_cluster(mm,qm);
  }

  write_mpqc(fr,qm,mm);

  if(NFILES<=0 || NFILES>MAXFILES){ /* This will only happen if the clustering has a bug */
    gmx_fatal(FARGS,"CLUSTERED MPQC FAILED");
  }

  /* BEGIN Big loop */
  for(i=0;i<NFILES;i++){
    RV=do_mpqc(qm,i);
    fflush(stdout); fflush(stdin); 
    QMener += read_mpqc(QMgrad,MMgrad,qm,mm,i);
  }
  /* END Big loop */
  /* put the QMMM forces in the force array and do the fshift */
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
  return(QMener);
}

int do_mpqc(t_QMrec *qm, int I)
{
  //NON PIPED IMPLEMENTATION OF CALL
  pid_t pid;
  int   rv,step_id;
  int	commpipe[2],stdcopy[2];		   /* This holds the fd for the in & out of the pipe */
  FILE  *in,*out;
  char  input_string[MPQCLEN],exename[200],execute[10],eo0[20],eo1[20],eo2[20],eo3[20],eo4[20],exe_list[100];

  //sprintf(eo0,"-o=run_mpqc%d.out",I);
  sprintf(eo1,"--nnodeperjob=%d",qm->nQMnodes); 
  sprintf(eo2,"--nthreadperproc=%d",qm->nQMcpus);
  sprintf(eo3,"--autoout");
  sprintf(eo4,"--simpout");
  sprintf(exename,"run_mpqc%d.inp",I);

  sprintf(input_string,"%s %s %s %s %s %s",qm->gauss_exe,eo1,eo2,eo3,eo4,exename);

#ifdef GMX_NO_SYSTEM
  printf("Warning-- No calls to system(3) supported on this platform.");
  gmx_fatal(FARGS,"Call to '%s' failed\n",input_string);
#else
  if ( system(input_string) != 0 )
    gmx_fatal(FARGS,"Call to '%s' failed\n",input_string);
#endif
  
  return(0);
}

void printClustering(t_QMrec *qm,FILE *fpo)
{
  int i,j,k,N,C=-1;

  N=qm->nrQMatoms;

  for(i=0;i<N;i++){
    if(!(C==(qm->QMcluster[i*2]))){
      C=qm->QMcluster[i*2];
      fprintf(fpo,"\n % POSITION::\n\%");
      for(j=XX;j<=ZZ;j++)
	fprintf(fpo,"% 9.5f ",qm->QMcluspos[qm->QMcluster[i*2]][j]);
    }
    fprintf(fpo,"\n %5d %5d",qm->QMcluster[i*2],qm->QMcluster[i*2+1]);
  }
}

void sortClusters(t_QMrec *qm)
{
 int i,j,k,l,m,nrEntries;

 int sa,sb;
 int ft,ft2;
 nrEntries=qm->nrQMatoms;

 for(sa = 1; sa < nrEntries; ++sa)
   for(sb = nrEntries-1; sb >= sa; --sb) {
     // compare adjacent elements 
     if(qm->QMcluster[ (sb - 1)*2] < qm->QMcluster[ sb*2 ]) {
       // exchange elements 
       ft = qm->QMcluster[ (sb - 1)*2 ];
       ft2= qm->QMcluster[ (sb - 1)*2+1];
       qm->QMcluster[ (sb - 1)*2 ] = qm->QMcluster[ sb*2 ];
       qm->QMcluster[ (sb - 1)*2 +1 ] = qm->QMcluster[ sb*2 + 1 ];
       qm->QMcluster[ sb*2 ] = ft;
       qm->QMcluster[ sb*2 + 1 ] = ft2;
     }
   }
}

void qmmmPBC(t_QMrec *qm, rvec dist[2])
{
  int i,j,k,l,m,n;

  //BEGIN PBC CLUSTER MODIFICATION
  dist[1][XX]=0.0; dist[1][YY]=0.0; dist[1][ZZ]=0.0;
  if(qm->ePBC==epbcXY){
    for(l=XX;l<=YY;l++){
      if(dist[0][l]>(qm->box[l][l])*0.5){
	dist[0][l]-=qm->box[l][l];
	dist[1][l]=-1;
      }
      else 
	if(abs(dist[0][l])>(qm->box[l][l])*0.5){
	  dist[0][l]+=qm->box[l][l];
	  dist[1][l]=1;
	}
    }
  }
  if(qm->ePBC==epbcXYZ){
    for(l=XX;l<=ZZ;l++){
      if(dist[0][l]>(qm->box[l][l])*0.5){
	dist[0][l]-=qm->box[l][l];
	dist[1][l]=-1;
      }
      else 
	if(abs(dist[0][l])>(qm->box[l][l])*0.5){
	  dist[0][l]+=qm->box[l][l];
	  dist[1][l]=1;
	}
    }	
  }
}

void clusterAlgorithm(t_QMrec *qm,FILE *fpo)
{
  int i,j,k,l,m,N;
  int flag,C,oC;

  N=qm->nrQMatoms;
  //	  qm->QMcluster[k]=C;
  //	  qm->QMcluster[k+1]=i;
  for(i=0;i<N;i++){
    qm->QMcluster[2*i]=i+1; //index
    qm->QMcV[i]=0;          //visited
  }
  C=-1;
  for(i=0;i<N;i++){
    flag=1;
    for(j=0;j<N;j++){
      if( !(qm->QMcV[j]) ){
	if( qm->QMcB[i*N+j] ){
	  qm->QMcV[j]=1;
	  qm->QMcV[i]=1;
	  if((qm->QMcluster[2*i])>0)
	    qm->QMcluster[2*j]=C;
	  else
	    qm->QMcluster[2*j]=(qm->QMcluster[2*i]);
	}
      }else{
	if( qm->QMcB[i*N+j] )
	  if(qm->QMcluster[2*j]<0)
	    qm->QMcluster[2*i]=qm->QMcluster[2*j];
      }
      if( !(qm->QMcB[i*N+j]) && flag){
	C--; flag=0;
      }
    }
  }
  for(j=(N-1);j>=0;j--){
    for(i=(N-1);i>=0;i--){
      if((qm->QMcluster[2*i])!=(qm->QMcluster[2*j])& (qm->QMcB[i*N+j])){
	(qm->QMcluster[2*i])=(qm->QMcluster[2*j]);
      }
    }
    (qm->QMcluster[2*j+1])=j;
  }

  //MAKE CLUSTERS SEQUENTIAL AND INIT CLUSTER POSITION VECTORS
  C=0; oC=(int)(qm->QMcluster[0]); qm->QMcV[0]=oC;
  for(i=0;i<N;i++){
    if(qm->QMcluster[2*i]<oC){
      for(j=XX;j<=ZZ;j++){
	qm->QMcluspos[C][j]=0.0;
      }
      C++;oC=qm->QMcluster[2*i];
      qm->QMcV[C]=oC;
      for(j=XX;j<=ZZ;j++){
	qm->QMcluspos[C][j]=0.0;
      }
    }
  }
  sortClusters(qm);

  fprintf(fpo,"\n\n");
  for(j=0;j<=C;j++){
    for(i=0;i<N;i++){
      if(qm->QMcluster[2*i] == qm->QMcV[j])
	qm->QMcluster[2*i]=j;
      if(j == C)
	fprintf(fpo,"Index::% 4d Cluster::% 4d \n",qm->QMcluster[2*i+1],(qm->QMcluster[2*i]));
    }
  }
  qm->nrQMc=C+1;
}

void makeDistance(t_QMrec *qm, FILE *fpo)
{
  int   i,j,k,l,N,C,oC,flag,NP;
  float RL=0.3;
  rvec  dist[2];

  N=qm->nrQMatoms;
  fprintf(fpo,"%Distance matrix from gromacs\nN=%d; A=[\n",N);
  for(i=0;i<N;i++){
    qm->QMcV[i]=i+1;
    for(j=0;j<N;j++){
      dist[0][XX]=qm->xQM[i][XX]-qm->xQM[j][XX];
      dist[0][YY]=qm->xQM[i][YY]-qm->xQM[j][YY];
      dist[0][ZZ]=qm->xQM[i][ZZ]-qm->xQM[j][ZZ];

      qmmmPBC(qm,dist);

      qm->QMcA[i*N+j]=sqrt( pow(dist[0][XX],2.0) + pow(dist[0][YY],2.0) + pow(dist[0][ZZ],2.0) );
      qm->QMcB[i*N+j]=(qm->QMcA[i*N+j]<=RL)?1.0:0.0;
      
      if(qm->QMcB[i*N+j]){
	for(k=XX;k<=ZZ;k++){
	  qm->rQM[i][k]=qm->xQM[i][k]+dist[1][k]*(qm->box[k][k]);
	}
      }
      fprintf(fpo," % 9.5f ",qm->QMcA[i*N+j]);
    }
    fprintf(fpo,"\n");
  }
  fprintf(fpo,"\n];\n");
  
  fprintf(fpo,"%Boolean matrix from gromacs\n B=[\n",N);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      fprintf(fpo," %d",qm->QMcB[i*N+j]);
    }
    fprintf(fpo,"\n");
  }
  fprintf(fpo,"];\n\n");
  
  for(i=0;i<N;i++){
    fprintf(fpo,"% 5d ",qm->QMcV[i]);
  }

  fprintf(fpo,"\n\n");
}

void clusterPositions(t_QMrec *qm)
{
  int   i,j,k,l,N,C,oC,flag,NP;
  N=qm->nrQMatoms;

  oC=0; NP=0;
  for(i=0;i<N;i++){
    if(oC-qm->QMcluster[i*2]==0){
      NP++;
      for(j=XX;j<=ZZ;j++){
	qm->QMcluspos[oC][j]+=qm->rQM[qm->QMcluster[2*i+1]][j];
      }
    }else{
      for(j=XX;j<=ZZ;j++){
	qm->QMcluspos[oC][j]=(qm->QMcluspos[oC][j])/((float)NP);
      }
      NP=1;
      oC=qm->QMcluster[i*2];
      if( i== (N-1) ){
	for(j=XX;j<=ZZ;j++){
	  qm->QMcluspos[oC][j]+=qm->rQM[qm->QMcluster[2*i+1]][j];
	}
      }
    }
  }
}

void clustering(t_QMrec *qm)
{
  int   i,j,k,l,N,C,oC,flag,NP;
  float RL=0.3;
  rvec  dist[2];
  FILE  *fpo;

  fpo=fopen("outpclust.m","w");

  N=qm->nrQMatoms;

  makeDistance(qm,fpo);     //CREATE DISTANCE MATRIX

  clusterAlgorithm(qm,fpo); //DO CLUSTER SEARCH

  C=qm->nrQMc;

  if(C>N){
    gmx_fatal(FARGS,"IMPOSSIBLE CLUSTERING SITUATION");
  }

  fprintf(fpo,"\nFOUND %6d CLUSTER(S)\n\n",C);
  
  clusterPositions(qm);
  printClustering(qm,fpo);

  fclose(fpo);
}

void one_cluster(t_MMrec *mm, t_QMrec *qm)
{
  int i,j,k;

  for(i=0;i<qm->nrQMatoms;i++){ //QM
    qm->QMcluster[i*2]=0;
    qm->QMcluster[i*2+1]=i;
  }
  for(i=0;i<mm->nrMMatoms;i++){ //MM
    mm->MMcluster[i]=1;
  }
  for(i=0;i<qm->nrQMatoms;i++){ //QA
    mm->MMcluster[i+mm->nrMMatoms]=0;
  }
}

void clusterMM(t_MMrec *mm, t_QMrec *qm)
{
  //SIMPLE DISTANCE SEARCH FOR CLOSEST CLUSTER
  int a,b,nrQMc=0,i,j,k,kk,l,nrQMa=0,nrMM,nrCOL,aa,bb;
  rvec R[2];
  double D[3],D0=100000000.0,RC=1.0,R0,D1=0.5;
  /*
    Motivation for RC choice 
    RC=0.8 is ca 200 atoms in spherical volume around an atom
    Strong polarisation effects occur within ca 1nm of an ionic center 
    If the charge is not extreme
  */

  nrQMc=qm->nrQMc; 
  nrQMa=qm->nrQMatoms;
  nrMM=mm->nrMMatoms; 
  nrCOL=nrMM+nrQMa;

  for(j=0;j<nrMM;j++){
    D0=100000000.0; k=0;
    for(i=0;i<nrQMa;i++){  /* CALCULATE DISTANCE BETWEEN POINT AND ALL QM ATOM CLUSTERS */
      a=qm->QMcluster[i*2];
      b=qm->QMcluster[i*2+1];
      R[0][XX]=mm->xMM[j][XX]-qm->rQM[i][XX];
      R[0][YY]=mm->xMM[j][YY]-qm->rQM[i][YY];
      R[0][ZZ]=mm->xMM[j][ZZ]-qm->rQM[i][ZZ];
      qmmmPBC(mm,R);
      R0=sqrt(pow(R[0][XX],2.0)+pow(R[0][YY],2.0)+pow(R[0][ZZ],2.0));
      for(kk=XX;kk<=ZZ;kk++){
	mm->rMM[j+a*nrCOL][kk]=mm->xMM[j][kk]+R[1][kk]*(qm->box[kk][kk]);
      }
      l=(R0<RC)?1:0;                      /* CUTOFF SEARCH */
      mm->MMcluster[nrCOL*a+j]     = l;   /* * */
      mm->MMcluster[nrCOL*nrQMc+j]+= l;   /* SUMS UP THE NUMBER OF TIMES AN ATOM IS CALLED */
      if(R0<D0){                          /* SEARCH FOR CLOSEST CLUSTER */
	k=i; D0=R0;
      }
    }
    //mm->MMcluster[nrCOL*nrQMc+j]+=(mm->MMcluster[nrCOL*nrQMc+j]==0)?1:0; 
    //mm->MMcluster[nrCOL*k+j] = 1;/* ENSURES THAT ALL ATOMS ARE PRESENT IN ATLEAST ONE CLUSTER */
  }
  //NOW ASSIGN Q CHARGES FROM THE OTHER CLUSTERS QM ATOMS (DONE)
  for(j=0;j<nrQMa;j++){         /* ADD Q ATOMS */
    a=qm->QMcluster[j*2];       /* CLUSTER INDEX a:(0 -> NC-1) */
    b=qm->QMcluster[j*2+1];     /* QM INDEX      b:(0 -> nrQM-1) */
    for(i=0; i<nrQMa ;i++){     /* DISTANCE BETWEEN POINT AND THE CLUSTER THIS QM ATOM BELONGS TO */
      aa=qm->QMcluster[i*2];
      bb=qm->QMcluster[i*2+1];
      l=((qm->QMcA[b*nrQMa+bb]<=RC) && (a!=aa) )?0:0; //TURNED OFF	
      for(kk=XX;kk<=ZZ;kk++){
	mm->rMM[b+nrMM+aa*nrCOL][kk]=qm->rQM[b][kk];
      }
      mm->MMcluster[(nrCOL*aa+(b+nrMM))] = l;
    }
  }
}

double mround(double T)
{
   double RV=0;
   RV=sgn(T)*ceil(max(T,-1*T));
   return(RV);
}

