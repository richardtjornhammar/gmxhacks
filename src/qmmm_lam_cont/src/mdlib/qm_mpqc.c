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
#include <stdlib.h>
#include "cluster.h"
#include <errno.h>

#define MPQCLEN 10000

void init_mpqc(t_commrec *cr, t_QMrec *qm, t_MMrec *mm)
{
  char *buf,*p;
  int i,n;
  FILE *npalog;

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

  if(qm->bPrint_npa){
    npalog=fopen("npa.log","w");
    fclose(npalog);
  }

  fprintf(stderr,"nr of background atoms = %d\nMPQC initialised...\n",mm->nrMMatoms);
}

void write_mpqc(t_forcerec *fr, t_QMrec *qm, t_MMrec *mm)
{
  // NOTE:: MPQC HAS OWN INPUT CHECKING FUNCTIONS
  char  wstring[MAXEL],dstring[MAXEL],tstring[MAXEL],estring[MAXEL];
  int   i,i0,j,qcflag=0,NRC=0,k=0,l=0,T,q=0,qq=0,zf=1,qi,qj,qk;
  int   nrZ[MAXFILES],nrQ[MAXFILES],nrCOL=0,nrMM=0,nrQM=0,a,b;
  double  tchrg[MAXFILES];
  FILE  *clustF[MAXFILES]; /* MAXFILES=1000 */ 
  FILE  *indexF[MAXFILES];  
  FILE  *allMMndx;

  strcpy(tstring,"run_mpqc");
  strcpy(dstring,"run_mpqc");
  strcpy(estring,".inp");
  NRC=qm->nrQMc;   nrQM=qm->nrQMatoms;
  nrMM=mm->nrMMatoms; nrCOL=nrMM+nrQM;

  if(qm->dI>0){
    NRC=qm->dI;
    qm->dI-=1;
  }

  if(NRC>=MAXFILES)
  {
    gmx_fatal(FARGS,"To many clusters. Your QM system is probably exploding");
  }

  for(i=qm->dI;i<NRC;i++){
    tchrg[i]=0;
    sprintf(dstring,"run_mpqc%d",i);
    clustF[i]=fopen(strcat(dstring,estring),"w");
    for(j=0;j<MAXEL;j++){
       if(dstring[j]=='.')
         break;
    }
    dstring[j+1]='n';dstring[j+2]='d';dstring[j+3]='x';
    indexF[i]=fopen(dstring,"w");
    fprintf(clustF[i],"%s","molecule<Molecule>:(\n unit = bohr\n");
    nrZ[i]=0;nrQ[i]=0;
  }

  if(nrMM >= 1 && (qm->includeQ>eIncludeNONE)) {
    allMMndx=fopen("mpqcMMatoms.ndx","w");
    fprintf(allMMndx,"[ MMINDEX ]\n");
    for(i=0;i<nrMM;i++){
      fprintf(allMMndx,"% 10d", mm->indexMM[i]+1 );  
      nrQ[0]++;
      if(nrQ[0]==10){
        fprintf(allMMndx,"\n");
        nrQ[0]=0;
      }
    }
  }
  /*
    Some of the atoms can be specified as Q and given a customizable charge. Such atoms are a point charge that do not have basis functions.
 ======>includeYES::
 This is the default QMMM option within MPQC
 ======>includeQ::
 If this option is true, then the Q atoms are included when computing the nuclear charge and the electric field due to the nuclear charge.
 ======>includeQQ::
 If this option is true, then the Q atoms are included when computing the nuclear repulsion energy and its derivatives.
 ======>includeQQQ::
 If this option is true, then the Q atoms are included when computing the nuclear repulsion energy and its derivatives and are also are included when computing the nuclear charge and the electric field due to the nuclear charge.
*/

  if(nrMM >= 1 && (qm->includeQ>eIncludeNONE)) {
    if(qm->includeQ==eIncludeQ||qm->includeQ==eIncludeQQQ)
      q=1;
    if(qm->includeQ==eIncludeQQ||qm->includeQ==eIncludeQQQ)
      qq=1;
    for(i=qm->dI;i<NRC;i++){
      if(!(qm->includeQ==eIncludeYES)){ //THIS IS THE MPQC DEFAULT FOR QMMM
        fprintf(clustF[i]," include_q=%d\n include_qq=%d\n",q,qq);
      }
      fprintf(clustF[i]," charge = [");
      fprintf(indexF[i],"[ MMrec_c%d ]\n",i);
    }
    for(i=0;i<nrMM;i++){
      for(j=qm->dI;j<NRC;j++){
	if(mm->MMcluster[nrCOL*j + i]){
	  fprintf(clustF[j],"% 4.2f ",mm->MMcharges[i]);
          fprintf(indexF[j],"% 10d", mm->indexMM[i]+1 );
          nrQ[j]++;
          if(nrQ[j]==10){
            fprintf(indexF[j],"\n");
            nrQ[j]=0;
          }
        }
      }
    }
    for(i=0;i<nrQM;i++){
      a=qm->QMcluster[i*2];      /* CLUSTER INDEX a:(0 -> NC-1)   */
      b=qm->QMcluster[i*2+1];    /* QM INDEX      b:(0 -> nrQM-1) */
      for(j=qm->dI;j<NRC;j++){
	if(mm->MMcluster[nrCOL*j+b+nrMM]){
	  fprintf(clustF[j],"% 4.2f ",qm->pcharges[b]);
          fprintf(indexF[j],"% 10d", qm->indexQM[i]+1 );  //MDRUN INDEX AND NDX INDEX DIFFER WITH 1
          nrQ[j]++;
          if(nrQ[j]==10){
            fprintf(indexF[j],"\n");
            nrQ[j]=0;
          }
        }
      }      
    }
    for(i=qm->dI;i<NRC;i++){
      fprintf(clustF[i],"%s","]");
    }
    qcflag=1;
  }

  for(i=qm->dI;i<NRC;i++){
    nrQ[i]=0;  
    fprintf(indexF[i],"\n\n[ QMrec_c%d ]\n",i);
    fprintf(clustF[i],"%s","\n{ atoms geometry } = {");
  }
  if(qcflag) {
    for(i=0;i<nrMM;i++){
      for(j=qm->dI;j<NRC;j++){
	if(mm->MMcluster[ nrCOL*j + i ]){
	  fprintf(clustF[j],"\n Q [% 6.3f % 6.3f % 6.3f]",mm->rMM[ nrCOL*j + i ][XX]/BORH2NM,mm->rMM[ nrCOL*j + i ][YY]/BORH2NM,mm->rMM[ nrCOL*j + i ][ZZ]/BORH2NM);
	}
      }
    }
    for(i=0;i<nrQM;i++){
      a=qm->QMcluster[i*2];      /* CLUSTER INDEX a:(0 -> NC-1)    */
      b=qm->QMcluster[i*2+1];    /*      QM INDEX b:(0 -> nrQM-1)  */
      for(j=qm->dI;j<NRC;j++){
	if(mm->MMcluster[b+nrMM+nrCOL*j]){
	  fprintf(clustF[j],"\n Q [% 6.3f % 6.3f % 6.3f]",mm->rMM[ b+nrMM+nrCOL*j ][XX]/BORH2NM,mm->rMM[ b+nrMM+nrCOL*j ][YY]/BORH2NM,mm->rMM[ b+nrMM+nrCOL*j ][ZZ]/BORH2NM);
	}
      }
    } 
  }
  for(i=0;i<nrQM;i++){
    j=qm->QMcluster[i*2];k=qm->QMcluster[i*2+1];
    if( !(j >= qm->dI && j < NRC)){
      break;
    }
    tchrg[j]+=qm->pcharges[k];                   //Charge contributions are summed up from the topology definition
    if(qm->atomicnumberQM[k]==0){                 //LA
        fprintf(clustF[j],"\n H [% 6.3f % 6.3f % 6.3f]",qm->rQM[k][XX]/BORH2NM,qm->rQM[k][YY]/BORH2NM,qm->rQM[k][ZZ]/BORH2NM);
    }
    else{
        fprintf(clustF[j],"\n%s [% 6.3f % 6.3f % 6.3f]",getAtomNameMPQC(qm->atomicnumberQM[k]),qm->rQM[k][XX]/BORH2NM,qm->rQM[k][YY]/BORH2NM,qm->rQM[k][ZZ]/BORH2NM);
    }
    if(zf==1){
      zf=qm->atomicnumberQM[k]>=20?0:1;
    }
    nrZ[j]++;

    fprintf(indexF[j],"% 10d", qm->indexQM[k]+1 );
    nrQ[j]++;
    if(nrQ[j]==10){
      fprintf(indexF[j],"\n");
      nrQ[j]=0;
    }
  }
  for(i=qm->dI;i<NRC;i++){
    fprintf(indexF[i],"\n");
    fprintf(clustF[i],"%s","\n }\n )\n basis<GaussianBasisSet>: (");
    fprintf(clustF[i],"\n name = \"%s\"\n",eQMbasis_names[qm->QMbasis]);
    fprintf(clustF[i]," molecule = $:molecule\n ) \n mpqc: (\n do_gradient = 1\n optimize = %d\n restart = %d\n checkpoint = %d\n savestate = %d\n  mole<",qm->bOPT,qm->bQMrestart,qm->bQMcheckpoint,qm->bQMsavestate);
    fprintf(clustF[i],"%s",eQMmethod_names[qm->QMmethod]);
    fprintf(clustF[i],">: (\n molecule = $:molecule \n basis = $:basis \n maxiter = 200 \n level_shift = %10.6f",qm->level_shift);
    if(qm->guessWF==10){
      fprintf(clustF[i],"\n gradient_accuracy=1e-5"); //WHICH IS EQUIVALENT TO AN ENERGY ACCURACY OF 1e-7
    }

    switch(qm->QMmethod) {
      case eQMmethodUKS:fprintf(clustF[i],"%s","\n functional<StdDenFunctional>: ( name = \"B3LYP\" )"); break;
      case eQMmethodCLKS:fprintf(clustF[i],"%s","\n functional<StdDenFunctional>: ( name = \"B3LYP\" )"); break;
      case eQMmethodHSOSKS:fprintf(clustF[i],"%s","\n functional<StdDenFunctional>: ( name = \"B3LYP\" )"); break;
      case eQMmethodMBPT2:fprintf(clustF[i],"%s","\n reference<CLHF>: (\n molecule = $:molecule\n basis = $:basis\n)"); break;
    }
    if(qm->guessWF>=1){
//      if((qm->guessWF+zf)>1){
      fprintf(clustF[i],"\n keep_guess_wavefunction = 1 \n");
      fprintf(clustF[i],"\n guess_wavefunction<CLHF>: (\n total_charge = %d\n",anint002(tchrg[i]));
      fprintf(clustF[i],"%s","molecule<Molecule>:(\n unit = bohr\n");
      fprintf(clustF[i],"%s","{ atoms geometry } = {");
      for(qi=0;qi<nrQM;qi++){
        qj=qm->QMcluster[qi*2];qk=qm->QMcluster[qi*2+1];
        if(qj==i){
          if(qm->atomicnumberQM[qk]==0){                 //LA
              fprintf(clustF[i],"\n H [% 6.3f % 6.3f % 6.3f]",qm->rQM[qk][XX]/BORH2NM,qm->rQM[qk][YY]/BORH2NM,qm->rQM[qk][ZZ]/BORH2NM);
          }
          else{
              fprintf(clustF[i],"\n%s [% 6.3f % 6.3f % 6.3f]",getAtomNameMPQC(qm->atomicnumberQM[qk]),qm->rQM[qk][XX]/BORH2NM,qm->rQM[qk][YY]/BORH2NM,qm->rQM[qk][ZZ]/BORH2NM);
          }
        }
      }
      fprintf(clustF[i],"\n }\n )\n");
      fprintf(clustF[i]," basis<GaussianBasisSet>: (\n");// BEGIN:: molecule = $:molecule
      fprintf(clustF[i],"%s","molecule<Molecule>:(\n unit = bohr\n");
      fprintf(clustF[i],"%s","{ atoms geometry } = {");
      for(qi=0;qi<nrQM;qi++){
        qj=qm->QMcluster[qi*2];qk=qm->QMcluster[qi*2+1];
        if(qj==i){
          if(qm->atomicnumberQM[qk]==0){                 //LA
              fprintf(clustF[i],"\n H [% 6.3f % 6.3f % 6.3f]",qm->rQM[qk][XX]/BORH2NM,qm->rQM[qk][YY]/BORH2NM,qm->rQM[qk][ZZ]/BORH2NM);
          }
          else{
              fprintf(clustF[i],"\n%s [% 6.3f % 6.3f % 6.3f]",getAtomNameMPQC(qm->atomicnumberQM[qk]),qm->rQM[qk][XX]/BORH2NM,qm->rQM[qk][YY]/BORH2NM,qm->rQM[qk][ZZ]/BORH2NM);
          }
        }
      }
      fprintf(clustF[i],"\n }\n )\n"); // END:: molecule = $:molecule 
      fprintf(clustF[i]," name = \"%s\"\n )\n )","6-31G");//eQMbasis_names[qm->QMbasis]);
/*      } else {
        fprintf(clustF[i], "\n     guess_wavefunction<CLHF>: (\n       total_charge = %d\n       molecule = $:molecule\n       basis<GaussianBasisSet>: (\n          molecule = $:molecule\n          name = \"6-31G\"\n       )\n     )",anint002(tchrg[i]));
      }
*/
    }
    fprintf(clustF[i], "\n print_npa = %d\n total_charge = %d \n)\n",qm->bPrint_npa,anint002(tchrg[i]));
    if(qm->bOPT){
      fprintf(clustF[i]," opt<QNewtonOpt>: (\n");
      fprintf(clustF[i]," function = $..:mole\n");
      fprintf(clustF[i]," update<BFGSUpdate>: ()\n");
      fprintf(clustF[i]," convergence<MolEnergyConvergence>: (\n");
      fprintf(clustF[i]," cartesian = yes\n");
      fprintf(clustF[i]," energy = $..:..:mole\n");
      fprintf(clustF[i],"\n)\n)\n");
    }
    fprintf(clustF[i],"%s",")\n"); 
  }
  if(nrMM >= 1 && (qm->includeQ>eIncludeNONE) ) {
    fclose(allMMndx);
  }
  for(i=qm->dI;i<NRC;i++){
    fclose(clustF[i]);
    fclose(indexF[i]);
    if(nrZ[i]==0){ 
	qm->nrQMc--;
    }
  }
}

real read_mpqc(rvec QMgrad[],rvec MMgrad[], t_QMrec *qm, t_MMrec *mm, int I)
{
  int  i,j,k,atnum,HL=0,nQMc,logi,nrHERE=0,nrMM,nrQM,nrCOL,leave=0,once=1,FAIL=0,npa; 
  char buf[MPQCLEN],atchar[2];
  real QMener=0.0,fei1,fei2,deltaE=0.0,tol=1;
  FILE *out,*mpqcerr,*mpqcMM,*mpqcQA;
  FILE *npalog;
  char txt1[20],txt2[20],txt3[20],txt4[20],txt5[20],filenmerr[50],filempqc[50],fileqa[50];
  char txt6[20],txt7[20],txt8[20];
  char filenm[100],gradtxt1[20],gradtxt2[20],enetxt1[20],enetxt2[20],enetxt3[20],enetxt4[20],hline[50],fline[MPQCLEN];
  char npatxt1[50],npatxt2[50],npatxt3[50],oline[MPQCLEN];
  char grdX[50],grdY[50],grdZ[50],eneQM[50];
  char dum1[50],dum2[50],dum3[50],dum4[50];
  rvec *Qgrad,*Mgrad,*QQgrad;
  
  sprintf(gradtxt1,"Total"); sprintf(gradtxt2,"Gradient:");
  sprintf(enetxt1,"Value"); sprintf(enetxt2,"of"); sprintf(enetxt3,"the"); sprintf(enetxt4,"MolecularEnergy:");
  sprintf(npatxt1,"Natural");sprintf(npatxt2,"Population");sprintf(npatxt3,"Analysis:");
  sprintf(filenm,"run_mpqc%d.out",I);
  
  out=fopen(filenm,"r");
  rewind(out);

  npa=(int)(qm->bPrint_npa);
  if(qm->bPrint_npa){
    npalog=fopen("npa.log","a");
  }  

  nQMc=qm->nrQMc;
  nrMM=mm->nrMMatoms;
  nrQM=qm->nrQMatoms;
  nrCOL=nrMM+nrQM;

  snew(Qgrad,nrQM+1);
  snew(QQgrad,nrQM+1);
  snew(Mgrad,nrMM+1);
  
  j=0;
  while(!feof(out)) {
    j++;
    if( NULL == fgets(fline,MPQCLEN,out)) {               /* Error handling */
      leave++;
      break;
    }
    else{
      sscanf(fline," %s %s %s %s %s %s %s %s\n",txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);

      if( !strcmp(txt1,"iter")) { /* Assigns energy from each interation */
	QMener=atof(txt5);
	deltaE=atof(txt8);
      }
      
      if( !strcmp(txt1,gradtxt1) && !strcmp(txt2,gradtxt2) && once) {
        nrHERE=0;once=0;
	if( nrMM>0 && qm->includeQ>eIncludeNONE){
          for(i=0;i<nrMM;i++){
	    logi=(int)mm->MMcluster[I*nrCOL+i];        /* BOOLEAN (IS THIS MM ATOM PART OF THIS CLUSTER) */
	    if(!((logi==0)||(logi==1))){
	      gmx_fatal(FARGS,"ERROR IN LOGIC (THIS SHOULD NEVER HAPPEN)");
            } 
	    if(logi){
	      HL++; nrHERE++;
	      if(NULL == fgets(fline,MPQCLEN,out))
		{
		  gmx_fatal(FARGS,"Error reading MM MPQC output (QM OUTPUT MIGHT NOT BE COMPLETE)");
		} 
              sscanf(fline,"%s %s %f %f %f\n",atchar,hline,&Mgrad[i][XX],&Mgrad[i][YY],&Mgrad[i][ZZ]);
	      /*sscanf(fline,"%s %s %s %s %s\n",atchar,hline,grdX,grdY,grdZ);
		Mgrad[i][XX]=atof(grdX);
		Mgrad[i][YY]=atof(grdY);
		Mgrad[i][ZZ]=atof(grdZ);*/
	    }
	  }
	  //REMOVE ALL EXTRA Q ATOMS HERE
	  nrHERE=0;
	  for(i=0;i<nrQM;i++){                
	    j=qm->QMcluster[2*i+1];
	    if(mm->MMcluster[(j+nrMM)+I*nrCOL]){ //IS THE QM ATOM A POINTCHARGE FOR THIS CLUSTER?
	      if(NULL == fgets(fline,MPQCLEN,out))
		{
		  gmx_fatal(FARGS,"Error reading QA MPQC output");
		}
	      HL++;
              sscanf(fline,"%s %s %f %f %f\n",atchar,hline,&QQgrad[j][XX],&QQgrad[j][YY],&QQgrad[j][ZZ]);
	      //sscanf(fline,"%s %s %s %s %s\n",atchar,hline,grdX,grdY,grdZ);
              //adds MM contribution from non-self clusters 
	      /*QQgrad[j][XX]=atof(grdX);
		QQgrad[j][YY]=atof(grdY);
		QQgrad[j][ZZ]=atof(grdZ);*/

	    }
	  } 
	  //END EXCESS QM Q
	}
	for(i=0;i<nrQM;i++){
	  if(qm->QMcluster[2*i]==I){
	    j=qm->QMcluster[2*i+1];
	    if(NULL == fgets(fline,MPQCLEN,out)){
	      gmx_fatal(FARGS,"Error reading QM MPQC output");
	    }
            sscanf(fline,"%s %s %f %f %f\n",atchar,hline,&Qgrad[j][XX],&Qgrad[j][YY],&Qgrad[j][ZZ]);
	    //sscanf(fline,"%s %s %s %s %s\n",atchar,hline,grdX,grdY,grdZ);
	    if(hline[1]=='\0'){
	      hline[1]=hline[0];
	      hline[0]=' ';
	      hline[2]='\0';
	    }
	    //CHECK THAT ATOMTYPE IS CORRECT
	    if( !(( !strcmp(" H",hline) & !strcmp("LA",getAtomNameMPQC(qm->atomicnumberQM[j])) ) | !strcmp(hline,getAtomNameMPQC(qm->atomicnumberQM[j])) )){ 
	      fprintf(stderr,"ATOM::%s ENTRY::%s",getAtomNameMPQC(qm->atomicnumberQM[j]),hline);
	      gmx_fatal(FARGS,"MISSMATCH OF QM ATOMS \nPOSSIBLE REASONS FOR ERROR:\n    1) LAUNCH environment variable not set\n    2) Non convergent run\n Else) True bug with code. (If mpqc printed timing information\n       then this is the case) ");
	    }
	    HL++;
	    /*Qgrad[j][XX]=atof(grdX);
	      Qgrad[j][YY]=atof(grdY);
	      Qgrad[j][ZZ]=atof(grdZ);*/
	  }
        }
      }
      
      if( (int)strspn(fline,"    Natural Population Analysis:")==32 && npa ){ // WE WILL ONLY GET NPA FROM CONVERGENT RUNS
        npa=0;
        if(NULL == fgets(fline,MPQCLEN,out)) {
	  gmx_fatal(FARGS,"Error reading MPQC NPA::0 output");
	}
        for(i=0;i<nrMM;i++){//MM
          if(mm->MMcluster[I*nrCOL+i]){
            if(NULL == fgets(fline,MPQCLEN,out)){
	      gmx_fatal(FARGS,"Error reading MPQC NPA output");
	    }
          }
        }
	for(i=0;i<nrQM;i++){//Q
          j=qm->QMcluster[2*i+1];
          if(mm->MMcluster[(j+nrMM)+I*nrCOL]){
            if(NULL == fgets(fline,MPQCLEN,out)){
	      gmx_fatal(FARGS,"Error reading MPQC NPA output");
	    }            
          }
        }
	for(i=0;i<nrQM;i++){//QM
	  if(qm->QMcluster[2*i]==I){
	    j=qm->QMcluster[2*i+1];
	    if(NULL == fgets(fline,MPQCLEN,out)){
	      gmx_fatal(FARGS,"Error reading MPQC QM NPA output");
	    }
            sprintf(oline,"%6d % 6.3f % 6.3f % 6.3f %s",I,qm->rQM[j][XX],qm->rQM[j][YY],qm->rQM[j][ZZ],fline); 
            fprintf(npalog,"%s",oline);
          }
        }  
      }
      
      if( strspn(fline,"  Value of the MolecularEnergy:") > 30  && leave==0) { //!strcmp(txt1,enetxt1) && !strcmp(txt2,enetxt2) && !strcmp(txt3,enetxt3) && !strcmp(txt4,enetxt4) ) { 
	sscanf(fline,"%s %s %s %s %f",dum1,dum2,dum3,dum4,&QMener);
        leave--;
      }
    }
  }
  
  if(leave){
    fprintf(stderr ,"\n!!QM DIDNT PROPERLY CONVERGE!!\nEXTRACTED ENERGY AND GRADIENT::\nENERGY=%10.6f TOL=%10.6e REM=%d\n",QMener,deltaE,HL);
    for(i=0;i<nrCOL;i++)
      if((mm->MMcluster[I*nrCOL+i]) &&  qm->includeQ>eIncludeNONE)
	HL--;
    for(i=0;i<qm->nrQMatoms;i++)
      if(qm->QMcluster[2*i]==I)
	HL--;
    if(deltaE>=1e-7){
      FAIL=1;
      qm->dI=I+1;
    }
  }

  if(!FAIL){
    for(i=0;i<nrQM;i++){
      j=qm->QMcluster[2*i+1];
      if(qm->QMcluster[2*i]==I || mm->MMcluster[(j+nrMM)+I*nrCOL]){
        for(k=XX;k<=ZZ;k++)
          QMgrad[j][k]+=Qgrad[j][k];
      }
      if(mm->MMcluster[(j+nrMM)+I*nrCOL]){
        for(k=XX;k<=ZZ;k++)
          QMgrad[j][k]+=QQgrad[j][k];
      }
    }
    if( nrMM>0 && qm->includeQ>eIncludeNONE){
      for(i=0;i<nrMM;i++){
        if((int)mm->MMcluster[I*nrCOL+i]){
          for(k=XX;k<=ZZ;k++)
            MMgrad[i][k]+=Mgrad[i][k];  
        }
      }
    } 
  }

  free(Qgrad);
  free(QQgrad);
  free(Mgrad);
  fclose(out);
  if(qm->bPrint_npa){
    fclose(npalog);
  }

  return(QMener); 
}

real call_mpqc(t_commrec *cr,  t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[])
{
  /* normal mpqc jobs */
  int i,I,j,RV=0,cl=1,NFILES=1,NFILER=1,GWF;
  real QMener=0.0,QME0=0.0;
  rvec *QMgrad,*MMgrad;
  char *exe;
  real *qmclusters,*mmclusters;
  FILE *stdver;

  /* Allocate */
  snew(QMgrad,qm->nrQMatoms);
  snew(MMgrad,mm->nrMMatoms);
  //fprintf(stderr,"nrMMatoms %d\n",mm->nrMMatoms);

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
  qm->dI=0;

  if(qm->bCLUSTER){
    clustering(qm);  
    if(mm->nrMMatoms>0){
      clusterMM(mm,qm);
    }
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

  for(i=0;i<NFILES;i++){
    RV=do_mpqc(qm,i);
    fflush(stdout); fflush(stdin); 
    qm->dI=0;
    QME0 = read_mpqc(QMgrad,MMgrad,qm,mm,i);
    if(qm->dI==0){            //SUCCESS
      QMener+=QME0;  
    }
    else{                
      fprintf(stderr,"FAILED ONCE: REDO\n");
      GWF=qm->guessWF;
      qm->guessWF=10;
      write_mpqc(fr,qm,mm);
      RV=do_mpqc(qm,i);
      fflush(stdout); fflush(stdin); 
      qm->dI=0;
      QME0 = read_mpqc(QMgrad,MMgrad,qm,mm,i);
      qm->guessWF=GWF;
      if(!(qm->dI==0)){          //FAILED AGAIN
        gmx_fatal(FARGS,"MPQC FAILED TWICE");
      }else{                     //SUCCESS
        QMener+=QME0;  
      }
    }
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
  char  input_string[MPQCLEN],exename[200],execute[10],eo0[20],eo1[20],eo2[20],eo3[20],eo4[20],tfi[100];
  char  options[1000],launch[1000],rmout[200];

  //sprintf(eo0,"-o=run_mpqc%d.out",I);
  if(qm->nQMnodes>0){
    sprintf(eo1,"--nnodeperjob=%d",qm->nQMnodes); 
  }
  else{
    sprintf(eo1," ");
  }
  if(qm->nQMcpus>0){
    sprintf(eo2,"--nthreadperproc=%d",qm->nQMcpus);
  }
  else{
    sprintf(eo2," ");
  }
  if(!(strstr(qm->gauss_exe,"mpqcrun")==NULL)){
    sprintf(eo3,"--autoout");
    sprintf(eo4,"--simpout");
  }
  else{
    sprintf(eo3," "); 
    sprintf(eo4,"-o=run_mpqc%d.out",I);
  }
  sprintf(options,"$OPTIONS");
  sprintf(launch,"--launch=\"$LAUNCH\"");

  sprintf(exename,"run_mpqc%d.inp",I);
  sprintf(input_string,"%s %s %s %s %s %s %s %s",qm->gauss_exe,options,eo1,eo2,eo3,eo4,launch,exename);

#ifdef GMX_NO_SYSTEM
  printf("Warning-- No calls to system(3) supported on this platform.");
  gmx_fatal(FARGS,"Call to '%s' failed\n",input_string);
#else
  sprintf(tfi,"touch run_mpqc%d.out",I);  
  sprintf(rmout,"rm run_mpqc%d.out",I);
  system(tfi);
  system(rmout); //REMOVES OLD OUTPUT FILE
  if ( system(input_string) != 0 )
    gmx_fatal(FARGS,"Call to '%s' failed\n",input_string);
#endif
  
  return(0);
}

