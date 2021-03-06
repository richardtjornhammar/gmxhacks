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
#define FEL     138.935485
#define nan2a0  18.89725989

char* getAtomNameMPQC(int i)
{
  char* atomName[120];

  atomName[0] = "LA";
  atomName[1] = " H"; 
  atomName[2] = "He";
  atomName[3] = "Li";
  atomName[4] = "Be";
  atomName[5] = " B";
  atomName[6] = " C";
  atomName[7] = " N";
  atomName[8] = " O";
  atomName[9] = " F";
  atomName[10] ="Ne";
  atomName[11] ="Na";
  atomName[12] ="Mg";
  atomName[13] ="Al";
  atomName[14] ="Si";
  atomName[15] =" P";
  atomName[16] =" S";
  atomName[17] ="Cl";
  atomName[18] ="Ar";
  atomName[19] =" K";
  atomName[20] ="Ca";
  atomName[21] ="Sc";
  atomName[22] ="Ti";
  atomName[23] =" V";
  atomName[24] ="Cr"; 
  atomName[25] ="Mn";
  atomName[26] ="Fe";
  atomName[27] ="Co";
  atomName[28] ="Ni";
  atomName[29] ="Cu";
  atomName[30] ="Zn";
  atomName[31] ="Ga";
  atomName[32] ="Ge";
  atomName[33] ="As";
  atomName[34] ="Se";
  atomName[35] ="Br";
  atomName[36] ="Kr";
  atomName[37] ="Rb";
  atomName[38] ="Sr";
  atomName[39] =" Y";
  atomName[40] ="Zr";
  atomName[41] ="Nb";
  atomName[42] ="Mo";
  atomName[43] ="Tc";
  atomName[44] ="Ru";
  atomName[45] ="Rh";
  atomName[46] ="Pd";
  atomName[47] ="Ag";
  atomName[48] ="Cd";
  atomName[49] ="In";
  atomName[50] ="Sn";
  atomName[51] ="Sb";
  atomName[52] ="Te";
  atomName[53] =" I";
  atomName[54] ="Xe";
  atomName[55] ="Cs";
  atomName[56] ="Ba";
  atomName[57] ="La";
  atomName[58] ="Ce";
  atomName[59] ="Pr";
  atomName[60] ="Nd"; 
  atomName[61] ="Pm";
  atomName[62] ="Sm";
  atomName[63] ="Eu";
  atomName[64] ="Gd";
  atomName[65] ="Tb";
  atomName[66] ="Dy";
  atomName[67] ="Ho";
  atomName[68] ="Er";
  atomName[69] ="Tm";
  atomName[70] ="Yb";
  atomName[71] ="Lu";
  atomName[72] ="Hf";
  atomName[73] ="Ta";
  atomName[74] =" W";
  atomName[75] ="Re";
  atomName[76] ="Os";
  atomName[77] ="Ir";
  atomName[78] ="Pt";
  atomName[79] ="Au";
  atomName[80] ="Hg";
  atomName[81] ="Tl";
  atomName[82] ="Pb";
  atomName[83] ="Bi";
  atomName[84] ="Po";
  atomName[85] ="At";
  atomName[86] ="Rn";
  atomName[87] ="Fr";
  atomName[88] ="Ra";
  atomName[89] ="Ac";
  atomName[90] ="Th";
  atomName[91] ="Pa";
  atomName[92] =" U";
  atomName[93] ="Np";
  atomName[94] ="Pu";
  atomName[95] ="Am";
  atomName[96] ="Cm";
  atomName[97] ="Bk";
  atomName[98] ="Cf";
  atomName[99] ="Es";
  atomName[100]="Fm";
  atomName[101]="Md";
  atomName[102]="No";
  atomName[103]="Lr";
  atomName[104]="Rf";
  atomName[105]="Db";
  atomName[106]="Sg";
  atomName[107]="Bh";
  atomName[108]="Hs";
  atomName[109]="Mt";
  atomName[110]="Ds";
  atomName[111]="Rg";

  return atomName[i];
}
/*
void get_append(t_QMrec *qm){
  FILE *fin;
  char  str[100],fn[100],*cptr;
  int   i,j,k;

  system("ls run_mpqc0.out* > blaj");
  fin=fopen("blaj","r");
  fgets(fn,100,fin);
  fclose(fin);

  for(i=0;i<100;i++)
    if(fn[i]=='o' && fn[i+1]=='u' && fn[i+2]=='t')
      break;
  i+=3;
  cptr=&fn[i];

  k=0;
  for(j=i;j<100;j++){
    if(k>=3 && !( (fn[j])=='0' || (fn[j])=='1' || (fn[j])=='2' || (fn[j])=='3'
 || (fn[j])=='4' || (fn[j])=='5' || (fn[j])=='6' || (fn[j])=='7' || (fn[j])=='8'
 || (fn[j])=='9' ))
      break;
    k+=fn[j]=='.'?1:0;
  }
  fn[j]='\0';
  strcpy(qm->mpqcappend,cptr);
}
*/
void init_mpqc(t_commrec *cr, t_QMrec *qm, t_MMrec *mm)
{
  char *buf,*p;
  int i,n;
  FILE *npalog;

  snew(buf,MPQCLEN);

  buf = getenv("MPQC_NNODEPERJOB");
  if (buf){
    sscanf(buf,"%d",&qm->nQMnodes);
  }
  else{
    qm->nQMnodes=1;
  }

  buf = getenv("MPQC_NPROCPERJOB");
  if (buf){
    sscanf(buf,"%d",&qm->nQMcpus);
  }
  else{
    qm->nQMcpus=1;
  }


  buf = getenv("MPQC_NPROCPERNODE");
  if (buf){
    sscanf(buf,"%d",&qm->nQMproc);
  }
  else{
    qm->nQMproc=1;
  }

  buf = getenv("MPQC_NTHREADPERPROC");
  if (buf){
    sscanf(buf,"%d",&qm->nQMthread);
  }
  else{
    qm->nQMthread=1;
  }

  i=0;
  buf=getenv("MPQC_EXE");
  if (buf){
    snew(qm->gauss_exe,MPQCLEN);
    sscanf(buf,"%s",qm->gauss_exe);
  }
  else{
    gmx_fatal(FARGS,"MUST >> export MPQC_EXE=[An executable with fully qualified path]\n");
  }

  if(qm->bPrint_npa){
    npalog=fopen("npa.log","w");
    fclose(npalog);
  }
  npalog=fopen("qmmm_ene.dat","w");
  fclose(npalog);

  fprintf(stderr,"nr of background atoms = %d\nMPQC initialised...\n",mm->nrMMatoms);
}

void write_mpqc(t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, int FAIL)
{
  // NOTE:: MPQC HAS OWN INPUT CHECKING FUNCTIONS
  char  wstring[MAXEL],dstring[MAXEL],tstring[MAXEL],estring[MAXEL];
  int   i,i0,j,qcflag=0,NRC=0,k=0,l=0,T,q=0,qq=0,zf=1,qi,qj,qk;
  int   nrZ[MAXFILES],nrQ[MAXFILES],nrCOL=0,nrMM=0,nrQM=0,a,b;
  double  tchrg[MAXFILES];
  FILE  *clustF[MAXFILES]; /* MAXFILES is set in cluster.h */ 
  FILE  *indexF[MAXFILES];  
  FILE  *allMMndx;

  strcpy(tstring,"run_mpqc");
  strcpy(dstring,"run_mpqc");
  strcpy(estring,".inp");
  NRC=qm->nrQMc;   nrQM=qm->nrQMatoms;
  nrMM=mm->nrMMatoms; nrCOL=nrMM+nrQM;

  if( NRC >= MAXFILES || qm->QMexec[0] >= MAXFILES)
  {
    gmx_fatal(FARGS,"To many clusters. Your QM system is probably exploding");
  }

  for(i=0;i<qm->QMexec[0];i++){
    tchrg[i]=0;
    sprintf(dstring,"run_mpqc%d",qm->QMexec[i+1]);
    clustF[i]=fopen(strcat(dstring,estring),"w");
    for(j=0;j<MAXEL;j++){
       if(dstring[j]=='.')
         break;
    }
    dstring[j+1]='n';dstring[j+2]='d';dstring[j+3]='x';
    indexF[i]=fopen(dstring,"w");
    fprintf(clustF[i],"%s","molecule<Molecule>:(\n  unit = bohr");
    nrZ[i]=0; nrQ[i]=0;
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
 ======>includeYES:: STABLE
 This is the default QMMM option within MPQC
 ======>includeQ::   NOT-STABLE
 If this option is true, then the Q atoms are included when computing the nuclear charge and the electric field due to the nuclear charge.
 ======>includeQQ::  NOT-STABLE
 If this option is true, then the Q atoms are included when computing the nuclear repulsion energy and its derivatives. Also all interaction energies between the different Q atoms are included.
 ======>includeQQQ:: NOT-STABLE
 If this option is true, then the Q atoms are included when computing the nuclear repulsion energy and its derivatives and are also are included when computing the nuclear charge and the electric field due to the nuclear charge.
*/
  if((qm->includeQ>eIncludeNONE)) {
    if(qm->includeQ==eIncludeQ||qm->includeQ==eIncludeQQQ)
      q=1;
    if(qm->includeQ==eIncludeQQ||qm->includeQ==eIncludeQQQ)
      qq=1;
    for(i=0;i<qm->QMexec[0];i++){
      if(!(qm->includeQ==eIncludeYES)){ //YES IS THE MPQC DEFAULT FOR QMMM
        fprintf(clustF[i],"\n  include_q=%d\n  include_qq=%d",q,qq);
      }
      fprintf(clustF[i],"\n  charge = [");
      fprintf(indexF[i],"[ MMrec_c%d ]\n",i);
    }
    if(nrMM>0){
      for(i=0;i<nrMM;i++){
        for(j=0;j<qm->QMexec[0];j++){
          if(mm->MMcluster[nrCOL*j + i]){
	    fprintf(clustF[j],"% 8.5f ",(mm->MMcharges[i])*(mm->fcharges[nrCOL*j + i]) );
            fprintf(indexF[j],"% 10d", mm->indexMM[i]+1 );
            nrQ[j]++;
            if(nrQ[j]==10){
              fprintf(indexF[j],"\n");
              nrQ[j]=0;
            }
          }
        }
      }
    }
    for(i=0;i<nrQM;i++){
//      a=qm->QMcluster[i*2];      // CLUSTER INDEX a:(0 -> NC-1)  
//      b=qm->QMcluster[i*2+1];    // QM INDEX      b:(0 -> nrQM-1)
      for(j=0;j<qm->QMexec[0];j++){
	if(mm->MMcluster[nrCOL*j+i+nrMM]){
	  fprintf(clustF[j],"% 8.5f ",(qm->pcharges[i])*(mm->fcharges[nrCOL*j+i+nrMM]));
          fprintf(indexF[j],"% 10d",qm->indexQM[i]+1);  //MDRUN INDEX AND NDX INDEX DIFFER WITH 1
          nrQ[j]++;
          if(nrQ[j]==10){
            fprintf(indexF[j],"\n");
            nrQ[j]=0;
          }
        }
      }      
    }
    for(i=0;i<qm->QMexec[0];i++){
      fprintf(clustF[i],"%s","]");
    }
    qcflag=1;
  }

  for(i=0;i<qm->QMexec[0];i++){
    nrQ[i]=0;  
    fprintf(indexF[i],"\n\n[ QMrec_c%d ]\n",i);
    fprintf(clustF[i],"%s","\n  { atoms ghost geometry } = {");
  }
  if(qcflag) {
    if(nrMM>0){
      for(i=0;i<nrMM;i++){
        for(j=0;j<qm->QMexec[0];j++){
	  if(mm->MMcluster[ nrCOL*j + i ]){
	    fprintf(clustF[j],"\n Q false [% 6.3f % 6.3f % 6.3f]",mm->rMM[ nrCOL*j + i ][XX]/BORH2NM,mm->rMM[ nrCOL*j + i ][YY]/BORH2NM,mm->rMM[ nrCOL*j + i ][ZZ]/BORH2NM);
	  }
        }
      }
    }
    for(i=0;i<nrQM;i++){
//      a=qm->QMcluster[i*2];      // CLUSTER INDEX a:(0 -> NC-1)    */
//      b=qm->QMcluster[i*2+1];    //      QM INDEX b:(0 -> nrQM-1)  */
    if(1){          // we always use Q atoms qm->useQMatom[i] needs corresponding at the charge declaration
      for(j=0;j<qm->QMexec[0];j++){
	if(mm->MMcluster[i+nrMM+nrCOL*j]){
	  fprintf(clustF[j],"\n Q false [% 6.3f % 6.3f % 6.3f]",mm->rMM[ i+nrMM+nrCOL*j ][XX]/BORH2NM,mm->rMM[ i+nrMM+nrCOL*j ][YY]/BORH2NM,mm->rMM[ i+nrMM+nrCOL*j ][ZZ]/BORH2NM);
	}
      }
    }
    } 
  }
  for(i=0;i<nrQM;i++){
    j=qm->QMcluster[i*2]; k=qm->QMcluster[i*2+1];
    if( !(j >= 0 && j < qm->QMexec[0])){
      break;
    }
    tchrg[j]+=qm->pcharges[k];                    //Charge contributions are summed up from the topology definition
    //if(qm->useQMatom[k]){
    if(qm->atomicnumberQM[k]==0){                 //LA a more clever way of doing this will be implemented
        fprintf(clustF[j],"\n H true  [% 6.3f % 6.3f % 6.3f]",qm->rQM[j*nrQM+k][XX]/BORH2NM,qm->rQM[j*nrQM+k][YY]/BORH2NM,qm->rQM[j*nrQM+k][ZZ]/BORH2NM);
    }
    else{
        fprintf(clustF[j],"\n%s false [% 6.3f % 6.3f % 6.3f]",getAtomNameMPQC(qm->atomicnumberQM[k]),qm->rQM[j*nrQM+k][XX]/BORH2NM,qm->rQM[j*nrQM+k][YY]/BORH2NM,qm->rQM[j*nrQM+k][ZZ]/BORH2NM);
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
    //}useqmatom
  }
  for(i=0;i<qm->QMexec[0];i++){
    fprintf(indexF[i],"\n");
    fprintf(clustF[i],"%s","\n  }\n)\nbasis<GaussianBasisSet>: (");
    fprintf(clustF[i],"\n  name = \"%s\"\n",eQMbasis_names[qm->QMbasis]);
    fprintf(clustF[i],"  molecule = $:molecule\n) \nmpqc: (\n  do_gradient = 1\n  optimize = %d\n  restart = %d\n  checkpoint = %d\n  savestate = %d\n  mole<",qm->bOPT,qm->bQMrestart,qm->bQMcheckpoint,qm->bQMsavestate);
    fprintf(clustF[i],"%s",eQMmethod_names[qm->QMmethod]);
    fprintf(clustF[i],">: (\n    molecule = $:molecule \n    basis = $:basis \n    maxiter = 300 \n    level_shift = %10.6f",(FAIL==0)?(qm->level_shift):(qm->level_shift+1.0));
    switch(qm->QMmethod) {
      case eQMmethodUKS:fprintf(clustF[i],"\n    functional<StdDenFunctional>: ( name = \"%s\" )",eQMfunctional_names[qm->QMfunctional]); break;
      case eQMmethodCLKS:fprintf(clustF[i],"\n    functional<StdDenFunctional>: ( name = \"%s\" )",eQMfunctional_names[qm->QMfunctional]); break;
      case eQMmethodHSOSKS:fprintf(clustF[i],"\n    functional<StdDenFunctional>: ( name = \"%s\" )",eQMfunctional_names[qm->QMfunctional]); break;
      case eQMmethodMBPT2:fprintf(clustF[i],"\n    reference<CLHF>: (\n      molecule = $:molecule\n      basis = $:basis\n)"); break;
    }
    if(qm->guessWF>=eGWFyes){
      fprintf(clustF[i],"\n    keep_guess_wavefunction = 1");
      fprintf(clustF[i],"\n    guess_wavefunction<CLHF>: (\n      total_charge = %d\n", (int)round(tchrg[i]) );// anint002(tchrg[i]));
      fprintf(clustF[i],"%s","      molecule = $:molecule\n");
      fprintf(clustF[i],"      basis<GaussianBasisSet>: (\n");
      fprintf(clustF[i],"%s","        molecule = $:molecule\n");
      fprintf(clustF[i],"        name = \"%s\"\n      )\n    )",(qm->guessWF==eGWFb)?("6-31G"):(  (qm->guessWF==eGWFs)?("3-21G"):("STO-3G") ) );
    }
    fprintf(clustF[i], "\n    print_npa = %d\n    total_charge = %d \n  )\n",qm->bPrint_npa,(int)round(tchrg[i]) );//anint002(tchrg[i]));
    if(qm->bOPT){
      fprintf(clustF[i],"  opt<QNewtonOpt>: (\n");
      fprintf(clustF[i],"    function = $..:mole\n");
      fprintf(clustF[i],"    update<BFGSUpdate>: ()\n");
      fprintf(clustF[i],"    convergence<MolEnergyConvergence>: (\n");
      fprintf(clustF[i],"    cartesian = yes\n");
      fprintf(clustF[i],"    energy = $..:..:mole\n");
      fprintf(clustF[i],"\n    )\n  )\n");
    }
    fprintf(clustF[i],"%s",")\n"); 
  }
  if(nrMM >= 1 && (qm->includeQ>eIncludeNONE) ) {
    fclose(allMMndx);
  }
  for(i=0;i<qm->QMexec[0];i++){
    fclose(clustF[i]);
    fclose(indexF[i]);
    if(nrZ[i]==0){ 
	qm->nrQMc--;
    }
  }
}

real read_mpqc(rvec QMgrad[],rvec MMgrad[], t_QMrec *qm, t_MMrec *mm, int I)
{
  int  i,j,k,atnum,HL=0,nQMc,logi,nrHERE=0,nrMM,nrQM,nrCOL,leave=0,once=1,FAIL=0,npa,iii=0; 
// BEGIN FOR NUC. REP.
  int  inrep=0;
  real Erep;
// END FOR NUC. REP.
  char buf[MAXEL],atchar[2];
  real QMener=0.0,fei1,fei2,deltaE=0.0,tol=1;
  FILE *out,*mpqcerr,*mpqcMM,*mpqcQA;
  FILE *npalog;
  char txt1[20],txt2[20],txt3[20],txt4[20],txt5[20],filenmerr[50],filempqc[50],fileqa[50];
  char txt6[20],txt7[20],txt8[20],str[5][100];
  char filenm[100],gradtxt1[20],gradtxt2[20],enetxt1[20],enetxt2[20],enetxt3[20],enetxt4[20],hline[50],fline[MAXEL];
  char npatxt1[50],npatxt2[50],npatxt3[50],oline[MAXEL];
  char grdX[50],grdY[50],grdZ[50],eneQM[50];
  char dum1[50],dum2[50],dum3[50],dum4[50];
  rvec *Qgrad,*Mgrad,*QQgrad;
  rvec  Ngrad,Fgrad;
  double DN;
  
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
  qm->Nrep[0]=0; qm->Nrep[1]=0;

  nQMc=qm->nrQMc;
  nrMM=mm->nrMMatoms;
  nrQM=qm->nrQMatoms;
  nrCOL=nrMM+nrQM;

  snew(Qgrad,nrQM+1);
  snew(QQgrad,nrQM+1);
  snew(Mgrad,nrMM+1);
  Ngrad[XX]=0.0; Ngrad[YY]=0.0; Ngrad[ZZ]=0.0;
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
	if( qm->includeQ>eIncludeNONE ){
          if(nrMM>0){
            for(i=0;i<nrMM;i++){
              logi=(int)mm->MMcluster[I*nrCOL+i];        // BOOLEAN (IS THIS MM ATOM PART OF THIS CLUSTER) 
              Mgrad[i][XX]=0.0; Mgrad[i][YY]=0.0; Mgrad[i][ZZ]=0.0;
              if(logi){
                qm->Nrep[0]++;
	        HL++; nrHERE++;
                if(NULL == fgets(fline,MPQCLEN,out))
		{
                  fprintf(stderr,"\nCURRENT FILENR::%d/%d NRQMATOMS::%d NRMMATOMS::%d ",I,nQMc,nrQM,nrMM);
		  gmx_fatal(FARGS,"\nError reading MM MPQC output (QM OUTPUT MIGHT NOT BE COMPLETE)");
		} 
//                sscanf(fline,"%s %s %f %f %f\n",atchar,hline,&Mgrad[i][XX],&Mgrad[i][YY],&Mgrad[i][ZZ]);
                sscanf(fline,"%s %s %s %s %s\n",atchar,hline,str[XX],str[YY],str[ZZ]);
	        if(hline[1]=='\0'){
	          hline[1]=hline[0];
	          hline[0]=' ';
	          hline[2]='\0';
	        }
                if(!(!strcmp(" Q",hline))){ //Check the type
                    fprintf(stderr,"\nGOT::%s\n",hline);
                    fprintf(stderr,"\nCURRENT FILENR::%d/%d NRQMATOMS::%d NRMMATOMS::%d ",I,nQMc,nrQM,nrMM);
                    fprintf(stderr,"\nError reading MM MPQC output");
                    FAIL=1;
                    break;
                }
                Mgrad[i][XX]=atof(str[XX]);
                Mgrad[i][YY]=atof(str[YY]);
                Mgrad[i][ZZ]=atof(str[ZZ]);  
              }
            }
	  }
	  //REMOVE ALL EXTRA Q ATOMS HERE
	  nrHERE=0;
	  for(i=0;i<nrQM;i++){                
	    // j=qm->QMcluster[2*i+1];
            QQgrad[i][XX]=0.0; QQgrad[i][YY]=0.0; QQgrad[i][ZZ]=0.0;
	    if(mm->MMcluster[(i+nrMM)+I*nrCOL]){     //IS THE QM ATOM A POINTCHARGE FOR THIS CLUSTER?
              qm->Nrep[0]++;
	      if(NULL == fgets(fline,MPQCLEN,out))
		{
                  fprintf(stderr,"\nCURRENT FILENR::%d/%d NRQMATOMS::%d NRMMATOMS::%d ",I,nQMc,nrQM,nrMM);
		  gmx_fatal(FARGS,"\nError reading QA MPQC output (QM OUTPUT MIGHT NOT BE COMPLETE)");
		}
	      HL++;
//              sscanf(fline,"%s %s %f %f %f\n",atchar,hline,&QQgrad[i][XX],&QQgrad[i][YY],&QQgrad[i][ZZ]);
              sscanf(fline,"%s %s %s %s %s\n",atchar,hline,str[XX],str[YY],str[ZZ]);
              QQgrad[i][XX]=atof(str[XX])*qm->fMMgradMask;
              QQgrad[i][YY]=atof(str[YY])*qm->fMMgradMask;
              QQgrad[i][ZZ]=atof(str[ZZ])*qm->fMMgradMask;
              Ngrad[XX]+=QQgrad[i][ZZ];
              Ngrad[YY]+=QQgrad[i][YY];
              Ngrad[ZZ]+=QQgrad[i][ZZ];
	      if(hline[1]=='\0'){
	        hline[1]=hline[0];
	        hline[0]=' ';
	        hline[2]='\0';
	      }
              if(!(!strcmp(" Q",hline))){ //Check the type
                fprintf(stderr,"\nGOT::%s\n",hline);
                fprintf(stderr,"\nCURRENT FILENR::%d/%d NRQMATOMS::%d NRMMATOMS::%d ",I,nQMc,nrQM,nrMM);
                fprintf(stderr,"\nError reading Q MPQC output");
                FAIL=1;
                break;
              }
	    }
	  } 
	  //END EXCESS QM Q
	}
	for(i=0;i<nrQM;i++){
          j=qm->QMcluster[2*i+1];
          Qgrad[j][XX]=0.0; Qgrad[j][YY]=0.0; Qgrad[j][ZZ]=0.0;
	  if(qm->QMcluster[2*i]==I){
          //  if(qm->useQMatom[j]){
              qm->Nrep[0]++;
              qm->Nrep[1]++;
	      if(NULL == fgets(fline,MPQCLEN,out)){
	        fprintf(stderr,"\nCURRENT FILENR::%d/%d NRQMATOMS::%d NRMMATOMS::%d ",I,nQMc,nrQM,nrMM);
	        fprintf(stderr,"\nError reading QM MPQC output");
                FAIL=1;
                break;
	      }
//            sscanf(fline,"%s %s %f %f %f\n",atchar,hline,&Qgrad[j][XX],&Qgrad[j][YY],&Qgrad[j][ZZ]);
              sscanf(fline,"%s %s %s %s %s\n",atchar,hline,str[XX],str[YY],str[ZZ]);
	      if(hline[1]=='\0'){
	        hline[1]=hline[0];
	        hline[0]=' ';
	        hline[2]='\0';
	      }
	      //CHECK THAT ATOMTYPE IS CORRECT
	      if( !(( !strcmp(" H",hline) & !strcmp("LA",getAtomNameMPQC(qm->atomicnumberQM[j])) ) | !strcmp(hline,getAtomNameMPQC(qm->atomicnumberQM[j])) )){ 
	        fprintf(stderr,"\nDUMPING DATA::\nATOM::%s ENTRY::%s",getAtomNameMPQC(qm->atomicnumberQM[j]),hline);
	        fprintf(stderr,"\nCURRENT FILENR::%d/%d NRQMATOMS::%d NRMMATOMS::%d ",I,nQMc,nrQM,nrMM);
                FAIL=1;
                break;
                //gmx_fatal(FARGS,"MISSMATCH OF QM ATOMS \nPOSSIBLE REASONS FOR ERROR:\n    1) LAUNCH environment variable not set\n    2) Non convergent run\n Else) True bug with code. (If mpqc printed timing information\n       then this is the case) ");
	      }
              Qgrad[j][XX]=atof(str[XX])*qm->fQMgradMask;
              Qgrad[j][YY]=atof(str[YY])*qm->fQMgradMask;
              Qgrad[j][ZZ]=atof(str[ZZ])*qm->fQMgradMask;
	      HL++;
          //  }useqmatom
	  }
        }
      }

      if( !(strstr(fline,"Natural Population Analysis")==NULL) && npa ) {
        npa=0;
        if(NULL == fgets(fline,MPQCLEN,out)) {
	  gmx_fatal(FARGS,"Error reading MPQC NPA::0 output");
	}
        dum2[0]='Q';
        while(!feof(out) && dum2[0]=='Q'){
          if(NULL == fgets(fline,MPQCLEN,out)) {
	    gmx_fatal(FARGS,"Error reading MPQC NPA::Q output");
	  }  
          sscanf(fline," %s %s %f",dum1,dum2,&fei1);
        }
	for(i=0;i<nrQM;i++){ // QM
	  if(qm->QMcluster[2*i]==I){
	    j=qm->QMcluster[2*i+1];
          //  if(qm->useQMatom[j]){
              if(!(npa==0)){
                if(NULL == fgets(fline,MPQCLEN,out)){
	          gmx_fatal(FARGS,"Error reading MPQC QM NPA output");
	        }
              }
              npa=1;
              sprintf(oline,"%6d % 6.3f % 6.3f % 6.3f %6d %s",I,qm->rQM[I*qm->nrQMatoms+j][XX],qm->rQM[I*qm->nrQMatoms+j][YY],qm->rQM[I*qm->nrQMatoms+j][ZZ],qm->indexQM[j],fline); 
              fprintf(npalog,"%s",oline);
          //  } useqmatom
          }
        }
        npa=0;
      }

      if( strspn(fline,"  Value of the MolecularEnergy:") > 30  && leave==0) { //!strcmp(txt1,enetxt1) && !strcmp(txt2,enetxt2) && !strcmp(txt3,enetxt3) && !strcmp(txt4,enetxt4) ) { 
	sscanf(fline,"%s %s %s %s %f",dum1,dum2,dum3,dum4,&QMener);
        leave--;
      }

      if( !( strstr(fline,"nuclear repulsion energy =")==NULL ) && (inrep < 2) ){
        sscanf(fline,"%s %s %s %s %s",str[0],str[1],str[2],str[3],str[4]);
        qm->Erep[10-inrep] = atof(str[4]); //this removed a nan error ...
        inrep++;
      }

    }
  }
  
  if(leave){
    fprintf(stderr ,"\n!!QM DIDNT PROPERLY CONVERGE!!\nEXTRACTED ENERGY AND GRADIENT::\nENERGY=%10.6f TOL=%10.6e REM=%d\n",QMener,deltaE,HL);
    FAIL=1;
    for(i=0;i<nrCOL;i++)
      if((mm->MMcluster[I*nrCOL+i]) &&  qm->includeQ>eIncludeNONE)
	HL--;
    for(i=0;i<qm->nrQMatoms;i++)
      if(qm->QMcluster[2*i]==I)
	HL--;
   /* if(deltaE>=1e-7){
      FAIL=1;
      qm->dI=I+1;
    }
   */
  }

  if(!FAIL){ //WE GOT OUR DATA FROM A TEXT FILE SO SANITY CHECKS MUST BE CONDUCTED
    for(i=0;i<nrQM;i++){
      j=qm->QMcluster[2*i+1];
      for(k=XX;k<=ZZ;k++){
        QMgrad[j][k]=0.0;   
        QMgrad[j+nrQM][k]=0.0;
      }  
      //if(qm->useQMatom[j]){
      if( qm->QMcluster[2*i]==I ) { // IS THIS ATOM PART OF THIS CLUSTER
        DN=(double)norm(Qgrad[j]);
        if( !(isnan(DN)) || !( isinf(DN)==0 ) ){ //SANITY CHECK
          for(k=XX;k<=ZZ;k++){
            if(qm->bQMexclusion)
              QMgrad[j][k]=Qgrad[j][k] + ( Qgrad[j][k] + (Ngrad[k])/(qm->Nrep[0]) ); // this looks odd but is correct
            else
              QMgrad[j][k]=Qgrad[j][k];
          }
        }else{
          fprintf(stderr,"\nWARNING FORCES ARE NaN OR INF IN QM ROUTINE\n");
        }
      }
      if(mm->MMcluster[(i+nrMM)+I*nrCOL]){ // IS THIS ATOM PART OF THIS CLUSTER
        DN=(double)norm(QQgrad[i]);
        if( !(isnan(DN)) || !( isinf(DN)==0 ) ){ //SANITY CHECK
          for(k=XX;k<=ZZ;k++)
            QMgrad[i+nrQM][k]=QQgrad[i][k];
        }else{
          fprintf(stderr,"\nWARNING FORCES ARE NaN OR INF IN QM ROUTINE\n");
        }
      }
      //}useqmatom
    }
    if( nrMM>0 && qm->includeQ>eIncludeNONE){
      for(i=0;i<nrMM;i++){
        for(k=XX;k<=ZZ;k++){
          MMgrad[i][k]=0.0;
        }
        if((int)mm->MMcluster[I*nrCOL+i]){ // IS THIS ATOM PART OF THIS CLUSTER
          DN=(double)norm(MMgrad[i]);
          if( !(isnan(DN)) || !( isinf(DN)==0 ) ){ //SANITY CHECK
            for(k=XX;k<=ZZ;k++)
              MMgrad[i][k]=Mgrad[i][k];  
          }else{
            fprintf(stderr,"\nWARNING FORCES ARE NaN OR INF IN QM ROUTINE\n");
          }
        }
      }
    } 
    qm->Erep[0]  += qm->Erep[10];
    qm->Erep[1]  += qm->Erep[9];
    qm->Erep[2]  += (qm->Erep[10]-qm->Erep[9])*( ((float)(qm->Nrep[1]))/((float)(qm->Nrep[0]))  );
    qm->Erep[3]  += ((float)(qm->Nrep[1]))/((float)(qm->Nrep[0]));
  }
  
  qm->dI=FAIL; 

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
  int i,I,J,j,RV=0,cl=1,NFILES=1,NFILER=1,GWF,k,l,a,aa,b,bb,qq,mnr;
  real QMener=0.0,QME0=0.0,rij,r12,r6,r2,C6ij,C12ij,DIener=0.0,theFORCE,DIe1=0.0,DIe2=0.0;
  rvec *QMgrad,*MMgrad,dist[2];
  char *exe;
  real *qmclusters,*mmclusters;
  FILE *stdver;
  FILE *npalog,*fptr;

  /* Allocate */
  //  fprintf(stderr,"allocating more memory\n");
  snew(QMgrad,2*(qm->nrQMatoms+1));
  snew(MMgrad,mm->nrMMatoms);
  snew(qm->QMexec,qm->nrQMc+1);

  qm->QMexec[0]=qm->nrQMc;
  for(i=0;i<qm->QMexec[0];i++){
    qm->QMexec[i+1]=i;
  }

  mnr=mm->nrMMatoms;
  qq=qm->nrQMatoms;
  aa=(qq*2+mnr);
  bb=(qq+mnr);

  //init
  for(i=0;i<11;i++){
    qm->Erep[i]=0.0;
  }
  for(i=0;i<qq;i++){
    for(j=XX;j<=ZZ;j++){
      QMgrad[i][j]=0.0;
      QMgrad[i+qq][j]=0.0;
    }
  }
  for(i=0;i<mnr;i++){
    MMgrad[i][XX]=0.0; MMgrad[i][YY]=0.0; MMgrad[i][ZZ]=0.0;
  }

  qm->dI=0;
  
  J=(qm->nrQMc==1)?(1):(qq);
  for(i=0;i<J*aa;i++){
    for(j=XX;j<=ZZ;j++){
      f[i][j]=0.0; 
      fshift[i][j]=0.0;
    }
  }

  write_mpqc(fr,qm,mm,0);

  NFILES=qm->nrQMc; 
  if(qm->bPrint_npa){
    npalog=fopen("npa.log","a");
    fprintf(npalog,"@%d\n",qm->nrQMc);
    fclose(npalog);
  }  

  RV=do_mpqc(qm);

  qm->QMexec[0]=0;
  for(I=0;I<qm->nrQMc;I++){
    qm->dI=0;
    QME0 = read_mpqc(QMgrad,MMgrad,qm,mm,I);
    if(qm->dI==0){            //SUCCESS
      QMener+=QME0;  
      //assign
      for(i=0;i<qm->nrQMatoms;i++){
        for(j=0;j<DIM;j++){
          f     [i+I*aa][j]    = HARTREE_BOHR2MD*QMgrad[i][j];
          fshift[i+I*aa][j]    = HARTREE_BOHR2MD*QMgrad[i][j];
          f     [i+I*aa+bb][j] = HARTREE_BOHR2MD*QMgrad[i+qq][j];
          fshift[i+I*aa+bb][j] = HARTREE_BOHR2MD*QMgrad[i+qq][j];
          QMgrad[i][j]=0.0; QMgrad[i+qq][j]=0.0;
        }
      }
      for(i=0;i<mm->nrMMatoms;i++){
        for(j=0;j<DIM;j++){
          f     [i+qq+I*aa][j] = HARTREE_BOHR2MD*MMgrad[i][j];      
          fshift[i+qq+I*aa][j] = HARTREE_BOHR2MD*MMgrad[i][j];
          MMgrad[i][j]=0.0;
        }
      }
    }
    else{
      //build new array with failed cases
      qm->QMexec[0]++;
      qm->QMexec[ qm->QMexec[0] ] = I;  
    }
  }
  if(qm->QMexec[0]>0){ //THIS MEANS WE HAVE FAILED CASES
    write_mpqc(fr,qm,mm,1);
    fprintf(stdout,"\nMPQC FAILED ONCE\nCASE FILES: ");
    for(i=0;i<qm->QMexec[0];i++){
      fprintf(stdout,"%d ",qm->QMexec[i+1]);
    }
    fprintf(stdout," |%d|\nRERUNNING FAILED CASES\n",qm->QMexec[0]);
    RV=do_mpqc(qm);
    qm->QMexec[0]=0;
    for(I=0;I<qm->QMexec[0];I++){
      qm->dI=0;
      QME0 = read_mpqc(QMgrad,MMgrad,qm,mm,qm->QMexec[I+1]);
      if(qm->dI==0){            //SUCCESS
        QMener+=QME0;  
        //assign
        for(i=0;i<qm->nrQMatoms;i++){
          for(j=0;j<DIM;j++){
            f     [i+qm->QMexec[I+1]*aa][j]    = HARTREE_BOHR2MD*QMgrad[i][j];
            fshift[i+qm->QMexec[I+1]*aa][j]    = HARTREE_BOHR2MD*QMgrad[i][j];
            f     [i+qm->QMexec[I+1]*aa+bb][j] = HARTREE_BOHR2MD*QMgrad[i+qq][j];
            fshift[i+qm->QMexec[I+1]*aa+bb][j] = HARTREE_BOHR2MD*QMgrad[i+qq][j];
            QMgrad[i][j]=0.0; QMgrad[i+qq][j]=0.0;
          }
        }
        for(i=0;i<mm->nrMMatoms;i++){
          for(j=0;j<DIM;j++){
            f     [i+qq+qm->QMexec[I+1]*aa][j] = HARTREE_BOHR2MD*MMgrad[i][j];      
            fshift[i+qq+qm->QMexec[I+1]*aa][j] = HARTREE_BOHR2MD*MMgrad[i][j];
            MMgrad[i][j]=0.0;
          }
        }
      }
      else{
        gmx_fatal(FARGS,"MPQC FAILED TWICE\n");   //THIS MEANS WE FAILED AGAIN  
      }
    }
  }

  qq=aa; mnr=bb;
  DIener=0.0;
  // NOW ADD DISPERSION ON Q ATOMS
  // REC.: DO NOT USE IT UNLESS YOUR QM MOLECULE IS ONE HUGE THING IN THE TOPOLOGY FILE
  //       IF THIS IS THE CASE THE BELOW SHOULD BE TESTED
  if(qm->bQMdispersion){
    for(k=0;k<qm->nrQMatoms;k++){         // ADD Q ATOMS                   
      a=qm->QMcluster[k*2];               // CLUSTER INDEX a:(0 -> NC-1)   
      b=qm->QMcluster[k*2+1];             // QM INDEX      b:(0 -> nrQM-1) 
      for(l=0;l<qm->nrQMatoms;l++){       // DISTANCE BETWEEN POINT AND THE CLUSTER THIS QM ATOM BELONGS TO 
        aa=qm->QMcluster[l*2];
        bb=qm->QMcluster[l*2+1];
        if((qm->QMcA[b*(qm->nrQMatoms)+bb]<=(qm->fQMelrc)) && (a!=aa) && (b!=bb) ) {	// IS THIS A QM ATOM THAT BELONG TO A DIFFERENT CLUSTER ?
          rvec_sub(qm->rQM[a*(qm->nrQMatoms)+b],mm->rMM[a*(mm->nrMMatoms+qm->nrQMatoms)+mm->nrMMatoms+bb],dist[0]);
          rij=qm->QMcA[b*(qm->nrQMatoms)+bb]; 
          r2 =rij*rij;
          r6 =r2*r2*r2;
          r12=r6*r6;
          C6ij=sqrt(qm->c6[b]*qm->c6[bb])/r6;      // DIVISION WITH LENGTH DONE HERE 
          C12ij=sqrt(qm->c12[b]*qm->c12[bb])/r12;  // AND HERE 
          DIe1+=C12ij;
          DIe2+=C6ij;
          for(j=XX;j<=ZZ;j++){
            theFORCE      = ( 12.0*C12ij - 6.0*C6ij  ) * ( dist[0][j] )/r2;
            f[b+a*qq+mnr][j]      += theFORCE;  //  test
            fshift[b+a*qq+mnr][j] += theFORCE;	//  FOR SR VIRIAL 
          }
          DIener += C12ij-C6ij;
        }
      }
    }
  }
  //END DISP.
  // BEGIN ELECTROEXCLUDE
  if(0){ //qm->bQMexclusion){
    for(k=0;k<qm->nrQMatoms;k++){
      a=qm->QMcluster[k*2];               // CLUSTER INDEX a:(0 -> NC-1)   
      b=qm->QMcluster[k*2+1]; 
      for(l=0;l<qm->nrQMatoms;l++){       // DISTANCE BETWEEN POINT AND THE CLUSTER THIS QM ATOM BELONGS TO 
        aa=qm->QMcluster[l*2];
        bb=qm->QMcluster[l*2+1];
        if( (qm->QMcA[b*(qm->nrQMatoms)+bb]<=(qm->fQMelrc)) && (a!=aa) && (b!=bb) ) { // b>bb
          rvec_sub(qm->rQM[a*(qm->nrQMatoms)+b],mm->rMM[a*(mm->nrMMatoms+qm->nrQMatoms)+mm->nrMMatoms+bb],dist[0]);
          rij=norm(dist[0]);
          if (rij!=qm->QMcA[b*(qm->nrQMatoms)+bb])
            fprintf(stderr,"\nWARNING::COORDINATE PANIC\n");
//          rij = (qm->QMcA[b*(qm->nrQMatoms)+bb]); // nanometer to bohr
          r2  = rij*rij;
          r6  = (qm->pcharges[b])*(qm->pcharges[bb]);
          DIener -= (r6/rij)*(qm->fQMexclusion)*FEL*0.5; // one sided contriubtion
          for(j=XX;j<=ZZ;j++){
            theFORCE      = -1.0*(qm->fQMexclusion)*( FEL*r6/r2  )*( dist[0][j] )/rij;
            f[b+a*qq+mnr][j]       += theFORCE;      //  test
            fshift[b+a*qq+mnr][j]  += theFORCE;	     //  FOR SR VIRIAL 
/*            f[bb+a*qq+mnr][j]      += -1.0*theFORCE; //  test
            fshift[bb+a*qq+mnr][j] += -1.0*theFORCE; //  FOR SR VIRIAL 
*/        }
        }
      } 
    }
  }
  // END ELECTROEXCLUDE
  // END Big loop 
  // put the QMMM forces in the force array and do the fshift 

  fptr=fopen("qmmm_ene.dat","a");
  fprintf(fptr,"% 20.10f % 20.10f % 20.10f % 20.10f % 20.10f % 6d % 20.10f\n",QMener,qm->Erep[2],qm->Erep[0],qm->Erep[1],qm->Erep[3],NFILES,DIener);
  fclose(fptr);

  free(QMgrad);
  free(MMgrad);
  free(qm->QMexec);

  QMener = QMener*HARTREE2KJ*AVOGADRO+DIener;
  return(QMener);
}

int do_mpqc(t_QMrec *qm)
{
  //NON PIPED IMPLEMENTATION OF CALL
  pid_t pid;
  int   rv,step_id;
  int	commpipe[2],stdcopy[2];		   /* This holds the fd for the in & out of the pipe */
  FILE  *in,*out;
  char  input_string[MPQCLEN*4],exename[MPQCLEN*2],execute[10],eo0[20],eo1[20],eo2[20],eo3[20],eo4[20],tfi[100],eo5[20];
  char  options[1000],launch[1000],rmout[200],eo6[20];

  //sprintf(eo0,"-o=run_mpqc%d.out",I);
  if(qm->nQMnodes>0){
    sprintf(eo1,"--nnodeperjob=%d",qm->nQMnodes); 
  }
  else{
    sprintf(eo1," ");
  }
  if(qm->nQMcpus>0){
    sprintf(eo2,"--nprocperjob=%d",qm->nQMcpus);
  }
  else{
    sprintf(eo2," ");
  }
  if(qm->nQMproc>0){
    sprintf(eo5,"--nprocpernode=%d",qm->nQMproc);
  }
  else{ 
    sprintf(eo5," ");
  }
  if(qm->nQMthread>0){
    sprintf(eo6,"--nthreadperproc=%d",qm->nQMthread);
  }
  else{
    sprintf(eo6," ");
  }

  if(!(strstr(qm->gauss_exe,"mpqcrun")==NULL)){
    sprintf(eo3,"--autoout");
    sprintf(eo4,"--simpout");
  }
  else{
    sprintf(eo3," "); 
    sprintf(eo4,"-o=run_mpqc0.out");
  }
  sprintf(options,"$OPTIONS");
  sprintf(launch,"--launch=\"$LAUNCH\"");

  sprintf(exename,"--rerun ");
  for(rv=0;rv<qm->QMexec[0];rv++){
    sprintf(rmout," run_mpqc%d.inp ",qm->QMexec[rv+1]);
    strcat(exename,rmout);
  }
  sprintf(input_string,"%s %s %s %s %s %s %s %s %s %s",qm->gauss_exe,options,eo1,eo2,eo5,eo6,eo3,eo4,launch,exename);

#ifdef GMX_NO_SYSTEM
  printf("Warning-- No calls to system(3) supported on this platform.");
  gmx_fatal(FARGS,"Call to '%s' failed\n",input_string);
#else
  sprintf(tfi,"touch run_mpqc*.out*");  
  sprintf(rmout,"rm run_mpqc*.out*");
  system(tfi);
  system(rmout); //REMOVES OLD OUTPUT FILE
  if ( system(input_string) != 0 )
    gmx_fatal(FARGS,"Call to '%s' failed\n",input_string);
#endif
  
  return(0);
}

