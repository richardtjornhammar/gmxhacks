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
#include <errno.h>
#include "qmstack.h"

#define   sgn(a)        ((a)<0?(-1):(1))
#define   max(a,b)     ((a)>(b)?(a):(b))

double mround(double T)
{
   double RV=0;
   RV=sgn(T)*ceil(max(T,-1*T));
   return(RV);
}

int anint000(double r)
{
   int sign,res;
   double frac,part;
   
   sign=(r>=0)?1:-1;
   frac=modf(r,&part);
   res=(part+((frac*sign>0.5)?1.0:0.0)*sign);

   return(res);   

}

int anint001(double r)
{
   int sign,res;
   double frac,part;
   
   sign=(r>=0)?1:-1;
   frac=modf(r,&part);

   res=(part+((frac*sign>0.5)?1.0:0.0)*sign);

   return(res);   

}

int anint002(double T)
{
   int RV=0;
   RV=(int)(sgn(T)*ceil(max(T,-1*T)));
   return(RV);
}

int anint003(double r)
{
   int sign,res;
   double frac,part;
   
   sign=(r>=0)?1:-1;
   frac=modf(r,&part);

   res=(((frac*sign>0.5)?1.0:0.0)*sign);

   return(res);  
}

double roundChrg(double chrg)
{
  double rchrg,schrg,cchrg;
  
  schrg=(chrg>=0.0)?(1.0):(-1.0);
  cchrg=(chrg>(-1*chrg))?(chrg):(-1*chrg);
  rchrg=ceil(cchrg)*schrg;
  
  return(rchrg);
}


void     qmsinit(qmstack *stk)
{
  stk -> cnt = 0;
  stk -> top = NULL;
  stk -> bottom = NULL;
}

void     qmspush(qmsdata d, qmstack *stk)
{
  elem *p,*q,*t;
  int i,N;

  p = malloc(sizeof(elem));
  p -> d = d;

  p -> next = stk -> top;
  stk -> top = p;
  stk -> cnt++;

  t = p;
  N = stk -> cnt;
  for(i=0;i<N-1;i++){
    q = t;
    t = t->next;
    t -> prev = q;
  }
  stk -> bottom = t;
}

void     qmsevac(qmstack *stk)
{
  elem *p;
  
  while( (stk -> cnt) >= 1){
    p = stk -> top;
    stk -> top = stk -> top -> next;
    stk -> cnt--;
    free(p);
  }
}

qmsdata     qmspop(qmstack *stk)
{
  qmsdata d;
  elem *p;
  
  d = stk -> top -> d;
  p = stk -> top;
  stk -> top = stk -> top -> next;
  stk -> cnt--;
  free(p);
  return d;
}

qmsdata     qmstop(qmstack *stk)
{
  return(stk -> top -> d);
}

boolean  qmsempty(const qmstack *stk)
{
  return((boolean)(stk -> cnt == EMPTYSTACK));
}

boolean  qmsfull(const qmstack *stk)
{
  return((boolean)(stk -> cnt == FULLSTACK));
}


void printClustering(t_QMrec *qm,FILE *fpo)
{
  int i,j,k,N,C=-1;

  N=qm->nrQMatoms;

  for(i=0;i<N;i++){
    if(!(C==(qm->QMcluster[i*2]))){
      C=qm->QMcluster[i*2];
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
     if(qm->QMcluster[ (sb - 1)*2] > qm->QMcluster[ sb*2 ]) {
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

void qmmmPBC(t_QMrec *qm, rvec R[2])
{
  int i,j,k,l,m,n,AA;
  real tmp=0.0;

  AA=( (qm->ePBC==epbcXY)?(YY):((qm->ePBC==epbcXYZ) ?(ZZ):(XX-1)) );

  for(j=XX;j<=AA;j++){
      R[1][j]=anint001((double)(R[0][j]/qm->box[j][j]))*(-1.0);
      R[0][j]+=qm->box[j][j]*R[1][j];
  }
}

void qmPBCsimple(t_QMrec *qm, rvec R[2])
{
  int i,j,k,l,m,n,AA;
  real tmp=0.0;

  AA=( (qm->ePBC==epbcXY)?(YY):((qm->ePBC==epbcXYZ) ?(ZZ):(XX-1)) );

  for(j=XX;j<=AA;j++){
      R[1][j]=(R[0][j]>0.5*qm->box[j][j])?(-1.0):( (R[0][j]>-0.5*qm->box[j][j])?(1.0):(0.0) );
      R[0][j]+=qm->box[j][j]*R[1][j];
  }
}

void clusterAlgorithm(t_QMrec *qm,FILE *fpo)
{
  int i,j,k,l,m,N,q;
  int flag,C,oC;
  int *nvisi;

   qmstack real_s,*s;
   qmstack real_NN,*NN;
   qmstack real_ndx,*ndx;
   N=qm->nrQMatoms;
   nvisi=malloc(sizeof(int)*N);

   s=&real_s;NN=&real_NN;ndx=&real_ndx;
   qmsinit(s);qmsinit(NN);qmsinit(ndx);

   for(i=0;i<N;i++){
     qmspush(i,ndx);
     nvisi[i]=i+1;
   }
   
   while(!qmsempty(ndx)){
     i=qmspop(ndx); qmsevac(NN);
     if(nvisi[i]>0){
       C--;
       for(j=0;j<N;j++){
	 if(qm->QMcB[i*N+j]){
	   qmspush(j,NN);
	 }
       }
       while(!qmsempty(NN)){
	 k=qmspop(NN);
	 nvisi[k]=C;
	 for(j=0;j<N;j++){
	   if(qm->QMcB[k*N+j]){
	     for(q=0;q<N;q++){
	       if(nvisi[q]==j+1){
		 qmspush(q,NN);
	       }
	     }
	   }
	 }
       }
     }
   }
   for(i=0;i<N;i++){
     qm->QMcluster[i*2+1]=i;
     qm->QMcluster[i*2]=nvisi[i]-C;
   }
  sortClusters(qm);

  fprintf(fpo,"\n\n");
  C=0; oC=qm->QMcluster[0];
  for(i=0;i<N;i++){
    if( !(qm->QMcluster[2*i] == C) && !(qm->QMcluster[2*i] == oC) ){
      C++;
      oC=qm->QMcluster[2*i];
    }
    qm->QMcluster[2*i]=C;
    fprintf(fpo,"Index::% 4d Cluster::% 4d \n",qm->QMcluster[2*i+1],(qm->QMcluster[2*i]));
  }
  fprintf(fpo,"\n\n");  
  qm->nrQMc=C+1;

  free(nvisi);
}

void makeDistance(t_QMrec *qm, FILE *fpo)
{
  int   i,j,k,l,N,C,oC,flag,NP,S;
  float RL=0.3,R1;
  rvec  dist[2];
  t_pbc  pbc;

  N=qm->nrQMatoms;

  fprintf(fpo,"%Distance matrix from gromacs\nN=%d; A=[\n",N);
  for(i=0;i<N;i++){
    qm->QMcV[i]=i+1;
    for(j=0;j<N;j++){
      rvec_sub(qm->xQM[i],qm->xQM[j],dist[0]);     
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
  if(mm->nrMMatoms>0){
    for(i=0;i<mm->nrMMatoms;i++){ //MM
      mm->MMcluster[i]=1;
    }
  }
  for(i=0;i<qm->nrQMatoms;i++){ //QA
    mm->MMcluster[i+mm->nrMMatoms]=0;
  }
}

void clusterMM(t_MMrec *mm, t_QMrec *qm)
{
  // SIMPLE DISTANCE SEARCH FOR CLOSEST CLUSTER
  int    a,b,nrQMc=0,i,j,k,kk,l,nrQMa=0,nrMM,nrCOL,aa,bb,S;
  rvec   R[2],dx;
  double D[3],D0=100000000.0,RC=1,R0,R1;
  char   errormsg[200];

  /*
   * Motivation for RC choice 
   * RC=1.0 is ca 600 atoms in spherical volume around an atom
   * Strong polarisation effects occur within ca 0.6nm of an ionic center 
   * If the charge is not extreme
   */

  nrQMc=qm->nrQMc; 
  nrQMa=qm->nrQMatoms;
  nrMM=mm->nrMMatoms; 
  nrCOL=nrMM+nrQMa;

  for(j=0;j<nrMM;j++){
    D0=100000000.0; k=0;
    for(i=0;i<nrQMa;i++){  			/* CALCULATE DISTANCE BETWEEN POINT AND ALL QM ATOM CLUSTERS */
      a=qm->QMcluster[i*2];			/* CLUSTER INDEX a:(0 -> NC-1) */
      b=qm->QMcluster[i*2+1];			/* QM INDEX      b:(0 -> nrQMa-1) */
      rvec_sub(qm->xQM[b],mm->xMM[j],R[0]);
      R1=sqrt(pow(R[0][XX],2.0)+pow(R[0][YY],2.0)+pow(R[0][ZZ],2.0));
      qmmmPBC(qm,R);                            
      R0=sqrt(pow(R[0][XX],2.0)+pow(R[0][YY],2.0)+pow(R[0][ZZ],2.0));

      for(kk=XX;kk<=ZZ;kk++){
	mm->rMM[j+a*nrCOL][kk]=mm->xMM[j][kk]+R[1][kk]*(qm->box[kk][kk]); //CREATES A SHIFTED COORDINATE
//        mm->rMM[j+a*nrCOL][kk]=qm->rQM[b][kk]-R[0][kk];
//	mm->rMM[j+a*nrCOL][kk]=qm->xQM[b][kk]+dx[kk];
      }
      l = (R0<RC) ? 1:0 ;                 /* CUTOFF SEARCH */
      mm->MMcluster[nrCOL*a+j]     = l;   /* * */
      mm->MMcluster[nrCOL*nrQMc+j]+= l;   /* SUMS UP THE NUMBER OF TIMES AN ATOM IS CALLED */
      if(R0<0.1){
        sprintf(errormsg,"BUG IN MM COORD::CLUSTER=%6d,QM=%6d::MM=%6d,NUM=%6d::NDX=[%6d,%6d]",a,b,j,mm->MMcluster[nrCOL*nrQMc+j],mm->indexMM[j]+1,qm->indexQM[b]+1);
        gmx_fatal(FARGS,errormsg);
      }
      if(R0<D0){                          /* SEARCH FOR CLOSEST CLUSTER AND ASSIGN IT */
	k=a; D0=R0;
      }
    }
    mm->MMcluster[nrCOL*k+j] = 1; /* ENSURES THAT ALL ATOMS ARE PRESENT IN THE CLOSEST CLUSTER */
  }
  //NOW ASSIGN Q CHARGES FROM THE OTHER CLUSTERS QM ATOMS (DONE)
  for(j=0;j<nrQMa;j++){         /* ADD Q ATOMS                   */
    a=qm->QMcluster[j*2];       /* CLUSTER INDEX a:(0 -> NC-1)   */
    b=qm->QMcluster[j*2+1];     /* QM INDEX      b:(0 -> nrQM-1) */
    for(i=0;i<nrQMa;i++){       /* DISTANCE BETWEEN POINT AND THE CLUSTER THIS QM ATOM BELONGS TO */
      aa=qm->QMcluster[i*2];
      bb=qm->QMcluster[i*2+1];
      l=((qm->QMcA[b*nrQMa+bb]<=RC) && (a!=aa) )?1:0; //TURNED ON	
      for(kk=XX;kk<=ZZ;kk++){
	mm->rMM[b+nrMM+aa*nrCOL][kk]=qm->rQM[b][kk];
      }
      mm->MMcluster[(nrCOL*aa+(b+nrMM))] = l;
    }
  }
}

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
