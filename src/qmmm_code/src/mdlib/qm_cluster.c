/*
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
    fprintf(fpo,"\n@ %5d %5d",qm->QMcluster[i*2],qm->QMcluster[i*2+1]);
  }
  fprintf(fpo,"\n\n"); 
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

void clusterAlgorithm(t_QMrec *qm,FILE *fpo)
{
  int i,j,k,l,m,N,q;
  int flag,C=0,oC;
  int *nvisi;

   qmstack real_s,*s;
   qmstack real_NN,*NN;
   qmstack real_ndx,*ndx;
   N=qm->nrQMatoms;
   nvisi=malloc(sizeof(int)*(N+1));

   s=&real_s;NN=&real_NN;ndx=&real_ndx;
   qmsinit(s);qmsinit(NN);qmsinit(ndx);

   for(i=0;i<N;i++){
     qmspush(i,ndx);
     nvisi[i]=i+1;
   }
   C=0;
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
    fprintf(fpo,"@ Index::% 4d Cluster::% 4d \n",qm->QMcluster[2*i+1],(qm->QMcluster[2*i]));
  }

  qm->nrQMc=C+1;

  free(nvisi);
}

void makeDistance(t_QMrec *qm, FILE *fpo, t_pbc *pbc)
{
  int    i,j,k,l,N,C,oC,flag,NP,S;
  real   RL=0.3,R1,scale=1.0;
  rvec   dist[2];

  N = qm->nrQMatoms;
  RL= qm->fQMcutoff;

  fprintf(fpo,"\n\n@ Distance matrix from gromacs\n# N=%d; A=[ \n # ",N);
  for(i=0;i<N;i++){
    qm->QMcV[i]=i+1;
    for(j=0;j<=i;j++){
      qm->QMcIS[i*N+j]=pbc_dx_aiuc(pbc,qm->xQM[i],qm->xQM[j],dist[0]);
      qm->QMcA[i*N+j]=norm(dist[0]);
      scale=(qm->atomicnumberQM[i]==6||qm->atomicnumberQM[i]==7||qm->atomicnumberQM[j]==6||qm->atomicnumberQM[j]==7)?1.5:1.0;
      scale=(qm->atomicnumberQM[i]==16||qm->atomicnumberQM[i]==30||qm->atomicnumberQM[j]==16||qm->atomicnumberQM[j]==30)?1.7:scale;
      qm->QMcB[i*N+j]=(qm->QMcA[i*N+j]<=(RL*scale))?(1.0):(0.0);
      qm->QMcA[j*N+i]=qm->QMcA[i*N+j];
      qm->QMcB[j*N+i]=qm->QMcB[i*N+j];
      qm->QMcIS[i*N+j]=qm->QMcIS[j*N+i];
      fprintf(fpo," % 9.5f ",qm->QMcA[i*N+j]);
    }
    fprintf(fpo,"\n # ");
  }
  fprintf(fpo,"\n# ];\n");
  fprintf(fpo,"@ Boolean matrix from gromacs N=%d \n B=[\n",N);
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

void clustering(t_QMrec *qm, t_pbc *pbc)
{
  int   i,j,k,l,N,C,oC,flag,NP;
  float RL=0.3;
  rvec  dist[2];
  FILE  *fpo;

  fpo=fopen("outpclust.m","w");
  fprintf(fpo,"%s","\n@ OPENED CLUSTERING\n\n");
  N=qm->nrQMatoms;
  makeDistance(qm, fpo, pbc);     // CREATE DISTANCE MATRIX

  if( (qm->nQMclusterupdate != -1) && ( (qm->nQMclusterupdate == qm->QMstep) || (qm->QMstep==0) ) ){
    fprintf(fpo,"%s","\n@ PERFORMING C.ALG\n");
    clusterAlgorithm(qm,fpo);     // DO CLUSTER SEARCH
    if(qm->nQMclusterupdate==0){
      qm->nQMclusterupdate=-1;
    }
    if( qm->nQMclusterupdate == qm->QMstep ){
      qm->QMstep=0;
    }
    fprintf(fpo,"\n@ FOUND %6d CLUSTER(S)\n\n",qm->nrQMc);
    printClustering(qm,fpo);
  }else{
    if( (qm->nQMclusterupdate == -1) && qm->QMstep == 0) {
      qm->nrQMc=qm->QMcluster[2*(qm->nrQMatoms-1)]+1; 
      fprintf(fpo,"\n@ WAS GIVEN %6d CLUSTER(S)\n\n",qm->nrQMc);
      sortClusters(qm);
      printClustering(qm,fpo);
    }
  }
  qm->QMstep+=1;

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
      mm->MMcluster[i]=0;
    }
  }
  for(i=0;i<qm->nrQMatoms;i++){ //QA
    mm->MMcluster[i+mm->nrMMatoms]=0;
  }
}

void clusterMM(t_MMrec *mm, t_QMrec *qm, t_pbc *pbc)
{
  // SIMPLE DISTANCE SEARCH FOR CLOSEST CLUSTER
  int    a,b,nrQMc=0,i,j,k,kk,l,nrQMa=0,nrMM,nrCOL,aa,bb,S,a0=-1;
  rvec   R[2],dx,*rcen;
  double D[3],D0=100000000.0,RC=1,R0,R1;
  char   errormsg[200];
  int    itmp[10];
  double dtmp[10];
  rvec   rtmp[2];

  nrQMc=qm->nrQMc; 
  nrQMa=qm->nrQMatoms;
  nrMM=mm->nrMMatoms; 
  nrCOL=nrMM+nrQMa;
  RC=qm->fQMelrc;

  snew( rcen,(nrQMc+1) );

  //init
  for(i=0;i<nrQMc;i++){
    for(j=0;j<nrQMa;j++){
      qm->QMcIS[i*nrQMa+j]=qm->shiftQM[j];
      mm->MMcIS[i*nrCOL+j+nrMM]=qm->shiftQM[j];
      mm->MMcluster[(nrCOL*i+j+nrMM)]=0;
    }
    for(j=0;j<nrMM;j++){
      mm->MMcIS[i*nrCOL+j]=mm->shiftMM[j];
      mm->MMcluster[nrCOL*i+j] = 0;
    }
  }

  a0=-1;
  for(i=0;i<nrQMa;i++){  	 // CALCULATE DISTANCE BETWEEN POINT AND ALL QM ATOM CLUSTERS 
    a=qm->QMcluster[i*2];	 // CLUSTER INDEX a:(0 -> NC-1)    monotonically increasing
    b=qm->QMcluster[i*2+1];	 // QM INDEX      b:(0 -> nrQMa-1)        not -"-
    if(!(a==a0)){                // If we found a new cluster. This means that the first atom of a cluster should be the center
      a0++;
      for(kk=XX;kk<=ZZ;kk++){
        rcen[a0][kk]=qm->xQM[b][kk];
      }
      qm->QMcIS[a0*nrQMa] = XYZ2IS(rcen[a0][XX],rcen[a0][YY],rcen[a0][ZZ]);
    }
  }

  for(a=0;a<nrQMc;a++){                         // CALCULATE A QM ATOM POSITION RELATIVE A CLUSTER CENTER                   
    for(i=0;i<nrQMa;i++){                       // DISTANCE BETWEEN POINT AND THE CLUSTER 
      qm->QMcIS[a*nrQMa+i]=pbc_dx_aiuc(pbc,rcen[a],qm->xQM[i],R[0]); 
      for(kk=XX;kk<=ZZ;kk++){
        qm->rQM[a*nrQMa+i][kk]=rcen[a][kk]-R[0][kk];
      }
    }
  }

  if(nrMM>0){
    for(j=0;j<nrMM;j++){                        // MM ATOMS
      D0=100000000.0; k=0; a0=-1;
      for(i=0;i<nrQMa;i++){  			// CALCULATE DISTANCE BETWEEN POINT AND ALL QM ATOMS 
        a=qm->QMcluster[i*2];			// CLUSTER INDEX a:(0 -> NC-1)    
        b=qm->QMcluster[i*2+1];			// QM INDEX      b:(0 -> nrQMa-1)  
        itmp[0] = pbc_dx_aiuc(pbc, rcen[a], mm->xMM[j], R[0]);
        for(kk=XX;kk<=ZZ;kk++)   
          rtmp[0][kk]= rcen[a][kk]-R[0][kk];
        kk = pbc_dx_aiuc(pbc, qm->xQM[b], mm->xMM[j], R[0]); 
        R0 = norm(R[0]);
        l = (R0<RC) ; 
        itmp[1] = (mm->MMcluster[nrCOL*a+j] == 0) ;
        if( itmp[1] && l ) { // if we found a new entry then we add it //R0
          mm->fcharges[nrCOL*a+j] = qm->fQMscreenS-qm->fQMscreenF*exp(-1.0*qm->fQMscreenK*pow( R0 - qm->fQMscreenP , qm->fQMscreenPow) );
          mm->MMcluster[nrCOL*a+j] = 1 ;
          mm->MMcIS[a*nrCOL+j]= itmp[0];  
          for(kk=XX;kk<=ZZ;kk++)                  // calculate MM atom position with respect to the cluster center a
            mm->rMM[a*nrCOL+j][kk] = rtmp[0][kk];
        }
        if(R0<D0){                              // SEARCH FOR CLOSEST CLUSTER AND ASSIGN IT 
          k=a; D0=R0;
        }
      }
//      if(qm->includeQ>eIncludeQ)
//        mm->MMcluster[nrCOL*k+j] = 0;     // NULLS OUT ENTRIES
    }
  }

  if(qm->nrQMc > 1){
  a0=-1;
  for(j=0;j<nrQMa;j++){                         // ADD Q ATOMS                   
    a=qm->QMcluster[j*2];                       // CLUSTER INDEX a:(0 -> NC-1)   (we sorted them earlier)
    b=qm->QMcluster[j*2+1];                     // QM INDEX      b:(0 -> nrQM-1) 
    mm->MMcIS[a*nrCOL+b+nrMM] = qm->QMcIS[a*nrQMa+b]; 
    for(i=0;i<nrQMa;i++){
      aa=qm->QMcluster[i*2];
      bb=qm->QMcluster[i*2+1];
      l=((qm->QMcA[b*nrQMa+bb] <= RC) && (a != aa) )?1:0;                         // TURNED ON
      itmp[1]=(mm->MMcluster[(nrCOL*a+(bb+nrMM))] == 0);
      if( l && itmp[1] ){                                                           // ADD THE Q ATOM once qm->QMcA[b*nrQMa+bb]
        mm->fcharges[(nrCOL*a+(bb+nrMM))]  = qm->fQMscreenS-qm->fQMscreenF*exp(-1.0*qm->fQMscreenK*pow( qm->QMcA[b*nrQMa+bb] - qm->fQMscreenP , qm->fQMscreenPow) );
        mm->MMcluster[(nrCOL*a+(bb+nrMM))] = 1;
        mm->MMcIS[bb+nrMM+a*nrCOL]=pbc_dx_aiuc(pbc,rcen[a],qm->xQM[bb],R[0]);
        for(kk=XX;kk<=ZZ;kk++)
          mm->rMM[bb+nrMM+a*nrCOL][kk]=rcen[a][kk]-R[0][kk];  
      }
/* else{
        mm->MMcIS[bb+nrMM+a*nrCOL] = qm->QMcIS[a*nrQMa]; //qm->QMcIS[a*nrQMa]; //these will never be used
      }         	*/						
    }
  }
  }

  free(rcen);
}

