#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "hlist.h"

void mMatVec(float R[3][3], float X[3], float T[3])
{
  T[0]=R[0][0]*X[0]+R[0][1]*X[1]+R[0][2]*X[2];
  T[1]=R[1][0]*X[0]+R[1][1]*X[1]+R[1][2]*X[2];
  T[2]=R[2][0]*X[0]+R[2][1]*X[1]+R[2][2]*X[2];
}

void rotate_top(LINK atom, LINK solv, float RP[3], int CM, float PSI[3])
{
  //ROTATES MOLECULE AROUND A POINT RP OR ITS CENTER OF MASS
  //PSI[0] ROT AROUND X, PSI[1] ROT AROUND Y, PSI[2] ROT AROUND Z
  int   i,j,k,l,m,n;
  int   lna=0,lns=0;
  float R[3],X[3],XX[3],ROT[3][3],T[3];
  float theta,phi;
  LINK  new;
 
  R[0]=0.0; R[1]=0.0; R[2]=0.0;
  X[0]=0.0; X[1]=0.0; X[2]=0.0;
  XX[0]=0.0; XX[1]=0.0; XX[2]=0.0;

  new=atom;
  lna=listlength(atom);
  for(i=0;i<lna;i++){
    R[0]+=new->data.pos[0];
    R[1]+=new->data.pos[1];
    R[2]+=new->data.pos[2];
    new=new->next;
    if(new==NULL)
      break;
  }
  lns=listlength(solv);
  if(CM){
    X[0]=R[0]/(float)lna; 
    X[1]=R[1]/(float)lna; 
    X[2]=R[2]/(float)lna; 
  }
  else{
    X[0]=RP[0]; 
    X[1]=RP[1]; 
    X[2]=RP[2]; 
  }
  
  //ASSEMBLE ROT MATRIX FROM EULER ANGLES
  
  ROT[0][0]=cos(PSI[0])*cos(PSI[2])-sin(PSI[0])*cos(PSI[1])*sin(PSI[2]);
  ROT[0][1]=sin(PSI[0])*cos(PSI[2])+cos(PSI[0])*cos(PSI[1])*sin(PSI[2]);
  ROT[0][2]=sin(PSI[1])*sin(PSI[2]);
  ROT[1][0]=-cos(PSI[0])*sin(PSI[2])-sin(PSI[0])*cos(PSI[1])*cos(PSI[2]);
  ROT[1][1]=-sin(PSI[0])*sin(PSI[2])+cos(PSI[0])*cos(PSI[1])*cos(PSI[2]);
  ROT[1][2]=sin(PSI[1])*cos(PSI[2]);
  ROT[2][0]=sin(PSI[1])*sin(PSI[0]);
  ROT[2][1]=-cos(PSI[0])*sin(PSI[1]);
  ROT[2][2]=cos(PSI[1]);
  new=atom;
  for(i=1;i<lna;i++){
    for(j=0;j<3;j++){
      XX[j]=new->data.pos[j]-X[j];
    }
    mMatVec(ROT,XX,T);
    new->data.pos[0]=T[0];
    new->data.pos[1]=T[1];
    new->data.pos[2]=T[2];
    new=new->next;
  }
  new=solv;
  for(i=1;i<lns;i++){
    for(j=0;j<3;j++){
      XX[j]=new->data.pos[j]-X[j];
    }
    mMatVec(ROT,XX,T);
    new->data.pos[0]=T[0];
    new->data.pos[1]=T[1];
    new->data.pos[2]=T[2];
    new=new->next;
  }
}
