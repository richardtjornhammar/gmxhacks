#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "hlist.h"

char* getAtomName(int i)
{
  char* atomName[120];

  atomName[0]="LA";
  atomName[1]="H "; 
  atomName[2]="He";
  atomName[3]="Li";
  atomName[4]="Be";
  atomName[5]="B ";
  atomName[6]="C ";
  atomName[7]="N ";
  atomName[8]="O ";
  atomName[9]="F ";
  atomName[10]="Ne";
  atomName[11]="Na";
  atomName[12]="Mg";
  atomName[13]="Al";
  atomName[14]="Si";
  atomName[15]="P ";
  atomName[16]="S ";
  atomName[17]="Cl";
  atomName[18]="Ar";
  atomName[19]="K ";
  atomName[20]="Ca";
  atomName[21]="Sc";
  atomName[22]="Ti";
  atomName[23]="V ";
  atomName[24]="Cr"; 
  atomName[25]="Mn";
  atomName[26]="Fe";
  atomName[27]="Co";
  atomName[28]="Ni";
  atomName[29]="Cu";
  atomName[30]="Zn";
  atomName[31]="Ga";
  atomName[32]="Ge";
  atomName[33]="As";
  atomName[34]="Se";
  atomName[35]="Br";
  atomName[36]="Kr";
  atomName[37]="Rb";
  atomName[38]="Sr";
  atomName[39]="Y ";
  atomName[40]="Zr";
  atomName[41]="Nb";
  atomName[42]="Mo";
  atomName[43]="Tc";
  atomName[44]="Ru";
  atomName[45]="Rh";
  atomName[46]="Pd";
  atomName[47]="Ag";
  atomName[48]="Cd";
  atomName[49]="In";
  atomName[50]="Sn";
  atomName[51]="Sb";
  atomName[52]="Te";
  atomName[53]="I ";
  atomName[54]="Xe";
  atomName[55]="Cs";
  atomName[56]="Ba";
  atomName[57]="La";
  atomName[58]="Ce";
  atomName[59]="Pr";
  atomName[60]="Nd"; 
  atomName[61]="Pm";
  atomName[62]="Sm";
  atomName[63]="Eu";
  atomName[64]="Gd";
  atomName[65]="Tb";
  atomName[66]="Dy";
  atomName[67]="Ho";
  atomName[68]="Er";
  atomName[69]="Tm";
  atomName[70]="Yb";
  atomName[71]="Lu";
  atomName[72]="Hf";
  atomName[73]="Ta";
  atomName[74]="W ";
  atomName[75]="Re";
  atomName[76]="Os";
  atomName[77]="Ir";
  atomName[78]="Pt";
  atomName[79]="Au";
  atomName[80]="Hg";
  atomName[81]="Tl";
  atomName[82]="Pb";
  atomName[83]="Bi";
  atomName[84]="Po";
  atomName[85]="At";
  atomName[86]="Rn";
  atomName[87]="Fr";
  atomName[88]="Ra";
  atomName[89]="Ac";
  atomName[90]="Th";
  atomName[91]="Pa";
  atomName[92]="U ";
  atomName[93]="Np";
  atomName[94]="Pu";
  atomName[95]="Am";
  atomName[96]="Cm";
  atomName[97]="Bk";
  atomName[98]="Cf";
  atomName[99]="Es";
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

int getAtomNr(char* Z){
  int i=0;
  char Zp[2],Zb[2];
  size_t n=2;

  for(i=0;i<=111;i++){
    sprintf(Zb,"%s", getAtomName(i));
    strncpy(Zp,Zb,n);
    if(Z[0] == Zp[0])
      if(i==1 | i==5 | i==6 | i==7 | i==8 | i==9 | i==15 | i==16 | i==19 | i==23 | i==39 | i==53 | i==92)
	break;
      else
	if(toupper(Z[1])==toupper(Zp[1]))
	  break;
  }
  return i;
}


void search_top(LINK hash[], int hlen, float CENTER[], int T[],
		float CRYST[], float DISPR[], MAT RES, float r, PAR CAL)
{
  int i,bb,j,m,sa,y,k,znr,z;
  int h=0,c=0,o=0,n=0,s=0;
  char Z;
  LINK at;
  float S[DIM3];
  float rpres;

  j=CAL->dat.B[0];
  m=CAL->dat.B[1];
  bb=CAL->dat.B[3];
  znr=CAL->dat.B[4];
  h=CAL->dat.B[5];
  c=CAL->dat.B[6];
  n=CAL->dat.B[7];
  o=CAL->dat.B[8];
  s=CAL->dat.B[10];

  for(i=1 ; i <= hlen; i++) {

    at=hash[i-1];
    for(sa=0;sa<3;sa++) {
      S[sa]=at->data.pos[sa]-CENTER[sa];
      switch (T[sa]) {
      case -1://shoot over box-dim
	if(at->data.pos[sa]<=DISPR[sa])
	  S[sa]+=CRYST[sa];
	break;
      case 1://shoot under box-dim
	if(at->data.pos[sa]>=(CRYST[sa]-DISPR[sa]))
	  S[sa]-=CRYST[sa];
      case 0:break;
      default: printf("ERROR2");
      }
    }
    rpres=sqrt( pow(S[0],2) + pow(S[1],2) + pow(S[2],2) );

    Z=at->data.name[0];
    switch(Z){
    case 'H':z=5; y=1;
      break;
    case 'C':z=6; y=1;
      break;
    case 'N':z=7; y=1;
      break;
    case 'O':z=8; y=1;
      break;
    case 'S':z=10; y=1;
      break;
    default: y=0;
    }
    //REGISTERS ALL HITS FOR HISTOGRAMS
    if( rpres <= r & y){
      RES->dat.H[0][i]+=1;
      if(!(RES->dat.H[11][i])){
	//ONLY ONE CONTRIBUTION PER FRAME
        RES->dat.H[11][i]=1;
        RES->dat.H[1][at->data.resnr]+=1;
        if(!(z==5)){ 
	  //IS NOT HYDROGEN
	  RES->dat.H[2][i]+=1;
	  if(rpres<=0.3)
	    RES->dat.H[9][i]+=1;
        }
      }
    }
    if( rpres <= r & !(at->data.edge) & y) {
      for(k=0;k<=j;k++){ 
	//is it a new atom?
	if(i==((RES->dat.H[z][k]))){ //NO
	  at->data.qmmm=1;
	}
      }
      if(!(at->data.qmmm)){
	RES->dat.H[z][j]=i;
	at -> data.edge=1;
	j++;
	switch(z){
	case 6:
	  //DONT CALCULATE C COORDINATION IN HISTOGRAMS
	  c++;
	  break;
	case 7:
	  n++;
	  break;
	case 8:
	  o++;
	  break;
	case 5:
	  h++;
	  break;
	case 10:
	  s++;
	  break;
	}
      }
      m++;
    }
  } 
  CAL->dat.B[0]=j;
  CAL->dat.B[1]=m;
  CAL->dat.B[5]=h;
  CAL->dat.B[6]=c;
  CAL->dat.B[7]=n;
  CAL->dat.B[8]=o;
  CAL->dat.B[10]=s;
}
