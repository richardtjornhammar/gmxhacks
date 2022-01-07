#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "hlist.h"
#include "define.h"

void sort (float *AR, int nrEntries) {
  int sa,sb;
  float ft,ft2;

  for(sa = 1; sa < nrEntries; ++sa)
    for(sb = nrEntries-1; sb >= sa; --sb) {
      // compare adjacent elements 
      if(AR[ (sb - 1)*2] > AR[ sb*2 ]) {
	// exchange elements 
	ft = AR[ (sb - 1)*2 ];
	ft2= AR[ (sb - 1)*2+1];
	AR[ (sb - 1)*2 ] = AR[ sb*2 ];
	AR[ (sb - 1)*2 +1 ] = AR[ sb*2 + 1 ];
	AR[ sb*2 ] = ft;
	AR[ sb*2 + 1 ] = ft2;
      }
    }
}

void sort_sq (float *AR, int nrEntries) {
  int sa,sb;
  float ft,ft2;

  for(sa = 1; sa < nrEntries; ++sa)
    for(sb = nrEntries-1; sb >= sa; --sb) {
      // compare adjacent elements 
      if(SQ(AR[ (sb - 1)*2]) > SQ(AR[ sb*2 ])) {
	// exchange elements 
	ft = AR[ (sb - 1)*2 ];
	ft2= AR[ (sb - 1)*2+1];
	AR[ (sb - 1)*2 ] = AR[ sb*2 ];
	AR[ (sb - 1)*2 +1 ] = AR[ sb*2 + 1 ];
	AR[ sb*2 ] = ft;
	AR[ sb*2 + 1 ] = ft2;
      }
    }
}

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

float dotP(float A[3], float B[3]){
  int i;
  float P=0.0;

  for(i=0;i<3;i++)
    P+=A[i]*B[i];
  
  return P;
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

//at,rpr,dr,Md,CELL,I,tframe
void calc_dipole_rad(LINK at,LINK cat, pvec mdens, float CELL[6], int I, int tframe)
{
  int i,j;
  float e[3],E[3],B[3],ff;
  float R[3],rlen=0.0,rtem=0.0,elen=0.0;
  LINK P,H;

  P=at->next->next;
  H=P->next; //printAtom(H);
  
  for(j=0;j<3;j++) //SMALLEST DISTANCE BETWEEN ATOMS INCLUDING MIRROR IMAGES SO MUST USE CUBIC BOX
    {
      R[j]=at->data.pos[j]-cat->data.pos[j];
      if(R[j]>CELL[j]*0.5)
	R[j]-=CELL[j];
      else
	if(abs(R[j])>CELL[j]*0.5)
	  R[j]+=CELL[j];
      rtem+=pow(R[j],2.0);
    }
  rlen=sqrt(rtem);rtem=0.0;
  
  for(j=0;j<3;j++) //SMALLEST DISTANCE BETWEEN ATOMS INCLUDING MIRROR IMAGES SO MUST USE CUBIC BOX
    {
      E[j]=H->data.pos[j]-at->data.pos[j]; //DIPOLE DIRECTION
      if(E[j]>CELL[j]*0.5)
	E[j]-=CELL[j];
      else
	if(abs(E[j])>CELL[j]*0.5)
	  E[j]+=CELL[j];
      rtem+=pow(E[j],2.0);
    }
  elen=sqrt(rtem);rtem=0.0;

  ff = dotP(R,E)/rlen/elen;
  mdens->dat[I] += ff; 
}

void calc_gID(pvec dens, float r, float dr, float rho,pvec gID, int verb){
  gID->dat[gID->i]=1.0/pi4/rho/r/r/dr*dens->dat[dens->i];
}

void calc_cID(pvec cD, pvec Md, pvec Nd, int i)
{
  if(Md->dat[i] == 0.0)
    cD->dat[i]=0.0;
  else
    cD->dat[i]=(Md->dat[i]/Nd->dat[i]);
}

void calc_Polarisation(pvec P,pvec cID,pvec gID, float rho,int verb)
{
  float mu=0.0473;
  int i;

  i=P->i;
  P->dat[i]=gID->dat[i]*cID->dat[i]*mu*rho;
}

void calc_simpPolarisation(LINK at, LINK cat, float CELL[6], pvec P, float r, float dr, int I, int do_charge)
{
  //N contains the amount of tip4p waters in this shell
  //M contains the orientation of tip4p waters in this shell projected on radius
  int i,j;
  float e[3],E[3],B[3],ff,p[3];
  float R[3],rlen=0.0,rtem=0.0,elen=0.0;
  LINK H1,H2,H;

  switch(do_charge){ 
  case 2://MOLECULE AND SOLVENT CHARGE HAS BEEN SPECIFIED; WE ARE WORKING WITH WATER
    H=at;
    for(j=0;j<3;j++)
      p[j]=H->data.mass_ch[0]*(H->data.pos[j]);
    H=H->next;
    if(!(H->data.name[0]=='O')){
      for(j=0;j<3;j++)
	p[j]+=H->data.mass_ch[0]*(H->data.pos[j]);
      
    }
    H=H->next;
    if(!(H->data.name[0]=='O')){
      for(j=0;j<3;j++)
	p[j]+=H->data.mass_ch[0]*(H->data.pos[j]);
      
    }
    H=H->next;
    if(!(H->data.name[0]=='O')){
      for(j=0;j<3;j++)
	p[j]+=H->data.mass_ch[0]*(H->data.pos[j]);
      
    }
    break;
  default://ASSUME TIP4P
    H1=at->next; 
    for(j=0;j<3;j++)
      p[j]=0.52*(H1->data.pos[j]);
    H2=H1->next;
    for(j=0;j<3;j++)
      p[j]+=0.52*(H2->data.pos[j]);
    H=H2->next; //printAtom(H);
    for(j=0;j<3;j++)
      p[j]-=1.04*(H->data.pos[j]);
    break;
  }

  for(j=0;j<3;j++) //SMALLEST DISTANCE BETWEEN ATOMS INCLUDING MIRROR IMAGES SO MUST USE CUBIC BOX
    {
      R[j]=at->data.pos[j]-cat->data.pos[j];
      if(R[j]>CELL[j]*0.5)
	R[j]-=CELL[j];
      else
	if(abs(R[j])>CELL[j]*0.5)
	  R[j]+=CELL[j];
      rtem+=pow(R[j],2.0);
    }
  rlen=sqrt(rtem);rtem=0.0;
  
  for(j=0;j<3;j++) //SMALLEST DISTANCE BETWEEN ATOMS INCLUDING MIRROR IMAGES SO MUST USE CUBIC BOX
    {
      E[j]=H->data.pos[j]-at->data.pos[j]; //DIPOLE DIRECTION
      if(E[j]>CELL[j]*0.5)
	E[j]-=CELL[j];
      else
	if(abs(E[j])>CELL[j]*0.5)
	  E[j]+=CELL[j];
      rtem+=pow(E[j],2.0);
    }
  elen=sqrt(rtem);rtem=0.0;

  ff = dotP(p,R)/rlen/pi4/rlen/rlen/dr;
  P->dat[I] += ff; 

}

void calc_polarisation_xyz(LINK at, LINK cat, float CELL[6], pvec Pp, pvec Px, pvec Py, pvec Pz, float r, float dr, int I, int do_charge, float nhat[3])
{
  //N contains the amount of tip4p waters in this shell
  //M contains the orientation of tip4p waters in this shell projected on radius
  int i,j;
  float e[3],E[3],B[3],ff,p[3];
  float R[3],rlen=0.0,rtem=0.0,elen=0.0;
  LINK H1,H2,H;

  switch(do_charge){ 
  case 2://MOLECULE AND SOLVENT CHARGE HAS BEEN SPECIFIED; WE ARE WORKING WITH WATER
    H=at;
    for(j=0;j<3;j++)
      p[j]=H->data.mass_ch[0]*(H->data.pos[j]);
    H=H->next;
    if(!(H->data.name[0]=='O')){
      for(j=0;j<3;j++)
	p[j]+=H->data.mass_ch[0]*(H->data.pos[j]);
    }
    H=H->next;
    if(!(H->data.name[0]=='O')){
      for(j=0;j<3;j++)
	p[j]+=H->data.mass_ch[0]*(H->data.pos[j]);
      
    }
    H=H->next;
    if(!(H->data.name[0]=='O')){
      for(j=0;j<3;j++)
	p[j]+=H->data.mass_ch[0]*(H->data.pos[j]);
      
    }
    break;
  default://ASSUME TIP4P
    H1=at->next; 
    for(j=0;j<3;j++)
      p[j]=0.52*(H1->data.pos[j]);
    H2=H1->next;
    for(j=0;j<3;j++)
      p[j]+=0.52*(H2->data.pos[j]);
    H=H2->next; //printAtom(H);
    for(j=0;j<3;j++)
      p[j]-=1.04*(H->data.pos[j]);
    break;
  }

  for(j=0;j<3;j++) //SMALLEST DISTANCE BETWEEN ATOMS INCLUDING MIRROR IMAGES SO MUST USE CUBIC BOX
    {
      R[j]=at->data.pos[j]-cat->data.pos[j];
      if(R[j]>CELL[j]*0.5)
	R[j]-=CELL[j];
      else
	if(abs(R[j])>CELL[j]*0.5)
	  R[j]+=CELL[j];
      R[j]*=nhat[j];       //PROJECT ONTO NHAT
      rtem+=pow(R[j],2.0);
    }
  rlen=sqrt(rtem);rtem=0.0;
  
  ff = (dotP(p,nhat))/dr;
  Pp->dat[I] += ff; 
  Px->dat[I] += p[0]/dr;
  Py->dat[I] += p[1]/dr;  
  Pz->dat[I] += p[2]/dr;
}

void calc_ndens(LINK solvh[], LINK cat,float r, float dr, pvec Ndens, pvec xi, int slen, float CELL[6],pvec Mdens, int verb)
{
  //THIS CURRENTLY DOESNT DISCERN ATOM TYPE
  int i,I,j,sa;
  LINK at;
  float S[3],rprs,rpr,B[6],DR[3],rtmp;
  int   T[3];
  
  I=Ndens->i;

  for(i=1;i<=slen;i++){

    at=solvh[i-1];
    rprs=0.0;

    for(j=0;j<3;j++) //SMALLEST DISTANCE BETWEEN ATOMS INCLUDING MIRROR IMAGES SO MUST USE CUBIC BOX
      {
	B[j]=at->data.pos[j]-cat->data.pos[j];
	if(B[j]>CELL[j]*0.5)
	  B[j]-=CELL[j];
	else
	  if(abs(B[j])>CELL[j]*0.5)
	    B[j]+=CELL[j];
	rprs+=pow(B[j],2.0);
      }
    rpr=sqrt(rprs);
 
    if(at->data.name[0]=='O' & at->data.name[1]=='W'){ //PICKS TIP4P OXYGEN
      rtmp=(r-rpr)>0?2.0*(r-rpr):-2.0*(r-rpr);
      if( (r-dr*0.5 < rpr) & (rpr <= r+dr*0.5) ){ //FOUND A FITTING ATOM fprintf(stdout,"* %f %f %f %f %s||",r,rpr,rtmp,dr,at->data.name);
	Ndens->dat[I]+=1.0;
	//	calc_dipole(at,cat,r,dr,Mdens,CELL,I);
      }
      xi->dat[I]+=(float)( rpr <= r );
    }
  }
}
float length(float R[3])
{
  float result=0.0;
  int i;
    
  for(i=0;i<3;i++)
    result+=SQ(R[i]);
  return sqrt(result);
}

void vcopy(float RC[3],float R[3])
{
  int i;
    
  for(i=0;i<3;i++)
    RC[i]=R[i];
}

void calc_atomPOL(LINK atomh[] ,int lna, LINK cat,float r0, float dr, float CELL[6], float nhat[3], pvec Pa, pvec Ef,float RESULTS[RESX][RESY], int tframe, int verb, int xyz, int dSIZE, int do_charge)
{
  int i,j,k,cgnr0=1,cgnr_o=0,I=0,ix;
  double ID;
  LINK at,tmp;
  float S[3],ff,rprc,rprf,rpc,rpf,rp,B[6],DR[3],rtmp=0.0,p[3],chrg=0.0;
  double rest;
  int   T[3];

  //THIS ROUTINE SHOULD CACLULATE A POLARISATION CONTRIBUTION FROM THE ATOMLIST
  //PROJECTED IN THE RADIAL DIRECTION FROM THE CENTRAL ATOM

  if( (!xyz) & ( do_charge > 1 )){
    for(i=1;i<=lna;i++){
      at=atomh[i-1];
      if(at==NULL)
	fatal("ATOM CANNOT BE NULL");
      if(!(at->data.cgnr==cgnr_o)){ //NEW CHARGE GROUP + ASSIGN OLD
	if(i>1){
     	  ff = dotP(p,DR)/rpc/pi4/rpc/rpc/dr;
	  I=(int)floor((rpc-r0)/dr);
	  if(I>=dSIZE){
	    fatal("INDEX EXCEEDED");
	  }
	  if(SQ(chrg)==0.0){//ONLY ASSIGN IF CHARGEGROUP IS ZERO
	    Pa->dat[I]+=ff;
	  }
	}
	p[0]=0.0;  p[1]=0.0; p[2]=0.0; rpc=0.0;  rpf=100.0;
	DR[0]=0.0;DR[1]=0.0;DR[2]=0.0; chrg=0.0; k=0;
      }
 
      rprc=0.0;
      for(j=0;j<3;j++) 
	{
	  //SMALLEST DISTANCE BETWEEN ATOMS INCLUDING 
	  //MIRROR IMAGES SO MUST USE CUBIC BOX
	  B[j]=at->data.pos[j]-cat->data.pos[j];
	  if(B[j]>CELL[j]*0.5)
	    B[j]-=CELL[j];
	  else
	    if(abs(B[j])>CELL[j]*0.5)
	      B[j]+=CELL[j];
	  //rprc+=pow(B[j],2.0);
	  p[j]+=at->data.mass_ch[0]*(at->data.pos[j]);
	  rprc+=SQ(B[j]);
	}
      rp=sqrt(rprc);
      I=(int)floor((rp-r0)/dr);
      Ef -> dat[I] += 1/pi4*at->data.mass_ch[0]/rp/rp;
      if(rp<rpf){
	vcopy(DR,B);
	rpc=rp;rpf=rp;
      }
      chrg+=at->data.mass_ch[0];
      //rpc=sqrt(rprc); //SHOULD TAKE CLOSEST 
      //AND FARTHEST DISTANCE IN ONE GROUP
      cgnr_o=at->data.cgnr;
      k++;
    }
  }
  fprintf(stdout,"COLLECTING RESULTS\n");
  for(ix=0;ix<dSIZE;ix++){
    RESULTS[7][ix]+=(Pa->dat[ix])/((float)tframe);
    RESULTS[8][ix]+=(Ef -> dat[ix])/((float)tframe);
  }
}

void calc_densDIR(LINK solvh[] ,int lns, LINK CAT,float r0, float dr, float CELL[6], float nhat[3], pvec Nd, pvec Md, pvec xi, pvec cPD, pvec gPD, pvec PHK, pvec Ps,pvec Ef, float RESULTS[RESX][RESY], int tframe, int verb, int xyz, int dSIZE, int do_charge)
{
  //THIS CURRENTLY DOESNT DISCERN ATOM TYPE
  int i,I,j,sa,ix;
  double ID;
  LINK at;
  float S[3],rprs,rpr,B[6],DR[3],rtmp=0.0;
  double rest;
  int   T[3];

  switch(xyz){
  case 0:{        // radial search
    for(i=1;i<=lns;i++){
      at=solvh[i-1];
      rprs=0.0;
      for(j=0;j<3;j++) //SMALLEST DISTANCE BETWEEN ATOMS INCLUDING MIRROR IMAGES SO MUST USE CUBIC BOX
	{
	  B[j]=at->data.pos[j]-CAT->data.pos[j];
	  if(B[j]>CELL[j]*0.5)
	    B[j]-=CELL[j];
	  else
	    if(abs(B[j])>CELL[j]*0.5)
	      B[j]+=CELL[j];
	  rprs+=pow(B[j],2.0);
	}
      rpr=sqrt(rprs);
      if(do_charge>0){
	I=(int)floor((rpr-r0)/dr);
	Ef -> dat[I] += 1/pi4*at->data.mass_ch[0]/rpr/rpr;
      }
      if(at->data.name[0]=='O' & at->data.name[1]=='W'){ //PICKS TIP4P OXYGEN
	I=(int)floor((rpr-r0)/dr);
	Nd->dat[I]+=1.0;
	calc_dipole_rad(at,CAT,Md,CELL,I,tframe);
	calc_simpPolarisation(at,CAT,CELL,Ps,r0,dr,I,do_charge);
	gPD->i=I;cPD->i=I;Nd->i=I;PHK->i=I;
	//if(r0+I*dr<){
	  calc_cID(cPD,Md,Nd,I);
	  calc_gID(Nd,(r0+(I-1)*dr),dr,RHOw,gPD,verb);
	  calc_Polarisation(PHK,cPD,gPD,(float)RHOw,verb);
	  //}
      }
    }
  }
    break;
  case 1:{ //PRODUCE P(x,y,z) IN NHAT DIRECTION
    for(i=1;i<=lns;i++){
      at=solvh[i-1];
      rprs=0.0;
      for(j=0;j<3;j++) //SMALLEST DISTANCE BETWEEN ATOMS INCLUDING MIRROR IMAGES SO MUST USE CUBIC BOX
	{
	  B[j]=(at->data.pos[j]-CAT->data.pos[j]);
	  if(B[j]>CELL[j]*0.5)
	    B[j]-=CELL[j];
	  else
	    if(abs(B[j])>CELL[j]*0.5)
	      B[j]+=CELL[j];
	  B[j]*=nhat[j];       //PROJECT ONTO NHAT
	  rprs+=pow(B[j],2.0); //CALCULATE DISTANCE OF PROJECTION
	}
      rpr=sqrt(rprs);
      
      if(at->data.name[0]=='O' & at->data.name[1]=='W'){ //PICKS TIP4P OXYGEN
	I=(int)floor((rpr-r0)/dr);
	gPD->i=I;cPD->i=I;Nd->i=I;PHK->i=I;
	Nd->dat[I]+=1.0;
	//calc_dipole_xyz(at,CAT,Md,CELL,I,tframe,nhat);
	calc_polarisation_xyz(at,CAT,CELL,Md,cPD,PHK,Ps,r0,dr,I,do_charge,nhat);
	//NOW PUT POLARISATION X IN CPD, Y IN PHK AND Z IN PS AND NHAT PROJ IN Md
      }
    } 
  }
    break;
  default:break;
  }
  
  for(ix=0;ix<dSIZE;ix++){
    rtmp+=Nd->dat[ix];
    rpr=r0+ix*dr;
    RESULTS[0][ix]=rpr;
    RESULTS[1][ix]+=Nd->dat[ix]/((float)tframe);
    RESULTS[2][ix]+=rtmp/((float)tframe);
    RESULTS[3][ix]+=Md->dat[ix]/((float)tframe);
    RESULTS[4][ix]+=cPD->dat[ix]/((float)tframe);
    RESULTS[5][ix]+=PHK->dat[ix]/((float)tframe);
    RESULTS[6][ix]+=(Ps->dat[ix])/((float)tframe);
    //RESULTS[8][ix]+=(Ef -> dat[ix])/((float)tframe);
  }
}


void corner_search(LINK SYS[],float CELL[6])
{
  int i,j,k;
  int len;
  float x=100.0,y=100.0,z=100.0;
  LINK at,ta;

  len=listlength(SYS[0]);

  //FIND MINIMAL COORDINATES
  for(i=1;i<=len;i++){
    at=SYS[i-1];
    x=(at->data.pos[0]<x?at->data.pos[0]:x);
    y=(at->data.pos[1]<y?at->data.pos[1]:y);
    z=(at->data.pos[2]<z?at->data.pos[2]:z);
  }
  CELL[3]=x;
  CELL[4]=y;
  CELL[5]=z;
}


void initPvec(pvec array, int N)
{
  int i;

  for(i=0;i<N;i++){
    //fprintf(stdout,"\n%d\n",i);
    array->dat[i]=0.0;
  }
}

void update_topology(LINK atom,LINK solv,LINK cat,FILE *in,int typ, float CELL[6],float diam[2],int tframe)
{
  if(!(typ==etGRO))
    fatal("NOT A SUPPORTED TYPE");
  else //PUT IN NEW COORDINATES
    gro_update(atom,solv,cat,in,CELL,diam,tframe);
}

void printResults(float R[RESX][RESY],FILE *ofile,int NR)
{
  int i;
  //P_w in 6 P_p in 7 E_f in 8
  for(i=0;i<NR;i++)
    fprintf(ofile,"\n%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f",R[0][i],R[1][i],R[2][i],R[3][i],R[4][i],R[5][i],R[8][i],R[7][i],R[6][i],R[7][i]+R[6][i],SGN(R[7][i]+R[6][i])*(R[7][i]+R[6][i]),(R[8][i]+R[7][i]+R[6][i])*R[8][i]);
  /*
    C12::D dot E
    C11::|Pol_tot|
    C10::Pol_tot
    C9 ::Pol_w
    C8 ::Pol_c
    C7 ::E
    C6 ::PHK
    C5 ::cID
    C4 ::M
    C3 ::Xi
    C2 ::N   
    C1 ::r 
  */
  fprintf(ofile,"\n");

}
