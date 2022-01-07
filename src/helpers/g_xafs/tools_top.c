#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "hlist.h"
#include "define.h"

void pdb_create(LINK atom, LINK solv, FILE *fpi_g, FILE *outs, float CELL[])
{
  int  i,j,k,l,m;
  int  resq,ike,nratoms,count;
  char line[SLEN],cont1[SLEN],cont2[SLEN],cont3[SLEN],cont4[SLEN],cont5[SLEN];
  char cont6[SLEN],cont7[SLEN],cont8[SLEN],cont9[SLEN],cont10[SLEN],cont11[SLEN];
  char fcryst[SLEN];
  LINK new;
  /* FOR TWO DIFFERENT PDB FORMATS GIVEN BELOW
  
 1      2    3   4   5   6      7       8       9       10   11     12 
ATOM     34 CL2A DPP     1       7.770  56.790  19.830  1.00  0.00
     X    X    X   X   X    X        X       X       X     X     X           X
     6    5    5   4   3    5        8       8       8     5     5          11 
ATOM    161  CG  LEU A  25      21.502 -16.911   9.118  1.00 49.63           C 

   */
  while(!feof(fpi_g)) {
    i=0;
    fgets(line,SLEN,fpi_g);
    if(!strncmp(line,"TITLE ",6)) {
    } //if
    if(!strncmp(line,"CRYST1",6)) {
      sscanf(line,"%6s%9f%9f%9f%9f%9f%9f",cont1,&CELL[0],&CELL[1],&CELL[2],&CELL[3],&CELL[4],&CELL[5]);
    } //if
    if(!strncmp(line,"ATOM  ",6)| !strncmp(line,"HETATM",6) ){// | !strncmp(line,"TER   ",6)) {
      i++;
      sscanf(line,"%6s%5s%5s%4s%3s%6s%8s%8s%8s%5s%5s%11s",cont1,cont2,cont3,cont4,cont5,cont6,cont7,cont8,cont9,cont10,cont11,fcryst);
      new = getAtom(atom,listlength(atom));
      strncpy(new -> data.name,cont3,DIME-1);
      new -> data.name[DIME-1]='\0';
      strncpy(new -> data.resi,cont4,DIME);
      new -> data.resi[DIME]='\0';
      new -> data.flag = 0;
      new -> data.c='X';
      if( isalpha(cont5[0])){
	new->data.resnr=atoi(cont6);
	new->data.pos[0]=atof(cont7);
	new->data.pos[1]=atof(cont8);
	new->data.pos[2]=atof(cont9);
	new->data.pos[3]=atof(cont10);
	new->data.pos[4]=atof(cont11);
      }
      else{
	new->data.resnr=atoi(cont5);
	new->data.pos[0]=atof(cont6);
	new->data.pos[1]=atof(cont7);
	new->data.pos[2]=atof(cont8);
	new->data.pos[3]=atof(cont9);
	new->data.pos[4]=atof(cont10);
      }
      new->data.ftyp[0]=cont3[0];
      new->data.flag=1;
      addAtom(atom);
    } //if
  } //while

}

void write_gro(LINK atom, LINK solv, FILE *outs, float CELL[], int typ)
{
  int lna,lns,i,j,k,l,m,nr=0,ornr=0;
  LINK at;
  lna=listlength(atom);
  lns=listlength(solv);

  fprintf(outs,"%s\n","TITLE     GENERATED STRUCTURE FILE");
  fprintf(outs,"%10d\n",lna+lns);
  for(i=1;i<=lna;i++){
    at=getAtom(atom,i);
    if(at==NULL)
      break;
    if(at->data.flag){
      l=((at->data.resnr-ornr)==0?0:1);
      if(l)
	nr++;
      fprintf(outs,"%5i%-5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",nr,at->data.resi,at->data.name,i,at->data.pos[0]*(typ==etGRO?1.0:0.1),at->data.pos[1]*(typ==etGRO?1.0:0.1),at->data.pos[2]*(typ==etGRO?1.0:0.1),at->data.pos[3],at->data.pos[4],at->data.pos[5]);
      ornr=at->data.resnr;
    }
  }
  j=i;
  if(lns>1)
    for(i=1;i<lns;i++){
      at=getAtom(solv,i);
      if(at==NULL)
	break;
      if(at->data.flag){
	l=((at->data.resnr-ornr)==0?0:1);
	if(l)
	  nr++;
        fprintf(outs,"%5i%-5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",nr,at->data.resi,at->data.name,i+j-1,at->data.pos[0]*(typ==etGRO?1.0:0.1),at->data.pos[1]*(typ==etGRO?1.0:0.1),at->data.pos[2]*(typ==etGRO?1.0:0.1),at->data.pos[3],at->data.pos[4],at->data.pos[5]);
	ornr=at->data.resnr;
      }
    }
    fprintf(outs,"     % 9.5f% 9.5f% 9.5f",CELL[0],CELL[1],CELL[2]);
}


void write_pdb_tot(LINK atom, FILE *outs, float CELL[], int typ, int NR)
{
  int lna,lns,i,j,k,l,m,nr=0,ornr=0;
  LINK at;

  lna=listlength(atom);
  if(NR>lna){
    fatal("COULDNT CREATE FINAL OUTPUT\n");
  }
	
  fprintf(outs,"%s\n","TITLE     GENERATED STRUCTURE FILE");
  fprintf(outs,"CRYST1%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f\n",CELL[0],CELL[1],CELL[2],CELL[3],CELL[4],CELL[5]);

  for(i=1;i<=NR;i++){
    at=getAtom(atom,i);
    if(at==NULL)
      break;
    if(at->data.flag){
      l=((at->data.resnr-ornr)==0?0:1);
      if(l)
	nr++;
      fprintf(outs,"ATOM  %5d %-4s%4s %c %5d  %8.3f%8.3f%8.3f%6.2f%6.2f           %c\n",i,at->data.name,at->data.resi,at->data.c,nr,at->data.pos[0]*(typ==etGRO?10.0:1.0),at->data.pos[1]*(typ==etGRO?10.0:1.0),at->data.pos[2]*(typ==etGRO?10.0:1.0),at->data.pos[3],at->data.pos[4],at->data.ftyp[0]);
      ornr=at->data.resnr;
    }
  }
  
  fprintf(outs,"END\n");
}

void write_pdb(LINK atom, LINK solv, FILE *fpi_g, FILE *outs, float CELL[], int typ)
{
  int lna,lns,i,j,k,l,m,nr=0,ornr=0;
  LINK at;

  lna=listlength(atom);
  lns=listlength(solv);

  fprintf(outs,"%s\n","TITLE     GENERATED STRUCTURE FILE");
  fprintf(outs,"CRYST1%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f\n",CELL[0],CELL[1],CELL[2],CELL[3],CELL[4],CELL[5]);

  for(i=1;i<=lna;i++){
    at=getAtom(atom,i);
    if(at==NULL)
      break;
    if(at->data.flag){
      l=((at->data.resnr-ornr)==0?0:1);
      if(l)
	nr++;
      fprintf(outs,"ATOM  %5d %-4s%4s %c %5d  %8.3f%8.3f%8.3f%6.2f%6.2f           %c\n",i,at->data.name,at->data.resi,at->data.c,nr,at->data.pos[0]*(typ==etGRO?10.0:1.0),at->data.pos[1]*(typ==etGRO?10.0:1.0),at->data.pos[2]*(typ==etGRO?10.0:1.0),at->data.pos[3],at->data.pos[4],at->data.ftyp[0]);
      ornr=at->data.resnr;
    }
  }
  j=i;
  if(lns>1)
    for(i=1;i<lns;i++){
      at=getAtom(solv,i);
      if(at==NULL)
	break;
      if(at->data.flag){
	l=((at->data.resnr-ornr)==0?0:1);
	if(l)
	  nr++;
	fprintf(outs,"ATOM  %5d %-4s%4s %c %5d  %8.3f%8.3f%8.3f%6.2f%6.2f           %c\n",i+j-1,at->data.name,at->data.resi,at->data.c,nr,at->data.pos[0]*(typ==etGRO?10.0:1.0),at->data.pos[1]*(typ==etGRO?10.0:1.0),at->data.pos[2]*(typ==etGRO?10.0:1.0),at->data.pos[3],at->data.pos[4],at->data.ftyp[0]);
	ornr=at->data.resnr;
      }
    }
  fprintf(outs,"END\n");
}

void gro_create(LINK atom, LINK solv, FILE *fpi_g, float CELL[])
{
  int   j,q,resq,ike,m,nratoms,a=0;
  char  line[SLEN],cont1[SLEN],cont2[SLEN],cont3[SLEN],cont4[SLEN],cont5[SLEN];
  char  cont6[SLEN],cont7[SLEN],cont8[SLEN],cont9[SLEN],cont10[SLEN];
  char  title[SLEN],fcryst[SLEN],blank[SLEN],dump1[10],dump2[5],alphab[26];
  char* chrptr=NULL;
  char* chrptr2=NULL;
  int   trasher=0,ratsher=0,lna=0,lns=0;
  LINK  new;

  strcpy(alphab,"ABCDEFGHIJKLMNOPQRSTUVWXYZ");
  fgets(title,SLEN,fpi_g);
  if(!strlen(title)) {
    fatal("GRO CREATE");
  }
  fgets(line,SLEN,fpi_g); sscanf(line,"%s",cont1);
  nratoms=atoi(cont1);
  j=0; q=0; resq=0; ike=1;
  while(j<nratoms) {
    fgets(line,SLEN,fpi_g);
    sscanf(line,"%s %s %s %s %s %s",cont1,cont2,cont3,cont4,cont5,cont6);    
    m=!(int)fabs((long unsigned int)strstr(cont1,"SOL"));
    if(strlen(cont2)>5 && j<nratoms){
      strcpy(cont6,cont5);strcpy(cont5,cont4);strcpy(cont4,cont3);
      chrptr2=strpbrk(cont2,"1234567890");
      if(strlen(chrptr2)>5)
	chrptr2++;
      strcpy(cont7,chrptr2);
      chrptr2[0]='\0';
      strcpy(cont8,cont2);
    }
    else
      strcpy(cont8,cont2);
    if(m) {
      if(q==0){
	new=atom;q=1;
      }
      else
	new=new->next;
      lna=j+1;
      trasher=atoi(cont1);
      chrptr=strpbrk(cont1,alphab);
      new->data.resnr=trasher;
      strncpy(new -> data.resi,chrptr,DIME);
      strncpy(new -> data.name,cont8, DIME-1); 
      new -> data.name[DIME-1]='\0';
      new -> data.resi[DIME]='\0';
      new -> data.edge=0;new->data.qmmm=0;
      new -> data.pos[0]=atof(cont4);//new->data.pos[0]*=10;
      new -> data.pos[1]=atof(cont5);//new->data.pos[1]*=10;
      new -> data.pos[2]=atof(cont6);//new->data.pos[2]*=10;
      new -> data.ftyp[0]=new->data.name[0];
      new -> data.ftyp[1]='\0';
      new -> data.flag=1;
      new -> data.c=' ';
      addAtom(atom);
    }
    else {
      if(resq==0){
	new=solv;resq=1;
      }
      else
	new=new->next;
      trasher=atoi(cont1);
      chrptr=strpbrk(cont1,alphab);
      new->data.resnr=trasher;
      strncpy(new -> data.resi,chrptr,DIME);
      strncpy(new -> data.name,cont8, DIME-1); 
      new -> data.resi[DIME]='\0';
      new -> data.name[DIME-1]='\0';
      //strncpy(new -> data.resi,"SOL",DIME);
      new->data.edge=0;new->data.qmmm=0;
      new->data.pos[0]=atof(cont4);//new->data.pos[0]*=10;
      new->data.pos[1]=atof(cont5);//new->data.pos[1]*=10;
      new->data.pos[2]=atof(cont6);//new->data.pos[2]*=10;
      new->data.ftyp[0]=new->data.name[0];
      new->data.ftyp[1]='\0';
      new->data.flag=1;
      new->data.c=' ';
      addAtom(solv);
    }
    strcpy(cont7,blank);
    strcpy(cont8,blank);
    j++;
  }
  fgets(fcryst,SLEN,fpi_g);
  sscanf(fcryst,"     %s %s %s",cont1,cont2,cont3);
  CELL[0]=atof(cont1);//*10;
  CELL[1]=atof(cont2);//*10;
  CELL[2]=atof(cont3);//*10;
  delAtom(atom,listlength(atom));
  delAtom(solv,listlength(solv));
}

void trans_top(LINK atom, LINK solv, float T[3], int lna, int lns)
{
  int i;
  LINK new;

  new=atom;
  for(i=0;i<lna;i++){
    new->data.pos[0]+=T[0];
    new->data.pos[1]+=T[1];
    new->data.pos[2]+=T[2];
    new=new->next;
  }
  new=solv;
  for(i=0;i<lns;i++){
    new->data.pos[0]+=T[0];
    new->data.pos[1]+=T[1];
    new->data.pos[2]+=T[2];
    new=new->next;
    if(new==NULL)
      break;
  }
}


void merge_top(LINK atom, LINK solv, LINK atomi, LINK solvi, float RP[3])
{
  //INSERTS IT TRANSLATED BY RP
  float CM[3],R[3],X[3];
  int   i,j,k,l,m,n;
  int   lna,lns,lia,lis;
  LINK new;

  new=atomi;

  lns=listlength(solv);
  lna=listlength(atom);

  lia=listlength(atomi);
  lis=listlength(solvi); 
  
  trans_top(atomi,solvi,RP,lia,lis);

  //append all data to the list
  new=getAtom(atom,lna-1);
  new->next=atomi;
  atomi->prev=new;

  new=getAtom(solv,lns);
  new->next=solvi;
  solvi->prev=new;
}

void assign_topology(LINK atom, LINK solv,LINK itpc, FILE *intop, int do_charge)
{
  //EXAMPLE LINE FROM TOPOLOGY FILE::
  //     1       Zn2+      1     ZN     ZN      1          2      65.37   ; qtot 2.0
  int i=0,j,k;
  int N[2];
  char line[SLEN],cont1[SLEN],cont2[SLEN],cont3[SLEN],cont4[SLEN],cont5[SLEN];
  char cont6[SLEN],cont7[SLEN],cont8[SLEN],cont9[SLEN],cont10[SLEN];
  LINK tmp,tmpi;

  N[0]=listlength(atom);
  N[1]=listlength(solv);
  tmp=atom;

  while(!feof(intop)) {
    fgets(line,SLEN,intop);
    if(!strncmp(line,"[ atoms ]",9)) {
      while(i<N[0]){
	fgets(line,SLEN,intop);
	if(line[0]==';')
	  fgets(line,SLEN,intop);
	sscanf(line,"%6s%11s%7s%7s%7s%7s%11s%12s",cont1,cont2,cont3,cont4,cont5,cont6,cont7,cont8);
	
	if(!strcmp(tmp->data.resi,cont4) & !strcmp(tmp->data.name,cont5)){
	  //printf("assigning\n\n");
	  //printAtom(tmp);
	  tmp->data.cgnr=atoi(cont6);
	  tmp->data.mass_ch[0]=atof(cont7); //chrg
	  tmp->data.mass_ch[1]=atof(cont8); //mass
	  strcpy(tmp->data.ftyp,cont2);
	  //printAtom(tmp);
	}
	else{
	  printAtom(tmp);
	  fatal("WRONG TOPOLOGY FILE");
	}
	tmp=tmp->next;
	i++;
      }
    }
  }//THIS MEANS THAT ALL NON SOL ATOMS MUST HAVE BEEN INITIALIZED
  
  if(do_charge>1){
    //SOL LOOP
    tmp=solv;tmpi=itpc;
    for(i=0;i<N[1];i++){
      strcpy(tmp->data.ftyp,tmpi->data.ftyp);
      tmp->data.mass_ch[0]=tmpi->data.mass_ch[0];
      tmp=tmp->next;
      tmpi=tmpi->next;
    }
  }
}

LINK centerAtom(LINK atom, char *label)
{
  int i=0,j,k,N;
  LINK this;

  this=atom;
  while(!(this==NULL)){
    if(( toupper(label[0]) == toupper(this->data.name[0]) & toupper(label[1]) == toupper(this->data.name[1]) ) ) //!strcmp(label,this->data.name)
      break;
    this=this->next;
    i++;
  }
  this->data.edge=i;
  return this;
}

void init_itp(LINK itpc, FILE *initp){
  int   i,j,k,clear=0;
  char  line[SLEN],cont1[SLEN],cont2[SLEN],cont3[SLEN],cont4[SLEN],cont5[SLEN];
  char  cont6[SLEN],cont7[SLEN],cont8[SLEN],cont9[SLEN],cont10[SLEN];
  LINK  tmp;

  //EXAMPLE::     1	 opls_135  1	SCH      C      1      -0.25
  //;   nr   type  resnr residue  atom   cgnr     charge       mass
  //THIS ROUTINE IS USED AS THE CYCLIC DEFINITION OF THE SOL
  j=0;
  while(!feof(initp)) {
    fgets(line,SLEN,initp);
    if(!strncmp(line,"[ atoms ]",9)) {
      fgets(line,SLEN,initp);
      while(!clear)
	if(line[0]==';' | line[0]=='#'){
	  sscanf(line,"%s %s",cont1,cont2);
	  if(!strncmp(cont1,"#endif",6))
	    break;
	  else
	    fgets(line,SLEN,initp);
	}
	else
	  clear=1;
      while(!(line[0]=='#')){
	sscanf(line,"%6s%11s%3s%7s%7s%7s%11s%12s",cont1,cont2,cont3,cont4,cont5,cont6,cont7,cont8);
	tmp = getAtom(itpc,listlength(itpc));
	tmp->data.mass_ch[0]=atof(cont7); //chrg
	tmp->data.mass_ch[1]=atof(cont8); //mass
	strcpy(tmp->data.ftyp,cont2);
	addAtom(itpc);
	j++;
	//fprintf(stdout,"\n%d\n",j);
	fgets(line,SLEN,initp);
      }
    }
  }
  //delAtom(itpc,listlength(itpc));
  tmp->next=itpc;
}

void gro_update(LINK atom, LINK solv, LINK cat, FILE *fpi_g, float CELL[], float DIA[],int NFRAMES)
{
  //update atoms from gro frame
  int j,q,resq,ike,m,nratoms;
  char line[SLEN],cont1[SLEN],cont2[SLEN],cont3[SLEN],cont4[SLEN],cont5[SLEN];
  char cont6[SLEN],cont7[SLEN],cont8[SLEN],cont9[SLEN],cont10[SLEN];
  char fcryst[SLEN],title[SLEN],blank[SLEN];
  char *chrptr2;
  LINK new;
  float x=100.0,y=100.0,z=100.0;
  float dia[2],dib[2];

  dia[0]=0.0;dia[1]=0.0;

  fgets(title,SLEN,fpi_g);
  if(!strlen(title)) { //the line should be "frame t= (float)"
     fprintf(stdout,"%s",title);
     fatal("GRO_UPDATE::DIDNT GET TITLE FROM GRO FILE\n");
  }
  if(strncmp(title,"frame",5)){
    fprintf(stdout,"GOT::%s",title);
    fatal("GRO_UPDATE::INDEX DRIFT");
  }
  fgets(line,SLEN,fpi_g); sscanf(line,"%s",cont1);
  nratoms=atoi(cont1);

  //fprintf(stdout,"GETTING: NRATOMS = %d\n\n",nratoms);

  j=0; q=0; resq=0; ike=1;
  while(j<nratoms) {
    fgets(line,SLEN,fpi_g);
    sscanf(line,"%s %s %s %s %s %s",cont1,cont2,cont3,cont4,cont5,cont6);    
    m=!(int)fabs((long unsigned int)strstr(cont1,"SOL"));
    if(strlen(cont2)>5 && j<nratoms){
      strcpy(cont6,cont5);strcpy(cont5,cont4);strcpy(cont4,cont3);
      chrptr2=strpbrk(cont2,"1234567890");
      if(strlen(chrptr2)>5)
	chrptr2++;
      strcpy(cont7,chrptr2);
      chrptr2[0]='\0';
      strcpy(cont8,cont2);
    }
    else
      strcpy(cont8,cont2);
    if(m) {
      if(q==0){
	new=atom;q=1;
      }
      else
	new=new->next;
      new->data.pos[0]=atof(cont4);
      new->data.pos[1]=atof(cont5);
      new->data.pos[2]=atof(cont6);
      /*  x=(new->data.pos[0]>x)?x:new->data.pos[0];
      y=(new->data.pos[1]>y)?y:new->data.pos[1];
      z=(new->data.pos[2]>z)?z:new->data.pos[2];
      dib[0]=sqrt(SQ(new->data.pos[0])+SQ(new->data.pos[1])+SQ(new->data.pos[2]));
      dia[0]=dia[0]<dib[0]?dib[0]:dia[0];
      dib[1]=sqrt(SQ(new->data.pos[0]-cat->data.pos[0])+SQ(new->data.pos[1]-cat->data.pos[1])+SQ(new->data.pos[2]-cat->data.pos[2]));
      dia[1]=dia[1]<dib[1]?dib[1]:dia[1];*/
    }
    else {
      if(resq==0){
	new=solv;resq=1;
      }
      else
	new=new->next;
      if(strncmp(new -> data.name,cont8, DIME-1)){ //CHECKING FOR MATCH IN SOLVENT
	fprintf(stdout,"DRIFT %s / %s\n",new -> data.name,cont8);
	fprintf(stdout,"DRIFT INFO: %d resnr: %d\n",j,atoi(cont1));
	fatal("ATOMS NOT MATCHING");
      }
      new->data.pos[0]=atof(cont4);
      new->data.pos[1]=atof(cont5);
      new->data.pos[2]=atof(cont6);
      /* x=(new->data.pos[0]>x)?x:new->data.pos[0];
      y=(new->data.pos[1]>y)?y:new->data.pos[1];
      z=(new->data.pos[2]>z)?z:new->data.pos[2];*/
    }
    strcpy(cont7,blank);
    strcpy(cont8,blank);
    j++;
  }// while nratoms 
  fgets(fcryst,SLEN,fpi_g);
  sscanf(fcryst,"     %s %s %s",cont1,cont2,cont3);
  CELL[0]=atof(cont1); CELL[1]=atof(cont2); CELL[2]=atof(cont3);
  /*  DIA[0]+=dia[0]/((float)(NFRAMES-1));
  DIA[1]+=dia[1]/((float)(NFRAMES-1));
  if(cat->data.edge){
    cat->data.pos[0]=x;
    cat->data.pos[1]=y;
    cat->data.pos[2]=z;
  }
  */
}

