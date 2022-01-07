/*	
*	C CODE WRITTEN BY 
*       RICHARD TJÖRNHAMMAR 
*       LATEST DATE: 2009-11-05
*	COMPILE: gcc -lm *.c
*
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "list.h"

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

struct ZNAME{
  char *name[120];
}ex_zn={
  "LA","H ","He","Li","Be","B ","C ","N ","O ","F ","Ne","Na","Mg",
  "Al","Si","P ","S ","Cl","Ar","K ","Ca","Sc","Ti","V ","Cr","Mn",
  "Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr",
  "Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb",
  "Te","I ","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd",
  "Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W ","Re","Os","Ir",
  "Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
  "Pa","U ","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
  "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg"};

int main (int argc, char *argv[]) {

  LINK		atom,new,temp,at,solv;
  LINK		hashAtom[MAXL],hashSolv[MAXL];
  FILE		*fpi_g,*fpo,*fpi_h,*fpi_hs,*frho,*fpos;
  int		a=0,b=1,bb=1,nb=0,c,nc,fret=0,nl=0,x=0,nw=0,qml=0,la=0,rso=1,y=0;
  int           TYPES[400],POTEN[400],ZPOT[400],AMOUNT[400]; //LINES FOR XAFS DATA GENERATION
  int		sa,sb,st,count,interr,solrsnr=0,version=4,trail=0;
  int		NFRAMES=1,frame=1,ifr=0,startfr=1,N=1,z=1,xanes=0;
  float		dr,r,qtot,rpast=0.0,rpres=0.0,r0=0.0,rx=0;
  float		radius[DIM3],Rd[2000],Rs[10000];
  int		Cndx=0,o;
  float		CRYST[DIM3+2],QMcenter[DIM3],DISPR[DIM3],S[DIM3];
  int           T[DIM3];
  float         rho_tot[12][200];
  char          *chrptr=NULL;
  char*         chrptr2=NULL;
  char*		fname_i=NULL;
  char*		fnend3=NULL;
  char*		fnend2=NULL;
  char*		fnend1=NULL;
  char*		fname_o=NULL;
  char*		atcen=NULL;				
  char*		fname_j=NULL;
  char*		fname_k=NULL;
  char*		lcont;
  int		resSEQ,resq,steps,keep=1,ike=1;
  float         charge=0.0;
  char		content[SLEN],Z[2];
  char		cont1[SLEN],cont2[SLEN],cont3[SLEN],cont4[SLEN] ,cont5[SLEN];
  char		cont6[SLEN],cont7[SLEN],cont8[SLEN],cont9[SLEN],cont10[SLEN];
  char 		title[SLEN],file1[SLEN],file2[SLEN],file3[SLEN],cont11[SLEN];
  char		prefix[SLEN];
  char          fcryst[SLEN];
  char		alab[5],cnt2[3],line[SLEN],qmline[SLEN];
  int		i,j,k,l,m,n,q,nratoms,sol=0,trasher=0,ratsher=0,lna=0,lns=0;
  char          blank[40]="                                        ",bokstav;
  char          frmat[10],dump1[10],dump2[5],alphab[26];
  MAT           RES,RESs;
  PAR           CAL,CALs;
  DEN           RHO,RHOs;
  struct ZNAME  ztags = ex_zn;

	for(i=0;i<DIME+3;i++)
		CRYST[i]=0.0;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
		  switch (argv[i][1]) {
		    // SETS QM CENTER ATOM LABEL
		  case 'c':	atcen = argv[++i];
		    a = 1;
		    break;
		  case 'x':     xanes=atoi(argv[++i]);
		    break;
		    // COUNT ALL Z
		  case 'z':     z=1; b=1;
		    break;
		    // THE INPUT FILE
		  case 'i':	fname_i = argv[++i]; 
		    break;
		    // THE OUTPUT FILE
		  case 'o':	fname_o = argv[++i];
		    break;
		  case 'r':     r=atof(argv[++i]);
		                dr=atof(argv[++i]);
		                steps=atoi(argv[++i]);
			        x=1; 
		    break;
		    // USE UNIQUE MARKERS N,C,O,H
		  case 'a':     b=1;
		    break;
		    // DONT SEARCH SOLVENT BUT ATOMS "BACKBONE"
		  case 'b':     bb=1;
		    break;
		    // SEARCH FOR TOTAL COORDINATION WITH TYPE TRACKING
		  case 'h':     fret=1;
		    break;
		    // SETS MAX NUMBER OF FRAMES IN FILE
		  case 'f':	NFRAMES=atoi(argv[++i]);
		    break;
		    // STARTS FRAME READ AT THIS FRAME
		  case 's':     startfr=atoi(argv[++i]);
		    break;
		    // TELLS US WHICH CENTER IN ORDER TO USE
		  case 'k':     keep=atoi(argv[++i]);
		    break;
		  default:	fprintf(stderr,"Unknown switch %s\n", argv[i]);
		  }
		}
	}
	
	if (fname_i != NULL) { 
	  printf("Input name is     : %s \n", fname_i);
	  nw=strlen(fname_i);
	}
    	else {
	  printf("\nInput options:\n");
	  printf("\n-i input file name\n-o output prefix generated files:: prefix{_s _b _histb _hists}\n-c specifies center atom label (-c ZN)\n-f NFRAMES specifies the number of frames to count\n-h produces histogram files\n-x [num] produces (1)xanes/(2)xafs input form data\n\nEx.: ./g_coord -i coord.gro -o test -c ZN -f 1 -x 2\n");
		return -1;
    	}
    	if (fname_o != NULL) {
	  printf("Ouput prefix is   : %s \n", fname_o);
	}
	else {
	  printf("Output prefix needed \n");
	  return -1;
	}
	
	if(!x) {
	  r=0.15; r0=r;
	  steps=75;
	  dr=0.01;
	}
	else
	  r0=r;

	printf("\nOpening:  %s\n",fname_i);	
	fpi_g = fopen(fname_i,"r");
	strcpy(prefix,fname_o);
	fpo = fopen(strcat(prefix,"_b"),"w+");
	strcpy(prefix,fname_o);
	fpos = fopen(strcat(prefix,"_s"),"w+");
	frho = fopen("rho_bin","wb");
	printf("\nDone opening");
	
	atom=initAtom();
	solv=initAtom();
	
	RES = malloc(sizeof(CONTAINED));
	CAL = malloc(sizeof(PARTITION));
	RHO = malloc(sizeof(DENSITY));
	RESs = malloc(sizeof(CONTAINED));
	CALs = malloc(sizeof(PARTITION));
	RHOs = malloc(sizeof(DENSITY));

	printf(": atom init\n");
	j=0; l=0; m=0; k=0; resq=0; count=0;
	strcpy(alphab,"ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	for(i=0;i<10;i++)
	  RHO->dat.D[i]=0.0;
	i=0;

	CAL->dat.B[3]=bb;
	CAL->dat.B[4]=z;
	CAL->dat.B[5]=0;
	CAL->dat.B[6]=0;
	CAL->dat.B[7]=0;
	CAL->dat.B[8]=0;
	CAL->dat.B[10]=0;

	for(q=0;q<MAXL;q++){
	  RES->dat.H[6][q]=0;
	  RES->dat.H[7][q]=0;
	  RES->dat.H[8][q]=0;
	}

	printf("...getting data\n\n");
	
	// GET ATOMS AND POSITIONS IN GRO FORMAT
	frame=1;
	while(frame<=NFRAMES) {
	  fgets(title,SLEN,fpi_g);
	  if(!strlen(title)) {
	    printf("\nERROR!\n");
	  }
	  fgets(line,SLEN,fpi_g); sscanf(line,"%s",cont1);
	  if(frame==1){
	    nratoms=atoi(cont1);
	  }
	  j=0; q=0; resq=0; ike=1;
	  while(j<nratoms) {
	    /*	    RES->dat.H[0][j]=0;
	    RES->dat.H[1][j]=0;
	    RES->dat.H[2][j]=0;
	    RES->dat.H[9][j]=0;
	    RES->dat.H[10][j]=0;
	    RES->dat.H[11][j]=0;*/
	    fgets(line,SLEN,fpi_g);
	    sscanf(line,"%s %s %s %s %s %s",cont1,cont2,cont3,cont4,cont5,cont6);    
	    if(1) {
	      m=!(int)fabs((long unsigned int)strstr(cont1,"SOL"));
	      if(strlen(cont2)>5 && j<nratoms){
		strcpy(cont6,cont5);strcpy(cont5,cont4);strcpy(cont4,cont3);
		if(frame==1){
		  chrptr2=strpbrk(cont2,"1234567890");
		  if(strlen(chrptr2)>5)
		    chrptr2++;
		  strcpy(cont7,chrptr2);
		  chrptr2[0]='\0';
		  strcpy(cont8,cont2);
		}
	      }
	      else
		strcpy(cont8,cont2);
	      if(m) {
		if(q==0){
		  new=atom;q=1;
		}
		else
		  new=new->next;
		if(frame==1){
		  lna=j+1;
		  trasher=atoi(cont1);
		  chrptr=strpbrk(cont1,alphab);
		  new->data.resnr=trasher;
		  strncpy(new -> data.resi,chrptr,DIME-1);
		  strncpy(new -> data.name,cont8, DIME-1); 
		  new -> data.name[DIME-1]='\0';
		  new -> data.resi[DIME-1]='\0';
		  new->data.edge=0;new->data.qmmm=0;
		}
		new->data.pos[0]=atof(cont4);
		new->data.pos[1]=atof(cont5);
		new->data.pos[2]=atof(cont6);
		if(a){
		  if(!strcmp(atcen,new->data.name)){
		    if(ike==keep&&frame==1){
		      Cndx=listlength(atom);
		    }
		    if(ike==keep){
		      for(i=0;i<3;i++)
			QMcenter[i]=new->data.pos[i];
		    }
		    ike++;
		  }
		  
		}
		if(frame==1)
		  addAtom(atom);
	      }
	      else {
		if(resq==0){
		  new=solv;resq=1;
		}
		else
		  new=new->next;
		if(frame==1){
		  trasher=atoi(cont1);
		  chrptr=strpbrk(cont1,alphab);
		  new->data.resnr=trasher;
		  strncpy(new -> data.resi,chrptr,DIME-1);
		  strncpy(new -> data.name,cont8, DIME-1); 
		  new -> data.resi[DIME-1]='\0';
		  new -> data.name[DIME-1]='\0';
		  strncpy(new -> data.resi,dump2,DIME);
		  new->data.edge=0;new->data.qmmm=0;
		}
		new->data.pos[0]=atof(cont4);
		new->data.pos[1]=atof(cont5);
		new->data.pos[2]=atof(cont6);
		if(frame==1)
		  addAtom(solv);
	      }
	      strcpy(cont7,blank);
	      strcpy(cont8,blank);
	    }// if 1
	    j++;
	  }// while nratoms 
	
	  fgets(fcryst,SLEN,fpi_g);
    	  sscanf(fcryst,"     %s %s %s",cont1,cont2,cont3);
	  CRYST[0]=atof(cont1); CRYST[1]=atof(cont2); CRYST[2]=atof(cont3);
 
	  if(frame==1) {
	    delAtom(atom,listlength(atom));delAtom(solv,listlength(solv));
	    lna=listlength(atom);
	    lns=listlength(solv);
	    printf("\nATOMLIST IS %d / %d FILLED\n",lna,FULL);
	    printf("\nSOLVLIST IS %d / %d FILLED\n",lns,FULL);
	    fprintf(stderr,"\nplease be patient...\nWORKING...\n");
	    mkListHash(hashAtom,atom);
	    mkListHash(hashSolv,solv);
	  }

	  if(a && frame>=startfr){ 
	    for(l=0;l<=steps;l++){ //R loop
	      rx=r;     
	      //IS QMCENTER CLOSE TO THE BOUNDARY
	      for(sa=0;sa<3;sa++){
		DISPR[sa]=0.0;T[sa]=0;
	      }
	      for(sa=0;sa<3;sa++) {
		if(QMcenter[sa]+rx>CRYST[sa]){
		  DISPR[sa]=QMcenter[sa]+rx-CRYST[sa];
		  T[sa]=-1;
		}
		if(QMcenter[sa]<rx){
		  DISPR[sa]=abs(CRYST[sa]-rx);
		  T[sa]=1;
		}
	      }
	      //INITS
	      j=0;m=0;
	      for(q=0;q<=12;q++){
		CAL->dat.B[q]=0;
		CALs->dat.B[q]=0;
	      }
	      
	      search_top(hashAtom,lna,QMcenter,T,CRYST,DISPR, RES,r,CAL);
	      search_top(hashSolv,lns,QMcenter,T,CRYST,DISPR,RESs,r,CALs);

	      RHO->dat.D[0]+=CAL->dat.B[0]; RHO->dat.D[1]+=CAL->dat.B[1];
	      RHO->dat.D[2]+=CAL->dat.B[6]; RHO->dat.D[3]+=CAL->dat.B[7];
	      RHO->dat.D[4]+=CAL->dat.B[8]; RHO->dat.D[5]+=CAL->dat.B[5];
	      RHO->dat.D[10]+=CAL->dat.B[10];

	      RHOs->dat.D[0]+=CALs->dat.B[0]; RHOs->dat.D[1]+=CALs->dat.B[1];
	      RHOs->dat.D[2]+=CALs->dat.B[6]; RHOs->dat.D[3]+=CALs->dat.B[7];
	      RHOs->dat.D[4]+=CALs->dat.B[8]; RHOs->dat.D[5]+=CALs->dat.B[5];
	      RHOs->dat.D[10]+=CALs->dat.B[10];

	      fwrite(RHO,sizeof(DENSITY),1,frho);
	      fwrite(RHOs,sizeof(DENSITY),1,frho);
	      r+=dr;
	    } //for l
	  } //if a
	  r=r0; //make ready for next frame
	  for(q=0;q<=12;q++){
	    RHO->dat.D[q]=0;	    
	    RHOs->dat.D[q]=0;
	  }
	  for(q=0;q<MAXL;q++){
	    RES->dat.H[5][q]=0;
	    RES->dat.H[6][q]=0;
	    RES->dat.H[7][q]=0;
	    RES->dat.H[8][q]=0;
	    RES->dat.H[10][q]=0;
	    RES->dat.H[11][q]=0;

	    RESs->dat.H[5][q]=0;
	    RESs->dat.H[6][q]=0;
	    RESs->dat.H[7][q]=0;
	    RESs->dat.H[8][q]=0;
            RESs->dat.H[10][q]=0;
            RESs->dat.H[11][q]=0;
	  }
	  for(q=0;q<lna;q++){
	    hashAtom[q]->data.qmmm=0;
	    hashAtom[q]->data.edge=0;
	  }
	  for(q=0;q<lns;q++){
	    hashSolv[q]->data.qmmm=0;
	    hashSolv[q]->data.edge=0;
	  }
	  frame++;
	} // while frames
	
	fclose(frho);
	frho=fopen("rho_bin","rb");
	rewind(frho);

	//CALC MEAN
	N=NFRAMES-startfr+1;
	for(i=(startfr-1);i<NFRAMES;i++){
	  for(j=0;j<=steps;j++){
	    fread(RHO,sizeof(DENSITY),1,frho);
	    fread(RHOs,sizeof(DENSITY),1,frho);   
	    rho_tot[0][j]+=RHO->dat.D[5]/N;
	    rho_tot[1][j]+=RHO->dat.D[2]/N;
	    rho_tot[2][j]+=RHO->dat.D[3]/N;
	    rho_tot[3][j]+=RHO->dat.D[4]/N;
	    rho_tot[8][j]+=RHO->dat.D[10]/N;
	    
	    rho_tot[4][j]+=RHOs->dat.D[5]/N;
	    rho_tot[5][j]+=RHOs->dat.D[2]/N;
	    rho_tot[6][j]+=RHOs->dat.D[3]/N;
	    rho_tot[7][j]+=RHOs->dat.D[4]/N;
	    rho_tot[9][j]+=RHOs->dat.D[10]/N;
	  }
	}
	for(j=0;j<=steps;j++){
	  fprintf(fpo,"%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",r0+j*dr,rho_tot[0][j],rho_tot[1][j],rho_tot[2][j],rho_tot[3][j],rho_tot[8][j]);
	  fprintf(fpos,"%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",r0+j*dr,rho_tot[4][j],rho_tot[5][j],rho_tot[6][j],rho_tot[7][j],rho_tot[9][j]);
	}
	printf("\n");
	if(fret){
	  strcpy(prefix,fname_o);
	  fpi_h = fopen(strcat(prefix,"_histb"),"w");
	  strcpy(prefix,fname_o);
	  fpi_hs= fopen(strcat(prefix,"_hists"),"w");
	  for(i=0;i<MAXL;i++){
	    j=(int)RES->dat.H[0][i];	//all hits
	    m=(int)RESs->dat.H[2][i];	//all non hydrogen
	    k=(int)RES->dat.H[2][i];	//all non hydrogen
	    n=(int)RES->dat.H[9][i];	//all non hydrogen <0.3
	    o=(int)RESs->dat.H[9][i];
	    if(j>0){
	      fprintf(fpi_h,"% 10d% 10d% 10d% 10d %4s %4s\n",i,j,k,n,hashAtom[i-1]->data.resi,hashAtom[i-1]->data.name);
	    }
	    if(m>0)
	      fprintf(fpi_hs,"% 10d% 10d% 10d %4s\n",i,m,o,hashSolv[i-1]->data.name);
	  }
	  fclose(fpi_h);
	  fclose(fpi_hs);
	}

	if(xanes>0){
	  for(i=0;i<400;i++)
	    AMOUNT[i]=0;
	  fpi_h=fopen("xafs.inp","w+");
	  if(xanes==1){
	    fprintf(fpi_h,"*SAMPEL FEFF FILE FOR XANES\n");
	    fprintf(fpi_h,"TITLE PROTEIN STRUCTURE\n\n");
	    fprintf(fpi_h,"XANES\n");
	    fprintf(fpi_h,"DEBYE   190  315     PROT AT 190K, DEBYE TEMP 315K\n");
	    fprintf(fpi_h,"\n* PROTEIN IS AMORPHOUS\n");
	  }
	  fprintf(fpi_h,"\nPOTENTIALS\n");
	  //COMPLETED XANES HEADER
	  if(xanes>1)
	    fprintf(fpi_h,"*    ipot   Z  element\n");
	  new=hashAtom[Cndx-1];
	  sprintf(Z,"%s",new->data.name);
	  q=getAtomNr(Z);
	  if(xanes==1)
	    fprintf(fpi_h,"%4d%5d\n",0,q);
	  else
	    fprintf(fpi_h,"   %4d%6d   %s\n",0,q,ztags.name[30]);
	  TYPES[0]=q;n=0;
	  TYPES[Cndx-1]=q;
	  ZPOT[q]=n;
	  POTEN[Cndx-1]=ZPOT[TYPES[Cndx-1]];
	  POTEN[0]=ZPOT[TYPES[Cndx-1]];
	  for(i=0;i<lna;i++){
	    if(!(i==Cndx-1)){
	      at=hashAtom[i];
	      sprintf(Z,"%s",at->data.name);
	      q=getAtomNr(Z);
	      m=1;
	      if(!(q==1)){
		for(j=0;j<i;j++)
		  if(q==TYPES[j])
		    m=0;
		TYPES[i]=q;
		if(m){
		  n++;
		  ZPOT[q]=n;
		  if(xanes==1)
		    fprintf(fpi_h,"%4d%5d\n",ZPOT[TYPES[i]],TYPES[i]);
		  else
		    fprintf(fpi_h,"   %4d%6d   %s\n",ZPOT[TYPES[i]],TYPES[i],ztags.name[q]);
		}
		POTEN[i]=ZPOT[TYPES[i]];
	      }
	    }
	  }
	  //POTENTIAL DONE
	  fprintf(stderr,"\nXANES ATOMS\n");
	  fprintf(fpi_h,"\nATOMS\n");
	  if(xanes>1)
	    fprintf(fpi_h,"*   x          y          z      ipot  tag              distance\n");
	  at=hashAtom[Cndx-1];
	  solrsnr=0;
	  for(i=0;i<lns;i++){
	    at=hashSolv[i];
	    r=sqrt( pow(at->data.pos[0]-new->data.pos[0],2) +
			  pow(at->data.pos[1]-new->data.pos[1],2) +
			  pow(at->data.pos[2]-new->data.pos[2],2))*10;
	    if(r<30){ // only search within 30A AND KEEP ALL VECTORS SAFE!
	      sprintf(Z,"%s",at->data.name);
	      q=getAtomNr(Z);
	      if(q>1){ //THIS MEANS THAT WE FOUND OXYGEN SINCE ONLY WATER HAS BEEN CONSIDERED HERE
                Rs[solrsnr*2]=r;
	        Rs[solrsnr*2+1]=i;
		solrsnr++;
	      }
	    }
	  }
	  steps=0;

	  fprintf(stderr,"need to caluclate distances and order the indices.\n");
	  for(i=0;i<lna;i++){
	    at=hashAtom[i];
	    sprintf(Z,"%s",at->data.name);
	    q=getAtomNr(Z);
	    if(q==16 & xanes==3)
		rx=1.1;
	    else
		rx=1.00;
	    Rd[i*2]=sqrt( pow(at->data.pos[0]-new->data.pos[0],2) +
			  pow(at->data.pos[1]-new->data.pos[1],2) +
			  pow(at->data.pos[2]-new->data.pos[2],2))*10*rx;
	    Rd[i*2+1]=i;
	  }
	  sort(Rd,lna);
	  sort(Rs,solrsnr);
	  fprintf(stderr,"sorted out\n");
	  if(xanes==1)
	    fprintf(fpi_h,"%9.4f  %9.4f  %9.4f %5d %9.4f \n",
		    (at->data.pos[0]-new->data.pos[0]),
		    (at->data.pos[1]-new->data.pos[1]),
		    (at->data.pos[2]-new->data.pos[2]),
		    ZPOT[TYPES[Cndx-1]],
		    sqrt( pow(at->data.pos[0]-new->data.pos[0],2) +
			  pow(at->data.pos[1]-new->data.pos[1],2) +
			  pow(at->data.pos[2]-new->data.pos[2],2)) );
	  else{
	    steps++;
	    q=TYPES[Cndx-1];
	    strcpy(cont1,ztags.name[TYPES[Cndx-1]]);
	    while((chrptr=strpbrk(cont1," "))!=NULL)
	      *chrptr='x';
	    cont2[0]=tolower(cont1[0]);cont1[0]=cont2[0];

	    if(q==16 & xanes==3)
		rx=1.08;
	    else
		rx=1.00;	    
	    fprintf(fpi_h,"% 9.4f  % 9.4f  % 9.4f %3d %s_%d %18.5f \n",
		    (at->data.pos[0]-new->data.pos[0])*10*rx ,
		    (at->data.pos[1]-new->data.pos[1])*10*rx ,
		    (at->data.pos[2]-new->data.pos[2])*10*rx ,
		    ZPOT[TYPES[Cndx-1]],cont1,++AMOUNT[TYPES[Cndx-1]],
		    sqrt( pow(at->data.pos[0]-new->data.pos[0],2) +
			  pow(at->data.pos[1]-new->data.pos[1],2) +
			  pow(at->data.pos[2]-new->data.pos[2],2))*10*rx); 
	  }
	  fprintf(stderr,"did some stuff\n");
	  for(i=0;i<lna;i++){
	    if(!( i==(Cndx-1) | i==0)){
	      j=Rd[i*2+1];
	      if(solrsnr>0)
              for(k=0;k<solrsnr;k++){
		if(Rs[k*2]<Rd[i*2]){
		  if(steps<500){
		    //FOUND ONE SO PRINT IT FIRST AND THEN NULL THE ENTRY
		    at=hashSolv[(int)Rs[k*2+1]];
		    fprintf(fpi_h,"%9.4f  %9.4f  %9.4f %3d %s_%d %18.5f \n",
			    (at->data.pos[0]-new->data.pos[0])*10,
			    (at->data.pos[1]-new->data.pos[1])*10,
			    (at->data.pos[2]-new->data.pos[2])*10,
			    ZPOT[8],"ox",++AMOUNT[8],
			    sqrt( pow(at->data.pos[0]-new->data.pos[0],2) +
				  pow(at->data.pos[1]-new->data.pos[1],2) +
				  pow(at->data.pos[2]-new->data.pos[2],2))*10 );
		    Rs[k*2]=1000000;
		    steps++;
		  }
		}
	      }
	      fprintf(stderr,"PERHAPS INTRODUCE THE SAME LIMIT ON PROTEIN\n");
	      at=hashAtom[j];
	      sprintf(Z,"%s",at->data.name);
	      q=getAtomNr(Z);
	      if(q==16 & xanes==3)
		  rx=1.08;
	      else
		  rx=1.00;
	      fprintf(stderr,"Working some more\n\n");	      
	      if(!(q==1)){
		if(xanes==1)
		  fprintf(fpi_h,"%9.4f %9.4f %9.4f %5d %9.4f \n",
			  (at->data.pos[0]-new->data.pos[0])*10,
			  (at->data.pos[1]-new->data.pos[1])*10,
			  (at->data.pos[2]-new->data.pos[2])*10,
			  ZPOT[TYPES[j]],
			  sqrt( pow(at->data.pos[0]-new->data.pos[0],2) +
				pow(at->data.pos[1]-new->data.pos[1],2) +
				pow(at->data.pos[2]-new->data.pos[2],2))*10 );
		else{
		  if(steps<500){
		    strcpy(cont1,ztags.name[TYPES[j]]);
		    while((chrptr=strpbrk(cont1," "))!=NULL)
		      *chrptr='x';
		    cont2[0]=tolower(cont1[0]);cont1[0]=cont2[0];
		    fprintf(fpi_h,"%9.4f  %9.4f  %9.4f %3d %s_%d %18.5f \n",
			    (at->data.pos[0]-new->data.pos[0])*10*rx,
			    (at->data.pos[1]-new->data.pos[1])*10*rx,
			    (at->data.pos[2]-new->data.pos[2])*10*rx,
			    ZPOT[TYPES[j]],cont1,++AMOUNT[TYPES[j]],
			    sqrt( pow(at->data.pos[0]-new->data.pos[0],2) +
				  pow(at->data.pos[1]-new->data.pos[1],2) +
				  pow(at->data.pos[2]-new->data.pos[2],2))*10*rx );
		    steps++;
		  } 
		}
	      }
	    }
	  }
	  if(xanes>1)
	    fprintf(fpi_h,"END\n");
	  fclose(fpi_h);
	}
	fclose(frho);
	fclose(fpi_g);
	fclose(fpo);
	fclose(fpos);
	remove("rho_bin");
	printf("\n");
	return 0;
}
