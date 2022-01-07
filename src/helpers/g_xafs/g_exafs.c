#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include "hlist.h"
#include "define.h"

const char *argp_program_version =
"g_exafs v1.0 by Richard Tjornhammar";

const char *argp_program_bug_address =
"<richard.tjornhammar@gmail.com>";

/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
  char *args[2];            /* ARG1 and ARG2 */
  int verbose;              /* The -v flag */
  char *outfile,*infile,*infile2;    /* Argument for -o and -i*/
  char *solitp;
  char *insert;             /* Argument for -j */
  int frames,start;
};

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

void fatal(char errstr[])
{
  fprintf(stderr,"ERROR::  %s\n",errstr);
  exit(0);
}

int file_ext(char filenm[], char filext[])
{
  char* ptr;
  
  ptr=strpbrk(filenm,".");
  ptr;
  if(ptr==NULL)
    return(-1);
  else{
    ptr++;
    return(strcmp(ptr,filext));
  }
}

void assign_vect(struct arguments arguments, float r[3], float n[3],int* xyz)
{
  int i;
  float len=0.0;
  
  for(i=0;i<3;i++){
    if(arguments.args[i+1][0]=='m')
      r[i]=-1.0*atof(&arguments.args[i+1][1]);
    else
      r[i]=atof(arguments.args[i+1]);
    if(arguments.args[i+4][0]=='m')
      n[i]=-1.0*atof(&arguments.args[i+4][1]);
    else
      n[i]=atof(arguments.args[i+4]);
    len+=SQ(n[i]);
  }
  if(len>0.0){
    for(i=0;i<3;i++){
      n[i]/=sqrt(len);
    }
    xyz[0]=1;
  }
}

pvec setup_pvec(int SIZE, struct arguments arguments)
{
  pvec PVEC=NULL;

  if(arguments.verbose)
    fprintf(stdout,"SETTING UP PVEC\n");
  PVEC = malloc(sizeof(double));
  if(PVEC == NULL)
    fatal("FAILED IN ALLOCATING MEMORY");
  if(arguments.verbose)
    fprintf(stdout,"SETTING UP PVEC -> dat\n");
  PVEC -> dat = malloc(sizeof(float)*(SIZE+3)); //padding
  if(PVEC -> dat == NULL)
    fatal("FAILED IN ALLOCATING MEMORY");

  return PVEC;
}

void kill_pvec(pvec P)
{
  free(P->dat);
  free(P);
}

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC}.
*/
static struct argp_option options[] =
{
  {"verbose", 'v', 0, 0, "Produce verbose output"},
  {"input",  'i', "INFILE", 0,
   "Input from INFILE instead"},
  {"itp",  's', "INFILE", 0,
   "Charge information for solvent topology"},
  {"topol",  'p', "INFILE", 0,
   "Charge information from topology"},
  {"output",  'o', "OUTFILE", 0,
   "Output to OUTFILE instead of to standard output"},
  {"insert",  'j', "INSERT", 0,
   "Appends molecule to OUTFILE"},
  {"frames",  'f', "#FRAMES", 0,
   "updates and performs ensemble average over #FRAMES"},
  {0}
};
/*
   PARSER. Field 2 in ARGP.
   Order of parameters: KEY, ARG, STATE.
*/
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'v':
      arguments->verbose = 1;
      break;
    case 'o':
      arguments->outfile = arg;
      break;
    case 'i':
      arguments->infile = arg;
      break;
    case 'p':
      arguments->infile2 = arg;
      break;
    case 's':
      arguments->solitp = arg;
      break;
    case 'j':
      arguments->insert = arg;
      break;
    case 'f':
      arguments->frames = atoi(arg);
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num >= 2)
	{
	  argp_usage(state);
	}
      arguments->args[state->arg_num] = arg;
      break;
    case ARGP_KEY_END:
      if (state->arg_num < 2)
	{
	  argp_usage (state);
	}
      break;
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/*
   ARGS_DOC. Field 3 in ARGP.
   A description of the non-option command-line arguments
     that we accept.
*/
static char args_doc[] = "LABEL CMD";//"ARG1 ARG2";

/*
  DOC.  Field 4 in ARGP.
  Program documentation.
*/
static char doc[] =
"g_exafs a small program for creating exafs input (feff6)";
/*
   The ARGP structure itself.
*/
static struct argp argp = {options, parse_opt, args_doc, doc};

/*
   The main function.
   Notice how now the only function call needed to process
   all command-line options and arguments nicely
   is argp_parse.
*/
int main (int argc, char **argv)
{
  //GENERAL
  int i,j,k,nr;

  //IO
  struct arguments arguments;
  FILE *out,*in,*ins,*intop,*initp,*frad;
  int typ,what,do_charge=0;
  char *exten[]={"GRO","PDB","XYZ","TOP","DAT"};
  char *commands[]={"XANES","EXAFS","PORTRUBED"};
  char *catm,*cmd;

  //COORDINATES
  LINK  atom,solv,itpc,solvh[MAXL],atomh[MAXL],dummy,at,new;
  float CELL[6],a,b,c;
  int   lna,lns;
  
  //ORIENT
  LINK CAT;
  pvec Nd,Md,xi,cPD,gPD,PHK,P,*PV,Pa,EF;
  int  datSIZE,ix,iy,iz,v=0,xyz=0;
  float dr=0.005,r=0.16,NSTEPS=0,mL,rp,nhat[3],center[3],rlen,rmax=0.0,rho=33.427;
  float diam[2];
  int  ifr,tframe,sframe;
  float RESULTS[RESX][RESY]; 
  
  //EXAFS
  int           TYPES[4000],POTEN[4000],ZPOT[200],AMOUNT[200]; //LINES FOR XAFS DATA GENERATION
  struct ZNAME  ztags = ex_zn;
  float		radius[3],Rd[4000],Rs[10000],rx=0;
  int           solrsnr,q,xanes=2,Cndx,n,m,steps;
  char          *chrptr;
  char          Z[2],cont1[MAXL],cont2[MAXL];

  /* DEFAULTS */
  arguments.outfile = NULL;
  arguments.infile = NULL;
  arguments.infile2 = NULL;
  arguments.solitp  = NULL;
  arguments.verbose = 0;
  arguments.frames = 1;
  nhat[0]=0.0; nhat[1]=0.0; nhat[2]=0.0;

  /* PARSER */
  argp_parse (&argp, argc, argv, 0, 0, &arguments);
  catm= arguments.args[0];
  cmd = arguments.args[1];

  i=0;
  while(!(cmd[i]=='\0')){
    cmd[i]=toupper(cmd[i]);
    i++;
  }
  for(i=1;i<=3;i++){
    if(!(strcmp(cmd,commands[i-1]))){
      xanes=i;
      fprintf(stderr,"%s :: %d \n",commands[i-1],i);
    }
  }
  i=0;
  while(!(catm[i]=='\0')){
    catm[i]=toupper(catm[i]);
    i++;
  }  
  /* OUTPUT */
  if (arguments.outfile){
    out = fopen (arguments.outfile, "w");
    if(out==NULL)
      fatal("COULD NOT OPEN OUTPUT FILE");
  }
  else
    out = stdout;

  if (arguments.verbose)
    fprintf(stdout,"INIT STRUCTURE\n");

  v=arguments.verbose;

  /* If in verbose mode */
  if (arguments.verbose)
    fprintf(stdout,"%s\n","VERBOSE");
  
  if(arguments.infile){
    if(!file_ext(arguments.infile,"gro"))
      typ=etGRO;
    if(!file_ext(arguments.infile,"pdb"))
      typ=etPDB;
    if(file_ext(arguments.infile,"gro") & file_ext(arguments.infile,"pdb"))
      fatal("INPUT FILE IS NOT A RECONGIZED FILETYPE");
    in=fopen(arguments.infile,"r");
    if(in==NULL)
      fatal("FAILED TO OPEN INPUT FILE");
  }
  else
    fatal("NO INPUT FILE GIVEN");
  
  if(arguments.infile2)
    if(!file_ext(arguments.infile2,"top")){
      do_charge=1;
      intop=fopen(arguments.infile2,"r");
      if(intop==NULL)
	fatal("FAILED TO OPEN INPUT TOP FILE");
    }
    else
      fatal("NOT A RECOGNIZED TOP FILE");

  if(arguments.solitp){
    if(!file_ext(arguments.solitp,"itp")){
      do_charge++;
      initp=fopen(arguments.solitp,"r");
      if(initp==NULL)
	fatal("FAILED TO OPEN INPUT ITP FILE");
      itpc=initAtom();
      init_itp(itpc,initp);
      fclose(initp);
    }
    else
      fatal("NOT A RECOGNIZED ITP FILE");
  }
    
  //GET COORDINATES
  atom=initAtom();
  solv=initAtom();
  dummy=initAtom();

  switch (typ) //FILLING OF TOPOLOGY MUST BE MADE SAFE!
    {
    case etGRO: 
      gro_create(atom,solv,in,CELL);
      lna=listlength(atom); lns=listlength(solv);
      if (arguments.verbose){
	fprintf(stdout,"ATOMLIST IS %d / %d FILLED\n",lna,FULL);
	fprintf(stdout,"SOLVLIST IS %d / %d FILLED\n",lns,FULL);
      }
      break;
    case etPDB:
      pdb_create(atom,solv,in,out,CELL);
      lna=listlength(atom); lns=listlength(solv);
      if (arguments.verbose){
	fprintf(stdout,"ATOMLIST IS %d / %d FILLED\n",lna,FULL);
	fprintf(stdout,"SOLVLIST IS %d / %d FILLED\n",lns,FULL);
      }
      break;
    default:
      fatal("DIED AT SWITCH\n\n");
    }

  //FILL EXISTING STRUCTURES WITH CHARGE AND FF DATA
  if(do_charge)
    assign_topology(atom,solv,itpc,intop,do_charge);

  //MUST BE CUBIC BOX; THE BINNING IS CUBICAL
  mL=(CELL[0]+CELL[1]+CELL[2])/3.0;

  rmax=mL;
  if(xyz==1){
    if(arguments.verbose)
      fprintf(stdout,"FINDING LARGEST DIMENSION\n");
    r=0.0;rmax=0.0;
    for(i=0;i<3;i++)
      if(CELL[i]>rmax)
	rmax=CELL[i];
  }

  if(arguments.verbose)
    fprintf(stdout,"\nSEARCHING FOR ATOM LABEL::%s\n",catm);
//MUST ASSIGN CENTER ATOM
  if(!(atom==NULL)){
    new=centerAtom(atom,catm);
    Cndx=new->data.edge;
    if(arguments.verbose)
      fprintf(stdout,"FOUND:: %s       NR:: %d\n",new->data.name,Cndx+1);
  }
  else{
    fatal("DIDNT FIND CENTRAL ATOM\n\n");
  }
//END
  mkListHash(solvh,solv);
  mkListHash(atomh,atom);

  //ATOM AND SOLV STRUCTURE NOW CONTAINS DATA
  if (arguments.verbose)
    fprintf(stdout,"\nINIT COMPLETE\n");
  tframe=arguments.frames;
  v=0;

//DO WHATEVER
	Cndx++;
	if(xanes>0){
	  for(i=0;i<4000;i++){
	    if(i<200){
	      ZPOT[i]=0;
	      AMOUNT[i]=0;
	    }
	    TYPES[i]=0;
	    POTEN[i]=0;
	  }
	  if(xanes==1){
	    fprintf(out,"*SAMPEL FEFF FILE FOR XANES\n");
	    fprintf(out,"TITLE PROTEIN STRUCTURE\n\n");
	    fprintf(out,"XANES\n");
	    fprintf(out,"DEBYE   190  315     PROT AT 190K, DEBYE TEMP 315K\n");
	    fprintf(out,"\n* PROTEIN IS AMORPHOUS\n");
	  }
	  fprintf(out,"POTENTIALS\n");
          if (arguments.verbose)
	    fprintf(stdout,"\nPOTENTIALS\n");
	  //COMPLETED XANES HEADER
	  if(xanes>1)
	    fprintf(out,"*    ipot   Z  element\n");
	  sprintf(Z,"%s",new->data.name);
	  q=getAtomNr(Z);
	  if(xanes==1)
	    fprintf(out,"%4d%5d\n",0,q);
	  else{
	    fprintf(out,"   %4d%6d   %s\n",0,q,ztags.name[q]);
	  }
	  TYPES[0]=0; 
	  n=0;
	  TYPES[Cndx-1]=q;
	  AMOUNT[TYPES[Cndx-1]]=0;
	  ZPOT[q]=n;
	  POTEN[Cndx-1]=ZPOT[TYPES[Cndx-1]];
	  POTEN[0]=ZPOT[TYPES[Cndx-1]];  
	  for(i=0;i<lna;i++){
	    if(!(i==Cndx-1)){
	      at=atomh[i];
	      sprintf(Z,"%s",at->data.name);
	      q=getAtomNr(Z);
	      m=1;
	      if(!(q==1) ){
		for(j=0;j<i;j++)
		  if( q == TYPES[j] )
		    m=0;
		TYPES[i]=q;
		AMOUNT[TYPES[Cndx-1]]=0;
		if(m){
		  n++;
		  ZPOT[q]=n;
		  if (arguments.verbose){
		    fprintf(stdout,"%d %d %d %d\n",TYPES[i],i,q,ZPOT[q]);
		  }
		  if(xanes==1)
		    fprintf(out,"%4d%5d\n",ZPOT[TYPES[i]],TYPES[i]);
		  else
		    fprintf(out,"   %4d%6d   %s\n",ZPOT[TYPES[i]],TYPES[i],ztags.name[q]);
		}
		POTEN[i]=ZPOT[TYPES[i]];
	      }
	    }
	  }
	  //POTENTIAL DONE
          if (arguments.verbose)
	    fprintf(stdout,"\nXANES ATOMS\n");
	  fprintf(out,"\nATOMS\n");
	  if(xanes>1)
	    fprintf(out,"*   x          y          z      ipot  tag              distance\n");
	  at=atomh[Cndx-1];
	  solrsnr=0;
	  for(i=0;i<lns;i++){
	    at=solvh[i];
	    r=sqrt( pow(at->data.pos[0]-new->data.pos[0],2) +
			  pow(at->data.pos[1]-new->data.pos[1],2) +
			  pow(at->data.pos[2]-new->data.pos[2],2))*10;
	    if(r<30){ // only search within 30A
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
          if (arguments.verbose)
	    fprintf(stdout,"\nDONE\n");
	  //fprintf(stderr,"need to caluclate distances and order the indices.\n");
	  for(i=0;i<lna;i++){
	    at=atomh[i];
	    sprintf(Z,"%s",at->data.name);
	    q=getAtomNr(Z);
	    if(q==16 & xanes==3)
		rx=1.08;
	    else
		rx=1.00;
	    Rd[i*2]=sqrt( pow(at->data.pos[0]-new->data.pos[0],2) +
			  pow(at->data.pos[1]-new->data.pos[1],2) +
			  pow(at->data.pos[2]-new->data.pos[2],2))*10*rx;
	    Rd[i*2+1]=i;
	  }
	  sort(Rd,lna);
	  sort(Rs,solrsnr);
	  //fprintf(stderr,"sorted out\n");
	  if(xanes==1)
	    fprintf(out,"%9.4f  %9.4f  %9.4f %5d %9.4f \n",
		    (at->data.pos[0]-new->data.pos[0]),
		    (at->data.pos[1]-new->data.pos[1]),
		    (at->data.pos[2]-new->data.pos[2]),
		    ZPOT[TYPES[Cndx-1]],
		    sqrt( pow(at->data.pos[0]-new->data.pos[0],2) +
			  pow(at->data.pos[1]-new->data.pos[1],2) +
			  pow(at->data.pos[2]-new->data.pos[2],2)) );
	  else{
	    at=new;
	    steps++;
	    q=TYPES[Cndx-1];
	    strcpy(cont1,ztags.name[TYPES[Cndx-1]]);
	    while((chrptr=strpbrk(cont1," "))!=NULL)
	      *chrptr='x';
	    cont2[0]=tolower(cont1[0]);cont1[0]=cont2[0];

	    if(q==16 & xanes==3)
		rx=1.06;
	    else
		rx=1.00;
	    fprintf(out,"% 9.4f  % 9.4f  % 9.4f %3d %s_%d %18.5f \n",
		    (at->data.pos[0]-new->data.pos[0])*10*rx ,
		    (at->data.pos[1]-new->data.pos[1])*10*rx ,
		    (at->data.pos[2]-new->data.pos[2])*10*rx ,
		    ZPOT[TYPES[Cndx-1]],cont1,++AMOUNT[TYPES[Cndx-1]],
		    sqrt( pow(at->data.pos[0]-new->data.pos[0],2) +
			  pow(at->data.pos[1]-new->data.pos[1],2) +
			  pow(at->data.pos[2]-new->data.pos[2],2))*10*rx); 
	  }
	  //fprintf(stderr,"did some stuff\n");
	  for(i=0;i<lna;i++){
	    if(!( i==(Cndx-1) | i==0)){
	      j=Rd[i*2+1];
	      if(solrsnr>0)
              for(k=0;k<solrsnr;k++){
		if(Rs[k*2]<Rd[i*2]){
		  if(steps<500){
		    //FOUND ONE SO PRINT IT FIRST AND THEN NULL THE ENTRY
		    at=solvh[(int)Rs[k*2+1]];
		    fprintf(out,"%9.4f  %9.4f  %9.4f %3d %s_%d %18.5f \n",
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
	      //fprintf(stderr,"PERHAPS INTRODUCE THE SAME LIMIT ON PROTEIN\n");
	      at=atomh[j];
	      sprintf(Z,"%s",at->data.name);
	      q=getAtomNr(Z);
	      if(q==16 & xanes==3)
		  rx=1.08;
	      else
		  rx=1.00;
	      //fprintf(stderr,"Working some more\n\n");	      
	      if(!(q==1)){
		if(xanes==1)
		  fprintf(out,"%9.4f %9.4f %9.4f %5d %9.4f \n",
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
		    fprintf(out,"%9.4f  %9.4f  %9.4f %3d %s_%d %18.5f \n",
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
	    fprintf(out,"END\n");
	}

//END WHATEVER
/*
  if(arguments.outfile){
    if(!file_ext(arguments.outfile,"gro")){
      write_gro(atom,solv,out,CELL,typ);
    }
    if(!file_ext(arguments.outfile,"pdb")){
      write_pdb(atom,solv,out,out,CELL,typ);
    }
    if(file_ext(arguments.outfile,"gro") & file_ext(arguments.outfile,"pdb"))
      fatal("OUTPUT FILE IS NOT A RECONGIZED FILETYPE");
    if(in==NULL)
      fatal("FAILED TO OPEN OUTPUT FILE");
  }
  else
      write_pdb(atom,solv,out,out,CELL,typ);
*/
  if(arguments.verbose)
    fprintf(stdout,"CLOSING STREAMS\n\n");

  if(do_charge)
    fclose(intop);
  fclose(in);
  fclose(out);
  return 0;
}
