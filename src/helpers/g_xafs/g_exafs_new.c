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
  int H;
  int pbcOFF;
  float r;
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

float r12(LINK at,LINK cat){
  float r;
  r=sqrt(pow(at->data.pos[0]-cat->data.pos[0],2.0) + pow(at->data.pos[1]-cat->data.pos[1],2.0)+pow(at->data.pos[2]-cat->data.pos[2],2.0) );
  return(r);
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
  {"includeH", 'H', 0, 0, "Produce file with H atoms "},
  {"center_off", 'c',0,0,"Dont center coordinates on atom LABEL"},
  {"input",  'i', "INFILE", 0,
   "Input from INFILE instead"},
  {"output",  'o', "OUTFILE", 0,
   "Output to OUTFILE instead of to standard output"},
  {"rad",  'r', "Portube factor for sulphurs", 1,
   "scales sulphur distances radially away from center"},
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
    case 'c':
      arguments->pbcOFF = 1;
      break;
    case 'H':
      arguments->H = 1;
      break;
    case 'o':
      arguments->outfile = arg;
      break;
    case 'i':
      arguments->infile = arg;
      break;
    case 'r':
      arguments->r = atof(arg);
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
"g_exafs a small program for creating exafs input (feff6)\nCMD can be any of \"XANES\" \"EXAFS\" \"PORTRUBED\"";
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
  FILE *out,*outpdb,*in,*intop,*initp,*frad;
  int typ,what,do_charge=0;
  char *exten[]={"GRO","PDB","XYZ","TOP","DAT"};
  char *commands[]={"XANES","EXAFS","PORTRUBED"};
  char *catm,*cmd;

  //COORDINATES
  LINK  atom,solv,itpc,solvh[MAXL],atomh[MAXL],toph[2*MAXL],dummy,at,new;
  float CELL[6],a,b,c,cscale;
  int   lna,lns;
  
  //ORIENT
  LINK CAT;
  pvec Nd,Md,xi,cPD,gPD,PHK,P,*PV,Pa,EF;
  int  datSIZE,ix,iy,iz,v=0,xyz=0,ii;
  float dr=0.005,r=0.16,mL,rp,nhat[3],center[3],rlen,rmax=0.0,rho=33.427;
  float diam[2];
  int  ifr,tframe,sframe;
  float RESULTS[RESX][RESY]; 
  
  //EXAFS
  int           TYPES[4000],POTEN[4000],ZPOT[200],AMOUNT[200]; //LINES FOR XAFS DATA GENERATION
  struct ZNAME  ztags = ex_zn;
  float		radius[3],Rd[4000],Rs[10000],rx=1.0;
  int           solrsnr,q,xanes=2,Cndx,n,m,steps,NSTEPS=500;
  char          *chrptr;
  char          Z[2],cont1[MAXL],cont2[MAXL];

  /* DEFAULTS */
  arguments.outfile = NULL;
  arguments.infile = NULL;
  arguments.infile2 = NULL;
  arguments.solitp  = NULL;
  arguments.verbose = 0;
  arguments.frames = 1;
  arguments.r=1.0;
  arguments.pbcOFF=0;
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
  if(xanes==3){
    rx=arguments.r;
    if(arguments.verbose)
      fprintf(stdout,"RADIAL SCALING CONSTANT:: %f\n",rx);
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
  if (arguments.verbose){
    fprintf(stdout,"%s\n","VERBOSE");
    if(arguments.H)
      fprintf(stdout,"INCLUDE H ATOMS\n");
  }

  if(arguments.infile){
    if(!file_ext(arguments.infile,"gro")){
      typ=etGRO;
      cscale=10.0;
    }
    if(!file_ext(arguments.infile,"pdb")){
      typ=etPDB;
      cscale=1.0;
    }
    if(file_ext(arguments.infile,"gro") & file_ext(arguments.infile,"pdb"))
      fatal("INPUT FILE IS NOT A RECONGIZED FILETYPE");
    in=fopen(arguments.infile,"r");
    if(in==NULL)
      fatal("FAILED TO OPEN INPUT FILE");
  }
  else
    fatal("NO INPUT FILE GIVEN");
  
    
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
   
  //MUST BE CUBIC BOX; THE BINNING IS CUBICAL
  mL=(CELL[0]+CELL[1]+CELL[2])/3.0;

  rmax=mL;
  for(i=0;i<3;i++){
    if(arguments.verbose)
      fprintf(stdout,"CELL[%d]=%f ",i,CELL[i]);
  }
  
  if(xyz==1){
    if(arguments.verbose)
      fprintf(stdout,"FINDING LARGEST DIMENSION\n\n");
    r=0.0;rmax=0.0;
    for(i=0;i<3;i++){
      if(CELL[i]>rmax)
	rmax=CELL[i];
    }
  }

  if(arguments.verbose)
    fprintf(stdout,"\n\nSEARCHING FOR ATOM LABEL::%s\n",catm);
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
  mkListHash(solvh,solv);
  mkListHash(atomh,atom);

  //CENTER TOPOLOGY ON CENTRAL ATOM
  if(!arguments.pbcOFF){
  for(i=0;i<lns;i++){
    mL=rmPBC(solvh[i],new,CELL);  
  }
  for(i=0;i<lna;i++){
    mL=rmPBC(atomh[i],new,CELL);
  }
  }

  for(i=0;i<lns;i++){
    solvh[i]->data.r=r12(solvh[i],new);
  }
  for(i=0;i<lna;i++){
    sprintf(Z,"%s",atomh[i]->data.name);
    q=getAtomNr(Z);
    if( (q == 16) && (xanes == 3) ){
      scale(atomh[i],new,rx);
      if(arguments.verbose)
	fprintf(stdout,"SCALED ONE S\n");
    }
    atomh[i]->data.r=r12(atomh[i],new);
  }

  atomh[lna-1]->next=solvh[0];
  solvh[0]->prev=atomh[lna-1];
  mkListHash(toph,atom);
  sortLL(toph,lna+lns-1);
  NSTEPS=(NSTEPS>(lna+lns-1))?(lna+lns-1):NSTEPS;

  if (arguments.verbose)
    fprintf(stdout,"\nNSTEPS=%d\n",NSTEPS+1);

  new=toph[0];
  printAtom(new);

  //THEY ARE NOW TOGETHER AND SORTED  
  //ATOM, SOLV AND TOPH STRUCTURES NOW CONTAINS DATA
  if (arguments.verbose)
    fprintf(stdout,"\nINIT COMPLETE\n");

  for(i=0;i<4000;i++){
    if(i<200){
      ZPOT[i]=0;
      AMOUNT[i]=0;
    }
    TYPES[i]=0;
    POTEN[i]=0;
  }
  if (arguments.verbose)
    fprintf(stdout,"POTENTIALS\n");  
  fprintf(out,"POTENTIALS\n");
  n=-1;
  if(xanes>1)
    fprintf(out,"*    ipot   Z  element\n");
  for(i=0;i<NSTEPS;i++){
    at=toph[i];
    sprintf(Z,"%s",at->data.name);
    q=getAtomNr(Z);
    if(q<40){
      m=1;
      if( !(q==1) || (arguments.H & q==1) ){
	for(j=0;j<i;j++)
	  if( q == TYPES[j] )
	    m=0;
	TYPES[i]=q;
	AMOUNT[TYPES[i]]=0;
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
  if (arguments.verbose)
    fprintf(stdout,"ATOMS\n");
  fprintf(out,"\nATOMS\n");
  if(xanes>1)
    fprintf(out,"*   x          y          z      ipot  tag              distance\n");
  j=0;
  for(i=0;i<NSTEPS;i++){
    at=toph[i];
    strcpy(cont1,ztags.name[TYPES[i]]);
    while((chrptr=strpbrk(cont1," "))!=NULL)
      *chrptr='x';
    cont2[0]=tolower(cont1[0]);cont1[0]=cont2[0];
    sprintf(Z,"%s",at->data.name);
    if(at->data.name[1]=='W')
      cont1[1]='w';
    q=getAtomNr(Z);
    if( q < 100 & j < NSTEPS & !(q==1) || (arguments.H & q==1) ){
      j++;
      fprintf(out,"% 9.4f  % 9.4f  % 9.4f %3d %s_%-3d %18.5f \n",
	      (at->data.pos[0]-new->data.pos[0])*cscale ,
	      (at->data.pos[1]-new->data.pos[1])*cscale ,
	      (at->data.pos[2]-new->data.pos[2])*cscale ,
	      ZPOT[TYPES[i]],cont1,++AMOUNT[TYPES[i]],
	      at->data.r*cscale);
    }
  }
  fprintf(out,"END\n");
  tframe=arguments.frames;
  v=0;
  if(arguments.verbose)
    fprintf(stdout,"CLOSING STREAMS\n\n");

  outpdb=fopen("out0.pdb","w");
  write_pdb_tot(toph[0],outpdb,CELL,typ,NSTEPS);
  fclose(outpdb);

  if(do_charge)
    fclose(intop);
  fclose(in);
  fclose(out);
  return 0;
}
