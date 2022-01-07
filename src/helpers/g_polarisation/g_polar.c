#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include "hlist.h"
#include "define.h"

const char *argp_program_version =
"g_tool v1.0 by Richard Tjornhammar";

const char *argp_program_bug_address =
"<richard.tjornhammar@gmail.com>";

/*typedef struct vect{
  float* dat ;
  int i;
} vector;

typedef vector *pvec;

enum {
  etGRO=0, etPDB=1, etXYZ=2, etTOP=3, etITP=4, etDAT=5
};

enum {
  cmdROT=0,cmdMROT=1,cmdPUT=2,cmdTRANS=3,cmdERROR=4
};
*/
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

/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
  char *args[1];            /* ARG1 and ARG2 */
  int verbose;              /* The -v flag */
  char *outfile,*infile,*infile2;    /* Argument for -o and -i*/
  char *solitp;
  char *insert;             /* Argument for -j */
  int frames,start;
};

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
      if (state->arg_num >= 1)
	{
	  argp_usage(state);
	}
      arguments->args[state->arg_num] = arg;
      break;
    case ARGP_KEY_END:
      if (state->arg_num < 1)
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
static char args_doc[] = "LABEL";//"ARG1 ARG2";

/*
  DOC.  Field 4 in ARGP.
  Program documentation.
*/
static char doc[] =
"g_gid -- A gromacs module to calculate polarisation";

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
  char *commands[]={"ROTATE","MROTATE","PUT","TRANSLATE","ERROR"};
  char *cmd;

  //COORDINATES
  LINK  atom,solv,itpc,solvh[MAXL];
  float CELL[6];
  int   lna,lns;
  
  //gID
  LINK CAT;
  pvec  gID,cID,Ndens,P,xi,Mdens;
  int  datSIZE,ir,v=0;
  float dr=0.02,r=0.10,NSTEPS=0,mL,rp,rmax,rho=33.427;
  int  ifr,tframe,sframe;
  float RESULTS[6][1000]; //THIS IS BAD CODING
  
  /* DEFAULTS */
  arguments.outfile = NULL;
  arguments.infile = NULL;
  arguments.infile2 = NULL;
  arguments.solitp  = NULL;
  arguments.verbose = 0;
  arguments.frames = 1;
  
  /* PARSER */
  argp_parse (&argp, argc, argv, 0, 0, &arguments);
  
  cmd=arguments.args[0];

  i=0;
  while(!(cmd[i]=='\0')){
    cmd[i]=toupper(cmd[i]);
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
    if(file_ext(arguments.infile,"gro")&file_ext(arguments.infile,"pdb"))
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
  //else
  //  fatal("MUST SPECIFY ITP FILE");
  
  //GET COORDINATES
  atom=initAtom();
  solv=initAtom();
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

  if(arguments.verbose)
    fprintf(stdout,"\nSETTING UPP FOR CORRELATION");

  mL=(CELL[0]+CELL[1]+CELL[2])/3.0;
  rmax=mL*0.5;
  datSIZE=floor((rmax-r)/(dr));
  NSTEPS=datSIZE;
  
  if(arguments.verbose)
    fprintf(stdout,"\nALLOCATING WORKSPACE\nSIZE:%d mL:%f a:%f r:%f rmax:%f dr:%f\n",datSIZE,mL,CELL[0],r,rmax,dr);
 
  if(arguments.verbose)
    fprintf(stdout,"\nSETTING UP gID");
  gID = malloc(sizeof(vector));
  if(gID==NULL)
    fatal("FAILED IN ALLOCATING MEMORY");
  if(arguments.verbose)
    fprintf(stdout,"\nSETTING UP gID -> dat\n"); 
  gID-> dat = malloc(sizeof(float)*datSIZE);
  if(gID->dat==NULL)
    fatal("FAILED IN ALLOCATING MEMORY");

  if(arguments.verbose)
    fprintf(stdout,"SETTING UP cID\n");
  cID=malloc(sizeof(double));
  if(cID==NULL)
    fatal("FAILED IN ALLOCATING MEMORY");
  if(arguments.verbose)
    fprintf(stdout,"SETTING UP cID->dat\n");
  cID-> dat = malloc(sizeof(float)*datSIZE);
  if(cID->dat==NULL)
    fatal("FAILED IN ALLOCATING MEMORY");

  if(arguments.verbose)
    fprintf(stdout,"SETTING UP Ndens\n");
  Ndens=malloc(sizeof(double));
  if(Ndens==NULL)
    fatal("FAILED IN ALLOCATING MEMORY");
  if(arguments.verbose)
    fprintf(stdout,"SETTING UP Ndens -> dat\n");
  Ndens-> dat = malloc(sizeof(float)*datSIZE);
  if(Ndens->dat==NULL)
    fatal("FAILED IN ALLOCATING MEMORY");

  if(arguments.verbose)
    fprintf(stdout,"SETTING UP Mdens\n");
  Mdens=malloc(sizeof(double));
  if(Mdens==NULL)
    fatal("FAILED IN ALLOCATING MEMORY");
  if(arguments.verbose)
    fprintf(stdout,"SETTING UP Mdens -> dat\n");
  Mdens-> dat = malloc(sizeof(float)*datSIZE);
  if(Mdens->dat==NULL)
    fatal("FAILED IN ALLOCATING MEMORY");

  if(arguments.verbose)
    fprintf(stdout,"SETTING UP P\n");
  P=malloc(sizeof(double)); 
  if(P==NULL)
    fatal("FAILED IN ALLOCATING MEMORY");
  if(arguments.verbose)
    fprintf(stdout,"SETTING UP P -> dat \n");
  P -> dat = malloc(sizeof(float)*datSIZE); 
  if(P->dat==NULL)
    fatal("FAILED IN ALLOCATING MEMORY");

  if(arguments.verbose)
    fprintf(stdout,"SETTING UP xi \n");
  xi = malloc(sizeof(double)); 
  if(xi==NULL)
    fatal("FAILED IN ALLOCATING MEMORY");
  if(arguments.verbose)
    fprintf(stdout,"SETTING UP xi -> dat\n");
  xi -> dat = malloc(sizeof(float)*datSIZE); 
  if(xi->dat==NULL)
    fatal("FAILED IN ALLOCATING MEMORY");
  
  if(arguments.verbose)
    fprintf(stdout,"\nSEARCHING FOR ATOM LABEL::%s\n",cmd);
  CAT=centerAtom(atom,cmd);
  printAtom(CAT);
  if(arguments.verbose)
    fprintf(stdout,"\n");

  //FILL EXISTING STRUCTURES WITH CHARGE AND FF DATA
  if(do_charge)
    assign_topology(atom,solv,itpc,intop,do_charge);

  //ATOM AND SOLV STRUCTURE NOW CONTAINS DATA
  if (arguments.verbose)
    fprintf(stdout,"\nINIT COMPLETE\n");

  mkListHash(solvh,solv);

  tframe=arguments.frames;
  v=0;
  for(ifr=0;ifr<tframe;ifr++){
    if (arguments.verbose)
      fprintf(stdout,"BEGIN INIT PVEC, frame::%d / %d \n",ifr+1,tframe);

    initPvec(gID,datSIZE);
    initPvec(cID,datSIZE);
    initPvec(P,datSIZE);
    initPvec(Ndens,datSIZE);
    initPvec(Mdens,datSIZE);
    initPvec(xi,datSIZE);
   
    if (arguments.verbose)
      fprintf(stdout,"END INIT PVEC \n");

    for(ir=0;ir<NSTEPS;ir++){
      Ndens->i=ir;
      gID->i=ir;
      cID->i=i;
      P->i=ir;
      xi->i=ir;

      rp=r+dr*ir;
      RESULTS[0][ir]=rp;
     
      calc_ndens(solvh,CAT,rp,dr,Ndens,xi,lns,CELL,Mdens,v);
      RESULTS[1][ir]+=Ndens->dat[ir]/tframe;
      RESULTS[2][ir]+=xi->dat[ir]/tframe;

      calc_gID(Ndens,rp,dr,rho,gID,v);
      RESULTS[3][ir]+=gID->dat[ir]/tframe;

      if(Ndens->dat[ir] > 0.0)
	cID->dat[ir]=(Mdens->dat[ir]/Ndens->dat[ir])/tframe;
      else
	cID->dat[ir]=0;
      RESULTS[4][ir]+=cID->dat[ir];

      calc_Polarisation(P,gID,cID,rho,v);
      RESULTS[5][ir]+=P->dat[ir]/tframe;

    }
    update_topology(atom,solv,CAT,in,typ,CELL);
  }

  printResults(RESULTS,out,NSTEPS);

  free(gID->dat);
  free(gID);
  free(cID->dat);
  free(cID);
  free(Ndens->dat);
  free(Ndens);
  free(Mdens->dat);
  free(Mdens);
  free(P->dat);
  free(P);
  free(xi->dat);
  free(xi);
  fclose(in);
  fclose(out);
  return 0;
}
