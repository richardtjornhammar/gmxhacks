#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include "hlist.h"
#include "define.h"

const char *argp_program_version =
"g_polar v1.0 by Richard Tjornhammar";

const char *argp_program_bug_address =
"<richard.tjornhammar@gmail.com>";

/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
  char *args[7];            /* ARG1 and ARG2 */
  int verbose;              /* The -v flag */
  char *outfile,*infile,*infile2;    /* Argument for -o and -i*/
  char *solitp;
  char *insert;             /* Argument for -j */
  int frames,start;
};

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
      if (state->arg_num >= 7)
	{
	  argp_usage(state);
	}
      arguments->args[state->arg_num] = arg;
      break;
    case ARGP_KEY_END:
      if (state->arg_num < 7)
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
static char args_doc[] = "LABEL X Y Z e1 e2 e3";//"ARG1 ARG2";

/*
  DOC.  Field 4 in ARGP.
  Program documentation.
*/
static char doc[] =
"g_orient -- A gromacs module to calculate volume polarisation\n   C12::D dot E\n    C11::|Pol_tot|\n    C10::Pol_tot\n    C9 ::Pol_w\n    C8 ::Pol_c\n    C7 ::E\n    C6 ::PHK\n    C5 ::cID\n    C4 ::M\n    C3 ::Xi\n    C2 ::N\n       C1 ::r\n ";
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
  LINK  atom,solv,itpc,solvh[MAXL],atomh[MAXL],dummy;
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

  assign_vect(arguments,center,nhat,&xyz);
  if (arguments.verbose)
    fprintf(stdout,"NHAT:: %f %f %f | xyz::%d\n",nhat[0],nhat[1],nhat[2],xyz);    

  for(i=0;i<3;i++)
    rmax+=SQ(center[i]);
  rlen=sqrt(rmax);

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

  if(arguments.verbose)
    fprintf(stdout,"\nSETTING UPP FOR CORRELATION CALCULATIONS\n");

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

  //CHANGE dr FOR SMALLER BINS
  datSIZE=floor((rmax-r)/(dr));
  NSTEPS=datSIZE;

  if(datSIZE>RESY)
    fatal("MAKE RESY IN DEFINE.H LARGER AND RECOMPILE");
  
  if(arguments.verbose)
    fprintf(stdout,"\nALLOCATING WORKSPACE\nSIZE:%d mL:%f a:%f r:%f rmax:%f dr:%f\n",datSIZE,mL,CELL[0],r,rmax,dr);

  // VOL || Nd,Md,xi,cPD,gPD,PHK,P;

  Nd = setup_pvec(datSIZE,arguments);
  Md = setup_pvec(datSIZE,arguments);
  xi = setup_pvec(datSIZE,arguments);
  cPD = setup_pvec(datSIZE,arguments);
  gPD = setup_pvec(datSIZE,arguments); 
  PHK = setup_pvec(datSIZE,arguments); 
  P = setup_pvec(datSIZE,arguments); 
  Pa = setup_pvec(datSIZE,arguments);
  EF = setup_pvec(datSIZE,arguments);

  if(arguments.verbose)
    fprintf(stdout,"\nSEARCHING FOR ATOM LABEL::%s\n",cmd);

  mkListHash(solvh,solv);
  mkListHash(atomh,atom);

  if(!(atom==NULL)){
    CAT=centerAtom(atom,cmd);
    if(CAT==NULL){
      fprintf(stdout,"WARNING:: DEPRECATED ATOM LABEL\n");
      CAT = dummy;
      CAT->data.edge=0; //INIT
      if(rlen>0.0){
	CAT->data.pos[0]=center[0];
	CAT->data.pos[1]=center[1];
	CAT->data.pos[2]=center[2];
	CAT->data.edge=0;
      }
      else {
	corner_search(solvh,CELL);
	CAT->data.edge=1;
	CAT->data.pos[0]=CELL[3];
	CAT->data.pos[1]=CELL[4];
	CAT->data.pos[2]=CELL[5];	
      }
    }
    else {
      CAT->data.edge=0;
    }
  }
  else{
    fprintf(stdout,"WARNING:: EMPTY ATOM LIST\n");
    CAT=dummy;
    if(rlen>0.0){
      CAT->data.pos[0]=center[0];
      CAT->data.pos[1]=center[1];
      CAT->data.pos[2]=center[2];
      CAT->data.edge=0;
    }
    else{
      corner_search(solvh,CELL);
      CAT->data.edge=1;
      CAT->data.pos[0]=CELL[3];
      CAT->data.pos[1]=CELL[4];
      CAT->data.pos[2]=CELL[5];
    }
  }

  if(arguments.verbose & !(CAT==NULL)){
    printAtom(CAT);
    printAtom(solv);
    printAtom(solv->next);
    printAtom(solv->next->next);
    printAtom(solv->next->next->next);
    printAtom(solv->next->next->next->next->next);
  }
  if(arguments.verbose){
    fprintf(stdout,"ASSIGNED CENTER COORDINATE: ");
    for(j=0;j<3;j++)
      fprintf(stdout,"%f ",CAT->data.pos[j]);
    fprintf(stdout,"\n");
  }

  //ATOM AND SOLV STRUCTURE NOW CONTAINS DATA
  if (arguments.verbose)
    fprintf(stdout,"\nINIT COMPLETE\n");

  tframe=arguments.frames;
  v=0;
  for(ifr=0;ifr<tframe;ifr++){
    if (arguments.verbose)
      fprintf(stdout,"BEGIN INIT FRAME:: %d / %d \n",ifr+1,tframe);
    
    initPvec(Nd,datSIZE);  
    if (arguments.verbose)
      fprintf(stdout,"INIT FIRST COMPLETE \n");

    initPvec(Md,datSIZE);  initPvec(EF,datSIZE);
    initPvec(cPD,datSIZE); initPvec(xi,datSIZE);
    initPvec(PHK,datSIZE); initPvec(P,datSIZE); 
    initPvec(gPD,datSIZE); initPvec(Pa,datSIZE);
    
    if (arguments.verbose)
      fprintf(stdout,"END INIT FRAME \n");

    if (arguments.verbose)
      fprintf(stdout,"CALCULATE VOLUME DENSITIES FOR SOLVENT\n");
    calc_densDIR(solvh,lns,CAT,r,dr,CELL,nhat,Nd,Md,xi,cPD,gPD,PHK,P,EF,RESULTS,tframe,v,xyz,datSIZE,do_charge);
    if (arguments.verbose)
      fprintf(stdout,"CALCULATE MOLECULE POLARISATION\n");
    calc_atomPOL(atomh,lna,CAT,r,dr,CELL,nhat,Pa,EF,RESULTS,tframe,v,xyz,datSIZE,do_charge);
      
    if(tframe>1){
      update_topology(atom,solv,CAT,in,typ,CELL,diam,tframe);
      if (arguments.verbose)
	fprintf(stdout,"UPDATED TOPOLOGY\n");
    }
  }

  if(arguments.verbose)
    fprintf(stdout,"\n\nPRINTING RESULTS\n");  
  printResults(RESULTS,out,NSTEPS);
  
  fprintf(stdout,"PROTEIN D::%f\nATO-PRO D::%f\n",diam[0],diam[1]);

  if(arguments.verbose)
    fprintf(stdout,"CLOSING STREAMS\n\n");

  if(do_charge)
    fclose(intop);
  fclose(in);
  fclose(out);
  return 0;
}
