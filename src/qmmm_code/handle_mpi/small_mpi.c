#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "argp.h"

#define CLEN 1000
#define NATOMS 20

int fatal( char str[]){
  fprintf(stderr,"%s\n",str);
  exit(1);
}

/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
  char *args[8];            /* ARG1 */
  int verbose;              /* The -v flag */  
};

static struct argp_option options[] =
  {
    {"verbose", 'v', 0, 0, "Be verbose"},
    {0}
  };

static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  struct arguments *arguments = (struct arguments*)state->input;
  
  switch (key)
    {
    case 'v':
      arguments->verbose = 1;
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num >= 8)
	{
	  argp_usage(state);
	}
      arguments->args[state->arg_num] = arg;
      break;
    case ARGP_KEY_END:
      if (state->arg_num < 8)
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
static char args_doc[] = "NCLUSTERS EXE INM INF INEXT OUTM OUTF OUTEXT";  // "ARG1 ARG2 ARG3 ARG4 ARG5 ARG6";

static char doc[] =
"THIS IS A SMALL MPI PROGRAM FOR DISHING OUT DATAJOBS TO GAUSSIAN \nEx.: mpirun -np 4 ./smpi  6 /lap/gaussian/03.e.01/g03/g03 \"<\" input com \">\" output out ";

static struct argp argp = {options, parse_opt, args_doc, doc};

int main ( int argc, char *argv[] )
{
  struct arguments arguments;
  int dest;
  int i,j,NCLUST;
  double res;
  int ierr;
  int master=0;
  int flag=1;
  int my_id;
  int num_procs;
  int source;
  MPI_Status status;
  int tag;
  char *p,*exe,*in_op,*out_op,*iname,*oname;
  char *in_ext,*out_ext;
  int  n;
  FILE *npalog;
  char line[CLEN],straeng[CLEN],exeline[CLEN];
  char st1[20],st2[50],buf[255];
  char enter_str[255],exit_str[255];
  char inname[255],outname[255],in2name[255],ut2namn[255];
  FILE *log,*com,*newcom;
  int  nrcase=0;

  /* DEFAULTS */
  arguments.verbose = 0;
  
  /* PARSER */
  argp_parse (&argp, argc, argv, 0, 0, &arguments);
  NCLUST  = atoi(arguments.args[0]);
  exe     = arguments.args[1];
  in_op   = arguments.args[2];
  iname   = arguments.args[3];
  in_ext  = arguments.args[4];
  out_op  = arguments.args[5];
  oname   = arguments.args[6];
  out_ext = arguments.args[7];
  
  arguments.verbose=1;
  /* If in verbose mode */
  if (arguments.verbose>1)
    fprintf(stdout,"%s\n","VERBOSE MODE");
  
  //  Initialize MPI.
  ierr = MPI_Init ( &argc, &argv );
  
  //  Get the number of processes.
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );
  
  //  Get the individual process ID.
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &my_id );
  
  //  Print a message.
  if ( my_id == master && arguments.verbose > 1) 
    {
      printf ( "\n" );
      printf ( "TYPE - Master process:\n" );
      printf ( "  C version\n" );
      printf ( "  An simple MPI program that executes a binary\n" );
      printf ( "\n" );
      printf ( "  Compiled on %s at %s.\n", __DATE__, __TIME__ );
      printf ( "\n" );
      printf ( "  The number of processes is %d\n", num_procs );
    }
  
  sprintf(exit_str,"cd ../..\n");
  for(j=0;j<NCLUST;j++){
    i=j%num_procs;
    if(my_id==i){
      sprintf(exeline,"cd  clusters/c%d && %s %s %s%d.%s %s %s%d.%s && cd ../.. \n",j,exe,in_op,iname,j,in_ext,out_op,oname,j,out_ext);  
      if(arguments.verbose)
        printf("proc(%3d)::%s",my_id,exeline);
      if ( system(exeline) != 0 ){
	fprintf(stderr,"FAILED::TESTING NEW CASE\n");
	sprintf(outname,"clusters/c%d/%s%d.%s",j,oname,j,out_ext);
	log=fopen(outname,"r+");
        nrcase=10;
	while(!feof(log)){
	  fgets(line,100,log);
	  if(!strncmp(line," NtrErr Called from FileIO.",26))
	    nrcase=1;
	  if(!strncmp(line," The SCF is confused.",20))
	    nrcase=2;
	  if(!strncmp(line," Convergence failure -- run terminated.",38))
	    nrcase=3;
	  if(!strncmp(line," Initial guess read from the checkpoint file:",45))
	    nrcase=4;
	  if(!strncmp(line," RdChkP:  Unable to locate IRWF",31))
	    nrcase=5;
	  if(!strncmp(line," Rare condition: small coef for last iteration:",48))
	    nrcase=6;
	  if(!strncmp(line," Error translating basis functions",35))
	    nrcase=7;
	}
	fclose(log);
        if(arguments.verbose>1){
  	  fprintf(stderr,"GOT CASENR:%d\n",nrcase);
        }
	sprintf(inname,"clusters/c%d/%s%d.%s",j,iname,j,in_ext);
	sprintf(in2name,"clusters/c%d/%s2%d.%s",j,iname,j,in_ext);
	if(nrcase>0){
	  com=fopen(inname,"r");
	  newcom=fopen(in2name,"w");
	  switch(nrcase){
	  case 1:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s Guess=Permute Charge \n",st1,st2); //after fail we dont do PA
	      }
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  case 2: 
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s/Auto Guess=Permute SCF=(QC,MaxCyc=500) Charge \n",st1,st2);
	      }
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  case 3:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s/Auto Guess=Permute Charge SCF=(QC,MaxCyc=500) \n",st1,st2); //SCF=(QC,MaxCycle=500)
	      }
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  case 4:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s/Auto Guess=Permute Charge SCF=(QC,MaxCyc=500) \n",st1,st2);
	      }
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  case 5:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s/Auto Guess=Permute Charge SCF=(QC,MaxCyc=500) \n",st1,st2);
	      }
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  case 6:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s/Auto Guess=Permute Charge SCF=(QC,MaxCyc=500)\n",st1,st2);
	      }
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  case 7:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s/Auto Guess=Permute Charge SCF=(QC,MaxCyc=500)\n",st1,st2);
	      }
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  default:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s/Auto Guess=Permute Charge SCF=(QC,MaxCyc=500)\n",st1,st2);
	      }
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  }
	  fclose(com);
	  fclose(newcom);
	  sprintf(buf,"cat %s | grep -v chk > %s ",in2name,inname);
	  system(buf);
          printf("proc(%3d)##%s",my_id,exeline);
	  if ( system(exeline) != 0 ){
	    fprintf(stderr,"Call to '%s' failed twice\n",exeline);
	    exit(1);
	  }
	  //	  system(exit_str);
	}
      }
      
      /*
	sprintf(line,"cd clusters/c%d && %s %s %s%d.%s %s %s%d.%s && cd ../.. \n",j,exe,in_op,iname,j,in_ext,out_op,oname,j,out_ext);  
      if(arguments.verbose)
        printf("proc(%3d)::%s",my_id,line);
      if ( system(line) != 0 ){
	fprintf(stderr,"FAILED::TESTING NEW CASE\n");
	sprintf(outname,"%s%d.%s",oname,j,out_ext);
	log=fopen(outname,"r+");
	while(!feof(log)){
	  fgets(line,100,log);
	  if(!strncmp(line," NtrErr Called from FileIO.",26))
	    nrcase=1;
	  if(!strncmp(line," The SCF is confused.",20))
	  nrcase=2;
	  if(!strncmp(line," Convergence failure -- run terminated.",38))
	    nrcase=3;
	  if(!strncmp(line," Initial guess read from the checkpoint file:",45))
	    nrcase=4;
	  if(!strncmp(line," RdChkP:  Unable to locate IRWF",31))
	    nrcase=5;
	  if(!strncmp(line," Rare condition: small coef for last iteration:",48))
	    nrcase=6;
	  if(!strncmp(line," Error translating basis functions",35))
	    nrcase=7;
	}
	fclose(log);
	
	sprintf(inname,"%s%d.%s",iname,j,in_ext);
	sprintf(in2name,"%s2%d.%s",iname,j,in_ext);
	if(nrcase>0){
	  com=fopen(inname,"r");
	  newcom=fopen(in2name,"w");
	  switch(nrcase){
	  case 1:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
	      }
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  case 2: 
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s Guess=Mix Pop=MK\n",st1,st2);
	      }
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  case 3:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
		}
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  case 4:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
	      }
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  case 5:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
	  }
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  case 6:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
		}
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  case 7:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
		}
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  default:
	    while(!feof(com)){
	      fgets(line,100,com);
	      if(line[0]=='#'){
		sscanf(line,"%s %s",st1,st2);
		fprintf(newcom,"%s %s Guess=Mix Charge Pop=MK\n",st1,st2);
		}
	      else
		fprintf(newcom,"%s",line);
	    }
	    break;
	  }
	  fclose(com);
	  fclose(newcom);
	  sprintf(buf,"mv clusters/c%d/%s clusters/c%d/%s ",j,in2name,j,inname);
	  system(buf);
	  if ( system(line) != 0 ){
              fprintf(stderr,"Call to '%s' failed twice\n",line);
              exit(1);
	  }
//	  system(exit_str);
	}
      }*/
    }
  }
  
  //  Shut down MPI.
  ierr = MPI_Finalize ( );
  
  if ( my_id == master && arguments.verbose >1) 
    {
      printf ( "\n" );
      printf ( "TYPE - Master process:\n" );
      printf ( "  Normal end of execution.\n" );
    }
  
  return 0;
}

