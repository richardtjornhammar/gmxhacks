#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#define STRLEN 100

void fatal(char str[])
{
  fprintf(stderr,"%s\n",str);
  exit(1);
}

int main(void){

  float dl=0.02,l0=0.0,l1=1.00,lsum=0.0;  
  int i,N;

  pid_t pid;
  int   rv,step_id;
  int   commpipe[2],stdcopy[2];            /* This holds the fd for the in & out of the pipe */
  FILE  *in,*out;
  char  *buf;
  char cmd[STRLEN],dircmd[STRLEN],direxe[STRLEN],flag[STRLEN],flag2[STRLEN];
  char selection[STRLEN],file[STRLEN],ldir[STRLEN],de[STRLEN],de2[STRLEN],dout[STRLEN];
  char file2[STRLEN],flag3[STRLEN],ldir2[STRLEN],fileo[STRLEN];
  char lout[STRLEN];

  buf=getenv("GMXBIN");
  if (buf){
    sscanf(buf,"%s",dircmd);
  }
  else{
    fatal("MUST HAVE THE BINPATH FROM $GMXBIN\n");
  }
  buf=getenv("PWD");
  if (buf){
    sscanf(buf,"%s",direxe);
  }
  else{
    fatal("MUST HAVE THE PATH FROM $PWD\n");
  }

  strcpy(flag,"-f"); strcpy(flag2,"-s"); strcpy(flag3,"-od");

  strcpy(cmd,"g_mindist");
  strcpy(file,"traj.xtc");
  strcpy(file2,"wrk-qm.tpr");
  strcpy(fileo,"mindist.xvg");
  strcat(dircmd,"/");strcat(direxe,"/");
  strcat(dircmd,cmd);
  sprintf(selection,"1\n 2\n");  

  lsum=l0;
  while(lsum<=l1){
    sprintf(ldir,"l%4.2f/%s",lsum,file);
    strcpy(de,direxe);
    strcat(de,ldir);

    sprintf(ldir2,"l%4.2f/%s",lsum,file2);
    strcpy(de2,direxe);
    strcat(de2,ldir2);

    sprintf(lout,"l%4.2f/%s",lsum,fileo);

    strcpy(dout,direxe);
    strcat(dout,lout);

    /* Setup communication pipeline first */
    if( pipe(commpipe) | pipe(stdcopy)){
      fprintf(stderr,"Pipe error!\n");
      exit(1);
    }
    dup2(STDOUT_FILENO, stdcopy[0]);
    dup2(STDIN_FILENO, stdcopy[1]);

    /* Attempt to fork and check for errors */
    if( (pid=fork()) == -1){
      fatal("Fork error. Exiting.\n");  /* something went wrong */
    }
   
    if(pid){
                                              /* A positive (non-negative) PID indicates the parent process */
      dup2(commpipe[1],STDOUT_FILENO);        /* Replace stdout with out side of the pipe */
      close(commpipe[0]);                     /* Close unused side of pipe (in side)      */
      setvbuf(stdout,(char*)NULL,_IONBF,0);   /* Set non-buffered output on stdout        */
      printf("%s\n",selection);               /*		EXECUTE		          */
      wait(&rv);			      /* Wait for child process to end            */
    }
    else{
                                            /* A zero PID indicates that this is the child process */
      dup2(commpipe[0],STDIN_FILENO);	    /* Replace stdin with the in side of the pipe */
      close(commpipe[1]);                   /* Close unused side of pipe (out side)       */
                                            /* Replace the child fork with a new process  */
      if(execl(dircmd,cmd,flag,de,flag2,de2,flag3,dout,NULL) == -1){ 
        fatal("execl Error!\n");
      }
    }
    dup2(stdcopy[0],STDOUT_FILENO);
    dup2(stdcopy[1],STDIN_FILENO);
    close(stdcopy[0]); close(stdcopy[1]);

    lsum+=dl;
  }

  return(0);

}
