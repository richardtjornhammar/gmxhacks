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

int do_whatever(char selection[],char dircmd[],char cmd[],char opt0[],char opt1[],char opt2[],char opt3[],char opt4[],char opt5[],char opt6[]){
  pid_t pid;
  int   rv,step_id;
  int   commpipe[2],stdcopy[2];            /* This holds the fd for the in & out of the pipe */

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
      if(execl(dircmd,cmd,opt0,opt1,opt2,opt3,opt4,opt5,NULL) == -1){ 
        fatal("execl Error!\n");
      }
    }
    dup2(stdcopy[0],STDOUT_FILENO);
    dup2(stdcopy[1],STDIN_FILENO);
    close(stdcopy[0]); close(stdcopy[1]);

    return(0);
}

int main(void){

  float dl=0.02,l0=0.00,l1=1.00,lsum=0.0;  
  int i,N;

  pid_t pid;
  int   rv,step_id;
  int   commpipe[2],stdcopy[2];            /* This holds the fd for the in & out of the pipe */
  FILE  *in,*out;
  char  *buf;
  char cmd[STRLEN],dircmd[STRLEN],direxe[STRLEN],flag[STRLEN],flag2[STRLEN];
  char selection[STRLEN],file[STRLEN],ldir[STRLEN],de[STRLEN],de2[STRLEN],dout[STRLEN];
  char file2[STRLEN],flag3[STRLEN],ldir2[STRLEN],fileo[STRLEN];
  char q1[10],q2[10],q3[10];
  char lout[STRLEN];

  buf=getenv("PWD");
  if (buf){
    sscanf(buf,"%s",direxe);
  }
  else{
    fatal("MUST HAVE THE PATH FROM $PWD\n");
  }
  sprintf(selection," ");  
  strcpy(de2,direxe);
  strcat(de2,"lamdqmcdl");
  do_whatever(selection,"lamdqmcdl",NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);

  strcpy(de2,direxe);
  strcat(de2,"binit"); 

  sprintf(flag,"-i");
  sprintf(file,"dqmcdl.dat");
  
  sprintf(q1,"-b");
  sprintf(q2,"20");
  sprintf(q3,"2");

  sprintf(flag2,"--output=");
  sprintf(fileo,"dqmcdl.hist");

  lsum=l0;
  while(lsum<=l1){
    sprintf(flag2,"--output=");
    sprintf(fileo,"dqmcdl.hist");
    sprintf(selection," "); sprintf(ldir," "); sprintf(de," ");
    sprintf(de2," "); sprintf(dout," ");

    sprintf(ldir,"/l%4.2f/%s",lsum,file);
    strcpy(de,direxe);
    strcat(de,ldir); //input file

    //./binit -i dqmcdl.dat -b 20 2 > dqmcdl.hist

    sprintf(lout,"/l%4.2f/%s",lsum,fileo);
    strcpy(dout,direxe);
    strcat(dout,lout);  //output file
    strcat(flag2,dout);

    //do_whatever(selection,"/usr/bin/touch","touch",de2,NULL,NULL,NULL,NULL,NULL,NULL);
    //printf("%s %s %s %s %s %s %s %s %s\n",selection,"./binit","binit",flag,de,q1,q2,q3,flag2);
    do_whatever(selection,"./binit","binit",flag,de,q1,q2,q3,flag2,NULL);

    sprintf(fileo,"./l%4.2f/.",lsum);
    //printf("\n\n%s\n\n",fileo);
    do_whatever(selection,"/bin/cp","/bin/cp","./gauB_qmc.m",fileo,NULL,NULL,NULL,NULL,NULL);

    sprintf(selection,"cd l%4.2f\n gauB_qmc \n cd ..\n exit \n",lsum);
    do_whatever(selection,"/usr/bin/octave","octave",NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    sprintf(selection," ");  

    lsum+=dl;
  }

  return(0);

}
