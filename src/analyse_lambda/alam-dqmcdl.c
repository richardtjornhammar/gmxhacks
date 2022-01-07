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

  float dl=0.02,l0=0.00,l1=1.00,lsum=0.0;  
  int i,N;

  pid_t pid;
  int   rv,step_id=-1,T=1000;
  int   commpipe[2],stdcopy[2];            /* This holds the fd for the in & out of the pipe */
  FILE  *ing,*inq,*out;
  char buf[STRLEN];
  char buf1[STRLEN],buf2[STRLEN];
  char file1[STRLEN],file2[STRLEN],ldir1[STRLEN],ldir2[STRLEN];
  char lout[STRLEN];
  float qm,dg,t;

  sprintf(file1,"dgdl.xvg");
  sprintf(file2,"qmener.dat");

  lsum=l0;
  while(lsum<=l1){
    sprintf(ldir1,"l%4.2f/%s",lsum,file1);
    sprintf(ldir2,"l%4.2f/%s",lsum,file2);
    sprintf(lout,"l%4.2f/%s",lsum,"dqmcdl.dat");

    ing=fopen(ldir1,"r");
    inq=fopen(ldir2,"r");
    out=fopen(lout,"w");

    while((!feof(ing) || !feof(inq)) || step_id<T) {
      
      if( NULL == fgets(buf1,STRLEN,ing) && step_id<T) {
	fatal(" Error handling ");
      }
      if(buf1[0]=='#'||buf1[0]=='@')
	while( !feof(ing) && (buf1[0]=='#'||buf1[0]=='@') ){
	  if( NULL == fgets(buf1,STRLEN,ing)  && step_id<T) {
	    fatal(" Error handling ");
	  }  
	}
      if( NULL == fgets(buf2,STRLEN,inq) && step_id<T) {               
        fatal(" Error handling ");
      }
      if(buf2[0]=='#'||buf2[0]=='@')
	while( !feof(inq) && (buf2[0]=='#'||buf2[0]=='@') ){
	  if( NULL == fgets(buf2,STRLEN,inq)  && step_id<T) {
	    fatal(" Error handling ");
	  }  
	}
      sscanf(buf1,"%f %f",&t,&dg);
      sscanf(buf2,"%f",&qm);      
      if(step_id<T)
	fprintf(out,"%f %f %f %f \n",t,dg+qm,dg,qm);
      step_id++;
    }

    fclose(out);
    fclose(ing);
    fclose(inq);
    step_id=-1;
    lsum+=dl;
  }

  return(0);

}
