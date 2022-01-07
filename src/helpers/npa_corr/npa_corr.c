#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define SLEN 1000
#define XX 0
#define YY 1
#define ZZ 2
#define MAX 100000

typedef double rvec[3];

void fatal(char str[])
{
  fprintf(stderr,"\nERROR:%s\n",str);
  exit(1);
}

int main(void)
{
  int i,j,k;
  int u,v,w;
  int resid,atoid,oresid;
  rvec pos,MU,mu;
  rvec store[MAX];
  char line[SLEN],typ[2];
  double charge;
  FILE *fin,*fout;

  fprintf(stderr,"MAIN\n");

  fin=fopen("npa.log","r");
  if(fin==NULL)
    fatal("COULD NOT OPEN npa.log");
  fout=fopen("mt.log","w");
  if(fout==NULL)
    fatal("COULD NOT OPEN mt.log");
  w=0;resid=0;oresid=0;v=0;u=0;

  fprintf(stderr,"FINISHED INIT\n");

  for(i=XX;i<=ZZ;i++){
    mu[i]=0.0;
    pos[i]=0.0;
    MU[i]=0.0;
  }
  while(!feof(fin) && w<MAX){
    if( NULL == fgets(line,SLEN,fin) ) {    
      ;
    }
    sscanf(line,"%d %lf %lf %lf %d %s %lf ",&resid,&pos[XX],&pos[YY],&pos[ZZ],&atoid,typ,&charge);
    if(resid==oresid){
      for(i=XX;i<=ZZ;i++)
        mu[i]+=pos[i]*charge;
    }
    else{
      u++;
      for(i=XX;i<=ZZ;i++){
        MU[i]+=mu[i];
        mu[i]=pos[i]*charge;
      }
      if(!(oresid==0) && resid==0){
        fprintf(fout,"%d  %f %f %f\n",v,MU[XX],MU[YY],MU[ZZ]);
        for(i=XX;i<=ZZ;i++){
          store[v][i]=MU[i]/u;
          MU[i]=0.0;
        }
        u=0;v++;
      }
      oresid=resid;
    }
    w++;
  }

  return(0);
}
