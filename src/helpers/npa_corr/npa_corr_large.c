#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "confio.h"
#include "gmx_fatal.h"
#include "futil.h"
#include "gstat.h"
#include "macros.h"
#include "maths.h"
#include "physics.h"
#include "index.h"
#include "string.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "strdb.h"
#include "xvgr.h"

#define SLEN 1000
#define XX 0
#define YY 1
#define ZZ 2
#define MAX 20000
#define DIM 3


void fatal(char str[])
{
  fprintf(stderr,"\nERROR:%s\n",str);
  exit(1);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "this is a small program meant to calculate dipole dependent properties",
  };
  static int n=1;
  static int tick=0;
  static real deltat=0.001;
  t_pargs pa[] = {
    { "-n", FALSE, etINT, {&n},
      "Plot data for atom number n (starting on 1)"
    },
    { "-tick", FALSE, etINT, {&tick},
      "one residue conists of this many atoms (default is off)"
    },
    { "-deltat", FALSE, etREAL, {&deltat},
      "the timestep used. default i 0.001"
    }
  };

  t_topology top;
  int        ePBC;

  int a;
  int i,j,k,ticker=1,jid;
  int u,v,w,umax=0,ustop=100,v3,n_alloc=0;
  int resid,atoid,oresid;
  rvec pos,MU,mu;
  real **store;
  char line[SLEN],typ[2];
  double charge;
  FILE *fin,*fout;

  t_filenm fnm[] = {
    { efXVG, "-o", "PAC",  ffWRITE }
  };
  
#define NFILE asize(fnm)
  int     npargs;
  t_pargs *ppa;

  npargs = asize(pa);
  ppa    = add_acf_pargs(&npargs,pa);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL);
  fprintf(stderr,"MAIN\n");

  snew(store,DIM*MAX);
  for(u=0;u<ustop;u++)
    srenew(store[u],DIM*MAX);
  n_alloc=DIM*MAX;

  fin=fopen("npa.log","r");
  if(fin==NULL)
    fatal("COULD NOT OPEN npa.log");
  fout=fopen("mt.log","w");
  if(fout==NULL)
    fatal("COULD NOT OPEN mt.log");

//INITS
  w=0;resid=0;oresid=0;v=0;u=0;
  for(i=XX;i<=ZZ;i++){
    mu[i]=0.0;
    pos[i]=0.0;
    MU[i]=0.0;
  }
  fprintf(stderr,"FINISHED INIT\n");

  ticker=0;u=1;umax=0;
  while(!feof(fin) && v<MAX){

    if(v >= n_alloc) {
      n_alloc+=100;
      for(i=0; i<ustop; i++){
        srenew(store[i],DIM*n_alloc);
     }
   }
   v3=3*v;

   if( ticker<tick ){
     if( NULL == fgets(line,SLEN,fin) ) { 
        ;
      }
      if(line[0]=='@')
        continue;
      sscanf(line,"%d %f %f %f %d %d %s %lf ",&resid,&pos[XX],&pos[YY],&pos[ZZ],&atoid,&jid,typ,&charge);
      ticker++; u=resid+1;
      umax=umax<u?u:umax;
      if(umax>ustop)
        fatal("EXCEEDED AMOUNT OF RESIDUES");
      for(i=XX;i<=ZZ;i++){
        mu[i]+=pos[i]*charge;
        store[(u-1)][(v)*3+i]=mu[i];
        MU[i]+=mu[i];
      }
      if(norm(mu)==0.0)
        fprintf(stderr,"\nWarning\n");
    }
    else{
      v++; ticker=0;
      fprintf(fout,"%f  %f %f %f\n",v*deltat,MU[XX],MU[YY],MU[ZZ]);
      for(i=XX;i<=ZZ;i++){
        MU[i]=0.0; 
        mu[i]=0.0;
      }
    }
  }
  v--;
  fprintf(stdout,"\nv - %d           um - %d         cgs - %d\n",v,umax,top.cgs.nr);
  
  fprintf(stderr,"DO ACORR\n");
  do_autocorr(ftp2fn(efXVG,NFILE,fnm),
     "Dipole-Dipole Autocorrelation Function",
     v,umax,store,deltat,eacVector,TRUE);  

/*
dt=0.001;
n=length(P);
dw=2*pi/n/dt;
om=(0:dw:(n-1)*dw);
F=real(fft(P(:,2)))'.*tanh((0:dw:(n-1)*dw)*pi*0.01);
plot(om(1:n/2),F(1:n/2)),axis([0 2500 0 2])
*/
  return(0);
}
