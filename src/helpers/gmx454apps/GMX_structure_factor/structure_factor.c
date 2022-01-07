/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
#include <gromacs/copyrite.h>
#include <gromacs/filenm.h>
#include <gromacs/macros.h>
#include <gromacs/pbc.h>
#include <gromacs/smalloc.h>
#include <gromacs/statutil.h>
#include <gromacs/vec.h>
#include <gromacs/xvgr.h>

#include <gromacs/nbsearch.h>
#include <gromacs/trajana.h>
#include <fftw3.h>

/*! \brief
 * Template analysis data structure.
 */
typedef struct
{
    double                 dQx;
    double                 dQy;
    double                 dQz;
    int                  bFlat;
    ivec                    nQ;
    rvec                   max;
    rvec                  qres;
    rvec                   min;
    rvec                    ZP;
    rvec                  mQZP;
    rvec                  mdZP;
    gmx_bool             bZoom,bRes;
    real                  count;
    gmx_ana_selection_t  *refsel;
    FILE                *fp;
    real                *ave;
    real                *n;
    real                *Sc;
    real                *Ss;
    real                *N;
    rvec                *dr;
    gmx_ana_nbsearch_t  *nb;
} t_analysisdata;


static int
analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
    t_analysisdata     *d = (t_analysisdata *)data;
    int                 g, i,j,k;
    int                 si,sj,sk,si0,sj0,sk0,nbx,nby,nbz,DX,DY,NX,NY;
    int                 pi,pj,pk,pi0,pj0,pk0,npx,npy,npz;
    real                frave;
    int                 rc,ret_code[2];
    rvec                r_ik,*x,*dr,dQ,Qzero,dQz;
    real                *Ss,*Sc,*R,qx,qy,qz,omega,*N,Re=0.0,Im=0.0;
    int                 *nQ;
    real               rhoc=0.0,rhos=0.0;
/* 

THIS FUNCTION SUPPLIES FUNCTIONALITY TO CALCULATE THE STRUCTURE FACTOR AT A SPECIFIC TIME
CURRENTLY THE RESULT IS JUST SUMMED UP INTO THE CORRECT S(Q) BIN

----( STATIC STRUCTURE FACTOR )----
NOTES:
current algorithm runs by using block summation

*/
    g=0;

    dQ[XX] = 2*M_PI/fr->box[XX][XX];
    dQ[YY] = 2*M_PI/fr->box[YY][YY];
    dQ[ZZ] = 2*M_PI/fr->box[ZZ][ZZ];

    Sc   = d->Sc;
    Ss   = d->Ss;
    N    = d->N;
    nQ   = d->nQ;
    x    = sel[g]->p.x;
    nr   = sel[g]->p.nr;
 
    rc = gmx_ana_nbsearch_pos_init(d->nb, pbc, &d->refsel->p);
    if (rc != 0){
      gmx_fatal(FARGS, "Neighborhood search initialization failed");
    }

    d->dQx   += dQ[XX];
    d->dQy   += dQ[YY];
    d->dQz   += dQ[ZZ];
    d->count += 1.0;
  
    if(!(d->bZoom)){
      dQ[XX] = (d->max[XX])/(d->nQ[XX])>dQ[XX]?(d->max[XX])/(d->nQ[XX]):dQ[XX]; //2*M_PI/(d->max[XX]);//>dQ[XX]?2*M_PI/(d->max[XX]):dQ[XX];
      dQ[YY] = (d->max[XX])/(d->nQ[YY])>dQ[YY]?(d->max[YY])/(d->nQ[YY]):dQ[YY]; //>dQ[YY]?2*M_PI/(d->max[XX]):dQ[YY];
      dQ[ZZ] = (d->max[XX])/(d->nQ[ZZ])>dQ[ZZ]?(d->max[ZZ])/(d->nQ[ZZ]):dQ[ZZ]; //>dQ[ZZ]?2*M_PI/(d->max[XX]):dQ[ZZ];
    }

    if(d->bZoom){
      if(!d->bRes){
        for(i=XX;i<=ZZ;i++){ 
          Qzero[i] = (d->ZP[ i ]) - dQ[ i ] * ( nearbyint(d->nQ[ i ]/2.0) + 2 );
          dQz[i]   = (dQ[ i ] * ( nearbyint(d->nQ[ i ]/2.0) + 2 ) * 2.0)/((float)d -> nQ[ i ]);
          d -> mQZP[i] += Qzero[i];
          d -> mdZP[i] += dQz[i];
        }  
      }else{
        for(i=XX;i<=ZZ;i++){
          Qzero[i]      = (d->ZP[ i ]) - d->qres[i] * ( nearbyint(d->nQ[ i ]/2.0) );
          dQz[i]        = d->qres[i];  
          d -> mQZP[i] += Qzero[i];
          d -> mdZP[i] += dQz[i];
        }
      }

    }else{
      for(i=XX;i<=ZZ;i++){
        dQz[i]   = dQ[ i ];
        Qzero[i] = 0.0;
      }
    }

    if(d->bFlat){
      for(i=0;i<nQ[XX];i++){
        for(j=0;j<nQ[YY];j++){  
          qx=Qzero[XX]+(i+1)*dQz[XX];
          qy=Qzero[YY]+(j+1)*dQz[YY];
          sk0=( i*nQ[YY]+ j );
          rhoc=0.0; rhos=0.0;
          for(si0=0;si0<nr;si0++){
            omega = x[si0][XX]*qx+x[si0][YY]*qy;
            rhoc += cos(omega);
            rhos += sin(omega);
          }
          N[sk0]  += 1.0;
          Ss[sk0] += (rhoc*rhoc+rhos*rhos)/nr;//nr;
        }
      }    
    }else{
      for(k=0;k<nQ[ZZ];k++){
        for(i=0;i<nQ[XX];i++){
          for(j=0;j<nQ[YY];j++){  
            qx=dQ[XX]*(i+1);
            qy=dQ[YY]*(j+1);
            qz=dQ[ZZ]*(k+1);
            sk0=( i*nQ[YY]+ j )*nQ[ZZ] + k;
            rhoc=0.0; rhos=0.0;
            for(si0=0;si0<nr;si0++){
              omega = x[si0][XX]*qx+x[si0][YY]*qy+x[si0][ZZ]*qz;
              rhoc += cos(omega);
              rhos += sin(omega);
            }
            N[sk0]  += 1.0;
            Ss[sk0] += (rhoc*rhoc+rhos*rhos)/nr;//nr;
          }
        }
      }
    }

    return 0;
}

int
gmx_template(int argc, char *argv[])
{
    const char         *desc[] = {
        "analysis groups.",
    };

    /* Command-line arguments */
    real                cutoff = 0;
    real                BOXS=10;
    real                BINDIM=50;
    rvec                zoom_point,resvect;
    gmx_bool            bArg   = FALSE, bZoom = FALSE, bRes = FALSE;
    t_pargs             pa[] = {
        {"-maxq", FALSE, etREAL, {&BOXS},
         "Maximum Q value"},
        {"-bindim",   FALSE, etREAL, {&BINDIM},
         "Bin dimension"},
        {"-b2d",   FALSE, etBOOL, {&bArg},
         "force pure 2d calculation"},
        {"-bZoom",   FALSE, etBOOL, {&bZoom},
         "will try to zoom in on a specific region of S(Q)"},
        {"-zoom_point",FALSE, etRVEC, {zoom_point},
         "will try to zoom in on this specific point with the best resolution possible"},
        {"-bRes",   FALSE, etBOOL, {&bRes},
         "Manually override resolution"},
        {"-resvect",FALSE, etRVEC, {resvect},
         "Assign the resolution"},

    };

    t_filenm            fnm[] = {
        {efXVG, "-o", "avedist", ffOPTWR},
    };
#define NFILE asize(fnm)

    gmx_ana_traj_t        *trj;
    output_env_t          oenv;
    t_analysisdata        d;
    int                   ngrps;
    gmx_ana_selection_t   **sel;
    int                   g;
    int                   rc;
    int                   pi,pj,pk,I,i,j,k,iqr,nR=250;
    rvec                  dq,dims,dv;
    real                  *Sr,dqR,Area;
    FILE                  *fout;

    zoom_point[XX]= 0.0;
    zoom_point[YY]= 0.0;
    zoom_point[ZZ]= 0.0;
    resvect[XX]   = 0.0;
    resvect[YY]   = 0.0;
    resvect[ZZ]   = 0.0;

    CopyRight(stderr, argv[0]);

    /*
     * Here, we can use flags to specify requirements for the selections and/or
     * other features of the library.
     */

    gmx_ana_traj_create(&trj, ANA_REQUIRE_TOP);
    gmx_ana_set_nrefgrps(trj, 1);
    gmx_ana_set_nanagrps(trj,-1);

    /*
     * If required, other functions can also be used to configure the library
     * before calling parse_trjana_args().
     */

    parse_trjana_args(trj, &argc, argv, PCA_CAN_VIEW,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL,
                      &oenv);

    /* First, we get some selection information from the structure */
    gmx_ana_get_refsel(trj, 0, &d.refsel);
    gmx_ana_get_nanagrps(trj, &ngrps);
    gmx_ana_get_anagrps(trj, &sel);

    if (ngrps != 1)
    {
        gmx_fatal(FARGS, "only one group permitted");
    }

    /* First, we initialize the neighborhood search for the first index
     * group.
     */
    rc = gmx_ana_nbsearch_create(&d.nb, cutoff, d.refsel->p.nr);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "neighborhood search initialization failed");
    }

    if(d.fp==NULL)
      d.fp=stderr;

    // GRID ON WHICH TO CALCULATE S(Q)
    dims[XX]  = BINDIM; dims[YY]  = BINDIM; dims[ZZ]  = BINDIM;
    d.max[XX] = BOXS;   d.max[YY] = BOXS;   d.max[ZZ] = BOXS;

    d.bZoom=bZoom;
    if(bZoom){
      for(i=XX;i<=ZZ;i++){
        d.ZP[i]=zoom_point[i];
        d.mQZP[i]=0.0;
        d.mdZP[i]=0.0;
      }
    }
    d.bRes=bRes;
    if(bRes){
      for(i=XX;i<=ZZ;i++){
        d.qres[i]=resvect[i];      
      }
    }

    d.bFlat=1;

    if(bArg){
      nR = BINDIM;
      dims[ZZ]=1;
    }else{
      nR = nearbyint( dims[XX]>dims[YY]?(dims[XX]>dims[ZZ]?dims[XX]:dims[ZZ]):(dims[YY]>dims[ZZ]?dims[YY]:dims[ZZ]));
    }

    snew(d.Ss,dims[XX]*dims[YY]*dims[ZZ]);
    snew(d.Sc,dims[XX]*dims[YY]*dims[ZZ]);
    snew(d.N ,dims[XX]*dims[YY]*dims[ZZ]);
    snew(Sr,nR*4);
    d.nQ[XX]  = dims[XX]; d.nQ[YY]  = dims[YY]; d.nQ[ZZ]  = dims[ZZ];
    d.count=0;

    /* We also open the output file if the user provided it */
    d.fp = NULL;
    if (opt2bSet("-o", NFILE, fnm))
    {
        d.fp = xvgropen(opt2fn("-o", NFILE, fnm), "Structure factor",
                        "Q [ nm -1 ]","S(Q)", oenv);
        xvgr_selections(d.fp, trj);
    }
    gmx_ana_do(trj, 0, &analyze_frame, &d);
    dq[XX]= (d.max[XX])/(d.nQ[XX]);
    dq[YY]= (d.max[YY])/(d.nQ[YY]);
    dq[ZZ]= (d.max[ZZ])/(d.nQ[ZZ]);

    if (d.fp) {
      fprintf(d.fp, "\n");

      if(d.bFlat){
        dqR   = dq[XX]>dq[YY]?dq[XX]:dq[YY]; 				
      }else{
        dqR   = dq[XX]>dq[YY]?(dq[XX]>dq[ZZ]?dq[XX]:dq[ZZ]):(dq[YY]>dq[ZZ]?dq[YY]:dq[ZZ]); 
      }
      for(i=0;i<d.nQ[XX];i++){
        for(j=0;j<d.nQ[YY];j++){
          if(d.bFlat){
            I= i*d.nQ[YY]+j ;
            dv[0] = sqrt((i+1)*(i+1)*dq[XX]*dq[XX]+(j+1)*(j+1)*dq[YY]*dq[YY]);
            dv[1] = d.Ss[I]/d.N[I]*sel[0]->p.nr;
            fprintf(d.fp,"% 10f % 10.6f\n ",dv[0],dv[1]);
          }else{
            for(k=0;k<d.nQ[ZZ];k++){
              I=( i*d.nQ[YY]+j )*d.nQ[ZZ]+k;
              dv[0] = sqrt((i+1)*(i+1)*dq[XX]*dq[XX]+(j+1)*(j+1)*dq[YY]*dq[YY]+(k+1)*(k+1)*dq[ZZ]*dq[ZZ]);
              dv[1] = d.Ss[I]/d.N[I]*sel[0]->p.nr;
              fprintf(d.fp,"% 10f % 10.6f\n ",dv[0],dv[1]);
            }
          }
        }
      }
      ffclose(d.fp);    
    }
    
    fprintf(stderr,"< dq_x > = %10.6f  [ nm^{-1} ]   |  max Qx = %10.6f  [ nm^{-1} ] \n",d.dQx/d.count,dq[XX]);
    fprintf(stderr,"< dq_y > = %10.6f  [ nm^{-1} ]   |  max Qy = %10.6f  [ nm^{-1} ] \n",d.dQy/d.count,dq[YY]);  

    if(bZoom){
      if(!bRes){
        d.qres[XX]=d.dQx/d.count;
        d.ZP[XX]=0.0;
      }
      fout=fopen("inp.dat","w"); 
      fprintf(fout," %10.6f %10.6f ",d.qres[XX],d.ZP[XX]);
      fprintf(fout," %10.6f %10.6f \n",d.qres[YY],d.ZP[YY]);
      fclose(fout);
    }
    if(!d.bFlat){
      fprintf(stderr,"< dq_z > = %10.6f  [ nm^{-1} ]\n",d.dQz/d.count);
    }

    fout=fopen("sq.dat","w"); 
    k=0;

    for(i=0;i<d.nQ[XX];i++){
      for(j=0;j<d.nQ[YY];j++){
        I=( i*d.nQ[YY]+j );
	fprintf(fout,"%10.6f ", d.Ss[I]/d.N[I]);
      }
      fprintf(fout,"\n");
    }       
    fclose(fout);

    thanx(stderr);
    return 0;
}

int
main(int argc, char *argv[])
{
    gmx_template(argc, argv);
    return 0;
}
