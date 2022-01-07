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
    double               dQx;
    double               dQy;
    double               dQz;
    ivec                  nQ;
    rvec                 max;
    rvec                 min;
    real                count;
    gmx_ana_selection_t *refsel;
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
    rvec                r_ik,*x,*dr,dQ;
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
    // BELOW SHOULD BE OPTIMIZED FOR SPECIFIC HARDWARE

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

    d->dQx += dQ[XX];
    d->dQy += dQ[YY];
    d->dQz += dQ[ZZ];
    d->count  += 1.0;
  
    dQ[XX] = (d->max[XX])/(d->nQ[XX]+1)>dQ[XX]?(d->max[XX])/(d->nQ[XX]+1):dQ[XX];
    dQ[YY] = (d->max[YY])/(d->nQ[YY]+1)>dQ[YY]?(d->max[YY])/(d->nQ[YY]+1):dQ[YY];
    dQ[ZZ] = (d->max[ZZ])/(d->nQ[ZZ]+1)>dQ[ZZ]?(d->max[ZZ])/(d->nQ[ZZ]+1):dQ[ZZ];

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
          Ss[sk0] += (rhoc*rhoc+rhos*rhos)/nr/nr;
        }
      }
    }

/* projection code
    if(0){ 
      for(si = 0; si<nr; si++){
        for(sk = si+1; sk<nr; sk++){
          ret_code[0] = pbc_dx_aiuc(pbc,x[si],x[sk],r_ik);              //any box
          pi=nearbyint(abs(r_ik[XX]+d->max[XX]*0.5)/d->max[XX]*nQ[XX]); //projected onto cubic grid
          pj=nearbyint(abs(r_ik[YY]+d->max[YY]*0.5)/d->max[YY]*nQ[YY]);
          pk=nearbyint(abs(r_ik[ZZ]+d->max[ZZ]*0.5)/d->max[ZZ]*nQ[ZZ]);
          sk0=(pi*nQ[YY]+pj)*nQ[ZZ]+pk;
          S[sk0<DY?sk0:DY]+=1; //sk0<DY?sk0:DY
          NX++;
        }
      }
    }
*/
    return 0;
}

/*! \brief
 * Function that implements the analysis tool.
 *
 * Following the style of Gromacs analysis tools, this function is called
 * \p gmx_something.
 */
int
gmx_template(int argc, char *argv[])
{
    const char         *desc[] = {
        "analysis groups.",
    };

    /* Command-line arguments */
    real                cutoff = 0;
    gmx_bool                bArg   = FALSE;
    t_pargs             pa[] = {
        {"-cutoff", FALSE, etREAL, {&cutoff},
         "Cutoff for distance calculation (0 = no cutoff)"},
        {"-arg2",   FALSE, etBOOL, {&bArg},
         "Example argument 2"},
    };
    /* The second argument is for demonstration purposes only */

    /* Output files */
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
    real                  *Sr,dqR;

    CopyRight(stderr, argv[0]);
    /* Here, we can use flags to specify requirements for the selections and/or
     * other features of the library. */

    gmx_ana_traj_create(&trj, ANA_REQUIRE_TOP);
    gmx_ana_set_nrefgrps(trj, 1);
    gmx_ana_set_nanagrps(trj,-1);

    /* If required, other functions can also be used to configure the library
     * before calling parse_trjana_args(). */
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
     * group. */
    rc = gmx_ana_nbsearch_create(&d.nb, cutoff, d.refsel->p.nr);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "neighborhood search initialization failed");
    }

    if(d.fp==NULL)
      d.fp=stderr;

    // GRID ON WHICH TO CALCULATE S(Q)
    dims[XX]  = 50; dims[YY]  = 50; dims[ZZ]  = 50;
    d.max[XX] = 30; d.max[YY] = 30; d.max[ZZ] = 30;
 
    nR = nearbyint( dims[XX]>dims[YY]?(dims[XX]>dims[ZZ]?dims[XX]:dims[ZZ]):(dims[YY]>dims[ZZ]?dims[YY]:dims[ZZ]));
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

    if (1) {
     fprintf(d.fp, "\n");

     dq[XX]= (d.max[XX])/(d.nQ[XX]+1);
     dq[YY]= (d.max[YY])/(d.nQ[YY]+1);
     dq[ZZ]= (d.max[ZZ])/(d.nQ[ZZ]+1);

     dqR   = dq[XX]>dq[YY]?(dq[XX]>dq[ZZ]?dq[XX]:dq[ZZ]):(dq[YY]>dq[ZZ]?dq[YY]:dq[ZZ]); //small incremental stepsize
     for(i=0;i<d.nQ[XX];i++)
       for(j=0;j<d.nQ[YY];j++)
         for(k=0;k<d.nQ[ZZ];k++){
           I=( i*d.nQ[YY]+j )*d.nQ[ZZ]+k;
           dv[0] = sqrt((i+1)*(i+1)*dq[XX]*dq[XX]+(j+1)*(j+1)*dq[YY]*dq[YY]+(k+1)*(k+1)*dq[ZZ]*dq[ZZ]);
           dv[1] = d.Ss[I]/d.N[I]*sel[0]->p.nr;
           fprintf(d.fp,"% 10f % 10.6f\n ",dv[0],dv[1]);
         }
    }
    /* For the template, we close the output file if one was opened */
    if (d.fp) {
      ffclose(d.fp);
    }

    thanx(stderr);
    return 0;
}

/*! \brief
 * The main function.
 */
int
main(int argc, char *argv[])
{
    gmx_template(argc, argv);
    return 0;
}
