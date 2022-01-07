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

/*! \brief
 * Template analysis data structure.
 */
typedef struct
{
    gmx_ana_selection_t  *refsel;
    FILE                 *fp;
    real                 *ave;
    real                 *n;
    real                 *R;
    int                   D;
    real                 idr;
    real                max;
    real                tol;
    rvec               BOXL;
    real              count;
    real               *svol;
    real                num;
    real               ivol;
    gmx_ana_nbsearch_t  *nb;
} t_analysisdata;

static int
analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
    t_analysisdata      *d = (t_analysisdata *)data;
    int                 g,i,j,k;
    real                frave,dz,r,theta=0,omega;
    int                 rc,skip,si,sj,ti,tj,ret_code;
    ivec                pos[2];
    rvec                *x,r_ik,r_ki;
    real                *svol,prev_vol,pres_vol,segvol,iphi,irad;

    rc = gmx_ana_nbsearch_pos_init(d->nb, pbc, &d->refsel->p);

    if (rc != 0)
    {
        gmx_fatal(FARGS, "Neighborhood search initialization failed");
    }
    g=0;
    frave=0;

    x    = sel[g]->p.x;
    nr   = sel[g]->p.nr;
    d->num=nr;
    g=d->D;

    d->BOXL[XX] += fr->box[XX][XX];
    d->BOXL[YY] += fr->box[YY][YY];
    d->BOXL[ZZ] += fr->box[ZZ][ZZ];
    d->count+=1.0;
    iphi = (d->D)/2/M_PI;
    irad = d->idr;

    d->R[0] = 0.0;
    if(1){ 
      for(i = 0; i<nr; i++){
        for(k = i+1; k<nr; k++){
          ret_code = pbc_dx_aiuc(pbc,x[i],x[k],r_ik);
          r = sqrt(r_ik[XX]*r_ik[XX]+r_ik[YY]*r_ik[YY]); 
          omega = r_ik[YY]/r_ik[XX];
          theta = atan(omega);

          if(r_ik[XX]< 0 && r_ik[YY]>=0)
            theta+=M_PI;
          if(r_ik[XX]< 0 && r_ik[YY]<0)
            theta-=M_PI;
          if(r_ik[XX]==0 && r_ik[YY]>0)
            theta=M_PI*0.5;
          if(r_ik[XX]==0 && r_ik[YY]<0)
            theta=-0.5*M_PI;
          if(r_ik[XX]==0 && r_ik[YY]==0)
            theta=0.0;

          theta+=M_PI;
          omega=(theta+M_PI)>2.0*M_PI?theta-M_PI:theta+M_PI;

	  // SUM UP 2D
          si=nearbyint( r*irad );
          sj=nearbyint( theta*iphi );
          if(si>=0 && si<=g && sj>=0 && sj<=g ){
            d->R[(g*g*5+si*g + sj  )]+=1.0;
            d->R[0]+=1.0; 
          }

          ti=si;
          tj=nearbyint( omega*iphi );
          if(ti>=0 && ti<=g && tj>=0 && tj<=g ){
            d->R[g*g*5+ti*g + tj ]+=1.0;
            d->R[0]+=1.0;
          }
        }
      }
    }
    for(i=0; i<g; i++) {
       d->R[g*g*5+g*i+g-1]+=d->R[g*g*5+i*g];
       d->R[g*g*5+g*i+g-1]*=0.5;
       d->R[g*g*5+g*i]=d->R[g*g*5+i*g+g-1];
    }
    prev_vol=0.0;
    for(i=0; i<g; i++) {
      r          = (i+0.5)*((float)d->max)*0.5/((float)g);
      pres_vol   = M_PI*r*r;
      segvol     = pres_vol-prev_vol;
      d->svol[i] = 1.0/segvol;
      prev_vol   = pres_vol;
      for(j=0; j<g; j++){
        d->R[g*g*5 + i*g + j]*=(d->svol[i]*( (float)g ));
      }
    }
    theta=0.0; si = g-2;
    for(i=0;i<g;i++){ 
      d->R[ g*g*5 + i ]=0.0;
      theta+=d->R[g*g*5 + si*g + i  ];
    }
    omega=(float)g/theta;
    for(i=0; i<g; i++){
      for(j=0;j<g;j++){
        d->R[g*g*6 + i*g + j]+= d->R[g*g*5 + i*g + j]*omega;
        d->R[g*g*5 + i*g + j] = 0.0;
        d->R[g*g*6 + j] = 0.0;
      }
    }
    /* We need to return 0 to tell that everything went OK */
    return 0;
}

int
gmx_template(int argc, char *argv[])
{
    const char         *desc[] = {
        "This is a template for writing your own analysis tools for",
        "analysis groups.",
    };

    real                BOXS = 10.0;
    real                BINDIM = 100;
    gmx_bool                bArg   = FALSE;
    t_pargs             pa[] = {
        {"-boxside", FALSE, etREAL, {&BOXS},
         "box side length"},
        {"-bindim",  FALSE, etREAL, {&BINDIM},
         "Amount of bins"},
    };

    t_filenm            fnm[] = {
        {efXVG, "-o", "avedist", ffOPTWR},
    };
#define NFILE asize(fnm)

    gmx_ana_traj_t        *trj;
    output_env_t          oenv;
    t_analysisdata        d;
    int                   ngrps;
    real                  *R,*svol,prev_vol,pres_vol;
    gmx_ana_selection_t   **sel;
    int                   g;
    int                   rc;
    int                   D;
    FILE                  *fout;
    int                   i,j;
    real                  r,segvol;
    real                  cutoff=0.0;

    CopyRight(stderr, argv[0]);

    gmx_ana_traj_create(&trj, ANA_REQUIRE_TOP);
    gmx_ana_set_nrefgrps(trj, 1);
    gmx_ana_set_nanagrps(trj, -1);

    parse_trjana_args(trj, &argc, argv, PCA_CAN_VIEW,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL,
                      &oenv);

    gmx_ana_get_refsel(trj, 0, &d.refsel);
    gmx_ana_get_nanagrps(trj, &ngrps);
    gmx_ana_get_anagrps(trj, &sel);

    rc = gmx_ana_nbsearch_create(&d.nb, cutoff, d.refsel->p.nr);
    if (rc != 0 || ngrps>1)
    {
        gmx_fatal(FARGS, "neighborhood search initialization failed");
    }

    snew(d.ave, ngrps);
    snew(d.n,   ngrps);
    d.D=(int)BINDIM; d.max=BOXS;
    d.idr=((float)d.D)/(((float)d.max))*2.0; d.tol=d.idr*0.5; //defaults

    snew(d.R,d.D*d.D*10);
    snew(d.svol,d.D);

    gmx_ana_do(trj, 0, &analyze_frame, &d);

    fout=fopen("g2d.dat","w");
    d.R[(int)((d.D)*(d.D)*6.0 )] = d.max*0.5*d.count;
    for(i=0;i<d.D;i++){
      for(j=0;j<d.D;j++){
        fprintf(fout," % 6.3f",d.R[(int)((d.D)*(d.D)*6.0 + i*d.D + j)]/d.count);
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
