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
    real                 dr;
    real                max;
    real                tol;
    rvec               BOXL;
    real              count;
    gmx_ana_nbsearch_t  *nb;
} t_analysisdata;

/*! \brief
 * Function that does the analysis for a single frame.
 *
 * It is called once for each frame.
 */
static int
analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
    t_analysisdata      *d = (t_analysisdata *)data;
    int                 g,i,j,k;
    real                frave,dz;
    int                 rc,skip;
    ivec                pos[2];
    rvec                *x;

    /* Here, you can do whatever analysis your program requires for a frame. */
    if (d->fp)
    {
        fprintf(d->fp, "%10.3f", fr->time);
    }
    rc = gmx_ana_nbsearch_pos_init(d->nb, pbc, &d->refsel->p);

    if (rc != 0)
    {
        gmx_fatal(FARGS, "Neighborhood search initialization failed");
    }
    g=0;
    frave=0;

    x    = sel[g]->p.x;
    nr   = sel[g]->p.nr;

    d->BOXL[XX] += fr->box[XX][XX];
    d->BOXL[YY] += fr->box[YY][YY];
    d->BOXL[ZZ] += fr->box[ZZ][ZZ];
    d->count+=1.0;

    for(i=0;i<nr;i+=2){
      skip=0;
      for(j=XX;j<=YY;j++){ //doing 2D
        pos[0][j]=nearbyint(x[i][j]/d->dr);
        pos[1][j]=nearbyint(x[i+1][j]/d->dr);
        if( !(pos[0][j]>0 && pos[0][j]<d->D) || !(pos[1][j]>0 && pos[1][j]<d->D) )
          skip++;
      }
      if(skip){ //must deduct skipped nr
        ;
      }
      dz=x[i][ZZ]-x[i+1][ZZ];
      if( dz < d->tol ){
        frave += 1.0;
      }
      else{
        if(dz>0 && pos[0][XX]>0 && pos[0][YY]>0){
          d->R[ (pos[0][XX]*(d->D) + pos[0][YY]) ]+=1.0;
        }
        if(dz<0 && pos[0][XX]>0 && pos[0][YY]>0){
          d->R[ (pos[0][XX]*(d->D) + pos[0][YY]) + d->D*d->D ]+=1.0;
        }
      }
    }

/* projection code
    if(0){ 
      for(si = 0; si<nr; si++){
        for(sk = si+1; sk<nr; sk++){
          ret_code[0] = pbc_dx_aiuc(pbc,x[si],x[sk],r_ik);              //any box
          pi=nearbyint(abs(r_ik[XX]+d->max[XX]*0.5)/d->max[XX]*nQ[XX]); //projected onto 3 square grids
          pj=nearbyint(abs(r_ik[YY]+d->max[YY]*0.5)/d->max[YY]*nQ[YY]);
          pk=nearbyint(abs(r_ik[ZZ]+d->max[ZZ]*0.5)/d->max[ZZ]*nQ[ZZ]);
          sk0=(pi*nQ[YY]+pj)*nQ[ZZ]+pk;
          S[sk0<DY?sk0:DY]+=1; //sk0<DY?sk0:DY
          NX++;
        }
      }
    }
*/

    if (d->fp)
    {
        fprintf(d->fp, " %.3f", frave);
        fprintf(d->fp, "\n");
    }

    /* We need to return 0 to tell that everything went OK */
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
        "This is a template for writing your own analysis tools for",
        "Gromacs. The advantage of using Gromacs for this is that you",
        "have access to all information in the topology, and your",
        "program will be able to handle all types of coordinates and",
        "trajectory files supported by Gromacs. In addition,",
        "you get a lot of functionality for free from the trajectory",
        "analysis library, including support for flexible dynamic",
        "selections. Go ahead an try it![PAR]",
        "To get started with implementing your own analysis program,",
        "follow the instructions in the README file provided.",
        "This template implements a simple analysis programs that calculates",
        "average distances from the a reference group to one or more",
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
    real                  *R;
    gmx_ana_selection_t   **sel;
    int                   g;
    int                   rc;
    int                   D;
    FILE                  *fout;
    int                   i,j;

    CopyRight(stderr, argv[0]);
    /* Here, we can use flags to specify requirements for the selections and/or
     * other features of the library. */
    gmx_ana_traj_create(&trj, ANA_REQUIRE_TOP);
    gmx_ana_set_nrefgrps(trj, 1);
    gmx_ana_set_nanagrps(trj, -1);
    /* If required, other functions can also be used to configure the library
     * before calling parse_trjana_args(). */
    parse_trjana_args(trj, &argc, argv, PCA_CAN_VIEW,
                      NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL,
                      &oenv);

    /* You can now do any initialization you wish, using the information
     * from the trj structure.
     * In particular, you should store any command-line parameter values that
     * the analysis part requires into d for them to be accessible. */

    /* First, we get some selection information from the structure */
    gmx_ana_get_refsel(trj, 0, &d.refsel);
    gmx_ana_get_nanagrps(trj, &ngrps);
    gmx_ana_get_anagrps(trj, &sel);

    /* First, we initialize the neighborhood search for the first index
     * group. */
    rc = gmx_ana_nbsearch_create(&d.nb, cutoff, d.refsel->p.nr);
    if (rc != 0 || ngrps>1)
    {
        gmx_fatal(FARGS, "neighborhood search initialization failed");
    }

    /* We then allocate memory in d to store intermediate data
     * and initialize the counters to zero */
    snew(d.ave, ngrps);
    snew(d.n,   ngrps);
    d.D=200; d.max=18.0; d.dr=d.max/d.D; d.tol=d.dr; //defaults

    snew(d.R,d.D*d.D*6);
    
    /* We also open the output file if the user provided it */
    d.fp = NULL;
    if (opt2bSet("-o", NFILE, fnm))
    {
        d.fp = xvgropen(opt2fn("-o", NFILE, fnm), "Average distance",
                        "Time [ps]", "Nr ", oenv);
        xvgr_selections(d.fp, trj);
    }

    /* Now, we do the actual analysis */
    gmx_ana_do(trj, 0, &analyze_frame, &d);

    fout=fopen("plaque_u.dat","w");
    for(i=0;i<d.D;i++){
      for(j=0;j<d.D;j++){
        fprintf(fout," % 6.3f",d.R[i*d.D+j]/d.count);
      }
      fprintf(fout,"\n");
    }
    fclose(fout);
    fout=fopen("plaque_b.dat","w");
    for(i=0;i<d.D;i++){
      for(j=0;j<d.D;j++){
        fprintf(fout," % 6.3f",d.R[i*d.D+j+d.D*d.D]/d.count);
      }
      fprintf(fout,"\n");
    }
    fclose(fout);

    if (d.fp)
    {
        ffclose(d.fp);
    }

    thanx(stderr);
    return 0;
}

/*! \brief
 * The main function.
 *
 * In Gromacs, most analysis programs are implemented such that the \p main
 * function is only a wrapper for a \p gmx_something function that does all
 * the work, and that convention is also followed here. */
int
main(int argc, char *argv[])
{
    gmx_template(argc, argv);
    return 0;
}
