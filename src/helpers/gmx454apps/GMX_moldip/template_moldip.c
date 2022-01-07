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

#include <fftw3.h>
#include <gromacs/nbsearch.h>
#include <gromacs/trajana.h>

/*! \brief
 * Template analysis data structure.
 */
// some usefull constants
#define SQ(X)     (X)*(X)
#define CB(X)     (X)*SQ(X)
#define GMX2DEBYE 48.032028118608
#define SQG2D     SQ(GMX2DEBYE)
#define c_0       299792458.00
#define nm        0.000000001
#define k_b       1.3806503*nm*nm*nm*10000
#define e_q       1.60217646*nm*nm*0.1
#define CMI1      5.30883742345939

typedef struct
{
    gmx_ana_selection_t *refsel;
    FILE                *fp;			// the output file
    FILE                *acfp;			// autocorrelation file
    real                *acorr_mu;		// for keeping all the mu vectors
    real                *mu_h;			// for histograms of the intramolecular properties
    rvec		*r_h;		     	// for histograms of the intramolecular properties
    real                **mu;  	 		// all the mu of the current frame ( or averages )
    real                *M;
    rvec                MU,MU2;			// the selection total MU and MU^2
    real                counts;			// a floating point counter	
    int                 n_alloc;                // for autocorrelation calculation
    int                 teller;		        // for autocorrelation calculation
    int                 tel3;			// for autocorrelation calculation
    real                t0,t1;			// for autocorrelation calculation
    int                 nr_mu;			// the number of mu's that are in one frame
    real                dx,xmin,xmax;		// spatial histogram binning
    real                dmu,mumin,mumax;	// mu histograms binning            
    int                 NBINS;                  // the amount of bins
    real                *n;			// not really used
    gmx_ana_nbsearch_t  *nb;			// select code requirement
    real                Q[3];			// net CG charges
} t_analysisdata;

static int
analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
    t_analysisdata     *d = (t_analysisdata *)data;
    int                 g, i, j;
    real                frave;
    int                 rc;
    real                Lx,dr;
    rvec                mu_temp,r_d[3],rnorm,MU_loc;
    real                amu;
    int                 ret_code;

    Lx=fr->box[XX][XX];

    if(d->t0==0.0){
      d->t0=fr->time;
    }
    d->t1=fr->time;

    rc = gmx_ana_nbsearch_pos_init(d->nb, pbc, &d->refsel->p);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "Neighborhood search initialization failed");
    }

    if (d->teller >= d->n_alloc){
      d->n_alloc+=2*d->nr_mu;
      for(i=0;i<d->nr_mu;i++){
         srenew(d->mu[i],d->n_alloc*3);
         srenew(d->M,d->n_alloc*3);
      }
    }

    j=0;
    d->tel3=3*(d->teller);
    MU_loc[XX]=0.0; MU_loc[YY]=0.0; MU_loc[ZZ]=0.0;
//  HARD CODED 3 ATOMS PER CHARGE GROUP
    for (i = 0; i < d->refsel->p.nr; i+=3)
    {
//	CURRENT CHARGEGROUP CONTRIBUTION
        mu_temp[XX]=0.0; mu_temp[YY]=0.0; mu_temp[ZZ]=0.0;
        mu_temp[XX]+=d->refsel->p.x[i][XX]*d->Q[0]+d->refsel->p.x[i+1][XX]*d->Q[1]+d->refsel->p.x[i+2][XX]*d->Q[2];
        mu_temp[YY]+=d->refsel->p.x[i][YY]*d->Q[0]+d->refsel->p.x[i+1][YY]*d->Q[1]+d->refsel->p.x[i+2][YY]*d->Q[2];
        mu_temp[ZZ]+=d->refsel->p.x[i][ZZ]*d->Q[0]+d->refsel->p.x[i+1][ZZ]*d->Q[1]+d->refsel->p.x[i+2][ZZ]*d->Q[2];

        d->mu[j][d->tel3+XX] = mu_temp[XX];
        d->mu[j][d->tel3+YY] = mu_temp[YY];
        d->mu[j][d->tel3+ZZ] = mu_temp[ZZ];

//	NORM OF CONTRIBUTION
        amu = sqrt(SQ(mu_temp[XX])+SQ(mu_temp[YY])+SQ(mu_temp[ZZ]))*GMX2DEBYE;

//	THE TOTAL BOX/SELECTION FLUCTUATION
        d->MU[XX]+=mu_temp[XX];
        d->MU[YY]+=mu_temp[YY];
        d->MU[ZZ]+=mu_temp[ZZ];
        MU_loc[XX]+=mu_temp[XX];
        MU_loc[YY]+=mu_temp[YY];
        MU_loc[ZZ]+=mu_temp[ZZ];

        d->MU2[XX]+=mu_temp[XX]*mu_temp[XX];
        d->MU2[YY]+=mu_temp[YY]*mu_temp[YY];
        d->MU2[ZZ]+=mu_temp[ZZ]*mu_temp[ZZ];

//	THE INTRAMOLECULAR LENGTHS
        ret_code = pbc_dx_aiuc(pbc,d->refsel->p.x[i],d->refsel->p.x[i+1],r_d[0]);
        ret_code = pbc_dx_aiuc(pbc,d->refsel->p.x[i],d->refsel->p.x[i+2],r_d[1]);
        ret_code = pbc_dx_aiuc(pbc,d->refsel->p.x[i+1],d->refsel->p.x[i+2],r_d[2]);

        rnorm[0] = sqrt( r_d[0][XX]*r_d[0][XX]+r_d[0][YY]*r_d[0][YY]+r_d[0][ZZ]*r_d[0][ZZ] );
        rnorm[1] = sqrt( r_d[1][XX]*r_d[1][XX]+r_d[1][YY]*r_d[1][YY]+r_d[1][ZZ]*r_d[1][ZZ] );
        rnorm[2] = sqrt( r_d[2][XX]*r_d[2][XX]+r_d[2][YY]*r_d[2][YY]+r_d[2][ZZ]*r_d[2][ZZ] );

//      COLLECT HISTOGRAMS
        if( amu > d->mumin && amu < d->mumax ){
           d->mu_h[ (int)nearbyint( (amu-d->mumin)/d->dmu ) ]+=1;
        }

        if( rnorm[0] > d->xmin && rnorm[0] < d->xmax ){
           d->r_h[ (int)nearbyint( (rnorm[0]-d->xmin)/d->dx ) ][0]+=1;
        }

        if( rnorm[1] > d->xmin && rnorm[1] < d->xmax ){
           d->r_h[ (int)nearbyint( (rnorm[1]-d->xmin)/d->dx ) ][1]+=1;
        }

        if( rnorm[2] > d->xmin && rnorm[2] < d->xmax ){
           d->r_h[ (int)nearbyint( (rnorm[2]-d->xmin)/d->dx ) ][2]+=1;
        }
        j++;
    }
    d->counts+=1.0;
    d->teller++;
    d->M[d->tel3+XX]=MU_loc[XX];
    d->M[d->tel3+YY]=MU_loc[YY];
    d->M[d->tel3+ZZ]=MU_loc[ZZ];

    return 0;
}

int
gmx_template(int argc, char *argv[])
{
    const char         *desc[] = {
        "This is a template for writing your own analysis tools",
    };

    // Command-line arguments
    real                **korre;
    real                *Z;
    real                Q     =  -1.0;
    int                 NBINS =   100;
    int                 n,m,M;
    gmx_bool            bArg  = FALSE;
    t_pargs             pa[] = {
        {"-NBINS", FALSE, etINT, {&NBINS},
         "set the number of bins"},
        {"-Q", FALSE, etREAL, {&Q},
         "The charge"},
        {"-spectrum",   FALSE, etBOOL, {&bArg},
         "Generate spectrum instead of AC"},
    };

    t_filenm            fnm[] = {
        {efXVG, "-o", "histograms", ffOPTWR},
        {efXVG, "-oac", "acorr",  ffOPTWR},
    };

#define NFILE asize(fnm)
    int     npargs;
    t_pargs *ppa;

    gmx_ana_traj_t       *trj;
    output_env_t          oenv;
    t_analysisdata        d;
    int                   ngrps;
    gmx_ana_selection_t **sel;
    int                   g=0,i,j;
    int                   rc;
    const char            *acfile;
//  FFTW3
    fftw_complex          *in1,*in2,*in3,*out,*im1,*im2,*im3,*ut;
    fftw_plan             p1,p2,p3,po,pm1,pm2,pm3,put;


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
    npargs = asize(pa);

    rc = gmx_ana_nbsearch_create(&d.nb, 0.0, d.refsel->p.nr);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "neighborhood search initialization failed");
    }
    d.nr_mu   = (int)nearbyint(d.refsel->p.nr/3.0);
    fprintf(stderr,"NGRPS=%d\n% 10d % 10d\nCGS= %10d\n",ngrps, d.refsel->p.nr, sel[g]->p.nr,d.nr_mu);


//  INIT HISTOGRAMS
    d.xmin = 0.05;
    d.xmax = 0.20;
    d.dx   = (d.xmax-d.xmin)/(float)NBINS;

    d.mumin = 1.9;
    d.mumax = 3.2;
    d.dmu   = (d.mumax-d.mumin)/(float)NBINS;
    d.MU[XX]=0.0; d.MU[YY]=0.0; d.MU[ZZ]=0.0;
    d.MU2[XX]=0.0; d.MU2[YY]=0.0; d.MU2[ZZ]=0.0;

    d.NBINS=NBINS;
    d.Q[0]=Q; d.Q[1]=-0.5*Q; d.Q[2]=d.Q[1];
    fprintf(stderr," % f % f % f \n",d.Q[0],d.Q[1],d.Q[2]);

    snew(d.n,   ngrps);
    snew(d.mu_h,NBINS);
    snew(d.r_h,NBINS);   // oh1, oh2, hh
    snew(d.mu,d.nr_mu);
    snew(d.M ,d.nr_mu);

    d.n_alloc = 10;
    for(i=0;i<d.nr_mu;i++)
       srenew(d.mu[i],d.n_alloc*3);

    d.t0      = 0.0;
    d.teller  = 0;
    d.tel3    = 0;

    for(i=0;i<ngrps+1;i++){
      if(i==0){
        fprintf(stderr,"[REF]");
        gmx_ana_selection_print_info(d.refsel);
      }else{
        fprintf(stderr,"[STD]");
        gmx_ana_selection_print_info(sel[i-1]);
      }
    }

//  ALSO OPEN OUTPUT FOR RESULTS IF PRESENT
    d.fp = NULL;
    if (opt2bSet("-o", NFILE, fnm))
    {
        d.fp = xvgropen(opt2fn("-o", NFILE, fnm), "\\x m\\f{}\\sT\\N",
                        "\\x m\\f{} [ D ]", "COUNTS", oenv);
        xvgr_selections(d.fp, trj);
    }

    d.acfp = NULL;
    if (opt2bSet("-oac", NFILE, fnm))
    {
        d.acfp = xvgropen(opt2fn("-oac", NFILE, fnm), "\\x < m\\f{}(t)\\xm\\f{}(0) >\\f{}\\sT\\N",
                        "t [ ps ]", "acorr", oenv);
        xvgr_selections(d.acfp, trj);
    }

    gmx_ana_do(trj, 0, &analyze_frame, &d);

    if (d.fp)
    {
        fprintf(d.fp,"@target G0.S0\n@type xy\n");
        for(i=0;i<NBINS;i++)
          fprintf(d.fp," % f % f \n",d.mumin+i*d.dmu,d.mu_h[i]/d.counts);
        fprintf(d.fp,"&\n");
        for(g=1;g<=3;g++){
          fprintf(d.fp,"@target G1.S%d\n@type xy\n",g);
          for(i=0;i<NBINS;i++)
            fprintf(d.fp," % f % f \n",d.xmin+i*d.dx,d.r_h[i][g-1]/d.counts);
          fprintf(d.fp,"&\n");
        }
        ffclose(d.fp);
    }
    fprintf(stderr,"<  MU  > = (% 20.10f ,% 20.10f ,% 20.10f) [ D ] \n< MU^2 > = (% 20.10f ,% 20.10f ,% 20.10f) [ D^2 ]\n",
    d.MU[XX]/d.counts*GMX2DEBYE,
    d.MU[YY]/d.counts*GMX2DEBYE,
    d.MU[ZZ]/d.counts*GMX2DEBYE,
    d.MU2[XX]/d.counts*GMX2DEBYE*GMX2DEBYE,
    d.MU2[YY]/d.counts*GMX2DEBYE*GMX2DEBYE,
    d.MU2[ZZ]/d.counts*GMX2DEBYE*GMX2DEBYE);

    fprintf(stderr,"READ %10d FRAMES\n",d.teller);

    if(d.acfp)
    {
      fprintf(stderr,"DOING AUTOCORRELATIONS (%d %d %f ) \n",d.tel3,d.teller,d.mu[d.nr_mu-1][d.tel3]);
      M=d.teller;
      in1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
      in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
      in3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
      out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);

      im1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
      im2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
      im3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
      ut  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);

      for(j=0;j<d.nr_mu;j++){
        if(j==0){

          for(i=0;i<d.teller;i++){
            im1[i][0]=d.M[i*3 + 0];
            im1[i][1]=0.0;
            im2[i][0]=d.M[i*3 + 1];
            im2[i][1]=0.0;
            im3[i][0]=d.M[i*3 + 2];
            im3[i][1]=0.0;
          }
          pm1 = fftw_plan_dft_1d(M,im1,im1, FFTW_FORWARD, FFTW_ESTIMATE);
          pm2 = fftw_plan_dft_1d(M,im2,im2, FFTW_FORWARD, FFTW_ESTIMATE);
          pm3 = fftw_plan_dft_1d(M,im3,im3, FFTW_FORWARD, FFTW_ESTIMATE);
          fftw_execute(pm1);
          fftw_execute(pm2);
          fftw_execute(pm3);
          for(i=0;i<M;i++){
            ut[i][0] = (im1[i][0]+im2[i][0]+im3[i][0])*(im1[i][0]+im2[i][0]+im3[i][0])+(im1[i][1]+im2[i][1]+im3[i][1])*(im1[i][1]+im2[i][1]+im3[i][1]);
            ut[i][1] = 0.0;
          }
          put = fftw_plan_dft_1d(M,ut,ut, FFTW_BACKWARD, FFTW_ESTIMATE);
          if(!bArg){
            fftw_execute(put);
          }

        }
        for(i=0;i<d.teller;i++){
          in1[i][0]=d.mu[j][i*3 + 0];
          in1[i][1]=0.0;
          in2[i][0]=d.mu[j][i*3 + 1];
          in2[i][1]=0.0;
          in3[i][0]=d.mu[j][i*3 + 2];
          in3[i][1]=0.0;
        }
        p1 = fftw_plan_dft_1d(M,in1,in1, FFTW_FORWARD, FFTW_ESTIMATE);
        p2 = fftw_plan_dft_1d(M,in2,in2, FFTW_FORWARD, FFTW_ESTIMATE);
        p3 = fftw_plan_dft_1d(M,in3,in3, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p1);
        fftw_execute(p2);
        fftw_execute(p3);

        for(i=0;i<M;i++){

          if(j==0)
            out[i][0] = (in1[i][0]+in2[i][0]+in3[i][0])*(in1[i][0]+in2[i][0]+in3[i][0])+(in1[i][1]+in2[i][1]+in3[i][1])*(in1[i][1]+in2[i][1]+in3[i][1]);
          else
            out[i][0]+= (in1[i][0]+in2[i][0]+in3[i][0])*(in1[i][0]+in2[i][0]+in3[i][0])+(in1[i][1]+in2[i][1]+in3[i][1])*(in1[i][1]+in2[i][1]+in3[i][1]);
	  /*
          if(j==0)                                                                                                                                                                             
            out[i][0] = SQ(in1[i][0]*in1[i][1]+in2[i][0]*in2[i][1]+in3[i][0]*in3[i][1]);                                        
          else                                                                                                                                                                                 
            out[i][0]+= SQ(in1[i][0]*in1[i][1]+in2[i][0]*in2[i][1]+in3[i][0]*in3[i][1]); 
	  */
          out[i][1] = 0.0 ;
        }
      }
      po = fftw_plan_dft_1d(M,out,out, FFTW_BACKWARD, FFTW_ESTIMATE);
      if(!bArg){
        fftw_execute(po);
      }
      M=nearbyint(M*0.5);
      if(bArg){
        for(i=0;i<M;i++)
          fprintf(d.acfp," % 20.10f % 20.10f % 20.10f\n",M_PI*i/((d.t1-d.t0)/d.teller)*CMI1/M,(out[i][0]-out[M][0])/out[0][0],(ut[i][0]-ut[M][0])/ut[0][0] );
      }else{
        for(i=0;i<M;i++)
          fprintf(d.acfp," % 20.10f % 20.10f % 20.10f\n",i*((d.t1-d.t0)/d.teller),out[i][0]/out[0][0],ut[i][0]/ut[0][0]);
      }
      fprintf(stderr,"DONE AUTOCORRELATIONS\n");
      ffclose(d.acfp);
      fftw_destroy_plan(p1);
      fftw_free(in1);
      fftw_destroy_plan(p2);
      fftw_free(in2);
      fftw_destroy_plan(p3);
      fftw_free(in3);
      fftw_destroy_plan(po);
      fftw_free(out);
      fftw_destroy_plan(pm1);
      fftw_free(im1);
      fftw_destroy_plan(pm2);
      fftw_free(im2);
      fftw_destroy_plan(pm3);
      fftw_free(im3);
      fftw_destroy_plan(put);
      fftw_free(ut);
    }

    thanx(stderr);
    return 0;
}

int
main(int argc, char *argv[])
{
    gmx_template(argc, argv);
    return 0;
}
