
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

typedef struct
{
    gmx_ana_selection_t *refsel;
    FILE                *fp;
    double              *ave;
    real                *n;
    real                 dr;
    gmx_bool             b2D;
    int                  NBINS;
    real                 counts;
    gmx_ana_nbsearch_t  *nb;
} t_analysisdata;

static int
analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
    t_analysisdata     *d = (t_analysisdata *)data;
    int                 g,i,j,k,l;
    real                dr;
    int                 rc;
    rvec                drp;
    d->counts+=1.0;

    rc = gmx_ana_nbsearch_pos_init(d->nb, pbc, &d->refsel->p);
    if (rc != 0)
    {
        gmx_fatal(FARGS, "Neighborhood search initialization failed");
    }
    for (j = 0; j < d->refsel->p.nr; ++j)
    {
        for (g = 0; g < nr; ++g)
        {
            for (i = 0; i < sel[g]->p.nr; ++i)
            {
                rc = pbc_dx_aiuc(pbc,d->refsel->p.x[j],sel[g]->p.x[i],drp);
                dr = sqrt(drp[XX]*drp[XX]+drp[YY]*drp[YY]+drp[ZZ]*drp[ZZ]);
                if(dr<=d->dr*(d->NBINS+1)){
                    k = nearbyint(  dr/(d->dr) );
                    for(l=k;l<=d->NBINS;l++)
                        d->ave[l]+=1.0/(d->refsel->p.nr);
                    d->n[k]+=1.0/(d->refsel->p.nr);
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

    real                cutoff = 0;
    gmx_bool            bArg   = TRUE;
    real                L = 1.0;
    int                 NBINS=100;
    t_pargs             pa[] = {
        {"-cutoff", FALSE, etREAL, {&cutoff},
         "Cutoff for distance calculation (0 = no cutoff)"},
        {"-NBINS", FALSE, etINT, {&NBINS},
         "set the number of bins in the chi plot"},
        {"-L", FALSE, etREAL, {&L},
         "Maximum bin value"},
        {"-arg2",   FALSE, etBOOL, {&bArg},
         "2D proj"},
    };

    t_filenm            fnm[] = {
        {efXVG, "-o", "avedist", ffOPTWR},
    };
#define NFILE asize(fnm)

    gmx_ana_traj_t       *trj;
    output_env_t          oenv;
    t_analysisdata        d;
    int                   ngrps;
    gmx_ana_selection_t **sel;
    int                   g,i;
    int                   rc;
    real                  sz[2];

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
    if (rc != 0)
    {
        gmx_fatal(FARGS, "neighborhood search initialization failed");
    }

    d.dr=L/(float)NBINS;
    d.NBINS=NBINS;
    d.b2D=bArg;

    snew(d.ave, ngrps*NBINS*2);
    snew(d.n,   ngrps*NBINS*2);

    d.fp = NULL;
    if (opt2bSet("-o", NFILE, fnm))
    {
        d.fp = xvgropen(opt2fn("-o", NFILE, fnm), "\\x x\\f{}",
                        "r [nm]", "\\x x", oenv);
        xvgr_selections(d.fp, trj);
    }
    gmx_ana_do(trj, 0, &analyze_frame, &d);

    fprintf(stderr,"%f %d %d %f\n",d.dr,d.NBINS,d.b2D,d.counts);
    
    if (d.fp) {
        d.n[0]=0.0;
        for(i=0;i<NBINS;i++){
            sz[0]=i*d.dr*i*d.dr;
            sz[1]=(i+1)*d.dr*(i+1)*d.dr;
            fprintf(d.fp,"% 10f % 10.6f % 10.6f\n ",i*d.dr,d.ave[i]/d.counts-d.ave[0]/d.counts,d.n[i]/d.counts);
        }
        ffclose(d.fp);
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
