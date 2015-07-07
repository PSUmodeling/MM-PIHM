#include "pihm.h"
#include "noah.h"

void SFlx (grid_struct * grid)
{
/*----------------------------------------------------------------------
* subroutine sflx - unified noahlsm version 1.0 july 2007
* ----------------------------------------------------------------------
* sub-driver for "noah lsm" family of physics subroutines for a
* soil/veg/snowpack land-surface model to update soil moisture, soil
* ice, soil temperature, skin temperature, snowpack water content,
* snowdepth, and all terms of the surface energy balance and surface
* water balance (excluding input atmospheric forcings of downward
* radiation and precip)
* --------------------------------------------------------------------*/
    int            *frzgra, *snowng;

    int            *nsoil, *vegtyp;
    int            *isurban;
    int            *nroot;
    int             kz, k, iout;

    int             rdlai2d;
    int             usemonalb;

    double         *shdmin, *shdmax, *dt, *dqsdt2, *lwdn, *prcp, *prcprain,
       *q2, *q2sat, *sfcprs, *sfcspd, *sfctmp, *snoalb, *soldn, *solnet,
       *tbot, *th2, *zlvl, *ffrozp;
    double         *embrd;
    double         *albedo;
    double         *cosz, *solardirect, *ch, *cm, *cmc, *sneqv, *sncovr,
       *snowh, *t1, *xlai, *shdfac, *z0brd, *emissi, *alb;
    double         *snotime1;
    double         *ribb;
    double         *sldpth;
    double         *et;
    double         *smav;
    double         *sh2o, *smc, *stc;
    double         *rtdis;
    double          zsoil[grid->nsoil];

    double         *eta_kinematic, *beta, *dew, *drip, *ec, *edir, *esnow, *eta, *etp, *flx1, *flx2, *flx3, *sheat, *pc, *runoff2, *runoff3, *rc, *rsmin, *rcq, *rcs, *rcsoil, *rct, *ssoil, *smcdry, *smcmax, *smcref, *smcwlt, *snomlt, *soilm, *soilw, *fdown, *q1;
    double         *cfactr, *cmcmax, *csoil, *czil, *df1, *dksat, *ett, *epsca, *f1, *fxexp, *frzx, *hs, *quartz, *rch, *rr, *rgl, *rsmax, *sndens, *sncond, *sbeta, *sn_new, *snup, *salp, *t24, *t2v, *topt, *zbot, *z0, *prcpf, *etns, *ptu;
    double          df1h, df1a, dsoil, dtot, frcsno, frcsoi, soilwm, soilww;
    double         *lvcoef;
    double          interp_fraction;
    double         *laimin, *laimax;
    double         *albedomin, *albedomax;
    double         *emissmin, *emissmax;
    double         *z0min, *z0max;

#ifdef _NOAH_
    double         *vgalpha, *vgbeta, *smcmin, *macksat, *areaf, *infil;
    double         *pcpdrp;
    int            *nmacd, *nwtbl;
    int            *mac_status;
#else
    double         *runoff1, *bexp, *dwsat, *kdt, *psisat, *refkdt, *slope;
#endif

/*----------------------------------------------------------------------
*   initialization
* --------------------------------------------------------------------*/

#ifndef _NOAH_
    grid->runoff1 = 0.0;
    grid->runoff2 = 0.0;
    grid->runoff3 = 0.0;
#endif
    grid->snomlt = 0.0;

/*----------------------------------------------------------------------
* calculate depth (negative) below ground from top skin sfc to bottom of
* each soil layer.  note:  sign of zsoil is negative (denoting below
* ground)
* --------------------------------------------------------------------*/

    zsoil[0] = -grid->sldpth[0];
    for (kz = 1; kz < grid->nsoil; kz++)
        zsoil[kz] = -grid->sldpth[kz] + zsoil[kz - 1];

/*----------------------------------------------------------------------
* next is crucial call to set the land-surface parameters, including
* soil-type and veg-type dependent parameters.
* --------------------------------------------------------------------*/

    //  RedPrm(grid, lsm, zsoil);           /* ys: RedPrm is now called in driver */

    frzgra = (int *)malloc (sizeof (int));
    snowng = (int *)malloc (sizeof (int));

    nsoil = &(grid->nsoil);
#ifndef _NOAH_
    slopetyp = &(grid->slopetyp);
    soiltyp = &(grid->soiltyp);
#endif
    vegtyp = &(grid->vegtyp);
    isurban = &(grid->isurban);
    nroot = &(grid->nroot);

    rdlai2d = grid->rdlai2d;
    usemonalb = grid->usemonalb;

    shdmin = &(grid->shdmin);
    shdmax = &(grid->shdmax);
    dt = &(grid->dt);
    dqsdt2 = &(grid->dqsdt2);
    lwdn = &(grid->lwdn);
    prcp = &(grid->prcp);
    prcprain = &(grid->prcprain);
    q2 = &(grid->q2);
    q2sat = &(grid->q2sat);
    sfcprs = &(grid->sfcprs);
    sfcspd = &(grid->sfcspd);
    sfctmp = &(grid->sfctmp);
    snoalb = &(grid->snoalb);
    soldn = &(grid->soldn);
    solnet = &(grid->solnet);
    tbot = &(grid->tbot);
    th2 = &(grid->th2);
    zlvl = &(grid->zlvl);
    ffrozp = &(grid->ffrozp);

    embrd = &(grid->embrd);
    albedo = &(grid->albedo);
    cosz = &(grid->cosz);
    solardirect = &(grid->solardirect);
    ch = &(grid->ch);
    cm = &(grid->cm);
    cmc = &(grid->cmc);
    sneqv = &(grid->sneqv);
    sncovr = &(grid->sncovr);
    snowh = &(grid->snowh);
    t1 = &(grid->t1);
    xlai = &(grid->xlai);
    shdfac = &(grid->shdfac);
    z0brd = &(grid->z0brd);
    emissi = &(grid->emissi);
    alb = &(grid->alb);

    snotime1 = &(grid->snotime1);

    ribb = &(grid->ribb);

    sldpth = &(grid->sldpth[0]);
    et = &(grid->et[0]);
    smav = &(grid->smav[0]);
    sh2o = &(grid->sh2o[0]);
    smc = &(grid->smc[0]);
    stc = &(grid->stc[0]);

    rtdis = &(grid->rtdis[0]);

    eta_kinematic = &(grid->eta_kinematic);
    beta = &(grid->beta);
    dew = &(grid->dew);
    drip = &(grid->drip);
    ec = &(grid->ec);
    edir = &(grid->edir);
    esnow = &(grid->esnow);
    eta = &(grid->eta);
    etp = &(grid->etp);
    flx1 = &(grid->flx1);
    flx2 = &(grid->flx2);
    flx3 = &(grid->flx3);
    sheat = &(grid->sheat);
    pc = &(grid->pc);
    runoff2 = &(grid->runoff2);
    runoff3 = &(grid->runoff3);
    rc = &(grid->rc);
    rsmin = &(grid->rsmin);
    rcq = &(grid->rcq);
    rcs = &(grid->rcs);
    rcsoil = &(grid->rcsoil);
    rct = &(grid->rct);
    ssoil = &(grid->ssoil);
    smcdry = &(grid->smcdry);
    smcmax = &(grid->smcmax);
    smcref = &(grid->smcref);
    smcwlt = &(grid->smcwlt);
    snomlt = &(grid->snomlt);
    soilm = &(grid->soilm);
    soilw = &(grid->soilw);
    fdown = &(grid->fdown);
    q1 = &(grid->q1);
#ifdef _NOAH_
    smcmin = &(grid->smcmin);
    vgalpha = &(grid->vgalpha);
    vgbeta = &(grid->vgbeta);
    macksat = &(grid->macksat);
    areaf = &(grid->areaf);
    infil = &(grid->infil);
    nmacd = &(grid->nmacd);
    mac_status = &(grid->mac_status);
    nwtbl = &(grid->nwtbl);
    pcpdrp = &(grid->pcpdrp);
#else
    runoff1 = &(grid->runoff1);
    bexp = &(grid->bexp);
    dwsat = &(grid->dwsat);
    kdt = &(grid->kdt);
    psisat = &(grid->psisat);
    slope = &(grid->slope);
    refkdt = &(grid->refkdt);
#endif
    cfactr = &(grid->cfactr);
    cmcmax = &(grid->cmcmax);
    csoil = &(grid->csoil);
    czil = &(grid->czil);
    df1 = (double *)malloc (sizeof (double));
    dksat = &(grid->dksat);
    ett = &(grid->ett);
    epsca = (double *)malloc (sizeof (double));
    f1 = &(grid->f1);
    fxexp = &(grid->fxexp);
    frzx = &(grid->frzx);
    hs = &(grid->hs);
    quartz = &(grid->quartz);
    rch = (double *)malloc (sizeof (double));
    rr = (double *)malloc (sizeof (double));
    rgl = &(grid->rgl);
    rsmax = &(grid->rsmax);
    sndens = (double *)malloc (sizeof (double));
    sncond = (double *)malloc (sizeof (double));
    sbeta = &(grid->sbeta);
    sn_new = (double *)malloc (sizeof (double));
    snup = &(grid->snup);
    salp = &(grid->salp);
    t24 = (double *)malloc (sizeof (double));
    t2v = (double *)malloc (sizeof (double));
    topt = &(grid->topt);
    zbot = &(grid->zbot);
    z0 = &(grid->z0);
    prcpf = (double *)malloc (sizeof (double));
    etns = (double *)malloc (sizeof (double));
    ptu = &(grid->ptu);

    lvcoef = &(grid->lvcoef);
    laimin = &(grid->laimin);
    laimax = &(grid->laimax);
    albedomin = &(grid->albedomin);
    albedomax = &(grid->albedomax);
    emissmin = &(grid->emissmin);
    emissmax = &(grid->emissmax);
    z0min = &(grid->z0min);
    z0max = &(grid->z0max);

#ifdef _NOAH_
    *pcpdrp = 0.;
#endif


    /*
     * urban 
     */
    if (*vegtyp == *isurban)
    {
        *shdfac = 0.05;
        *rsmin = 400.0;
        *smcmax = 0.45;
#ifdef _NOAH_
        *smcmin = 0.0;
#endif
        *smcref = 0.42;
        *smcwlt = 0.40;
        *smcdry = 0.40;
    }

#ifdef _NOAH_

/*----------------------------------------------------------------------
* ys: flux-pihm uses lai as a forcing variable
* vegetation fraction is calculated from lai following noah-mp
* --------------------------------------------------------------------*/

    if (*xlai >= *laimax)
    {
        *embrd = *emissmax;
        *alb = *albedomin;
        *z0brd = *z0max;
    }
    else if (*xlai <= *laimin)
    {
        *embrd = *emissmin;
        *alb = *albedomax;
        *z0brd = *z0min;
    }
    else
    {
        if (*laimax > *laimin)
        {
            interp_fraction = (*xlai - *laimin) / (*laimax - *laimin);
            /*
             * bound interp_fraction between 0 and 1 
             */
            interp_fraction = interp_fraction < 1.0 ? interp_fraction : 1.0;
            interp_fraction = interp_fraction > 0.0 ? interp_fraction : 0.0;
            /*
             * scale emissivity and lai between emissmin and emissmax by interp_fraction 
             */
            *embrd = ((1.0 - interp_fraction) * *emissmin) + (interp_fraction * *emissmax);
            *alb = ((1.0 - interp_fraction) * *albedomax) + (interp_fraction * *albedomin);
            *z0brd = ((1.0 - interp_fraction) * *z0min) + (interp_fraction * *z0max);
        }
        else
        {
            *embrd = 0.5 * *emissmin + 0.5 * *emissmax;
            *alb = 0.5 * *albedomin + 0.5 * *albedomax;
            *z0brd = 0.5 * *z0min + 0.5 * *z0max;
        }
    }

//    *shdfac = 1. - exp (-0.52 * (*xlai));
      *shdfac = 1. - exp (-0.75 * (*xlai));
#else
    if (*shdfac >= *shdmax)
    {
        *embrd = *emissmax;
        if (!rdlai2d)
            *xlai = *laimax;
        if (!usemonalb)
            *alb = *albedomin;
        *z0brd = *z0max;
    }
    else if (*shdfac <= *shdmin)
    {
        *embrd = *emissmin;
        if (!rdlai2d)
            *xlai = *laimin;
        if (!usemonalb)
            *alb = *albedomax;
        *z0brd = *z0min;
    }
    else
    {
        if (*shdmax > *shdmin)
        {
            interp_fraction = (*shdfac - *shdmin) / (*shdmax - *shdmin);
            /*
             * bound interp_fraction between 0 and 1 
             */
            interp_fraction = interp_fraction < 1.0 ? interp_fraction : 1.0;
            interp_fraction = interp_fraction > 0.0 ? interp_fraction : 0.0;
            /*
             * scale emissivity and lai between emissmin and emissmax by interp_fraction 
             */
            *embrd =
               ((1.0 - interp_fraction) * *emissmin) +
               (interp_fraction * *emissmax);
            if (!rdlai2d)
                *xlai =
                   ((1.0 - interp_fraction) * *laimin) +
                   (interp_fraction * *laimax);
            if (!usemonalb)
                *alb =
                   ((1.0 - interp_fraction) * *albedomax) +
                   (interp_fraction * *albedomin);
            *z0brd =
               ((1.0 - interp_fraction) * *z0min) +
               (interp_fraction * *z0max);
        }
        else
        {
            *embrd = 0.5 * *emissmin + 0.5 * *emissmax;
            if (!rdlai2d)
                *xlai = 0.5 * *laimin + 0.5 * *laimax;
            if (!usemonalb)
                *alb = 0.5 * *albedomin + 0.5 * *albedomax;
            *z0brd = 0.5 * *z0min + 0.5 * *z0max;
        }
    }
#endif

    /*
     * initialize precipitation logicals. 
     */
    *snowng = 0;
    *frzgra = 0;

/*----------------------------------------------------------------------
* if input snowpack is nonzero, then compute snow density "sndens" and
* snow thermal conductivity "sncond" (note that CSnow is a function
* subroutine)
* --------------------------------------------------------------------*/
    if (*sneqv <= 1.0e-7)       /* safer if kmh (2008/03/25) */
    {
        *sneqv = 0.0;
        *sndens = 0.0;
        *snowh = 0.0;
        *sncond = 1.0;
    }
    else
    {
        *sndens = *sneqv / *snowh;
        if (*sndens > 1.0)
        {
            printf ("physical snow depth is less than snow water equiv.\n");
            exit (0);
        }
        CSnow (sncond, sndens);
    }

#ifdef _DEBUG_
    //printf ("nsoil = %d, isurban = %d, nroot = %d\n", *nsoil, *isurban, *nroot);
    //printf ("rdlai2d = %d, usemonalb = %d\n", rdlai2d, usemonalb);
    //printf ("shdmin = %f, shdmax = %f, dt = %f, dqsdt2 = %f, lwdn = %f, prcp = %f, prcprain = %f, q2 = %f, q2sat = %f, sfcprs = %f, sfcspd = %f, sfctmp = %f, snoalb = %f, soldn = %f, solnet = %f, tbot = %f, th2, = %f, zlvl = %f, ffrozp = %f\n", *shdmin, *shdmax, *dt, *dqsdt2, *lwdn, *prcp, *prcprain, *q2, *q2sat, *sfcprs, *sfcspd, *sfctmp, *snoalb, *soldn, *solnet, *tbot, *th2, *zlvl, *ffrozp);
    //printf ("ch = %f, cm = %f, cmc = %f, sneqv = %f, sncovr = %f, snowh = %f, t1 = %f, xlai = %f, shdfac = %f, z0brd = %f, emissi = %f, alb = %f\n", *ch, *cm, *cmc, *sneqv, *sncovr, *snowh, *t1, *xlai, *shdfac, *z0brd, *emissi, *alb);
    //printf ("snotime1 = %f\n", *snotime1);
    //printf ("ribb = %f\n", *ribb);
    //for (kz = 0; kz < *nsoil; kz++) printf ("sldpth[%d] = %f ", kz, sldpth[kz]);
    //printf ("\n");
    //for (kz = 0; kz < *nsoil; kz++)
    //    printf ("sh2o[%d] = %f ", kz, sh2o[kz]);
    //printf ("\n");
    //for (kz = 0; kz < *nsoil; kz++)
    //    printf ("smc[%d] = %f ", kz, smc[kz]);
    //printf ("\n");
    //for (kz = 0; kz < *nsoil; kz++)
    //    printf ("stc[%d] = %f ", kz, stc[kz]);
    //printf ("\n");
    //for (kz = 0; kz < *nsoil; kz++)
    //    printf ("rtdis[%d] = %f ", kz, rtdis[kz]);
    //printf ("\n");
    //for (kz = 0; kz < *nsoil; kz++)
    //    printf ("zsoil[%d] = %f ", kz, zsoil[kz]);
    //printf ("\n");
    //printf ("vgalpha = %f, vgbeta = %f, cfactr = %f, cmcmax = %f, csoil = %f, czil = %f, df1 = %f, dksat = %f, fxexp = %f, frzx = %f, LVH2O = %f, rgl = %f, rsmax = %f, sbeta = %f, topt = %f, hs = %f, zbot = %f, snup = %f, salp = %f, smcmax = %f, smcwlt = %f, smcref = %f, smcdry = %f, quartz = %f, laimin = %f, laimax = %f, emissmin = %f, emissimax = %f, albedomin = %f, albedomax = %f, z0min = %f, z0max = %f\n", *vgalpha, *vgbeta, *cfactr, *cmcmax, *csoil, *czil, *df1, *dksat, *fxexp, *frzx, LVH2O, *rgl, *rsmax, *sbeta, *topt, *hs, *zbot, *snup, *salp, *smcmax, *smcwlt, *smcref, *smcdry, *quartz, *laimin, *laimax, *emissmin, *emissmax, *albedomin, *albedomax, *z0min, *z0max);
#endif

/*----------------------------------------------------------------------
* determine if it's precipitating and what kind of precip it is.
* if it's prcping and the air temp is colder than 0 c, it's snowing!
* if it's prcping and the air temp is warmer than 0 c, but the grnd
* temp is colder than 0 c, freezing rain is presumed to be falling.
* --------------------------------------------------------------------*/

    if (*prcp > 0.0)

        /*
         * snow defined when fraction of frozen precip (ffrozp) > 0.5,
         * passed in from model microphysics. 
         */
    {
        if (*ffrozp > 0.5)
            *snowng = 1;
        else
        {
            if (*t1 <= TFREEZ)
                *frzgra = 1;
        }
    }

/*----------------------------------------------------------------------
* if either prcp flag is set, determine new snowfall (converting prcp
* rate from kg m-2 s-1 to a liquid equiv snow depth in meters) and add
* it to the existing snowpack.
* note that since all precip is added to snowpack, no precip infiltrates
* into the soil so that prcp1 is set to zero.
* --------------------------------------------------------------------*/
    if (*snowng || *frzgra)
    {
        *sn_new = *prcp * *dt * 0.001;
        *sneqv = *sneqv + *sn_new;
        *prcpf = 0.0;

/*----------------------------------------------------------------------
* update snow density based on new snowfall, using old and new snow.
* update snow thermal conductivity
* --------------------------------------------------------------------*/
        SnowNew (sfctmp, sn_new, snowh, sndens);
        CSnow (sncond, sndens);
    }

/*----------------------------------------------------------------------
* precip is liquid (rain), hence save in the precip variable that
* later can wholely or partially infiltrate the soil (along with
* any canopy "drip" added to this later)
* --------------------------------------------------------------------*/
    else
        *prcpf = *prcp;

/*----------------------------------------------------------------------
* determine snowcover and albedo over land.
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* if snow depth=0, set snow fraction=0, albedo=snow free albedo.
* --------------------------------------------------------------------*/
    if (*sneqv == 0.0)
    {
        *sncovr = 0.0;
        *albedo = *alb;
        *emissi = *embrd;
    }
    else
    {

/*----------------------------------------------------------------------
* determine snow fractional coverage.
* determine surface albedo modification due to snowdepth state.
* --------------------------------------------------------------------*/
        SnFrac (sneqv, snup, salp, snowh, sncovr);

        *sncovr = *sncovr < 0.98 ? *sncovr : 0.98;

        AlCalc (alb, snoalb, embrd, shdfac, shdmin, sncovr, t1, albedo,
           emissi, dt, snowng, snotime1, lvcoef);
    }

/*----------------------------------------------------------------------
* next calculate the subsurface heat flux, which first requires
* calculation of the thermal diffusivity.  treatment of the
* latter follows that on pages 148-149 from "heat transfer in
* cold climates", by v. j. lunardini (published in 1981
* by van nostrand reinhold co.) i.e. treatment of two contiguous
* "plane parallel" mediums (namely here the first soil layer
* and the snowpack layer, if any). this diffusivity treatment
* behaves well for both zero and nonzero snowpack, including the
* limit of very thin snowpack.  this treatment also eliminates
* the need to impose an arbitrary upper bound on subsurface
* heat flux when the snowpack becomes extremely thin.
* ----------------------------------------------------------------------
* first calculate thermal diffusivity of top soil layer, using
* both the frozen and liquid soil moisture, following the
* soil thermal diffusivity function of peters-lidard et al.
* (1998,jas, vol 55, 1209-1224), which requires the specifying
* the quartz content of the given soil class (see routine RedPrm)
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* next add subsurface heat flux reduction effect from the
* overlying green canopy, adapted from section 2.1.2 of
* peters-lidard et al. (1997, jgr, vol 102(d4))
* --------------------------------------------------------------------*/
#ifdef _NOAH_
    TDfCnd (df1, smc, quartz, smcmax, smcmin, sh2o);
#else
    TDfCnd (df1, smc, quartz, smcmax, sh2o);
#endif

    /*
     * urban 
     */
    if (*vegtyp == *isurban)
        *df1 = 3.24;

    *df1 = *df1 * exp (*sbeta * *shdfac);

    /*
     * kmh 09/03/2006
     * kmh 03/25/2008  change sncovr threshold to 0.97
     */
    if (*sncovr > 0.97)
        *df1 = *sncond;

/*----------------------------------------------------------------------
* finally "plane parallel" snowpack effect following
* v.j. linardini reference cited above. note that dtot is
* combined depth of snowdepth and thickness of first soil layer
* --------------------------------------------------------------------*/

    dsoil = -(0.5 * zsoil[0]);
    if (*sneqv == 0.)
        *ssoil = *df1 * (*t1 - stc[0]) / dsoil;
    else
    {
        dtot = *snowh + dsoil;
        frcsno = *snowh / dtot;

        /*
         * 1. harmonic mean (series flow) 
         */
        //      df1 = (sncond*df1)/(frcsoi*sncond+frcsno*df1)
        frcsoi = dsoil / dtot;

        /*
         * 2. arithmetic mean (parallel flow) 
         */
        //      df1 = frcsno*sncond + frcsoi*df1
        df1h = (*sncond * *df1) / (frcsoi * *sncond + frcsno * *df1);

        /*
         * 3. geometric mean (intermediate between harmonic and arithmetic mean) 
         */
        //      df1 = (sncond**frcsno)*(df1**frcsoi)
        /*
         * weigh df by snow fraction 
         */
        //      df1 = df1h*sncovr + df1a*(1.0-sncovr)
        //      df1 = df1h*sncovr + df1*(1.0-sncovr)
        df1a = frcsno * *sncond + frcsoi * *df1;

/*----------------------------------------------------------------------
* calculate subsurface heat flux, ssoil, from final thermal diffusivity
* of surface mediums, df1 above, and skin temperature and top
* mid-layer soil temperature
* --------------------------------------------------------------------*/
        *df1 = df1a * *sncovr + *df1 * (1.0 - *sncovr);
        *ssoil = *df1 * (*t1 - stc[0]) / dtot;
    }

/*----------------------------------------------------------------------
* determine surface roughness over snowpack using snow condition from
* the previous timestep.
* --------------------------------------------------------------------*/
    if (*sncovr > 0.)
        Snowz0 (sncovr, z0, z0brd, snowh);
    else
        *z0 = *z0brd;

/*----------------------------------------------------------------------
* next call routine sfcdif to calculate the sfc exchange coef (ch) for
* heat and moisture.

* note !!!
* do not call sfcdif until after above call to RedPrm, in case
* alternative values of roughness length (z0) and zilintinkevich coef
* (czil) are set there via namelist i/o.

* note !!!
* routine sfcdif returns a ch that represents the wind spd times the
* "original" nondimensional "ch" typical in literature.  hence the ch
* returned from sfcdif has units of m/s.  the important companion
* coefficient of ch, carried here as "rch", is the ch from sfcdif times
* air density and parameter "CP".  "rch" is computed in "call Penman".
* rch rather than ch is the coeff usually invoked later in eqns.

* note !!!
*-----------------------------------------------------------------------
* sfcdif also returns the surface exchange coefficient for momentum, cm,
* also known as the surface drage coefficient. needed as a state variable
* for iterative/implicit solution of ch in sfcdif
* --------------------------------------------------------------------*/

    /*
     * if (!lch)
     * {
     * *t1v = *t1 * (1.0 + 0.61 * *q2);
     * *th2v = *th2 * (1.0 + 0.61 * *q2);
     * SfcDifOff(zlvl, z0, t1v, th2v, sfcspd, czil, cm, ch);
     * }
     */

/*----------------------------------------------------------------------
* call Penman subroutine to calculate potential evaporation (etp), and
* other partial products and sums save in common/rite for later
* calculations.
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* calculate total downward radiation (solar plus longwave) needed in
* Penman ep subroutine that follows
* --------------------------------------------------------------------*/
    //  *fdown = *soldn * (1.0- *albedo) + *lwdn;
    *fdown = *solnet + *lwdn;

/*----------------------------------------------------------------------
* calc virtual temps and virtual potential temps needed by subroutines
* Penman.
* --------------------------------------------------------------------*/
    *t2v = *sfctmp * (1.0 + 0.61 * *q2);

#ifdef _DEBUG_
    iout = 1;
#else
    iout = 0;
#endif
    //if (iout == 1)
    //{
    //    printf ("before Penman\n");
    //    printf ("sfctmp = %lf sfcprs = %lf ch = %lf t2v = %lf th2 = %lf prcp = %lf fdown = %lf t24 = %lf ssoil = %lf q2 = %f q2sat = %lf etp = %lf rch = %lf epsca = %lf rr = %lf snowng = %d frzgra = %d dqsdt2 = %lf flx2 = %lf snowh = %lf sneqv = %lf dsoil = %lf frcsno = %lf sncovr = %lf dtot = %lf zsoil(1) = %lf df1 = %lf t1 = %lf stc1 = %lf albedo = %lf smc = %lf stc = %lf sh2o = %lf\n", *sfctmp, *sfcprs, *ch, *t2v, *th2, *prcp, *fdown, *t24, *ssoil, *q2, *q2sat, *etp, *rch, *epsca, *rr, *snowng, *frzgra, *dqsdt2, *flx2, *snowh, *sneqv, dsoil, frcsno, *sncovr, dtot, zsoil[0], *df1, *t1, stc[0], *albedo, smc[0], stc[1], sh2o[0]);
    //    //      for (k = 0; k < *nsoil; k++)
    //    //          printf("sh2o = %f, smc = %f\t", sh2o[k], smc[k]);
    //    //      printf("\n");
    //}

    Penman (sfctmp, sfcprs, ch, t2v, th2, prcp, fdown, t24, ssoil, q2, q2sat, etp, rch, epsca, rr, snowng, frzgra, dqsdt2, flx2, emissi, sneqv, t1, sncovr);

/*----------------------------------------------------------------------
* call CanRes to calculate the canopy resistance and convert it into pc
* if nonzero greenness fraction
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
*  frozen ground extension: total soil water "smc" was replaced
*  by unfrozen soil water "sh2o" in call to CanRes below
* --------------------------------------------------------------------*/
    if (*shdfac > 0.)
        CanRes (soldn, ch, sfctmp, q2, sfcprs, sh2o, zsoil, nsoil, smcwlt,
           smcref, rsmin, rc, pc, nroot, q2sat, dqsdt2, topt, rsmax, rgl, hs,
           xlai, rcs, rct, rcq, rcsoil, emissi);
    else
        *rc = 0.0;

/*----------------------------------------------------------------------
* now decide major pathway branch to take depending on whether snowpack
* exists or not:
* --------------------------------------------------------------------*/
    *esnow = 0.0;
    if (*sneqv == 0.0)
    {
#ifdef _NOAH_
        NoPac (etp, eta, prcp, pcpdrp, smc, smcmax, smcmin, smcwlt, smcref,
           smcdry, cmc, cmcmax, nsoil, dt, shdfac, sbeta, q2, t1, sfctmp, t24,
           th2, fdown, f1, emissi, ssoil, stc, epsca, vgalpha, vgbeta,
           macksat, areaf, nmacd, mac_status, nwtbl, pc, rch, rr, cfactr, sh2o, frzx,
           zsoil, dksat, tbot, zbot, infil, runoff2, runoff3, edir, ec, et, ett,
           nroot, rtdis, quartz, fxexp, csoil, beta, drip, dew, flx1, flx3,
           vegtyp, isurban);
#else
        NoPac (etp, eta, prcp, smc, smcmax, smcwlt, smcref, smcdry, cmc,
           cmcmax, nsoil, dt, shdfac, sbeta, q2, t1, sfctmp, t24, th2, fdown,
           f1, emissi, ssoil, stc, epsca, bexp, pc, rch, rr, cfactr, sh2o,
           slope, kdt, frzx, psisat, zsoil, dksat, dwsat, tbot, zbot, runoff1,
           runoff2, runoff3, edir, ec, et, ett, nroot, rtdis, quartz, fxexp,
           csoil, beta, drip, dew, flx1, flx3, vegtyp, isurban);
#endif
        *eta_kinematic = *eta;
    }
    else
    {
#ifdef _NOAH_
        SnoPac (etp, eta, prcp, prcpf, pcpdrp, snowng, smc, smcmax, smcmin,
           smcwlt, smcref, smcdry, cmc, cmcmax, nsoil, dt, sbeta, df1, q2, t1,
           sfctmp, t24, th2, fdown, f1, ssoil, stc, epsca, sfcprs, vgalpha,
           vgbeta, macksat, areaf, nmacd, mac_status, nwtbl, pc, rch, rr, cfactr, sncovr,
           sneqv, sndens, snowh, sh2o, frzx, zsoil, dksat, tbot, zbot, shdfac,
           infil, runoff2, runoff3, edir, ec, et, ett, nroot, snomlt, rtdis,
           quartz, fxexp, csoil, beta, drip, dew, flx1, flx2, flx3, esnow,
           etns, emissi, ribb, soldn, isurban, vegtyp);
#else
        SnoPac (etp, eta, prcp, prcpf, snowng, smc, smcmax, smcwlt, smcref,
           smcdry, cmc, cmcmax, nsoil, dt, sbeta, df1, q2, t1, sfctmp, t24,
           th2, fdown, f1, ssoil, stc, epsca, sfcprs, bexp, pc, rch, rr,
           cfactr, sncovr, sneqv, sndens, snowh, sh2o, slope, kdt, frzx,
           psisat, zsoil, dwsat, dksat, tbot, zbot, shdfac, runoff1, runoff2,
           runoff3, edir, ec, et, ett, nroot, snomlt, rtdis, quartz, fxexp,
           csoil, beta, drip, dew, flx1, flx2, flx3, esnow, etns, emissi,
           ribb, soldn, isurban, vegtyp);
#endif
        *eta_kinematic = *esnow + *etns;
    }

    /*
     * calculate effective mixing ratio at grnd level (skin) 
     */

    //  *q1 = *q2 + *eta * CP / *rch;
    *q1 = *q2 + *eta_kinematic * CP / *rch;

/*----------------------------------------------------------------------
* determine sensible heat (h) in energy units (w m-2)
* --------------------------------------------------------------------*/

    *sheat = -(*ch * CP * *sfcprs) / (RD * *t2v) * (*th2 - *t1);

/*----------------------------------------------------------------------
* convert evap terms from kinematic (kg m-2 s-1) to energy units (w m-2)
* --------------------------------------------------------------------*/
    *edir = *edir * LVH2O;
    *ec = *ec * LVH2O;
    for (k = 0; k < *nsoil; k++)
        et[k] = et[k] * LVH2O;
    *ett = *ett * LVH2O;
    *esnow = *esnow * LSUBS;
    *etp = *etp * ((1. - *sncovr) * LVH2O + *sncovr * LSUBS);
    if (*etp > 0.)
        *eta = *edir + *ec + *ett + *esnow;
    else
        *eta = *etp;

/*----------------------------------------------------------------------
* determine beta (ratio of actual to potential evap)
* --------------------------------------------------------------------*/
    if (*etp == 0.0)
        *beta = 0.0;
    else
        *beta = *eta / *etp;

/*----------------------------------------------------------------------
* convert the sign of soil heat flux so that:
*   ssoil>0: warm the surface  (night time)
*   ssoil<0: cool the surface  (day time)
* --------------------------------------------------------------------*/
    *ssoil = -1.0 * *ssoil;

/*----------------------------------------------------------------------
*  for the case of land:
*  convert runoff3 (internal layer runoff from supersat) from m to m s-1
*  and add to subsurface runoff/drainage/baseflow.  runoff2 is already
*  a rate at this point
* --------------------------------------------------------------------*/

    *runoff3 = *runoff3 / *dt;
#ifndef _NOAH_
    *runoff2 = *runoff2 + *runoff3;
#endif
    /* definitions of soilm and soilw have been changed in flux-pihm for
     * coupling purpose */
    *soilm = -1.0 * sh2o[0] * zsoil[0];
    for (k = 1; k < *nsoil; k++)
        *soilm += sh2o[k] * (zsoil[k - 1] - zsoil[k]);

    *soilw = -1.0 * sh2o[0] * zsoil[0];
    for (k = 1; k < *nroot; k++)
        *soilw += sh2o[k] * (zsoil[k - 1] - zsoil[k]);
    *soilw /= -zsoil[*nroot - 1];

    //*soilm = -1.0 * smc[0] * zsoil[0];
    //for (k = 1; k < *nsoil; k++)
    //    *soilm = *soilm + smc[k] * (zsoil[k - 1] - zsoil[k]);
    //soilwm = -1.0 * (*smcmax - *smcwlt) * zsoil[0];
    //soilww = -1.0 * (smc[0] - *smcwlt) * zsoil[0];

    //for (k = 0; k < *nsoil; k++)
    //    smav[k] = (smc[k] - *smcwlt) / (*smcmax - *smcwlt);

    //if (*nroot > 0)
    //{
    //    for (k = 1; k < *nroot; k++)
    //    {
    //        soilwm = soilwm + (*smcmax - *smcwlt) * (zsoil[k - 1] - zsoil[k]);
    //        soilww = soilww + (smc[k] - *smcwlt) * (zsoil[k - 1] - zsoil[k]);
    //    }
    //}
    //if (soilwm < 1.e-6)
    //{
    //    soilwm = 0.0;
    //    *soilw = 0.0;
    //    *soilm = 0.0;
    //}
    //else
    //    *soilw = soilww / soilwm;

    /*
     * #ifdef _debug_
     * printf("prcp = %f, q2 = %f, q2sat = %f, sfcprs = %f, sfcspd = %f, sfctmp = %f, solnet = %f, zlvl = %f\n", *prcp, *q2, *q2sat, *sfcprs, *sfcspd, *sfctmp, *solnet, *zlvl);
     * printf("czil = %f, ch = %f, cm = %f\n", *czil, *ch, *cm);
     * printf("xlai = %f, shdfac = %f, z0 = %f\n", *xlai, *shdfac, *z0);
     * printf("ec = %f, edir = %f, ett = %f, etp = %f\n", *ec, *edir, *ett, *etp);
     * printf("pc = %f\n", *pc);
     * #endif
     */
    free (frzgra);
    free (snowng);
    free (df1);
    free (epsca);
    free (rch);
    free (rr);
    free (sndens);
    free (sncond);
    free (sn_new);
    free (t24);
    free (t2v);
    free (prcpf);
    free (etns);

/*----------------------------------------------------------------------
  end subroutine sflx
* --------------------------------------------------------------------*/
}

void
AlCalc (double *alb, double *snoalb, double *embrd, double *shdfac,
   double *shdmin, double *sncovr, double *tsnow, double *albedo,
   double *emissi, double *dt, int *snowng, double *snotime1, double *lvcoef)
{

/*----------------------------------------------------------------------
* calculate albedo including snow effect (0 -> 1)
*   alb     snowfree albedo
*   snoalb  maximum (deep) snow albedo
*   shdfac    areal fractional coverage of green vegetation
*   shdmin    minimum areal fractional coverage of green vegetation
*   sncovr  fractional snow cover
*   albedo  surface albedo including snow effect
*   tsnow   snow surface temperature (k)
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* snoalb is argument representing maximum albedo over deep snow,
* as passed into sflx, and adapted from the satellite-based maximum
* snow albedo fields provided by d. robinson and g. kukla
* (1985, jcam, vol 24, 402-411)
* --------------------------------------------------------------------*/
    double          snoalb2;
    double          snoalb1;
    double          snacca = 0.94, snaccb = 0.58;

    /*
     * turn of vegetation effect 
     */
    //      albedo = alb + (1.0- (shdfac - shdmin))* sncovr * (snoalb - alb)
    //      albedo = (1.0-sncovr)*alb + sncovr*snoalb !this is equivalent to below

    *albedo = *alb + *sncovr * (*snoalb - *alb);
    *emissi = *embrd + *sncovr * (EMISSI_S - *embrd);

    /*
     * base formulation (dickinson et al., 1986, cogley et al., 1990) 
     */

    /*
     * if (tsnow.le.263.16) then
     * albedo=snoalb
     * else
     * if (tsnow.lt.273.16) then
     * tm=0.1*(tsnow-263.16)
     * snoalb1=0.5*((0.9-0.2*(tm**3))+(0.8-0.16*(tm**3)))
     * else
     * snoalb1=0.67
     * if(sncovr.gt.0.95) snoalb1= 0.6
     * snoalb1 = alb + sncovr*(snoalb-alb)
     * endif
     * endif
     * albedo = alb + sncovr*(snoalb1-alb)
     * 
     * isba formulation (verseghy, 1991; baker et al., 1990)
     * snoalb1 = snoalb+coef*(0.85-snoalb)
     * snoalb2=snoalb1
     * m          lstsnw=lstsnw+1
     * snotime1 = snotime1 + dt
     * if (snowng) then
     * snoalb2=snoalb
     * *m             lstsnw=0
     * snotime1 = 0.0
     * else
     * if (tsnow.lt.273.16) then
     * *              snoalb2=snoalb-0.008*lstsnw*dt/86400
     * *m              snoalb2=snoalb-0.008*snotime1/86400
     * snoalb2=(snoalb2-0.65)*exp(-0.05*dt/3600)+0.65
     * *              snoalb2=(albedo-0.65)*exp(-0.01*dt/3600)+0.65
     * else
     * snoalb2=(snoalb2-0.5)*exp(-0.0005*dt/3600)+0.5
     * *              snoalb2=(snoalb-0.5)*exp(-0.24*lstsnw*dt/86400)+0.5
     * *m              snoalb2=(snoalb-0.5)*exp(-0.24*snotime1/86400)+0.5
     * endif
     * endif
     * 
     * *               print*,'snoalb2',snoalb2,'albedo',albedo,'dt',dt
     * albedo = alb + sncovr*(snoalb2-alb)
     * if (albedo .gt. snoalb2) albedo=snoalb2
     * *m          lstsnw1=lstsnw
     * *          snotime = snotime1
     */

    /*
     * formulation by livneh
     * * ----------------------------------------------------------------------
     * * snoalb is considered as the maximum snow albedo for new snow, at
     * * a value of 85%. snow albedo curve defaults are from bras p.263. should
     * * not be changed except for serious problems with snow melt.
     * * to implement accumulatin parameters, snacca and snaccb, assert that it
     * * is indeed accumulation season. i.e. that snow surface temp is below
     * * zero and the date falls between october and february
     * * --------------------------------------------------------------------
     */
    snoalb1 = *snoalb + *lvcoef * (0.85 - *snoalb);
    snoalb2 = snoalb1;

/*---------------- initial lstsnw ------------------------------------*/
    if (*snowng)
        *snotime1 = 0.0;
    else
    {
        *snotime1 = *snotime1 + *dt;
        //      if (tsnow.lt.273.16) then
        snoalb2 = snoalb1 * pow (snacca, pow (*snotime1 / 86400.0, snaccb));
        //      else
        //          snoalb2 =snoalb1*(snthwa**((snotime1/86400.0)**snthwb))
        //               endif
    }

    snoalb2 = snoalb2 > *alb ? snoalb2 : *alb;
    *albedo = *alb + *sncovr * (snoalb2 - *alb);
    if (*albedo > snoalb2)
        *albedo = snoalb2;

    //          if (tsnow.lt.273.16) then
    //            albedo=snoalb-0.008*dt/86400
    //          else
    //            albedo=(snoalb-0.5)*exp(-0.24*dt/86400)+0.5
    //          endif

    //      if (albedo > snoalb) albedo = snoalb

/*----------------------------------------------------------------------
  end subroutine AlCalc
* --------------------------------------------------------------------*/
}

void
CanRes (double *solar, double *ch, double *sfctmp, double *q2, double *sfcprs,
   double *smc, double *zsoil, int *nsoil, double *smcwlt, double *smcref,
   double *rsmin, double *rc, double *pc, int *nroot, double *q2sat,
   double *dqsdt2, double *topt, double *rsmax, double *rgl, double *hs,
   double *xlai, double *rcs, double *rct, double *rcq, double *rcsoil,
   double *emissi)
{

/*----------------------------------------------------------------------
* subroutine CanRes
* ----------------------------------------------------------------------
* calculate canopy resistance which depends on incoming solar radiation,
* air temperature, atmospheric water vapor pressure deficit at the
* lowest model level, and soil moisture (preferably unfrozen soil
* moisture rather than total)
* ----------------------------------------------------------------------
* source:  jarvis (1976), noilhan and planton (1989, mwr), jacquemin and
* noilhan (1990, blm)
* see also:  chen et al (1996, jgr, vol 101(d3), 7251-7268), eqns 12-14
* and table 2 of sec. 3.1.2
* ----------------------------------------------------------------------
* input:
*   solar   incoming solar radiation
*   ch      surface exchange coefficient for heat and moisture
*   sfctmp  air temperature at 1st level above ground
*   q2      air humidity at 1st level above ground
*   q2sat   saturation air humidity at 1st level above ground
*   dqsdt2  slope of saturation humidity function wrt temp
*   sfcprs  surface pressure
*   smc     volumetric soil moisture
*   zsoil   soil depth (negative sign, as it is below ground)
*   nsoil   no. of soil layers
*   nroot   no. of soil layers in root zone (1.le.nroot.le.nsoil)
*   xlai    leaf area index
*   smcwlt  wilting point
*   smcref  reference soil moisture (where soil water deficit stress
*             sets in)
* rsmin, rsmax, topt, rgl, hs are canopy stress parameters set in
*   surboutine RedPrm
* output:
*   pc  plant coefficient
*   rc  canopy resistance
* --------------------------------------------------------------------*/

    int             k;
    double          delta, ff, gx, rr;
    double          part[*nsoil];
    double          slv = 2.501000e6;

/*----------------------------------------------------------------------
* initialize canopy resistance multiplier terms.
* --------------------------------------------------------------------*/
    *rcs = 0.0;
    *rct = 0.0;
    *rcq = 0.0;
    *rcsoil = 0.0;

/*----------------------------------------------------------------------
* contribution due to incoming solar radiation
* --------------------------------------------------------------------*/
    *rc = 0.0;
    ff = 0.55 * 2.0 * *solar / (*rgl * *xlai);
    *rcs = (ff + *rsmin / *rsmax) / (1.0 + ff);

//    printf("solar = %lf, rgl = %lf, xlai = %lf, ff = %lf, rsmin = %lf, rsmax = %lf, rcs = %lf\n", *solar, *rgl, *xlai, ff, *rsmin, *rsmax, *rcs);

/*----------------------------------------------------------------------
* contribution due to air temperature at first model level above ground
* rct expression from noilhan and planton (1989, mwr).
* --------------------------------------------------------------------*/
    *rcs = *rcs > 0.0001 ? *rcs : 0.0001;
    *rct = 1.0 - 0.0016 * pow (*topt - *sfctmp, 2.0);

/*----------------------------------------------------------------------
* contribution due to vapor pressure deficit at first model level.
* rcq expression from ssib
* --------------------------------------------------------------------*/
    *rct = *rct > 0.0001 ? *rct : 0.0001;
    *rcq = 1.0 / (1.0 + *hs * (*q2sat - *q2));

/*----------------------------------------------------------------------
* contribution due to soil moisture availability.
* determine contribution from each soil layer, then add them up.
* --------------------------------------------------------------------*/
    *rcq = *rcq > 0.01 ? *rcq : 0.01;
    gx = (smc[0] - *smcwlt) / (*smcref - *smcwlt);
    if (gx > 1.)
        gx = 1.;
    if (gx < 0.)
        gx = 0.;

/*----------------------------------------------------------------------
* use soil depth as weighting factor
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* use root distribution as weighting factor
*      part(1) = rtdis(1) * gx
* --------------------------------------------------------------------*/
    part[0] = (zsoil[0] / zsoil[*nroot - 1]) * gx;
    for (k = 1; k < *nroot; k++)
    {
        gx = (smc[k] - *smcwlt) / (*smcref - *smcwlt);
        if (gx > 1.)
            gx = 1.;
        if (gx < 0.)
            gx = 0.;

/*----------------------------------------------------------------------
* use soil depth as weighting factor
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* use root distribution as weighting factor
*        part(k) = rtdis(k) * gx
* --------------------------------------------------------------------*/
        part[k] = ((zsoil[k] - zsoil[k - 1]) / zsoil[*nroot - 1]) * gx;
    }
    for (k = 0; k < *nroot; k++)
        *rcsoil = *rcsoil + part[k];

/*----------------------------------------------------------------------
* determine canopy resistance due to all factors.  convert canopy
* resistance (rc) to plant coefficient (pc) to be used with potential
* evap in determining actual evap.  pc is determined by:
*   pc * linerized Penman potential evap =
*   Penman-monteith actual evaporation (containing rc term).
* --------------------------------------------------------------------*/
    *rcsoil = *rcsoil > 0.0001 ? *rcsoil : 0.0001;

    *rc = *rsmin / (*xlai * *rcs * *rct * *rcq * *rcsoil);
    //  rr = (4.* SIGMA * RD / CP)* (sfctmp **4.)/ (sfcprs * ch) + 1.0;
    rr =
       (4. * *emissi * SIGMA * RD / CP) * pow (*sfctmp,
       4.) / (*sfcprs * *ch) + 1.0;

    delta = (slv / CP) * *dqsdt2;

    *pc = (rr + delta) / (rr * (1. + *rc * *ch) + delta);

/*----------------------------------------------------------------------
  end subroutine CanRes
* --------------------------------------------------------------------*/
}

void CSnow (double *sncond, double *dsnow)
{

/*----------------------------------------------------------------------
* subroutine CSnow
* function CSnow
* ----------------------------------------------------------------------
* calculate snow thermal conductivity
* --------------------------------------------------------------------*/
    double          c;
    double          unit = 0.11631;

/*----------------------------------------------------------------------
* sncond in units of cal/(cm*hr*c), returned in w/(m*c)
* csnow in units of cal/(cm*hr*c), returned in w/(m*c)
* basic version is dyachkova equation (1960), for range 0.1-0.4
* --------------------------------------------------------------------*/
    c = 0.328 * pow (10, 2.25 * *dsnow);
    //  csnow=unit*c

/*----------------------------------------------------------------------
* de vaux equation (1933), in range 0.1-0.6
* ----------------------------------------------------------------------
*      sncond=0.0293*(1.+100.*dsnow**2)
*      csnow=0.0293*(1.+100.*dsnow**2)

* ----------------------------------------------------------------------
* e. andersen from flerchinger
* ----------------------------------------------------------------------
*      sncond=0.021+2.51*dsnow**2
*      csnow=0.021+2.51*dsnow**2

*      sncond = unit * c
* double snow thermal conductivity
*/
    *sncond = 2.0 * unit * c;

/*----------------------------------------------------------------------
  end subroutine CSnow
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
DEvap (double *edir, double *etp1, double *smc, double *zsoil, double *shdfac,
   double *smcmax, double *dksat, double *smcdry, double *smcref,
   double *smcwlt, double *fxexp)
#else
void
DEvap (double *edir, double *etp1, double *smc, double *zsoil, double *shdfac,
   double *smcmax, double *bexp, double *dksat, double *dwsat, double *smcdry,
   double *smcref, double *smcwlt, double *fxexp)
#endif
{

/*----------------------------------------------------------------------
* subroutine DEvap
* function DEvap
* ----------------------------------------------------------------------
* calculate direct soil evaporation
* --------------------------------------------------------------------*/
    double          fx, sratio;

/*----------------------------------------------------------------------
* direct evap a function of relative soil moisture availability, linear
* when fxexp=1.
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* fx > 1 represents demand control
* fx < 1 represents flux control
* --------------------------------------------------------------------*/

    sratio = (*smc - *smcdry) / (*smcmax - *smcdry);
    if (sratio > 0.)
    {
        fx = pow (sratio, *fxexp);
        fx = fx > 1. ? 1. : fx;
        fx = fx < 0. ? 0. : fx;
    }
    else
        fx = 0.;

/*----------------------------------------------------------------------
* allow for the direct-evap-reducing effect of shade
* --------------------------------------------------------------------*/
    *edir = fx * (1.0 - *shdfac) * *etp1;

/*----------------------------------------------------------------------
  end subroutine DEvap
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
Evapo (double *eta1, double *smc, int *nsoil, double *cmc, double *etp1,
   double *dt, double *zsoil, double *sh2o, double *smcmax, double *pc,
   double *smcwlt, double *dksat, double *smcref, double *shdfac,
   double *cmcmax, double *smcdry, double *cfactr, double *edir, double *ec,
   double *et, double *ett, double *sfctmp, double *q2, int *nroot,
   double *rtdis, double *fxexp)
#else
void
Evapo (double *eta1, double *smc, int *nsoil, double *cmc, double *etp1,
   double *dt, double *zsoil, double *sh2o, double *smcmax, double *bexp,
   double *pc, double *smcwlt, double *dksat, double *dwsat, double *smcref,
   double *shdfac, double *cmcmax, double *smcdry, double *cfactr,
   double *edir, double *ec, double *et, double *ett, double *sfctmp,
   double *q2, int *nroot, double *rtdis, double *fxexp)
#endif
{

/*----------------------------------------------------------------------
* subroutine Evapo
* ----------------------------------------------------------------------
* calculate soil moisture flux.  the soil moisture content (smc - a per
* unit volume measurement) is a dependent variable that is updated with
* prognostic eqns. the canopy moisture content (cmc) is also updated.
* frozen ground version:  new states added: sh2o, and frozen ground
* correction factor, frzfact and parameter slope.
* --------------------------------------------------------------------*/
    int             k;
    double          cmc2ms;

/*----------------------------------------------------------------------
* executable code begins here if the potential evapotranspiration is
* greater than zero.
* --------------------------------------------------------------------*/
    *edir = 0.;
    *ec = 0.;
    *ett = 0.;
    for (k = 0; k < *nsoil; k++)
        et[k] = 0.;

/*----------------------------------------------------------------------
* retrieve direct evaporation from soil surface.  call this function
* only if veg cover not complete.
* frozen ground version:  sh2o states replace smc states.
* --------------------------------------------------------------------*/
    if (*etp1 > 0.0)
    {
        if (*shdfac < 1.)
#ifdef _NOAH_
            DEvap (edir, etp1, smc, zsoil, shdfac, smcmax, dksat, smcdry,
               smcref, smcwlt, fxexp);
#else
            DEvap (edir, etp1, smc, zsoil, shdfac, smcmax, bexp, dksat, dwsat,
               smcdry, smcref, smcwlt, fxexp);
#endif

/*----------------------------------------------------------------------
* initialize plant total transpiration, retrieve plant transpiration,
* and accumulate it for all soil layers.
* --------------------------------------------------------------------*/

        if (*shdfac > 0.0)
        {
            Transp (et, nsoil, etp1, sh2o, cmc, zsoil, shdfac, smcwlt, cmcmax,
               pc, cfactr, smcref, sfctmp, q2, nroot, rtdis);
            for (k = 0; k < *nsoil; k++)
                *ett = *ett + et[k];

/*----------------------------------------------------------------------
* calculate canopy evaporation.
* if statements to avoid tangent linear problems near cmc=0.0.
* --------------------------------------------------------------------*/
            if (*cmc > 0.0)
                *ec =
                   *shdfac * pow ((*cmc / *cmcmax > 1. ? 1. : *cmc / *cmcmax),
                   *cfactr) * *etp1;
            else
                *ec = 0.0;

/*----------------------------------------------------------------------
* ec should be limited by the total amount of available water on the
* canopy.  -f.chen, 18-oct-1994
* --------------------------------------------------------------------*/
            cmc2ms = *cmc / *dt;
            *ec = cmc2ms < *ec ? cmc2ms : *ec;
        }
    }

/*----------------------------------------------------------------------
* total up evap and transp types to obtain actual evapotransp
* --------------------------------------------------------------------*/
    *eta1 = *edir + *ett + *ec;

/*----------------------------------------------------------------------
  end subroutine Evapo
* --------------------------------------------------------------------*/
}

void Fac2Mit (double *smcmax, double *flimit)
{
    *flimit = 0.90;

    if (*smcmax == 0.395)
        *flimit = 0.59;
    else if ((*smcmax == 0.434) || (*smcmax == 0.404))
        *flimit = 0.85;
    else if ((*smcmax == 0.465) || (*smcmax == 0.406))
        *flimit = 0.86;
    else if ((*smcmax == 0.476) || (*smcmax == 0.439))
        *flimit = 0.74;
    else if ((*smcmax == 0.200) || (*smcmax == 0.464))
        *flimit = 0.80;

/*----------------------------------------------------------------------
  end subroutine Fac2Mit
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
FrH2O (double *freew, double *tkelv, double *smc, double *sh2o, double *smcmax,
   double *smcmin, double *vgalpha, double *vgbeta)
#else
void
FrH2O (double *freew, double *tkelv, double *smc, double *sh2o, double *smcmax,
   double *bexp, double *psis)
#endif
{

/*----------------------------------------------------------------------
* subroutine FrH2O
* ----------------------------------------------------------------------
* calculate amount of supercooled liquid soil water content if
* temperature is below 273.15k (t0).  requires newton-type iteration to
* solve the nonlinear implicit equation given in eqn 17 of koren et al
* (1999, jgr, vol 104(d16), 19569-19585).
* ----------------------------------------------------------------------
* new version (june 2001): much faster and more accurate newton
* iteration achieved by first taking log of eqn cited above -- less than
* 4 (typically 1 or 2) iterations achieves convergence.  also, explicit
* 1-step solution option for special case of parameter ck=0, which
* reduces the original implicit equation to a simpler explicit form,
* known as the "flerchinger eqn". improved handling of solution in the
* limit of freezing point temperature t0.
* ----------------------------------------------------------------------
* input:

*   tkelv.........temperature (kelvin)
*   smc...........total soil moisture content (volumetric)
*   sh2o..........liquid soil moisture content (volumetric)
*   smcmax........saturation soil moisture content (from RedPrm)
*   b.............soil type "b" parameter (from RedPrm)
*   psis..........saturated soil matric potential (from RedPrm)

* output:
*   freew..........supercooled liquid water content
* --------------------------------------------------------------------*/
    double          denom, df, dswl, fk, swl, swlk;
    int             nlog, kcount;
    //      parameter(ck = 0.0)
    double          ck = 8.0, error = 0.005, hlice = 3.335e5, gs = 9.81, t0 = 273.15;

#ifdef _NOAH_
    double          mx;
    mx = *vgbeta / (1 - *vgbeta);
#else

/*----------------------------------------------------------------------
* limits on parameter b: b < 5.5  (use parameter blim)
* simulations showed if b > 5.5 unfrozen water content is
* non-realistically high at very low temperatures.
* --------------------------------------------------------------------*/
    bx = *bexp;

/*----------------------------------------------------------------------
* initializing iterations counter and iterative solution flag.
* --------------------------------------------------------------------*/
    if (*bexp > blim)
        bx = blim;
#endif
    nlog = 0;

/*----------------------------------------------------------------------
*  if temperature not significantly below freezing (t0), sh2o = smc
* --------------------------------------------------------------------*/
    kcount = 0;
    //      frh2o = smc
    if (*tkelv > (t0 - 1.e-3))
        *freew = *smc;
    else
    {

/*----------------------------------------------------------------------
* option 1: iterated solution for nonzero ck
* in koren et al, jgr, 1999, eqn 17
* ----------------------------------------------------------------------
* initial guess for swl (frozen content)
* --------------------------------------------------------------------*/
        if (ck != 0.0)
        {
            swl = *smc - *sh2o;

/*----------------------------------------------------------------------
* keep within bounds.
* --------------------------------------------------------------------*/
#ifdef _NOAH_
            if (swl > (*smc - *smcmin - 0.02))
                swl = *smc - *smcmin - 0.02;
#else
            if (swl > (*smc - 0.02))
                swl = *smc - 0.02;
#endif

/*----------------------------------------------------------------------
*  start of iterations
* --------------------------------------------------------------------*/
            if (swl < 0.)
                swl = 0.;
          c1001:
            if (!((nlog < 10) && (kcount == 0)))
                goto c1002;
            nlog = nlog + 1;
#ifdef _NOAH_
            df =
               log ((gs / *vgalpha / hlice) * pow (1. + ck * swl,
                  2.) * pow (pow ((*smc - swl - *smcmin) / (*smcmax -
                        *smcmin), mx) - 1.,
                  1. / *vgbeta)) - log (-(*tkelv - t0) / *tkelv);
            denom =
               2. * ck / (1. + ck * swl) - 1. / (1 - *vgbeta) / (*smcmax -
               *smcmin) * pow ((*smc - swl - *smcmin) / (*smcmax - *smcmin),
               mx - 1.) / (pow ((*smc - swl - *smcmin) / (*smcmax - *smcmin),
                  mx) - 1.);
#else
            df =
               log ((*psis * gs / hlice) * pow (1. + ck * swl,
                  2.) * pow (*smcmax / (*smc - swl),
                  bx)) - log (-(*tkelv - t0) / *tkelv);
            denom = 2. * ck / (1. + ck * swl) + bx / (*smc - swl);
#endif
            swlk = swl - df / denom;

/*----------------------------------------------------------------------
* bounds useful for mathematical solution.
* --------------------------------------------------------------------*/
#ifdef _NOAH_
            if (swlk > (*smc - *smcmin - 0.02))
                swlk = *smc - *smcmin - 0.02;
#else
            if (swlk > (*smc - 0.02))
                swlk = *smc - 0.02;
#endif
            if (swlk < 0.)
                swlk = 0.;

/*----------------------------------------------------------------------
* mathematical solution bounds applied.
* --------------------------------------------------------------------*/
            dswl = fabs (swlk - swl);

/*----------------------------------------------------------------------
* if more than 10 iterations, use explicit method (ck=0 approx.)
* when dswl less or eq. error, no more iterations required.
* --------------------------------------------------------------------*/
            swl = swlk;
            if (dswl <= error)
                kcount = kcount + 1;

/*----------------------------------------------------------------------
*  end of iterations
* ----------------------------------------------------------------------
* bounds applied within do-block are valid for physical solution.
* --------------------------------------------------------------------*/
            //          frh2o = smc - swl
            goto c1001;
          c1002:
            *freew = *smc - swl;
        }

/*----------------------------------------------------------------------
* end option 1
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* option 2: explicit solution for flerchinger eq. i.e. ck=0
* in koren et al., jgr, 1999, eqn 17
* apply physical bounds to flerchinger solution
* --------------------------------------------------------------------*/
        if (kcount == 0)
        {
            //          printf ("flerchinger used in new version. iterations= %d\n",
            //             nlog);
#ifdef _NOAH_
            fk =
               pow (pow (-(*tkelv - t0) / *tkelv * *vgalpha * hlice / gs,
                  *vgbeta) + 1., 1. / mx) * (*smcmax - *smcmin);
#else
            fk =
               pow ((hlice / (gs * (-*psis))) * ((*tkelv - t0) / *tkelv),
               -1 / bx) * *smcmax;
#endif
            //          frh2o = min (fk, smc)
            if (fk < 0.02)
                fk = 0.02;
            *freew = fk < *smc ? fk : *smc;

/*----------------------------------------------------------------------
* end option 2
* --------------------------------------------------------------------*/
        }
    }

/*----------------------------------------------------------------------
  end subroutine FrH2O
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
HRT (double *rhsts, double *stc, double *smc, double *smcmax, double *smcmin,
   int *nsoil, double *zsoil, double *yy, double *zz1, double *tbot,
   double *zbot, double *sh2o, double *dt, double *vgalpha, double *vgbeta,
   double *f1, double *df1, double *quartz, double *csoil, double *ai,
   double *bi, double *ci, int *vegtyp, int *isurban)
#else
void
HRT (double *rhsts, double *stc, double *smc, double *smcmax, int *nsoil,
   double *zsoil, double *yy, double *zz1, double *tbot, double *zbot,
   double *psisat, double *sh2o, double *dt, double *bexp, double *f1,
   double *df1, double *quartz, double *csoil, double *ai, double *bi,
   double *ci, int *vegtyp, int *isurban)
#endif
{

/*----------------------------------------------------------------------
* subroutine HRT
* ----------------------------------------------------------------------
* calculate the right hand side of the time tendency term of the soil
* thermal diffusion equation.  also to compute ( prepare ) the matrix
* coefficients for the tri-diagonal matrix of the implicit time scheme.
* --------------------------------------------------------------------*/
    int             itavg;
    int             k;

    double          ddz, ddz2, denom, df1k, dtsdz, dtsdz2, hcpct, ssoil, sice,
       csoil_loc;
    double         *df1n, *qtot, *tavg, *tbk, *tbk1, *tsnsr, *tsurf;
    double          t0 = 273.15, cair = 1004.0, cice = 2.106e6, ch2o = 4.2e6;

    df1n = (double *)malloc (sizeof (double));
    qtot = (double *)malloc (sizeof (double));
    tavg = (double *)malloc (sizeof (double));
    tbk = (double *)malloc (sizeof (double));
    tbk1 = (double *)malloc (sizeof (double));
    tsnsr = (double *)malloc (sizeof (double));
    tsurf = (double *)malloc (sizeof (double));

    /*
     * urban 
     */
    if (*vegtyp == *isurban)
        csoil_loc = 3.0e6;
    else
        csoil_loc = *csoil;

/*----------------------------------------------------------------------
* initialize logical for soil layer temperature averaging.
* --------------------------------------------------------------------*/
    itavg = 1;

/*----------------------------------------------------------------------
* begin section for top soil layer
* ----------------------------------------------------------------------
* calc the heat capacity of the top soil layer
* --------------------------------------------------------------------*/
    hcpct =
       sh2o[0] * ch2o + (1.0 - *smcmax) * csoil_loc + (*smcmax -
       smc[0]) * cair + (smc[0] - sh2o[0]) * cice;

/*----------------------------------------------------------------------
* calc the matrix coefficients ai, bi, and ci for the top layer
* --------------------------------------------------------------------*/
    ddz = 1.0 / (-0.5 * zsoil[1]);
    ai[0] = 0.0;
    ci[0] = (*df1 * ddz) / (zsoil[0] * hcpct);

/*----------------------------------------------------------------------
* calculate the vertical soil temp gradient btwn the 1st and 2nd soil
* layers.  then calculate the subsurface heat flux. use the temp
* gradient and subsfc heat flux to calc "right-hand side tendency
* terms", or "rhsts", for top soil layer.
* --------------------------------------------------------------------*/
    bi[0] = -ci[0] + *df1 / (0.5 * zsoil[0] * zsoil[0] * hcpct * *zz1);
    dtsdz = (stc[0] - stc[1]) / (-0.5 * zsoil[1]);
    ssoil = *df1 * (stc[0] - *yy) / (0.5 * zsoil[0] * *zz1);
    //      rhsts[0] = (*df1 * dtsdz - ssoil) / (zsoil[0] * hcpct);
    denom = (zsoil[0] * hcpct);

/*----------------------------------------------------------------------
* next capture the vertical difference of the heat flux at top and
* bottom of first soil layer for use in heat flux constraint applied to
* potential soil freezing/thawing in routine SnkSrc.
* --------------------------------------------------------------------*/
    //  *qtot = ssoil - *df1*dtsdz;
    rhsts[0] = (*df1 * dtsdz - ssoil) / denom;

/*----------------------------------------------------------------------
* calculate frozen water content in 1st soil layer.
* --------------------------------------------------------------------*/
    *qtot = -1.0 * rhsts[0] * denom;

/*----------------------------------------------------------------------
* if temperature averaging invoked (itavg=true; else skip):
* set temp "tsurf" at top of soil column (for use in freezing soil
* physics later in function subroutine SnkSrc).  if snowpack content is
* zero, then tsurf expression below gives tsurf = skin temp.  if
* snowpack is nonzero (hence argument zz1=1), then tsurf expression
* below yields soil column top temperature under snowpack.  then
* calculate temperature at bottom interface of 1st soil layer for use
* later in function subroutine SnkSrc
* --------------------------------------------------------------------*/
    sice = smc[0] - sh2o[0];
    if (itavg)
    {
        *tsurf = (*yy + (*zz1 - 1) * stc[0]) / *zz1;

/*----------------------------------------------------------------------
* if frozen water present or any of layer-1 mid-point or bounding
* interface temperatures below freezing, then call SnkSrc to
* compute heat source/sink (and change in frozen water content)
* due to possible soil water phase change
* --------------------------------------------------------------------*/
        TBnd (stc, stc + 1, zsoil, zbot, 0, nsoil, tbk);
        if ((sice > 0.) || (stc[0] < t0) || (*tsurf < t0) || (*tbk < t0))
        {
            //          *tsnsr = SnkSrc (tavg,smc(1),sh2o(1),
            TmpAvg (tavg, tsurf, stc, tbk, zsoil, nsoil, 0);
#ifdef _NOAH_
            SnkSrc (tsnsr, tavg, smc, sh2o, zsoil, nsoil, smcmax, smcmin,
               vgalpha, vgbeta, dt, 0, qtot);
#else
            SnkSrc (tsnsr, tavg, smc, sh2o, zsoil, nsoil, smcmax, psisat,
               bexp, dt, 0, qtot);
#endif
            //          rhsts(1) = rhsts(1) - *tsnsr / ( zsoil(1) * hcpct )
            rhsts[0] = rhsts[0] - *tsnsr / denom;
        }
    }
    else
    {
        //          *tsnsr = SnkSrc (stc(1),smc(1),sh2o(1),
        if ((sice > 0.) || (stc[0] < t0))
        {
#ifdef _NOAH_
            SnkSrc (tsnsr, stc, smc, sh2o, zsoil, nsoil, smcmax, smcmin,
               vgalpha, vgbeta, dt, 0, qtot);
#else
            SnkSrc (tsnsr, stc, smc, sh2o, zsoil, nsoil, smcmax, psisat, bexp,
               dt, 0, qtot);
#endif
            //          rhsts(1) = rhsts(1) - *tsnsr / ( zsoil(1) * hcpct )
            rhsts[0] = rhsts[0] - *tsnsr / denom;
        }

/*----------------------------------------------------------------------
* this ends section for top soil layer.
* --------------------------------------------------------------------*/
    }

    /*
     * initialize ddz2 
     */

/*--------------------------------------------------------------------*/

    ddz2 = 0.0;
    df1k = *df1;

/*----------------------------------------------------------------------
* loop thru the remaining soil layers, repeating the above process
* (except subsfc or "ground" heat flux not repeated in lower layers)
* ----------------------------------------------------------------------
* calculate heat capacity for this soil layer.
* --------------------------------------------------------------------*/
    for (k = 1; k < *nsoil; k++)
    {
        hcpct =
           sh2o[k] * ch2o + (1.0 - *smcmax) * csoil_loc + (*smcmax -
           smc[k]) * cair + (smc[k] - sh2o[k]) * cice;

/*----------------------------------------------------------------------
* this section for layer 2 or greater, but not last layer.
* ----------------------------------------------------------------------
* calculate thermal diffusivity for this layer.
* --------------------------------------------------------------------*/
        if (k != *nsoil - 1)
        {

/*----------------------------------------------------------------------
* calc the vertical soil temp gradient thru this layer
* --------------------------------------------------------------------*/
#ifdef _NOAH_
            TDfCnd (df1n, smc + k, quartz, smcmax, smcmin, sh2o + k);
#else
            TDfCnd (df1n, smc + k, quartz, smcmax, sh2o + k);
#endif

            /*
             * urban 
             */
            if (*vegtyp == *isurban)
                *df1n = 3.24;

            denom = 0.5 * (zsoil[k - 1] - zsoil[k + 1]);

/*----------------------------------------------------------------------
* calc the matrix coef, ci, after calc'ng its partial product
* --------------------------------------------------------------------*/
            dtsdz2 = (stc[k] - stc[k + 1]) / denom;
            ddz2 = 2. / (zsoil[k - 1] - zsoil[k + 1]);

/*----------------------------------------------------------------------
* if temperature averaging invoked (itavg=true; else skip):  calculate
* temp at bottom of layer.
* --------------------------------------------------------------------*/
            ci[k] = -*df1n * ddz2 / ((zsoil[k - 1] - zsoil[k]) * hcpct);
            if (itavg)
                TBnd (stc + k, stc + k + 1, zsoil, zbot, k, nsoil, tbk1);
        }
        else
        {

/*----------------------------------------------------------------------
* special case of bottom soil layer:  calculate thermal diffusivity for
* bottom layer.
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* calc the vertical soil temp gradient thru bottom layer.
* --------------------------------------------------------------------*/
#ifdef _NOAH_
            TDfCnd (df1n, smc + k, quartz, smcmax, smcmin, sh2o + k);
#else
            TDfCnd (df1n, smc + k, quartz, smcmax, sh2o + k);
#endif

            /*
             * urban 
             */
            if (*vegtyp == *isurban)
                *df1n = 3.24;

            denom = 0.5 * (zsoil[k - 1] + zsoil[k]) - *zbot;

/*----------------------------------------------------------------------
* set matrix coef, ci to zero if bottom layer.
* --------------------------------------------------------------------*/
            dtsdz2 = (stc[k] - *tbot) / denom;

/*----------------------------------------------------------------------
* if temperature averaging invoked (itavg=true; else skip):  calculate
* temp at bottom of last layer.
* --------------------------------------------------------------------*/
            ci[k] = 0.;
            if (itavg)
                TBnd (stc + k, tbot, zsoil, zbot, k, nsoil, tbk1);

/*----------------------------------------------------------------------
* this ends special loop for bottom layer.
* --------------------------------------------------------------------*/
        }

/*----------------------------------------------------------------------
* calculate rhsts for this layer after calc'ng a partial product.
* --------------------------------------------------------------------*/
        denom = (zsoil[k] - zsoil[k - 1]) * hcpct;
        rhsts[k] = (*df1n * dtsdz2 - df1k * dtsdz) / denom;
        *qtot = -1.0 * denom * rhsts[k];

        sice = smc[k] - sh2o[k];

        if (itavg)
        {
            TmpAvg (tavg, tbk, stc + k, tbk1, zsoil, nsoil, k);
            //                  *tsnsr = SnkSrc(tavg,smc(k),sh2o(k),zsoil,nsoil,
            if ((sice > 0.) || (stc[k] < t0) || (*tbk < t0) || (*tbk1 < t0))
            {
#ifdef _NOAH_
                SnkSrc (tsnsr, tavg, smc + k, sh2o + k, zsoil, nsoil, smcmax,
                   smcmin, vgalpha, vgbeta, dt, k, qtot);
#else
                SnkSrc (tsnsr, tavg, smc + k, sh2o + k, zsoil, nsoil, smcmax,
                   psisat, bexp, dt, k, qtot);
#endif
                rhsts[k] = rhsts[k] - *tsnsr / denom;
            }
        }
        else
        {
            //            tsnsr = SnkSrc(stc(k),smc(k),sh2o(k),zsoil,nsoil,
            if ((sice > 0.) || (stc[k] < t0))
            {
#ifdef _NOAH_
                SnkSrc (tsnsr, stc + k, smc + k, sh2o + k, zsoil, nsoil,
                   smcmax, smcmin, vgalpha, vgbeta, dt, k, qtot);
#else
                SnkSrc (tsnsr, stc + k, smc + k, sh2o + k, zsoil, nsoil,
                   smcmax, psisat, bexp, dt, k, qtot);
#endif
                rhsts[k] = rhsts[k] - *tsnsr / denom;
            }
        }

/*----------------------------------------------------------------------
* calc matrix coefs, ai, and bi for this layer.
* --------------------------------------------------------------------*/
        ai[k] = -df1k * ddz / ((zsoil[k - 1] - zsoil[k]) * hcpct);

/*----------------------------------------------------------------------
* reset values of df1, dtsdz, ddz, and tbk for loop to next soil layer.
* --------------------------------------------------------------------*/
        bi[k] = -(ai[k] + ci[k]);
        *tbk = *tbk1;
        df1k = *df1n;
        dtsdz = dtsdz2;
        ddz = ddz2;
    }

    free (df1n);
    free (qtot);
    free (tavg);
    free (tbk);
    free (tbk1);
    free (tsnsr);
    free (tsurf);

/*----------------------------------------------------------------------
  end subroutine HRT
* --------------------------------------------------------------------*/
}


void
HStep (double *stcout, double *stcin, double *rhsts, double *dt, int *nsoil,
   double *ai, double *bi, double *ci)
{

/*----------------------------------------------------------------------
* subroutine HStep
* ----------------------------------------------------------------------
* calculate/update the soil temperature field.
* --------------------------------------------------------------------*/
    int             k;

    double          rhstsin[*nsoil];
    double          ciin[*nsoil];

/*----------------------------------------------------------------------
* create finite difference values for use in Rosr12 routine
* --------------------------------------------------------------------*/
    for (k = 0; k < *nsoil; k++)
    {
        rhsts[k] = rhsts[k] * *dt;
        ai[k] = ai[k] * *dt;
        bi[k] = 1. + bi[k] * *dt;
        ci[k] = ci[k] * *dt;
    }

/*----------------------------------------------------------------------
* copy values for input variables before call to Rosr12
* --------------------------------------------------------------------*/
    for (k = 0; k < *nsoil; k++)
        rhstsin[k] = rhsts[k];
    for (k = 0; k < *nsoil; k++)
        ciin[k] = ci[k];

/*----------------------------------------------------------------------
* solve the tri-diagonal matrix equation
* --------------------------------------------------------------------*/
    Rosr12 (ci, ai, bi, ciin, rhstsin, rhsts, nsoil);

/*----------------------------------------------------------------------
* calc/update the soil temps using matrix solution
* --------------------------------------------------------------------*/
    for (k = 0; k < *nsoil; k++)
        stcout[k] = stcin[k] + ci[k];

/*----------------------------------------------------------------------
* end subroutine HStep
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
NoPac (double *etp, double *eta, double *prcp, double *pcpdrp, double *smc,
   double *smcmax, double *smcmin, double *smcwlt, double *smcref,
   double *smcdry, double *cmc, double *cmcmax, int *nsoil, double *dt,
   double *shdfac, double *sbeta, double *q2, double *t1, double *sfctmp,
   double *t24, double *th2, double *fdown, double *f1, double *emissi,
   double *ssoil, double *stc, double *epsca, double *vgalpha, double *vgbeta,
   double *macksat, double *areaf, int *nmacd, int *mac_status, int *nwtbl, double *pc,
   double *rch, double *rr, double *cfactr, double *sh2o, double *frzfact,
   double *zsoil, double *dksat, double *tbot, double *zbot, double *infil,
   double *runoff2, double *runoff3, double *edir, double *ec, double *et,
   double *ett, int *nroot, double *rtdis, double *quartz, double *fxexp,
   double *csoil, double *beta, double *drip, double *dew, double *flx1,
   double *flx3, int *vegtyp, int *isurban)
#else
void
NoPac (double *etp, double *eta, double *prcp, double *smc, double *smcmax,
   double *smcwlt, double *smcref, double *smcdry, double *cmc,
   double *cmcmax, int *nsoil, double *dt, double *shdfac, double *sbeta,
   double *q2, double *t1, double *sfctmp, double *t24, double *th2,
   double *fdown, double *f1, double *emissi, double *ssoil, double *stc,
   double *epsca, double *bexp, double *pc, double *rch, double *rr,
   double *cfactr, double *sh2o, double *slope, double *kdt, double *frzfact,
   double *psisat, double *zsoil, double *dksat, double *dwsat, double *tbot,
   double *zbot, double *runoff1, double *runoff2, double *runoff3,
   double *edir, double *ec, double *et, double *ett, int *nroot,
   double *rtdis, double *quartz, double *fxexp, double *csoil, double *beta,
   double *drip, double *dew, double *flx1, double *flx3, int *vegtyp,
   int *isurban)
#endif
{

/*----------------------------------------------------------------------
* subroutine NoPac
* ----------------------------------------------------------------------
* calculate soil moisture and heat flux values and update soil moisture
* content and soil heat content values for the case when no snow pack is
* present.
* --------------------------------------------------------------------*/

    int             k;

    double          et1[*nsoil];
    double         *ec1, *edir1, *ett1, *df1, *eta1, *etp1, *prcp1, *yy, *zz1;
    double          yynum;

    ec1 = (double *)malloc (sizeof (double));
    edir1 = (double *)malloc (sizeof (double));
    ett1 = (double *)malloc (sizeof (double));
    df1 = (double *)malloc (sizeof (double));
    eta1 = (double *)malloc (sizeof (double));
    etp1 = (double *)malloc (sizeof (double));
    prcp1 = (double *)malloc (sizeof (double));
    yy = (double *)malloc (sizeof (double));
    zz1 = (double *)malloc (sizeof (double));

/*----------------------------------------------------------------------
* executable code begins here:
* convert etp fnd prcp from kg m-2 s-1 to m s-1 and initialize dew.
* --------------------------------------------------------------------*/
    *prcp1 = *prcp * 0.001;
    *etp1 = *etp * 0.001;
    *dew = 0.0;

/*----------------------------------------------------------------------
* initialize evap terms.
* --------------------------------------------------------------------*/
    *edir = 0.;
    *edir1 = 0.;
    *ec1 = 0.;
    *ec = 0.;
    for (k = 0; k < *nsoil; k++)
    {
        et[k] = 0.;
        et1[k] = 0.;
    }
    *ett = 0.;
    *ett1 = 0.;

    if (*etp > 0.0)
    {
#ifdef _NOAH_
        Evapo (eta1, smc, nsoil, cmc, etp1, dt, zsoil, sh2o, smcmax, pc,
           smcwlt, dksat, smcref, shdfac, cmcmax, smcdry, cfactr, edir1, ec1,
           et1, ett1, sfctmp, q2, nroot, rtdis, fxexp);
        SmFlx (smc, nsoil, cmc, dt, prcp1, pcpdrp, zsoil, sh2o, frzfact,
           smcmax, smcmin, vgalpha, vgbeta, macksat, areaf, nmacd, mac_status, nwtbl,
           smcwlt, dksat, shdfac, cmcmax, infil, runoff2, runoff3, edir1, ec1,
           et1, drip);
#else
        Evapo (eta1, smc, nsoil, cmc, etp1, dt, zsoil, sh2o, smcmax, bexp, pc,
           smcwlt, dksat, dwsat, smcref, shdfac, cmcmax, smcdry, cfactr,
           edir1, ec1, et1, ett1, sfctmp, q2, nroot, rtdis, fxexp);
        SmFlx (smc, nsoil, cmc, dt, prcp1, zsoil, sh2o, slope, kdt, frzfact,
           smcmax, bexp, smcwlt, dksat, dwsat, shdfac, cmcmax, runoff1,
           runoff2, runoff3, edir1, ec1, et1, drip);
#endif

/*----------------------------------------------------------------------
* convert modeled evapotranspiration from  m s-1  to  kg m-2 s-1.
* --------------------------------------------------------------------*/

        *eta = *eta1 * 1000.0;
    }

/*----------------------------------------------------------------------
* if etp < 0, assume dew forms (transform etp1 into dew and reinitialize
* etp1 to zero).
* --------------------------------------------------------------------*/
    else
    {
        *dew = -*etp1;

/*----------------------------------------------------------------------
* convert prcp from 'kg m-2 s-1' to 'm s-1' and add dew amount.
* --------------------------------------------------------------------*/

        *prcp1 = *prcp1 + *dew;
#ifdef _NOAH_
        SmFlx (smc, nsoil, cmc, dt, prcp1, pcpdrp, zsoil, sh2o, frzfact,
           smcmax, smcmin, vgalpha, vgbeta, macksat, areaf, nmacd, mac_status, nwtbl,
           smcwlt, dksat, shdfac, cmcmax, infil, runoff2, runoff3, edir1, ec1,
           et1, drip);
#else
        SmFlx (smc, nsoil, cmc, dt, prcp1, zsoil, sh2o, slope, kdt, frzfact,
           smcmax, bexp, smcwlt, dksat, dwsat, shdfac, cmcmax, runoff1,
           runoff2, runoff3, edir1, ec1, et1, drip);
#endif

/*----------------------------------------------------------------------
* convert modeled evapotranspiration from 'm s-1' to 'kg m-2 s-1'.
* --------------------------------------------------------------------*/
        //      *eta = *eta1 * 1000.0
    }

/*----------------------------------------------------------------------
* based on etp and e values, determine beta
* --------------------------------------------------------------------*/

    if (*etp <= 0.0)
    {
        *beta = 0.0;
        *eta = *etp;
        if (*etp < 0.0)
            *beta = 1.0;
    }
    else
        *beta = *eta / *etp;

/*----------------------------------------------------------------------
* convert modeled evapotranspiration components 'm s-1' to 'kg m-2 s-1'.
* --------------------------------------------------------------------*/
    *edir = *edir1 * 1000.;
    *ec = *ec1 * 1000.;
    for (k = 0; k < *nsoil; k++)
        et[k] = et1[k] * 1000.;
    *ett = *ett1 * 1000.;

/*----------------------------------------------------------------------
* get soil thermal diffuxivity/conductivity for top soil lyr,
* calc. adjusted top lyr soil temp and adjusted soil flux, then
* call ShFlx to compute/update soil heat flux and soil temps.
* --------------------------------------------------------------------*/

#ifdef _NOAH_
    TDfCnd (df1, smc, quartz, smcmax, smcmin, sh2o);
#else
    TDfCnd (df1, smc, quartz, smcmax, sh2o);
#endif

    /*
     * urban 
     */
    if (*vegtyp == *isurban)
        *df1 = 3.24;

/*----------------------------------------------------------------------
* vegetation greenness fraction reduction in subsurface heat flux
* via reduction factor, which is convenient to apply here to thermal
* diffusivity that is later used in HRT to compute sub sfc heat flux
* (see additional comments on veg effect sub-sfc heat flx in
* routine sflx)
* --------------------------------------------------------------------*/
    *df1 = *df1 * exp (*sbeta * *shdfac);

/*----------------------------------------------------------------------
* compute intermediate terms passed to routine HRT (via routine
* ShFlx below) for use in computing subsurface heat flux in HRT
* --------------------------------------------------------------------*/
    yynum = *fdown - *emissi * SIGMA * *t24;
    *yy = *sfctmp + (yynum / *rch + *th2 - *sfctmp - *beta * *epsca) / *rr;

    *zz1 = *df1 / (-0.5 * zsoil[0] * *rch * *rr) + 1.0;

    /*
     * urban 
     */
#ifdef _NOAH_
    ShFlx (ssoil, stc, smc, smcmax, smcmin, nsoil, t1, dt, yy, zz1, zsoil,
       tbot, zbot, smcwlt, sh2o, vgalpha, vgbeta, f1, df1, quartz, csoil,
       vegtyp, isurban);
#else
    ShFlx (ssoil, stc, smc, smcmax, nsoil, t1, dt, yy, zz1, zsoil, tbot, zbot,
       smcwlt, psisat, sh2o, bexp, f1, df1, quartz, csoil, vegtyp, isurban);
#endif

/*----------------------------------------------------------------------
* set flx1 and flx3 (snopack phase change heat fluxes) to zero since
* they are not used here in SnoPac.  flx2 (freezing rain heat flux) was
* similarly initialized in the Penman routine.
* --------------------------------------------------------------------*/
    *flx1 = CPH2O * *prcp * (*t1 - *sfctmp);
    *flx3 = 0.0;

    free (ec1);
    free (edir1);
    free (ett1);
    free (df1);
    free (eta1);
    free (etp1);
    free (prcp1);
    free (yy);
    free (zz1);

/*----------------------------------------------------------------------
  end subroutine NoPac
* --------------------------------------------------------------------*/
}

void Penman (double *sfctmp, double *sfcprs, double *ch, double *t2v, double *th2, double *prcp, double *fdown, double *t24, double *ssoil, double *q2, double *q2sat, double *etp, double *rch, double *epsca, double *rr, int *snowng, int *frzgra, double *dqsdt2, double *flx2, double *emissi_in, double *sneqv, double *t1, double *sncovr)
{

/*----------------------------------------------------------------------
* subroutine Penman
* ----------------------------------------------------------------------
* calculate potential evaporation for the current point.  various
* partial sums/products are also calculated and passed back to the
* calling routine for later use.
* --------------------------------------------------------------------*/

    double          a, delta, fnet, rad, rho, emissi, elcp1, lvs;

    double          elcp = 2.4888e+3, lsubc = 2.501000e+6;

/*----------------------------------------------------------------------
* executable code begins here:
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* prepare partial quantities for Penman equation.
* --------------------------------------------------------------------*/
    emissi = *emissi_in;
    elcp1 = (1.0 - *sncovr) * elcp + *sncovr * elcp * LSUBS / lsubc;
    lvs = (1.0 - *sncovr) * lsubc + *sncovr * LSUBS;

    *flx2 = 0.0;
    //delta = elcp * dqsdt2
    delta = elcp1 * *dqsdt2;
    *t24 = *sfctmp * *sfctmp * *sfctmp * *sfctmp;
    //rr = t24 * 6.48e-8 / (sfcprs * ch) + 1.0
    *rr = emissi * *t24 * 6.48e-8 / (*sfcprs * *ch) + 1.0;
    rho = *sfcprs / (RD * *t2v);

/*----------------------------------------------------------------------
* adjust the partial sums / products with the latent heat
* effects caused by falling precipitation.
* --------------------------------------------------------------------*/
    *rch = rho * CP * *ch;
    if (!*snowng)
    {
        if (*prcp > 0.0)
            *rr = *rr + CPH2O * *prcp / *rch;
    }
    else
        *rr = *rr + CPICE * *prcp / *rch;

/*----------------------------------------------------------------------
* include the latent heat effects of frzng rain converting to ice on
* impact in the calculation of flx2 and fnet.
* --------------------------------------------------------------------*/
    //      fnet = fdown - SIGMA * t24- ssoil
    fnet = *fdown - emissi * SIGMA * *t24 - *ssoil;
    if (*frzgra)
    {
        *flx2 = -LSUBF * (*prcp);
        fnet = fnet - *flx2;

/*----------------------------------------------------------------------
* finish Penman equation calculations.
* --------------------------------------------------------------------*/
    }
    rad = fnet / *rch + *th2 - *sfctmp;
    //  a = elcp * (q2sat - q2)
    a = elcp1 * (*q2sat - *q2);
    *epsca = (a * *rr + rad * delta) / (delta + *rr);
    //  etp = epsca * rch / lsubc;
    *etp = *epsca * *rch / lvs;
#ifdef _DEBUG_
    printf("rr = %lf, rad = %lf, delta = %lf, ch = %lf, epsca = %f, etp = %lg, dqsdt2 = %lf\n", *rr, rad, delta, *ch, *epsca, *etp, *dqsdt2);
#endif

/*----------------------------------------------------------------------
  end subroutine Penman
* --------------------------------------------------------------------*/
}

//void RedPrm (grid_struct * grid, lsm_struct lsm, double *zsoil)
//{
//
///*----------------------------------------------------------------------
//* internally set (default valuess)
//* all soil and vegetation parameters required for the execusion of
//* the noah lsm are defined in vegparm.tbl, soilparm.tb, and genparm.tbl.
//* ----------------------------------------------------------------------
//*     vegetation parameters:
//*             albbrd: sfc background snow-free albedo
//*             cmxtbl: max cnpy capacity
//*              z0brd: background roughness length
//*             shdfac: green vegetation fraction
//*              nroot: rooting depth
//*              rsmin: mimimum stomatal resistance
//*              rsmax: max. stomatal resistance
//*                rgl: parameters used in radiation stress function
//*                 hs: parameter used in vapor pressure deficit functio
//*               topt: optimum transpiration air temperature.
//*             cmcmax: maximum canopy water capacity
//*             cfactr: parameter used in the canopy inteception calculation
//*               snup: threshold snow depth (in water equivalent m) that
//*                     implies 100 percent snow cover
//*                lai: leaf area index
//*
//* ----------------------------------------------------------------------
//*      soil parameters:
//*        smcmax: max soil moisture content (porosity)
//*        smcref: reference soil moisture  (field capacity)
//*        smcwlt: wilting point soil moisture
//*        smcwlt: air dry soil moist content limits
//*       ssatpsi: sat (saturation) soil potential
//*         dksat: sat soil conductivity
//*          bexp: b parameter
//*        ssatdw: sat soil diffusivity
//*           f1: soil thermal diffusivity/conductivity coef.
//*        quartz: soil quartz content
//*  modified by f. chen (12/22/97)  to use the statsgo soil map
//*  modified by f. chen (01/22/00)  to include playa, lava, and white san
//*  modified by f. chen (08/05/02)  to include additional parameters for the noah
//* note: satdw = bb*satdk*(satpsi/maxsmc)
//*         f11 = alog10(satpsi) + bb*alog10(maxsmc) + 2.0
//*       refsmc1=maxsmc*(5.79e-9/satdk)**(1/(2*bb+3)) 5.79e-9 m/s= 0.5 mm
//*       refsmc=refsmc1+1./3.(maxsmc-refsmc1)
//*       wltsmc1=maxsmc*(200./satpsi)**(-1./bb)    (wetzel and chang, 198
//*       wltsmc=wltsmc1-0.5*wltsmc1
//* note: the values for playa is set for it to have a thermal conductivit
//* as sand and to have a hydrulic conductivity as clay
//*
//* ----------------------------------------------------------------------
//* class parameter 'slopetyp' was included to estimate linear reservoir
//* coefficient 'slope' to the baseflow runoff out of the bottom layer.
//* lowest class (slopetyp=0) means highest slope parameter = 1.
//* definition of slopetyp from 'zobler' slope type:
//* slope class  percent slope
//* 1            0-8
//* 2            8-30
//* 3            > 30
//* 4            0-30
//* 5            0-8 & > 30
//* 6            8-30 & > 30
//* 7            0-8, 8-30, > 30
//* 9            glacial ice
//* blank        ocean/sea
//*       slope_data: linear reservoir coefficient
//*       sbeta_data: parameter used to caluculate vegetation effect on soil heat
//*       fxexp_dat:  soil evaporation exponent used in DEvap
//*       csoil_data: soil heat capacity [j m-3 k-1]
//*       salp_data: shape parameter of  distribution function of snow cover
//*       refdk_data and refkdt_data: parameters in the surface runoff parameteriz
//*       frzk_data: frozen ground parameter
//*       zbot_data: depth[m] of lower boundary soil temperature
//*       czil_data: calculate roughness length of heat
//*       smlow_data and mhigh_data: two soil moisture wilt, soil moisture referen
//*                 parameters
//* set maximum number of soil-, veg-, and slopetyp in data statement.
//* --------------------------------------------------------------------*/
//
//    int             i;
//
//    double          frzfact;
//
//    /*
//     * save
//     * * ----------------------------------------------------------------------
//     */
//    if (grid->soiltyp > lsm->soiltbl.slcats)
//    {
//        printf ("warning: too many input soil types\n");
//        exit (0);
//    }
//    if (grid->vegtyp > lsm->vegtbl.lucats)
//    {
//        printf ("warning: too many input landuse types\n");
//        exit (0);
//    }
//    if (grid->slopetyp > lsm->genprmt.slpcats)
//    {
//        printf ("warning: too many input slope types\n");
//        exit (0);
//    }
//
///*----------------------------------------------------------------------
//*  set-up soil parameters
//* --------------------------------------------------------------------*/
//    grid->csoil = lsm->genprmt.csoil_data;
//#ifdef _NOAH_
//    grid->vgalpha = lsm->soiltbl.vga[grid->soiltyp - 1];
//    grid->vgbeta = lsm->soiltbl.vgb[grid->soiltyp - 1];
//    grid->smcmin = lsm->soiltbl.minsmc[grid->soiltyp - 1];
//    grid->macksat = lsm->soiltbl.macksat[grid->soiltyp - 1];
//    grid->areaf = lsm->soiltbl.areaf[grid->soiltyp - 1];
//    grid->nmacd = lsm->soiltbl.nmacd[grid->soiltyp - 1];
//#else
//    grid->bexp = lsm->soiltbl.bb[grid->soiltyp - 1];
//    grid->psisat = lsm->soiltbl.satpsi[grid->soiltyp - 1];
//    grid->dwsat = lsm->soiltbl.satdw[grid->soiltyp - 1];
//#endif
//    grid->dksat = lsm->soiltbl.satdk[grid->soiltyp - 1];
//    grid->f1 = lsm->soiltbl.f11[grid->soiltyp - 1];
//    grid->quartz = lsm->soiltbl.qtz[grid->soiltyp - 1];
//    grid->smcdry = lsm->soiltbl.drysmc[grid->soiltyp - 1];
//    grid->smcmax = lsm->soiltbl.maxsmc[grid->soiltyp - 1];
//    grid->smcref = lsm->soiltbl.refsmc[grid->soiltyp - 1];
//    grid->smcwlt = lsm->soiltbl.wltsmc[grid->soiltyp - 1];
//
///*----------------------------------------------------------------------
//* set-up universal parameters (not dependent on soiltyp, vegtyp or
//* slopetyp)
//* --------------------------------------------------------------------*/
//    grid->zbot = lsm->genprmt.zbot_data;
//    grid->salp = lsm->genprmt.salp_data;
//    grid->sbeta = lsm->genprmt.sbeta_data;
//    grid->frzk = lsm->genprmt.frzk_data;
//    grid->fxexp = lsm->genprmt.fxexp_data;
//    grid->ptu = 0.;             /* (not used yet) to satisify intent(out) */
//    grid->czil = lsm->genprmt.czil_data;
//    grid->lvcoef = lsm->genprmt.lvcoef_data;
//#ifndef _NOAH_
//    grid->slope = lsm->genprmt.slope_data[grid->slopetyp - 1];
//    grid->refkdt = lsm->genprmt.refkdt_data;
//    grid->refdk = lsm->genprmt.refdk_data;
//    grid->kdt = grid->refkdt * grid->dksat / grid->refdk;
//#endif
//
///*----------------------------------------------------------------------
//* to adjust frzk parameter to actual soil type: frzk * frzfact
//* --------------------------------------------------------------------*/
//    frzfact = (grid->smcmax / grid->smcref) * (0.412 / 0.468);
//    grid->frzx = grid->frzk * frzfact;
//
///*----------------------------------------------------------------------
//* set-up vegetation parameters
//* --------------------------------------------------------------------*/
//    grid->topt = lsm->vegtbl.topt_data;
//    grid->cfactr = lsm->vegtbl.cfactr_data;
//    grid->rsmax = lsm->vegtbl.rsmax_data;
//#ifdef _NOAH_
//    grid->cmcmax = lsm->vegtbl.cmcfactrtbl[grid->vegtyp - 1] * grid->xlai;
//#else
//    grid->cmcmax = lsm->vegtbl.cmcmax_data;
//#endif
//    grid->nroot = lsm->vegtbl.nrotbl[grid->vegtyp - 1];
//    grid->snup = lsm->vegtbl.snuptbl[grid->vegtyp - 1];
//    grid->rsmin = lsm->vegtbl.rstbl[grid->vegtyp - 1];
//    grid->rgl = lsm->vegtbl.rgltbl[grid->vegtyp - 1];
//    grid->hs = lsm->vegtbl.hstbl[grid->vegtyp - 1];
//    grid->emissmin = lsm->vegtbl.emissmintbl[grid->vegtyp - 1];
//    grid->emissmax = lsm->vegtbl.emissmaxtbl[grid->vegtyp - 1];
//    grid->laimin = lsm->vegtbl.laimintbl[grid->vegtyp - 1];
//    grid->laimax = lsm->vegtbl.laimaxtbl[grid->vegtyp - 1];
//    grid->z0min = lsm->vegtbl.z0mintbl[grid->vegtyp - 1];
//    grid->z0max = lsm->vegtbl.z0maxtbl[grid->vegtyp - 1];
//    grid->albedomin = lsm->vegtbl.albedomintbl[grid->vegtyp - 1];
//    grid->albedomax = lsm->vegtbl.albedomaxtbl[grid->vegtyp - 1];
//
//    grid->isurban = lsm->vegtbl.isurban;
//
//    if (grid->vegtyp == lsm->vegtbl.bare)
//        grid->shdfac = 0.0;
//
//    if (grid->nroot > grid->nsoil)
//    {
//        printf ("error: too many root layers %d, %d\n", grid->nsoil,
//           grid->nroot);
//        exit (0);
//    }
//
///*----------------------------------------------------------------------
//* calculate root distribution.  present version assumes uniform
//* distribution based on soil layer depths.
//* --------------------------------------------------------------------*/
//    for (i = 0; i < grid->nroot; i++)
//        grid->rtdis[i] = -grid->sldpth[i] / zsoil[grid->nroot - 1];
//
///*----------------------------------------------------------------------
//*  set-up slope parameter
//* --------------------------------------------------------------------*/
//
//    /*
//     * print*,'end of prmred'
//     * *       print*,'vegtyp',vegtyp,'soiltyp',soiltyp,'slopetyp',slopetyp,    &
//     * *    & 'cfactr',cfactr,'cmcmax',cmcmax,'rsmax',rsmax,'topt',topt,        &
//     * *    & 'refkdt',refkdt,'kdt',kdt,'sbeta',sbeta, 'shdfac',shdfac,         &
//     * *    &  'rsmin',rsmin,'rgl',rgl,'hs',hs,'zbot',zbot,'frzx',frzx,         &
//     * *    &  'psisat',psisat,'slope',slope,'snup',snup,'salp',salp,'bexp',    &
//     * *    &   bexp,                                                           &
//     * *    &  'dksat',dksat,'dwsat',dwsat,                                     &
//     * *    &  'smcmax',smcmax,'smcwlt',smcwlt,'smcref',smcref,'smcdry',smcdry, &
//     * *    &  'f1',f1,'quartz',quartz,'fxexp',fxexp,                           &
//     * *    &  'rtdis',rtdis,'sldpth',sldpth,'zsoil',zsoil, 'nroot',nroot,      &
//     * *    &  'nsoil',nsoil,'z0',z0,'czil',czil,'lai',lai,                     &
//     * *    &  'csoil',csoil,'ptu',ptu,                                         &
//     * *    &  'local', local
//     */
//
//}

void
Rosr12 (double *p, double *a, double *b, double *c, double *d, double *delta,
   int *nsoil)
{

/*----------------------------------------------------------------------
* subroutine Rosr12
* ----------------------------------------------------------------------
* invert (solve) the tri-diagonal matrix problem shown below:
* ###                                            ### ###  ###   ###  ###
* #b(1), c(1),  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
* #a(2), b(2), c(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #
* # 0  , a(3), b(3), c(3),  0  ,   . . .  ,    0   # #      #   # d(3) #
* # 0  ,  0  , a(4), b(4), c(4),   . . .  ,    0   # # p(4) #   # d(4) #
* # 0  ,  0  ,  0  , a(5), b(5),   . . .  ,    0   # # p(5) #   # d(5) #
* # .                                          .   # #  .   # = #   .  #
* # .                                          .   # #  .   #   #   .  #
* # .                                          .   # #  .   #   #   .  #
* # 0  , . . . , 0 , a(m-2), b(m-2), c(m-2),   0   # #p(m-2)#   #d(m-2)#
* # 0  , . . . , 0 ,   0   , a(m-1), b(m-1), c(m-1)# #p(m-1)#   #d(m-1)#
* # 0  , . . . , 0 ,   0   ,   0   ,  a(m) ,  b(m) # # p(m) #   # d(m) #
* ###                                            ### ###  ###   ###  ###
* --------------------------------------------------------------------*/

    int             k, kk;

/*----------------------------------------------------------------------
* initialize eqn coef c for the lowest soil layer
* --------------------------------------------------------------------*/
    c[*nsoil - 1] = 0.0;
    p[0] = -c[0] / b[0];

/*----------------------------------------------------------------------
* solve the coefs for the 1st soil layer
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* solve the coefs for soil layers 2 thru nsoil
* --------------------------------------------------------------------*/
    delta[0] = d[0] / b[0];
    for (k = 1; k < *nsoil; k++)
    {
        p[k] = -c[k] * (1.0 / (b[k] + a[k] * p[k - 1]));
        delta[k] =
           (d[k] - a[k] * delta[k - 1]) * (1.0 / (b[k] + a[k] * p[k - 1]));
    }

/*----------------------------------------------------------------------
* set p to delta for lowest soil layer
* --------------------------------------------------------------------*/
    p[*nsoil - 1] = delta[*nsoil - 1];

/*----------------------------------------------------------------------
* adjust p for soil layers 2 thru nsoil
* --------------------------------------------------------------------*/
    for (k = 1; k < *nsoil; k++)
    {
        kk = *nsoil - k - 1;
        p[kk] = p[kk] * p[kk + 1] + delta[kk];
    }

/*----------------------------------------------------------------------
  end subroutine Rosr12
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
ShFlx (double *ssoil, double *stc, double *smc, double *smcmax,
   double *smcmin, int *nsoil, double *t1, double *dt, double *yy,
   double *zz1, double *zsoil, double *tbot, double *zbot, double *smcwlt,
   double *sh2o, double *vgalpha, double *vgbeta, double *f1, double *df1,
   double *quartz, double *csoil, int *vegtyp, int *isurban)
#else
void
ShFlx (double *ssoil, double *stc, double *smc, double *smcmax, int *nsoil,
   double *t1, double *dt, double *yy, double *zz1, double *zsoil,
   double *tbot, double *zbot, double *smcwlt, double *psisat, double *sh2o,
   double *bexp, double *f1, double *df1, double *quartz, double *csoil,
   int *vegtyp, int *isurban)
#endif
{

/*----------------------------------------------------------------------
* subroutine ShFlx
* ----------------------------------------------------------------------
* update the temperature state of the soil column based on the thermal
* diffusion equation and update the frozen soil moisture content based
* on the temperature.
* --------------------------------------------------------------------*/
    int             i;

    double          ai[*nsoil], bi[*nsoil], ci[*nsoil], stcf[*nsoil],
       rhsts[*nsoil];

/*----------------------------------------------------------------------
* HRT routine calcs the right hand side of the soil temp dif eqn
* --------------------------------------------------------------------*/

    /*
     * land case 
     */

#ifdef _NOAH_
    HRT (rhsts, stc, smc, smcmax, smcmin, nsoil, zsoil, yy, zz1, tbot, zbot,
       sh2o, dt, vgalpha, vgbeta, f1, df1, quartz, csoil, ai, bi, ci, vegtyp,
       isurban);
#else
    HRT (rhsts, stc, smc, smcmax, nsoil, zsoil, yy, zz1, tbot, zbot, psisat,
       sh2o, dt, bexp, f1, df1, quartz, csoil, ai, bi, ci, vegtyp, isurban);
#endif

    HStep (stcf, stc, rhsts, dt, nsoil, ai, bi, ci);

    for (i = 0; i < *nsoil; i++)
        stc[i] = stcf[i];

/*----------------------------------------------------------------------
* in the no snowpack case (via routine NoPac branch,) update the grnd
* (skin) temperature here in response to the updated soil temperature
* profile above.  (note: inspection of routine SnoPac shows that t1
* below is a dummy variable only, as skin temperature is updated
* differently in routine SnoPac)
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* calculate surface soil heat flux
* --------------------------------------------------------------------*/
    *t1 = (*yy + (*zz1 - 1.0) * stc[0]) / *zz1;
    *ssoil = *df1 * (stc[0] - *t1) / (0.5 * zsoil[0]);

/*----------------------------------------------------------------------
  end subroutine ShFlx
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
SmFlx (double *smc, int *nsoil, double *cmc, double *dt, double *prcp1,
   double *pcpdrp, double *zsoil, double *sh2o, double *frzfact,
   double *smcmax, double *smcmin, double *vgalpha, double *vgbeta,
   double *macksat, double *areaf, int *nmacd, int *mac_status, int *nwtbl, double *smcwlt,
   double *dksat, double *shdfac, double *cmcmax, double *infil,
   double *runoff2, double *runoff3, double *edir, double *ec, double *et,
   double *drip)
#else
void
SmFlx (double *smc, int *nsoil, double *cmc, double *dt, double *prcp1,
   double *zsoil, double *sh2o, double *slope, double *kdt, double *frzfact,
   double *smcmax, double *bexp, double *smcwlt, double *dksat, double *dwsat,
   double *shdfac, double *cmcmax, double *runoff1, double *runoff2,
   double *runoff3, double *edir, double *ec, double *et, double *drip)
#endif
{

/*----------------------------------------------------------------------
* subroutine SmFlx
* ----------------------------------------------------------------------
* calculate soil moisture flux.  the soil moisture content (smc - a per
* unit volume measurement) is a dependent variable that is updated with
* prognostic eqns. the canopy moisture content (cmc) is also updated.
* frozen ground version:  new states added: sh2o, and frozen ground
* correction factor, frzfact and parameter slope.
* --------------------------------------------------------------------*/

    int             i;
    double          ai[*nsoil], bi[*nsoil], ci[*nsoil];
    double          rhstt[*nsoil];
    double          sice[*nsoil];
    double         *dummy;
    double          excess;
    double         *rhsct;
    double          trhsct;
#ifndef _NOAH_
    double         *pcpdrp;
#endif
    //double          fac2;
    double         *flimit;
#ifdef _NOAH_
    double          kd = 6.54e-7, bfactr = 3.89;
#endif
    dummy = (double *)malloc (sizeof (double));
#ifndef _NOAH_
    pcpdrp = (double *)malloc (sizeof (double));
#endif
    rhsct = (double *)malloc (sizeof (double));
    flimit = (double *)malloc (sizeof (double));

/*----------------------------------------------------------------------
* executable code begins here.
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* compute the right hand side of the canopy eqn term ( rhsct )
* --------------------------------------------------------------------*/
    *dummy = 0.;

/*----------------------------------------------------------------------
* convert rhsct (a rate) to trhsct (an amount) and add it to existing
* cmc.  if resulting amt exceeds max capacity, it becomes drip and will
* fall to the grnd.
* --------------------------------------------------------------------*/
    *rhsct = *shdfac * *prcp1 - *ec;
    *drip = 0.;
    trhsct = *dt * *rhsct;
    excess = *cmc + trhsct;

/*----------------------------------------------------------------------
* pcpdrp is the combined prcp1 and drip (from cmc) that goes into the
* soil
* --------------------------------------------------------------------*/
#ifdef _NOAH_

/*----------------------------------------------------------------------
* pihm drip calculation following rutter and mortan (1977 jae)
* --------------------------------------------------------------------*/
    if (excess > 0)
    {
        if (excess >= *cmcmax)
        {
            *drip = (kd * *cmcmax * exp (bfactr)) * *dt + excess - *cmcmax;
            *rhsct = *rhsct - kd * *cmcmax * exp (bfactr);
        }
        else
        {
            *drip = (kd * *cmcmax * exp (bfactr * excess / *cmcmax)) * *dt;
            *rhsct = *rhsct - kd * *cmcmax * exp (bfactr * excess / *cmcmax);
        }
    }
#else
    if (excess > *cmcmax)
        *drip = excess - *cmcmax;
#endif
    *pcpdrp = (1. - *shdfac) * *prcp1 + *drip / *dt;

/*----------------------------------------------------------------------
* store ice content at each soil layer before calling SRT and SStep
* --------------------------------------------------------------------*/
    for (i = 0; i < *nsoil; i++)
        sice[i] = smc[i] - sh2o[i];

/*----------------------------------------------------------------------
* call subroutines SRT and SStep to solve the soil moisture
* tendency equations.
* if the infiltrating precip rate is nontrivial,
*   (we consider nontrivial to be a precip total over the time step
*    exceeding one one-thousandth of the water holding capacity of
*    the first soil layer)
* then call the SRT/SStep subroutine pair twice in the manner of
*   time scheme "f" (implicit state, averaged coefficient)
*   of section 2 of kalnay and kanamitsu (1988, mwr, vol 116,
*   pages 1945-1958)to minimize 2-delta-t oscillations in the
*   soil moisture value of the top soil layer that can arise because
*   of the extreme nonlinear dependence of the soil hydraulic
*   diffusivity coefficient and the hydraulic conductivity on the
*   soil moisture state
* otherwise call the SRT/SStep subroutine pair once in the manner of
*   time scheme "d" (implicit state, explicit coefficient)
*   of section 2 of kalnay and kanamitsu
* pcpdrp is units of kg/m**2/s or mm/s, zsoil is negative depth in m
* ----------------------------------------------------------------------
*  according to dr. ken mitchell's suggestion, add the second contraint
*  to remove numerical instability of runoff and soil moisture
*  flimit is a limit value for fac2
* --------------------------------------------------------------------*/
#ifdef _NOAH_

/*----------------------------------------------------------------------
* frozen ground version:
* smc states replaced by sh2o states in SRT subr.  sh2o & sice states
* inc&uded in SStep subr.  frozen ground correction factor, frzfact
* added.  all water balance calculations using unfrozen water
* --------------------------------------------------------------------*/
    if (*nwtbl == 0)
    {
        for (i = 0; i < *nsoil; i++)
        {
            smc[i] = *smcmax;
            sh2o[i] = smc[i] - sice[i];
        }
    }
    else
    {

        //      if ((*pcpdrp * *dt) > (0.0001 * 1000.0 * (-zsoil[0]) * *smcmax))
        //      {
        //          SRT(rhstt, edir, et, sh2o, sh2o, nwtbl, pcpdrp, zsoil, dksat, smcmax, smcmin, vgalpha, vgbeta, macksat, areaf, nmacd, infil, runoff2, dt, smcwlt, frzfact, sice, ai, bi, ci);
        //          SStep(sh2ofg, sh2o, dummy, rhstt, rhsct, dt, nwtbl, smcmax, smcmin, cmcmax, runoff3, zsoil, smc, sice, ai, bi, ci);
        //          for (k = 0; k < *nsoil; k++)
        //              sh2oa[k] = (sh2o[k] + sh2ofg[k]) * 0.5;
        //          SRT (rhstt, edir, et, sh2o, sh2oa, nwtbl, pcpdrp, zsoil, dksat, smcmax, smcmin, vgalpha, vgbeta, macksat, areaf, nmacd, infil, runoff2, dt, smcwlt, frzfact, sice, ai, bi, ci);
        //          SStep (sh2o, sh2o, cmc, rhstt, rhsct, dt, nwtbl, smcmax, cmcmax, smcmin, runoff3, zsoil, smc, sice, ai, bi, ci);
        //      }
        //      else
        //      {
        //          SRT(rhstt, edir, et, sh2o, sh2o, nwtbl, pcpdrp, zsoil, dksat, smcmax, smcmin, vgalpha, vgbeta, macksat, areaf, nmacd, infil, runoff2, dt, smcwlt, frzfact, sice, ai, bi, ci);
        //          SStep (sh2o, sh2o, cmc, rhstt, rhsct, dt, nwtbl, smcmax, smcmin, cmcmax, runoff3, zsoil, smc, sice, ai, bi, ci);
        //      }
        SRT (rhstt, edir, et, sh2o, sh2o, nsoil, nwtbl, pcpdrp, zsoil, dksat, smcmax,
           smcmin, vgalpha, vgbeta, macksat, areaf, nmacd, mac_status, infil, runoff2, dt,
           smcwlt, frzfact, sice, ai, bi, ci);
        SStep (sh2o, sh2o, cmc, rhstt, rhsct, dt, nsoil, smcmax, smcmin,
           cmcmax, runoff3, zsoil, smc, sice, ai, bi, ci);
    }

#else
    fac2 = 0.0;
    for (i = 0; i < *nsoil; i++)
    {
        fac2 = fac2 > (sh2o[i] / *smcmax) ? fac2 : (sh2o[i] / *smcmax);
    }
    Fac2Mit (smcmax, flimit);

/*----------------------------------------------------------------------
* frozen ground version:
* smc states replaced by sh2o states in SRT subr.  sh2o & sice states
* inc&uded in SStep subr.  frozen ground correction factor, frzfact
* added.  all water balance calculations using unfrozen water
* --------------------------------------------------------------------*/

    if (((*pcpdrp * *dt) > (0.0001 * 1000.0 * (-zsoil[0]) * *smcmax))
       || (fac2 > *flimit))
    {
        SRT (rhstt, edir, et, sh2o, sh2o, nsoil, pcpdrp, zsoil, dwsat, dksat,
           smcmax, bexp, runoff1, runoff2, dt, smcwlt, slope, kdt, frzfact,
           sice, ai, bi, ci);
        SStep (sh2ofg, sh2o, dummy, rhstt, rhsct, dt, nsoil, smcmax, cmcmax,
           runoff3, zsoil, smc, sice, ai, bi, ci);
        for (k = 0; k < *nsoil; k++)
            sh2oa[k] = (sh2o[k] + sh2ofg[k]) * 0.5;
        SRT (rhstt, edir, et, sh2o, sh2oa, nsoil, pcpdrp, zsoil, dwsat, dksat,
           smcmax, bexp, runoff1, runoff2, dt, smcwlt, slope, kdt, frzfact,
           sice, ai, bi, ci);
        SStep (sh2o, sh2o, cmc, rhstt, rhsct, dt, nsoil, smcmax, cmcmax,
           runoff3, zsoil, smc, sice, ai, bi, ci);
    }
    else
    {
        SRT (rhstt, edir, et, sh2o, sh2o, nsoil, pcpdrp, zsoil, dwsat, dksat,
           smcmax, bexp, runoff1, runoff2, dt, smcwlt, slope, kdt, frzfact,
           sice, ai, bi, ci);
        SStep (sh2o, sh2o, cmc, rhstt, rhsct, dt, nsoil, smcmax, cmcmax,
           runoff3, zsoil, smc, sice, ai, bi, ci);
        //      runof = runoff

    }
#endif
    free (dummy);
#ifndef _NOAH_
    free (pcpdrp);
#endif
    free (rhsct);
    free (flimit);

/*----------------------------------------------------------------------
  end subroutine SmFlx
* --------------------------------------------------------------------*/
}

void
SnFrac (double *sneqv, double *snup, double *salp, double *snowh,
   double *sncovr)
{

/*----------------------------------------------------------------------
* subroutine SnFrac
* ----------------------------------------------------------------------
* calculate snow fraction (0 -> 1)
* sneqv   snow water equivalent (m)
* snup    threshold sneqv depth above which sncovr=1
* salp    tuning parameter
* sncovr  fractional snow cover
* --------------------------------------------------------------------*/

    double          rsnow;

/*----------------------------------------------------------------------
* snup is veg-class dependent snowdepth threshhold (set in routine
* RedPrm) above which snocvr=1.
* --------------------------------------------------------------------*/
    if (*sneqv < *snup)
    {
        rsnow = *sneqv / *snup;
        *sncovr = 1. - (exp (-*salp * rsnow) - rsnow * exp (-*salp));
    }
    else
        *sncovr = 1.0;

    /*
     * formulation of dickinson et al. 1986 
     */
    //  z0n = 0.035;

    //  *sncovr = *snowh / (*snowh + 5. * z0n);

    /*
     * formulation of marshall et al. 1994 
     */
    //  *sncovr = *sneqv / (*sneqv + 2. * z0n);

/*----------------------------------------------------------------------
  end subroutine SnFrac
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
SnkSrc (double *tsnsr, double *tavg, double *smc, double *sh2o, double *zsoil,
   int *nsoil, double *smcmax, double *smcmin, double *vgalpha,
   double *vgbeta, double *dt, int k, double *qtot)
#else
void
SnkSrc (double *tsnsr, double *tavg, double *smc, double *sh2o, double *zsoil,
   int *nsoil, double *smcmax, double *psisat, double *bexp, double *dt,
   int k, double *qtot)
#endif
{

/*----------------------------------------------------------------------
* subroutine SnkSrc
* ----------------------------------------------------------------------
* calculate sink/source term of the thermal diffusion equation. (sh2o)
* is available liqued water.
* --------------------------------------------------------------------*/

    double          dz, *freew, xh2o;

    double          dh2o = 1.0000e3, hlice = 3.3350e5;

    freew = (double *)malloc (sizeof (double));

    if (k == 0)
        dz = -zsoil[0];
    else
        dz = zsoil[k - 1] - zsoil[k];

/*----------------------------------------------------------------------
* via function FrH2O, compute potential or 'equilibrium' unfrozen
* supercooled free water for given soil type and soil layer temperature.
* function frh20 invokes eqn (17) from v. koren et al (1999, jgr, vol.
* 104, pg 19573).  (aside:  latter eqn in journal in centigrade units.
* routine FrH2O use form of eqn in kelvin units.)
* --------------------------------------------------------------------*/

    /*
     * freew = FrH2O(tavg,smc,sh2o,smcmax,bexp,psisat)   
     */

/*----------------------------------------------------------------------
* in next block of code, invoke eqn 18 of v. koren et al (1999, jgr,
* vol. 104, pg 19573.)  that is, first estimate the new amountof liquid
* water, 'xh2o', implied by the sum of (1) the liquid water at the begin
* of current time step, and (2) the freeze of thaw change in liquid
* water implied by the heat flux 'qtot' passed in from routine HRT.
* second, determine if xh2o needs to be bounded by 'free' (equil amt) or
* if 'free' needs to be bounded by xh2o.
* --------------------------------------------------------------------*/
#ifdef _NOAH_
    FrH2O (freew, tavg, smc, sh2o, smcmax, smcmin, vgalpha, vgbeta);
#else
    FrH2O (freew, tavg, smc, sh2o, smcmax, bexp, psisat);
#endif

/*----------------------------------------------------------------------
* first, if freezing and remaining liquid less than lower bound, then
* reduce extent of freezing, thereby letting some or all of heat flux
* qtot cool the soil temp later in routine HRT.
* --------------------------------------------------------------------*/
    xh2o = *sh2o + *qtot * *dt / (dh2o * hlice * dz);

    if (xh2o < *sh2o && xh2o < *freew)
    {
        if (*freew > *sh2o)
            xh2o = *sh2o;
        else
            xh2o = *freew;
    }

/*----------------------------------------------------------------------
* second, if thawing and the increase in liquid water greater than upper
* bound, then reduce extent of thaw, thereby letting some or all of heat
* flux qtot warm the soil temp later in routine HRT.
* --------------------------------------------------------------------*/

    if (xh2o > *sh2o && xh2o > *freew)
    {
        if (*freew < *sh2o)
            xh2o = *sh2o;
        else
            xh2o = *freew;
    }

/*----------------------------------------------------------------------
* calculate phase-change heat source/sink term for use in routine HRT
* and update liquid water to reflcet final freeze/thaw increment.
* --------------------------------------------------------------------*/
    //      SnkSrc = -dh2o*hlice*dz*(xh2o-sh2o)/ *dt
    if (xh2o < 0.)
        xh2o = 0.;
    if (xh2o > *smc)
        xh2o = *smc;
    *tsnsr = -dh2o * hlice * dz * (xh2o - *sh2o) / *dt;

    *sh2o = xh2o;

    free (freew);

/*----------------------------------------------------------------------
  end subroutine SnkSrc
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
SnoPac (double *etp, double *eta, double *prcp, double *prcpf, double *pcpdrp,
   int *snowng, double *smc, double *smcmax, double *smcmin, double *smcwlt,
   double *smcref, double *smcdry, double *cmc, double *cmcmax, int *nsoil,
   double *dt, double *sbeta, double *df1, double *q2, double *t1,
   double *sfctmp, double *t24, double *th2, double *fdown, double *f1,
   double *ssoil, double *stc, double *epsca, double *sfcprs, double *vgalpha,
   double *vgbeta, double *macksat, double *areaf, int *nmacd, int *mac_status, int *nwtbl,
   double *pc, double *rch, double *rr, double *cfactr, double *sncovr,
   double *esd, double *sndens, double *snowh, double *sh2o, double *frzfact,
   double *zsoil, double *dksat, double *tbot, double *zbot, double *shdfac,
   double *infil, double *runoff2, double *runoff3, double *edir, double *ec,
   double *et, double *ett, int *nroot, double *snomlt, double *rtdis,
   double *quartz, double *fxexp, double *csoil, double *beta, double *drip,
   double *dew, double *flx1, double *flx2, double *flx3, double *esnow,
   double *etns, double *emissi, double *ribb, double *soldn, int *isurban,
   int *vegtyp)
#else
void
SnoPac (double *etp, double *eta, double *prcp, double *prcpf, int *snowng,
   double *smc, double *smcmax, double *smcwlt, double *smcref,
   double *smcdry, double *cmc, double *cmcmax, int *nsoil, double *dt,
   double *sbeta, double *df1, double *q2, double *t1, double *sfctmp,
   double *t24, double *th2, double *fdown, double *f1, double *ssoil,
   double *stc, double *epsca, double *sfcprs, double *bexp, double *pc,
   double *rch, double *rr, double *cfactr, double *sncovr, double *esd,
   double *sndens, double *snowh, double *sh2o, double *slope, double *kdt,
   double *frzfact, double *psisat, double *zsoil, double *dwsat,
   double *dksat, double *tbot, double *zbot, double *shdfac, double *runoff1,
   double *runoff2, double *runoff3, double *edir, double *ec, double *et,
   double *ett, int *nroot, double *snomlt, double *rtdis, double *quartz,
   double *fxexp, double *csoil, double *beta, double *drip, double *dew,
   double *flx1, double *flx2, double *flx3, double *esnow, double *etns,
   double *emissi, double *ribb, double *soldn, int *isurban, int *vegtyp)
#endif
{

/*----------------------------------------------------------------------
* subroutine SnoPac
* ----------------------------------------------------------------------
* calculate soil moisture and heat flux values & update soil moisture
* content and soil heat content values for the case when a snow pack is
* present.
* --------------------------------------------------------------------*/

    int             k;

    double          et1[*nsoil];
    double          denom;
    double          dsoil;
    double          dtot;
    double          esnow1;
    double          esnow2;
    //double          etp2;
    double          etp3;
    double          etanrg;
    double          ex;
    double          seh;
    double          sncond;
    double          t12;
    double          t12a;
    double          t12b;
    double          t14;
    double         *ec1;
    double         *edir1;
    double         *ett1;
    double         *etp1;
    double         *etns1;
    double         *prcp1;
    double         *ssoil1;
    double         *t11;
    double         *yy;
    double         *zz1;
    double          esdmin = 1.0e-6;
    double          lsubc = 2.501e6;
    double          snoexp = 2.0;

    ec1 = (double *)malloc (sizeof (double));
    edir1 = (double *)malloc (sizeof (double));
    ett1 = (double *)malloc (sizeof (double));
    etp1 = (double *)malloc (sizeof (double));
    etns1 = (double *)malloc (sizeof (double));
    prcp1 = (double *)malloc (sizeof (double));
    ssoil1 = (double *)malloc (sizeof (double));
    t11 = (double *)malloc (sizeof (double));
    yy = (double *)malloc (sizeof (double));
    zz1 = (double *)malloc (sizeof (double));

/*----------------------------------------------------------------------
* executable code begins here:
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* initialize evap terms.
* ----------------------------------------------------------------------
* conversions:
* esnow [kg m-2 s-1]
* esdflx [kg m-2 s-1] .le. esnow
* esnow1 [m s-1]
* esnow2 [m]
* etp [kg m-2 s-1]
* etp1 [m s-1]
* etp2 [m]
* --------------------------------------------------------------------*/
    *dew = 0.;
    *edir = 0.;
    *edir1 = 0.;
    *ec1 = 0.;
    *ec = 0.;
    //      EMISSI_S=0.95    ! for snow

    for (k = 0; k < *nsoil; k++)
    {
        et[k] = 0.;
        et1[k] = 0.;
    }
    *ett = 0.;
    *ett1 = 0.;
    *etns = 0.;
    *etns1 = 0.;
    *esnow = 0.;
    esnow1 = 0.;
    esnow2 = 0.;

/*----------------------------------------------------------------------
* convert potential evap (etp) from kg m-2 s-1 to etp1 in m s-1
* --------------------------------------------------------------------*/
    *prcp1 = *prcpf * 0.001;

/*----------------------------------------------------------------------
* if etp<0 (downward) then dewfall (=frostfall in this case).
* --------------------------------------------------------------------*/
    *beta = 1.0;
    if (*etp <= 0.0)
    {
        if ((*ribb >= 0.1) && (*fdown > 150.0))
        {
            *etp =
               ((*etp * (1.0 - *ribb) <
                  0. ? *etp * (1.0 - *ribb) : 0.0) * *sncovr / 0.980 +
               *etp * (0.980 - *sncovr)) / 0.980;
        }
        if (*etp == 0.)
            *beta = 0.0;
        *etp1 = *etp * 0.001;
        *dew = -*etp1;
        esnow2 = *etp1 * *dt;
        etanrg = *etp * ((1. - *sncovr) * lsubc + *sncovr * LSUBS);
    }
    else
    {
        *etp1 = *etp * 0.001;
        /*
         * land case 
         */
        if (*sncovr < 1.)
        {
#ifdef _NOAH_
            Evapo (etns1, smc, nsoil, cmc, etp1, dt, zsoil, sh2o, smcmax, pc,
               smcwlt, dksat, smcref, shdfac, cmcmax, smcdry, cfactr, edir1,
               ec1, et1, ett1, sfctmp, q2, nroot, rtdis, fxexp);
#else
            Evapo (etns1, smc, nsoil, cmc, etp1, dt, zsoil, sh2o, smcmax,
               bexp, pc, smcwlt, dksat, dwsat, smcref, shdfac, cmcmax, smcdry,
               cfactr, edir1, ec1, et1, ett1, sfctmp, q2, nroot, rtdis,
               fxexp);
#endif

            /*
             * --------------------------------------------------------------------------
             */
            *edir1 = *edir1 * (1. - *sncovr);
            *ec1 = *ec1 * (1. - *sncovr);
            for (k = 0; k < *nsoil; k++)
                et1[k] = et1[k] * (1. - *sncovr);
            *ett1 = *ett1 * (1. - *sncovr);
            //          *etns1 = *edir1+ *ec1+ *ett1;
            *etns1 = *etns1 * (1. - *sncovr);

/*--------------------------------------------------------------------------*/
            *edir = *edir1 * 1000.;
            *ec = *ec1 * 1000.;
            for (k = 0; k < *nsoil; k++)
                et[k] = et1[k] * 1000.;
            *ett = *ett1 * 1000.;
            *etns = *etns1 * 1000.;

/*--------------------------------------------------------------------*/

        }
        *esnow = *etp * *sncovr;
        esnow1 = *esnow * 0.001;
        esnow2 = esnow1 * *dt;
        etanrg = *esnow * LSUBS + *etns * lsubc;
    }

/*----------------------------------------------------------------------
* if precip is falling, calculate heat flux from snow sfc to newly
* accumulating precip.  note that this reflects the flux appropriate for
* the not-yet-updated skin temperature (t1).  assumes temperature of the
* snowfall striking the ground is =sfctmp (lowest model level air temp).
* --------------------------------------------------------------------*/
    *flx1 = 0.0;
    if (*snowng)
        *flx1 = CPICE * *prcp * (*t1 - *sfctmp);
    else
    {
        if (*prcp > 0.0)
            *flx1 = CPH2O * *prcp * (*t1 - *sfctmp);

/*----------------------------------------------------------------------
* calculate an 'effective snow-grnd sfc temp' (t12) based on heat fluxes
* between the snow pack and the soil and on net radiation.
* include flx1 (precip-snow sfc) and flx2 (freezing rain latent heat)
* fluxes.  flx1 from above, flx2 brought in via commom block rite.
* flx2 reflects freezing rain latent heat flux using t1 calculated in
* Penman.
* --------------------------------------------------------------------*/
    }
    dsoil = -(0.5 * zsoil[0]);
    dtot = *snowh + dsoil;
    denom = 1.0 + *df1 / (dtot * *rr * *rch);

    /*
     * surface emissivity weighted by snow cover fraction
     * *      t12a = ( (*fdown - *flx1 - *flx2 -                                   &
     * *     &       ((*sncovr*EMISSI_S)+*emissi*(1.0-*sncovr))*SIGMA **t24)/ *rch    &
     * *     &       + *th2 - *sfctmp - etanrg / *rch ) / *rr
     */
    t12a =
       ((*fdown - *flx1 - *flx2 - *emissi * SIGMA * *t24) / *rch + *th2 -
       *sfctmp - etanrg / *rch) / *rr;
    t12b = *df1 * stc[0] / (dtot * *rr * *rch);

/*----------------------------------------------------------------------
* if the 'effective snow-grnd sfc temp' is at or below freezing, no snow
* melt will occur.  set the skin temp to this effective temp.  reduce
* (by sublimination ) or increase (by frost) the depth of the snowpack,
* depending on sign of etp.
* update soil heat flux (ssoil) using new skin temperature (t1)
* since no snowmelt, set accumulated snowmelt to zero, set 'effective'
* precip from snowmelt to zero, set phase-change heat flux from snowmelt
* to zero.
* ----------------------------------------------------------------------
* sub-freezing block
* --------------------------------------------------------------------*/
    t12 = (*sfctmp + t12a + t12b) / denom;

    if (t12 <= TFREEZ)
    {
        *t1 = t12;
        *ssoil = *df1 * (*t1 - stc[0]) / dtot;

        //      *esd = max (0.0, *esd- etp2)
        *esd = *esd - esnow2 > 0. ? (*esd - esnow2) : 0.;
        *flx3 = 0.0;
        ex = 0.0;

        *snomlt = 0.0;
    }

/*----------------------------------------------------------------------
* if the 'effective snow-grnd sfc temp' is above freezing, snow melt
* will occur.  call the snow melt rate,ex and amt, snomlt.  revise the
* effective snow depth.  revise the skin temp because it would have chgd
* due to the latent heat released by the melting. calc the latent heat
* released, flx3. set the effective precip, prcp1 to the snow melt rate,
* ex for use in SmFlx.  adjustment to t1 to account for snow patches.
* calculate qsat valid at freezing point.  note that esat (saturation
* vapor pressure) value of 6.11e+2 used here is that valid at frzzing
* point.  note that etp from call Penman in sflx is ignored here in
* favor of bulk etp over 'open water' at freezing temp.
* update soil heat flux (s) using new skin temperature (t1)
* ----------------------------------------------------------------------
* above freezing block
* --------------------------------------------------------------------*/
    else
    {
        *t1 =
           TFREEZ * pow (*sncovr, snoexp) + t12 * (1.0 - pow (*sncovr,
              snoexp));
        *beta = 1.0;

/*----------------------------------------------------------------------
* if potential evap (sublimation) greater than depth of snowpack.
* beta<1
* snowpack has sublimated away, set depth to zero.
* --------------------------------------------------------------------*/
        *ssoil = *df1 * (*t1 - stc[0]) / dtot;

        if (*esd - esnow2 <= esdmin)
        {
            *esd = 0.0;
            ex = 0.0;
            *snomlt = 0.0;
            *flx3 = 0.0;

/*----------------------------------------------------------------------
* sublimation less than depth of snowpack
* snowpack (esd) reduced by esnow2 (depth of sublimated snow)
* --------------------------------------------------------------------*/
        }
        else
        {
            *esd = *esd - esnow2;
            etp3 = *etp * lsubc;
            seh = *rch * (*t1 - *th2);
            t14 = *t1 * *t1;
            t14 = t14 * t14;
            //          *flx3 = *fdown - *flx1 - *flx2 - ((*sncovr*EMISSI_S)+*emissi*(1-*sncovr))*SIGMA*t14 - *ssoil - seh - etanrg;
            *flx3 =
               *fdown - *flx1 - *flx2 - *emissi * SIGMA * t14 - *ssoil - seh -
               etanrg;
            if (*flx3 <= 0.0)
                *flx3 = 0.0;

/*----------------------------------------------------------------------
* snowmelt reduction depending on snow cover
* --------------------------------------------------------------------*/
            ex = *flx3 * 0.001 / LSUBF;

/*----------------------------------------------------------------------
* esdmin represents a snowpack depth threshold value below which we
* choose not to retain any snowpack, and instead include it in snowmelt.
* --------------------------------------------------------------------*/
            *snomlt = ex * *dt;
            if (*esd - *snomlt >= esdmin)
            {
                *esd = *esd - *snomlt;

/*----------------------------------------------------------------------
* snowmelt exceeds snow depth
* --------------------------------------------------------------------*/
            }
            else
            {
                ex = *esd / *dt;
                *flx3 = ex * 1000.0 * LSUBF;
                *snomlt = *esd;

                *esd = 0.0;

/*----------------------------------------------------------------------
* end of 'esd .le. etp2' if-block
----------------------------------------------------------------------*/
            }
        }

/*----------------------------------------------------------------------
* end of 't12 .le. TFREEZ' if-block
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* if non-glacial land, add snowmelt rate (ex) to precip rate to be used
* in subroutine SmFlx (soil moisture evolution) via infiltration.
*
* runoff/baseflow later near the end of sflx (after return from call to
* subroutine SnoPac)
* --------------------------------------------------------------------*/
        *prcp1 = *prcp1 + ex;

/*----------------------------------------------------------------------
* set the effective potnl evapotransp (etp1) to zero since this is snow
* case, so surface evap not calculated from edir, ec, or ett in SmFlx
* (below).
* SmFlx returns updated soil moisture values for non-glacial land.
* --------------------------------------------------------------------*/
    }
#ifdef _NOAH_
    SmFlx (smc, nsoil, cmc, dt, prcp1, pcpdrp, zsoil, sh2o, frzfact, smcmax,
       smcmin, vgalpha, vgbeta, macksat, areaf, nmacd, mac_status, nwtbl, smcwlt, dksat,
       shdfac, cmcmax, infil, runoff2, runoff3, edir1, ec1, et1, drip);
#else
    SmFlx (smc, nsoil, cmc, dt, prcp1, zsoil, sh2o, slope, kdt, frzfact,
       smcmax, bexp, smcwlt, dksat, dwsat, shdfac, cmcmax, runoff1, runoff2,
       runoff3, edir1, ec1, et1, drip);
#endif

/*----------------------------------------------------------------------
* before call ShFlx in this snowpack case, set zz1 and yy arguments to
* special values that ensure that ground heat flux calculated in ShFlx
* matches that already computer for below the snowpack, thus the sfc
* heat flux to be computed in ShFlx will effectively be the flux at the
* snow top surface.  t11 is a dummy arguement so we will not use the
* skin temp value as revised by ShFlx.
* --------------------------------------------------------------------*/
    *zz1 = 1.0;
    *yy = stc[0] - 0.5 * *ssoil * zsoil[0] * *zz1 / *df1;

/*----------------------------------------------------------------------
* ShFlx will calc/update the soil temps.  note:  the sub-sfc heat flux
* (ssoil1) and the skin temp (t11) output from this ShFlx call are not
* used  in any subsequent calculations. rather, they are dummy variables
* here in the SnoPac case, since the skin temp and sub-sfc heat flux are
* updated instead near the beginning of the call to SnoPac.
* --------------------------------------------------------------------*/
    *t11 = *t1;
#ifdef _NOAH_
    ShFlx (ssoil1, stc, smc, smcmax, smcmin, nsoil, t11, dt, yy, zz1, zsoil,
       tbot, zbot, smcwlt, sh2o, vgalpha, vgbeta, f1, df1, quartz, csoil,
       vegtyp, isurban);
#else
    ShFlx (ssoil1, stc, smc, smcmax, nsoil, t11, dt, yy, zz1, zsoil, tbot,
       zbot, smcwlt, psisat, sh2o, bexp, f1, df1, quartz, csoil, vegtyp,
       isurban);
#endif

/*----------------------------------------------------------------------
* snow depth and density adjustment based on snow compaction.  yy is
* assumed to be the soil temperture at the top of the soil column.
* --------------------------------------------------------------------*/
    /*
     * land 
     */
    if (*esd > 0.)
        SnowPack (esd, dt, snowh, sndens, t1, yy);
    else
    {
        *esd = 0.;
        *snowh = 0.;
        *sndens = 0.;
        sncond = 1.;
        *sncovr = 0.;
    }

    free (ec1);
    free (edir1);
    free (ett1);
    free (etp1);
    free (etns1);
    free (prcp1);
    free (ssoil1);
    free (t11);
    free (yy);
    free (zz1);

/*----------------------------------------------------------------------
  end subroutine SnoPac
* --------------------------------------------------------------------*/
}

void SnowPack (double *esd, double *dtsec, double *snowh, double *sndens, double *tsnow, double *tsoil)
{
/*----------------------------------------------------------------------
* subroutine SnowPack
* ----------------------------------------------------------------------
* calculate compaction of snowpack under conditions of increasing snow
* density, as obtained from an approximate solution of e. anderson's
* differential equation (3.29), noaa technical report nws 19, by victor
* koren, 03/25/95.
* ----------------------------------------------------------------------
* esd     water equivalent of snow (m)
* dtsec   time step (sec)
* snowh   snow depth (m)
* sndens  snow density (g/cm3=dimensionless fraction of h2o density)
* tsnow   snow surface temperature (k)
* tsoil   soil surface temperature (k)

* subroutine will return new values of snowh and sndens
* --------------------------------------------------------------------*/

    int             ipol;
    int             j;
    double          bfac;
    double          dsx;
    double          dthr;
    double          dw;
    double          snowhc;
    double          pexp;
    double          tavgc;
    double          tsnowc;
    double          tsoilc;
    double          esdc;
    double          esdcx;
    double          c1 = 0.01;
    double          c2 = 21.0;

/*----------------------------------------------------------------------
* conversion into simulation units
* --------------------------------------------------------------------*/
    snowhc = *snowh * 100.0;
    esdc = *esd * 100.0;
    dthr = *dtsec / 3600.0;
    tsnowc = *tsnow - 273.15;
    tsoilc = *tsoil - 273.15;

/*----------------------------------------------------------------------
* calculating of average temperature of snow pack
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* calculating of snow depth and density as a result of compaction
*  sndens=ds0*(exp(bfac*esd)-1.)/(bfac*esd)
*  bfac=dthr*c1*exp(0.08*tavgc-c2*ds0)
* note: bfac*esd in sndens eqn above has to be carefully treated
* numerically below:
*   c1 is the fractional increase in density (1/(cm*hr))
*   c2 is a constant (cm3/g) kojima estimated as 21 cms/g
* --------------------------------------------------------------------*/
    tavgc = 0.5 * (tsnowc + tsoilc);
    if (esdc > 1.e-2)
        esdcx = esdc;
    else
        esdcx = 1.e-2;

    //      dsx = *sndens*((dexp(bfac*esdc)-1.)/(bfac*esdc))

/*----------------------------------------------------------------------
* the function of the form (e**x-1)/x embedded in above expression
* for dsx was causing numerical difficulties when the denominator "x"
* (i.e. bfac*esdc) became zero or approached zero (despite the fact that
* the analytical function (e**x-1)/x has a well defined limit as
* "x" approaches zero), hence below we replace the (e**x-1)/x
* expression with an equivalent, numerically well-behaved
* polynomial expansion.

* number of terms of polynomial expansion, and hence its accuracy,
* is governed by iteration limit "ipol".
*      ipol greater than 9 only makes a difference on double
*            precision (relative errors given in percent %).
*       ipol=9, for rel.error <~ 1.6 e-6 % (8 significant digits)
*       ipol=8, for rel.error <~ 1.8 e-5 % (7 significant digits)
*       ipol=7, for rel.error <~ 1.8 e-4 % ...
* --------------------------------------------------------------------*/
    bfac = dthr * c1 * exp (0.08 * tavgc - c2 * *sndens);
    ipol = 4;
    pexp = 0.;
    //  pexp = (1. + pexp)*bfac*esdc/real(j+1);
    for (j = ipol; j > 0; j--)
    {
        pexp = (1. + pexp) * bfac * esdcx / (double)(j + 1);
    }

    pexp = pexp + 1.;

/*----------------------------------------------------------------------
* above line ends polynomial substitution
* ----------------------------------------------------------------------
*     end of korean formulation
* --------------------------------------------------------------------*/

    /*
     *     base formulation (cogley et al., 1990)
     *     convert density from g/cm3 to kg/m3
     *       dsm=*sndens*1000.0
     
     *       dsx=dsm+*dtsec*0.5*dsm*g* *esd/
     *    &      (1e7*exp(-0.02*dsm+kn/(tavgc+273.16)-14.643))
     
     *  &   convert density from kg/m3 to g/cm3
     *       dsx=dsx/1000.0
     
     *     end of cogley et al. formulation
     */

/*----------------------------------------------------------------------
* set upper/lower limit on snow density
* --------------------------------------------------------------------*/
    dsx = *sndens * (pexp);
    if (dsx > 0.40)
        dsx = 0.40;
    if (dsx < 0.05)
        dsx = 0.05;

/*----------------------------------------------------------------------
* update of snow depth and density depending on liquid water during
* snowmelt.  assumed that 13% of liquid water can be stored in snow per
* day during snowmelt till snow density 0.40.
* --------------------------------------------------------------------*/
    *sndens = dsx;
    if (tsnowc >= 0.)
    {
        dw = 0.13 * dthr / 24.;
        *sndens = *sndens * (1. - dw) + dw;
        if (*sndens >= 0.40)
            *sndens = 0.40;

/*----------------------------------------------------------------------
* calculate snow depth (cm) from snow water equivalent and snow density.
* change snow depth units to meters
* --------------------------------------------------------------------*/
    }
    snowhc = esdc / *sndens;
    *snowh = snowhc * 0.01;

/*----------------------------------------------------------------------
  end subroutine SnowPack
* --------------------------------------------------------------------*/
}

void Snowz0 (double *sncovr, double *z0, double *z0brd, double *snowh)
{

/*----------------------------------------------------------------------
* subroutine Snowz0
* ----------------------------------------------------------------------
* calculate total roughness length over snow
* sncovr  fractional snow cover
* z0      roughness length (m)
* z0s     snow roughness length:=0.001 (m)
* --------------------------------------------------------------------*/
    double          z0s = 0.001;
    double          burial;
    double          z0eff;

    //m z0 = (1.- sncovr)* z0brd + sncovr * z0s
    burial = 7.0 * *z0brd - *snowh;
    if (burial < 0.0007)
        z0eff = z0s;
    else
        z0eff = burial / 7.0;

    *z0 = (1. - *sncovr) * *z0brd + *sncovr * z0eff;

/*----------------------------------------------------------------------
  end subroutine Snowz0
* --------------------------------------------------------------------*/
}

void SnowNew (double *temp, double *newsn, double *snowh, double *sndens)
{

/*----------------------------------------------------------------------
* subroutine SnowNew
* ----------------------------------------------------------------------
* calculate snow depth and density to account for the new snowfall.
* new values of snow depth & density returned.

* temp    air temperature (k)
* newsn   new snowfall (m)
* snowh   snow depth (m)
* sndens  snow density (g/cm3=dimensionless fraction of h2o density)
* --------------------------------------------------------------------*/

    double          dsnew, hnewc, snowhc, newsnc, tempc;

/*----------------------------------------------------------------------
* conversion into simulation units
* --------------------------------------------------------------------*/
    snowhc = *snowh * 100.0;
    newsnc = *newsn * 100.0;

/*----------------------------------------------------------------------
* calculating new snowfall density depending on temperature
* equation from gottlib l. 'a general runoff model for snowcovered
* and glacierized basin', 6th nordic hydrological conference,
* vemadolen, sweden, 1980, 172-177pp.
*---------------------------------------------------------------------*/
    tempc = *temp - 273.15;
    if (tempc <= -15.)
        dsnew = 0.05;
    else
        dsnew = 0.05 + 0.0017 * pow (tempc + 15., 1.5);

/*----------------------------------------------------------------------
* adjustment of snow density depending on new snowfall
* --------------------------------------------------------------------*/
    hnewc = newsnc / dsnew;
    if (snowhc + hnewc < 0.001)
        *sndens = dsnew > *sndens ? dsnew : *sndens;
    else
        *sndens = (snowhc * *sndens + hnewc * dsnew) / (snowhc + hnewc);
    snowhc = snowhc + hnewc;
    *snowh = snowhc * 0.01;

/*----------------------------------------------------------------------
  end subroutine SnowNew
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
SRT (double *rhstt, double *edir, double *et, double *sh2o, double *sh2oa,
   int *nsoil, int *nwtbl, double *pcpdrp, double *zsoil, double *dksat, double *smcmax,
   double *smcmin, double *vgalpha, double *vgbeta, double *macksat,
   double *areaf, int *nmacd, int *mac_status, double *infil, double *runoff2, double *dt,
   double *smcwlt, double *frzx, double *sice, double *ai, double *bi,
   double *ci)
#else
void
SRT (double *rhstt, double *edir, double *et, double *sh2o, double *sh2oa,
   int *nsoil, double *pcpdrp, double *zsoil, double *dwsat, double *dksat,
   double *smcmax, double *bexp, double *runoff1, double *runoff2, double *dt,
   double *smcwlt, double *slope, double *kdt, double *frzx, double *sice,
   double *ai, double *bi, double *ci)
#endif
{

/*----------------------------------------------------------------------
* subroutine SRT
* ----------------------------------------------------------------------
* calculate the right hand side of the time tendency term of the soil
* water diffusion equation.  also to compute ( prepare ) the matrix
* coefficients for the tri-diagonal matrix of the implicit time scheme.
* --------------------------------------------------------------------*/
    int             iohinf;
    int             k, ks;
    double          ddz;
    double          ddz2;
    double          denom;
    double          denom2;
    double          numer;
    double          pddum;
    double          sstt;
    double         *mxsmc, *mxsmc2;
    double         *sicemax;
    double         *wcnd;
    double         *wcnd2;
    double         *wdf;
    double         *wdf2;
#ifdef _NOAH_
    double         *dsmdz, *dsmdz2;
    int             macpore[*nsoil];
#else
    double          dsmdz, dsmdz2;
#endif

/*----------------------------------------------------------------------
* frozen ground version:
* reference frozen ground parameter, cvfrz, is a shape parameter of
* areal distribution function of soil ice content which equals 1/cv.
* cv is a coefficient of spatial variation of soil ice content.  based
* on field data cv depends on areal mean of frozen depth, and it close
* to constant = 0.6 if areal mean frozen depth is above 20 cm.  that is
* why parameter cvfrz = 3 (int{1/0.6*0.6}).
* current logic doesn't allow cvfrz be bigger than 3
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* determine rainfall infiltration rate and runoff.  include the
* infiltration formule from schaake and koren model.
* modified by q duan
* ----------------------------------------------------------------------
* ----------------------------------------------------------------------
* let sicemax be the greatest, if any, frozen water content within soil
* layers.
* --------------------------------------------------------------------*/
    sicemax = (double *)malloc (sizeof (double));
    wcnd = (double *)malloc (sizeof (double));
    wcnd2 = (double *)malloc (sizeof (double));
    wdf = (double *)malloc (sizeof (double));
    wdf2 = (double *)malloc (sizeof (double));
    iohinf = 1;
    *sicemax = 0.0;

    for (ks = 0; ks < *nsoil; ks++)
    {
        if (sice[ks] > *sicemax)
            *sicemax = sice[ks];
/*----------------------------------------------------------------------
* determine rainfall infiltration rate and runoff
* --------------------------------------------------------------------*/
    }

#ifdef _NOAH_
    dsmdz = (double *)malloc (sizeof (double));
    dsmdz2 = (double *)malloc (sizeof (double));

    pddum = *infil;

    for (k = 0; k < *nsoil; k++)
        macpore[k] = 0;
    for (k = 0; k < *nmacd - 1; k++)
        macpore[k] = 1;

///*----------------------------------------------------------------------
//* ys: if lateral runoff (runoff2) is negative (grid is a sink) and the
//* interface layer is close to saturation, lateral runoff is added to the
//* layer above
//* --------------------------------------------------------------------*/
//    denom2 = zsoil[*nsoil - 2] - zsoil[*nsoil - 1];
//    if (*nsoil == 2)
//        denom = -zsoil[*nsoil - 1];
//    else
//        denom = zsoil[*nsoil - 3] - zsoil[*nsoil - 1];
//
//    *dsmdz = (sh2o[*nsoil - 2] - sh2o[*nsoil - 1]) / (denom * 0.5);
//
//    WDfCnd (wdf, wcnd, sh2oa + *nsoil - 2, smcmax, smcmin, vgalpha, vgbeta,
//       dksat, macksat, areaf, mac_status, sicemax, dsmdz, macpore + *nsoil - 2);
//
//    if (*runoff2 < 0
//       && (*smcmax - sh2o[*nsoil - 1]) / *dt * denom2 <
//       *wdf * *dsmdz + *wcnd - et[*nsoil - 1] - *runoff2)
//        rflag = 1;
//    else
//        rflag = 0;

    mxsmc = sh2oa;

    *dsmdz = (sh2o[0] - sh2o[1]) / (-0.5 * zsoil[1]);
    WDfCnd (wdf, wcnd, mxsmc, smcmax, smcmin, vgalpha, vgbeta, dksat, macksat,
       areaf, mac_status, sicemax, dsmdz, macpore);

/*----------------------------------------------------------------------
* calc the matrix coefficients ai, bi, and ci for the top layer
* --------------------------------------------------------------------*/
    ddz = 1. / (-.5 * zsoil[1]);
    ai[0] = 0.0;
    bi[0] = *wdf * ddz / (-zsoil[0]);

/*----------------------------------------------------------------------
* calc rhstt for the top layer after calc'ng the vertical soil moisture
* gradient btwn the top and next to top layers.
* --------------------------------------------------------------------*/
    ci[0] = -bi[0];
    rhstt[0] = (*wdf * *dsmdz + *wcnd - pddum + *edir + et[0]) / zsoil[0];

    if (*nwtbl == 1)
        rhstt[0] += *runoff2 / zsoil[0];

/*----------------------------------------------------------------------
* initialize ddz2
* --------------------------------------------------------------------*/
    sstt = *wdf * *dsmdz + *wcnd + *edir + et[0];

/*----------------------------------------------------------------------
* loop thru the remaining soil layers, repeating the abv process
* --------------------------------------------------------------------*/
    ddz2 = 0.0;
    for (k = 1; k < *nsoil; k++)
    {
        denom2 = (zsoil[k - 1] - zsoil[k]);
        if (k < *nsoil - 1)
        {

/*----------------------------------------------------------------------
* again, to avoid spurious drainage behavior, 'upstream differencing' in
* line below replaced with new approach in 2nd line:
* 'mxsmc2 = max (sh2oa(k), sh2oa(k+1))'
* --------------------------------------------------------------------*/
            mxsmc2 = sh2oa + k;
            denom = (zsoil[k - 1] - zsoil[k + 1]);
            *dsmdz2 = (sh2o[k] - sh2o[k + 1]) / (denom * 0.5);
            WDfCnd (wdf2, wcnd2, mxsmc2, smcmax, smcmin, vgalpha, vgbeta,
               dksat, macksat, areaf, mac_status, sicemax, dsmdz2, macpore + k);

/*-----------------------------------------------------------------------
* calc some partial products for later use in calc'ng rhstt
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* calc the matrix coef, ci, after calc'ng its partial product
* --------------------------------------------------------------------*/
            ddz2 = 2.0 / denom;
            ci[k] = -*wdf2 * ddz2 / denom2;
        }
        else
        {

/*----------------------------------------------------------------------
* slope of bottom layer is introduced
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* retrieve the soil water diffusivity and hydraulic conductivity for
* this layer
* --------------------------------------------------------------------*/
            *wdf2 = 0;
            *wcnd2 = 0;

/*----------------------------------------------------------------------
* calc a partial product for later use in calc'ng rhstt
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* set matrix coef ci to zero
* --------------------------------------------------------------------*/
            *dsmdz2 = 0.0;
            ci[k] = 0.0;

/*----------------------------------------------------------------------
* calc rhstt for this layer after calc'ng its numerator
* --------------------------------------------------------------------*/
        }

        numer = (*wdf2 * *dsmdz2) + *wcnd2 - (*wdf * *dsmdz) - *wcnd + et[k];
        if (k == *nwtbl - 1)
                numer = numer + *runoff2;

/*----------------------------------------------------------------------
* calc matrix coefs, ai, and bi for this layer
* --------------------------------------------------------------------*/
        rhstt[k] = numer / (-denom2);
        ai[k] = -*wdf * ddz / denom2;

/*----------------------------------------------------------------------
* reset values of wdf, wcnd, dsmdz, and ddz for loop to next lyr
* runoff2:  sub-surface or baseflow runoff
* --------------------------------------------------------------------*/
        bi[k] = -(ai[k] + ci[k]);

        if (k != *nsoil - 1)
        {
            *wdf = *wdf2;
            *wcnd = *wcnd2;
            *dsmdz = *dsmdz2;
            ddz = ddz2;
        }
    }

    free (dsmdz);
    free (dsmdz2);
#else
    pddum = *pcpdrp;
    *runoff1 = 0.0;

/*----------------------------------------------------------------------
* modified by q. duan, 5/16/94
* --------------------------------------------------------------------*/
    //        if (iohinf == 1) then

    if (*pcpdrp != 0.0)
    {
        dt1 = *dt / 86400.;
        smcav = *smcmax - *smcwlt;

/*----------------------------------------------------------------------
* frozen ground version:
* --------------------------------------------------------------------*/
        dmax[0] = -zsoil[0] * smcav;

        dice = -zsoil[0] * sice[0];
        dmax[0] = dmax[0] * (1.0 - (sh2oa[0] + sice[0] - *smcwlt) / smcav);

        dd = dmax[0];

/*----------------------------------------------------------------------
* frozen ground version:
* --------------------------------------------------------------------*/
        for (ks = 1; ks < *nsoil; ks++)
        {
            dice = dice + (zsoil[ks - 1] - zsoil[ks]) * sice[ks];
            dmax[ks] = (zsoil[ks - 1] - zsoil[ks]) * smcav;
            dmax[ks] =
               dmax[ks] * (1.0 - (sh2oa[ks] + sice[ks] - *smcwlt) / smcav);
            dd = dd + dmax[ks];

/*----------------------------------------------------------------------
* val = (1.-exp(-kdt*sqrt(dt1)))
* in below, remove the sqrt in above
* --------------------------------------------------------------------*/
        }
        val = (1. - exp (-*kdt * dt1));
        ddt = dd * val;
        px = *pcpdrp * *dt;
        if (px < 0.0)
            px = 0.0;

/*----------------------------------------------------------------------
* frozen ground version:
* reduction of infiltration based on frozen ground parameters
* --------------------------------------------------------------------*/
        infmax = (px * (ddt / (px + ddt))) / *dt;
        fcr = 1.;
        if (dice > 1.e-2)
        {
            acrt = (double)cvfrz **frzx / dice;
            sum = 1.;
            ialp1 = cvfrz - 1;
            for (j = 1; j < ialp1 + 1; j++)
            {
                k = 1;
                for (jj = j + 1; jj < ialp1; jj++)
                {
                    k = k * jj;
                }
                sum = sum + pow (acrt, (double)(cvfrz - j)) / (double)k;
            }
            fcr = 1. - exp (-acrt) * sum;
        }

/*----------------------------------------------------------------------
* correction of infiltration limitation:
* if infmax .le. hydrolic conductivity assign infmax the value of
* hydrolic conductivity
* ---------------------------------------------------------------------*/
        //         mxsmc = max ( sh2oa(1), sh2oa(2) )
        infmax = infmax * fcr;
        mxsmc = sh2oa;
        WDfCnd (wdf, wcnd, mxsmc, smcmax, bexp, dksat, dwsat, sicemax);
        infmax = infmax > *wcnd ? infmax : *wcnd;

        infmax = infmax < px / *dt ? infmax : px / *dt;
        if (*pcpdrp > infmax)
        {
            *runoff1 = *pcpdrp - infmax;
            pddum = infmax;
        }

/*----------------------------------------------------------------------
* to avoid spurious drainage behavior, 'upstream differencing' in line
* below replaced with new approach in 2nd line:
* 'mxsmc = max(sh2oa(1), sh2oa(2))'
* --------------------------------------------------------------------*/
    }

    mxsmc = sh2oa;
    WDfCnd (wdf, wcnd, mxsmc, smcmax, bexp, dksat, dwsat, sicemax);

/*----------------------------------------------------------------------
* calc the matrix coefficients ai, bi, and ci for the top layer
* --------------------------------------------------------------------*/
    ddz = 1. / (-.5 * zsoil[1]);
    ai[0] = 0.0;
    bi[0] = *wdf * ddz / (-zsoil[0]);

/*----------------------------------------------------------------------
* calc rhstt for the top layer after calc'ng the vertical soil moisture
* gradient btwn the top and next to top layers.
* --------------------------------------------------------------------*/
    ci[0] = -bi[0];
    dsmdz = (sh2o[0] - sh2o[1]) / (-0.5 * zsoil[1]);
    rhstt[0] = (*wdf * dsmdz + *wcnd - pddum + *edir + et[0]) / zsoil[0];

/*----------------------------------------------------------------------
* initialize ddz2
* --------------------------------------------------------------------*/
    sstt = *wdf * dsmdz + *wcnd + *edir + et[0];

/*----------------------------------------------------------------------
* loop thru the remaining soil layers, repeating the abv process
* --------------------------------------------------------------------*/
    ddz2 = 0.0;
    for (k = 1; k < *nsoil; k++)
    {
        denom2 = (zsoil[k - 1] - zsoil[k]);
        if (k != *nsoil - 1)
        {

/*----------------------------------------------------------------------
* again, to avoid spurious drainage behavior, 'upstream differencing' in
* line below replaced with new approach in 2nd line:
* 'mxsmc2 = max (sh2oa(k), sh2oa(k+1))'
* --------------------------------------------------------------------*/
            slopx = 1.;

            mxsmc2 = sh2oa + k;
            WDfCnd (wdf2, wcnd2, mxsmc2, smcmax, bexp, dksat, dwsat, sicemax);

/*-----------------------------------------------------------------------
* calc some partial products for later use in calc'ng rhstt
* --------------------------------------------------------------------*/
            denom = (zsoil[k - 1] - zsoil[k + 1]);

/*----------------------------------------------------------------------
* calc the matrix coef, ci, after calc'ng its partial product
* --------------------------------------------------------------------*/
            dsmdz2 = (sh2o[k] - sh2o[k + 1]) / (denom * 0.5);
            ddz2 = 2.0 / denom;
            ci[k] = -*wdf2 * ddz2 / denom2;
        }
        else
        {

/*----------------------------------------------------------------------
* slope of bottom layer is introduced
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* retrieve the soil water diffusivity and hydraulic conductivity for
* this layer
* --------------------------------------------------------------------*/
            slopx = *slope;
            WDfCnd (wdf2, wcnd2, sh2oa + *nsoil - 1, smcmax, bexp, dksat,
               dwsat, sicemax);

/*----------------------------------------------------------------------
* calc a partial product for later use in calc'ng rhstt
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* set matrix coef ci to zero
* --------------------------------------------------------------------*/
            dsmdz2 = 0.0;
            ci[k] = 0.0;

/*----------------------------------------------------------------------
* calc rhstt for this layer after calc'ng its numerator
* --------------------------------------------------------------------*/
        }
        numer =
           (*wdf2 * dsmdz2) + slopx * *wcnd2 - (*wdf * dsmdz) - *wcnd + et[k];

/*----------------------------------------------------------------------
* calc matrix coefs, ai, and bi for this layer
* --------------------------------------------------------------------*/
        rhstt[k] = numer / (-denom2);
        ai[k] = -*wdf * ddz / denom2;

/*----------------------------------------------------------------------
* reset values of wdf, wcnd, dsmdz, and ddz for loop to next lyr
* runoff2:  sub-surface or baseflow runoff
* --------------------------------------------------------------------*/
        bi[k] = -(ai[k] + ci[k]);
        if (k == *nsoil - 1)
            *runoff2 = slopx * *wcnd2;
        if (k != *nsoil - 1)
        {
            *wdf = *wdf2;
            *wcnd = *wcnd2;
            dsmdz = dsmdz2;
            ddz = ddz2;
        }
    }
#endif

    free (sicemax);
    free (wcnd);
    free (wcnd2);
    free (wdf);
    free (wdf2);

/*----------------------------------------------------------------------
  end subroutine SRT
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
SStep (double *sh2oout, double *sh2oin, double *cmc, double *rhstt,
   double *rhsct, double *dt, int *nsoil, double *smcmax, double *smcmin,
   double *cmcmax, double *runoff3, double *zsoil, double *smc, double *sice,
   double *ai, double *bi, double *ci)
#else
void
SStep (double *sh2oout, double *sh2oin, double *cmc, double *rhstt,
   double *rhsct, double *dt, int *nsoil, double *smcmax, double *cmcmax,
   double *runoff3, double *zsoil, double *smc, double *sice, double *ai,
   double *bi, double *ci)
#endif
{

/*----------------------------------------------------------------------
* subroutine SStep
* ----------------------------------------------------------------------
* calculate/update soil moisture content values and canopy moisture
* content values.
* --------------------------------------------------------------------*/
    int             k, kk11;

    double          rhsttin[*nsoil], ciin[*nsoil];
    double          sh2omid[*nsoil];
    double          ddz, stot, wplus;

/*----------------------------------------------------------------------
* create 'amount' values of variables to be input to the
* tri-diagonal matrix routine.
* --------------------------------------------------------------------*/
    for (k = 0; k < *nsoil; k++)
    {
        rhstt[k] = rhstt[k] * *dt;
        ai[k] = ai[k] * *dt;
        bi[k] = 1. + bi[k] * *dt;
        ci[k] = ci[k] * *dt;
    }

/*----------------------------------------------------------------------
* copy values for input variables before call to Rosr12
* --------------------------------------------------------------------*/
    for (k = 0; k < *nsoil; k++)
        rhsttin[k] = rhstt[k];
    for (k = 0; k < *nsoil; k++)
        ciin[k] = ci[k];

/*----------------------------------------------------------------------
* call Rosr12 to solve the tri-diagonal matrix
* --------------------------------------------------------------------*/
    Rosr12 (ci, ai, bi, ciin, rhsttin, rhstt, nsoil);

/*----------------------------------------------------------------------
* sum the previous smc value and the matrix solution to get a
* new value.  min allowable value of smc will be 0.02.
* runoff3: runoff within soil layers
* --------------------------------------------------------------------*/
    wplus = 0.0;
    *runoff3 = 0.;

    for (k = *nsoil - 1; k >= 0; k--)
    {
        if (k != 0)
            ddz = zsoil[k - 1] - zsoil[k];
        else
            ddz = -zsoil[0];

        sh2omid[k] = sh2oin[k] + ci[k] + wplus / ddz;
        stot = sh2omid[k] + sice[k];

        if (stot > *smcmax)
        {
            if (k == 0)
                ddz = -zsoil[0];
            else
            {
                kk11 = k - 1;
                ddz = -zsoil[k] + zsoil[kk11];
            }
            wplus = (stot - *smcmax) * ddz;
        }
        else
            wplus = 0.0;

        if (stot < *smcmax)
            smc[k] = stot;
        else
            smc[k] = *smcmax;

        sh2omid[k] = smc[k] - sice[k];
    }

    ddz = -zsoil[0];
    for (k = 0; k < *nsoil; k++)
    {
        if (k != 0)
            ddz = zsoil[k - 1] - zsoil[k];
        sh2oout[k] = sh2omid[k] + wplus / ddz;
        stot = sh2oout[k] + sice[k];
        if (stot > *smcmax)
        {
            if (k == 0)
                ddz = -zsoil[0];
            else
            {
                kk11 = k - 1;
                ddz = -zsoil[k] + zsoil[kk11];
            }
            wplus = (stot - *smcmax) * ddz;
        }
        else
            wplus = 0.;

        smc[k] = stot < *smcmax ? stot : *smcmax;
#ifdef _NOAH_
        smc[k] = smc[k] > *smcmin + 0.02 ? smc[k] : *smcmin + 0.02;
#else
        smc[k] = smc[k] > 0.02 ? smc[k] : 0.02;
#endif
        sh2oout[k] = smc[k] - sice[k];
        sh2oout[k] = sh2oout[k] > 0 ? sh2oout[k] : 0;
    }

/*----------------------------------------------------------------------
* update canopy water content/interception (cmc).  convert rhsct to
* an 'amount' value and add to previous cmc value to get new cmc.
* --------------------------------------------------------------------*/
    *runoff3 = wplus;
    *cmc = *cmc + *dt * *rhsct;
    if (*cmc < 1.e-20)
        *cmc = 0.0;
    *cmc = *cmc < *cmcmax ? *cmc : *cmcmax;

/*----------------------------------------------------------------------
  end subroutine SStep
* --------------------------------------------------------------------*/
}

void
TBnd (double *tu, double *tb, double *zsoil, double *zbot, int k, int *nsoil,
   double *tbnd1)
{

/*----------------------------------------------------------------------
* subroutine TBnd
* ----------------------------------------------------------------------
* calculate temperature on the boundary of the layer by interpolation of
* the middle layer temperatures
* --------------------------------------------------------------------*/
    double          zb, zup;

/*----------------------------------------------------------------------
* use surface temperature on the top of the first layer
* --------------------------------------------------------------------*/
    if (k == 0)
        zup = 0.;
    else
        zup = zsoil[k - 1];

/*----------------------------------------------------------------------
* use depth of the constant bottom temperature when interpolate
* temperature into the last layer boundary
* --------------------------------------------------------------------*/
    if (k == *nsoil - 1)
        zb = 2. * *zbot - zsoil[k];
    else
        zb = zsoil[k + 1];

/*----------------------------------------------------------------------
* linear interpolation between the average layer temperatures
* --------------------------------------------------------------------*/

    *tbnd1 = *tu + (*tb - *tu) * (zup - zsoil[k]) / (zup - zb);

/*----------------------------------------------------------------------
  end subroutine TBnd
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
TDfCnd (double *df, double *smc, double *qz, double *smcmax, double *smcmin,
   double *sh2o)
#else
void
TDfCnd (double *df, double *smc, double *qz, double *smcmax, double *sh2o)
#endif
{

/*----------------------------------------------------------------------
* subroutine TDfCnd
* ----------------------------------------------------------------------
* calculate thermal diffusivity and conductivity of the soil for a given
* point and time.
* ----------------------------------------------------------------------
* peters-lidard approach (peters-lidard et al., 1998)
* june 2001 changes: frozen soil condition.
* --------------------------------------------------------------------*/
    double          ake, gammd, thkdry, thkice, thko, thkqtz, thksat, thks,
       thkw, satratio, xu, xunfroz;

/*----------------------------------------------------------------------
* we now get quartz as an input argument (set in routine RedPrm):
*      data quartz /0.82, 0.10, 0.25, 0.60, 0.52,
*     &             0.35, 0.60, 0.40, 0.82/
* ----------------------------------------------------------------------
* if the soil has any moisture content compute a partial sum/product
* otherwise use a constant value which works well with most soils
* ----------------------------------------------------------------------
*  thkw ......water thermal conductivity
*  thkqtz ....thermal conductivity for quartz
*  thko ......thermal conductivity for other soil components
*  thks ......thermal conductivity for the solids combined(quartz+other)
*  thkice ....ice thermal conductivity
*  smcmax ....porosity (= smcmax)
*  qz .........quartz content (soil type dependent)
* ----------------------------------------------------------------------
* use as in peters-lidard, 1998 (modif. from johansen, 1975).

*                                  pablo grunmann, 08/17/98
* refs.:
*      farouki, o.t.,1986: thermal properties of soils. series on rock
*              and soil mechanics, vol. 11, trans tech, 136 pp.
*      johansen, o., 1975: thermal conductivity of soils. ph.d. thesis,
*              university of trondheim,
*      peters-lidard, c. d., et al., 1998: the effect of soil thermal
*              conductivity parameterization on surface energy fluxes
*              and temperatures. journal of the atmospheric sciences,
*              vol. 55, pp. 1209-1224.
* --------------------------------------------------------------------*/
    // needs parameters
    // porosity(soil type):
    //      poros = smcmax
    // saturation ratio:
    // parameters  w/(m.k)
#ifdef _NOAH_
    satratio = (*smc - *smcmin) / (*smcmax - *smcmin);
#else
    satratio = *smc / *smcmax;
#endif

    /*
     * ice conductivity: 
     */
    thkice = 2.2;

    /*
     * water conductivity: 
     */
    thkw = 0.57;

    /*
     * thermal conductivity of "other" soil components 
     */
    //      if (qz .le. 0.2) thko = 3.0
    thko = 2.0;

    /*
     * quartz' conductivity 
     */
    thkqtz = 7.7;

    /*
     * solids' conductivity 
     */
    thks = pow (thkqtz, *qz) * pow (thko, 1. - *qz);

    /*
     * unfrozen fraction (from 1., i.e., 100%liquid, to 0. (100% frozen)) 
     */
    xunfroz = *sh2o / *smc;

    /*
     * unfrozen volume for saturation (porosity*xunfroz) 
     */
    xu = xunfroz * (*smcmax);

    /*
     * saturated thermal conductivity 
     */
    thksat =
       pow (thks, 1. - *smcmax) * pow (thkice, *smcmax - xu) * pow (thkw, xu);

    /*
     * dry density in kg/m3 
     */
    gammd = (1. - *smcmax) * 2700.;

    /*
     * dry thermal conductivity in w.m-1.k-1 
     */
    thkdry = (0.135 * gammd + 64.7) / (2700. - 0.947 * gammd);

    /*
     * frozen 
     */
    if ((*sh2o + 0.0005) < *smc)
        ake = satratio;

    /*
     * unfrozen
     * range of validity for the kersten number (ake)
     */
    else
    {

        /*
         * kersten number (using "fine" formula, valid for soils containing at
         * least 5% of particles with diameter less than 2.e-6 meters.)
         * (for "coarse" formula, see peters-lidard et al., 1998).
         */

        if (satratio > 0.1)
            ake = log10 (satratio) + 1.0;

        /*
         * use k = kdry 
         */
        else
            ake = 0.0;
    }

    /*
     * thermal conductivity 
     */

    *df = ake * (thksat - thkdry) + thkdry;

/*----------------------------------------------------------------------
  end subroutine TDfCnd
* --------------------------------------------------------------------*/
}

void
TmpAvg (double *tavg, double *tup, double *tm, double *tdn, double *zsoil,
   int *nsoil, int k)
{

/*----------------------------------------------------------------------
* subroutine TmpAvg
* ----------------------------------------------------------------------
* calculate soil layer average temperature (tavg) in freezing/thawing
* layer using up, down, and middle layer temperatures (tup, tdn, tm),
* where tup is at top boundary of layer, tdn is at bottom boundary of
* layer.  tm is layer prognostic state temperature.
* --------------------------------------------------------------------*/
    double          dz, dzh, x0, xdn, xup;
    double          t0 = 2.7315e2;

/*--------------------------------------------------------------------*/

    if (k == 0)
        dz = -zsoil[0];
    else
        dz = zsoil[k - 1] - zsoil[k];

    dzh = dz * 0.5;
    if (*tup < t0)
    {
        if (*tm < t0)
        {

/*----------------------------------------------------------------------
* tup, tm, tdn < t0
* --------------------------------------------------------------------*/
            if (*tdn < t0)
                *tavg = (*tup + 2.0 * *tm + *tdn) / 4.0;

/*----------------------------------------------------------------------
* tup & tm < t0,  tdn .ge. t0
* --------------------------------------------------------------------*/
            else
            {
                x0 = (t0 - *tm) * dzh / (*tdn - *tm);
                *tavg =
                   0.5 * (*tup * dzh + *tm * (dzh + x0) + t0 * (2. * dzh -
                      x0)) / dz;
            }
        }
        else
        {

/*----------------------------------------------------------------------
* tup < t0, tm .ge. t0, tdn < t0
* --------------------------------------------------------------------*/
            if (*tdn < t0)
            {
                xup = (t0 - *tup) * dzh / (*tm - *tup);
                xdn = dzh - (t0 - *tm) * dzh / (*tdn - *tm);
                *tavg =
                   0.5 * (*tup * xup + t0 * (2. * dz - xup - xdn) +
                   *tdn * xdn) / dz;
            }

/*----------------------------------------------------------------------
* tup < t0, tm .ge. t0, tdn .ge. t0
* --------------------------------------------------------------------*/
            else
            {
                xup = (t0 - *tup) * dzh / (*tm - *tup);
                *tavg = 0.5 * (*tup * xup + t0 * (2. * dz - xup)) / dz;
            }
        }
    }
    else
    {
        if (*tm < t0)
        {

/*----------------------------------------------------------------------
* tup .ge. t0, tm < t0, tdn < t0
* --------------------------------------------------------------------*/
            if (*tdn < t0)
            {
                xup = dzh - (t0 - *tup) * dzh / (*tm - *tup);
                *tavg =
                   0.5 * (t0 * (dz - xup) + *tm * (dzh + xup) +
                   *tdn * dzh) / dz;
            }

/*----------------------------------------------------------------------
* tup .ge. t0, tm < t0, tdn .ge. t0
* --------------------------------------------------------------------*/
            else
            {
                xup = dzh - (t0 - *tup) * dzh / (*tm - *tup);
                xdn = (t0 - *tm) * dzh / (*tdn - *tm);
                *tavg =
                   0.5 * (t0 * (2. * dz - xup - xdn) + *tm * (xup +
                      xdn)) / dz;
            }
        }
        else
        {

/*----------------------------------------------------------------------
* tup .ge. t0, tm .ge. t0, tdn < t0
* --------------------------------------------------------------------*/
            if (*tdn < t0)
            {
                xdn = dzh - (t0 - *tm) * dzh / (*tdn - *tm);
                *tavg = (t0 * (dz - xdn) + 0.5 * (t0 + *tdn) * xdn) / dz;
            }

/*----------------------------------------------------------------------
* tup .ge. t0, tm .ge. t0, tdn .ge. t0
* --------------------------------------------------------------------*/
            else
                *tavg = (*tup + 2.0 * *tm + *tdn) / 4.0;
        }
    }

/*----------------------------------------------------------------------
  end subroutine TmpAvg
* --------------------------------------------------------------------*/
}

void
Transp (double *et, int *nsoil, double *etp1, double *smc, double *cmc,
   double *zsoil, double *shdfac, double *smcwlt, double *cmcmax, double *pc,
   double *cfactr, double *smcref, double *sfctmp, double *q2, int *nroot,
   double *rtdis)
{

/*----------------------------------------------------------------------
* subroutine Transp
* ----------------------------------------------------------------------
* calculate transpiration for the veg class.
* --------------------------------------------------------------------*/
    int             i, k;
    double          denom;
    double          etp1a;
    //.....real part(nsoil)
    double          gx[*nroot];
    double          rtx, sgx;

/*----------------------------------------------------------------------
* initialize plant transp to zero for all soil layers.
* --------------------------------------------------------------------*/
    for (k = 0; k < *nsoil; k++)
        et[k] = 0.;

/*----------------------------------------------------------------------
* calculate an 'adjusted' potential transpiration
* if statement below to avoid tangent linear problems near zero
* note: gx and other terms below redistribute transpiration by layer,
* et(k), as a function of soil moisture availability, while preserving
* total etp1a.
* --------------------------------------------------------------------*/
    if (*cmc != 0.0)
        etp1a = *shdfac * *pc * *etp1 * (1.0 - pow (*cmc / *cmcmax, *cfactr));
    else
        etp1a = *shdfac * *pc * *etp1;
    sgx = 0.0;
    for (i = 0; i < *nroot; i++)
    {
        gx[i] = (smc[i] - *smcwlt) / (*smcref - *smcwlt);
        if (gx[i] < 0)
            gx[i] = 0;
        if (gx[i] > 1)
            gx[i] = 1.0;
        sgx = sgx + gx[i];
    }

    sgx = sgx / (double)*nroot;
    denom = 0.;
    for (i = 0; i < *nroot; i++)
    {
        rtx = rtdis[i] + gx[i] - sgx;
        gx[i] = gx[i] * (rtx > 0. ? rtx : 0.);
        denom = denom + gx[i];
    }

    if (denom <= 0.0)
        denom = 1.;
    for (i = 0; i < *nroot; i++)
    {
        et[i] = etp1a * gx[i] / denom;

/*----------------------------------------------------------------------
* above code assumes a vertically uniform root distribution
* code below tests a variable root distribution
* ----------------------------------------------------------------------
*		et[0] = (zsoil[0] / zsoil(*nroot) ) * gx * etp1a;
*		et[0] = (zsoil[0] / zsoil(*nroot) ) * etp1a;
* ----------------------------------------------------------------------
* using root distribution as weighting factor
* ----------------------------------------------------------------------
*		et[0] = rtdis[0] * etp1a;
*		et[0] = etp1a * part[0];
* ----------------------------------------------------------------------
* loop down thru the soil layers repeating the operation above,
* but using the thickness of the soil layer (rather than the
* absolute depth of each layer) in the final calculation.
* ----------------------------------------------------------------------
*		for (k = 1; k < *nroot; k++)
		{
*			gx = (smc[k] - *smcwlt ) / (*smcref - *smcwlt);
*			gx = gx < 1.0 ? gx : 1.0;
*			gx = gx > 0.0 ? gx : 0.0;
* test canopy resistance
*			gx = 1.0;
*			et[k] = ((zsoil[k] - zsoil[k-1]) / zsoil[*nroot - 1]) * gx * etp1a;
*			et[k] = ((zsoil[k] - zsoil[k-1]) / zsoil[*nroot - 1]) * etp1a;
* ----------------------------------------------------------------------
* using root distribution as weighting factor
* ----------------------------------------------------------------------
*			et[k] = rtdis[k] * etp1a;
*			et[k] = etp1a * part[k];
*		}
*/
    }

/*----------------------------------------------------------------------
  end subroutine Transp
* --------------------------------------------------------------------*/
}

#ifdef _NOAH_
void
WDfCnd (double *wdf, double *wcnd, double *smc, double *smcmax, double *smcmin, double *vgalpha, double *vgbeta, double *dksat, double *macksat, double *areaf, int *mac_status, double *sicemax, double *dsmdz, int *macpore)
#else
void
WDfCnd (double *wdf, double *wcnd, double *smc, double *smcmax, double *bexp, double *dksat, double *dwsat, double *sicemax)
#endif
{

/*----------------------------------------------------------------------
* subroutine WDfCnd
* ----------------------------------------------------------------------
* calculate soil water diffusivity and soil hydraulic conductivity.
* --------------------------------------------------------------------*/
    double          expon, factr1, factr2, vkwgt;
#ifdef _NOAH_
    double          satkfunc, dpsidsm;
#endif

/*----------------------------------------------------------------------
*     calc the ratio of the actual to the max psbl soil h2o content
* --------------------------------------------------------------------*/

#ifdef _NOAH_
    factr1 = 0.05 / (*smcmax - *smcmin);
    factr2 = (*smc - *smcmin) / (*smcmax - *smcmin);

/*----------------------------------------------------------------------
* factr2 should avoid to be 0 or 1
* --------------------------------------------------------------------*/
    if (factr2 > 1. - .0005)
        factr2 = 1. - .0005;
    if (factr2 < 0. + .0005)
        factr2 = .0005;
    factr1 = factr1 < factr2 ? factr1 : factr2;
    expon = 1.0 - 1. / *vgbeta;

    satkfunc =
       pow (factr2, 0.5) * pow (1. - pow (1. - pow (factr2, 1. / expon),
          expon), 2.);
    dpsidsm =
       (1. - expon) / *vgalpha / expon / (*smcmax -
       *smcmin) * pow (pow (factr2, -1. / expon) - 1.,
       0. - expon) * pow (factr2, -(1. / expon + 1.));

    if (*macpore == 1)
        *wcnd = EFFKV (satkfunc, factr2, *mac_status, *macksat, *dksat, *areaf);
    else
        *wcnd = *dksat * satkfunc;

    *wdf = *wcnd * dpsidsm;

    //  *wdf = (1. - expon) * (*dksat * (1. - *areaf) + *macksat * *areaf) / *vgalpha / expon / (*smcmax - *smcmin) * pow(factr2, 0.5 - 1. / expon) * (pow(1. - pow(factr2, 1. / expon) ,-expon) + pow(1. - pow(factr2, 1. / expon) ,expon) - 2.);
    //  *wcnd = sqrt(factr2) * pow(1. - pow(1. - pow(factr2, 1. / expon), expon), 2.) * (*dksat * (1. - *areaf) + *macksat * *areaf);

    if (*sicemax > 0.0)
    {
        vkwgt = 1. / (1. + pow (500. * *sicemax, 3.));
        satkfunc =
           pow (factr1, 0.5) * pow (1. - pow (1. - pow (factr1, 1. / expon),
              expon), 2.);
        dpsidsm =
           (1. - expon) / *vgalpha / expon / (*smcmax -
           *smcmin) * pow (pow (factr1, -1. / expon) - 1.,
           0. - expon) * pow (factr1, -(1. / expon + 1.));
        if (*macpore == 1)
            *wdf =
               vkwgt * *wdf + (1. - vkwgt) * dpsidsm * EFFKV (satkfunc, factr1, *mac_status, *macksat, *dksat, *areaf);
        else
            *wdf = vkwgt * *wdf + (1. - vkwgt) * dpsidsm * satkfunc * *dksat;
    }


#else

/*----------------------------------------------------------------------
* prep an expntl coef and calc the soil water diffusivity
* --------------------------------------------------------------------*/
    factr1 = 0.05 / *smcmax;
    factr2 = *smc / *smcmax;
    factr1 = factr1 < factr2 ? factr1 : factr2;
    expon = *bexp + 2.0;

/*----------------------------------------------------------------------
* frozen soil hydraulic diffusivity.  very sensitive to the vertical
* gradient of unfrozen water. the latter gradient can become very
* extreme in freezing/thawing situations, and given the relatively
* few and thick soil layers, this gradient sufferes serious
* trunction errors yielding erroneously high vertical transports of
* unfrozen water in both directions from huge hydraulic diffusivity.
* therefore, we found we had to arbitrarily constrain wdf
* --
* version d_10cm: ........  factr1 = 0.2/smcmax
* weighted approach...................... pablo grunmann, 28_sep_1999.
* --------------------------------------------------------------------*/
    *wdf = *dwsat * pow (factr2, expon);
    if (*sicemax > 0.0)
    {
        vkwgt = 1. / (1. + pow (500. * *sicemax, 3.));
        *wdf = vkwgt * *wdf + (1. - vkwgt) * *dwsat * pow (factr1, expon);

/*----------------------------------------------------------------------
* reset the expntl coef and calc the hydraulic conductivity
* --------------------------------------------------------------------*/
    }
    expon = (2.0 * *bexp) + 3.0;
    *wcnd = *dksat * pow (factr2, expon);
#endif

/*----------------------------------------------------------------------
  end subroutine WDfCnd
* --------------------------------------------------------------------*/
}


void SfcDifOff (double *zlm, double *zlm_wind, double *z0, double *thz0, double *thlm, double *sfcspd, double *czil, double *akms, double *akhs, int *vegtyp, int *isurban, int *iz0tlnd)
{

/*----------------------------------------------------------------------
* subroutine sfcdif (renamed SfcDifOff to avoid clash with eta pbl)
* ----------------------------------------------------------------------
* calculate surface layer exchange coefficients via iterative process.
* see chen et al (1997, blm)
* --------------------------------------------------------------------*/

    double          zilfc, zu, zt, rdz, cxch;
    double          dthv, du2, btgh, wstar2, ustar, zslu, zslt, rlogu, rlogt;
    double          rlmo, zetalt, zetalu, zetau, zetat, xlu4, xlt4, xu4, xt4;
    //!cc   ......real ztfc

    double          xlu, xlt, xu, xt, psmz, simm, pshz, simh, ustark, rlmn,
       rlma;

    int             ilech, itr;

    double          wwst = 1.2;
    double          wwst2;
    double          g = 9.8, vkrm = 0.40, excm = 0.001, beta = 1. / 270.;
    double          btg, elfc;
    double          wold = .15;
    double          wnew;
    int             itrmx = 5;

    double          epsu2 = 1.e-4;
    double          epsust = 0.07;
    //  double epsa = 1.e-8;
    double          ztmin = -5.;
    double          ztmax = 1.;
    double          hpbl = 1000.0;
    double          sqvisc = 258.2;

    wwst2 = wwst * wwst;
    btg = beta * g;
    elfc = vkrm * btg;
    wnew = 1. - wold;

/*----------------------------------------------------------------------
*     ztfc: ratio of zoh/zom  less or equal than 1
*     c......ztfc=0.1
*     czil: constant c in zilitinkevich, s. s.1995,:note about zt
* --------------------------------------------------------------------*/
    ilech = 0;

/*--------------------------------------------------------------------*/
    if ((*iz0tlnd == 0) || (*vegtyp == *isurban))
    {
        /* just use the original czil value. */
        zilfc = -*czil * vkrm * sqvisc;
    }
    else
    {
        /* modify czil according to chen & zhang, 2009
         * czil = 10 ** -0.40 h, ( where h = 10*zo ) */
        *czil = pow (10.0, -0.4 * (*z0 / 0.07));
        zilfc = -*czil * vkrm * sqvisc;
    }
    //     c.......zt=z0*ztfc
    zu = *z0;
    rdz = 1. / *zlm_wind;
    cxch = excm * rdz;
    dthv = *thlm - *thz0;

/*----------------------------------------------------------------------
* beljars correction of ustar
* --------------------------------------------------------------------*/
    du2 = (*sfcspd * *sfcspd) > epsu2 ? (*sfcspd * *sfcspd) : epsu2;

    /*
     * cc   if statements to avoid tangent linear problems near zero 
     */
    btgh = btg * hpbl;
    if (btgh * *akhs * dthv != 0.0)
        wstar2 = wwst2 * pow (fabs (btgh * *akhs * dthv), 2. / 3.);
    else
        wstar2 = 0.0;

/*----------------------------------------------------------------------
* zilitinkevitch approach for zt
* --------------------------------------------------------------------*/
    ustar = sqrt (*akms * sqrt (du2 + wstar2));
    ustar = ustar > epsust ? ustar : epsust;

/*--------------------------------------------------------------------*/
    zt = exp (zilfc * sqrt (ustar * *z0)) * *z0;
    zslu = *zlm_wind + zu;
    //     print*,'zslt=',zslt
    //     print*,'zlm=',zlm
    //     print*,'zt=',zt

    zslt = *zlm + zt;
//    zslt = *zlm_wind + zt;
    rlogu = log (zslu / zu);

    rlogt = log (zslt / zt);
    //     print*,'rlmo=',rlmo
    //     print*,'elfc=',elfc
    //     print*,'akhs=',akhs
    //     print*,'dthv=',dthv
    //     print*,'ustar=',ustar

    rlmo = elfc * *akhs * dthv / pow (ustar, 3);

/*----------------------------------------------------------------------
* 1./monin-obukkhov length-scale
* --------------------------------------------------------------------*/
    for (itr = 0; itr < itrmx; itr++)
    {
        zetalt = zslt * rlmo;
        zetalt = zetalt > ztmin ? zetalt : ztmin;
        rlmo = zetalt / zslt;
        zetalu = zslu * rlmo;
        zetau = zu * rlmo;

        zetat = zt * rlmo;
        if (ilech == 0)
        {
            if (rlmo < 0.)
            {
                xlu4 = 1. - 16. * zetalu;
                xlt4 = 1. - 16. * zetalt;
                xu4 = 1. - 16. * zetau;
                xt4 = 1. - 16. * zetat;
                xlu = sqrt (sqrt (xlu4));
                xlt = sqrt (sqrt (xlt4));
                xu = sqrt (sqrt (xu4));

                xt = sqrt (sqrt (xt4));
                //     print*,'-----------1------------'
                //     print*,'psmz=',psmz
                //     print*,'Pspmu(zetau)=',Pspmu(zetau)
                //     print*,'xu=',xu
                //     print*,'------------------------'
                psmz = Pspmu (xu);
                simm = Pspmu (xlu) - psmz + rlogu;
                pshz = Psphu (xt);
                simh = Psphu (xlt) - pshz + rlogt;
            }
            else
            {
                zetalu = zetalu < ztmax ? zetalu : ztmax;
                zetalt = zetalt < ztmax ? zetalt : ztmax;
                //     print*,'-----------2------------'
                //     print*,'psmz=',psmz
                //     print*,'Pspms(zetau)=',Pspms(zetau)
                //     print*,'zetau=',zetau
                //     print*,'------------------------'
                psmz = Pspms (zetau);
                simm = Pspms (zetalu) - psmz + rlogu;
                pshz = Psphs (zetat);
                simh = Psphs (zetalt) - pshz + rlogt;
            }
        }

/*----------------------------------------------------------------------
* lech's functions
* --------------------------------------------------------------------*/
        else
        {
            if (rlmo < 0.)
            {
                //     print*,'-----------3------------'
                //     print*,'psmz=',psmz
                //     print*,'Pslmu(zetau)=',Pslmu(zetau)
                //     print*,'zetau=',zetau
                //     print*,'------------------------'
                psmz = Pslmu (zetau);
                simm = Pslmu (zetalu) - psmz + rlogu;
                pshz = Pslhu (zetat);
                simh = Pslhu (zetalt) - pshz + rlogt;
            }
            else
            {
                zetalu = zetalu < ztmax ? zetalu : ztmax;
                zetalt = zetalt < ztmax ? zetalt : ztmax;
                //     print*,'-----------4------------'
                //     print*,'psmz=',psmz
                //     print*,'Pslms(zetau)=',Pslms(zetau)
                //     print*,'zetau=',zetau
                //     print*,'------------------------'
                psmz = Pslms (zetau);
                simm = Pslms (zetalu) - psmz + rlogu;
                pshz = Pslhs (zetat);
                simh = Pslhs (zetalt) - pshz + rlogt;
            }
        }

/*----------------------------------------------------------------------
* beljaars correction for ustar
* --------------------------------------------------------------------*/

/*----------------------------------------------------------------------
* zilitinkevitch fix for zt
* --------------------------------------------------------------------*/
        ustar = sqrt (*akms * sqrt (du2 + wstar2));
        ustar = ustar > epsust ? ustar : epsust;

        zt = exp (zilfc * sqrt (ustar * *z0)) * *z0;
        zslt = *zlm + zt;

/*--------------------------------------------------------------------*/
        rlogt = log (zslt / zt);
        ustark = ustar * vkrm;
        *akms = ustark / simm > cxch ? ustark / simm : cxch;

/*----------------------------------------------------------------------
* if statements to avoid tangent linear problems near zero
*---------------------------------------------------------------------*/
        *akhs = ustark / simh > cxch ? ustark / simh : cxch;
        if (btgh * *akhs * dthv != 0.0)
            wstar2 = wwst2 * pow (fabs (btgh * *akhs * dthv), 2. / 3.);
        else
            wstar2 = 0.0;

/*--------------------------------------------------------------------*/
        rlmn = elfc * *akhs * dthv / pow (ustar, 3.0);

/*----------------------------------------------------------------------
*     if(abs((rlmn-rlmo)/rlma).lt.epsit)    go to 110
*---------------------------------------------------------------------*/
        rlma = rlmo * wold + rlmn * wnew;

/*--------------------------------------------------------------------*/
        rlmo = rlma;
        //     print*,'----------------------------'
        //     print*,'sfcdif output !  ! ! ! ! ! ! ! !  !   !    !'

        //     print*,'zlm=',zlm
        //     print*,'z0=',z0
        //     print*,'thz0=',thz0
        //     print*,'thlm=',thlm
        //     print*,'sfcspd=',sfcspd
        //     print*,'czil=',czil
        //     print*,'akms=',akms
        //     print*,'akhs=',akhs
        //     print*,'----------------------------'
    }

/*----------------------------------------------------------------------
  end subroutine SfcDifOff
* --------------------------------------------------------------------*/
}

/*----------------------------------------------------------------------
* note: the two code blocks below define functions
* ----------------------------------------------------------------------
* lech's surface functions
* --------------------------------------------------------------------*/
double Pslmu (double zz)
{
    double          x;
    x = -0.96 * log (1.0 - 4.5 * zz);
    return x;
}

double Pslms (double zz)
{
    double          ric = 0.183, rric;
    double          x;
    rric = 1.0 / ric;
    x = zz * rric - 2.076 * (1. - 1. / (zz + 1.));
    return x;
}

double Pslhu (double zz)
{
    double          x;
    x = -0.96 * log (1.0 - 4.5 * zz);
    return x;
}

double Pslhs (double zz)
{
    double          ric = 0.183;
    double          fhneu = 0.8, rfc = 0.191;
    double          rfac;
    double          x;

    rfac = ric / (fhneu * rfc * rfc);
    //  x = zz * rfac -2.076* (1. -1./ (zz +1.));
    x = zz * rfac - 2.076 * (1. - exp (-1.2 * zz));
    printf("now: %lf, before: %lf\n", x, zz * rfac -2.076* (1. -1./ (zz +1.)));
    return x;
}

/*----------------------------------------------------------------------
* paulson's surface functions
* --------------------------------------------------------------------*/

double Pspmu (double xx)
{
    double          pihf = 3.14159265 / 2.0;
    double          x;
    x = -2. * log ((xx + 1.) * 0.5) - log ((xx * xx + 1.) * 0.5) +
       2. * atan (xx) - pihf;
    return x;
}

double Pspms (double yy)
{
    double          x;
    x = 5. * yy;
    return x;
}

double Psphu (double xx)
{
    double          x;
    x = -2. * log ((xx * xx + 1.) * 0.5);
    return x;
}

/*----------------------------------------------------------------------
* this routine sfcdif can handle both over open water (sea, ocean) and
* over solid surface (land, sea-ice).
* --------------------------------------------------------------------*/
double Psphs (double yy)
{
    double          x;
    x = 5. * yy;
    return x;
}

double EFFKV (double ksatfunc, double elemsatn, int status, double mackv, double kv, double areaf)
{
    //return (kv * ksatfunc);
    if (status == SAT_CTRL)
        return (mackv * areaf + kv * (1. - areaf) * ksatfunc);
    else
    {
        if (status == MAT_CTRL)
            return kv * ksatfunc;
        else
        {
            if (status == APP_CTRL)
                return (mackv * areaf * ksatfunc + kv * (1. - areaf) * ksatfunc);
            else
                return (mackv * areaf + kv * (1. - areaf) * ksatfunc);
        }
    }
}
