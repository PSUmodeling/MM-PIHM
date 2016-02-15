#include "pihm.h"

void Noah (int t, pihm_struct pihm)
{
    int             i, j;
    double          zenith, azimuth;
    double          sdir, sdif;
    elem_struct    *elem;

    /*
     * Calculate Sun position for topographic solar radiation
     */
    SunPos (t, pihm->latitude, pihm->longitude, pihm->elevation,
        pihm->noahtbl.tbot, &zenith, &azimuth);

    /*
     * Calculate avearge water fluxes from PIHM
     */
    AvgFlux (pihm->elem, pihm->numele, AVG);

    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        /* Read forcing */
        elem->ps.sfcspd = *elem->forc.meteo[SFCSPD_TS];
        elem->es.sfctmp = *elem->forc.meteo[SFCTMP_TS];
        elem->ps.rh = *elem->forc.meteo[RH_TS];
        elem->ps.sfcprs = *elem->forc.meteo[PRES_TS];
        elem->ef.longwave = *elem->forc.meteo[LONGWAVE_TS];
        elem->wf.prcp = *elem->forc.meteo[PRCP_TS] * 0.001;
        /* Calculate solar radiation */
        if (pihm->ctrl.rad_mode > 0)
        {
            sdir = *elem->forc.rad[SOLAR_DIR_TS];
            sdif = *elem->forc.rad[SOLAR_DIF_TS];

            elem->ef.soldn =
                TopoRadn (sdir, sdif, zenith, azimuth, elem->topo.slope,
                elem->topo.aspect, elem->topo.h_phi, elem->topo.svf);
        }
        else
        {
            elem->ef.soldn = *elem->forc.meteo[SOLAR_TS];
        }

        CalHum (&elem->ps, &elem->es);

        elem->ps.ffrozp = FrozRain (elem->wf.prcp, elem->es.sfctmp);

        elem->ps.alb = BADVAL;

        if (elem->forc.lai_type > 0)
        {
            elem->ps.xlai = *elem->forc.lai;
        }
        else
        {
            elem->ps.xlai = MonthlyLAI (t, elem->lc.type);
        }

        elem->ws.cmcmax = elem->lc.cmcfactr * elem->ps.xlai;

        if (elem->ps.q1 == BADVAL)
        {
            elem->ps.q1 = elem->ps.q2;
        }

        elem->ef.solnet = elem->ef.soldn * (1.0 - elem->ps.albedo);
        elem->ef.lwdn = elem->ef.longwave * elem->ps.emissi;

        elem->ps.nwtbl = FindLayer (elem->ps.sldpth, elem->ps.nsoil,
            elem->soil.depth - elem->ws.gw);
        elem->ps.nwtbl = (elem->ps.nwtbl > elem->ps.nsoil) ?
            elem->ps.nsoil : elem->ps.nwtbl;

        for (j = 0; j < elem->ps.nsoil; j++)
        {
            elem->ws.smc[j] = (elem->ws.smc[j] > elem->soil.smcmin + 0.02) ?
                elem->ws.smc[j] : elem->soil.smcmin + 0.02;
            elem->ws.smc[j] = (elem->ws.smc[j] < elem->soil.smcmax) ?
                elem->ws.smc[j] : elem->soil.smcmax;
            elem->ws.sh2o[j] = (elem->ws.sh2o[j] < elem->ws.smc[j]) ?
                elem->ws.sh2o[j] : elem->ws.smc[j];
        }

        /*
         * Run Noah LSM
         */
        SFlx (&elem->ws, &elem->wf, &elem->avgwf, &elem->es, &elem->ef,
            &elem->ps, &elem->lc, &elem->soil, pihm->ctrl.etstep);

        /*
         * Transfer Noah variables to PIHM
         */
        elem->wf.netprcp = elem->wf.pcpdrp;
        /* ET: convert from w m-2 to m s-1 */
        elem->wf.ec = elem->ef.ec / LVH2O / 1000.0;
        elem->wf.ett = elem->ef.ett / LVH2O / 1000.0;
        elem->wf.edir = elem->ef.edir / LVH2O / 1000.0;

        ZeroWaterFlux (&elem->avgwf);
    }
}

void SFlx (ws_struct *ws, wf_struct *wf, const wf_struct *avgwf,
    es_struct *es, ef_struct *ef, ps_struct *ps, lc_struct *lc,
    soil_struct *soil, int etstep)
{
    /*
     * subroutine SFlx - unified noahlsm version 1.0 july 2007
     *
     * Sub-driver for "Noah LSM" family of physics subroutines for a
     * soil/veg/snowpack land-surface model to update soil moisture, soil ice,
     * soil temperature, skin temperature, snowpack water content, snowdepth,
     * and all terms of the surface energy balance and surface water balance
     * (excluding input atmospheric forcings of downward radiation and precip)
     */
    int             frzgra, snowng;
    const int       IZ0TLND = 0;
    int             kz;
    double          zsoil[MAXLYR];
    double          df1;
    double          df1h, df1a;
    double          dsoil;
    double          dtot;
    double          frcsno, frcsoi;
    double          t1v;
    double          th2v;
    double          t2v;
    double          t24;
    double          interp_fraction;
    double          sn_new;
    double          dt;
    double          prcpf;
    double          soilwm;
    double          soilww;
    double          smav[MAXLYR];
    int             k;

    dt = (double)etstep;

    /*
     * Initialization
     */
    wf->snomlt = 0.0;

    /*
     * Calculate depth (negative) below ground from top skin sfc to bottom of
     * each soil layer.  note:  sign of zsoil is negative (denoting below
     * ground)
     */

    zsoil[0] = - ps->sldpth[0];
    for (kz = 1; kz < ps->nsoil; kz++)
    {
        zsoil[kz] = - ps->sldpth[kz] + zsoil[kz - 1];
    }

    /*
     * Next is crucial call to set the land-surface parameters, including
     * soil-type and veg-type dependent parameters.
     */
    //RedPrm(grid, lsm, zsoil);           /* ys: RedPrm is now called in driver */

    wf->pcpdrp = 0.0;

    /*
     * Urban 
     */
    if (lc->isurban)
    {
        lc->shdfac = 0.05;
        lc->rsmin = 400.0;
        soil->smcmax = 0.45;
        soil->smcmin = 0.0;
        soil->smcref = 0.42;
        soil->smcwlt = 0.40;
        soil->smcdry = 0.40;
    }

    /* YS: Flux-PIHM uses LAI as a forcing variable.
     * Vegetation fraction is calculated from LAI following Noah-MP */
    if (ps->xlai >= lc->laimax)
    {
        ps->embrd = lc->emissmax;
        ps->alb = lc->albedomin;
        ps->z0brd = lc->z0max;
    }
    else if (ps->xlai <= lc->laimin)
    {
        ps->embrd = lc->emissmin;
        ps->alb = lc->albedomax;
        ps->z0brd = lc->z0min;
    }
    else
    {
        if (lc->laimax > lc->laimin)
        {
            interp_fraction = (ps->xlai - lc->laimin) / (lc->laimax - lc->laimin);

            /* Bound interp_fraction between 0 and 1 */
            interp_fraction = (interp_fraction < 1.0) ? interp_fraction : 1.0;
            interp_fraction = (interp_fraction > 0.0) ? interp_fraction : 0.0;

            /* Scale emissivity and LAI between emissmin and emissmax by
             * interp_fraction */
            ps->embrd =
                ((1.0 - interp_fraction) * lc->emissmin) +
                (interp_fraction * lc->emissmax);
            ps->alb =
                ((1.0 - interp_fraction) * lc->albedomax) +
                (interp_fraction * lc->albedomin);
            ps->z0brd =
                ((1.0 - interp_fraction) * lc->z0min) +
                (interp_fraction * lc->z0max);
        }
        else
        {
            ps->embrd = 0.5 * lc->emissmin + 0.5 * lc->emissmax;
            ps->alb = 0.5 * lc->albedomin + 0.5 * lc->albedomax;
            ps->z0brd = 0.5 * lc->z0min + 0.5 * lc->z0max;
        }
    }

    //lc->shdfac = 1.0 - exp (-0.52 * (ps->xlai));
    lc->shdfac = 1.0 - exp (-0.75 * (ps->xlai));

    /*
     * Initialize precipitation logicals. 
     */
    snowng = 0;
    frzgra = 0;

    /* If input snowpack is nonzero, then compute snow density "sndens" and
     * snow thermal conductivity "SnCond" (note that CSnow is a function
     * subroutine */
    if (ws->sneqv <= 1.0e-7)       /* Safer if KMH (2008/03/25) */
    {
        ws->sneqv = 0.0;
        ps->sndens = 0.0;
        ps->snowh = 0.0;
        ps->sncond = 1.0;
    }
    else
    {
        ps->sndens = ws->sneqv / ps->snowh;
        if (ps->sndens > 1.0)
        {
            printf ("ERROR: "
                "Physical snow depth is less than snow water equiv.\n");
            PihmExit (1);
        }
        ps->sncond = CSnow (ps->sndens);
    }

    /* Determine if it's precipitating and what kind of precip it is.
     * If it's prcping and the air temp is colder than 0 C, it's snowing!
     * If it's prcping and the air temp is warmer than 0 C, but the grnd
     * temp is colder than 0 C, freezing rain is presumed to be falling. */
    if (wf->prcp > 0.0)
    {
        /* Snow defined when fraction of frozen precip (ffrozp) > 0.5, passed
         * in from model microphysics.  */
        if (ps->ffrozp > 0.5)
        {
            snowng = 1;
        }
        else
        {
            if (es->t1 <= TFREEZ)
            {
                frzgra = 1;
            }
        }
    }

    /* If either prcp flag is set, determine new snowfall (converting prcp
     * rate from kg m-2 s-1 to a liquid equiv snow depth in meters) and add
     * it to the existing snowpack.
     * Note that since all precip is added to snowpack, no precip infiltrates
     * into the soil so that prcp1 is set to zero. */
    if (snowng || frzgra)
    {
        //sn_new = wf->prcp * dt * 0.001;
        sn_new = wf->prcp * dt;
        ws->sneqv += sn_new;
        prcpf = 0.0;

        /* Update snow density based on new snowfall, using old and new snow.
         * Update snow thermal conductivity */
        SnowNew (es, sn_new, ps);
        ps->sncond = CSnow (ps->sndens);
    }
    else
    {
        /* Precip is liquid (rain), hence save in the precip variable that later
         * can wholely or partially infiltrate the soil (along with any canopy
         * "drip" added to this later) */
        prcpf = wf->prcp;
    }

    /* 
     * Determine snowcover and albedo over land.
     */
    if (ws->sneqv == 0.0)
    {
        /* If snow depth = 0, set snow fraction = 0, albedo = snow free
         * albedo. */
        ps->sncovr = 0.0;
        ps->albedo = ps->alb;
        ps->emissi = ps->embrd;
    }
    else
    {
        /* Determine snow fractional coverage.
         * Determine surface albedo modification due to snowdepth state. */
        ps->sncovr = SnFrac (ws->sneqv, lc->snup, ps->salp, ps->snowh);

        ps->sncovr = (ps->sncovr < 0.98) ? ps->sncovr : 0.98;

        AlCalc (ps, dt, snowng);
    }

    /*
     * Next calculate the subsurface heat flux, which first requires
     * calculation of the thermal diffusivity. Treatment of the latter
     * follows that on Pages 148-149 from "Heat transfer in cold climates", by
     * V. J. Lunardini (published in 1981 by van Nostrand Reinhold Co.) i.e.
     * treatment of two contiguous "plane parallel" mediums (namely here the
     * first soil layer and the snowpack layer, if any). This diffusivity
     * treatment behaves well for both zero and nonzero snowpack, including
     * the limit of very thin snowpack.  this treatment also eliminates the
     * need to impose an arbitrary upper bound on subsurface heat flux when
     * the snowpack becomes extremely thin.
     *
     * First calculate thermal diffusivity of top soil layer, using both the
     * frozen and liquid soil moisture, following the soil thermal diffusivity
     * function of Peters-Lidard et al. (1998, JAS, Vol 55, 1209-1224), which
     * requires the specifying the quartz content of the given soil class (see
     * routine RedPrm)
     *
     * Next add subsurface heat flux reduction effect from the overlying green
     * canopy, adapted from Section 2.1.2 of Peters-Lidard et al. (1997, JGR,
     * Vol 102(D4))
     */
    df1 = TDfCnd (ws->smc[0], soil->quartz, soil->smcmax, soil->smcmin,
        ws->sh2o[0]);

    /* Urban */
    if (lc->isurban)
    {
        df1 = 3.24;
    }

    df1 *= exp (ps->sbeta * lc->shdfac);

    /* KMH 09/03/2006
     * KMH 03/25/2008  Change sncovr threshold to 0.97 */
    if (ps->sncovr > 0.97)
    {
        df1 = ps->sncond;
    }

    /* Finally "plane parallel" snowpack effect following V. J. Linardini
     * reference cited above. Note that dtot is combined depth of snowdepth
     * and thickness of first soil layer */
    dsoil = - (0.5 * zsoil[0]);
    if (ws->sneqv == 0.0)
    {
        ef->ssoil = df1 * (es->t1 - es->stc[0]) / dsoil;
    }
    else
    {
        dtot = ps->snowh + dsoil;
        frcsno = ps->snowh / dtot;

        /* 1. harmonic mean (series flow) */
        //df1 = (sncond * df1) / (frcsoi * sncond + frcsno * df1);
        frcsoi = dsoil / dtot;

        /* 2. arithmetic mean (parallel flow) */
        //df1 = frcsno * sncond + frcsoi * df1;
        df1h =
            (ps->sncond * df1) / (frcsoi * ps->sncond + frcsno * df1);

        /* 3. geometric mean (intermediate between harmonic and arithmetic
         * mean) */
        //df1 = pow (sncond, frcsno) * pow(df1, frcsoi);
        /* weigh df by snow fraction */
        //df1 = df1h * sncovr + df1a * (1.0-sncovr);
        //df1 = df1h * sncovr + df1 * (1.0-sncovr);
        df1a = frcsno * ps->sncond + frcsoi * df1;

        /* Calculate subsurface heat flux, ssoil, from final thermal
         * diffusivity of surface mediums, df1 above, and skin temperature and
         * top mid-layer soil temperature */
        df1 = df1a * ps->sncovr + df1 * (1.0 - ps->sncovr);
        ef->ssoil = df1 * (es->t1 - es->stc[0]) / dtot;
    }

    /*
     * Determine surface roughness over snowpack using snow condition from the
     * previous timestep.
     */
    if (ps->sncovr > 0.0)
    {
        ps->z0 = Snowz0 (ps->sncovr, ps->z0brd, ps->snowh);
    }
    else
    {
        ps->z0 = ps->z0brd;
    }

    /*
     * Next call function SfcDif to calculate the sfc exchange coef (ch) for
     * heat and moisture.
     *
     * Note !!!
     * Do not call SfcDif until after above call to RedPrm, in case
     * alternative values of roughness length (z0) and Zilintinkevich coef
     * (czil) are set there via namelist i/o.
     *
     * Note !!!
     * Function SfcDif returns a ch that represents the wind spd times the
     * "original" nondimensional "ch" typical in literature. Hence the ch
     * returned from SfcDif has units of m/s. The important companion
     * coefficient of ch, carried here as "rch", is the ch from sfcdif times
     * air density and parameter "CP". "rch" is computed in "Penman".
     * rch rather than ch is the coeff usually invoked later in eqns.
     *
     * Note !!!
     * SfcDif also returns the surface exchange coefficient for momentum, cm,
     * also known as the surface drage coefficient. Needed as a state variable
     * for iterative/implicit solution of ch in SfcDif
     */

     t1v = es->t1 * (1.0 + 0.61 * ps->q2);
     th2v = es->th2 * (1.0 + 0.61 * ps->q2);
     SfcDifOff (ps, lc, t1v, th2v, IZ0TLND);

    /*
     * Call Penman function to calculate potential evaporation (ETP), and
     * other partial products and sums save in common/rite for later
     * calculations.
     */

    /* Calculate total downward radiation (solar plus longwave) needed in
     * Penman ep subroutine that follows */
    //es->fdown = es->soldn * (1.0- ps->albedo) + ef->lwdn;
    ef->fdown = ef->solnet + ef->lwdn;

    /* Calc virtual temps and virtual potential temps needed by Penman. */
    t2v = es->sfctmp * (1.0 + 0.61 * ps->q2);
    
    Penman (wf, es, ef, ps, &t24, t2v, snowng, frzgra);

    /*
     * Call CanRes to calculate the canopy resistance and convert it into pc
     * if nonzero greenness fraction
     * Frozen ground extension: total soil water "smc" was replaced by
     * unfrozen soil water "sh2o" in call to CanRes below
     */
    if (lc->shdfac > 0.0)
    {
        CanRes (ws, es, ef, ps, zsoil, soil, lc);
    }
    else
    {
        ps->rc = 0.0;
    }

    /*
     * Now decide major pathway branch to take depending on whether snowpack
     * exists or not
     */
    wf->esnow = 0.0;

    if (ws->sneqv == 0.0)
    {
        NoPac (ws, wf, avgwf, es, ef, ps, lc, soil, zsoil, dt, t24);
        ps->eta_kinematic = wf->eta * 1000.0;
    }
    else
    {
        SnoPac (ws, wf, avgwf, es, ef, ps, lc, soil, snowng, zsoil, dt, t24,
            prcpf, df1);
        ps->eta_kinematic = (wf->esnow + wf->etns) * 1000.0;
    }

    /* Calculate effective mixing ratio at grnd level (skin) */

    ps->q1 = ps->q2 + ps->eta_kinematic * CP / ps->rch;

    /* Determine sensible heat (H) in energy units (W m-2) */
    ef->sheat =
        - (ps->ch * CP * ps->sfcprs) / (RD * t2v) * (es->th2 - es->t1);

    /* Convert evap terms from rate (m s-1) to energy units (w m-2) */
    ef->edir = wf->edir * 1000.0 * LVH2O;
    ef->ec = wf->ec * 1000.0 * LVH2O;
    for (k = 0; k < ps->nsoil; k++)
    {
        ef->et[k] = wf->et[k] * 1000.0 * LVH2O;
    }
    ef->ett = wf->ett * 1000.0 * LVH2O;
    ef->esnow = wf->esnow * 1000.0 * LSUBS;
    ef->etp = wf->etp * 1000.0 *
        ((1.0 - ps->sncovr) * LVH2O + ps->sncovr * LSUBS);
    if (ef->etp > 0.0)
    {
        ef->eta = ef->edir + ef->ec + ef->ett + ef->esnow;
    }
    else
    {
        ef->eta = ef->etp;
    }

    /* Determine beta (ratio of actual to potential evap) */
    if (ef->etp == 0.0)
    {
        ps->beta = 0.0;
    }
    else
    {
        ps->beta = ef->eta / ef->etp;
    }

    /* Convert the sign of soil heat flux so that:
     *   ssoil>0: warm the surface  (night time)
     *   ssoil<0: cool the surface  (day time) */
    ef->ssoil *= -1.0;

    ///* Definitions of soilm and soilw have been changed in Flux-PIHM for
    // * coupling purpose */
    //ws->soilm = -1.0 * ws->sh2o[0] * zsoil[0];
    //for (k = 1; k < ps->nsoil; k++)
    //{
    //    ws->soilm += ws->sh2o[k] * (zsoil[k - 1] - zsoil[k]);
    //}

    //ws->soilw = -1.0 * ws->sh2o[0] * zsoil[0];
    //for (k = 1; k < lc->nroot; k++)
    //{
    //    ws->soilw += ws->sh2o[k] * (zsoil[k - 1] - zsoil[k]);
    //}
    //ws->soilw /= - zsoil[lc->nroot - 1];

    ws->soilm = -1.0 * ws->smc[0] * zsoil[0];
    for (k = 1; k < ps->nsoil; k++)
    {
        ws->soilm += ws->smc[k] * (zsoil[k - 1] - zsoil[k]);
    }

    soilwm = -1.0 * (soil->smcmax - soil->smcwlt) * zsoil[0];
    soilww = -1.0 * (ws->smc[0] - soil->smcwlt) * zsoil[0];

    for (k = 0; k < ps->nsoil; k++)
    {
        smav[k] = (ws->smc[k] - soil->smcwlt) / (soil->smcmax - soil->smcwlt);
    }

    if (lc->nroot > 0)
    {
        for (k = 1; k < lc->nroot; k++)
        {
            soilwm += (soil->smcmax - soil->smcwlt) * (zsoil[k - 1] - zsoil[k]);
            soilww += (ws->smc[k] - soil->smcwlt) * (zsoil[k - 1] - zsoil[k]);
        }
    }
    if (soilwm < 1.0e-6)
    {
        soilwm = 0.0;
        ws->soilw = 0.0;
        ws->soilm = 0.0;
    }
    else
    {
        ws->soilw = soilww / soilwm;
    }
}

void AlCalc (ps_struct *ps, double dt, int snowng)
{
    /*
     * Calculate albedo including snow effect (0 -> 1)
     */
    /* snoalb is argument representing maximum albedo over deep snow, as
     * passed into SFlx, and adapted from the satellite-based maximum snow
     * albedo fields provided by D. Robinson and G. Kukla (1985, JCAM, Vol 24,
     * 402-411) */
    double          snoalb2;
    double          snoalb1;
    const double    SNACCA = 0.94;
    const double    SNACCB = 0.58;

    /* Turn of vegetation effect */
    //ps->albedo = ps->alb + (1.0 - (lc->shdfac - lc->shdmin)) * ps->sncovr *
    //  (ps->snoalb - ps->alb);
    //ps->albedo = (1.0 - ps->sncovr) * ps->alb +
    //  ps->sncovr * ps->snoalb;    /* this is equivalent to below */

    ps->albedo = ps->alb + ps->sncovr * (ps->snoalb - ps->alb);
    ps->emissi = ps->embrd + ps->sncovr * (EMISSI_S - ps->embrd);

    /* Formulation by livneh
     * snoalb is considered as the maximum snow albedo for new snow, at a
     * value of 85%. Snow albedo curve defaults are from Bras P.263. should
     * not be changed except for serious problems with snow melt.
     * To implement accumulating parameters, SNACCA and SNACCB, assert that it
     * is indeed accumulation season. i.e. that snow surface temp is below
     * zero and the date falls between october and february */
    snoalb1 = ps->snoalb + ps->lvcoef * (0.85 - ps->snoalb);
    snoalb2 = snoalb1;

    /* Initial lstsnw */
    if (snowng)
    {
        ps->snotime1 = 0.0;
    }
    else
    {
        ps->snotime1 += dt;
        snoalb2 = snoalb1 * pow (SNACCA, pow (ps->snotime1 / 86400.0, SNACCB));
    }

    snoalb2 = (snoalb2 > ps->alb) ? snoalb2 : ps->alb;
    ps->albedo = ps->alb + ps->sncovr * (snoalb2 - ps->alb);
    ps->albedo = (ps->albedo > snoalb2) ? snoalb2 : ps->albedo;
}

void CanRes (ws_struct *ws, es_struct *es, ef_struct *ef, ps_struct *ps,
    const double *zsoil, const soil_struct *soil, const lc_struct *lc)
{
    /*
     * Function CanRes
     *
     * Calculate canopy resistance which depends on incoming solar radiation,
     * air temperature, atmospheric water vapor pressure deficit at the
     * lowest model level, and soil moisture (preferably unfrozen soil
     * moisture rather than total)
     *
     * Source:  Jarvis (1976), Noilhan and Planton (1989, MWR), Jacquemin and
     * Noilhan (1990, BLM)
     * See also:  Chen et al. (1996, JGR, Vol 101(D3), 7251-7268), Eqns 12-14
     * and Table 2 of Sec. 3.1.2
     */
    int             k;
    double          delta, ff, gx, rr;
    double          part[MAXLYR];
    const double    SLV = 2.501000e6;

    /* Initialize canopy resistance multiplier terms. */
    ps->rcs = 0.0;
    ps->rct = 0.0;
    ps->rcq = 0.0;
    ps->rcsoil = 0.0;

    ps->rc = 0.0;

    /* Contribution due to incoming solar radiation */
    ff = 0.55 * 2.0 * ef->soldn / (lc->rgl * ps->xlai);
    ps->rcs = (ff + lc->rsmin / lc->rsmax) / (1.0 + ff);
    ps->rcs = (ps->rcs > 0.0001) ? ps->rcs : 0.0001;

    /* Contribution due to air temperature at first model level above ground
     * rct expression from Noilhan and Planton (1989, MWR). */
    ps->rct = 1.0 - 0.0016 * pow (lc->topt - es->sfctmp, 2.0);
    ps->rct = (ps->rct > 0.0001) ? ps->rct : 0.0001;

    /* Contribution due to vapor pressure deficit at first model level.
     * rcq expression from ssib */
    ps->rcq = 1.0 / (1.0 + lc->hs * (ps->q2sat - ps->q2));
    ps->rcq = (ps->rcq > 0.01) ? ps->rcq : 0.01;

    /* Contribution due to soil moisture availability.
     * Determine contribution from each soil layer, then add them up. */
    gx = (ws->sh2o[0] - soil->smcwlt) / (soil->smcref - soil->smcwlt);
    gx = (gx > 1.0) ? 1.0 : gx;
    gx = (gx < 0.0) ? 0.0 : gx;

    /* Use root distribution as weighting factor */
    //part[0]1 = rtdis[0] * gx;
    /* Use soil depth as weighting factor */
    part[0] = (zsoil[0] / zsoil[lc->nroot - 1]) * gx;
    for (k = 1; k < lc->nroot; k++)
    {
        gx = (ws->sh2o[k] - soil->smcwlt) / (soil->smcref - soil->smcwlt);
        gx = (gx > 1.0) ? 1.0 : gx;
        gx = (gx < 0.0) ? 0.0 : gx;

        /* Use root distribution as weighting factor */
        //part[k] = rtdis[k] * gx;
        /* Use soil depth as weighting factor */
        part[k] = ((zsoil[k] - zsoil[k - 1]) / zsoil[lc->nroot - 1]) * gx;
    }

    for (k = 0; k < lc->nroot; k++)
    {
        ps->rcsoil += part[k];
    }
    ps->rcsoil = (ps->rcsoil > 0.0001) ? ps->rcsoil : 0.0001;

    /* Determine canopy resistance due to all factors.  convert canopy
     * resistance (rc) to plant coefficient (pc) to be used with potential
     * evap in determining actual evap. pc is determined by:
     *   pc * linerized Penman potential evap =
     *   Penman-monteith actual evaporation (containing rc term). */
    ps->rc = lc->rsmin / (ps->xlai * ps->rcs * ps->rct * ps->rcq * ps->rcsoil);
    //rr = (4.0 * SIGMA * RD / CP) * pow(es->sfctmp, 4.0) /
    //  (ps->sfcprs * ps->ch) + 1.0;
    rr = (4.0 * ps->emissi * SIGMA * RD / CP) *
        pow (es->sfctmp, 4.0) / (ps->sfcprs * ps->ch) + 1.0;

    delta = (SLV / CP) * ps->dqsdt2;

    ps->pc = (rr + delta) / (rr * (1.0 + ps->rc * ps->ch) + delta);
}

double CSnow (double dsnow)
{
    /*
     * Function CSnow
     *
     * Calculate snow thermal conductivity
     */
    double          c;
    double          sncond;
    const double    UNIT = 0.11631;

    /* sncond in units of cal/(cm*hr*c), returned in W/(m*C)
     * csnow in units of cal/(cm*hr*c), returned in W/(m*C)
     * Basic version is Dyachkova equation (1960), for range 0.1-0.4 */

    c = 0.328 * pow (10, 2.25 * dsnow);

    /* De Vaux equation (1933), in range 0.1-0.6 */
    //sncond = 0.0293 * (1.0 + 100.0 * dsnow * dsnow);
    //csnow = 0.0293 * (1.0 + 100.0 * dsnow * dsnow);

    /* E. Andersen from Flerchinger */
    //sncond = 0.021 + 2.51 * dsnow * dsnow;
    //csnow = 0.021 + 2.51 * dsnow * dsnow;

    //sncond = UNIT * c;
    
    /* Double snow thermal conductivity */
    sncond = 2.0 * UNIT * c;

    return (sncond);
}

void DEvap (const ws_struct *ws, wf_struct *wf, const ps_struct *ps,
    const lc_struct *lc, const soil_struct *soil)
{
    /*
     * Function DEvap
     *
     * Calculate direct soil evaporation
     */
    double          fx, sratio;

    /* Direct evap a function of relative soil moisture availability, linear
     * when fxexp = 1.
     * fx > 1 represents demand control
     * fx < 1 represents flux control */

    sratio = (ws->sh2o[0] - soil->smcdry) / (soil->smcmax - soil->smcdry);
    if (sratio > 0.0)
    {
        fx = pow (sratio, ps->fxexp);
        fx = (fx > 1.0) ? 1.0 : fx;
        fx = (fx < 0.0) ? 0.0 : fx;
    }
    else
    {
        fx = 0.0;
    }

    /* Allow for the direct-evap-reducing effect of shade */
    wf->edir = fx * (1.0 - lc->shdfac) * wf->etp;
}

void Evapo (ws_struct *ws, wf_struct *wf, ps_struct *ps, const lc_struct *lc,
    const soil_struct *soil, const double *zsoil, double dt)
{
    /*
     * Function Evapo
     *
     * Calculate soil moisture flux. The soil moisture content (smc - a per
     * unit volume measurement) is a dependent variable that is updated with
     * prognostic eqns. The canopy moisture content (cmc) is also updated.
     * Frozen ground version: new states added: sh2o, and frozen ground
     * correction factor, frzfact and parameter slope.
     */
    int             k;
    double          cmc2ms;

    /* Executable code begins here if the potential evapotranspiration is
     * greater than zero. */
    wf->edir = 0.0;
    wf->ec = 0.0;
    wf->ett = 0.0;
    for (k = 0; k < ps->nsoil; k++)
    {
        wf->et[k] = 0.0;
    }

    if (wf->etp > 0.0)
    {
        if (lc->shdfac < 1.0)
        {
            /* Retrieve direct evaporation from soil surface. Call this
             * function only if veg cover not complete.
             * Frozen ground version:  sh2o states replace smc states. */
            DEvap (ws, wf, ps, lc, soil);
        }

        if (lc->shdfac > 0.0)
        {
            /* Initialize plant total transpiration, retrieve plant transpiration,
             * and accumulate it for all soil layers. */
            Transp (ws, wf, ps, lc, soil, zsoil);
            for (k = 0; k < ps->nsoil; k++)
            {
                wf->ett += wf->et[k];
            }

            /* Calculate canopy evaporation.
             * If statements to avoid tangent linear problems near
             * cmc = 0.0. */
            if (ws->cmc > 0.0)
            {
                wf->ec =
                    lc->shdfac * pow ((ws->cmc / ws->cmcmax > 1.0) ?
                        1.0 : ws->cmc / ws->cmcmax, lc->cfactr) * wf->etp;
            }
            else
            {
                wf->ec = 0.0;
            }

            /* ec should be limited by the total amount of available water on
             * the canopy.  F.Chen, 18-oct-1994 */
            cmc2ms = ws->cmc / dt;
            wf->ec = (cmc2ms < wf->ec) ? cmc2ms : wf->ec;
        }
    }

    /* Total up evap and transp types to obtain actual evapotransp */ 
    wf->etns = wf->edir + wf->ett + wf->ec;
}
//
//void Fac2Mit (double *smcmax, double *flimit)
//{
//    *flimit = 0.90;
//
//    if (*smcmax == 0.395)
//        *flimit = 0.59;
//    else if ((*smcmax == 0.434) || (*smcmax == 0.404))
//        *flimit = 0.85;
//    else if ((*smcmax == 0.465) || (*smcmax == 0.406))
//        *flimit = 0.86;
//    else if ((*smcmax == 0.476) || (*smcmax == 0.439))
//        *flimit = 0.74;
//    else if ((*smcmax == 0.200) || (*smcmax == 0.464))
//        *flimit = 0.80;
//
///*----------------------------------------------------------------------
//  end subroutine Fac2Mit
//* --------------------------------------------------------------------*/
//}
//
double FrH2O (double tkelv, double smc, double sh2o, const soil_struct *soil)
{
    /*
     * Function FrH2O
     *
     * Calculate amount of supercooled liquid soil water content if
     * temperature is below 273.15k (t0). Requires Newton-type iteration to
     * solve the nonlinear implicit equation given in Eqn 17 of Koren et al
     * (1999, JGR, Vol 104(D16), 19569-19585).
     *
     * New version (June 2001): much faster and more accurate Newton
     * iteration achieved by first taking log of eqn cited above -- less than
     * 4 (typically 1 or 2) iterations achieves convergence. Also, explicit
     * 1-step solution option for special case of parameter ck = 0, which
     * reduces the original implicit equation to a simpler explicit form,
     * known as the "Flerchinger eqn". Improved handling of solution in the
     * limit of freezing point temperature t0.
     */
    double          denom;
    double          df;
    double          dswl;
    double          fk;
    double          swl;
    double          swlk;
    double          satn;
    int             nlog, kcount;
    const double    CK = 8.0;
    const double    ERROR = 0.005;
    double          mx;
    double          freew;

    nlog = 0;

    /* If temperature not significantly below freezing (t0), sh2o = smc */
    kcount = 0;

    if (tkelv > (TFREEZ - 1.0e-3))
    {
        freew = smc;
    }
    else
    {
        if (CK != 0.0)
        {
            /* Option 1: iterated solution for nonzero CK in Koren et al, JGR,
             * 1999, Eqn 17
             * Initial guess for swl (frozen content) */
            swl = smc - sh2o;

            /* Keep within bounds. */
            if (swl > (smc - 0.02))
            {
                swl = smc - 0.02;
            }
            swl = (swl < 0.0) ? 0.0 : swl;

            /* Start of iterations */
            while (nlog < 10 && kcount == 0)
            {
                nlog++;

                satn = (smc - swl - soil->smcmin) / (soil->smcmax - soil->smcmin);
                mx = soil->beta / (1.0 - soil->beta);

                df = log(GRAV / soil->alpha / LSUBF) +
                    1.0 / soil->beta * log (pow (satn, mx) - 1.0) +
                    2.0 * log (1.0 + CK * swl) -
                    log (- (tkelv - TFREEZ) / tkelv);

                denom = 1.0 / (soil->beta - 1.0) / (soil->smcmax - soil->smcmin) *
                    pow (satn, mx - 1.0) / (pow (satn, mx) - 1.0) +
                    2.0 * CK / (1.0 + CK * swl);

                swlk = swl - df / denom;

                /* Bounds useful for mathematical solution. */
                swlk = (swlk > smc - 0.02) ? smc - 0.02 : swlk;
                swlk = (swlk < 0.0) ? 0.0 : swlk;

                /* Mathematical solution bounds applied. */
                dswl = fabs (swlk - swl);

                /* If more than 10 iterations, use explicit method (ck = 0
                 * approx.) when dswl less or eq. error, no more iterations
                 * required. */
                swl = swlk;
                if (dswl <= ERROR)
                {
                    kcount++;
                }
                /* End of iterations
                 * Bounds applied within do-block are valid for physical
                 * solution */
            }

            freew = smc - swl;
        }

        /* Option 2: explicit solution for flerchinger eq. i.e. ck = 0
         * in Koren et al., JGR, 1999, Eqn 17
         * Apply physical bounds to flerchinger solution */
        if (kcount == 0)
        {
            fk = pow (pow (- (tkelv - TFREEZ) / tkelv * soil->alpha * LSUBF / GRAV,
                soil->beta), 1.0 / mx) * (soil->smcmax - soil->smcmin) - soil->smcmin;
            fk = (fk < 0.02) ? 0.02 : fk;

            freew = (fk < smc) ? fk : smc;
        }
    }

    return (freew);
}

void HRT (ws_struct *ws, es_struct *es, ef_struct *ef, ps_struct *ps,
    const lc_struct *lc, const soil_struct *soil, double *rhsts,
    const double *zsoil, double yy, double zz1, double dt, double df1,
    double *ai, double *bi, double *ci)
{
    /*
     * Function HRT
     *
     * Calculate the right hand side of the time tendency term of the soil
     * thermal diffusion equation.  also to compute ( prepare ) the matrix
     * coefficients for the tri-diagonal matrix of the implicit time scheme.
     */
    int             itavg;
    int             k;

    double          ddz, ddz2;
    double          denom;
    double          df1k;
    double          dtsdz;
    double          dtsdz2;
    double          hcpct;
    double          sice;
    double          csoil_loc;
    double          df1n;
    double          qtot;
    double          tavg;
    double          tbk;
    double          tbk1;
    double          tsnsr;
    double          tsurf;
    const double    CH2O = 4.2e6;

    /* Urban */
    if (lc->isurban)
    {
        csoil_loc = 3.0e6;
    }
    else
    {
        csoil_loc = soil->csoil;
    }

    /* Initialize logical for soil layer temperature averaging. */
    itavg = 1;

    /* Begin section for top soil layer
     * Calc the heat capacity of the top soil layer */
    hcpct = ws->sh2o[0] * CH2O + (1.0 - soil->smcmax) * csoil_loc +
        (soil->smcmax - ws->smc[0]) * CP + (ws->smc[0] - ws->sh2o[0]) * CPICE;

    /* Calc the matrix coefficients ai, bi, and ci for the top layer */
    ddz = 1.0 / (- 0.5 * zsoil[1]);
    ai[0] = 0.0;
    ci[0] = (df1 * ddz) / (zsoil[0] * hcpct);

    /* Calculate the vertical soil temp gradient btwn the 1st and 2nd soil
     * layers. Then calculate the subsurface heat flux. use the temp gradient
     * and subsfc heat flux to calc "right-hand side tendency terms", or
     * "rhsts", for top soil layer. */
    bi[0] = -ci[0] + df1 / (0.5 * zsoil[0] * zsoil[0] * hcpct * zz1);
    dtsdz = (es->stc[0] - es->stc[1]) / (- 0.5 * zsoil[1]);
    ef->ssoil = df1 * (es->stc[0] - yy) / (0.5 * zsoil[0] * zz1);
    //rhsts[0] = (*df1 * dtsdz - ssoil) / (zsoil[0] * hcpct);
    denom = (zsoil[0] * hcpct);
    rhsts[0] = (df1 * dtsdz - ef->ssoil) / denom;

    /* Next capture the vertical difference of the heat flux at top and bottom
     * of first soil layer for use in heat flux constraint applied to
     * potential soil freezing/thawing in routine SnkSrc. */
    //  *qtot = ssoil - *df1*dtsdz;
    qtot = - 1.0 * rhsts[0] * denom;

    /* Calculate frozen water content in 1st soil layer. */
    sice = ws->smc[0] - ws->sh2o[0];

    /* If temperature averaging invoked (itavg = true; else skip):
     * Set temp "tsurf" at top of soil column (for use in freezing soil
     * physics later in function subroutine SnkSrc). If snowpack content is
     * zero, then tsurf expression below gives tsurf = skin temp. If
     * snowpack is nonzero (hence argument zz1 = 1), then tsurf expression
     * below yields soil column top temperature under snowpack. Then
     * calculate temperature at bottom interface of 1st soil layer for use
     * later in function subroutine SnkSrc */
    if (itavg)
    {
        tsurf = (yy + (zz1 - 1.0) * es->stc[0]) / zz1;

        /* If frozen water present or any of layer-1 mid-point or bounding
         * interface temperatures below freezing, then call SnkSrc to
         * compute heat source/sink (and change in frozen water content)
         * due to possible soil water phase change */
        tbk = TBnd (es->stc[0], es->stc[1], zsoil, ps->zbot, 0, ps->nsoil);

        if ((sice > 0.0) || (es->stc[0] < TFREEZ) ||
            (tsurf < TFREEZ) || (tbk < TFREEZ))
        {
            tavg = TmpAvg (tsurf, es->stc[0], tbk, zsoil, ps->nsoil, 0);
            SnkSrc (&tsnsr, tavg, ws->smc[0], &ws->sh2o[0], soil, zsoil, ps->nsoil,
                dt, 0, qtot);
            rhsts[0] = rhsts[0] - tsnsr / denom;
        }
    }
    else
    {
        if ((sice > 0.0) || (es->stc[0] < TFREEZ))
        {
            SnkSrc (&tsnsr, es->stc[0], ws->smc[0], &ws->sh2o[0], soil, zsoil, ps->nsoil,
                dt, 0, qtot);
            rhsts[0] = rhsts[0] - tsnsr / denom;
        }
        /* This ends section for top soil layer. */
    }

    /* Initialize ddz2 */
    ddz2 = 0.0;
    df1k = df1;

    /* Loop thru the remaining soil layers, repeating the above process
     * (except subsfc or "ground" heat flux not repeated in lower layers)
     * Calculate heat capacity for this soil layer. */
    for (k = 1; k < ps->nsoil; k++)
    {
        hcpct = ws->sh2o[k] * CH2O + (1.0 - soil->smcmax) * csoil_loc +
            (soil->smcmax - ws->smc[k]) * CP +
            (ws->smc[k] - ws->sh2o[k]) * CPICE;

        /* This section for layer 2 or greater, but not last layer.
         * Calculate thermal diffusivity for this layer. */
        if (k != ps->nsoil - 1)
        {
            /* Calc the vertical soil temp gradient thru this layer */
            df1n = TDfCnd (ws->smc[k], soil->quartz, soil->smcmax,
                soil->smcmin, ws->sh2o[k]);

            /* Urban */
            if (lc->isurban)
            {
                df1n = 3.24;
            }

            denom = 0.5 * (zsoil[k - 1] - zsoil[k + 1]);

            /* Calc the matrix coef, ci, after calc'ng its partial product */
            dtsdz2 = (es->stc[k] - es->stc[k + 1]) / denom;
            ddz2 = 2.0 / (zsoil[k - 1] - zsoil[k + 1]);

            /* If temperature averaging invoked (itavg = true; else skip)
             * Calculate temp at bottom of layer. */
            ci[k] = - df1n * ddz2 / ((zsoil[k - 1] - zsoil[k]) * hcpct);
            if (itavg)
            {
                tbk1 = TBnd (es->stc[k], es->stc[k + 1], zsoil, ps->zbot, k,
                    ps->nsoil);
            }
        }
        else
        {
            /* Special case of bottom soil layer
             * Calculate thermal diffusivity for bottom layer. */
            df1n = TDfCnd (ws->smc[k], soil->quartz, soil->smcmax,
                soil->smcmin, ws->sh2o[k]);

            /* Urban */
            if (lc->isurban)
            {
                df1n = 3.24;
            }

            /* Calc the vertical soil temp gradient thru bottom layer. */
            denom = 0.5 * (zsoil[k - 1] + zsoil[k]) - ps->zbot;
            dtsdz2 = (es->stc[k] - ps->tbot) / denom;

            /* Set matrix coef, ci to zero if bottom layer. */
            ci[k] = 0.0;

            /* If temperature averaging invoked (itavg = true; else skip)
             * Calculate temp at bottom of last layer. */
            if (itavg)
            {
                tbk1 = TBnd (es->stc[k], ps->tbot, zsoil, ps->zbot, k,
                    ps->nsoil);
            }

            /* This ends special loop for bottom layer. */
        }

        /* Calculate rhsts for this layer after calc'ng a partial product. */
        denom = (zsoil[k] - zsoil[k - 1]) * hcpct;
        rhsts[k] = (df1n * dtsdz2 - df1k * dtsdz) / denom;
        qtot = - 1.0 * denom * rhsts[k];

        sice = ws->smc[k] - ws->sh2o[k];

        if (itavg)
        {
            tavg = TmpAvg (tbk, es->stc[k], tbk1, zsoil, ps->nsoil, k);
            if ((sice > 0.0) || (es->stc[k] < TFREEZ) ||
                (tbk < TFREEZ) || (tbk1 < TFREEZ))
            {
                SnkSrc (&tsnsr, tavg, ws->smc[k], &ws->sh2o[k], soil, zsoil,
                    ps->nsoil, dt, k, qtot);
                rhsts[k] = rhsts[k] - tsnsr / denom;
            }
        }
        else
        {
            if ((sice > 0.0) || (es->stc[k] < TFREEZ))
            {
                SnkSrc (&tsnsr, es->stc[k], ws->smc[k], &ws->sh2o[k], soil,
                    zsoil, ps->nsoil, dt, k, qtot);
                rhsts[k] = rhsts[k] - tsnsr / denom;
            }
        }

        /* Calc matrix coefs, ai, and bi for this layer. */
        ai[k] = - df1k * ddz / ((zsoil[k - 1] - zsoil[k]) * hcpct);
        bi[k] = -(ai[k] + ci[k]);

        /* Reset values of df1, dtsdz, ddz, and tbk for loop to next soil
         * layer. */
        tbk = tbk1;
        df1k = df1n;
        dtsdz = dtsdz2;
        ddz = ddz2;
    }
}

void HStep (es_struct *es, double *rhsts, double dt, int nsoil, double *ai,
    double *bi, double *ci)
{
    /*
     * Subroutine HStep
     *
     * Calculate/update the soil temperature field.
     */
    int             k;
    double          rhstsin[MAXLYR];
    double          ciin[MAXLYR];

    /* Create finite difference values for use in Rosr12 routine */
    for (k = 0; k < nsoil; k++)
    {
        rhsts[k] *= dt;
        ai[k] *= dt;
        bi[k] *= dt;
        bi[k] += 1.0;
        ci[k] *= dt;
    }

    /* Copy values for input variables before call to Rosr12 */
    for (k = 0; k < nsoil; k++)
    {
        rhstsin[k] = rhsts[k];
    }
    for (k = 0; k < nsoil; k++)
    {
        ciin[k] = ci[k];
    }

    /* Solve the tri-diagonal matrix equation */
    Rosr12 (ci, ai, bi, ciin, rhstsin, rhsts, nsoil);

    /* Calc/update the soil temps using matrix solution */
    for (k = 0; k < nsoil; k++)
    {
        es->stc[k] += ci[k];
    }
}

void NoPac (ws_struct *ws, wf_struct *wf, const wf_struct *avgwf,
    es_struct *es, ef_struct *ef, ps_struct *ps, lc_struct *lc,
    soil_struct *soil, const double *zsoil, double dt, double t24)
{
    /*
     * Function NoPac
     *
     * Calculate soil moisture and heat flux values and update soil moisture
     * content and soil heat content values for the case when no snow pack is
     * present.
     */
    int             k;
    double          df1;
    double          yy;
    double          zz1;
    double          yynum;
    double          prcpf;

    //*prcp1 = *prcp * 0.001;
    //*etp1 = *etp * 0.001;
    prcpf = wf->prcp;
    wf->dew = 0.0;

    /* Initialize evap terms */
    //*edir = 0.;
    //*edir1 = 0.;
    //*ec1 = 0.;
    //*ec = 0.;
    wf->edir = 0.0;
    wf->ec = 0.0;

    for (k = 0; k < ps->nsoil; k++)
    {
        wf->et[k] = 0.0;
        //et1[k] = 0.;
    }

    wf->ett = 0.0;
    //*ett1 = 0.;

    if (wf->etp > 0.0)
    {
        Evapo (ws, wf, ps, lc, soil, zsoil, dt);

        wf->eta = wf->etns;

        SmFlx (ws, wf, avgwf, ps, lc, soil, zsoil, prcpf, dt);

        /* Convert modeled evapotranspiration from m s-1 to kg m-2 s-1. */
        //*eta = *eta1 * 1000.0;
    }
    else
    {
        /* If etp < 0, assume dew forms (transform etp into dew and
         * reinitialize etp to zero). */
        wf->dew = - wf->etp;

        /* Add dew amount to prcp */
        prcpf += wf->dew;
        SmFlx (ws, wf, avgwf, ps, lc, soil, zsoil, prcpf, dt);

        /* Convert modeled evapotranspiration from 'm s-1' to 'kg m-2 s-1'. */
        //*eta = *eta1 * 1000.0
    }

    /* Based on etp and e values, determine beta */
    if (wf->etp <= 0.0)
    {
        ps->beta = 0.0;
        wf->eta = wf->etp;
        if (wf->etp < 0.0)
        {
            ps->beta = 1.0;
        }
    }
    else
    {
        ps->beta = wf->eta / wf->etp;
    }

    /* Convert modeled evapotranspiration components 'm s-1' to 'kg m-2 s-1'. */
    //*edir = *edir1 * 1000.;
    //*ec = *ec1 * 1000.;
    //for (k = 0; k < *nsoil; k++)
    //    et[k] = et1[k] * 1000.;
    //*ett = *ett1 * 1000.;

    /* Get soil thermal diffuxivity/conductivity for top soil lyr, Calc.
     * adjusted top lyr soil temp and adjusted soil flux, then call ShFlx to
     * compute/update soil heat flux and soil temps. */
    df1 = TDfCnd (ws->smc[0], soil->quartz, soil->smcmax, soil->smcmin,
        ws->sh2o[0]);

    /* Urban */
    if (lc->isurban)
    {
        df1 = 3.24;
    }

    /* Vegetation greenness fraction reduction in subsurface heat flux via
     * reduction factor, which is convenient to apply here to thermal
     * diffusivity that is later used in HRT to compute sub sfc heat flux
     * (see additional comments on veg effect sub-sfc heat flx in function
     * SFlx) */
    df1 *= exp (ps->sbeta * lc->shdfac);

    /* Compute intermediate terms passed to routine HRT (via routine ShFlx
     * below) for use in computing subsurface heat flux in HRT */
    yynum = ef->fdown - ps->emissi * SIGMA * t24;
    yy = es->sfctmp +
        (yynum / ps->rch + es->th2 - es->sfctmp - ps->beta * ef->epsca) /
            ps->rr;

    zz1 = df1 / (- 0.5 * zsoil[0] * ps->rch * ps->rr) + 1.0;

    ShFlx (ws, es, ef, ps, lc, soil, dt, yy, zz1, zsoil, df1);

    /* Set flx1 and flx3 (snopack phase change heat fluxes) to zero since
     * they are not used here in SnoPac. flx2 (freezing rain heat flux) was
     * similarly initialized in the Penman routine. */
    ef->flx1 = CPH2O * wf->prcp * 1000.0 * (es->t1 - es->sfctmp);
    ef->flx3 = 0.0;
}

void Penman (wf_struct *wf, es_struct *es, ef_struct *ef, ps_struct *ps,
    double *t24, double t2v, int snowng, int frzgra)
{
    /*
     * Function Penman
     *
     * Calculate potential evaporation for the current point. Various partial
     * sums/products are also calculated and passed back to the calling
     * routine for later use.
     */
    double          a;
    double          delta;
    double          fnet;
    double          rad;
    double          rho;
    double          emissi;
    double          elcp1;
    double          lvs;
    const double    ELCP = 2.4888e+3;
    const double    LSUBC = 2.501000e+6;

    /* executable code begins here:
     * Prepare partial quantities for Penman equation. */
    emissi = ps->emissi;
    elcp1 = (1.0 - ps->sncovr) * ELCP + ps->sncovr * ELCP * LSUBS / LSUBC;
    lvs = (1.0 - ps->sncovr) * LSUBC + ps->sncovr * LSUBS;

    ef->flx2 = 0.0;
    //delta = ELCP * dqsdt2;
    delta = elcp1 * ps->dqsdt2;
    *t24 = es->sfctmp * es->sfctmp * es->sfctmp * es->sfctmp;
    //ps->rr = *t24 * 6.48e-8 / (ps->sfcprs * ps->ch) + 1.0;
    ps->rr = emissi * *t24 * 6.48e-8 / (ps->sfcprs * ps->ch) + 1.0;
    rho = ps->sfcprs / (RD * t2v);

    /* Adjust the partial sums / products with the latent heat effects caused
     * by falling precipitation. */
    ps->rch = rho * CP * ps->ch;
    if (!snowng)
    {
        if (wf->prcp > 0.0)
        {
            //ps->rr += CPH2O * wf->prcp / ps->rch;
            ps->rr += CPH2O * wf->prcp * 1000.0 / ps->rch;
        }
    }
    else
    {
        //ps->rr += CPICE * wf->prcp / ps->rch;
        ps->rr += CPICE * wf->prcp * 1000.0 / ps->rch;
    }

    /* Include the latent heat effects of frzng rain converting to ice on
     * impact in the calculation of flx2 and fnet. */
    //fnet = ef->fdown - SIGMA * t24- ef->ssoil;
    fnet = ef->fdown - emissi * SIGMA * *t24 - ef->ssoil;
    if (frzgra)
    {
        ef->flx2 = - LSUBF * wf->prcp * 1000.0;
        fnet -= ef->flx2;

        /* Finish Penman equation calculations */
    }

    rad = fnet / ps->rch + es->th2 - es->sfctmp;
    //a = ELCP * (ps->q2sat - ps->q2);
    a = elcp1 * (ps->q2sat - ps->q2);
    ef->epsca = (a * ps->rr + rad * delta) / (delta + ps->rr);
    //etp = epsca * rch / LSUBC;
    wf->etp = ef->epsca * ps->rch / lvs / 1000.0;
}
////void RedPrm (grid_struct * grid, lsm_struct lsm, double *zsoil)
////{
////
/////*----------------------------------------------------------------------
////* internally set (default valuess)
////* all soil and vegetation parameters required for the execusion of
////* the noah lsm are defined in vegparm.tbl, soilparm.tb, and genparm.tbl.
////* ----------------------------------------------------------------------
////*     vegetation parameters:
////*             albbrd: sfc background snow-free albedo
////*             cmxtbl: max cnpy capacity
////*              z0brd: background roughness length
////*             shdfac: green vegetation fraction
////*              nroot: rooting depth
////*              rsmin: mimimum stomatal resistance
////*              rsmax: max. stomatal resistance
////*                rgl: parameters used in radiation stress function
////*                 hs: parameter used in vapor pressure deficit functio
////*               topt: optimum transpiration air temperature.
////*             cmcmax: maximum canopy water capacity
////*             cfactr: parameter used in the canopy inteception calculation
////*               snup: threshold snow depth (in water equivalent m) that
////*                     implies 100 percent snow cover
////*                lai: leaf area index
////*
////* ----------------------------------------------------------------------
////*      soil parameters:
////*        smcmax: max soil moisture content (porosity)
////*        smcref: reference soil moisture  (field capacity)
////*        smcwlt: wilting point soil moisture
////*        smcwlt: air dry soil moist content limits
////*       ssatpsi: sat (saturation) soil potential
////*         dksat: sat soil conductivity
////*          bexp: b parameter
////*        ssatdw: sat soil diffusivity
////*           f1: soil thermal diffusivity/conductivity coef.
////*        quartz: soil quartz content
////*  modified by f. chen (12/22/97)  to use the statsgo soil map
////*  modified by f. chen (01/22/00)  to include playa, lava, and white san
////*  modified by f. chen (08/05/02)  to include additional parameters for the noah
////* note: satdw = bb*satdk*(satpsi/maxsmc)
////*         f11 = alog10(satpsi) + bb*alog10(maxsmc) + 2.0
////*       refsmc1=maxsmc*(5.79e-9/satdk)**(1/(2*bb+3)) 5.79e-9 m/s= 0.5 mm
////*       refsmc=refsmc1+1./3.(maxsmc-refsmc1)
////*       wltsmc1=maxsmc*(200./satpsi)**(-1./bb)    (wetzel and chang, 198
////*       wltsmc=wltsmc1-0.5*wltsmc1
////* note: the values for playa is set for it to have a thermal conductivit
////* as sand and to have a hydrulic conductivity as clay
////*
////* ----------------------------------------------------------------------
////* class parameter 'slopetyp' was included to estimate linear reservoir
////* coefficient 'slope' to the baseflow runoff out of the bottom layer.
////* lowest class (slopetyp=0) means highest slope parameter = 1.
////* definition of slopetyp from 'zobler' slope type:
////* slope class  percent slope
////* 1            0-8
////* 2            8-30
////* 3            > 30
////* 4            0-30
////* 5            0-8 & > 30
////* 6            8-30 & > 30
////* 7            0-8, 8-30, > 30
////* 9            glacial ice
////* blank        ocean/sea
////*       slope_data: linear reservoir coefficient
////*       sbeta_data: parameter used to caluculate vegetation effect on soil heat
////*       fxexp_dat:  soil evaporation exponent used in DEvap
////*       csoil_data: soil heat capacity [j m-3 k-1]
////*       salp_data: shape parameter of  distribution function of snow cover
////*       refdk_data and refkdt_data: parameters in the surface runoff parameteriz
////*       frzk_data: frozen ground parameter
////*       zbot_data: depth[m] of lower boundary soil temperature
////*       czil_data: calculate roughness length of heat
////*       smlow_data and mhigh_data: two soil moisture wilt, soil moisture referen
////*                 parameters
////* set maximum number of soil-, veg-, and slopetyp in data statement.
////* --------------------------------------------------------------------*/
////
////    int             i;
////
////    double          frzfact;
////
////    /*
////     * save
////     * * ----------------------------------------------------------------------
////     */
////    if (grid->soiltyp > lsm->soiltbl.slcats)
////    {
////        printf ("warning: too many input soil types\n");
////        PihmExit (0);
////    }
////    if (grid->vegtyp > lsm->vegtbl.lucats)
////    {
////        printf ("warning: too many input landuse types\n");
////        PihmExit (0);
////    }
////    if (grid->slopetyp > lsm->genprmt.slpcats)
////    {
////        printf ("warning: too many input slope types\n");
////        PihmExit (0);
////    }
////
/////*----------------------------------------------------------------------
////*  set-up soil parameters
////* --------------------------------------------------------------------*/
////    grid->csoil = lsm->genprmt.csoil_data;
////#ifdef _NOAH_
////    grid->vgalpha = lsm->soiltbl.vga[grid->soiltyp - 1];
////    grid->vgbeta = lsm->soiltbl.vgb[grid->soiltyp - 1];
////    grid->smcmin = lsm->soiltbl.minsmc[grid->soiltyp - 1];
////    grid->macksat = lsm->soiltbl.macksat[grid->soiltyp - 1];
////    grid->areaf = lsm->soiltbl.areaf[grid->soiltyp - 1];
////    grid->nmacd = lsm->soiltbl.nmacd[grid->soiltyp - 1];
////#else
////    grid->bexp = lsm->soiltbl.bb[grid->soiltyp - 1];
////    grid->psisat = lsm->soiltbl.satpsi[grid->soiltyp - 1];
////    grid->dwsat = lsm->soiltbl.satdw[grid->soiltyp - 1];
////#endif
////    grid->dksat = lsm->soiltbl.satdk[grid->soiltyp - 1];
////    grid->f1 = lsm->soiltbl.f11[grid->soiltyp - 1];
////    grid->quartz = lsm->soiltbl.qtz[grid->soiltyp - 1];
////    grid->smcdry = lsm->soiltbl.drysmc[grid->soiltyp - 1];
////    grid->smcmax = lsm->soiltbl.maxsmc[grid->soiltyp - 1];
////    grid->smcref = lsm->soiltbl.refsmc[grid->soiltyp - 1];
////    grid->smcwlt = lsm->soiltbl.wltsmc[grid->soiltyp - 1];
////
/////*----------------------------------------------------------------------
////* set-up universal parameters (not dependent on soiltyp, vegtyp or
////* slopetyp)
////* --------------------------------------------------------------------*/
////    grid->zbot = lsm->genprmt.zbot_data;
////    grid->salp = lsm->genprmt.salp_data;
////    grid->sbeta = lsm->genprmt.sbeta_data;
////    grid->frzk = lsm->genprmt.frzk_data;
////    grid->fxexp = lsm->genprmt.fxexp_data;
////    grid->ptu = 0.;             /* (not used yet) to satisify intent(out) */
////    grid->czil = lsm->genprmt.czil_data;
////    grid->lvcoef = lsm->genprmt.lvcoef_data;
////#ifndef _NOAH_
////    grid->slope = lsm->genprmt.slope_data[grid->slopetyp - 1];
////    grid->refkdt = lsm->genprmt.refkdt_data;
////    grid->refdk = lsm->genprmt.refdk_data;
////    grid->kdt = grid->refkdt * grid->dksat / grid->refdk;
////#endif
////
/////*----------------------------------------------------------------------
////* to adjust frzk parameter to actual soil type: frzk * frzfact
////* --------------------------------------------------------------------*/
////    frzfact = (grid->smcmax / grid->smcref) * (0.412 / 0.468);
////    grid->frzx = grid->frzk * frzfact;
////
/////*----------------------------------------------------------------------
////* set-up vegetation parameters
////* --------------------------------------------------------------------*/
////    grid->topt = lsm->vegtbl.topt_data;
////    grid->cfactr = lsm->vegtbl.cfactr_data;
////    grid->rsmax = lsm->vegtbl.rsmax_data;
////#ifdef _NOAH_
////    grid->cmcmax = lsm->vegtbl.cmcfactrtbl[grid->vegtyp - 1] * grid->xlai;
////#else
////    grid->cmcmax = lsm->vegtbl.cmcmax_data;
////#endif
////    grid->nroot = lsm->vegtbl.nrotbl[grid->vegtyp - 1];
////    grid->snup = lsm->vegtbl.snuptbl[grid->vegtyp - 1];
////    grid->rsmin = lsm->vegtbl.rstbl[grid->vegtyp - 1];
////    grid->rgl = lsm->vegtbl.rgltbl[grid->vegtyp - 1];
////    grid->hs = lsm->vegtbl.hstbl[grid->vegtyp - 1];
////    grid->emissmin = lsm->vegtbl.emissmintbl[grid->vegtyp - 1];
////    grid->emissmax = lsm->vegtbl.emissmaxtbl[grid->vegtyp - 1];
////    grid->laimin = lsm->vegtbl.laimintbl[grid->vegtyp - 1];
////    grid->laimax = lsm->vegtbl.laimaxtbl[grid->vegtyp - 1];
////    grid->z0min = lsm->vegtbl.z0mintbl[grid->vegtyp - 1];
////    grid->z0max = lsm->vegtbl.z0maxtbl[grid->vegtyp - 1];
////    grid->albedomin = lsm->vegtbl.albedomintbl[grid->vegtyp - 1];
////    grid->albedomax = lsm->vegtbl.albedomaxtbl[grid->vegtyp - 1];
////
////    grid->isurban = lsm->vegtbl.isurban;
////
////    if (grid->vegtyp == lsm->vegtbl.bare)
////        grid->shdfac = 0.0;
////
////    if (grid->nroot > grid->nsoil)
////    {
////        printf ("error: too many root layers %d, %d\n", grid->nsoil,
////           grid->nroot);
////        PihmExit (0);
////    }
////
/////*----------------------------------------------------------------------
////* calculate root distribution.  present version assumes uniform
////* distribution based on soil layer depths.
////* --------------------------------------------------------------------*/
////    for (i = 0; i < grid->nroot; i++)
////        grid->rtdis[i] = -grid->sldpth[i] / zsoil[grid->nroot - 1];
////
/////*----------------------------------------------------------------------
////*  set-up slope parameter
////* --------------------------------------------------------------------*/
////
////    /*
////     * print*,'end of prmred'
////     * *       print*,'vegtyp',vegtyp,'soiltyp',soiltyp,'slopetyp',slopetyp,    &
////     * *    & 'cfactr',cfactr,'cmcmax',cmcmax,'rsmax',rsmax,'topt',topt,        &
////     * *    & 'refkdt',refkdt,'kdt',kdt,'sbeta',sbeta, 'shdfac',shdfac,         &
////     * *    &  'rsmin',rsmin,'rgl',rgl,'hs',hs,'zbot',zbot,'frzx',frzx,         &
////     * *    &  'psisat',psisat,'slope',slope,'snup',snup,'salp',salp,'bexp',    &
////     * *    &   bexp,                                                           &
////     * *    &  'dksat',dksat,'dwsat',dwsat,                                     &
////     * *    &  'smcmax',smcmax,'smcwlt',smcwlt,'smcref',smcref,'smcdry',smcdry, &
////     * *    &  'f1',f1,'quartz',quartz,'fxexp',fxexp,                           &
////     * *    &  'rtdis',rtdis,'sldpth',sldpth,'zsoil',zsoil, 'nroot',nroot,      &
////     * *    &  'nsoil',nsoil,'z0',z0,'czil',czil,'lai',lai,                     &
////     * *    &  'csoil',csoil,'ptu',ptu,                                         &
////     * *    &  'local', local
////     */
////
////}
//
void Rosr12 (double *p, double *a, double *b, double *c, double *d,
    double *delta, int nsoil)
{
    /*
     * Function Rosr12
     *
     * Invert (solve) the tri-diagonal matrix problem shown below:
     * ###                                            ### ###  ###   ###  ###
     * #b[0], c[0],  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
     * #a[1], b[1], c[1],  0  ,  0  ,   . . .  ,    0   # #      #   #      #
     * # 0  , a[2], b[2], c[2],  0  ,   . . .  ,    0   # #      #   # d[2] #
     * # 0  ,  0  , a[3], b[3], c[3],   . . .  ,    0   # # p[3] #   # d[3] #
     * # 0  ,  0  ,  0  , a[4], b[4],   . . .  ,    0   # # p[4] #   # d[4] #
     * # .                                          .   # #  .   # = #   .  #
     * # .                                          .   # #  .   #   #   .  #
     * # .                                          .   # #  .   #   #   .  #
     * # 0  , . . . , 0 , a[m-2], b[m-2], c[m-2],   0   # #p(m-2)#   #d[m-2]#
     * # 0  , . . . , 0 ,   0   , a[m-1], b[m-1], c[m-1]# #p(m-1)#   #d[m-1]#
     * # 0  , . . . , 0 ,   0   ,   0   ,  a(m) ,  b[m] # # p(m) #   # d[m] #
     * ###                                            ### ###  ###   ###  ###
     * --------------------------------------------------------------------*/

    int             k, kk;

    /* Initialize eqn coef c for the lowest soil layer */
    c[nsoil - 1] = 0.0;
    p[0] = -c[0] / b[0];

    /* Solve the coefs for the 1st soil layer */
    delta[0] = d[0] / b[0];

    /* Solve the coefs for soil layers 2 thru nsoil */
    for (k = 1; k < nsoil; k++)
    {
        p[k] = -c[k] * (1.0 / (b[k] + a[k] * p[k - 1]));
        delta[k] =
            (d[k] - a[k] * delta[k - 1]) * (1.0 / (b[k] + a[k] * p[k - 1]));
    }

    /* Set p to delta for lowest soil layer */
    p[nsoil - 1] = delta[nsoil - 1];

    /* Adjust p for soil layers 2 thru nsoil */
    for (k = 1; k < nsoil; k++)
    {
        kk = nsoil - k - 1;
        p[kk] = p[kk] * p[kk + 1] + delta[kk];
    }
}

void ShFlx (ws_struct *ws, es_struct *es, ef_struct *ef, ps_struct *ps,
    const lc_struct *lc, const soil_struct *soil, double dt, double yy,
    double zz1, const double *zsoil, double df1)
{
    /*
     * Function ShFlx
     *
     * Update the temperature state of the soil column based on the thermal
     * diffusion equation and update the frozen soil moisture content based
     * on the temperature.
     */
    double          ai[MAXLYR], bi[MAXLYR], ci[MAXLYR];
    double          rhsts[MAXLYR];

    /* HRT routine calcs the right hand side of the soil temp dif eqn
     * Land case */
    HRT (ws, es, ef, ps, lc, soil, rhsts, zsoil, yy, zz1, dt, df1, ai, bi,
        ci);

    HStep (es, rhsts, dt, ps->nsoil, ai, bi, ci);

    /* In the no snowpack case (via routine NoPac branch,) update the grnd
     * (skin) temperature here in response to the updated soil temperature
     * profile above.  (note: inspection of routine SnoPac shows that t1
     * below is a dummy variable only, as skin temperature is updated
     * differently in routine SnoPac)
     *
     * Calculate surface soil heat flux */
    es->t1 = (yy + (zz1 - 1.0) * es->stc[0]) / zz1;
    ef->ssoil = df1 * (es->stc[0] - es->t1) / (0.5 * zsoil[0]);
}

void SmFlx (ws_struct *ws, wf_struct *wf, const wf_struct *avgwf,
    ps_struct *ps, const lc_struct *lc, const soil_struct *soil,
    const double *zsoil, double prcp, double dt)
{
    /*
     * Function SmFlx
     *
     * Calculate soil moisture flux. The soil moisture content (smc - a per
     * unit volume measurement) is a dependent variable that is updated with
     * prognostic eqns. The canopy moisture content (cmc) is also updated.
     * Frozen ground version: new states added: sh2o, and frozen ground
     * correction factor, frzfact and parameter slope.
     */
    int             i;
    double          ai[MAXLYR], bi[MAXLYR], ci[MAXLYR];
    double          rhstt[MAXLYR];
    double          sice[MAXLYR];
    double          excess;
    double          rhsct;
    double          trhsct;
    const double    KD = 6.54e-7;
    const double    BFACTR = 3.89;

    /* Executable code begins here.
     * Compute the right hand side of the canopy eqn term (rhsct) */

    /* Convert rhsct (a rate) to trhsct (an amount) and add it to existing
     * cmc. If resulting amt exceeds max capacity, it becomes drip and will
     * fall to the grnd. */
    rhsct = lc->shdfac * prcp - wf->ec;
    wf->drip = 0.0;
    trhsct = dt * rhsct;
    excess = ws->cmc + trhsct;

    /* pcpdrp is the combined prcp and drip (from cmc) that goes into the
     * soil */

    /* PIHM drip calculation following rutter and mortan (1977 jae) */
    if (excess > 0.0)
    {
        if (excess >= ws->cmcmax)
        {
            wf->drip =
                (KD * ws->cmcmax * exp (BFACTR)) + (excess - ws->cmcmax) / dt;
            rhsct -= KD * ws->cmcmax * exp (BFACTR);
        }
        else
        {
            wf->drip =
                (KD * ws->cmcmax * exp (BFACTR * excess / ws->cmcmax));
            rhsct -= KD * ws->cmcmax * exp (BFACTR * excess / ws->cmcmax);
        }
    }
    wf->pcpdrp = (1.0 - lc->shdfac) * prcp + wf->drip;

    /* Store ice content at each soil layer before calling SRT and SStep */
    for (i = 0; i < ps->nsoil; i++)
    {
        sice[i] = ws->smc[i] - ws->sh2o[i];
    }

    /* Call subroutines SRT and SStep to solve the soil moisture tendency
     * equations.
     * Call the SRT/SStep function in the manner of time scheme "d" (implicit
     * state, explicit coefficient) of Section 2 of Kalnay and Kanamitsu
     * pcpdrp is units of m/s, zsoil is negative depth in m
     * According to Dr. Ken Mitchell's suggestion, add the second contraint
     * to remove numerical instability of runoff and soil moisture
     * flimit is a limit value for fac2
     * Frozen ground version:
     * smc states replaced by sh2o states in SRT subr. sh2o & sice states
     * included in SStep subr. Frozen ground correction factor, frzfact added.
     * All water balance calculations using unfrozen water */
    if (ps->nwtbl == 0)
    {
        for (i = 0; i < ps->nsoil; i++)
        {
            ws->smc[i] = soil->smcmax;
            ws->sh2o[i] = ws->smc[i] - sice[i];
        }
    }
    else
    {
        SRT (ws, wf, avgwf, ps, soil, rhstt, sice, ai, bi, ci, zsoil, dt);
        SStep (ws, wf, ps, soil, rhstt, rhsct, zsoil, sice, ai, bi, ci, dt);
    }
}

double SnFrac (double sneqv, double snup, double salp, double snowh)
{
    /*
     * Function SnFrac
     *
     * Calculate snow fraction (0 -> 1)
     */
    double          rsnow;
    double          sncovr;

    /* snup is veg-class dependent snowdepth threshhold (set in routine
     * RedPrm) above which snocvr = 1. */
    if (sneqv < snup)
    {
        rsnow = sneqv / snup;
        sncovr = 1.0 - (exp (- salp * rsnow) - rsnow * exp (- salp));
    }
    else
    {
        sncovr = 1.0;
    }

    /* Formulation of Dickinson et al. 1986 */
    //z0n = 0.035;

    //sncovr = snowh / (snowh + 5.0 * z0n);

    /* Formulation of Marshall et al. 1994 */
    //sncovr = sneqv / (sneqv + 2.0 * z0n);

    return (sncovr);
}

void SnkSrc (double *tsnsr, double tavg, double smc, double *sh2o,
    const soil_struct *soil, const double *zsoil, int nsoil, double dt, int k,
    double qtot)
{
    /*
     * Subroutine SnkSrc
     *
     * Calculate sink/source term of the thermal diffusion equation. (sh2o)
     * is available liqued water.
     */
    double          dz;
    double          freew;
    double          xh2o;

    double          DH2O = 1.0000e3;

    if (k == 0)
    {
        dz = - zsoil[0];
    }
    else
    {
        dz = zsoil[k - 1] - zsoil[k];
    }

    /* Via function FrH2O, compute potential or 'equilibrium' unfrozen
     * supercooled free water for given soil type and soil layer temperature.
     * Function FrH2O invokes Eqn (17) from V. Koren et al (1999, JGR, Vol.
     * 104, Pg 19573). (Aside: latter eqn in journal in centigrade units.
     * routine FrH2O use form of eqn in kelvin units.) */
    freew = FrH2O (tavg, smc, *sh2o, soil);

    /* In next block of code, invoke Eqn 18 of V. Koren et al (1999, JGR,
     * Vol. 104, Pg 19573.) that is, first estimate the new amountof liquid
     * water, 'xh2o', implied by the sum of (1) the liquid water at the begin
     * of current time step, and (2) the freeze of thaw change in liquid water
     * implied by the heat flux 'qtot' passed in from routine HRT. Second,
     * determine if xh2o needs to be bounded by 'freew' (equil amt) or if
     * 'freew' needs to be bounded by xh2o. */

    xh2o = *sh2o + qtot * dt / (DH2O * LSUBF * dz);

    /* First, if freezing and remaining liquid less than lower bound, then
     * reduce extent of freezing, thereby letting some or all of heat flux
     * qtot cool the soil temp later in routine HRT. */
    if (xh2o < *sh2o && xh2o < freew)
    {
        if (freew > *sh2o)
        {
            xh2o = *sh2o;
        }
        else
        {
            xh2o = freew;
        }
    }

    /* Second, if thawing and the increase in liquid water greater than upper
     * bound, then reduce extent of thaw, thereby letting some or all of heat
     * flux qtot warm the soil temp later in routine HRT. */
    if (xh2o > *sh2o && xh2o > freew)
    {
        if (freew < *sh2o)
        {
            xh2o = *sh2o;
        }
        else
        {
            xh2o = freew;
        }
    }

    /* Calculate phase-change heat source/sink term for use in routine HRT
     * and update liquid water to reflcet final freeze/thaw increment. */
    if (xh2o < 0.0)
    {
        xh2o = 0.0;
    }
    if (xh2o > smc)
    {
        xh2o = smc;
    }

    *tsnsr = - DH2O * LSUBF * dz * (xh2o - *sh2o) / dt;

    *sh2o = xh2o;
}

void SnoPac (ws_struct *ws, wf_struct *wf, const wf_struct *avgwf,
    es_struct *es, ef_struct *ef, ps_struct *ps, lc_struct *lc,
    const soil_struct *soil, int snowng, const double *zsoil, double dt,
    double t24, double prcpf, double df1)
{
    /*
     * Function SnoPac
     *
     * Calculate soil moisture and heat flux values & update soil moisture
     * content and soil heat content values for the case when a snow pack is
     * present.
     */
    int             k;

    double          denom;
    double          dsoil;
    double          dtot;
    double          esnow1;
    double          esnow2;
    double          etp3;
    double          etanrg;
    double          ex;
    double          seh;
    double          sncond;
    double          t12;
    double          t12a;
    double          t12b;
    double          t14;
    double         ssoil1;
    double         t11;
    double         yy;
    double         zz1;
    const double    ESDMIN = 1.0e-6;
    const double    SNOEXP = 2.0;

    /* Executable code begins here:
     * Initialize evap terms.
     * conversions:
     * esnow [kg m-2 s-1]
     * esdflx [kg m-2 s-1] .le. esnow
     * esnow1 [m s-1]
     * esnow2 [m]
     * etp [kg m-2 s-1]
     * etp1 [m s-1]
     * etp2 [m]*/
    wf->dew = 0.0;
    wf->edir = 0.0;
    wf->ec = 0.0;

    for (k = 0; k < ps->nsoil; k++)
    {
        wf->et[k] = 0.0;
    }
    wf->ett = 0.0;
    wf->esnow = 0.0;
    esnow1 = 0.0;
    esnow2 = 0.0;

    /* Convert potential evap (etp) from kg m-2 s-1 to etp1 in m s-1
     * if etp < 0 (downward) then dewfall (= frostfall in this case). */
    ps->beta = 1.0;

    if (wf->etp <= 0.0)
    {
        if ((ps->ribb >= 0.1) && (ef->fdown > 150.0))
        {
            wf->etp =
                ((wf->etp * (1.0 - ps->ribb) < 0.0 ?
                    wf->etp * (1.0 - ps->ribb) : 0.0)
                    * ps->sncovr / 0.980 + wf->etp * (0.980 - ps->sncovr)) / 0.980;
        }
        if (wf->etp == 0.0)
        {
            ps->beta = 0.0;
        }
        //etp1 = *etp * 0.001;
        wf->dew = -wf->etp;
        esnow2 = wf->etp * dt;
        etanrg = wf->etp * 1000.0 * ((1.0 - ps->sncovr) * LVH2O +
            ps->sncovr * LSUBS);
    }
    else
    {
        //*etp1 = *etp * 0.001;
        /* Land case */
        if (ps->sncovr < 1.0)
        {
            Evapo (ws, wf, ps, lc, soil, zsoil, dt);

            wf->edir *= (1.0 - ps->sncovr);
            wf->ec *= (1.0 - ps->sncovr);
            for (k = 0; k < ps->nsoil; k++)
            {
                wf->et[k] = wf->et[k] * (1.0 - ps->sncovr);
            }
            wf->ett *= (1.0 - ps->sncovr);
            wf->etns *= (1.0 - ps->sncovr);
        }
        wf->esnow = wf->etp * ps->sncovr;
        //esnow1 = *esnow * 0.001;
        esnow2 = wf->esnow * dt;
        etanrg = wf->esnow * 1000.0 * LSUBS + wf->etns * 1000.0 * LVH2O;
    }

    /* If precip is falling, calculate heat flux from snow sfc to newly
     * accumulating precip.  note that this reflects the flux appropriate for
     * the not-yet-updated skin temperature (t1).  assumes temperature of the
     * snowfall striking the ground is = sfctmp (lowest model level air
     * temp). */
    ef->flx1 = 0.0;
    if (snowng)
    {
        ef->flx1 = CPICE * wf->prcp * 1000.0 * (es->t1 - es->sfctmp);
    }
    else
    {
        if (wf->prcp > 0.0)
        {
            ef->flx1 = CPH2O * wf->prcp * 1000.0 * (es->t1 - es->sfctmp);
        }
    }

    /* Calculate an 'effective snow-grnd sfc temp' (t12) based on heat fluxes
     * between the snow pack and the soil and on net radiation.
     * Include flx1 (precip-snow sfc) and flx2 (freezing rain latent heat)
     * fluxes. flx1 from above, flx2 brought in via commom block rite.
     * flx2 reflects freezing rain latent heat flux using t1 calculated in
     * Penman. */
    dsoil = - (0.5 * zsoil[0]);
    dtot = ps->snowh + dsoil;
    denom = 1.0 + df1 / (dtot * ps->rr * ps->rch);

    /* Surface emissivity weighted by snow cover fraction */
    t12a =
        ((ef->fdown - ef->flx1 - ef->flx2 - ps->emissi * SIGMA * t24) / ps->rch +
        es->th2 - es->sfctmp - etanrg / ps->rch) / ps->rr;
    t12b = df1 * es->stc[0] / (dtot * ps->rr * ps->rch);

    /* If the 'effective snow-grnd sfc temp' is at or below freezing, no snow
     * melt will occur. Set the skin temp to this effective temp. Reduce (by
     * sublimination ) or increase (by frost) the depth of the snowpack,
     * depending on sign of etp.
     * Update soil heat flux (ssoil) using new skin temperature (t1) since no
     * snowmelt, set accumulated snowmelt to zero, set 'effective' precip from
     * snowmelt to zero, set phase-change heat flux from snowmelt to zero. */
    t12 = (es->sfctmp + t12a + t12b) / denom;

    if (t12 <= TFREEZ)
    {
        /* Sub-freezing block */
        es->t1 = t12;
        ef->ssoil = df1 * (es->t1 - es->stc[0]) / dtot;

        ws->sneqv = (ws->sneqv - esnow2 > 0.0) ? ws->sneqv - esnow2 : 0.0;
        ef->flx3 = 0.0;
        ex = 0.0;

        wf->snomlt = 0.0;
    }
    /* If the 'effective snow-grnd sfc temp' is above freezing, snow melt
     * will occur.  call the snow melt rate, ex and amt, snomlt. Revise the
     * effective snow depth. Revise the skin temp because it would have chgd
     * due to the latent heat released by the melting. Calc the latent heat
     * released, flx3. Set the effective precip, prcp1 to the snow melt rate,
     * ex for use in SmFlx. Adjustment to t1 to account for snow patches.
     * Calculate qsat valid at freezing point. Note that esat (saturation
     * vapor pressure) value of 6.11e+2 used here is that valid at frzzing
     * point. Note that etp from call Penman in sflx is ignored here in
     * favor of bulk etp over 'open water' at freezing temp. Update soil heat
     * flux (s) using new skin temperature (t1) */
    else
    {
        /* Above freezing block */
        es->t1 = TFREEZ * pow (ps->sncovr, SNOEXP) +
            t12 * (1.0 - pow (ps->sncovr, SNOEXP));
        ps->beta = 1.0;

        /* If potential evap (sublimation) greater than depth of snowpack.
         * beta < 1
         * snowpack has sublimated away, set depth to zero. */
        ef->ssoil = df1 * (es->t1 - es->stc[0]) / dtot;

        if (ws->sneqv - esnow2 <= ESDMIN)
        {
            ws->sneqv = 0.0;
            ex = 0.0;
            wf->snomlt = 0.0;
            ef->flx3 = 0.0;

        }
        else
        {
            /* Sublimation less than depth of snowpack
             * snowpack (esd) reduced by esnow2 (depth of sublimated snow) */
            ws->sneqv -= esnow2;
            etp3 = wf->etp * 1000.0 * LVH2O;
            seh = ps->rch * (es->t1 - es->th2);
            t14 = es->t1 * es->t1;
            t14 = t14 * t14;
            ef->flx3 =
                ef->fdown - ef->flx1 - ef->flx2 - ps->emissi * SIGMA * t14 - ef->ssoil -
                seh - etanrg;
            if (ef->flx3 <= 0.0)
            {
                ef->flx3 = 0.0;
            }

            /* Snowmelt reduction depending on snow cover */
            ex = ef->flx3 * 0.001 / LSUBF;

            /* ESDMIN represents a snowpack depth threshold value below which
             * we choose not to retain any snowpack, and instead include it in
             * snowmelt. */
            //snomlt = ex * *dt;
            wf->snomlt = ex;
            if (ws->sneqv - wf->snomlt * dt >= ESDMIN)
            {
                ws->sneqv -= wf->snomlt * dt;
            }
            else
            {
                /* Snowmelt exceeds snow depth */
                ex = ws->sneqv / dt;
                ef->flx3 = ex * 1000.0 * LSUBF;
                wf->snomlt = ws->sneqv / dt;

                ws->sneqv = 0.0;
            }
        }

        /* If non-glacial land, add snowmelt rate (ex) to precip rate to be
         * used in subroutine SmFlx (soil moisture evolution) via
         * infiltration.
         * Runoff/baseflow later near the end of sflx (after return from call
         * to subroutine SnoPac) */
        prcpf += ex;

        /* Set the effective potnl evapotransp (etp1) to zero since this is
         * snow case, so surface evap not calculated from edir, ec, or ett in
         * SmFlx (below).
         * SmFlx returns updated soil moisture values for non-glacial land. */
    }

    SmFlx (ws, wf, avgwf, ps, lc, soil, zsoil, prcpf, dt);

    /* Before call ShFlx in this snowpack case, set zz1 and yy arguments to
     * special values that ensure that ground heat flux calculated in ShFlx
     * matches that already computer for below the snowpack, thus the sfc
     * heat flux to be computed in ShFlx will effectively be the flux at the
     * snow top surface.  t11 is a dummy arguement so we will not use the
     * skin temp value as revised by ShFlx. */
    zz1 = 1.0;
    yy = es->stc[0] - 0.5 * ef->ssoil * zsoil[0] * zz1 / df1;

    /* ShFlx will calc/update the soil temps. Note: the sub-sfc heat flux
     * (ssoil1) and the skin temp (t11) output from this ShFlx call are not
     * used  in any subsequent calculations. Rather, they are dummy variables
     * here in the SnoPac case, since the skin temp and sub-sfc heat flux are
     * updated instead near the beginning of the call to SnoPac. */
    t11 = es->t1;
    ssoil1 = ef->ssoil;

    ShFlx (ws, es, ef, ps, lc, soil, dt, yy, zz1, zsoil, df1);

    es->t1 = t11;
    ef->ssoil = ssoil1;

    /* Snow depth and density adjustment based on snow compaction. yy is
     * assumed to be the soil temperture at the top of the soil column. */
    if (ws->sneqv > 0.0)
    {
        SnowPack (ws->sneqv, dt, &ps->snowh, &ps->sndens, es->t1, yy);
    }
    else
    {
        ws->sneqv = 0.0;
        ps->snowh = 0.0;
        ps->sndens = 0.0;
        sncond = 1.0;
        ps->sncovr = 0.0;
    }
}

void SnowPack (double esd, double dtsec, double *snowh, double *sndens,
    double tsnow, double tsoil)
{
    /*
     * Function SnowPack
     *
     * Calculate compaction of snowpack under conditions of increasing snow
     * density, as obtained from an approximate solution of E. Anderson's
     * differential equation (3.29), NOAA technical report NWS 19, by Victor
     * Koren, 03/25/95.
     *
     * esd     water equivalent of snow (m)
     * dtsec   time step (sec)
     * snowh   snow depth (m)
     * sndens  snow density (g/cm3=dimensionless fraction of h2o density)
     * tsnow   snow surface temperature (k)
     * tsoil   soil surface temperature (k)
     * Function will return new values of snowh and sndens
     */
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
    const double    C1 = 0.01;
    const double    C2 = 21.0;

    /* Conversion into simulation units */
    snowhc = *snowh * 100.0;
    esdc = esd * 100.0;
    dthr = dtsec / 3600.0;
    tsnowc = tsnow - 273.15;
    tsoilc = tsoil - 273.15;

    /* Calculating of average temperature of snow pack */
    tavgc = 0.5 * (tsnowc + tsoilc);

    /* Calculating of snow depth and density as a result of compaction
     *  sndens = ds0 * (exp (bfac * esd) - 1.) / (bfac * esd)
     *  bfac = dthr * c1 * exp (0.08 * tavgc - c2 * ds0)
     * Note: bfac * esd in sndens eqn above has to be carefully treated
     * numerically below:
     *   c1 is the fractional increase in density (1/(cm*hr))
     *   c2 is a constant (cm3/g) kojima estimated as 21 cms/g */
    if (esdc > 1.0e-2)
    {
        esdcx = esdc;
    }
    else
    {
        esdcx = 1.0e-2;
    }

    /* The function of the form (e**x - 1) / x embedded in above expression
     * for dsx was causing numerical difficulties when the denominator "x"
     * (i.e. bfac*esdc) became zero or approached zero (despite the fact that
     * the analytical function (e**x-1)/x has a well defined limit as
     * "x" approaches zero), hence below we replace the (e**x-1)/x
     * expression with an equivalent, numerically well-behaved
     * polynomial expansion.
     * Number of terms of polynomial expansion, and hence its accuracy,
     * is governed by iteration limit "ipol".
     *      ipol greater than 9 only makes a difference on double
     *            precision (relative errors given in percent %).
     *       ipol=9, for rel.error <~ 1.6 e-6 % (8 significant digits)
     *       ipol=8, for rel.error <~ 1.8 e-5 % (7 significant digits)
     *       ipol=7, for rel.error <~ 1.8 e-4 % ... */
    bfac = dthr * C1 * exp (0.08 * tavgc - C2 * *sndens);
    ipol = 4;
    pexp = 0.0;

    for (j = ipol; j > 0; j--)
    {
        pexp = (1.0 + pexp) * bfac * esdcx / (double)(j + 1);
    }

    pexp += 1.0;

    /* Set upper/lower limit on snow density */
    dsx = *sndens * pexp;
    dsx = (dsx > 0.40) ? 0.40: dsx;
    dsx = (dsx < 0.05) ? 0.05 : dsx;

    /* Update of snow depth and density depending on liquid water during
     * snowmelt. Assumed that 13% of liquid water can be stored in snow per
     * day during snowmelt till snow density 0.40. */
    *sndens = dsx;
    if (tsnowc >= 0.0)
    {
        dw = 0.13 * dthr / 24.0;
        *sndens = *sndens * (1.0 - dw) + dw;
        if (*sndens >= 0.40)
        {
            *sndens = 0.40;
        }
    }
    /* Calculate snow depth (cm) from snow water equivalent and snow density.
     * Change snow depth units to meters */
    snowhc = esdc / *sndens;
    *snowh = snowhc * 0.01;
}

double Snowz0 (double sncovr, double z0brd, double snowh)
{
    /*
     * Function Snowz0
     *
     * Calculate total roughness length over snow
     * z0s     snow roughness length:=0.001 (m)
     */
    const double    Z0S = 0.001;
    double          burial;
    double          z0eff;
    double          z0;

    //z0 = (1.0 - sncovr)* z0brd + sncovr * Z0S;
    burial = 7.0 * z0brd - snowh;
    if (burial < 0.0007)
    {
        z0eff = Z0S;
    }
    else
    {
        z0eff = burial / 7.0;
    }

    z0 = (1.0 - sncovr) * z0brd + sncovr * z0eff;

    return (z0);
}

void SnowNew (const es_struct *es, double newsn, ps_struct *ps)
{
    /*
     * Function SnowNew
     *
     * Calculate snow depth and density to account for the new snowfall.
     * New values of snow depth & density returned.
     */

    double          dsnew;
    double          hnewc;
    double          snowhc;
    double          newsnc;
    double          tempc;

    /* Conversion into simulation units */
    snowhc = ps->snowh * 100.0;
    newsnc = newsn * 100.0;

    /* Calculating new snowfall density depending on temperature equation
     * from Gottlib L. 'A general runoff model for snowcovered and
     * glacierized basin', 6th Nordic Hydrological Conference, Vemadolen,
     * Sweden, 1980, 172-177pp. */
    tempc = es->sfctmp - 273.15;
    if (tempc <= -15.0)
    {
        dsnew = 0.05;
    }
    else
    {
        dsnew = 0.05 + 0.0017 * pow (tempc + 15.0, 1.5);
    }

    /* Adjustment of snow density depending on new snowfall */
    hnewc = newsnc / dsnew;
    if (snowhc + hnewc < 0.001)
    {
        ps->sndens = (dsnew > ps->sndens) ? dsnew : ps->sndens;
    }
    else
    {
        ps->sndens = (snowhc * ps->sndens + hnewc * dsnew) / (snowhc + hnewc);
    }
    snowhc = snowhc + hnewc;
    ps->snowh = snowhc * 0.01;
}

void SRT (ws_struct *ws, wf_struct *wf, const wf_struct *avgwf, ps_struct *ps,
    const soil_struct *soil, double *rhstt, double *sice, double *ai,
    double *bi, double *ci, const double *zsoil, double dt)
{
    /*
     * Function SRT
     * Calculate the right hand side of the time tendency term of the soil
     * water diffusion equation. Also to compute (prepare) the matrix
     * coefficients for the tri-diagonal matrix of the implicit time scheme.
     */
    int             iohinf;
    int             j, jj, k, ks;
    double          ddz;
    double          ddz2;
    double          denom;
    double          denom2;
    double          numer;
    double          pddum;
    double          mxsmc, mxsmc2;
    double          sicemax;
    double          wcnd;
    double          wcnd2;
    double          wdf;
    double          wdf2;
    double          dsmdz,  dsmdz2;
    double          dice;
    const double    CVFRZ = 3.0;
    double          acrt;
    double          sum;
    double          ialp1;
    double          weight[MAXLYR];
    double          smctot;
    int             macpore[MAXLYR];
    /* Frozen ground version:
     * Reference frozen ground parameter, cvfrz, is a shape parameter of areal
     * distribution function of soil ice content which equals 1/cv.
     * cv is a coefficient of spatial variation of soil ice content. Based on
     * field data cv depends on areal mean of frozen depth, and it close to
     * constant = 0.6 if areal mean frozen depth is above 20 cm. That is why
     * parameter cvfrz = 3 (int{1/0.6*0.6}). Current logic doesn't allow cvfrz
     * be bigger than 3 */

    /* Let sicemax be the greatest, if any, frozen water content within soil
     * layers. */
    iohinf = 1;
    sicemax = 0.0;

    for (ks = 0; ks < ps->nsoil; ks++)
    {
        sicemax = (sice[ks] > sicemax) ? sice[ks] : sicemax;
    }

    /* Calculate infiltration reduction factor due to frozen soil */
    dice = -zsoil[0] * sice[0];
    for (ks = 1; ks < ps->nsoil; ks++)
    {
        dice += (zsoil[ks - 1] - zsoil[ks]) * sice[ks];
    }

    ps->fcr = 1.0;

    if (dice > 1.0e-2)
    {
        acrt = CVFRZ * ps->frzx / dice;
        sum = 1.0;
        ialp1 = (int)CVFRZ - 1;
        for (j = 1; j < ialp1 + 1; j++)
        {
            k = 1;
            for (jj = j + 1; jj < ialp1; jj++)
            {
                k *= jj;
            }
            sum += pow (acrt, CVFRZ - (double)j) / (double)k;
        }

        ps->fcr = 1.0 - exp (-acrt) * sum;
    }

    /* Determine rainfall infiltration rate and runoff */
    pddum = avgwf->infil;

    for (k = 0; k < ps->nsoil; k++)
    {
        macpore[k] = 0;
    }
    for (k = 0; k < ps->nmacd - 1; k++)
    {
        macpore[k] = 1;
    }

    /* Determine runoff from each layer */
    smctot = (ps->nwtbl <= 1) ? -zsoil[0] * ws->sh2o[0] : 0.0;
    for (ks = 1; ks < ps->nsoil; ks++)
    {
        if (ks >= ps->nwtbl - 1)
        {
            smctot += (zsoil[ks - 1] - zsoil[ks]) * ws->sh2o[ks];
        }
    }

    for (ks = 0; ks < MAXLYR; ks++)
    {
        weight[ks] = 0.0;
    }

    for (ks = 0; ks < ps->nsoil; ks++)
    {
        if (ks >= ps->nwtbl - 1)
        {
            if (ks == 0)
            {
                weight[ks] = -zsoil[0] * ws->sh2o[0] / smctot;
            }
            else
            {
                weight[ks] = (zsoil[ks - 1] - zsoil[ks]) * ws->sh2o[ks] / smctot;
            }
        }
    }

    if (weight[0] + weight[1] + weight[2] + weight[3] + weight[4] + weight[5] + weight[6] + weight[7] + weight[8] + weight[9] + weight[10] < 0.99
        ||weight[0] + weight[1] + weight[2] + weight[3] + weight[4] + weight[5] + weight[6] + weight[7] + weight[8] + weight[9] + weight[10] > 1.01 )
    printf ("total weight %lf\n", weight[0] + weight[1] + weight[2] + weight[3] + weight[4] + weight[5] + weight[6] + weight[7] + weight[8] + weight[9] + weight[10]);

    mxsmc = ws->sh2o[0];

    dsmdz = (ws->sh2o[0] - ws->sh2o[1]) / (- 0.5 * zsoil[1]);
    WDfCnd (&wdf, &wcnd, mxsmc, sicemax, dsmdz, macpore[0], soil, ps);

    /* Calc the matrix coefficients ai, bi, and ci for the top layer */
    ddz = 1.0 / (- 0.5 * zsoil[1]);
    ai[0] = 0.0;
    bi[0] = wdf * ddz / (- zsoil[0]);
    ci[0] = -bi[0];
    /* Calc rhstt for the top layer after calc'ng the vertical soil moisture
     * gradient btwn the top and next to top layers. */
    rhstt[0] = (wdf * dsmdz + wcnd - pddum + wf->edir + wf->et[0]) / zsoil[0];

    rhstt[0] += avgwf->runoff2 * weight[0] / zsoil[0];

    /* Loop thru the remaining soil layers, repeating the abv process */
    /* Initialize ddz2 */
    ddz2 = 0.0;
    for (k = 1; k < ps->nsoil; k++)
    {
        denom2 = (zsoil[k - 1] - zsoil[k]);
        if (k < ps->nsoil - 1)
        {
            mxsmc2 = ws->sh2o[k];
            denom = zsoil[k - 1] - zsoil[k + 1];
            dsmdz2 = (ws->sh2o[k] - ws->sh2o[k + 1]) / (denom * 0.5);
            WDfCnd (&wdf2, &wcnd2, mxsmc2, sicemax, dsmdz2, macpore[k], soil,
                ps);
            /* Calc some partial products for later use in calc'ng rhstt
             * Calc the matrix coef, ci, after calc'ng its partial product */
            ddz2 = 2.0 / denom;
            ci[k] = - wdf2 * ddz2 / denom2;
        }
        else
        {
            /* Retrieve the soil water diffusivity and hydraulic conductivity
             * for this layer */
            wdf2 = 0.0;
            wcnd2 = 0.0;

            /* Calc a partial product for later use in calc'ng rhstt */
            dsmdz2 = 0.0;
            /* Set matrix coef ci to zero */
            ci[k] = 0.0;
        }

        /* Calc rhstt for this layer after calc'ng its numerator */
        numer = (wdf2 * dsmdz2) + wcnd2 - (wdf * dsmdz) - wcnd + wf->et[k];
        
        numer = numer + avgwf->runoff2 * weight[k];

        rhstt[k] = numer / (-denom2);

        /* Calc matrix coefs, ai, and bi for this layer */
        ai[k] = - wdf * ddz / denom2;
        bi[k] = -(ai[k] + ci[k]);

        /* Reset values of wdf, wcnd, dsmdz, and ddz for loop to next lyr */
        if (k != ps->nsoil - 1)
        {
            wdf = wdf2;
            wcnd = wcnd2;
            dsmdz = dsmdz2;
            ddz = ddz2;
        }
    }
}

void SStep (ws_struct *ws, wf_struct *wf, ps_struct *ps,
    const soil_struct *soil, double *rhstt, double rhsct, const double *zsoil,
    double *sice, double *ai, double *bi, double *ci, double dt)
{
    /*
     * Function SStep
     *
     * Calculate/update soil moisture content values and canopy moisture
     * content values.
     */
    int             k, kk11;

    double          rhsttin[MAXLYR], ciin[MAXLYR];
    double          sh2omid[MAXLYR];
    double          ddz, stot, wplus;

    /* Create 'amount' values of variables to be input to the tri-diagonal
     * matrix routine. */
    for (k = 0; k < ps->nsoil; k++)
    {
        rhstt[k] *= dt;
        ai[k] *= dt;
        bi[k] *= dt;
        bi[k] += 1.0;
        ci[k] *= dt;
    }

    /* Copy values for input variables before call to Rosr12 */
    for (k = 0; k < ps->nsoil; k++)
    {
        rhsttin[k] = rhstt[k];
    }
    for (k = 0; k < ps->nsoil; k++)
    {
        ciin[k] = ci[k];
    }

    /* Call Rosr12 to solve the tri-diagonal matrix */
    Rosr12 (ci, ai, bi, ciin, rhsttin, rhstt, ps->nsoil);

    /* Sum the previous smc value and the matrix solution to get a new value.
     * Min allowable value of smc will be 0.02.
     * Runoff3: runoff within soil layers */
    wplus = 0.0;
    wf->runoff3 = 0.0;

    for (k = ps->nsoil - 1; k >= 0; k--)
    {
        if (k != 0)
        {
            ddz = zsoil[k - 1] - zsoil[k];
        }
        else
        {
            ddz = -zsoil[0];
        }

        sh2omid[k] = ws->sh2o[k] + ci[k] + wplus / ddz;
        stot = sh2omid[k] + sice[k];

        if (stot > soil->smcmax)
        {
            ws->smc[k] = soil->smcmax;
            wplus = (stot - soil->smcmax) * ddz;
        }
        else
        {
            if (k > ps->nwtbl - 1)
            {
                ws->smc[k] = soil->smcmax;
                wplus = (stot - soil->smcmax) * ddz;
            }
            else
            {
                ws->smc[k] = stot;
                wplus = 0.0;
            }
        }

        sh2omid[k] = ws->smc[k] - sice[k];
    }

    ddz = -zsoil[0];
    for (k = 0; k < ps->nsoil; k++)
    {
        if (k != 0)
        {
            ddz = zsoil[k - 1] - zsoil[k];
        }
        ws->sh2o[k] = sh2omid[k] + wplus / ddz;
        stot = ws->sh2o[k] + sice[k];
        if (stot > soil->smcmax)
        {
            wplus = (stot - soil->smcmax) * ddz;
        }
        else
        {
            wplus = 0.0;
        }

        ws->smc[k] = (stot < soil->smcmax) ? stot : soil->smcmax;
        ws->smc[k] = (ws->smc[k] > soil->smcmin + 0.02) ?
            ws->smc[k] : soil->smcmin + 0.02;
        ws->sh2o[k] = ws->smc[k] - sice[k];
        ws->sh2o[k] = (ws->sh2o[k] > 0.0) ? ws->sh2o[k] : 0.0;
    }

    /* Update canopy water content/interception (cmc). Convert rhsct to an
     * "amount" value and add to previous cmc value to get new cmc. */
    wf->runoff3 = wplus / dt;
    ws->cmc += dt * rhsct;
    if (ws->cmc < 1.0e-20)
    {
        ws->cmc = 0.0;
    }
    ws->cmc = (ws->cmc < ws->cmcmax) ? ws->cmc : ws->cmcmax;
}

double TBnd (double tu, double tb, const double *zsoil, double zbot, int k,
    int nsoil)
{
    /*
     * Subroutine TBnd
     *
     * Calculate temperature on the boundary of the layer by interpolation of
     * the middle layer temperatures */
    double          zb, zup;
    double          tbnd1;

    /* Use surface temperature on the top of the first layer */
    if (k == 0)
    {
        zup = 0.0;
    }
    else
    {
        zup = zsoil[k - 1];
    }

    /* Use depth of the constant bottom temperature when interpolate
     * temperature into the last layer boundary */
    if (k == nsoil - 1)
    {
        zb = 2.0 * zbot - zsoil[k];
    }
    else
    {
        zb = zsoil[k + 1];
    }

    /* Linear interpolation between the average layer temperatures */
    tbnd1 = tu + (tb - tu) * (zup - zsoil[k]) / (zup - zb);

    return (tbnd1);
}
//
double TDfCnd (double smc, double qz, double smcmax, double smcmin,
    double sh2o)
{

    /*
     * Function TDfCnd
     *
     * Calculate thermal diffusivity and conductivity of the soil for a given
     * point and time.
     * Peters-Lidard approach (Peters-Lidard et al., 1998)
     * June 2001 changes: frozen soil condition.
     */
    double          df;
    double          ake;
    double          gammd;
    double          thkdry;
    const double    THKICE = 2.2;
    double          thko;
    const double    THKQTZ = 7.7;
    double          thksat;
    double          thks;
    const double    THKW = 0.57;
    double          satratio;
    double          xu;
    double          xunfroz;

    /* We now get quartz as an input argument
     *
     * If the soil has any moisture content compute a partial sum/product
     * otherwise use a constant value which works well with most soils
     *
     *  thkw ......water thermal conductivity
     *  thkqtz ....thermal conductivity for quartz
     *  thko ......thermal conductivity for other soil components
     *  thks ......thermal conductivity for the solids combined(quartz+other)
     *  THKICE ....ice thermal conductivity
     *  smcmax ....porosity (= smcmax)
     *  qz .........quartz content (soil type dependent)
     *
     * Use as in Peters-lidard, 1998 (Modif. from Johansen, 1975).
     * Pablo grunmann, 08/17/98
     * Refs.:
     *  Farouki, O. T.,1986: Thermal properties of soils. Series on rock and
     *      soil mechanics, Vol. 11, trans tech, 136 pp.
     *  Johansen, O., 1975: Thermal conductivity of soils. Ph.D. thesis,
     *      University of trondheim
     *  Peters-Lidard, C. D., et al., 1998: The effect of soil thermal
     *      conductivity parameterization on surface energy fluxes and
     *      temperatures. Journal of the Atmospheric Sciences, Vol. 55,
     *      pp. 1209-1224.
     */
    satratio = (smc - smcmin) / (smcmax - smcmin);

    /* Thermal conductivity of "other" soil components */
    //if (qz < 0.2) thko = 3.0;
    thko = 2.0;

    /* Solids' conductivity */
    thks = pow (THKQTZ, qz) * pow (thko, 1.0 - qz);

    /* Unfrozen fraction (from 1, i.e., 100% liquid, to 0 (100% frozen)) */
    xunfroz = sh2o / smc;

    /* Unfrozen volume for saturation (porosity * xunfroz) */
    xu = xunfroz * smcmax;

    /* Saturated thermal conductivity */
    thksat =
        pow (thks, 1.0 - smcmax) * pow (THKICE, smcmax - xu) * pow (THKW, xu);

    /* Dry density in kg/m3 */
    gammd = (1.0 - smcmax) * 2700.0;

    /* Dry thermal conductivity in W m-1 K-1 */
    thkdry = (0.135 * gammd + 64.7) / (2700.0 - 0.947 * gammd);

    if ((sh2o + 0.0005) < smc)
    {
        /* Frozen */
        ake = satratio;
    }
    else
    {
        /* Unfrozen
         * range of validity for the kersten number (ake)
         * Kersten number (using "fine" formula, valid for soils containing at
         * least 5% of particles with diameter less than 2.e-6 meters.)
         * (for "coarse" formula, see Peters-Lidard et al., 1998). */
        if (satratio > 0.1)
        {
            ake = log10 (satratio) + 1.0;
        }
        else
        {
            /* use k = kdry */
            ake = 0.0;
        }
    }

    /* Thermal conductivity */

    df = ake * (thksat - thkdry) + thkdry;

    return (df);
}

double TmpAvg (double tup, double tm, double tdn, const double *zsoil,
    int nsoil, int k)
{
    /*
     * Function TmpAvg
     *
     * Calculate soil layer average temperature (tavg) in freezing/thawing
     * layer using up, down, and middle layer temperatures (tup, tdn, tm),
     * where tup is at top boundary of layer, tdn is at bottom boundary of
     * layer. tm is layer prognostic state temperature.
     */
    double          dz;
    double          dzh;
    double          x0;
    double          xdn;
    double          xup;
    double          tavg;

    if (k == 0)
    {
        dz = - zsoil[0];
    }
    else
    {
        dz = zsoil[k - 1] - zsoil[k];
    }

    dzh = dz * 0.5;

    if (tup < TFREEZ)
    {
        if (tm < TFREEZ)
        {
            if (tdn < TFREEZ)
            {
                /* tup, tm, tdn < TFREEZ */
                tavg = (tup + 2.0 * tm + tdn) / 4.0;
            }
            else
            {
                /* tup & tm < TFREEZ, tdn > TFREEZ */
                x0 = (TFREEZ - tm) * dzh / (tdn - tm);
                tavg =
                    0.5 * (tup * dzh + tm * (dzh + x0) + TFREEZ * (2.0 * dzh -
                    x0)) / dz;
            }
        }
        else
        {
            if (tdn < TFREEZ)
            {
                /* tup < TFREEZ, tm > TFREEZ, tdn < TFREEZ */
                xup = (TFREEZ - tup) * dzh / (tm - tup);
                xdn = dzh - (TFREEZ - tm) * dzh / (tdn - tm);
                tavg =
                    0.5 * (tup * xup + TFREEZ * (2.0 * dz - xup - xdn) +
                    tdn * xdn) / dz;
            }
            else
            {
                /* tup < TFREEZ, tm > TFREEZ, tdn > TFREEZ */
                xup = (TFREEZ - tup) * dzh / (tm - tup);
                tavg = 0.5 * (tup * xup + TFREEZ * (2.0 * dz - xup)) / dz;
            }
        }
    }
    else
    {
        if (tm < TFREEZ)
        {
            if (tdn < TFREEZ)
            {
                /* tup > TFREEZ, tm < TFREEZ, tdn < TFREEZ */
                xup = dzh - (TFREEZ - tup) * dzh / (tm - tup);
                tavg = 0.5 * (TFREEZ * (dz - xup) + tm * (dzh + xup) +
                    tdn * dzh) / dz;
            }
            else
            {
                /* tup > TFREEZ, tm < TFREEZ, tdn > TFREEZ */
                xup = dzh - (TFREEZ - tup) * dzh / (tm - tup);
                xdn = (TFREEZ - tm) * dzh / (tdn - tm);
                tavg = 0.5 * (TFREEZ * (2.0 * dz - xup - xdn) + tm * (xup +
                    xdn)) / dz;
            }
        }
        else
        {
            if (tdn < TFREEZ)
            {
                /* tup > TFREEZ, tm > TFREEZ, tdn < TFREEZ */
                xdn = dzh - (TFREEZ - tm) * dzh / (tdn - tm);
                tavg = (TFREEZ * (dz - xdn) + 0.5 * (TFREEZ + tdn) * xdn) / dz;
            }
            else
            {
                /* tup > TFREEZ, tm > TFREEZ, tdn > TFREEZ */
                tavg = (tup + 2.0 * tm + tdn) / 4.0;
            }
        }
    }

    return (tavg);
}

void Transp (const ws_struct *ws, wf_struct *wf, const ps_struct *ps,
    const lc_struct *lc, const soil_struct *soil, const double *zsoil)
{
    /*
     * Function Transp
     *
     * Calculate transpiration for the veg class.
     */
    int             i, k;
    double          denom;
    double          etpa;
    double          gx[MAXLYR];
    double          rtx, sgx;

    /* Initialize plant transp to zero for all soil layers. */
    for (k = 0; k < ps->nsoil; k++)
    {
        wf->et[k] = 0.0;
    }

    /* Calculate an 'adjusted' potential transpiration
     * If statement below to avoid tangent linear problems near zero
     * Note: gx and other terms below redistribute transpiration by layer,
     * et(k), as a function of soil moisture availability, while preserving
     * total etpa. */
    if (ws->cmc != 0.0)
    {
        etpa = lc->shdfac * ps->pc * wf->etp *
            (1.0 - pow (ws->cmc / ws->cmcmax, lc->cfactr));
    }
    else
    {
        etpa = lc->shdfac * ps->pc * wf->etp;
    }

    sgx = 0.0;
    for (i = 0; i < lc->nroot; i++)
    {
        gx[i] = (ws->smc[i] - soil->smcwlt) / (soil->smcref - soil->smcwlt);
        gx[i] = (gx[i] < 0.0) ? 0.0 : gx[i];
        gx[i] = (gx[i] > 1.0) ? 1.0 : gx[i];
        sgx += gx[i];
    }
    sgx = sgx / (double)lc->nroot;

    denom = 0.0;
    for (i = 0; i < lc->nroot; i++)
    {
        rtx = lc->rtdis[i] + gx[i] - sgx;
        gx[i] *= (rtx > 0.0 ? rtx : 0.0);
        denom += gx[i];
    }
    denom = (denom <= 0.0) ? 1.0 : denom;

    for (i = 0; i < lc->nroot; i++)
    {
        wf->et[i] = etpa * gx[i] / denom;

        /* Above code assumes a vertically uniform root distribution
         * Code below tests a variable root distribution */
        //wf->et[0] = (zsoil[0] / zsoil[ps->nroot - 1]) * gx * etpa;
        //wf->et[0] = (zsoil[0] / zsoil[ps->nroot - 1]) * etpa;
        ///* Using root distribution as weighting factor */
        //wf->et[0] = (lc->rtdis[0] * etpa);
        //wf->et[0] = etpa * part[0];
        ///* Loop down thru the soil layers repeating the operation above,
        // * but using the thickness of the soil layer (rather than the
        // * absolute depth of each layer) in the final calculation. */
        //for (k = 0; k < lc->nroot; k++)
        //{ 
        //    gx = (ws->smc[k] - soil->smcwlt ) / (soil->smcref - soil->smcwlt);
        //    gx = (gx < 0.0) ? 0.0 : gx;
        //    gx = (gx > 1.0) ? 1.0 : gx;
        //    /*test canopy resistance */
        //    gx = 1.0;
        //    wf->et[k] = ((zsoil[k] - zsoil[k-1]) /
        //    zsoil[lc->nroot - 1]) * gx * etpa;
        //    wf->et[k] = ((zsoil[k] - zsoil[k-1]) /
        //    zsoil[lc->nroot - 1]) * etpa;
        //    /*using root distribution as weighting factor */
        //    wf->et[k] = lc->rtdis[k] * etpa;
        //    wf->et[k] = etpa * part[k];
        //}
    }
}

void WDfCnd (double *wdf, double *wcnd, double smc, double sicemax,
    double dsmdz, int macpore, const soil_struct *soil, const ps_struct *ps)
{
    /*
     * Function WDfCnd
     *
     * Calculate soil water diffusivity and soil hydraulic conductivity.
     */
    double          expon;
    double          factr1;
    double          factr2;
    double          vkwgt;
    double          satkfunc;
    double          dpsidsm;

    /* Calc the ratio of the actual to the max psbl soil h2o content */
    factr1 = 0.05 / (soil->smcmax - soil->smcmin);
    factr2 = (smc - soil->smcmin) / (soil->smcmax - soil->smcmin);

    /* factr2 should avoid to be 0 or 1 */
    factr2 = (factr2 > 1.0 - 5.0e-4) ? 1.0 - 5.0e-4 : factr2;
    factr2 = (factr2 < 0.0 + 5.0e-4) ? 5.0e-4 : factr2;

    factr1 = (factr1 < factr2) ? factr1 : factr2;
    expon = 1.0 - 1.0 / soil->beta;

    satkfunc = KrFunc (soil->alpha, soil->beta, factr2);
    dpsidsm =
        (1.0 - expon) / soil->alpha / expon / (soil->smcmax - soil->smcmin) *
        pow (pow (factr2, - 1.0 / expon) - 1.0, - expon) *
        pow (factr2, - (1.0 / expon + 1.0));

    if (macpore == 1)
    {
        *wcnd =
            EFFKV (satkfunc, factr2, ps->macpore_status, soil->kmacv, soil->ksatv,
            soil->areafh);
    }
    else
    {
        *wcnd = soil->ksatv * satkfunc;
    }

    *wdf = *wcnd * dpsidsm;

    if (sicemax > 0.0)
    {
        vkwgt = 1.0 / (1.0 + pow (500.0 * sicemax, 3.0));
        satkfunc = KrFunc (soil->alpha, soil->beta, factr1);
        dpsidsm =
            (1.0 - expon) / soil->alpha / expon
            / (soil->smcmax - soil->smcmin) *
            pow (pow (factr1, - 1.0 / expon) - 1.0, - expon) *
            pow (factr1, - (1.0 / expon + 1.0));
        if (macpore == 1)
        {
            *wdf =
                vkwgt * *wdf + (1.0 - vkwgt) * dpsidsm * EFFKV (satkfunc,
                factr1, ps->macpore_status, soil->kmacv, soil->ksatv,
                soil->areafh);
        }
        else
        {
            *wdf = vkwgt * *wdf + (1.0 - vkwgt) * dpsidsm * satkfunc *
                soil->ksatv;
        }
    }
}

void SfcDifOff (ps_struct *ps, const lc_struct *lc, double t1v, double th2v,
    int iz0tlnd)
{
    /*
     * Calculate surface layer exchange coefficients via iterative process.
     * See Chen et al. (1997, BLM)
     */
    double          zilfc;
    double          zu;
    double          zt;
    double          rdz;
    double          cxch;
    double          dthv;
    double          du2;
    double          btgh;
    double          wstar2;
    double          ustar;
    double          zslu;
    double          zslt;
    double          rlogu;
    double          rlogt;
    double          rlmo;
    double          zetalt;
    double          zetalu;
    double          zetau;
    double          zetat;
    double          xlu4;
    double          xlt4;
    double          xu4;
    double          xt4;

    double          xlu;
    double          xlt;
    double          xu;
    double          xt;
    double          psmz;
    double          simm;
    double          pshz;
    double          simh;
    double          ustark;
    double          rlmn;
    double          rlma;

    int             ilech;
    int             itr;

    const double    WWST = 1.2;
    double          wwst2;
    const double    VKRM = 0.40;
    const double    EXCM = 0.001;
    const double    BETA = 1.0 / 270.0;
    double          btg;
    double          elfc;
    double          WOLD = 0.15;
    double          wnew;
    const int       ITRMX = 5;

    const double    EPSU2 = 1.0e-4;
    const double    EPSUST = 0.07;
    const double    ZTMIN = -5.0;
    const double    ZTMAX = 1.0;
    const double    HPBL = 1000.0;
    const double    SQVISC = 258.2;

    wwst2 = WWST * WWST;
    btg = BETA * GRAV;
    elfc = VKRM * btg;
    wnew = 1.0 - WOLD;

    /* czil: constant C in Zilitinkevich, S. S.1995 */

    ilech = 0;

    if ((iz0tlnd == 0) || lc->isurban)
    {
        /* Just use the original Czil value. */
        zilfc = -ps->czil * VKRM * SQVISC;
    }
    else
    {
        /* Modify czil according to Chen & Zhang, 2009
         * czil = 10 ** -0.40 h, ( where h = 10*zo ) */
        ps->czil = pow (10.0, -0.4 * (ps->z0 / 0.07));
        zilfc = -ps->czil * VKRM * SQVISC;
    }

    zu = ps->z0;
    rdz = 1.0 / ps->zlvl_wind;
    cxch = EXCM * rdz;
    dthv = th2v - t1v;

    /* Beljars correction of ustar */
    du2 = (ps->sfcspd * ps->sfcspd > EPSU2) ? ps->sfcspd * ps->sfcspd : EPSU2;

    btgh = btg * HPBL;
    /* If statements to avoid tangent linear problems near zero */
    if (btgh * ps->ch * dthv != 0.0)
    {
        wstar2 = wwst2 * pow (fabs (btgh * ps->ch * dthv), 2.0 / 3.0);
    }
    else
    {
        wstar2 = 0.0;
    }

    ustar = sqrt (ps->cm * sqrt (du2 + wstar2));
    ustar = (ustar > EPSUST) ? ustar : EPSUST;

    /* zilitinkevitch approach for zt */
    zt = exp (zilfc * sqrt (ustar * ps->z0)) * ps->z0;
    zslu = ps->zlvl_wind + zu;

    zslt = ps->zlvl + zt;

    rlogu = log (zslu / zu);

    rlogt = log (zslt / zt);

    rlmo = elfc * ps->ch * dthv / pow (ustar, 3);

    for (itr = 0; itr < ITRMX; itr++)
    {
        zetalt = zslt * rlmo;
        zetalt = zetalt > ZTMIN ? zetalt : ZTMIN;
        rlmo = zetalt / zslt;
        zetalu = zslu * rlmo;
        zetau = zu * rlmo;

        zetat = zt * rlmo;

        /* 1./Monin-Obukkhov length-scale */
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

                psmz = Pspmu (xu);
                simm = Pspmu (xlu) - psmz + rlogu;
                pshz = Psphu (xt);
                simh = Psphu (xlt) - pshz + rlogt;
            }
            else
            {
                zetalu = zetalu < ZTMAX ? zetalu : ZTMAX;
                zetalt = zetalt < ZTMAX ? zetalt : ZTMAX;

                psmz = Pspms (zetau);
                simm = Pspms (zetalu) - psmz + rlogu;
                pshz = Psphs (zetat);
                simh = Psphs (zetalt) - pshz + rlogt;
            }
        }
        /* Lech's functions */
        else
        {
            if (rlmo < 0.0)
            {
                psmz = Pslmu (zetau);
                simm = Pslmu (zetalu) - psmz + rlogu;
                pshz = Pslhu (zetat);
                simh = Pslhu (zetalt) - pshz + rlogt;
            }
            else
            {
                zetalu = zetalu < ZTMAX ? zetalu : ZTMAX;
                zetalt = zetalt < ZTMAX ? zetalt : ZTMAX;

                psmz = Pslms (zetau);
                simm = Pslms (zetalu) - psmz + rlogu;
                pshz = Pslhs (zetat);
                simh = Pslhs (zetalt) - pshz + rlogt;
            }
        }

        /* Beljaars correction for ustar */
        ustar = sqrt (ps->cm * sqrt (du2 + wstar2));
        ustar = ustar > EPSUST ? ustar : EPSUST;

        /* Zilitinkevitch fix for zt */
        zt = exp (zilfc * sqrt (ustar * ps->z0)) * ps->z0;
        zslt = ps->zlvl + zt;

        rlogt = log (zslt / zt);
        ustark = ustar * VKRM;
        ps->cm = (ustark / simm) > cxch ? ustark / simm : cxch;

        ps->ch = (ustark / simh) > cxch ? ustark / simh : cxch;

        /* if statements to avoid tangent linear problems near zero */
        if (btgh * ps->ch * dthv != 0.0)
        {
            wstar2 = wwst2 * pow (fabs (btgh * ps->ch * dthv), 2.0 / 3.0);
        }
        else
        {
            wstar2 = 0.0;
        }

        rlmn = elfc * ps->ch * dthv / pow (ustar, 3.0);

        rlma = rlmo * WOLD + rlmn * wnew;

        rlmo = rlma;
    }
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
    printf ("now: %lf, before: %lf\n", x,
        zz * rfac - 2.076 * (1. - 1. / (zz + 1.)));
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

double EFFKV (double ksatfunc, double elemsatn, int status, double mackv,
    double kv, double areaf)
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
                return (mackv * areaf * ksatfunc + kv * (1. -
                        areaf) * ksatfunc);
            else
                return (mackv * areaf + kv * (1. - areaf) * ksatfunc);
        }
    }
}
