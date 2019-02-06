#include "pihm.h"

void SFlxGlacial(wstate_struct *ws, wflux_struct *wf, estate_struct *es,
    eflux_struct *ef, pstate_struct *ps, lc_struct *lc, epconst_struct *epc,
    soil_struct *soil, double dt)
{
    /*
     * Sub-driver for "Noah LSM" family of physics subroutines for a
     * soil/veg/snowpack land-surface model to update ice temperature, skin
     * temperature, snowpack water content, snowdepth, and all terms of the
     * surface energy balance (excluding input atmospheric forcings of downward
     * radiation and precip
     */
    int             frzgra, snowng;
    const int       IZ0TLND = 0;
    double          df1;
    double          df1a;
    double          dsoil;
    double          dtot;
    double          frcsno, frcsoi;
    double          t1v;
    double          th2v;
    double          t2v;
    double          t24;
    double          interp_fraction;
    double          sn_new;
    double          prcpf;
    double          soilwm;
    double          soilww;
    double          smav[MAXLYR];
    int             k;

    /*
     * Initialization
     */
    wf->snomlt = 0.0;
    wf->pcpdrp = 0.0;

    /* Flux-PIHM uses LAI as a forcing variable.
     * Vegetation fraction is calculated from LAI following Noah-MP */
    if (ps->proj_lai >= lc->laimax)
    {
        ps->embrd = lc->emissmax;
        ps->alb = lc->albedomin;
        ps->z0brd = lc->z0max;
    }
    else if (ps->proj_lai <= lc->laimin)
    {
        ps->embrd = lc->emissmin;
        ps->alb = lc->albedomax;
        ps->z0brd = lc->z0min;
    }
    else
    {
        if (lc->laimax > lc->laimin)
        {
            interp_fraction =
                (ps->proj_lai - lc->laimin) / (lc->laimax - lc->laimin);

            /* Bound interp_fraction between 0 and 1 */
            interp_fraction = (interp_fraction < 1.0) ? interp_fraction : 1.0;
            interp_fraction = (interp_fraction > 0.0) ? interp_fraction : 0.0;

            /* Scale emissivity and LAI between emissmin and emissmax by
             * interp_fraction */
            ps->embrd = ((1.0 - interp_fraction) * lc->emissmin) +
                (interp_fraction * lc->emissmax);
            ps->alb = ((1.0 - interp_fraction) * lc->albedomax) +
                (interp_fraction * lc->albedomin);
            ps->z0brd = ((1.0 - interp_fraction) * lc->z0min) +
                (interp_fraction * lc->z0max);
        }
        else
        {
            ps->embrd = 0.5 * lc->emissmin + 0.5 * lc->emissmax;
            ps->alb = 0.5 * lc->albedomin + 0.5 * lc->albedomax;
            ps->z0brd = 0.5 * lc->z0min + 0.5 * lc->z0max;
        }
    }

    /* Initialize precipitation logicals. */
    snowng = 0;
    frzgra = 0;

    /* If input snowpack is nonzero, then compute snow density "sndens" and
     * snow thermal conductivity "sncond" subroutine */
    if (ws->sneqv <= 1.0e-7)    /* Safer if KMH (2008/03/25) */
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
            PIHMprintf(VL_ERROR,
                "Error: Physical snow depth is less than snow water equiv.\n");
            PIHMexit(EXIT_FAILURE);
        }

        ps->sncond = CSnow(ps->sndens);
    }

    /* Determine if it's precipitating and what kind of precip it is.
     * If it's prcping and the air temp is colder than 0 C, it's snowing!
     * If it's prcping and the air temp is warmer than 0 C, but the grnd temp is
     * colder than 0 C, freezing rain is presumed to be falling. */
    if (wf->prcp > 0.0)
    {
        /* Snow defined when fraction of frozen precip (ffrozp) > 0.5, passed in
         * from model microphysics.  */
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

    /* If either prcp flag is set, determine new snowfall and add it to the
     * existing snowpack.
     * Note that since all precip is added to snowpack, no precip infiltrates
     * into the soil so that prcpf is set to zero. */
    if (snowng || frzgra)
    {
        sn_new = wf->prcp * dt;
        ws->sneqv += sn_new;
        prcpf = 0.0;

        /* Update snow density based on new snowfall, using old and new snow.
         * Update snow thermal conductivity */
        SnowNew(es, sn_new, ps);
        ps->sncond = CSnow(ps->sndens);
    }
    else
    {
        /* Precip is liquid (rain), hence save in the precip variable (along
         * with any canopy "drip" added to this later) */
        prcpf = wf->prcp;
    }

    /*
     * Determine snowcover and albedo over land.
     */
    if (ws->sneqv == 0.0)
    {
        /* If snow depth = 0, set snow fraction = 0, albedo = snow free albedo.
         */
        ps->sncovr = 0.0;
        ps->albedo = ps->alb;
        ps->emissi = ps->embrd;
    }
    else
    {
        /* Determine snow fractional coverage.
         * Determine surface albedo modification due to snowdepth state. */
        ps->sncovr = SnFrac(ws->sneqv, lc->snup, ps->salp);
        ps->sncovr = (ps->sncovr < 0.98) ? ps->sncovr : 0.98;

        AlCalc(ps, dt, snowng);
    }

    ps->icecond = 2.4;

    /* Thermal conductivity */
    df1 = (ps->snowh * ps->sncond + ps->iceh * ps->icecond) /
        (ps->snowh + ps->iceh);

    /* Finally "plane parallel" snowpack effect following V. J. Linardini
     * reference cited above. Note that dtot is combined depth of snowdepth and
     * thickness of first soil layer */
    dsoil = -(0.5 * ps->zsoil[0]);

    dtot = ps->snowh + ps->iceh + dsoil;
    dtot = (dtot > 2.0 * dsoil) ? 2.0 * dsoil : dtot;

    /* Calculate subsurface heat flux, ssoil, from final thermal
     * diffusivity of surface mediums, df1 above, and skin temperature and
     * top mid-layer soil temperature */
    ef->ssoil = df1 * (es->t1 - es->stc[0]) / dtot;

    /*
     * Determine surface roughness over snowpack using snow condition from the
     * previous timestep.
     */
    ps->z0 = Snowz0(ps->sncovr, ps->z0brd, ps->snowh + ps->iceh);

    /*
     * Next call function SfcDif to calculate the sfc exchange coef (ch) for
     * heat and moisture.
     *
     * Note !!!
     * Do not call SfcDif until after above call to RedPrm, in case alternative
     * values of roughness length (z0) and Zilintinkevich coef (czil) are set
     * there via namelist i/o.
     *
     * Note !!!
     * Function SfcDif returns a ch that represents the wind spd times the
     * "original" nondimensional "ch" typical in literature. Hence the ch
     * returned from SfcDif has units of m/s. The important companion
     * coefficient of ch, carried here as "rch", is the ch from sfcdif times air
     * density and parameter "CP". "rch" is computed in "Penman".
     * rch rather than ch is the coeff usually invoked later in eqns.
     *
     * Note !!!
     * SfcDif also returns the surface exchange coefficient for momentum, cm,
     * also known as the surface drag coefficient. Needed as a state variable
     * for iterative/implicit solution of ch in SfcDif
     */
    t1v = es->t1 * (1.0 + 0.61 * ps->q2);
    th2v = es->th2 * (1.0 + 0.61 * ps->q2);

    SfcDifOff(ps, lc, t1v, th2v, IZ0TLND);

    /*
     * Call Penman function to calculate potential evaporation (ETP), and other
     * partial products and sums save in common/rite for later calculations.
     */

    /* Calculate total downward radiation (solar plus longwave) needed in
     * Penman ep subroutine that follows */
    ef->fdown = ef->solnet + ef->lwdn;

    /* Calc virtual temps and virtual potential temps needed by Penman. */
    t2v = es->sfctmp * (1.0 + 0.61 * ps->q2);

    PenmanGlacial(wf, es, ef, ps, &t24, t2v, snowng, frzgra);

    ps->rc = 0.0;

    IcePac(ws, wf, es, ef, ps, lc, soil, snowng, dt, t24, prcpf, df1);
    ps->eta_kinematic = wf->esnow * 1000.0;

    /* Calculate effective mixing ratio at grnd level (skin) */
    ps->q1 = ps->q2 + ps->eta_kinematic * CP / ps->rch;

    /* Determine sensible heat (H) in energy units (W m-2) */
    ef->sheat = -(ps->ch * CP * ps->sfcprs) / (RD * t2v) * (es->th2 - es->t1);

    /* Convert evap terms from rate (m s-1) to energy units (w m-2) */
    ef->edir = wf->edir * 1000.0 * LVH2O;
    ef->ec = wf->ec * 1000.0 * LVH2O;
    for (k = 0; k < ps->nsoil; k++)
    {
        ef->et[k] = wf->et[k] * 1000.0 * LVH2O;
    }
    ef->ett = wf->ett * 1000.0 * LVH2O;
    ef->esnow = wf->esnow * 1000.0 * LSUBS;
    ef->etp = wf->etp * 1000.0 * LSUBS;
    ef->eta = (ef->etp > 0.0) ? ef->esnow : ef->etp;

    /* Determine beta (ratio of actual to potential evap) */
    ps->beta = (ef->etp == 0.0) ? 0.0 : (ef->eta / ef->etp);

    /* Convert the sign of soil heat flux so that:
     *   ssoil>0: warm the surface  (night time)
     *   ssoil<0: cool the surface  (day time) */
    ef->ssoil *= -1.0;

    ws->soilm = -1.0 * ws->smc[0] * ps->zsoil[0];
    for (k = 1; k < ps->nsoil; k++)
    {
        ws->soilm += ws->smc[k] * (ps->zsoil[k - 1] - ps->zsoil[k]);
    }

    ps->soilw = 0.0;
}

void PenmanGlacial(wflux_struct *wf, const estate_struct *es, eflux_struct *ef,
    pstate_struct *ps, double *t24, double t2v, int snowng, int frzgra)
{
    /*
     * Function Penman
     *
     * Calculate potential evaporation for the current point. Various partial
     * sums/products are also calculated and passed back to the calling routine
     * for later use.
     */
    double          a;
    double          delta;
    double          fnet;
    double          rad;
    double          rho;
    double          elcp1;
    double          lvs;
    const double    ELCP = 2.4888e+3;
    const double    LSUBC = 2.501000e+6;

    /* Prepare partial quantities for Penman equation. */
    elcp1 = (es->t1 > TFREEZ) ? ELCP : ELCP * LSUBS / LSUBC;
    lvs = (es->t1 > TFREEZ) ? LSUBC : LSUBS;

    delta = elcp1 * ps->dqsdt2;
    a = elcp1 * (ps->q2sat - ps->q2);
    *t24 = es->sfctmp * es->sfctmp * es->sfctmp * es->sfctmp;
    ps->rr = ps->emissi * *t24 * 6.48e-8 / (ps->sfcprs * ps->ch) + 1.0;

    rho = ps->sfcprs / (RD * t2v);
    ps->rch = rho * CP * ps->ch;

    /* Adjust the partial sums/products with the latent heat effects caused by
     * falling precipitation. */
    if (!snowng)
    {
        if (wf->prcp > 0.0)
        {
            ps->rr += CPH2O * wf->prcp * 1000.0 / ps->rch;
        }
    }
    else
    {
        ps->rr += CPICE * wf->prcp * 1000.0 / ps->rch;
    }

    /* Include the latent heat effects of frzng rain converting to ice on impact
     * in the calculation of flx2 and fnet. */
    ef->flx2 = (frzgra) ? -LSUBF * wf->prcp * 1000.0 : 0.0;

    fnet = ef->fdown - ps->emissi * SIGMA * *t24 - ef->ssoil - ef->flx2;

    /* Finish Penman equation calculations */
    rad = fnet / ps->rch + es->th2 - es->sfctmp;
    ps->epsca = (a * ps->rr + rad * delta) / (delta + ps->rr);
    wf->etp = ps->epsca * ps->rch / lvs / 1000.0;
}

void IcePac(wstate_struct *ws, wflux_struct *wf, estate_struct *es,
    eflux_struct *ef, pstate_struct *ps, const lc_struct *lc,
    const soil_struct *soil, int snowng, double dt, double t24, double prcpf,
    double df1)
{
    /*
     * Function IcePac
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
    double          ssoil1;
    double          t11;
    double          yy;
    double          zz1;
    double          sniceeqv;
    const double    ESDMIN = 1.0e-6;
    const double    SNOEXP = 2.0;

    /* Initialize evap terms. */
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

    ps->beta = 1.0;

    /* If etp < 0 (downward) then dewfall (= frostfall in this case). */
    if (wf->etp <= 0.0)
    {
        if ((ps->ribb >= 0.1) && (ef->fdown > 150.0))
        {
            wf->etp = (((wf->etp * (1.0 - ps->ribb) < 0.0) ?
                wf->etp * (1.0 - ps->ribb) : 0.0) / 0.980 +
                wf->etp * (0.980 - 1.0)) / 0.980;
        }

        if (wf->etp == 0.0)
        {
            ps->beta = 0.0;
        }

        wf->dew = -wf->etp;
        esnow2 = wf->etp * dt;
        etanrg = wf->etp * 1000.0 * LSUBS;
    }
    else
    {
        wf->esnow = wf->etp;
        esnow2 = wf->esnow * dt;
        etanrg = wf->esnow * 1000.0 * LSUBS;
    }

    /* If precip is falling, calculate heat flux from snow sfc to newly
     * accumulating precip.  note that this reflects the flux appropriate for
     * the not-yet-updated skin temperature (t1).  assumes temperature of the
     * snowfall striking the ground is = sfctmp (lowest model level air temp) */
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
     * fluxes. flx1 from above, flx2 brought in via common block rite.
     * flx2 reflects freezing rain latent heat flux using t1 calculated in
     * Penman. */
    dsoil = -0.5 * ps->zsoil[0];
    dtot = ps->snowh + ps->iceh + dsoil;
    denom = 1.0 + df1 / (dtot * ps->rr * ps->rch);

    t12a = ((ef->fdown - ef->flx1 - ef->flx2 - ps->emissi * SIGMA * t24) /
        ps->rch + es->th2 - es->sfctmp - etanrg / ps->rch) / ps->rr;
    t12b = df1 * es->stc[0] / (dtot * ps->rr * ps->rch);

    t12 = (es->sfctmp + t12a + t12b) / denom;

    /* If the 'effective snow-grnd sfc temp' is at or below freezing, no snow
     * melt will occur. Set the skin temp to this effective temp. Reduce (by
     * sublimation) or increase (by frost) the depth of the snowpack, depending
     * on sign of etp.
     * Update soil heat flux (ssoil) using new skin temperature (t1) since no
     * snowmelt, set accumulated snowmelt to zero, set 'effective' precip from
     * snowmelt to zero, set phase-change heat flux from snowmelt to zero. */
    sniceeqv = (ws->sneqv > 0.0) ? ws->sneqv : ps->iceh * ICEDENS;

    if (t12 <= TFREEZ)
    {
        /* Sub-freezing block */
        es->t1 = t12;
        ef->ssoil = df1 * (es->t1 - es->stc[0]) / dtot;

        /* Snow melts before ice */
        if (sniceeqv - esnow2 > 0.0)
        {
            sniceeqv -= esnow2;
        }
        else
        {
            sniceeqv = 0.0;
            esnow2 = sniceeqv;
            wf->esnow = esnow2 / dt;
        }

        ef->flx3 = 0.0;
        ex = 0.0;

        wf->snomlt = 0.0;
    }
    else
    {
        /* If the 'effective snow-grnd sfc temp' is above freezing, snow melt
         * will occur. Call the snow melt rate, ex and amt, snomlt. Revise the
         * effective snow depth. Revise the skin temp because it would have chgd
         * due to the latent heat released by the melting. Calc the latent heat
         * released, flx3.
         * Set the effective precip, prcp1 to the snow melt rate, ex for use in
         * SmFlx. Adjustment to t1 to account for snow patches.
         * Calculate qsat valid at freezing point. Note that esat (saturation
         * vapor pressure) value of 6.11e+2 used here is that valid at freezing
         * point.
         * Note that etp from call Penman in sflx is ignored here in favor of
         * bulk etp over 'open water' at freezing temp. Update soil heat flux
         * (s) using new skin temperature (t1) */
        /* Above freezing block */
        es->t1 = TFREEZ;
        ps->beta = 1.0;
        dtot = (dtot > 2.0 * dsoil) ? 2.0 * dsoil : dtot;
        ef->ssoil = df1 * (es->t1 - es->stc[0]) / dtot;

        if (sniceeqv - esnow2 <= ESDMIN)
        {
            /* If potential evap (sublimation) greater than depth of snowpack.
             * beta < 1
             * snowpack has sublimated away, set depth to zero. */
            sniceeqv = 0.0;
            ex = 0.0;
            wf->snomlt = 0.0;
            ef->flx3 = 0.0;
        }
        else
        {
            /* Sublimation less than depth of snowpack
             * Snowpack (esd) reduced by esnow2 (depth of sublimated snow) */
            sniceeqv -= esnow2;
            etp3 = wf->etp * 1000.0 * LVH2O;
            seh = ps->rch * (es->t1 - es->th2);
            t14 = es->t1 * es->t1 * es->t1 * es->t1;
            ef->flx3 =
                ef->fdown - ef->flx1 - ef->flx2 - ps->emissi * SIGMA * t14 -
                ef->ssoil - seh - etanrg;
            ef->flx3 = (ef->flx3 <= 0.0) ? 0.0 : ef->flx3;

            ex = ef->flx3 * 0.001 / LSUBF;

            /* Snowmelt reduction depending on snow cover */
            wf->snomlt = ex;

            /* ESDMIN represents a snowpack depth threshold value below which
             * we choose not to retain any snowpack, and instead include it in
             * snowmelt. */
            if (sniceeqv - wf->snomlt * dt >= ESDMIN)
            {
                sniceeqv -= wf->snomlt * dt;
            }
            else
            {
                /* Snowmelt exceeds snow depth */
                ex = sniceeqv / dt;
                ef->flx3 = ex * 1000.0 * LSUBF;
                wf->snomlt = sniceeqv / dt;

                sniceeqv = 0.0;
            }
        }

        prcpf += ex;
    }

    if (ws->sneqv > 0.0)
    {
        ws->sneqv = sniceeqv;
    }
    else
    {
        ps->iceh = sniceeqv / ICEDENS;
    }

    wf->pcpdrp =  prcpf;

    /* Before call ShFlx in this snowpack case, set zz1 and yy arguments to
     * special values that ensure that ground heat flux calculated in ShFlx
     * matches that already computer for below the snowpack, thus the sfc heat
     * flux to be computed in ShFlx will effectively be the flux at the snow top
     * surface.  t11 is a dummy argument so we will not use the skin temp value
     * as revised by ShFlx. */
    zz1 = 1.0;
    yy = es->stc[0] - 0.5 * ef->ssoil * ps->zsoil[0] * zz1 / df1;

    t11 = es->t1;
    ssoil1 = ef->ssoil;

    /* ShFlx will calc/update the soil temps. Note: the sub-sfc heat flux
     * (ssoil1) and the skin temp (t11) output from this ShFlx call are not used
     * in any subsequent calculations. Rather, they are dummy variables here in
     * the SnoPac case, since the skin temp and sub-sfc heat flux are updated
     * instead near the beginning of the call to SnoPac. */
    ShFlx(ws, es, ps, lc, soil, dt, yy, zz1, df1);

    es->t1 = t11;
    ef->ssoil = ssoil1;

    /* Snow depth and density adjustment based on snow compaction. yy is assumed
     * to be the soil temperature at the top of the soil column. */
    if (ws->sneqv > 0.0)
    {
        SnowPack(ws->sneqv, dt, &ps->snowh, &ps->sndens, es->t1, yy);
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
