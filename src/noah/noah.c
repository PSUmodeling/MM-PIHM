#include "pihm.h"

void Noah(elem_struct *elem, const lctbl_struct *lctbl, const calib_struct *cal,
    double dt)
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        if (elem[i].lc.glacier == 1 && elem[i].ps.iceh <= 0.0)
        {
            elem[i].lc.glacier = 0;
            elem[i].attrib.lc_type = IGBP_BARREN;

            _InitLc(&elem[i], lctbl, cal);
        }

        CalHum(&elem[i].ps, &elem[i].es);

        elem[i].ps.ffrozp = FrozRain(elem[i].wf.prcp, elem[i].es.sfctmp);

        elem[i].ps.alb = BADVAL;

        if (elem[i].ps.q1 == BADVAL)
        {
            elem[i].ps.q1 = elem[i].ps.q2;
        }

        elem[i].ef.solnet = elem[i].ef.soldn * (1.0 - elem[i].ps.albedo);
        elem[i].ef.lwdn = elem[i].ef.longwave * elem[i].ps.emissi;

        for (j = 0; j < elem[i].ps.nlayers; j++)
        {
            elem[i].ws.smc[j] =
                (elem[i].ws.smc[j] > elem[i].soil.smcmin + SH2OMIN) ?
                elem[i].ws.smc[j] : elem[i].soil.smcmin + SH2OMIN;
            elem[i].ws.smc[j] = (elem[i].ws.smc[j] < elem[i].soil.smcmax) ?
                elem[i].ws.smc[j] : elem[i].soil.smcmax;
            elem[i].ws.swc[j] = (elem[i].ws.swc[j] < elem[i].ws.smc[j]) ?
                elem[i].ws.swc[j] : elem[i].ws.smc[j];
        }

        /*
         * Run Noah LSM
         */
        if (elem[i].lc.glacier == 1)
        {
            SFlxGlacial(&elem[i].ws, &elem[i].wf, &elem[i].es, &elem[i].ef,
                &elem[i].ps, &elem[i].lc, &elem[i].soil, dt);
        }
        else
        {
#if defined(_CYCLES_)
            SFlx(dt, &elem[i].weather, &elem[i].cs, &elem[i].soil, &elem[i].lc,
                elem[i].crop, &elem[i].ps, &elem[i].ws, &elem[i].wf,
                &elem[i].es, &elem[i].ef);
#else
            SFlx(&elem[i].ws, &elem[i].wf, &elem[i].es, &elem[i].ef,
                &elem[i].ps, &elem[i].lc, &elem[i].epc, &elem[i].soil, dt);
#endif
        }

        /* ET: convert from W m-2 to m s-1 */
        elem[i].wf.ec = elem[i].ef.ec / LVH2O / 1000.0;
        elem[i].wf.ett = elem[i].ef.ett / LVH2O / 1000.0;
        elem[i].wf.edir = elem[i].ef.edir / LVH2O / 1000.0;
    }
}

void NoahHydrol(elem_struct *elem, double dt)
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             k;

        /* Find water table position */
        elem[i].ps.nwtbl = FindWaterTable(elem[i].ps.soil_depth, elem[i].ps.nlayers,
            elem[i].ws.gw, elem[i].ps.satdpth);

        for (k = 0; k < elem[i].ps.nlayers; k++)
        {
            elem[i].ws.smc[k] =
                (elem[i].ws.smc[k] > elem[i].soil.smcmin + SH2OMIN) ?
                elem[i].ws.smc[k] : elem[i].soil.smcmin + SH2OMIN;
            elem[i].ws.smc[k] =
                (elem[i].ws.smc[k] < elem[i].soil.smcmax) ?
                elem[i].ws.smc[k] : elem[i].soil.smcmax;
            elem[i].ws.swc[k] =
                (elem[i].ws.swc[k] < elem[i].ws.smc[k]) ?
                elem[i].ws.swc[k] : elem[i].ws.smc[k];
        }

        SmFlx(&elem[i].ws, &elem[i].wf, &elem[i].ps, &elem[i].soil, dt);

#if defined(_CYCLES_)
        double          wflux[MAXLYR + 1];

        /*
         * Calcluate vertical transport of solute
         */
        for (k = 0; k < MAXLYR + 1; k++)
        {
            wflux[k] = 0.0;
        }

        for (k = 0; k < elem[i].ps.nlayers; k++)
        {
            /* Note in Noah and Cycles flux k represents flux from different
             * layers
             *
             *     Flux-PIHM               Cycles
             *
             *    smflx[k - 1]            wflux[k]
             * --------|--------     --------|--------
             *         V                     V
             *                 layer k
             *      smflx[k]            wflux[k + 1]
             * --------|--------     --------|--------
             *         V                     V
             */
            wflux[k + 1] = elem[i].wf.smflx[k] * RHOH2O * dt;
        }

        SoluteTransp(KD_NO3, 0.0, wflux, elem[i].ws.smc, &elem[i].soil,
            &elem[i].ps, elem[i].ns.no3);

        SoluteTransp(KD_NH4, 0.0, wflux, elem[i].ws.smc, &elem[i].soil,
            &elem[i].ps, elem[i].ns.nh4);

# if defined(_DEBUG_)
    for (k = 0; k < elem[i].ps.nlayers; k++)
    {
        if (isnan(elem[i].ns.no3[k]))
        {
            pihm_exit(EXIT_FAILURE);
        }
    }
# endif
#endif
    }

#if TEMP_DISABLED
# if defined(_DEBUG_)
    for (i = 0; i < nelem; i++)
    {
        int             k;
        double          pihm_soilm, noah_soilm;

        pihm_soilm = ((elem[i].ws.unsat + elem[i].ws.gw > elem[i].soil.depth) ?
            elem[i].soil.depth : elem[i].ws.unsat + elem[i].ws.gw) *
            elem[i].soil.porosity + elem[i].soil.depth * elem[i].soil.smcmin;
        noah_soilm = elem[i].ws.smc[0] * elem[i].ps.soil_depth[0];
        for (k = 1; k < elem[i].ps.nlayers; k++)
        {
            noah_soilm += elem[i].ws.smc[k] * elem[i].ps.soil_depth[k];
        }
        pihm_printf(VL_NORMAL, "%d: %lf, %lf\n", i + 1, pihm_soilm, noah_soilm);
    }
# endif
#endif
}

#if defined(_CYCLES_)
void SFlx(double dt, const weather_struct *weather, const cstate_struct *cs,
    soil_struct *soil, lc_struct *lc, crop_struct crop[], phystate_struct *ps,
    wstate_struct *ws, wflux_struct *wf, estate_struct *es, eflux_struct *ef)
#else
void SFlx(wstate_struct *ws, wflux_struct *wf, estate_struct *es,
    eflux_struct *ef, phystate_struct *ps, lc_struct *lc, epconst_struct *epc,
    soil_struct *soil, double dt)
#endif
{
    /*
     * Sub-driver for "Noah LSM" family of physics subroutines for a soil/veg/
     * snowpack land-surface model to update soil moisture, soil ice, soil
     * temperature, skin temperature, snowpack water content, snowdepth, and all
     * terms of the surface energy balance, modified for Flux-PIHM
     */
    int             frzgra;             /* Flag that indicates freezing rain (-)
                                         */
    int             snowng;             /* Flag that indicates snow (-) */
    const int       IZ0TLND = 0;        /* Option to turn on (IZ0TLND=1) or off
                                         * (IZ0TLND=0) the vegetation-category-
                                         * dependent calculation of the
                                         * Zilitinkevich coefficient CZIL in the
                                         * SfcDif subroutines */
    double          df1;                /* Thermal conductivity of surface
                                         * mediums (W m-1 K-1) */
    double          df1a;               /* Thermal conductivity of surface snow
                                         * mediums (W m-1 K-1) */
    double          dsoil;              /* Soil thickness for soil heat flux
                                         * calculation (m) */
    double          dtot;               /* Soil thickness plus snow depth for
                                         * soil heat flux calculation (m) */
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

    /* Urban */
    if (lc->isurban)
    {
        lc->shdfac = 0.05;
#if !defined(_CYCLES_)
        epc->rsmin = 400.0;
#endif
        /* In original Noah LSM, urban model changes urban soil porosity. In
         * Flux-PIHM, porosity should not be changed because it is also used in
         * PIHM hydrology calculation. In addition, field capacity and wilting
         * point are not changed to allow transpiration from the non-urban
         * fraction. Parameter smcdry is adjusted according to Noah urban smcdry
         * saturation level to limit direct evaporation from non-vegetated area,
         * which includes urban area */
        soil->smcdry = (0.40 / 0.45) * soil->porosity + soil->smcmin;
    }

#if defined(_CYCLES_)
    lc->shdfac = CommRadIntcp(crop);
#endif

    /* Set minimum LAI for non-barren land cover to improve performance */
    if (lc->shdfac > 0.0)
    {
        ps->proj_lai = MAX(ps->proj_lai, 0.5);
    }

#if !defined(_CYCLES_)
    /* Calculate maximum canopy moisture capacity */
    ws->cmcmax = lc->shdfac * lc->cmcfactr * ps->proj_lai;
#endif

    /* Flux-PIHM uses LAI as a forcing variable. Vegetation fraction is constant
     * unless coupled to Cycles */
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

    /* Initialize precipitation logicals */
    snowng = 0;
    frzgra = 0;

    /* If input snowpack is nonzero, then compute snow density "sndens" and
     * snow thermal conductivity "sncond" subroutine */
    if (ws->sneqv <= 1.0E-7)
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
            pihm_printf(VL_ERROR,
                "Error: Physical snow depth is less than snow water equiv.\n");
            pihm_exit(EXIT_FAILURE);
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
        /* If snow depth=0, set snow fraction=0, albedo=snow free albedo */
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

    /*
     * Next calculate the subsurface heat flux, which first requires calculation
     * of the thermal diffusivity. Treatment of the latter follows that on Pages
     * 148-149 from "Heat transfer in cold climates", by V. J. Lunardini
     * (published in 1981 by van Nostrand Reinhold Co.) i.e. treatment of two
     * contiguous "plane parallel" mediums (namely here the first soil layer and
     * the snowpack layer, if any). This diffusivity treatment behaves well for
     * both zero and nonzero snowpack, including the limit of very thin
     * snowpack. This treatment also eliminates the need to impose an arbitrary
     * upper bound on subsurface heat flux when the snowpack becomes extremely
     * thin.
     *
     * First calculate thermal diffusivity of top soil layer, using both the
     * frozen and liquid soil moisture, following the soil thermal diffusivity
     * function of Peters-Lidard et al. (1998, JAS, Vol 55, 1209-1224), which
     * requires the specifying the quartz content of the given soil class
     *
     * Next add subsurface heat flux reduction effect from the overlying green
     * canopy, adapted from Section 2.1.2 of Peters-Lidard et al. (1997, JGR,
     * Vol 102(D4))
     */
    df1 = TDfCnd(ws->smc[0], soil->quartz, soil->smcmax, soil->smcmin,
        ws->swc[0]);

    /* Urban */
    df1 = (lc->isurban) ? 3.24 : df1;

    df1 *= exp(ps->sbeta * lc->shdfac);

    df1 = (ps->sncovr > 0.97) ? ps->sncond : df1;

    /* Finally "plane parallel" snowpack effect following V. J. Linardini
     * reference cited above. Note that dtot is combined depth of snowdepth and
     * thickness of first soil layer */
    dsoil = -(0.5 * ps->zsoil[0]);
    if (ws->sneqv == 0.0)
    {
        ef->ssoil = df1 * (es->t1 - es->stc[0]) / dsoil;
    }
    else
    {
        dtot = ps->snowh + dsoil;
        frcsno = ps->snowh / dtot;
        frcsoi = dsoil / dtot;

        df1a = frcsno * ps->sncond + frcsoi * df1;

        df1 = df1a * ps->sncovr + df1 * (1.0 - ps->sncovr);

        /* Calculate subsurface heat flux, ssoil, from final thermal diffusivity
         * of surface mediums, df1 above, and skin temperature and top mid-layer
         * soil temperature */
        ef->ssoil = df1 * (es->t1 - es->stc[0]) / dtot;
    }

    /*
     * Determine surface roughness over snowpack using snow condition from the
     * previous timestep.
     */
    ps->z0 = (ps->sncovr > 0.0) ?
        Snowz0(ps->sncovr, ps->z0brd, ps->snowh) : ps->z0brd;

    /*
     * Next call function SfcDif to calculate the sfc exchange coef (ch) for
     * heat and moisture.
     *
     * Note !!!
     * Function SfcDif returns a ch that represents the wind spd times the
     * "original" nondimensional "ch" typical in literature. Hence the ch
     * returned from SfcDif has units of m/s. The important companion
     * coefficient of ch, carried here as "rch", is the ch from sfcdif times air
     * density and parameter "CP". "rch" is computed in "Penman". rch rather
     * than ch is the coeff usually invoked later in eqns.
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
     * Call Penman function to calculate potential evaporation (etp), and other
     * partial products and sums save in common/rite for later calculations.
     */

    /* Calculate total downward radiation (solar plus longwave) needed in Penman
     * eq subroutine that follows */
    ef->fdown = ef->solnet + ef->lwdn;

    /* Calc virtual temps and virtual potential temps needed by Penman */
    t2v = es->sfctmp * (1.0 + 0.61 * ps->q2);

    Penman(wf, es, ef, ps, &t24, t2v, snowng, frzgra);

    /*
     * Call CanRes to calculate the canopy resistance and convert it into pc if
     * nonzero greenness fraction
     */
    if (lc->shdfac > 0.0)
    {
#if defined(_CYCLES_)
        CanRes(es, ps);
#else
        CanRes(ws, es, ef, ps, soil, epc);
#endif
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
#if defined(_CYCLES_)
        NoPac(soil, lc, weather, cs, dt, t24, crop, ps, ws, wf, es, ef);
#else
        NoPac(ws, wf, es, ef, ps, lc, soil, dt, t24);
#endif
        ps->eta_kinematic = wf->eta * 1000.0;
    }
    else
    {
#if defined(_CYCLES_)
        SnoPac(soil, lc, weather, cs, snowng, dt, t24, prcpf, df1, crop, ps, ws,
        wf, es, ef);
#else
        SnoPac(ws, wf, es, ef, ps, lc, soil, snowng, dt, t24, prcpf, df1);
#endif
        ps->eta_kinematic = (wf->esnow + wf->etns) * 1000.0;
    }

    /* Calculate effective mixing ratio at grnd level (skin) */
    ps->q1 = ps->q2 + ps->eta_kinematic * CP / ps->rch;

    /* Determine sensible heat (H) in energy units (W m-2) */
    ef->sheat = -(ps->ch * CP * ps->sfcprs) / (RD * t2v) * (es->th2 - es->t1);

    /* Convert evap terms from rate (m s-1) to energy units (w m-2) */
    ef->edir = wf->edir * 1000.0 * LVH2O;
    ef->ec = wf->ec * 1000.0 * LVH2O;
    for (k = 0; k < ps->nlayers; k++)
    {
        ef->et[k] = wf->et[k] * 1000.0 * LVH2O;
    }
    ef->ett = wf->ett * 1000.0 * LVH2O;
    ef->esnow = wf->esnow * 1000.0 * LSUBS;
    ef->etp =
        wf->etp * 1000.0 * ((1.0 - ps->sncovr) * LVH2O + ps->sncovr * LSUBS);
    if (ef->etp > 0.0)
    {
        ef->eta = ef->edir + ef->ec + ef->ett + ef->esnow;
    }
    else
    {
        ef->eta = ef->etp;
    }

    /* Determine beta (ratio of actual to potential evap) */
    ps->beta = (ef->etp == 0.0) ? 0.0 : (ef->eta / ef->etp);

    /* Convert the sign of soil heat flux so that:
     *   ssoil>0: warm the surface  (night time)
     *   ssoil<0: cool the surface  (day time) */
    ef->ssoil *= -1.0;

    ws->soilm = -1.0 * ws->smc[0] * ps->zsoil[0];
    for (k = 1; k < ps->nlayers; k++)
    {
        ws->soilm += ws->smc[k] * (ps->zsoil[k - 1] - ps->zsoil[k]);
    }

    soilwm = -1.0 * (soil->smcmax - soil->smcwlt) * ps->zsoil[0];
    soilww = -1.0 * (ws->smc[0] - soil->smcwlt) * ps->zsoil[0];

    for (k = 0; k < ps->nlayers; k++)
    {
        smav[k] = (ws->smc[k] - soil->smcwlt) / (soil->smcmax - soil->smcwlt);
    }

    if (ps->nroot > 1)
    {
        for (k = 1; k < ps->nroot; k++)
        {
            soilwm += (soil->smcmax - soil->smcwlt) *
                (ps->zsoil[k - 1] - ps->zsoil[k]);
            soilww += (ws->smc[k] - soil->smcwlt) *
                (ps->zsoil[k - 1] - ps->zsoil[k]);
        }
    }

    if (soilwm < 1.0E-6)
    {
        soilwm = 0.0;
        ps->soilw = 0.0;
        ws->soilm = 0.0;
    }
    else
    {
        ps->soilw = soilww / soilwm;
    }
}

void AlCalc(phystate_struct *ps, double dt, int snowng)
{
    /*
     * Calculate albedo including snow effect (0 -> 1)
     * snoalb is argument representing maximum albedo over deep snow, as passed
     * into SFlx, and adapted from the satellite-based maximum snow albedo
     * fields provided by D. Robinson and G. Kukla (1985, JCAM, Vol 24, 402-411)
     */
    double          snoalb2;
    double          snoalb1;
    const double    SNACCA = 0.94;
    const double    SNACCB = 0.58;

    ps->albedo = ps->alb + ps->sncovr * (ps->snoalb - ps->alb);
    ps->emissi = ps->embrd + ps->sncovr * (EMISSI_S - ps->embrd);

    /* Formulation by livneh
     * snoalb is considered as the maximum snow albedo for new snow, at a value
     * of 85%. Snow albedo curve defaults are from Bras P.263. should not be
     * changed except for serious problems with snow melt.
     * To implement accumulating parameters, SNACCA and SNACCB, assert that it
     * is indeed accumulation season. i.e. that snow surface temp is below zero
     * and the date falls between October and February */
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
        snoalb2 *= pow(SNACCA, pow(ps->snotime1 / 86400.0, SNACCB));
    }

    snoalb2 = (snoalb2 > ps->alb) ? snoalb2 : ps->alb;
    ps->albedo = ps->alb + ps->sncovr * (snoalb2 - ps->alb);
    ps->albedo = (ps->albedo > snoalb2) ? snoalb2 : ps->albedo;
}

void CalHum(phystate_struct *ps, estate_struct *es)
{
    const double    A2 = 17.67;
    const double    A3 = 273.15;
    const double    A4 = 29.65;
    const double    ELWV = 2.501e6;
    double          a23m4;
    const double    E0 = 611.0;
    const double    RVV = 461.0;
    const double    EPSILON = 0.622;
    double          e;
    double          esat;
    double          svp;
    const double    SVP1 = 611.2;
    const double    SVP2 = 17.67;
    const double    SVP3 = 29.65;
    const double    SVPT0 = 273.15;

    double          t2v;
    double          rho;
    double          rh;

    rh = ps->rh / 100.0;

    svp = SVP1 * exp(SVP2 * (es->sfctmp - SVPT0) / (es->sfctmp - SVP3));
    e = rh * svp;

    ps->q2 = (0.622 * e) / (ps->sfcprs - (1.0 - 0.622) * e);

    es->th2 = es->sfctmp + (0.0098 * ps->zlvl);
    t2v = es->sfctmp * (1.0 + 0.61 * ps->q2);
    rho = ps->sfcprs / (RD * t2v);

    a23m4 = A2 * (A3 - A4);

    esat = E0 * exp(ELWV / RVV * (1.0 / A3 - 1.0 / es->sfctmp));

    ps->q2sat = EPSILON * esat / (ps->sfcprs - (1.0 - EPSILON) * esat);

    ps->dqsdt2 = ps->q2sat * a23m4 / ((es->sfctmp - A4) * (es->sfctmp - A4));
}

#if defined(_CYCLES_)
void CanRes(const estate_struct *es, phystate_struct *ps)
#else
void CanRes(const wstate_struct *ws, const estate_struct *es,
    const eflux_struct *ef, phystate_struct *ps, const soil_struct *soil,
    const epconst_struct *epc)
#endif
{
    /*
     * Calculate canopy resistance which depends on incoming solar radiation,
     * air temperature, atmospheric water vapor pressure deficit at the lowest
     * model level, and soil moisture (preferably unfrozen soil moisture rather
     * than total)
     * Source:  Jarvis (1976), Noilhan and Planton (1989, MWR), Jacquemin and
     * Noilhan (1990, BLM)
     * See also: Chen et al. (1996, JGR, Vol 101(D3), 7251-7268), Eqns 12-14
     * and Table 2 of Sec. 3.1.2
     */
    double          delta;
    double          rr;
#if !defined(_CYCLES_)
    double          ff;
    double          gx;
    int             k;
    double          part[MAXLYR];
#endif
    const double    SLV = 2.501000e6;
#if defined(_CYCLES_)
    const double    RC = 70.0;              /* Cycles resistance value
                                             * 0.00081 day m-1 = 70 s m-1*/
#endif

#if defined(_CYCLES_)
    ps->rc = RC;
#else
    /* Initialize canopy resistance multiplier terms. */
    ps->rcs = 0.0;
    ps->rct = 0.0;
    ps->rcq = 0.0;
    ps->rcsoil = 0.0;

    ps->rc = 0.0;

    /* Contribution due to incoming solar radiation */
    ff = 0.55 * 2.0 * ef->soldn / (epc->rgl * ps->proj_lai);
    ps->rcs = (ff + epc->rsmin / epc->rsmax) / (1.0 + ff);
    ps->rcs = (ps->rcs > 0.0001) ? ps->rcs : 0.0001;

    /* Contribution due to air temperature at first model level above ground rct
     * expression from Noilhan and Planton (1989, MWR). */
    ps->rct = 1.0 - 0.0016 * pow(epc->topt - es->sfctmp, 2.0);
    ps->rct = (ps->rct > 0.0001) ? ps->rct : 0.0001;

    /* Contribution due to vapor pressure deficit at first model level.
     * rcq expression from ssib */
    ps->rcq = 1.0 / (1.0 + epc->hs * (ps->q2sat - ps->q2));
    ps->rcq = (ps->rcq > 0.01) ? ps->rcq : 0.01;

    /* Contribution due to soil moisture availability.
     * Determine contribution from each soil layer, then add them up. */
    gx = (ws->swc[0] - soil->smcwlt) / (soil->smcref - soil->smcwlt);
    gx = (gx > 1.0) ? 1.0 : gx;
    gx = (gx < 0.0) ? 0.0 : gx;

# if NOT_YET_IMPLEMENTED
    /* Use root distribution as weighting factor */
    part[0] = rtdis[0] * gx;
# endif
    /* Use soil depth as weighting factor */
    part[0] = (ps->zsoil[0] / ps->zsoil[ps->nroot - 1]) * gx;
    for (k = 1; k < ps->nroot; k++)
    {
        /* Frozen ground extension: total soil water "smc" was replaced by
         * unfrozen soil water "swc" */
        gx = (ws->swc[k] - soil->smcwlt) / (soil->smcref - soil->smcwlt);
        gx = (gx > 1.0) ? 1.0 : gx;
        gx = (gx < 0.0) ? 0.0 : gx;

# if NOT_YET_IMPLEMENTED
        /* Use root distribution as weighting factor */
        part[k] = rtdis[k] * gx;
# endif
        /* Use soil depth as weighting factor */
        part[k] = ((ps->zsoil[k] - ps->zsoil[k - 1]) /
            ps->zsoil[ps->nroot - 1]) * gx;
    }

    for (k = 0; k < ps->nroot; k++)
    {
        ps->rcsoil += part[k];
    }
    ps->rcsoil = (ps->rcsoil > 0.0001) ? ps->rcsoil : 0.0001;

    /* Determine canopy resistance due to all factors. Convert canopy resistance
     * (rc) to plant coefficient (pc) to be used with potential evap in
     * determining actual evap. pc is determined by:
     *   pc * linearized Penman potential evap =
     *   Penman-Monteith actual evaporation (containing rc term). */
    ps->rc = epc->rsmin /
        (ps->proj_lai * ps->rcs * ps->rct * ps->rcq * ps->rcsoil);
#endif

    rr = (4.0 * ps->emissi * SIGMA * RD / CP) * pow(es->sfctmp, 4.0) /
        (ps->sfcprs * ps->ch) + 1.0;

    delta = (SLV / CP) * ps->dqsdt2;

    ps->pc = (rr + delta) / (rr * (1.0 + ps->rc * ps->ch) + delta);
}

double CSnow(double dsnow)
{
    /*
     * Calculate snow thermal conductivity
     */
    double          c;
    double          sncond;
    const double    UNIT = 0.11631;

    /* Sncond in units of cal/(cm*hr*c), returned in W/(m*C)
     * Csnow in units of cal/(cm*hr*c), returned in W/(m*C)
     * Basic version is Dyachkova equation (1960), for range 0.1-0.4 */
    c = 0.328 * pow(10, 2.25 * dsnow);

#if NOT_YET_IMPLEMENTED
    /* De Vaux equation (1933), in range 0.1-0.6 */
    sncond = 0.0293 * (1.0 + 100.0 * dsnow * dsnow);
    csnow = 0.0293 * (1.0 + 100.0 * dsnow * dsnow);

    /* E. Andersen from Flerchinger */
    sncond = 0.021 + 2.51 * dsnow * dsnow;
    csnow = 0.021 + 2.51 * dsnow * dsnow;

    sncond = UNIT * c;
#endif

    /* Double snow thermal conductivity */
    sncond = 2.0 * UNIT * c;

    return sncond;
}

void DEvap(const wstate_struct *ws, wflux_struct *wf, const phystate_struct *ps,
    const lc_struct *lc, const soil_struct *soil)
{
    /*
     * Calculate direct soil evaporation
     */
    double          fx, sratio;

    /* Direct evap a function of relative soil moisture availability, linear
     * when fxexp = 1.
     * fx > 1 represents demand control
     * fx < 1 represents flux control */
    sratio = (ws->swc[0] - soil->smcdry) / (soil->smcmax - soil->smcdry);
    if (sratio > 0.0)
    {
        fx = pow(sratio, ps->fxexp);
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

#if defined(_CYCLES_)
void Evapo(const soil_struct *soil, const lc_struct *lc,
    const weather_struct *weather, const phystate_struct *ps,
    const estate_struct *es, const cstate_struct *cs, crop_struct crop[],
    wstate_struct *ws, wflux_struct *wf)
#else
void Evapo(const wstate_struct *ws, wflux_struct *wf, const phystate_struct *ps,
    const lc_struct *lc, const soil_struct *soil, double dt)
#endif
{
    /*
     * Calculate soil moisture flux. The soil moisture content (smc - a per unit
     * volume measurement) is a dependent variable that is updated with
     * prognostic eqns. The canopy moisture content (cmc) is also updated.
     * Frozen ground version: new states added: swc, and frozen ground
     * correction factor, frzfact and parameter slope.
     */
    int             k;
    double          cmc2ms;

    /* Executable code begins here if the potential evapotranspiration is
     * greater than zero. */
    wf->edir = 0.0;
    wf->ec = 0.0;
    wf->ett = 0.0;
    for (k = 0; k < ps->nlayers; k++)
    {
        wf->et[k] = 0.0;
    }

    if (wf->etp > 0.0)
    {
#if defined(_CYCLES_)
        /* Evaporation from residue (Cycles function) */
        SoilEvap(wf->etp * RHOH2O * DAYINSEC, ps->sncovr, crop, soil, ps, ws,
            wf);

        ResidueEvap(wf->etp * RHOH2O * DAYINSEC, ps->sncovr, crop, cs, ps, ws,
            wf);

        if (NumActiveCrop(crop) > 0)
        {
            WaterUptake(ps->pc * wf->etp * RHOH2O * DAYINSEC, soil, weather, ps,
                crop, ws, wf);
        }
#else
        if (lc->shdfac < 1.0)
        {
            /* Retrieve direct evaporation from soil surface. Call this function
             * only if veg cover not complete.
             * Frozen ground version:  swc states replace smc states. */
            DEvap(ws, wf, ps, lc, soil);
        }

        if (lc->shdfac > 0.0)
        {
            /* Initialize plant total transpiration, retrieve plant
             * transpiration, and accumulate it for all soil layers */
            Transp(ws, wf, ps, lc, soil);

            /* When coupled to Cycles, canopy evaporation is replaced by residue
             * evaporation */
            /* Calculate canopy evaporation.
             * If statements to avoid tangent linear problems near cmc = 0.0. */
            wf->ec = (ws->cmc > 0.0) ?
                (lc->shdfac * pow((ws->cmc / ws->cmcmax > 1.0) ?
                1.0 : ws->cmc / ws->cmcmax, lc->cfactr) * wf->etp) : 0.0;

            /* Ec should be limited by the total amount of available water on
             * the canopy */
            cmc2ms = ws->cmc / dt;
            wf->ec = (cmc2ms < wf->ec) ? cmc2ms : wf->ec;
        }
#endif

        for (k = 0; k < ps->nlayers; k++)
        {
            wf->ett += wf->et[k];
        }
    }

    /* Total up evap and transp types to obtain actual evapotransp */
    wf->etns = wf->edir + wf->ett + wf->ec;
}

double FrozRain(double prcp, double sfctmp)
{
    double          ffrozp;

    ffrozp = (prcp > 0.0 && sfctmp < TFREEZ) ? 1.0 : 0.0;

    return ffrozp;
}

double FrH2O(double tkelv, double smc, double swc, const soil_struct *soil)
{
    /*
     * Calculate amount of supercooled liquid soil water content if temperature
     * is below 273.15K (t0). Requires Newton-type iteration to solve the
     * nonlinear implicit equation given in Eqn 17 of Koren et al
     * (1999, JGR, Vol 104(D16), 19569-19585).
     *
     * New version (June 2001): much faster and more accurate Newton iteration
     * achieved by first taking log of eqn cited above -- less than 4
     * (typically 1 or 2) iterations achieves convergence. Also, explicit
     * 1-step solution option for special case of parameter ck = 0, which
     * reduces the original implicit equation to a simpler explicit form, known
     * as the "Flerchinger eqn". Improved handling of solution in the limit of
     * freezing point temperature t0.
     *
     * In Flux-PIHM, van Genuchten parameters are used. See Technical
     * Documentation for details
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
    const double    TOL = 0.005;
    double          mx;
    double          freew;

    nlog = 0;
    kcount = 0;

    /* If temperature not significantly below freezing (t0), swc = smc */
    if (tkelv > (TFREEZ - 1.0e-3))
    {
        freew = smc;
    }
    else
    {
        /* Option 1: iterated solution for nonzero CK in Koren et al, JGR, 1999,
         * Eqn 17 */
        if (CK != 0.0)
        {
            /* Initial guess for swl (frozen content) */
            swl = smc - swc;

            /* Keep within bounds. */
            swl = (swl > smc - SH2OMIN) ? (smc - SH2OMIN) : swl;
            swl = MAX(swl, 1.0E-4);

            /* Start of iterations */
            while (nlog < 10 && kcount == 0)
            {
                nlog++;

                satn = (smc - swl - soil->smcmin) /
                    (soil->smcmax - soil->smcmin);
                satn = MAX(satn, SATMIN);
                satn = MIN(satn, 1.0 - 1.0E-4);

                mx = soil->beta / (1.0 - soil->beta);

                df = log(GRAV / soil->alpha / LSUBF) +
                    1.0 / soil->beta * log(pow(satn, mx) - 1.0) +
                    2.0 * log(1.0 + CK * swl) - log(-(tkelv - TFREEZ) / tkelv);

                denom = 1.0 / (soil->beta - 1.0) /
                    (soil->smcmax - soil->smcmin) * pow(satn, mx - 1.0) /
                    (pow(satn, mx) - 1.0) + 2.0 * CK / (1.0 + CK * swl);

                swlk = swl - df / denom;

                /* Bounds useful for mathematical solution. */
                swlk = (swlk > smc - SH2OMIN) ? (smc - SH2OMIN) : swlk;
                swlk = (swlk < 0.0) ? 0.0 : swlk;

                /* Mathematical solution bounds applied. */
                dswl = fabs(swlk - swl);

                /* If more than 10 iterations, use explicit method (ck = 0
                 * approx.) when dswl less or eq. error, no more iterations
                 * required. */
                swl = swlk;
                if (dswl <= TOL)
                {
                    kcount++;
                }
                /* End of iterations
                 * Bounds applied within do-block are valid for physical
                 * solution */
            }

            freew = smc - swl;
        }

        /* Option 2: explicit solution for flerchinger eq. i.e. ck = 0 in Koren
         * et al., JGR, 1999, Eqn 17
         * Apply physical bounds to flerchinger solution */
        if (kcount == 0)
        {
            fk = pow(pow(-(tkelv - TFREEZ) / tkelv * soil->alpha *
                    LSUBF / GRAV, soil->beta), 1.0 / mx) *
                (soil->smcmax - soil->smcmin) - soil->smcmin;
            fk = (fk < SH2OMIN) ? SH2OMIN : fk;

            freew = (fk < smc) ? fk : smc;
        }
    }

    return freew;
}

void HRT(wstate_struct *ws, const estate_struct *es, const phystate_struct *ps,
    const lc_struct *lc, const soil_struct *soil, double *rhsts, double yy,
    double zz1, double dt, double df1, double *ai, double *bi, double *ci)
{
    /*
     * Calculate the right hand side of the time tendency term of the soil
     * thermal diffusion equation. Also to compute (prepare) the matrix
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
    double          ssoil;
    const double    CH2O = 4.2E6;       /* Volumetric heat capacity of water
                                         * (J m-3 K-1) */
    const double    CICE = 1.26E6;      /* Volumetric heat capacity of ice
                                         * (J m-3 K-1) */

    /* Urban */
    csoil_loc = (lc->isurban) ? 3.0E6 : soil->csoil;

    /* Initialize logical for soil layer temperature averaging. */
    itavg = 1;

    /* Begin section for top soil layer
     * Calc the heat capacity of the top soil layer */
    hcpct = ws->swc[0] * CH2O + (1.0 - soil->smcmax) * csoil_loc +
        (soil->smcmax - ws->smc[0]) * CP + (ws->smc[0] - ws->swc[0]) * CICE;

    /* Calc the matrix coefficients ai, bi, and ci for the top layer */
    ddz = 1.0 / (-0.5 * ps->zsoil[1]);
    ai[0] = 0.0;
    ci[0] = (df1 * ddz) / (ps->zsoil[0] * hcpct);

    /* Calculate the vertical soil temp gradient btwn the 1st and 2nd soil
     * layers. Then calculate the subsurface heat flux. use the temp gradient
     * and subsfc heat flux to calc "right-hand side tendency terms", or
     * "rhsts," or top soil layer. */
    bi[0] = -ci[0] + df1 / (0.5 * ps->zsoil[0] * ps->zsoil[0] * hcpct * zz1);
    dtsdz = (es->stc[0] - es->stc[1]) / (-0.5 * ps->zsoil[1]);
    ssoil = df1 * (es->stc[0] - yy) / (0.5 * ps->zsoil[0] * zz1);
    denom = ps->zsoil[0] * hcpct;

    rhsts[0] = (df1 * dtsdz - ssoil) / denom;

    /* Next capture the vertical difference of the heat flux at top and bottom
     * of first soil layer for use in heat flux constraint applied to potential
     * soil freezing/thawing in routine SnkSrc. */
    qtot = -1.0 * rhsts[0] * denom;

    /* Calculate frozen water content in 1st soil layer. */
    sice = ws->smc[0] - ws->swc[0];

    /* If temperature averaging invoked (itavg = true; else skip):
     * Set temp "tsurf" at top of soil column (for use in freezing soil physics
     * later in function subroutine SnkSrc). If snowpack content is zero, then
     * tsurf expression below gives tsurf = skin temp. If snowpack is nonzero
     * (hence argument zz1 = 1), then tsurf expression below yields soil column
     * top temperature under snowpack. Then calculate temperature at bottom
     * interface of 1st soil layer for use later in function subroutine SnkSrc
     */
    if (itavg)
    {
        tsurf = (yy + (zz1 - 1.0) * es->stc[0]) / zz1;

        /* If frozen water present or any of layer-1 mid-point or bounding
         * interface temperatures below freezing, then call SnkSrc to compute
         * heat source/sink (and change in frozen water content) due to possible
         * soil water phase change */
        tbk = TBnd(es->stc[0], es->stc[1], ps->zsoil, ps->zbot, 0, ps->nlayers);

        if (sice > 0.0 || es->stc[0] < TFREEZ || tsurf < TFREEZ || tbk < TFREEZ)
        {
            tavg = TmpAvg(tsurf, es->stc[0], tbk, ps->zsoil, 0);
            SnkSrc(&tsnsr, tavg, ws->smc[0], &ws->swc[0], soil, ps->zsoil, dt,
                0, qtot);
            rhsts[0] -= tsnsr / denom;
        }
    }
    else
    {
        if (sice > 0.0 || es->stc[0] < TFREEZ)
        {
            SnkSrc(&tsnsr, es->stc[0], ws->smc[0], &ws->swc[0], soil,
                ps->zsoil, dt, 0, qtot);
            rhsts[0] -= tsnsr / denom;
        }    /* This ends section for top soil layer. */
    }

    /* Initialize ddz2 */
    ddz2 = 0.0;
    df1k = df1;

    /* Loop thru the remaining soil layers, repeating the above process (except
     * subsfc or "ground" heat flux not repeated in lower layers)
     * Calculate heat capacity for this soil layer. */
    for (k = 1; k < ps->nlayers; k++)
    {
        hcpct = ws->swc[k] * CH2O + (1.0 - soil->smcmax) * csoil_loc +
            (soil->smcmax - ws->smc[k]) * CP +
            (ws->smc[k] - ws->swc[k]) * CICE;

        /* This section for Layer 2 or greater, but not last layer.
         * Calculate thermal diffusivity for this layer. */
        if (k != ps->nlayers - 1)
        {
            /* Calc the vertical soil temp gradient thru this layer */
            df1n = TDfCnd(ws->smc[k], soil->quartz, soil->smcmax, soil->smcmin,
                ws->swc[k]);

            /* Urban */
            df1n = (lc->isurban) ? 3.24 : df1n;

            denom = 0.5 * (ps->zsoil[k - 1] - ps->zsoil[k + 1]);

            /* Calc the matrix coef, ci, after calc'ng its partial product */
            dtsdz2 = (es->stc[k] - es->stc[k + 1]) / denom;
            ddz2 = 2.0 / (ps->zsoil[k - 1] - ps->zsoil[k + 1]);

            /* If temperature averaging invoked (itavg = true; else skip)
             * Calculate temp at bottom of layer. */
            ci[k] = -df1n * ddz2 / ((ps->zsoil[k - 1] - ps->zsoil[k]) * hcpct);
            if (itavg)
            {
                tbk1 = TBnd(es->stc[k], es->stc[k + 1], ps->zsoil, ps->zbot, k,
                    ps->nlayers);
            }
        }
        else
        {
            /* Special case of bottom soil layer
             * Calculate thermal diffusivity for bottom layer. */
            df1n = TDfCnd(ws->smc[k], soil->quartz, soil->smcmax, soil->smcmin,
                ws->swc[k]);

            /* Urban */
            df1n = (lc->isurban) ? 3.24 : df1n;

            /* Calc the vertical soil temp gradient thru bottom layer. */
            denom = 0.5 * (ps->zsoil[k - 1] + ps->zsoil[k]) - ps->zbot;
            dtsdz2 = (es->stc[k] - ps->tbot) / denom;

            /* Set matrix coef, ci to zero if bottom layer. */
            ci[k] = 0.0;

            /* If temperature averaging invoked (itavg = true; else skip)
             * Calculate temp at bottom of last layer. */
            if (itavg)
            {
                tbk1 = TBnd(es->stc[k], ps->tbot, ps->zsoil, ps->zbot, k,
                    ps->nlayers);
            }    /* This ends special loop for bottom layer. */
        }

        /* Calculate rhsts for this layer after calc'ng a partial product. */
        denom = (ps->zsoil[k] - ps->zsoil[k - 1]) * hcpct;
        rhsts[k] = (df1n * dtsdz2 - df1k * dtsdz) / denom;
        qtot = -1.0 * denom * rhsts[k];

        sice = ws->smc[k] - ws->swc[k];

        if (itavg)
        {
            tavg = TmpAvg(tbk, es->stc[k], tbk1, ps->zsoil, k);
            if (sice > 0.0 || es->stc[k] < TFREEZ || tbk < TFREEZ ||
                tbk1 < TFREEZ)
            {
                SnkSrc(&tsnsr, tavg, ws->smc[k], &ws->swc[k], soil, ps->zsoil,
                    dt, k, qtot);
                rhsts[k] = rhsts[k] - tsnsr / denom;
            }
        }
        else
        {
            if (sice > 0.0 || es->stc[k] < TFREEZ)
            {
                SnkSrc(&tsnsr, es->stc[k], ws->smc[k], &ws->swc[k], soil,
                    ps->zsoil, dt, k, qtot);
                rhsts[k] = rhsts[k] - tsnsr / denom;
            }
        }

        /* Calc matrix coefs, ai, and bi for this layer. */
        ai[k] = -df1k * ddz / ((ps->zsoil[k - 1] - ps->zsoil[k]) * hcpct);
        bi[k] = -(ai[k] + ci[k]);

        /* Reset values of df1, dtsdz, ddz, and tbk for loop to next soil layer.
         */
        tbk = tbk1;
        df1k = df1n;
        dtsdz = dtsdz2;
        ddz = ddz2;
    }
}

void HStep(estate_struct *es, double *rhsts, double dt, int nlayers, double *ai,
    double *bi, double *ci)
{
    /*
     * Calculate/update the soil temperature field.
     */
    int             k;
    double          rhstsin[MAXLYR];
    double          ciin[MAXLYR];

    /* Create finite difference values for use in Rosr12 routine */
    for (k = 0; k < nlayers; k++)
    {
        rhsts[k] *= dt;
        ai[k] *= dt;
        bi[k] = 1.0 + bi[k] * dt;
        ci[k] *= dt;
    }

    /* Copy values for input variables before call to Rosr12 */
    for (k = 0; k < nlayers; k++)
    {
        rhstsin[k] = rhsts[k];
        ciin[k] = ci[k];
    }

    /* Solve the tri-diagonal matrix equation */
    Rosr12(ci, ai, bi, ciin, rhstsin, rhsts, nlayers);

    /* Calc/update the soil temps using matrix solution */
    for (k = 0; k < nlayers; k++)
    {
        es->stc[k] += ci[k];
    }
}

#if defined(_CYCLES_)
void NoPac(const soil_struct *soil, const lc_struct *lc,
    const weather_struct *weather, const cstate_struct *cs, double dt,
    double t24, crop_struct crop[], phystate_struct *ps, wstate_struct *ws,
    wflux_struct *wf, estate_struct *es, eflux_struct *ef)
#else
void NoPac(wstate_struct *ws, wflux_struct *wf, estate_struct *es,
    eflux_struct *ef, phystate_struct *ps, const lc_struct *lc,
    const soil_struct *soil, double dt, double t24)
#endif
{
    /*
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

    prcpf = wf->prcp;

    /* Initialize dew */
    wf->dew = 0.0;

    /* Initialize evap terms */
    wf->edir = 0.0;
    wf->ec = 0.0;
    for (k = 0; k < ps->nlayers; k++)
    {
        wf->et[k] = 0.0;
    }
    wf->ett = 0.0;

    if (wf->etp > 0.0)
    {
#if defined(_CYCLES_)
        Evapo(soil, lc, weather, ps, es, cs, crop, ws, wf);
#else
        Evapo(ws, wf, ps, lc, soil, dt);
#endif

        wf->eta = wf->etns;
    }
    else
    {
        /* If etp < 0, assume dew forms (transform etp into dew and reinitialize
         * etp to zero). */
        wf->dew = -wf->etp;

        /* Add dew amount to prcp */
        prcpf += wf->dew;
    }

#if defined(_CYCLES_)
    prcpf += wf->irrig / RHOH2O / DAYINSEC;
#endif

#if defined(_CYCLES_)
    ResidueWetting(prcpf * RHOH2O * DAYINSEC, cs, ps, ws, wf);
#else
    PcpDrp(ws, wf, lc, prcpf, dt);
#endif

    /* Based on etp and e values, determine beta */
    if (wf->etp <= 0.0)
    {
        ps->beta = (wf->etp < 0.0) ? 1.0 : 0.0;
        wf->eta = wf->etp;
    }
    else
    {
        ps->beta = wf->eta / wf->etp;
    }

    /* Get soil thermal diffusivity/conductivity for top soil lyr, Calc.
     * adjusted top lyr soil temp and adjusted soil flux, then call ShFlx to
     * compute/update soil heat flux and soil temps. */
    df1 = TDfCnd(ws->smc[0], soil->quartz, soil->smcmax, soil->smcmin,
        ws->swc[0]);

    /* Urban */
    df1 = (lc->isurban) ? 3.24 : df1;

    /* Vegetation greenness fraction reduction in subsurface heat flux via
     * reduction factor, which is convenient to apply here to thermal
     * diffusivity that is later used in HRT to compute sub sfc heat flux (see
     * additional comments on veg effect sub-sfc heat flx in function SFlx) */
    df1 *= exp(ps->sbeta * lc->shdfac);

    /* Compute intermediate terms passed to routine HRT (via routine ShFlx
     * below) for use in computing subsurface heat flux in HRT */
    yynum = ef->fdown - ps->emissi * SIGMA * t24;
    yy = es->sfctmp +
        (yynum / ps->rch + es->th2 - es->sfctmp - ps->beta * ps->epsca) /
        ps->rr;

    zz1 = df1 / (-0.5 * ps->zsoil[0] * ps->rch * ps->rr) + 1.0;

    ShFlx(ws, es, ps, lc, soil, dt, yy, zz1, df1);

    /* In the no snowpack case update the grnd (skin) temperature in response to
     * the updated soil temperature profile */
    es->t1 = (yy + (zz1 - 1.0) * es->stc[0]) / zz1;

    /* Calculate surface soil heat flux */
    ef->ssoil = df1 * (es->stc[0] - es->t1) / (0.5 * ps->zsoil[0]);

    /* Set flx1 and flx3 (snopack phase change heat fluxes) to zero since they
     * are not used here in SnoPac. flx2 (freezing rain heat flux) was similarly
     * initialized in the Penman routine. */
    ef->flx1 = CPH2O * wf->prcp * 1000.0 * (es->t1 - es->sfctmp);
    ef->flx3 = 0.0;
}

void PcpDrp(wstate_struct *ws, wflux_struct *wf, const lc_struct *lc,
    double prcp, double dt)
{
    /*
     * Separated from Noah SmFlx function
     * The canopy moisture content (cmc) is updated.
     */
    double          excess;
    double          rhsct;
    double          trhsct;
    const double    KD = 6.54E-7;
    const double    BFACTR = 3.89;

    /* Compute the right hand side of the canopy eqn term (rhsct)
     * Convert rhsct (a rate) to trhsct (an amount) and add it to existing cmc.
     * If resulting amt exceeds max capacity, it becomes drip and will fall to
     * the grnd. */
    rhsct = lc->shdfac * prcp - wf->ec;
    wf->drip = 0.0;
    trhsct = dt * rhsct;
    excess = ws->cmc + trhsct;

    /* Pcpdrp is the combined prcp and drip (from cmc) that goes into the soil
     */
    if (excess > 0.0)
    {
        /* PIHM drip calculation following Rutter and Mortan (1977 JAE) */
        wf->drip = (excess >= ws->cmcmax) ?
            (KD * ws->cmcmax * exp(BFACTR)) + (excess - ws->cmcmax) / dt :
            (KD * ws->cmcmax * exp(BFACTR * excess / ws->cmcmax));

        wf->drip = (wf->drip > excess / dt) ? excess / dt : wf->drip;
    }

    rhsct -= wf->drip;

    wf->pcpdrp = (1.0 - lc->shdfac) * prcp + wf->drip;

    /* Update canopy water content/interception (cmc). Convert rhsct to an
     * "amount" value and add to previous cmc value to get new cmc. */
    ws->cmc += dt * rhsct;

    ws->cmc = (ws->cmc < 1.0E-20) ? 0.0 : ws->cmc;
    ws->cmc = (ws->cmc < ws->cmcmax) ? ws->cmc : ws->cmcmax;
}

void Penman(wflux_struct *wf, const estate_struct *es, eflux_struct *ef,
    phystate_struct *ps, double *t24, double t2v, int snowng, int frzgra)
{
    /*
     * Calculate potential evaporation for the current point. Various partial
     * sums/products are also calculated and passed back to the calling routine
     * for later use.
     */
    double          a;
    double          delta;
    double          fnet;
    double          rad;
    double          rho;
    double          emissi;
    double          elcp1;
    double          lvs;
    const double    ELCP = 2.4888E3;

    /* Prepare partial quantities for Penman equation. */
    emissi = ps->emissi;
    elcp1 = (1.0 - ps->sncovr) * ELCP + ps->sncovr * ELCP * LSUBS / LVH2O;
    lvs = (1.0 - ps->sncovr) * LVH2O + ps->sncovr * LSUBS;

    ef->flx2 = 0.0;
    delta = elcp1 * ps->dqsdt2;
    *t24 = es->sfctmp * es->sfctmp * es->sfctmp * es->sfctmp;
    ps->rr = emissi * *t24 * 6.48E-8 / (ps->sfcprs * ps->ch) + 1.0;
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
    fnet = ef->fdown - emissi * SIGMA * *t24 - ef->ssoil;
    if (frzgra)
    {
        ef->flx2 = -LSUBF * wf->prcp * 1000.0;
        fnet -= ef->flx2;
    }

    /* Finish Penman equation calculations */
    rad = fnet / ps->rch + es->th2 - es->sfctmp;
    a = elcp1 * (ps->q2sat - ps->q2);
    ps->epsca = (a * ps->rr + rad * delta) / (delta + ps->rr);
    wf->etp = ps->epsca * ps->rch / lvs / 1000.0;
}

void Rosr12(double *p, const double *a, const double *b, double *c,
    const double *d, double *delta, int nlayers)
{
    /*
     * Invert (solve) the tri-diagonal matrix problem shown below:
     * ###                                            ### ###  ###   ###  ###
     * #b[0], c[0],  0  ,   0   ,   0   , . . . ,   0   # #      #   #      #
     * #a[1], b[1], c[1],   0   ,   0   , . . . ,   0   # #      #   #      #
     * # 0  , a[2], b[2],  c[2] ,   0   , . . . ,   0   # #      #   # d[2] #
     * # 0  ,  0  , a[3],  b[3] ,  c[3] , . . . ,   0   # # p[3] #   # d[3] #
     * # 0  ,  0  ,  0  ,  a[4] ,  b[4] , . . . ,   0   # # p[4] #   # d[4] #
     * # .                                          .   # #  .   # = #   .  #
     * # .                                          .   # #  .   #   #   .  #
     * # .                                          .   # #  .   #   #   .  #
     * # 0  , . . . , 0 , a[m-2], b[m-2], c[m-2],   0   # #p(m-2)#   #d[m-2]#
     * # 0  , . . . , 0 ,   0   , a[m-1], b[m-1], c[m-1]# #p(m-1)#   #d[m-1]#
     * # 0  , . . . , 0 ,   0   ,   0   ,  a(m) ,  b[m] # # p(m) #   # d[m] #
     * ###                                            ### ###  ###   ###  ###
     */
    int             k, kk;

    /* Initialize eqn coef c for the lowest soil layer */
    c[nlayers - 1] = 0.0;

    /* Solve the coefs for the 1st soil layer */
    p[0] = -c[0] / b[0];
    delta[0] = d[0] / b[0];

    /* Solve the coefs for soil layers 2 thru nlayers */
    for (k = 1; k < nlayers; k++)
    {
        p[k] = -c[k] * (1.0 / (b[k] + a[k] * p[k - 1]));
        delta[k] =
            (d[k] - a[k] * delta[k - 1]) * (1.0 / (b[k] + a[k] * p[k - 1]));
    }

    /* Set p to delta for lowest soil layer */
    p[nlayers - 1] = delta[nlayers - 1];

    /* Adjust p for soil layers 2 thru nlayers */
    for (k = 1; k < nlayers; k++)
    {
        kk = nlayers - k - 1;
        p[kk] = p[kk] * p[kk + 1] + delta[kk];
    }
}

void ShFlx(wstate_struct *ws, estate_struct *es, const phystate_struct *ps,
    const lc_struct *lc, const soil_struct *soil, double dt, double yy,
    double zz1, double df1)
{
    /*
     * Update the temperature state of the soil column based on the thermal
     * diffusion equation and update the frozen soil moisture content based on
     * the temperature.
     */
    double          ai[MAXLYR], bi[MAXLYR], ci[MAXLYR];
    double          rhsts[MAXLYR];

    /* HRT routine calcs the right hand side of the soil temp dif eqn */
    /* Land case */
    HRT(ws, es, ps, lc, soil, rhsts, yy, zz1, dt, df1, ai, bi, ci);

    HStep(es, rhsts, dt, ps->nlayers, ai, bi, ci);
}


void SmFlx(wstate_struct *ws, wflux_struct *wf, phystate_struct *ps,
    const soil_struct *soil, double dt)
{
    /*
     * Calculate soil moisture flux. The soil moisture content (smc - a per unit
     * volume measurement) is a dependent variable that is updated with
     * prognostic eqns.
     * Frozen ground version: new states added: swc, and frozen ground
     * correction factor, frzfact and parameter slope.
     */
    int             i;
    double          ai[MAXLYR], bi[MAXLYR], ci[MAXLYR];
    double          rhstt[MAXLYR];
    double          sice[MAXLYR];

    /* Store ice content at each soil layer before calling SRT and SStep */
    for (i = 0; i < ps->nlayers; i++)
    {
        sice[i] = ws->smc[i] - ws->swc[i];
    }

    /* Call subroutines SRT and SStep to solve the soil moisture tendency
     * equations.
     * Call the SRT/SStep function in the manner of time scheme "d" (implicit
     * state, explicit coefficient) of Section 2 of Kalnay and Kanamitsu pcpdrp
     * is units of m/s, zsoil is negative depth in m
     * According to Dr. Ken Mitchell's suggestion, add the second contraint to
     * remove numerical instability of runoff and soil moisture
     *
     * Frozen ground version:
     * smc states replaced by swc states in SRT subr. swc & sice states
     * included in SStep subr. Frozen ground correction factor, frzfact added.
     * All water balance calculations using unfrozen water */
    if (0 == ps->nwtbl)
    {
        /* Special case: all soil layers are saturated */
        for (i = 0; i < ps->nlayers; i++)
        {
            ws->smc[i] = soil->smcmax;
            ws->swc[i] = ws->smc[i] - sice[i];
        }
    }
    else
    {
        SRT(ws, wf, ps, soil, rhstt, sice, ai, bi, ci);
        SStep(ws, wf, ps, soil, rhstt, sice, ai, bi, ci, dt);
    }
}

double SnFrac(double sneqv, double snup, double salp)
{
    /*
     * Calculate snow fraction (0 -> 1)
     */
    double          rsnow;
    double          sncovr;

    /* Snup is veg-class dependent snowdepth threshold above which snocvr = 1 */
    if (sneqv < snup)
    {
        rsnow = sneqv / snup;
        sncovr = 1.0 - (exp(-salp * rsnow) - rsnow * exp(-salp));
    }
    else
    {
        sncovr = 1.0;
    }

#if NOT_YET_IMPLEMENTED
    /* Formulation of Dickinson et al. 1986 */
    z0n = 0.035;
    sncovr = snowh / (snowh + 5.0 * z0n);

    /* Formulation of Marshall et al. 1994 */
    sncovr = sneqv / (sneqv + 2.0 * z0n);
#endif

    return sncovr;
}

void SnkSrc(double *tsnsr, double tavg, double smc, double *swc,
    const soil_struct *soil, const double *zsoil, double dt, int k, double qtot)
{
    /*
     * Calculate sink/source term of the thermal diffusion equation. (swc) is
     * available liquid water.
     */
    double          dz;
    double          freew;
    double          xh2o;

    double          DH2O = 1.0000E3;

    dz = (0 == k) ? -zsoil[0] : zsoil[k - 1] - zsoil[k];

    /* Via function FrH2O, compute potential or 'equilibrium' unfrozen
     * supercooled free water for given soil type and soil layer temperature.
     * Function FrH2O invokes Eqn (17) from V. Koren et al (1999, JGR, Vol. 104,
     * Pg 19573). (Aside: latter eqn in journal in centigrade units.
     * routine FrH2O use form of eqn in kelvin units.) */
    freew = FrH2O(tavg, smc, *swc, soil);

    /* In next block of code, invoke Eqn 18 of V. Koren et al (1999, JGR,
     * Vol. 104, Pg 19573.) that is, first estimate the new amount of liquid
     * water, 'xh2o', implied by the sum of (1) the liquid water at the begin of
     * current time step, and (2) the freeze of thaw change in liquid water
     * implied by the heat flux 'qtot' passed in from routine HRT. Second,
     * determine if xh2o needs to be bounded by 'freew' (equil amt) or if
     * 'freew' needs to be bounded by xh2o. */
    xh2o = *swc + qtot * dt / (DH2O * LSUBF * dz);

    /* First, if freezing and remaining liquid less than lower bound, then
     * reduce extent of freezing, thereby letting some or all of heat flux qtot
     * cool the soil temp later in routine HRT. */
    if (xh2o < *swc && xh2o < freew)
    {
        xh2o = (freew > *swc) ? *swc : freew;
    }

    /* Second, if thawing and the increase in liquid water greater than upper
     * bound, then reduce extent of thaw, thereby letting some or all of heat
     * flux qtot warm the soil temp later in routine HRT. */
    if (xh2o > *swc && xh2o > freew)
    {
        xh2o = (freew < *swc) ? *swc : freew;
    }

    xh2o = (xh2o < 0.0) ? 0.0 : xh2o;
    xh2o = (xh2o > smc) ? smc : xh2o;

    /* Calculate phase-change heat source/sink term for use in routine HRT and
     * update liquid water to reflect final freeze/thaw increment. */
    *tsnsr = -DH2O * LSUBF * dz * (xh2o - *swc) / dt;

    *swc = xh2o;
}

#if defined(_CYCLES_)
void SnoPac(const soil_struct *soil, const lc_struct *lc,
    const weather_struct *weather, const cstate_struct *cs, int snowng,
    double dt, double t24, double prcpf, double df1, crop_struct crop[],
    phystate_struct *ps, wstate_struct *ws, wflux_struct *wf, estate_struct *es,
    eflux_struct *ef)
#else
void SnoPac(wstate_struct *ws, wflux_struct *wf, estate_struct *es,
    eflux_struct *ef, phystate_struct *ps, const lc_struct *lc,
    const soil_struct *soil, int snowng, double dt, double t24, double prcpf,
    double df1)
#endif
{
    /*
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
    double          yy;
    double          zz1;
    const double    ESDMIN = 1.0E-6;
    const double    SNOEXP = 2.0;

    /* Initialize evap terms. */
    wf->dew = 0.0;
    wf->edir = 0.0;
    wf->ec = 0.0;

    for (k = 0; k < ps->nlayers; k++)
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
        if (ps->ribb >= 0.1 && ef->fdown > 150.0)
        {
            wf->etp =
                (((wf->etp * (1.0 - ps->ribb) < 0.0) ?
                wf->etp * (1.0 - ps->ribb) : 0.0)
                * ps->sncovr / 0.980 + wf->etp * (0.980 - ps->sncovr)) / 0.980;
        }

        ps->beta = (wf->etp == 0.0) ? 0.0 : ps->beta;

        wf->dew = -wf->etp;
        esnow2 = wf->etp * dt;
        etanrg = wf->etp * 1000.0 *
            ((1.0 - ps->sncovr) * LVH2O + ps->sncovr * LSUBS);
    }
    else
    {
        /* Land case */
        if (ps->sncovr < 1.0)
        {
#if defined(_CYCLES_)
            Evapo(soil, lc, weather, ps, es, cs, crop, ws, wf);
#else
            Evapo(ws, wf, ps, lc, soil, dt);
#endif

            wf->edir *= (1.0 - ps->sncovr);
            wf->ec *= (1.0 - ps->sncovr);
            for (k = 0; k < ps->nlayers; k++)
            {
                wf->et[k] = wf->et[k] * (1.0 - ps->sncovr);
            }
            wf->ett *= (1.0 - ps->sncovr);
            wf->etns *= (1.0 - ps->sncovr);
        }
        wf->esnow = wf->etp * ps->sncovr;
        esnow2 = wf->esnow * dt;
        etanrg = wf->esnow * 1000.0 * LSUBS + wf->etns * 1000.0 * LVH2O;
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

    dsoil = -0.5 * ps->zsoil[0];
    dtot = ps->snowh + dsoil;
    denom = 1.0 + df1 / (dtot * ps->rr * ps->rch);

    /* Calculate an 'effective snow-grnd sfc temp' (t12) based on heat fluxes
     * between the snow pack and the soil and on net radiation.
     * Include flx1 (precip-snow sfc) and flx2 (freezing rain latent heat)
     * fluxes. flx1 from above, flx2 brought in via common block rite.
     * flx2 reflects freezing rain latent heat flux using t1 calculated in
     * Penman. */
    /* Surface emissivity weighted by snow cover fraction */
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
    if (t12 <= TFREEZ)
    {
        /* Sub-freezing block */
        es->t1 = t12;
        ef->ssoil = df1 * (es->t1 - es->stc[0]) / dtot;

        if (ws->sneqv - esnow2 > 0.0)
        {
            ws->sneqv -= esnow2;
        }
        else
        {
            ws->sneqv = 0.0;
            esnow2 = ws->sneqv;
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
        es->t1 = TFREEZ * pow(ps->sncovr, SNOEXP) +
            t12 * (1.0 - pow(ps->sncovr, SNOEXP));
        ps->beta = 1.0;

        ef->ssoil = df1 * (es->t1 - es->stc[0]) / dtot;

        if (ws->sneqv - esnow2 <= ESDMIN)
        {
            /* If potential evap (sublimation) greater than depth of snowpack.
             * beta < 1
             * snowpack has sublimated away, set depth to zero. */
            ws->sneqv = 0.0;
            ex = 0.0;
            wf->snomlt = 0.0;
            ef->flx3 = 0.0;
        }
        else
        {
            /* Sublimation less than depth of snowpack
             * Snowpack (esd) reduced by esnow2 (depth of sublimated snow) */
            ws->sneqv -= esnow2;
            etp3 = wf->etp * 1000.0 * LVH2O;
            seh = ps->rch * (es->t1 - es->th2);
            t14 = es->t1 * es->t1;
            t14 = t14 * t14;
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

        /* If non-glacial land, add snowmelt rate (ex) to precip rate to be used
         * in subroutine SmFlx (soil moisture evolution) via infiltration.
         * Runoff/baseflow later near the end of sflx (after return from call to
         * subroutine SnoPac) */
        prcpf += ex;

        /* Set the effective potnl evapotransp (etp1) to zero since this is snow
         * case, so surface evap not calculated from edir, ec, or ett in SmFlx
         * (below).
         * SmFlx returns updated soil moisture values for non-glacial land. */
    }

#if defined(_CYCLES_)
    prcpf += wf->irrig / RHOH2O / DAYINSEC;
#endif

#if defined(_CYCLES_)
    ResidueWetting(prcpf * RHOH2O * DAYINSEC, cs, ps, ws, wf);
#else
    PcpDrp(ws, wf, lc, prcpf, dt);
#endif

    /* Before call ShFlx in this snowpack case, set zz1 and yy arguments to
     * special values that ensure that ground heat flux calculated in ShFlx
     * matches that already computer for below the snowpack, thus the sfc heat
     * flux to be computed in ShFlx will effectively be the flux at the snow top
     * surface.  t11 is a dummy argument so we will not use the skin temp value
     * as revised by ShFlx. */
    zz1 = 1.0;
    yy = es->stc[0] - 0.5 * ef->ssoil * ps->zsoil[0] * zz1 / df1;

    /* ShFlx will calc/update the soil temps */
    ShFlx(ws, es, ps, lc, soil, dt, yy, zz1, df1);

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

void SnowPack(double esd, double dtsec, double *snowh, double *sndens,
    double tsnow, double tsoil)
{
    /*
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
    tsnowc = tsnow - TFREEZ;
    tsoilc = tsoil - TFREEZ;

    /* Calculating of average temperature of snow pack */
    tavgc = 0.5 * (tsnowc + tsoilc);

    esdcx = (esdc > 1.0E-2) ? esdc : 1.0E-2;

    /* Calculating of snow depth and density as a result of compaction
     *   sndens = ds0 * (exp (bfac * esd) - 1.) / (bfac * esd)
     *   bfac = dthr * c1 * exp (0.08 * tavgc - c2 * ds0)
     * Note: bfac * esd in sndens eqn above has to be carefully treated
     * numerically below:
     *   c1 is the fractional increase in density (1/(cm*hr))
     *   c2 is a constant (cm3/g) kojima estimated as 21 cms/g */
    bfac = dthr * C1 * exp(0.08 * tavgc - C2 * *sndens);

    /* The function of the form (e**x - 1) / x embedded in above expression for
     * dsx was causing numerical difficulties when the denominator "x"
     * (i.e. bfac*esdc) became zero or approached zero (despite the fact that
     * the analytical function (e**x-1)/x has a well defined limit as "x"
     * approaches zero), hence below we replace the (e**x-1)/x expression with
     * an equivalent, numerically well-behaved polynomial expansion.
     * Number of terms of polynomial expansion, and hence its accuracy, is
     * governed by iteration limit "ipol".
     *   ipol greater than 9 only makes a difference on double precision
     *   (relative errors given in percent %).
     *     ipol=9, for rel.error <~ 1.6 e-6 % (8 significant digits)
     *     ipol=8, for rel.error <~ 1.8 e-5 % (7 significant digits)
     *     ipol=7, for rel.error <~ 1.8 e-4 % ... */
    ipol = 4;
    pexp = 0.0;
    for (j = ipol; j > 0; j--)
    {
        pexp = (1.0 + pexp) * bfac * esdcx / (double)(j + 1);
    }

    pexp += 1.0;

    /* Set upper/lower limit on snow density */
    dsx = *sndens * pexp;
    dsx = (dsx > 0.40) ? 0.40 : dsx;
    dsx = (dsx < 0.05) ? 0.05 : dsx;

    /* Update of snow depth and density depending on liquid water during
     * snowmelt. Assumed that 13% of liquid water can be stored in snow per day
     * during snowmelt till snow density 0.40. */
    *sndens = dsx;
    if (tsnowc >= 0.0)
    {
        dw = 0.13 * dthr / 24.0;
        *sndens = *sndens * (1.0 - dw) + dw;
        *sndens = (*sndens >= 0.40) ? 0.40 : *sndens;
    }
    /* Calculate snow depth (cm) from snow water equivalent and snow density.
     * Change snow depth units to meters */
    snowhc = esdc / *sndens;
    *snowh = snowhc * 0.01;
}

double Snowz0(double sncovr, double z0brd, double snowh)
{
    /*
     * Calculate total roughness length over snow
     */
    const double    Z0S = 0.001;    /* Snow roughness length (m) */
    double          burial;
    double          z0eff;
    double          z0;

    burial = 7.0 * z0brd - snowh;
    z0eff = (burial <= 0.0007) ? Z0S : (burial / 7.0);

    z0 = (1.0 - sncovr) * z0brd + sncovr * z0eff;

    return z0;
}

void SnowNew(const estate_struct *es, double newsn, phystate_struct *ps)
{
    /*
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

    /* Calculating new snowfall density depending on temperature equation from
     * Gottlib L. 'A general runoff model for snowcovered and glacierized
     * basin', 6th Nordic Hydrological Conference, Vemadolen, Sweden, 1980,
     * 172-177pp. */
    tempc = es->sfctmp - 273.15;

    dsnew = (tempc <= -15.0) ? 0.05 : 0.05 + 0.0017 * pow(tempc + 15.0, 1.5);

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
    snowhc += hnewc;
    ps->snowh = snowhc * 0.01;
}

void SRT(wstate_struct *ws, wflux_struct *wf, phystate_struct *ps,
    const soil_struct *soil, double *rhstt, double *sice, double *ai,
    double *bi, double *ci)
{
    /*
     * Calculate the right hand side of the time tendency term of the soil water
     * diffusion equation. Also to compute (prepare) the matrix coefficients for
     * the tri-diagonal matrix of the implicit time scheme.
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
    double          dsmdz, dsmdz2;
    double          dice;
    const double    CVFRZ = 3.0;
    double          acrt;
    double          sum;
    int             ialp1;
    /* Frozen ground version:
     * Reference frozen ground parameter, cvfrz, is a shape parameter of areal
     * distribution function of soil ice content which equals 1/cv.
     * Cv is a coefficient of spatial variation of soil ice content. Based on
     * field data cv depends on areal mean of frozen depth, and it close to
     * constant = 0.6 if areal mean frozen depth is above 20 cm. That is why
     * parameter cvfrz = 3 (int{1/0.6*0.6}). Current logic doesn't allow cvfrz
     * be bigger than 3 */
    /* Let sicemax be the greatest, if any, frozen water content within soil
     * layers. */
    iohinf = 1;
    sicemax = 0.0;

    for (ks = 0; ks < ps->nlayers; ks++)
    {
        sicemax = (sice[ks] > sicemax) ? sice[ks] : sicemax;
    }

    /* Calculate infiltration reduction factor due to frozen soil */
    dice = -ps->zsoil[0] * sice[0];
    for (ks = 1; ks < ps->nlayers; ks++)
    {
        dice += (ps->zsoil[ks - 1] - ps->zsoil[ks]) * sice[ks];
    }

    ps->fcr = 1.0;

    if (dice > 1.0E-2)
    {
        acrt = CVFRZ * ps->frzx / dice;
        sum = 1.0;
        ialp1 = roundi(CVFRZ) - 1;
        for (j = 1; j < ialp1 + 1; j++)
        {
            k = 1;
            for (jj = j + 1; jj < ialp1 + 1; jj++)
            {
                k *= jj;
            }
            sum += pow(acrt, CVFRZ - (double)j) / (double)k;
        }

        ps->fcr = 1.0 - exp(-acrt) * sum;
    }

    /* Determine rainfall infiltration rate and runoff */
    pddum = wf->eqv_infil;

    mxsmc = ws->swc[0];

    dsmdz = (ws->swc[0] - ws->swc[1]) / (-0.5 * ps->zsoil[1]);
    WDfCnd(&wdf, &wcnd, mxsmc, sicemax, soil);

    /* Calc the matrix coefficients ai, bi, and ci for the top layer */
    ddz = 1.0 / (-0.5 * ps->zsoil[1]);
    ai[0] = 0.0;
    bi[0] = wdf * ddz / (-ps->zsoil[0]);
    ci[0] = -bi[0];
    /* Calc rhstt for the top layer after calc'ng the vertical soil moisture
     * gradient btwn the top and next to top layers. */
    rhstt[0] = (wdf * dsmdz + wcnd - pddum + wf->edir + wf->et[0]) /
        ps->zsoil[0];

    rhstt[0] += wf->runoff2_lyr[0] / ps->zsoil[0];

    /* Loop thru the remaining soil layers, repeating the abv process */
    /* Initialize ddz2 */
    ddz2 = 0.0;
    for (k = 1; k < ps->nlayers; k++)
    {
        denom2 = (ps->zsoil[k - 1] - ps->zsoil[k]);
        if (k < ps->nlayers - 1)
        {
            mxsmc2 = ws->swc[k];
            denom = ps->zsoil[k - 1] - ps->zsoil[k + 1];
            dsmdz2 = (ws->swc[k] - ws->swc[k + 1]) / (denom * 0.5);
            WDfCnd(&wdf2, &wcnd2, mxsmc2, sicemax, soil);
            /* Calc some partial products for later use in calc'ng rhstt
             * Calc the matrix coef, ci, after calc'ng its partial product */
            ddz2 = 2.0 / denom;
            ci[k] = -wdf2 * ddz2 / denom2;
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

        numer += wf->runoff2_lyr[k];

        rhstt[k] = numer / (-denom2);

        /* Calc matrix coefs, ai, and bi for this layer */
        ai[k] = -wdf * ddz / denom2;
        bi[k] = -(ai[k] + ci[k]);

        /* Reset values of wdf, wcnd, dsmdz, and ddz for loop to next lyr */
        if (k != ps->nlayers - 1)
        {
            wdf = wdf2;
            wcnd = wcnd2;
            dsmdz = dsmdz2;
            ddz = ddz2;
        }
    }
}

void SStep(wstate_struct *ws, wflux_struct *wf, phystate_struct *ps,
    const soil_struct *soil, double *rhstt, double *sice, double *ai,
    double *bi, double *ci, double dt)
{
    /*
     * Calculate/update soil moisture content values and canopy moisture content
     * values.
     */
    int             k;
    double          rhsttin[MAXLYR];
    double          ciin[MAXLYR];
    double          sh2o0[MAXLYR];

    /* Create 'amount' values of variables to be input to the tri-diagonal
     * matrix routine. */
    for (k = 0; k < ps->nlayers; k++)
    {
        rhstt[k] *= dt;
        ai[k] *= dt;
        bi[k] = 1.0 + bi[k] * dt;
        ci[k] *= dt;
        sh2o0[k] = ws->swc[k];
        wf->smflx[k] = 0.0;
    }

    /* Copy values for input variables before call to Rosr12 */
    for (k = 0; k < ps->nlayers; k++)
    {
        rhsttin[k] = rhstt[k];
        ciin[k] = ci[k];
    }

    /* Call Rosr12 to solve the tri-diagonal matrix */
    Rosr12(ci, ai, bi, ciin, rhsttin, rhstt, ps->nlayers);

    /* Sum the previous smc value and the matrix solution to get a new value. */
    for (k = 0; k < ps->nlayers; k++)
    {
        ws->swc[k] += ci[k];
    }

    AdjSmProf(soil, ps, sice, dt, wf, ws);

    /* Calculate soil moisture flux within soil layers */
    for (k = ps->nlayers - 1; k > 0; k--)
    {
        /* Positive smflx[k] is flux out of soil layer k */
        wf->smflx[k - 1] =
            (ws->swc[k] - sh2o0[k]) * ps->soil_depth[k] / dt +
            wf->runoff2_lyr[k] + wf->et[k] + wf->smflx[k];
    }
}

double TBnd(double tu, double tb, const double *zsoil, double zbot, int k,
    int nlayers)
{
    /*
     * Calculate temperature on the boundary of the layer by interpolation of
     * the middle layer temperatures */
    double          zb, zup;
    double          tbnd1;

    /* Use surface temperature on the top of the first layer */
    zup = (k == 0) ? 0.0 : zsoil[k - 1];

    /* Use depth of the constant bottom temperature when interpolate temperature
     * into the last layer boundary */
    zb = (k == nlayers - 1) ? 2.0 * zbot - zsoil[k] : zsoil[k + 1];

    /* Linear interpolation between the average layer temperatures */
    tbnd1 = tu + (tb - tu) * (zup - zsoil[k]) / (zup - zb);

    return tbnd1;
}

double TDfCnd(double smc, double qz, double smcmax, double smcmin, double swc)
{
    /*
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
     * Pablo Grunmann, 08/17/98
     * Refs.:
     *  Farouki, O. T.,1986: Thermal properties of soils. Series on rock and
     *    soil mechanics, Vol. 11, trans tech, 136 pp.
     *  Johansen, O., 1975: Thermal conductivity of soils. Ph.D. thesis,
     *    University of Trondheim
     *  Peters-Lidard, C. D., et al., 1998: The effect of soil thermal
     *    conductivity parameterization on surface energy fluxes and
     *    temperatures. Journal of the Atmospheric Sciences, Vol. 55,
     *    pp. 1209-1224.
     */
    satratio = (smc - smcmin) / (smcmax - smcmin);

    /* Thermal conductivity of "other" soil components */
    thko = 2.0;

    /* Solids' conductivity */
    thks = pow(THKQTZ, qz) * pow(thko, 1.0 - qz);

    /* Unfrozen fraction (from 1, i.e., 100% liquid, to 0 (100% frozen)) */
    xunfroz = swc / smc;

    /* Unfrozen volume for saturation (porosity * xunfroz) */
    xu = xunfroz * smcmax;

    /* Saturated thermal conductivity */
    thksat = pow(thks, 1.0 - smcmax) * pow(THKICE, smcmax - xu) * pow(THKW, xu);

    /* Dry density in kg/m3 */
    gammd = (1.0 - smcmax) * 2700.0;

    /* Dry thermal conductivity in W m-1 K-1 */
    thkdry = (0.135 * gammd + 64.7) / (2700.0 - 0.947 * gammd);

    if (swc + 0.0005 < smc)
    {
        /* Frozen */
        ake = satratio;
    }
    else
    {
        /* Unfrozen
         * range of validity for the kersten number (ake) Kersten number (using
         * "fine" formula, valid for soils containing at least 5% of particles
         * with diameter less than 2.e-6 meters.)
         * (for "coarse" formula, see Peters-Lidard et al., 1998). */
        ake = (satratio > 0.1) ? log10(satratio) + 1.0 : 0.0;
    }

    /* Thermal conductivity */
    df = ake * (thksat - thkdry) + thkdry;

    return df;
}

double TmpAvg(double tup, double tm, double tdn, const double *zsoil, int k)
{
    /*
     * Calculate soil layer average temperature (tavg) in freezing/thawing layer
     * using up, down, and middle layer temperatures (tup, tdn, tm), where tup
     * is at top boundary of layer, tdn is at bottom boundary of layer.
     * Tm is layer prognostic state temperature.
     */
    double          dz;
    double          dzh;
    double          x0;
    double          xdn;
    double          xup;
    double          tavg;

    dz = (k == 0) ? -zsoil[0] : zsoil[k - 1] - zsoil[k];

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
                tavg = 0.5 *
                    (tup * dzh + tm * (dzh + x0) + TFREEZ * (2.0 * dzh - x0)) /
                    dz;
            }
        }
        else
        {
            if (tdn < TFREEZ)
            {
                /* tup < TFREEZ, tm > TFREEZ, tdn < TFREEZ */
                xup = (TFREEZ - tup) * dzh / (tm - tup);
                xdn = dzh - (TFREEZ - tm) * dzh / (tdn - tm);
                tavg = 0.5 *
                    (tup * xup + TFREEZ * (2.0 * dz - xup - xdn) + tdn * xdn) /
                    dz;
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
                tavg = 0.5 *
                    (TFREEZ * (dz - xup) + tm * (dzh + xup) + tdn * dzh) / dz;
            }
            else
            {
                /* tup > TFREEZ, tm < TFREEZ, tdn > TFREEZ */
                xup = dzh - (TFREEZ - tup) * dzh / (tm - tup);
                xdn = (TFREEZ - tm) * dzh / (tdn - tm);
                tavg = 0.5 *
                    (TFREEZ * (2.0 * dz - xup - xdn) + tm * (xup + xdn)) / dz;
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

    return tavg;
}

void Transp(const wstate_struct *ws, wflux_struct *wf, const phystate_struct *ps,
    const lc_struct *lc, const soil_struct *soil)
{
    /*
     * Calculate transpiration for the veg class
     */
    int             i, k;
    double          denom;
    double          etpa;
    double          gx[MAXLYR];
    double          rtx, sgx;

    /* Initialize plant transp to zero for all soil layers. */
    for (k = 0; k < ps->nlayers; k++)
    {
        wf->et[k] = 0.0;
    }

    /* Calculate an 'adjusted' potential transpiration
     * If statement below to avoid tangent linear problems near zero
     * Note: gx and other terms below redistribute transpiration by layer,
     * et(k), as a function of soil moisture availability, while preserving
     * total etpa. */
    etpa = (ws->cmc != 0.0) ?
        lc->shdfac * ps->pc * wf->etp *
        (1.0 - pow(ws->cmc / ws->cmcmax, lc->cfactr)) :
        lc->shdfac * ps->pc * wf->etp;

    sgx = 0.0;
    for (i = 0; i < ps->nroot; i++)
    {
        gx[i] = (ws->smc[i] - soil->smcwlt) / (soil->smcref - soil->smcwlt);
        gx[i] = (gx[i] < 0.0) ? 0.0 : gx[i];
        gx[i] = (gx[i] > 1.0) ? 1.0 : gx[i];
        sgx += gx[i];
    }
    sgx /= (double)ps->nroot;

    denom = 0.0;
    for (i = 0; i < ps->nroot; i++)
    {
        rtx = ps->rtdis[i] + gx[i] - sgx;
        gx[i] *= ((rtx > 0.0) ? rtx : 0.0);
        denom += gx[i];
    }

    denom = (denom <= 0.0) ? 1.0 : denom;

    for (i = 0; i < ps->nroot; i++)
    {
        wf->et[i] = etpa * gx[i] / denom;
    }
# if NOT_YET_IMPLEMENTED
    /* Above code assumes a vertically uniform root distribution
     * Code below tests a variable root distribution */
    wf->et[0] = (zsoil[0] / zsoil[ps->nroot - 1]) * gx * etpa;
    wf->et[0] = (zsoil[0] / zsoil[ps->nroot - 1]) * etpa;
    /* Using root distribution as weighting factor */
    wf->et[0] = (lc->rtdis[0] * etpa);
    wf->et[0] = etpa * part[0];
    /* Loop down thru the soil layers repeating the operation above, but using
     * the thickness of the soil layer (rather than the absolute depth of each
     * layer) in the final calculation. */
    for (k = 0; k < ps->nroot; k++)
    {
        gx = (ws->smc[k] - soil->smcwlt ) / (soil->smcref - soil->smcwlt);
        gx = (gx < 0.0) ? 0.0 : gx;
        gx = (gx > 1.0) ? 1.0 : gx;
        /* Test canopy resistance */
        gx = 1.0;
        wf->et[k] = ((zsoil[k] - zsoil[k-1]) /
        zsoil[ps->nroot - 1]) * gx * etpa;
        wf->et[k] = ((zsoil[k] - zsoil[k-1]) /
        zsoil[ps->nroot - 1]) * etpa;
        /* Using root distribution as weighting factor */
        wf->et[k] = lc->rtdis[k] * etpa;
        wf->et[k] = etpa * part[k];
    }
# endif
}

void WDfCnd(double *wdf, double *wcnd, double smc, double sicemax,
    const soil_struct *soil)
{
    /*
     * Calculate soil water diffusivity and soil hydraulic conductivity.
     * Flux-PIHM: using van Genuchten parameters
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

    /* Factr2 should avoid to be 0 or 1 */
    factr2 = (factr2 > 1.0 - 5.0E-4) ? 1.0 - 5.0E-4 : factr2;
    factr2 = (factr2 < 0.0 + 5.0E-4) ? 5.0E-4 : factr2;

    factr1 = (factr1 < factr2) ? factr1 : factr2;
    expon = 1.0 - 1.0 / soil->beta;

    satkfunc = KrFunc(soil->beta, factr2);
    dpsidsm =
        (1.0 - expon) / soil->alpha / expon / (soil->smcmax - soil->smcmin) *
        pow(pow(factr2, -1.0 / expon) - 1.0, -expon) *
        pow(factr2, -(1.0 / expon + 1.0));

    *wcnd = soil->ksatv * satkfunc;

    *wdf = *wcnd * dpsidsm;

    if (sicemax > 0.0)
    {
        vkwgt = 1.0 / (1.0 + pow(500.0 * sicemax, 3.0));
        satkfunc = KrFunc(soil->beta, factr1);
        dpsidsm = (1.0 - expon) / soil->alpha / expon /
            (soil->smcmax - soil->smcmin) *
            pow(pow(factr1, -1.0 / expon) - 1.0, -expon) *
            pow(factr1, -(1.0 / expon + 1.0));
        *wdf = vkwgt * *wdf + (1.0 - vkwgt) * dpsidsm * satkfunc * soil->ksatv;
    }
}

void SfcDifOff(phystate_struct *ps, const lc_struct *lc, double t1v,
    double th2v, int iz0tlnd)
{
    /*
     * Calculate surface layer exchange coefficients via iterative process.
     * See Chen et al. (1997, BLM)
     *
     * This routine sfcdif can handle both over open water (sea, ocean) and
     * over solid surface (land, sea-ice).
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

    ilech = 0;

    if (iz0tlnd == 0 || lc->isurban)
    {
        /* czil: constant C in Zilitinkevich, S. S.1995 */
        zilfc = -ps->czil * VKRM * SQVISC;
    }
    else
    {
        /* Modify czil according to Chen & Zhang, 2009
         * czil = 10 ** -0.40 h, ( where h = 10*zo ) */
        ps->czil = pow(10.0, -0.4 * (ps->z0 / 0.07));
        zilfc = -ps->czil * VKRM * SQVISC;
    }

    zu = ps->z0;
    rdz = 1.0 / ps->zlvl_wind;
    cxch = EXCM * rdz;
    dthv = th2v - t1v;

    /* Beljars correction of ustar */
    du2 = ps->sfcspd * ps->sfcspd;
    du2 = (du2 > EPSU2) ? du2 : EPSU2;

    btgh = btg * HPBL;
    /* If statements to avoid tangent linear problems near zero */
    wstar2 = (btgh * ps->ch * dthv != 0.0) ?
        wwst2 * pow(fabs(btgh * ps->ch * dthv), 2.0 / 3.0) : 0.0;

    ustar = sqrt(ps->cm * sqrt(du2 + wstar2));
    ustar = (ustar > EPSUST) ? ustar : EPSUST;

    /* Zilitinkevitch approach for zt */
    zt = exp(zilfc * sqrt(ustar * ps->z0)) * ps->z0;
    zslu = ps->zlvl_wind + zu;

    zslt = ps->zlvl + zt;

    rlogu = log(zslu / zu);

    rlogt = log(zslt / zt);

    rlmo = elfc * ps->ch * dthv / pow(ustar, 3);

    for (itr = 0; itr < ITRMX; itr++)
    {
        zetalt = zslt * rlmo;
        zetalt = (zetalt > ZTMIN) ? zetalt : ZTMIN;
        rlmo = zetalt / zslt;
        zetalu = zslu * rlmo;
        zetau = zu * rlmo;
        zetat = zt * rlmo;

        /* 1. Monin-Obukkhov length-scale */
        if (ilech == 0)
        {
            if (rlmo < 0.0)
            {
                xlu4 = 1.0 - 16.0 * zetalu;
                xlt4 = 1.0 - 16.0 * zetalt;
                xu4 = 1.0 - 16.0 * zetau;
                xt4 = 1.0 - 16.0 * zetat;

                xlu = sqrt(sqrt(xlu4));
                xlt = sqrt(sqrt(xlt4));
                xu = sqrt(sqrt(xu4));
                xt = sqrt(sqrt(xt4));

                psmz = Pspmu(xu);
                simm = Pspmu(xlu) - psmz + rlogu;
                pshz = Psphu(xt);
                simh = Psphu(xlt) - pshz + rlogt;
            }
            else
            {
                zetalu = (zetalu < ZTMAX) ? zetalu : ZTMAX;
                zetalt = (zetalt < ZTMAX) ? zetalt : ZTMAX;

                psmz = Pspms(zetau);
                simm = Pspms(zetalu) - psmz + rlogu;
                pshz = Psphs(zetat);
                simh = Psphs(zetalt) - pshz + rlogt;
            }
        }
        else
        {
            /* Lech's functions */
            if (rlmo < 0.0)
            {
                psmz = Pslmu(zetau);
                simm = Pslmu(zetalu) - psmz + rlogu;
                pshz = Pslhu(zetat);
                simh = Pslhu(zetalt) - pshz + rlogt;
            }
            else
            {
                zetalu = zetalu < ZTMAX ? zetalu : ZTMAX;
                zetalt = zetalt < ZTMAX ? zetalt : ZTMAX;

                psmz = Pslms(zetau);
                simm = Pslms(zetalu) - psmz + rlogu;
                pshz = Pslhs(zetat);
                simh = Pslhs(zetalt) - pshz + rlogt;
            }
        }

        /* Beljaars correction for ustar */
        ustar = sqrt(ps->cm * sqrt(du2 + wstar2));
        ustar = (ustar > EPSUST) ? ustar : EPSUST;

        /* Zilitinkevitch fix for zt */
        zt = exp(zilfc * sqrt(ustar * ps->z0)) * ps->z0;
        zslt = ps->zlvl + zt;

        rlogt = log(zslt / zt);
        ustark = ustar * VKRM;

        ps->cm = (ustark / simm > cxch) ? ustark / simm : cxch;
        ps->ch = (ustark / simh > cxch) ? ustark / simh : cxch;

        /* If statements to avoid tangent linear problems near zero */
        wstar2 = (btgh * ps->ch * dthv != 0.0) ?
            wwst2 * pow(fabs(btgh * ps->ch * dthv), 2.0 / 3.0) : 0.0;

        rlmn = elfc * ps->ch * dthv / pow(ustar, 3.0);

        rlma = rlmo * WOLD + rlmn * wnew;

        rlmo = rlma;
    }
}

void AdjSmProf(const soil_struct *soil, const phystate_struct *ps,
    const double *sice, double dt, wflux_struct *wf, wstate_struct *ws)
{
    /* Min allowable value of smc will be SH2OMIN. */
    /* In Flux-PIHM, the soil layers are gone thru twice:
     * 1. From bottom to top, to make sure all layers below water table is
     * saturated;
     * 2. From top to bottom, to make sure soil moisture from all layers are
     * within plausible ranges */
    int             k;
    double          ddz;
    double          sh2omid[MAXLYR];
    double          wplus;
    double          stot;
    double          stotmin;

    /* Runoff3: runoff within soil layers */
    wplus = 0.0;
    wf->runoff3 = 0.0;

    for (k = ps->nlayers - 1; k >= 0; k--)
    {
        ddz = (k != 0) ? ps->zsoil[k - 1] - ps->zsoil[k] : -ps->zsoil[0];

        sh2omid[k] = ws->swc[k] + wplus / ddz;
        stot = sh2omid[k] + sice[k];

        if (stot > soil->smcmax)
        {
            ws->smc[k] = soil->smcmax;
            wplus = (stot - soil->smcmax) * ddz;
        }
        else
        {
            stotmin = (ps->satdpth[k] * soil->smcmax +
                (ps->soil_depth[k] - ps->satdpth[k]) * (soil->smcmin + SH2OMIN)) /
                ps->soil_depth[k];
            stotmin = (stotmin > soil->smcmax) ? soil->smcmax : stotmin;
            stotmin = (stotmin < soil->smcmin + SH2OMIN) ?
                (soil->smcmin + SH2OMIN) : stotmin;

            if (stot < stotmin)
            {
                ws->smc[k] = stotmin;
                wplus = (stot - stotmin) * ddz;
            }
            else
            {
                ws->smc[k] = stot;
                wplus = 0.0;
            }
        }

        sh2omid[k] = ws->smc[k] - sice[k];
    }

    for (k = 0; k < ps->nlayers; k++)
    {
        ddz = (k == 0) ? -ps->zsoil[0] : ps->zsoil[k - 1] - ps->zsoil[k];

        ws->swc[k] = sh2omid[k] + wplus / ddz;

        stot = ws->swc[k] + sice[k];
        wplus = (stot > soil->smcmax) ? (stot - soil->smcmax) * ddz : 0.0;

        ws->smc[k] = (stot < soil->smcmax) ? stot : soil->smcmax;
        ws->smc[k] = (ws->smc[k] > soil->smcmin + SH2OMIN) ?
            ws->smc[k] : (soil->smcmin + SH2OMIN);
        ws->swc[k] = ws->smc[k] - sice[k];
        ws->swc[k] = (ws->swc[k] > 0.0) ? ws->swc[k] : 0.0;
    }

    wf->runoff3 = wplus / dt;
}

/* Lech's surface functions */
double Pslmu(double zz)
{
    return -0.96 * log(1.0 - 4.5 * zz);
}

double Pslms(double zz)
{
    const double    RIC = 0.183;
    double          rric;

    rric = 1.0 / RIC;
    return zz * rric - 2.076 * (1.0 - 1.0 / (zz + 1.0));
}

double Pslhu(double zz)
{
    return -0.96 * log(1.0 - 4.5 * zz);
}

double Pslhs(double zz)
{
    const double    RIC = 0.183;
    const double    FHNEU = 0.8;
    const double    RFC = 0.191;
    double          rfac;

    rfac = RIC / (FHNEU * RFC * RFC);
    return zz * rfac - 2.076 * (1.0 - exp(-1.2 * zz));
}

/* Paulson's surface functions */
double Pspmu(double xx)
{
    return -2.0 * log((xx + 1.0) * 0.5) - log((xx * xx + 1.0) * 0.5) +
        2.0 * atan(xx) - PI * 0.5;
}

double Pspms(double yy)
{
    return 5.0 * yy;
}

double Psphu(double xx)
{
    return -2.0 * log((xx * xx + 1.0) * 0.5);
}

double Psphs(double yy)
{
    return 5.0 * yy;
}
