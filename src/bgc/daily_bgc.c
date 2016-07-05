/*
 * daily_bgc.c
 * Daily BGC model logic
 */

#include "pihm.h"

void DailyBgc (pihm_struct pihm, int t, int simstart, const double *naddfrac, int first_balance)
{
    //metvar_struct  *metv;
    co2control_struct *co2;
    ndepcontrol_struct *ndepctrl;
    ctrl_struct *ctrl;
    epconst_struct *epc;
    epvar_struct   *epv;
    psn_struct     *psn_sun, *psn_shade;
    wstate_struct  *ws;
    wflux_struct   *wf;
    estate_struct  *es;
    eflux_struct   *ef;
    cstate_struct  *cs;
    cflux_struct   *cf;
    nstate_struct  *ns;
    nflux_struct   *nf;
    ntemp_struct   *nt;
    pstate_struct  *ps;
    soil_struct    *soil;
    lc_struct      *lc;
    phenology_struct *phen;
    summary_struct *summary;
    int             i, k;
    int             simday;
    struct tm      *timestamp;
    time_t          rawtime;
    double          co2lvl;
    double          vwc;

    /* miscelaneous variables for program control in main */
    int             annual_alloc;

    double          daily_ndep, daily_nfix;

    rawtime = (int)t;
    timestamp = gmtime (&rawtime);

    co2 = &pihm->co2;
    ndepctrl = &pihm->ndepctrl;
    ctrl = &pihm->ctrl;

    /* Get co2 and ndep */
    if (ctrl->spinup)      /* Spinup mode */
    {
        co2lvl = co2->co2ppm;
        daily_ndep = ndepctrl->ndep / 365.0;
        daily_nfix = ndepctrl->nfix / 365.0;
    }
    else                        /* Model mode */
    {
        /* atmospheric CO2 and Ndep handling */
        if (!(co2->varco2))
        {
            /* constant CO2, constant Ndep */
            co2lvl = co2->co2ppm;
            daily_ndep = ndepctrl->ndep / 365.0;
            daily_nfix = ndepctrl->nfix / 365.0;
        }
        else
        {
            /* When varco2 = 1, use file for co2 */
            if (co2->varco2 == 1)
                co2lvl = GetCO2 (pihm->forc.co2[0], t);
            if (co2lvl < -999)
            {
                printf ("Error finding CO2 value on %4.4d-%2.2d-%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday);
                exit (1);
            }

            /* When varco2 = 2, use the constant CO2 value, but can vary
             * Ndep */
            if (co2->varco2 == 2)
                co2lvl = co2->co2ppm;

            if (ndepctrl->varndep == 0)
            {
                /* Increasing CO2, constant Ndep */
                daily_ndep = ndepctrl->ndep / 365.0;
                daily_nfix = ndepctrl->nfix / 365.0;
            }
            else
            {
                daily_ndep = GetNdep (pihm->forc.ndep[0], t);
                daily_nfix = ndepctrl->nfix / 365.0;
                if (daily_ndep < -999)
                {
                    printf ("Error finding NDEP %4.4d-%2.2d-%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday);
                    exit (1);
                }
                else
                {
                    daily_ndep = daily_ndep / 365.0;
                }
            }
        }
    }

    simday = (t - simstart) / 24 / 3600 - 1;

    if (ctrl->spinup)      /* Spinup mode */
    {
        for (i = 0; i < pihm->numele; i++)
        {
            daymet (&pihm->elem[i].stor, &pihm->elem[i].daily.ws,
                &pihm->elem[i].daily.wf, &pihm->elem[i].daily.es,
                &pihm->elem[i].daily.ef, &pihm->elem[i].daily.ps,
                simday);
        }

        for (i = 0; i < pihm->numriv; i++)
        {
            daymet (&pihm->riv[i].stor, &pihm->riv[i].daily.ws,
                &pihm->riv[i].daily.wf, &pihm->riv[i].daily.es,
                &pihm->riv[i].daily.ef, &pihm->riv[i].daily.ps,
                simday);
        }
    }

    for (i = 0; i < pihm->numele; i++)
    {
        epc = &pihm->elem[i].epc;
        epv = &pihm->elem[i].epv;
        soil = &pihm->elem[i].soil;
        ws = &pihm->elem[i].daily.ws;
        wf = &pihm->elem[i].daily.wf;
        es = &pihm->elem[i].daily.es;
        ef = &pihm->elem[i].daily.ef;
        ps = &pihm->elem[i].daily.ps;
        cs = &pihm->elem[i].cs;
        cf = &pihm->elem[i].cf;
        ns = &pihm->elem[i].ns;
        nf = &pihm->elem[i].nf;
        nt = &pihm->elem[i].nt;
        phen = &pihm->elem[i].phen;
        psn_sun = &pihm->elem[i].psn_sun;
        psn_shade = &pihm->elem[i].psn_shade;
        summary = &pihm->elem[i].summary;

        ps->co2 = co2lvl;

        PrecisionControl (cs, ns);

        /* Zero all the daily flux variables */
        MakeZeroFluxStruct (cf, nf);

        /* Phenology fluxes */
        Phenology (epc, ps, es, phen, epv, cs, cf, ns, nf);

        /* Test for the annual allocation day */
        if (phen->remdays_litfall == 1)
        {
            annual_alloc = 1;
        }
        else
        {
            annual_alloc = 0;
        }

        /* Calculate leaf area index, sun and shade fractions, and specific
         * leaf area for sun and shade canopy fractions, then calculate
         * canopy radiation interception and transmission */
        RadTrans (cs, ef, ps, epc, epv);

        /* Update the ann max LAI for annual diagnostic output */
        if (epv->proj_lai > epv->ytd_maxplai)
        {
            epv->ytd_maxplai = epv->proj_lai;
        }

        /* Soil water potential */
        vwc = (ws->unsat + ws->gw) / soil->depth * soil->smcmax;

        SoilPsi (soil, vwc, &epv->psi);

        /* daily maintenance respiration */
        MaintResp (cs, ns, epc, es, ps, cf, epv);

        /* Begin canopy bio-physical process simulation */
        if (cs->leafc && ps->dayl)
        {
            /* conductance */
            canopy_et (metv, epc, epv, wf);
        }
        /* Do photosynthesis only when it is part of the current growth season, as
         * defined by the remdays_curgrowth flag.  This keeps the occurrence of
         * new growth consistent with the treatment of litterfall and
         * allocation */

        ////printf ("leafc %lf dormant %lf, dayl %lf, soilc = %lf\n", cs->leafc, epv->dormant_flag, metv->dayl, summary->soilc);

        //if (cs->leafc && !epv->dormant_flag && metv->dayl)
        //    total_photosynthesis (metv, epc, epv, cf, psn_sun, psn_shade);
        //else
        //    epv->assim_sun = epv->assim_shade = 0.0;

        //nf->ndep_to_sminn = daily_ndep;
        //nf->nfix_to_sminn = daily_nfix;

        ///* daily litter and soil decomp and nitrogen fluxes */
        //decomp (metv->tsoil, epc, epv, cs, cf, ns, nf, nt);


        ///* Daily allocation gets called whether or not this is a current growth
        // * day, because the competition between decomp immobilization fluxes and
        // * plant growth N demand is resolved here.  On days with no growth, no
        // * allocation occurs, but immobilization fluxes are updated normally */
        //daily_allocation (cf, cs, nf, ns, epc, epv, nt, naddfrac[i], ctrl->spinup);

        ///* reassess the annual turnover rates for livewood --> deadwood, and for
        // * evergreen leaf and fine root litterfall. This happens once each year,
        // * on the annual_alloc day (the last litterfall day) */
        //if (annual_alloc)
        //    annual_rates (epc, epv);

        ///* daily growth respiration */
        //growth_resp (epc, cf);

        ///* daily update of carbon state variables */
        //daily_carbon_state_update (cf, cs, annual_alloc, epc->woody, epc->evergreen);

        ///* daily update of nitrogen state variables */
        //daily_nitrogen_state_update (nf, ns, annual_alloc, epc->woody, epc->evergreen);
    }



    ///* Calculate N leaching loss.  This is a special state variable update
    // * routine, done after the other fluxes and states are reconciled in order
    // * to avoid negative sminn under heavy leaching potential */
    //nleaching (pihm->elem, numele, pihm->riv, numriv);

    //for (i = 0; i < numele; i++)
    //{
    //    //sitec = &pihm->elem[i].sitec;
    //    metv = &pihm->elem[i].metv;
    //    epc = &pihm->elem[i].epc;
    //    epv = &pihm->elem[i].epv;
    //    ws = &pihm->elem[i].ws;
    //    wf = &pihm->elem[i].wf;
    //    cs = &pihm->elem[i].cs;
    //    cf = &pihm->elem[i].cf;
    //    ns = &pihm->elem[i].ns;
    //    nf = &pihm->elem[i].nf;
    //    nt = &pihm->elem[i].nt;
    //    phen = &pihm->elem[i].phen;
    //    psn_sun = &pihm->elem[i].psn_sun;
    //    psn_shade = &pihm->elem[i].psn_shade;
    //    summary = &pihm->elem[i].summary;

    //    /* Calculate daily mortality fluxes and update state variables */
    //    /* This is done last, with a special state update procedure, to insure
    //     * that pools don't go negative due to mortality fluxes conflicting with
    //     * other proportional fluxes */
    //    mortality (epc, cs, cf, ns, nf);

    //    /* Test for carbon balance */
    //    check_carbon_balance (cs, &epv->old_c_balance, first_balance);

    //    /* Test for nitrogen balance */
    //    check_nitrogen_balance (ns, &epv->old_n_balance, first_balance);

    //    /* Calculate carbon summary variables */
    //    csummary (cf, cs, summary);
    //}
}
