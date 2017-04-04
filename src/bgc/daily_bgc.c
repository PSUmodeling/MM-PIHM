
/*
 * daily_bgc.c
 * Daily BGC model logic
 */

#include "pihm.h"

void DailyBgc (pihm_struct pihm, int t, int simstart)
{
    //metvar_struct  *metv;
    co2control_struct *co2;
    ndepcontrol_struct *ndepctrl;
    ctrl_struct    *ctrl;
    epconst_struct *epc;
    epvar_struct   *epv;
    psn_struct     *psn_sun, *psn_shade;
    daily_struct   *daily;
    cstate_struct  *cs;
    cflux_struct   *cf;
    nstate_struct  *ns;
    nflux_struct   *nf;
    ntemp_struct   *nt;
    eflux_struct   *ef;
    pstate_struct  *ps;
    soil_struct    *soil;
    phenology_struct *phen;
    summary_struct *summary;
    int             i, k;
    int             simday;
    struct tm      *timestamp;
    time_t          rawtime;
    double          co2lvl;
    double          vwc;
    double          droot;

    /* miscelaneous variables for program control in main */
    int             annual_alloc;

    double          daily_ndep, daily_nfix;

    rawtime = (int)t;
    timestamp = gmtime (&rawtime);

    co2 = &pihm->co2;
    ndepctrl = &pihm->ndepctrl;
    ctrl = &pihm->ctrl;

    /* Get co2 and ndep */
    if (ctrl->bgc_spinup)       /* Spinup mode */
    {
        co2lvl = co2->co2ppm;
        daily_ndep = ndepctrl->ndep / 365.0;
        daily_nfix = ndepctrl->nfix / 365.0;
    }
    else                        /* Model mode */
    {
        /* Atmospheric CO2 handling */
        if (!(co2->varco2))
        {
            /* constant CO2, constant Ndep */
            co2lvl = co2->co2ppm;
            daily_ndep = ndepctrl->ndep / 365.0;
            daily_nfix = ndepctrl->nfix / 365.0;
        }
        else
        {
            co2lvl = GetCO2 (pihm->forc.co2[0], t);
            if (co2lvl < -999)
            {
                PIHMprintf (VL_ERROR,
                    "Error finding CO2 value on %4.4d-%2.2d-%2.2d\n",
                    timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                    timestamp->tm_mday);
                PIHMexit (EXIT_FAILURE);
            }
        }

        /* Ndep handling */
        if (!(ndepctrl->varndep))
        {
            /* Constant Ndep */
            daily_ndep = ndepctrl->ndep / 365.0;
            daily_nfix = ndepctrl->nfix / 365.0;
        }
        else
        {
            daily_ndep = GetNdep (pihm->forc.ndep[0], t);
            daily_nfix = ndepctrl->nfix / 365.0;
            if (daily_ndep < -999)
            {
                PIHMprintf (VL_ERROR,
                    "Error finding NDEP %4.4d-%2.2d-%2.2d\n",
                    timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                    timestamp->tm_mday);
                PIHMexit (EXIT_FAILURE);
            }
            else
            {
                daily_ndep = daily_ndep / 365.0;
            }
        }
    }

    simday = (t - simstart) / 24 / 3600 - 1;

    if (ctrl->bgc_spinup)       /* Spinup mode */
    {
        for (i = 0; i < pihm->numele; i++)
        {
            DayMet (&pihm->elem[i].stor, &pihm->elem[i].daily, simday);
        }

        for (i = 0; i < pihm->numriv; i++)
        {
            RiverDayMet (&pihm->riv[i].stor, &pihm->riv[i].daily, simday);
        }
    }

    for (i = 0; i < pihm->numele; i++)
    {
        epc = &pihm->elem[i].epc;
        epv = &pihm->elem[i].epv;
        soil = &pihm->elem[i].soil;
        daily = &pihm->elem[i].daily;
        cs = &pihm->elem[i].cs;
        cf = &pihm->elem[i].cf;
        ns = &pihm->elem[i].ns;
        nf = &pihm->elem[i].nf;
        nt = &pihm->elem[i].nt;
        ef = &pihm->elem[i].ef;
        ps = &pihm->elem[i].ps;
        phen = &pihm->elem[i].phen;
        psn_sun = &pihm->elem[i].psn_sun;
        psn_shade = &pihm->elem[i].psn_shade;
        summary = &pihm->elem[i].summary;

        ps->co2 = co2lvl;

        PrecisionControl (cs, ns);

        /* Zero all the daily flux variables */
        MakeZeroFluxStruct (cf, nf);

        /* Phenology fluxes */
        Phenology (epc, daily, phen, epv, cs, cf, ns, nf);

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
        RadTrans (cs, daily, ef, ps, epc, epv);

        /* Update the ann max LAI for annual diagnostic output */
        if (ps->proj_lai > epv->ytd_maxplai)
        {
            epv->ytd_maxplai = ps->proj_lai;
        }

        /* Soil water potential */
        //vwc = (ws->unsat + ws->gw) / soil->depth * soil->smcmax;
        vwc = daily->avg_sh2o[0] * ps->sldpth[0];
        droot = ps->sldpth[0];

        if (ps->nsoil > 1)
        {
            for (k = 1; k < ps->nsoil; k++)
            {
                vwc += daily->avg_sh2o[k] * ps->sldpth[k];
                droot += ps->sldpth[k];
            }
        }

        vwc /= droot;

        SoilPsi (soil, vwc, &epv->psi);

        /* daily maintenance respiration */
        MaintResp (cs, ns, epc, daily, cf, epv);

        /* Begin canopy bio-physical process simulation */
        if (cs->leafc && daily->dayl)
        {
            /* Conductance */
            CanopyCond (epc, daily, ps, soil, epv);
        }

        /* Do photosynthesis only when it is part of the current growth season, as
         * defined by the remdays_curgrowth flag.  This keeps the occurrence of
         * new growth consistent with the treatment of litterfall and
         * allocation */

        if (cs->leafc && !epv->dormant_flag && daily->dayl)
        {
            TotalPhotosynthesis (epc, daily, ps, epv, cf, psn_sun, psn_shade);
        }
        else
        {
            epv->assim_sun = epv->assim_shade = 0.0;
        }

        nf->ndep_to_sminn = daily_ndep;
        nf->nfix_to_sminn = daily_nfix;

        /* Daily litter and soil decomp and nitrogen fluxes */
        Decomp (daily->avg_stc[0] - TFREEZ, epc, epv, cs, cf, ns, nf, nt);

        /* Daily allocation gets called whether or not this is a current growth
         * day, because the competition between decomp immobilization fluxes and
         * plant growth N demand is resolved here.  On days with no growth, no
         * allocation occurs, but immobilization fluxes are updated normally */
        DailyAllocation (cf, cs, nf, ns, epc, epv, nt, ctrl->bgc_spinup);

        /* Reassess the annual turnover rates for livewood --> deadwood, and for
         * evergreen leaf and fine root litterfall. This happens once each year,
         * on the annual_alloc day (the last litterfall day) */
        if (annual_alloc)
        {
            AnnualRates (epc, epv);
        }

        /* Daily growth respiration */
        GrowthResp (epc, cf);

        /* Daily update of carbon state variables */
        DailyCarbonStateUpdate (cf, cs, annual_alloc, epc->woody,
            epc->evergreen);

        /* Daily update of nitrogen state variables */
        DailyNitrogenStateUpdate (nf, ns, annual_alloc, epc->woody,
            epc->evergreen);
    }

    /* Calculate N leaching loss.  This is a special state variable update
     * routine, done after the other fluxes and states are reconciled in order
     * to avoid negative sminn under heavy leaching potential */
    NLeaching (pihm->elem, pihm->numele, pihm->riv, pihm->numriv);

    for (i = 0; i < pihm->numele; i++)
    {
        epc = &pihm->elem[i].epc;
        epv = &pihm->elem[i].epv;
        cs = &pihm->elem[i].cs;
        cf = &pihm->elem[i].cf;
        ns = &pihm->elem[i].ns;
        nf = &pihm->elem[i].nf;
        summary = &pihm->elem[i].summary;

        /* Calculate daily mortality fluxes and update state variables */
        /* This is done last, with a special state update procedure, to insure
         * that pools don't go negative due to mortality fluxes conflicting with
         * other proportional fluxes */
        Mortality (epc, cs, cf, ns, nf);

        /* Test for carbon balance */
        CheckCarbonBalance (cs, &epv->old_c_balance, first_day);

        /* Test for nitrogen balance */
        CheckNitrogenBalance (ns, &epv->old_n_balance, first_day);

        /* Calculate carbon summary variables */
        CSummary (cf, cs, summary);
    }
}
