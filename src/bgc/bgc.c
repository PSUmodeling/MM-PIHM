#include "pihm.h"

void Bgc (pihm_struct pihm, int t, double dt, int first_balance)
{

    int             i;
    double          co2lvl;
    double          ndep, nfix;

    /* Get co2 and ndep */
    if (pihm->ctrl.bgc_spinup)      /* Spinup mode */
    {
        co2lvl = pihm->co2.co2ppm;
        ndep = pihm->ndepctrl.ndep / 365.0 / DAYINSEC;
        nfix = pihm->ndepctrl.nfix / 365.0 / DAYINSEC;
    }
    else                            /* Model mode */
    {
        /* Atmospheric CO2 handling */
        if (!(pihm->co2.varco2))
        {
            /* Constant CO2, constant Ndep */
            co2lvl = pihm->co2.co2ppm;
            ndep = pihm->ndepctrl.ndep / 365.0 / DAYINSEC;
            nfix = pihm->ndepctrl.nfix / 365.0 / DAYINSEC;
        }
        else
        {
            co2lvl = GetCO2 (pihm->forc.co2[0], t);
        }

        /* Ndep handling */
        if (!(pihm->ndepctrl.varndep))
        {
            /* Constant Ndep */
            ndep = pihm->ndepctrl.ndep / 365.0 / DAYINSEC;
            nfix = pihm->ndepctrl.nfix / 365.0 / DAYINSEC;
        }
        else
        {
            ndep = GetNdep (pihm->forc.ndep[0], t);
            ndep = ndep / 365.0 / DAYINSEC;
            nfix = pihm->ndepctrl.nfix / 365.0 / DAYINSEC;
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int         k;
        epconst_struct *epc;
        epvar_struct *epv;
        soil_struct *soil;
        wstate_struct *ws;
        estate_struct *es;
        eflux_struct *ef;
        pstate_struct *ps;
        cstate_struct *cs;
        cflux_struct *cf;
        nstate_struct *ns;
        nflux_struct *nf;
        ntemp_struct *nt;
        psn_struct *psn_sun, *psn_shade;
        summary_struct *summary;
        int         annual_alloc;

        epc = &pihm->elem[i].epc;
        epv = &pihm->elem[i].epv;
        soil = &pihm->elem[i].soil;
        ws = &pihm->elem[i].ws;
        es = &pihm->elem[i].es;
        ef = &pihm->elem[i].ef;
        ps = &pihm->elem[i].ps;
        cs = &pihm->elem[i].cs;
        cf = &pihm->elem[i].cf;
        ns = &pihm->elem[i].ns;
        nf = &pihm->elem[i].nf;
        nt = &pihm->elem[i].nt;
        psn_sun = &pihm->elem[i].psn_sun;
        psn_shade = &pihm->elem[i].psn_shade;
        summary = &pihm->elem[i].summary;

        ps->co2 = co2lvl;

        PrecisionControl (cs, ns);

        /* Zero all the flux variables */
        MakeZeroFluxStruct (cf, nf);

        /* Phenology fluxes */
        Phenology (epc, epv, es, cs, cf, ns, nf, dt);

        /* Test for the annual allocation day */
        if (epv->offset_flag == 1 && epv->offset_counter <= (int)round (dt))
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

        /* Soil water potential */
        SoilPsi (soil, ws->soilm / soil->depth, &epv->psi);

        /* Maintenance respiration */
        MaintResp (epc, epv, es, ef, cs, cf, ns);

        /* Begin canopy bio-physical process simulation */
        if (cs->leafc && ef->soldn)
        {
            /* Conductance */
            CanopyCond (epc, epv, ef, ps);
        }

        /* Do photosynthesis only when it is part of the current growth season, as
         * defined by the remdays_curgrowth flag.  This keeps the occurrence of
         * new growth consistent with the treatment of litterfall and
         * allocation */
        if (cs->leafc && !epv->dormant_flag && ef->soldn > 0.0)
        {
            TotalPhotosynthesis (epc, epv, es, ps, cf, psn_sun, psn_shade);
        }
        else
        {
            epv->assim_sun = epv->assim_shade = 0.0;
        }

        nf->ndep_to_sminn = ndep;
        nf->nfix_to_sminn = nfix;

        /* Daily litter and soil decomp and nitrogen fluxes */
        Decomp (es->stc[0] - TFREEZ, epc, epv, cs, cf, ns, nf, nt, dt);

        /* Allocation gets called whether or not this is a current growth
         * day, because the competition between decomp immobilization fluxes and
         * plant growth N demand is resolved here.  On days with no growth, no
         * allocation occurs, but immobilization fluxes are updated normally */
        Allocation (cf, cs, nf, ns, epc, epv, nt, dt);

        /* Growth respiration */
        GrowthResp (epc, cf);

        /* Update of carbon state variables */
        CarbonStateUpdate (cf, cs, annual_alloc, epc->woody, epc->evergreen,
            dt);

        /* Update of nitrogen state variables */
        NitrogenStateUpdate (nf, ns, annual_alloc, epc->woody, epc->evergreen,
            dt);
    }

    /* Calculate N leaching loss.  This is a special state variable update
     * routine, done after the other fluxes and states are reconciled in order
     * to avoid negative sminn under heavy leaching potential */
    //NTransport (pihm->elem, pihm->riv, DAYINSEC);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        epconst_struct *epc;
        epvar_struct *epv;
        cstate_struct *cs;
        cflux_struct *cf;
        nstate_struct *ns;
        nflux_struct *nf;
        summary_struct *summary;

        epc = &pihm->elem[i].epc;
        epv = &pihm->elem[i].epv;
        cs = &pihm->elem[i].cs;
        cf = &pihm->elem[i].cf;
        ns = &pihm->elem[i].ns;
        nf = &pihm->elem[i].nf;
        summary = &pihm->elem[i].summary;

        /* Calculate mortality fluxes and update state variables */
        /* This is done last, with a special state update procedure, to insure
         * that pools don't go negative due to mortality fluxes conflicting with
         * other proportional fluxes */
        Mortality (epc, cs, cf, ns, nf, dt);

        /* Test for carbon balance */
        CheckCarbonBalance (cs, &epv->old_c_balance, first_balance);

        /* Test for nitrogen balance */
        CheckNitrogenBalance (ns, &epv->old_n_balance, first_balance);

        /* Calculate carbon summary variables */
        CSummary (cf, cs, summary);
    }
}
