#include "pihm.h"

void DailyBgc(int t, pihm_struct pihm)
{

    int             i;
    double          co2lvl;
    double          dayl, prev_dayl;
    double          ndep, nfix;
    spa_data        spa, prev_spa;
    double         *vwc;
    elem_struct    *elem;
    co2control_struct *co2;
    ndepcontrol_struct *ndepctrl;
    forc_struct    *forc;
    siteinfo_struct *siteinfo;

    elem = &pihm->elem[0];
    co2 = &pihm->co2;
    ndepctrl = &pihm->ndepctrl;
    forc = &pihm->forc;
    siteinfo = &pihm->siteinfo;

# if defined(_OPENMP)
#  pragma omp parallel for
# endif
    for (i = 0; i < nelem; i++)
    {
        // First test for nitrogen balance from previous day
        CheckNitrogenBalance(&elem[i].ns, &elem[i].epv.old_n_balance);
    }

    // BGC module for the current day

    // Get co2 and ndep
    if (spinup_mode)    // Spinup mode
    {
        co2lvl = co2->co2ppm;
        ndep = ndepctrl->ndep / 365.0;
        nfix = ndepctrl->nfix / 365.0;
    }
    else    // Model mode
    {
        // Atmospheric CO2 handling
        if (!(co2->varco2))
        {
            // Constant CO2, constant Ndep
            co2lvl = co2->co2ppm;
            ndep = ndepctrl->ndep / 365.0;
            nfix = ndepctrl->nfix / 365.0;
        }
        else
        {
            co2lvl = GetCO2(t, &forc->co2[0]);
        }

        // Ndep handling
        if (!(ndepctrl->varndep))
        {
            // Constant Ndep
            ndep = ndepctrl->ndep / 365.0;
            nfix = ndepctrl->nfix / 365.0;
        }
        else
        {
            ndep = GetNdep(t, &forc->ndep[0]);
            ndep = ndep / 365.0;
            nfix = ndepctrl->nfix / 365.0;
        }
    }

    // Calculate daylengths
    SunPos(t, siteinfo, &spa);
    SunPos(t - DAYINSEC, siteinfo, &prev_spa);

    dayl = (spa.sunset - spa.sunrise) * 3600.0;
    dayl = (dayl < 0.0) ? (dayl + 24.0 * 3600.0) : dayl;

    prev_dayl = (prev_spa.sunset - prev_spa.sunrise) * 3600.0;
    prev_dayl = (prev_dayl < 0.0) ? (prev_dayl + 24.0 * 3600.0) : prev_dayl;

    // Calculate average soil water content for all model grids
    vwc = (double *)malloc(nelem * sizeof(double));

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             kz;

        vwc[i] = elem[i].daily.avg_sh2o[0] * elem[i].ps.soil_depth[0];
        if (elem[i].ps.nlayers > 1)
        {
            for (kz = 1; kz < elem[i].ps.nlayers; kz++)
            {
                vwc[i] += elem[i].daily.avg_sh2o[kz] *
                    elem[i].ps.soil_depth[kz];
            }
        }
        vwc[i] /= elem[i].soil.depth;
    }

# if defined(_OPENMP)
#  pragma omp parallel for
# endif
    for (i = 0; i < nelem; i++)
    {
        daily_struct   *daily;
        epconst_struct *epc;
        epvar_struct   *epv;
        soil_struct    *soil;
        eflux_struct   *ef;
        phystate_struct *ps;
        cstate_struct  *cs;
        cflux_struct   *cf;
        nstate_struct  *ns;
        nflux_struct   *nf;
        ntemp_struct   *nt;
        solute_struct  *solute;
        psn_struct     *psn_sun, *psn_shade;
        summary_struct *summary;
        int             annual_alloc;

        daily = &elem[i].daily;
        epc = &elem[i].epc;
        epv = &elem[i].epv;
        soil = &elem[i].soil;
        ef = &elem[i].ef;
        ps = &elem[i].ps;
        cs = &elem[i].cs;
        cf = &elem[i].cf;
        ns = &elem[i].ns;
        nf = &elem[i].nf;
        nt = &elem[i].nt;
        solute = &elem[i].solute[0];
        psn_sun = &elem[i].psn_sun;
        psn_shade = &elem[i].psn_shade;
        summary = &elem[i].summary;

        // Determine daylengths
        epv->dayl = dayl;
        epv->prev_dayl = prev_dayl;

        // Determine CO2 level
        ps->co2 = co2lvl;

        PrecisionControl(cs, ns);

        // Zero all the flux variables
        MakeZeroFluxStruct(cf, nf);

        // Phenology fluxes
        Phenology(epc, daily, cs, ns, epv, cf, nf);

        // Test for the annual allocation day
        annual_alloc = (epv->offset_flag == 1 && epv->offset_counter == 1) ? 1 : 0;

        // Calculate leaf area index, sun and shade fractions, and specific leaf area for sun and shade canopy
        // fractions, then calculate canopy radiation interception and transmission
        RadTrans(cs, epc, daily, epv, ps, ef);

        SoilPsi(soil, vwc[i], &epv->psi);

        // Maintenance respiration
        MaintResp(epc, daily, cs, ns, epv, cf);

        // Begin canopy bio-physical process simulation
        if (cs->leafc && epv->dayl)
        {
            // Conductance
            CanopyCond(soil, epc, daily, ps, ef, epv);
        }

        // Do photosynthesis only when it is part of the current growth season, as defined by the remdays_curgrowth
        // flag. This keeps the occurrence of new growth consistent with the treatment of litterfall and allocation
        if (cs->leafc && !epv->dormant_flag && epv->dayl)
        {
            TotalPhotosynthesis(epc, daily, ps, epv, cf, psn_sun, psn_shade);
        }
        else
        {
            epv->assim_sun = epv->assim_shade = 0.0;
        }

        nf->ndep_to_sminn = ndep;
        nf->nfix_to_sminn = nfix;

        // Daily litter and soil decomp and nitrogen fluxes
        Decomp(daily->avg_stc[0] - TFREEZ, epc, cs, ns, epv, cf, nf, nt);

        // Allocation gets called whether or not this is a current growth day, because the competition between decomp
        // immobilization fluxes and plant growth N demand is resolved here. On days with no growth, no allocation
        // occurs, but immobilization fluxes are updated normally.
        DailyAllocation(epc, cs, ns, epv, cf, nf, nt);

        // Growth respiration
        GrowthResp(epc, cf);

        // Update of carbon state variables
        DailyCarbonStateUpdate(annual_alloc, epc->woody, epc->evergreen, cs, cf);

        // Update of nitrogen state variables
        DailyNitrogenStateUpdate(annual_alloc, epc->woody, epc->evergreen, ns, nf, solute);

        // Calculate mortality fluxes and update state variables
        // This is done last, with a special state update procedure, to ensure that pools don't go negative due to
        // mortality fluxes conflicting with other proportional fluxes
        Mortality(epc, cs, cf, ns, nf);

        // Test for carbon balance
        CheckCarbonBalance(cs, &epv->old_c_balance);

#if OBSOLETE
        // Nitrogen balance is checked the next day at the beginning of DailyBgc function because a bgc cycle is not
        // finished until N state variables are solved by CVODE
        CheckNitrogenBalance (ns, &epv->old_n_balance);
#endif

        // Calculate carbon summary variables
        CSummary(cs, cf, summary);

        if (spinup_mode)
        {
            elem[i].spinup.soilc += summary->soilc;
            elem[i].spinup.totalc += summary->totalc;
        }
    }

    first_balance = 0;

    free(vwc);
}
