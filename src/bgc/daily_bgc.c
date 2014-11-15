/*
bgc.c
Core BGC model logic

Includes in-line output handling routines that write to daily and annual
output files. This is the only library module that has external
I/O connections, and so it is the only module that includes bgc_io.h.

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information

Revisions since 4.1.2
	Merged spinup_bgc.c with bgc.c to eliminate
	code duplication
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

void daily_bgc(bgc_struct BGCM, bgc_grid *grid, const double t, const double naddfrac, int first_balance)
{
    siteconst_struct *sitec;
    metvar_struct   *metv;
    co2control_struct *co2;
    ndepcontrol_struct *ndepctrl;
    control_struct  *ctrl;    
    epconst_struct  *epc;
    epvar_struct    *epv;
    psn_struct      *psn_sun, *psn_shade;
    wstate_struct *ws;
    wflux_struct *wf;
    cstate_struct *cs;
    cflux_struct *cf;
    nstate_struct *ns;
    nflux_struct *nf;
    ntemp_struct *nt;
    phenology_struct *phen;
    summary_struct *summary;
    struct tm      *timestamp;
    time_t          *rawtime;

    /* miscelaneous variables for program control in main */
    int simyr, yday, metyr, metday;
    int annual_alloc;
    int outv;
    int i, nmetdays;
    double tair_avg, tdiff;
    int dayout;

    double daily_ndep, daily_nfix, ndep_scalar, ndep_diff, ndep;
    int ind_simyr;

    sitec = &grid->sitec;
    metv = &grid->metv;
    co2 = &BGCM->co2;
    ndepctrl = &BGCM->ndepctrl;
    ctrl = &BGCM->ctrl;
    epc = &grid->epc;
    epv = &grid->epv;
    ws = &grid->ws;
    wf = &grid->wf;
    cs = &grid->cs;
    cf = &grid->cf;
    ns = &grid->ns;
    nf = &grid->nf;
    nt = &grid->nt;
    phen = &grid->phen;
    psn_sun = &grid->psn_sun;
    psn_shade = &grid->psn_shade;
    summary = &grid->summary;
//    printf ("BGC daily cycle\n");
    rawtime = (time_t *) malloc (sizeof (time_t));
    *rawtime = (int)t;
    timestamp = gmtime (rawtime);

//    printf ("BGC Time = %4.4d-%2.2d-%2.2d %2.2d:%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);

    /* Get co2 and ndep */
    if (ctrl->spinup == 1)      /* Spinup mode */
    {
        metv->co2 = co2->co2ppm;
        daily_ndep = ndepctrl->ndep / 365.0;
        daily_nfix = ndepctrl->nfix / 365.0;
    }
    else                /* Model mode */
    {
        /* atmospheric CO2 and Ndep handling */
        if (!(co2->varco2))
        {
            /* constant CO2, constant Ndep */
            metv->co2 = co2->co2ppm;
            daily_ndep = ndepctrl->ndep / 365.0;
            daily_nfix = ndepctrl->nfix / 365.0;
        }
        else 
        {
            /* when varco2 = 1, use file for co2 */
            if (co2->varco2 == 1)
                metv->co2 = get_co2(BGCM->Forcing[CO2_TS][0], t);
            if (metv->co2 < -999)
            {
                printf ("Error finding CO2 value on %4.4d-%2.2d-%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday);
                exit (1);
            }

            /* when varco2 = 2, use the constant CO2 value, but can vary Ndep */
            if (co2->varco2 == 2)
                metv->co2 = co2->co2ppm;

            if (ndepctrl->varndep == 0)
            {
                /* increasing CO2, constant Ndep */
                daily_ndep = ndepctrl->ndep / 365.0;
                daily_nfix = ndepctrl->nfix / 365.0;
            }
            else
            {
                daily_ndep = get_ndep(BGCM->Forcing[NDEP_TS][0], t);
                daily_nfix = ndepctrl->nfix / 365.0;
                if(daily_ndep < -999)
                {
                    printf("Error finding NDEP %4.4d-%2.2d-%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday);
                    exit (1);
                }
                else
                {
                    daily_ndep = daily_ndep / 365.0;	
                }
            }
        }
    }
//    printf ("co2 = %lf, ndep = %lf, nfix = %lf\n", metv->co2, daily_ndep, daily_nfix);

    precision_control (ws, cs, ns);

    /* zero all the daily flux variables */
    make_zero_flux_struct (wf, cf, nf);

    /* test for the annual allocation day */
    if (phen->remdays_litfall == 1)
    {
        annual_alloc = 1;
        printf ("Annual allocation\n");
    }
    else
        annual_alloc = 0;

    /* phenology fluxes */
    phenology (epc, metv, phen, epv, cs, cf, ns, nf);

    /* calculate leaf area index, sun and shade fractions, and specific
       leaf area for sun and shade canopy fractions, then calculate
       canopy radiation interception and transmission */
    radtrans(cs, epc, metv, epv, sitec->sw_alb);

    /* update the ann max LAI for annual diagnostic output */
    if (epv->proj_lai > epv->ytd_maxplai)
        epv->ytd_maxplai = epv->proj_lai;
    
    /* soil water potential */
    epv->vwc = metv->swc;
    soilpsi (sitec, epv->vwc, &epv->psi);

    /* daily maintenance respiration */
    maint_resp (cs, ns, epc, metv, cf, epv);

    /* begin canopy bio-physical process simulation */
    if (cs->leafc && metv->dayl)
    {
        /* conductance */
        canopy_et (metv, epc, epv, wf);
    }
    /* do photosynthesis only when it is part of the current growth season, as
     * defined by the remdays_curgrowth flag.  This keeps the occurrence of
     * new growth consistent with the treatment of litterfall and allocation */
    if (cs->leafc && !epv->dormant_flag && metv->dayl)
        total_photosynthesis (metv, epc, epv, cf, psn_sun, psn_shade);
    else
        epv->assim_sun = epv->assim_shade = 0.0;

    nf->ndep_to_sminn = daily_ndep;
    nf->nfix_to_sminn = daily_nfix;
    
    /* daily litter and soil decomp and nitrogen fluxes */
    decomp (metv->tsoil, epc, epv, cs, cf, ns, nf, nt);


    /* Daily allocation gets called whether or not this is a current growth
     * day, because the competition between decomp immobilization fluxes and
     * plant growth N demand is resolved here.  On days with no growth, no
     * allocation occurs, but immobilization fluxes are updated normally */
    daily_allocation (cf, cs, nf, ns, epc, epv, nt, naddfrac, ctrl->spinup);

    /* reassess the annual turnover rates for livewood --> deadwood, and for
     * evergreen leaf and fine root litterfall. This happens once each year,
     * on the annual_alloc day (the last litterfall day) */
    if (annual_alloc)
        annual_rates(epc, epv);

    /* daily growth respiration */
    growth_resp (epc, cf);

    /* daily update of carbon state variables */
    daily_carbon_state_update (cf, cs, annual_alloc, epc->woody, epc->evergreen);

    /* daily update of nitrogen state variables */
    daily_nitrogen_state_update (nf, ns, annual_alloc, epc->woody, epc->evergreen);

    /* calculate N leaching loss.  This is a special state variable update
     * routine, done after the other fluxes and states are reconciled in order
     * to avoid negative sminn under heavy leaching potential */
//    nleaching(ns, nf, ws, wf);

    /* calculate daily mortality fluxes and update state variables */
    /* this is done last, with a special state update procedure, to insure
     * that pools don't go negative due to mortality fluxes conflicting with
     * other proportional fluxes */
    mortality (epc, cs, cf, ns, nf);

    /* test for carbon balance */
    check_carbon_balance(cs, first_balance);

    /* test for nitrogen balance */
    check_nitrogen_balance(ns, first_balance);

    /* calculate carbon summary variables */
    csummary(cf, cs, summary);
}
