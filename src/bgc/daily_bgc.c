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

void daily_bgc(bgc_struct BGCM, bgc_grid *grid, const double t, const double naddfrac)
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
    struct tm      *timestamp;
    time_t          *rawtime;

    /* miscelaneous variables for program control in main */
    int simyr, yday, metyr, metday;
    int first_balance;
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
//    mortality (epc, cs, cf, ns, nf);
//
//    /* test for carbon balance */
//    if (ok && check_carbon_balance(&cs, first_balance))
//    {
//        bgc_printf(BV_ERROR, "Error in check_carbon_balance() from bgc()\n");
//        bgc_printf(BV_ERROR, "%d\n",metday);
//        ok=0;
//    }
//
//    bgc_printf(BV_DIAG, "%d\t%d\tdone carbon balance\n",simyr,yday);
//
//    /* test for nitrogen balance */
//    if (ok && check_nitrogen_balance(&ns, first_balance))
//    {
//        bgc_printf(BV_ERROR, "Error in check_nitrogen_balance() from bgc()\n");
//        bgc_printf(BV_ERROR, "%d\n",metday);
//        ok=0;
//    }
//
//    bgc_printf(BV_DIAG, "%d\t%d\tdone nitrogen balance\n",simyr,yday);
//
//    /* calculate carbon summary variables */
//    if (ok && csummary(&cf, &cs, &summary))
//    {
//        bgc_printf(BV_ERROR, "Error in csummary() from bgc()\n");
//        ok=0;
//    } 
//
//    bgc_printf(BV_DIAG, "%d\t%d\tdone carbon summary\n",simyr,yday);
//
//    /* calculate water summary variables */
//    if (ok && wsummary(&ws,&wf,&summary))
//    {
//        printf("Error in wsummary() from bgc()\n");
//        ok=0;
//    }
//
//    bgc_printf(BV_DIAG, "%d\t%d\tdone water summary\n", simyr,yday);
//
//    /* DAILY OUTPUT HANDLING */
//    /* fill the daily output array if daily output is requested,
//       or if the monthly or annual average of daily output variables
//       have been requested */
//    bgc_printf(BV_DIAG, "Number of daily outputs: %d\n", ctrl.ndayout);
//    if (ok && dayout)
//    {
//        /* fill the daily output array */
//        for (outv=0 ; outv<ctrl.ndayout ; outv++)
//        {
//            bgc_printf(BV_DIAG, "Outv: %d, ", outv);
//            bgc_printf(BV_DIAG, "DayCode: %d, ", ctrl.daycodes[outv]);
//            bgc_printf(BV_DIAG, "Output: %f\n", *output_map[ctrl.daycodes[outv]]);
//            dayarr[outv] = (float) *output_map[ctrl.daycodes[outv]];
//        }
//    }
//    /* only write daily outputs if requested */
//    if (ok && ctrl.dodaily)
//    {
//        /* write the daily output array to daily output file */
//        if (fwrite(dayarr, sizeof(float), ctrl.ndayout, bgcout->dayout.ptr)
//                != (size_t)ctrl.ndayout)
//        {
//            bgc_printf(BV_ERROR, "Error writing to %s: simyear = %d, simday = %d\n",
//                    bgcout->dayout.name,simyr,yday);
//            ok=0;
//        }
//
//        bgc_printf(BV_DIAG, "%d\t%d\tdone daily output\n",simyr,yday);
//        if(ok && bgcout->bgc_ascii)
//        {	
//
//            output_ascii(dayarr,ctrl.ndayout,bgcout->dayoutascii.ptr);
//
//        }
//
//    }
//    /*******************/
//    /* MONTHLY OUTPUTS */
//    /*******************/
//
//    /* MONTHLY AVERAGE OF DAILY OUTPUT VARIABLES */
//    if (ctrl.domonavg)
//    {
//        /* update the monthly average array */
//        for (outv=0 ; outv<ctrl.ndayout ; outv++)
//        {
//            monavgarr[outv] += dayarr[outv];
//
//            switch (ctrl.daycodes[outv])
//            {
//                /* Leaf area index */
//                case 545:   
//                    if(dayarr[outv] > monmaxlai) monmaxlai = dayarr[outv]; 
//                    break;
//            }
//        }
//
//        /* if this is the last day of the current month, output... */
//        if (yday == endday[curmonth])
//        {
//            /* finish the averages */
//            for (outv=0 ; outv<ctrl.ndayout ; outv++)
//            {
//                if (summary_sanity == SANE)
//                {
//                    switch (ctrl.daycodes[outv])
//                    {
//                        /* Leaf area index */
//                        /* Maximum monthly */
//                        case 545:
//                            monavgarr[outv] = monmaxlai; 
//                            break;
//                            /* Snow water */
//                        case 21:
//                            monavgarr[outv] = dayarr[outv] - eomsnoww; 
//                            eomsnoww = dayarr[outv]; 
//                            break;
//                            /* Soil water content */
//                        case 20:
//                            monavgarr[outv] = dayarr[outv] - eomsoilw;
//                            eomsoilw = dayarr[outv];
//                            break;
//                        default:
//                            monavgarr[outv] /= (float)mondays[curmonth];
//                            break;
//                    }
//                }
//                else 
//                {
//                    monavgarr[outv] /= (float)mondays[curmonth];
//                }
//            }
//
//            /* write to file */
//            if (fwrite(monavgarr, sizeof(float), ctrl.ndayout, bgcout->monavgout.ptr)
//                    != (size_t)ctrl.ndayout)
//            {
//                bgc_printf(BV_ERROR, "Error writing to %s: simyear = %d, simday = %d\n",
//                        bgcout->monavgout.name,simyr,yday);
//                ok=0;
//            }
//
//            if(ok && bgcout->bgc_ascii)
//            {
//                output_ascii(monavgarr,ctrl.ndayout, bgcout->monoutascii.ptr);
//
//            }
//
//            /* reset monthly average variables for next month */
//            for (outv=0 ; outv<ctrl.ndayout ; outv++)
//            {
//                monavgarr[outv] = 0.0;
//                monmaxlai = 0.0;
//                monmaxsnoww = 0.0;
//            }
//
//            /* increment current month counter */
//            curmonth++;
//
//            bgc_printf(BV_DIAG, "%d\t%d\tdone monavg output\n",simyr,yday);
//
//        }
//    }
//
//    /* ANNUAL AVERAGE OF DAILY OUTPUT VARIABLES */
//    if (ctrl.doannavg)
//    {
//        /* update the annual average array */
//        for (outv=0 ; outv<ctrl.ndayout ; outv++)
//        {
//            annavgarr[outv] += dayarr[outv];
//            switch (ctrl.daycodes[outv])
//            {
//                /* Leaf area index */
//                case 545:
//                    if(dayarr[outv] > annmaxplai) annmaxplai = dayarr[outv];
//                    break;
//            }
//        }
//
//        /* if this is the last day of the year, output... */
//        if (yday == 364)
//        {
//            /* finish averages */
//            for (outv=0 ; outv<ctrl.ndayout ; outv++)
//            {
//                if (summary_sanity == SANE)
//                {
//                    switch (ctrl.daycodes[outv])
//                    {
//                        /* Leaf area index*/ 
//                        case 545:
//                            annavgarr[outv] = (float)annmaxplai;
//                            break;
//                        default: 
//                            annavgarr[outv] /= 365.0;
//                            break;
//                    }
//                }
//                else
//                {
//                    annavgarr[outv] /= 365.0;
//                }
//            }
//
//            /* write to file */
//            if (fwrite(annavgarr, sizeof(float), ctrl.ndayout, bgcout->annavgout.ptr)
//                    != (size_t)ctrl.ndayout)
//            {
//                bgc_printf(BV_ERROR, "Error writing to %s: simyear = %d, simday = %d\n",
//                        bgcout->annavgout.name,simyr,yday);
//                ok=0;
//            }
//
//            /* reset annual average variables for next month */
//            for (outv=0 ; outv<ctrl.ndayout ; outv++)
//            {
//                annavgarr[outv] = 0.0;
//                annmaxplai = 0.0;
//            }
//
//            bgc_printf(BV_DIAG, "%d\t%d\tdone annavg output\n",simyr,yday);
//
//        }
//    }
//
//    if (mode == MODE_MODEL)
//    {
//        /* very simple annual summary variables for text file output */
//        if (epv.proj_lai > (double)annmaxlai) annmaxlai = (float)epv.proj_lai;
//        annet += wf.canopyw_evap + wf.snoww_subl + wf.soilw_evap +
//            wf.soilw_trans;
//        annoutflow += wf.soilw_outflow;
//        annnpp += summary.daily_npp * 1000.0;
//        annnbp += summary.daily_nee * 1000.0;
//        annprcp += metv.prcp;
//        anntavg += metv.tavg/365.0;
//    }
//    else if (mode == MODE_SPINUP)
//    {
//        /* spinup control */
//        /* keep a tally of total soil C during successive
//           met cycles for comparison */
//        if (metcycle == 1)
//        {
//            tally1 += summary.soilc;
//            tally1b += summary.totalc;
//        }
//        if (metcycle == 2)
//        {
//            tally2 += summary.soilc;
//            tally2b += summary.totalc;
//        }
//    }
//
//    /* at the end of first day of simulation, turn off the 
//       first_balance switch */
//    if (first_balance) first_balance = 0;
//
//
//    /* Make zero fluxes */
//
//    
//

}
