/* 
 * phenology.c
 * daily phenology fluxes
 * 
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 * Biome-BGC version 4.2 (final release)
 * See copyright.txt for Copyright information
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 */

#include "bgc.h"

void phenology (const epconst_struct * epc, const metvar_struct * metv, phenology_struct * phen, epvar_struct * epv, cstate_struct * cs, cflux_struct * cf, nstate_struct * ns, nflux_struct * nf)
{
    int             woody, evergreen;
    double          ndays;
    /* phenology model variables */
    int             ws_flag;
    double          onset_critsum;
    double          critdayl = 39300.0; /* seconds */

    double          leaflitfallc, frootlitfallc;
    double          livestemtovrc, livestemtovrn;
    double          livecroottovrc, livecroottovrn;
    double          drate;

    /* set some local flags to control the phenology model behavior */
    /* woody=1 --> woody veg type        woody=0 --> non-woody veg type */
    /* evergreen=1 --> evergreen type    evergreen=0 --> deciduous type */
    /* south=1 --> southern hemisphere   south=0 --> northern hemisphere */
    woody = epc->woody;
    evergreen = epc->evergreen;

    onset_critsum = exp (4.795 + 0.129 * epv->annavg_t2m);

    /* define the phenology signals for cases in which the phenology signals
     * are constant between years */
    if (evergreen)
    {
        epv->dormant_flag = 0.;
        epv->day_leafc_litfall_increment = cs->leafc * epc->leaf_turnover / 365.;
        epv->day_frootc_litfall_increment = cs->frootc * epc->froot_turnover / 365.;

        /* litterfall happens everyday.  To prevent litterfall from driving
         * pools negative in the case of a very high mortality, fluxes are
         * checked and set to zero when the pools get too small. */

        /* leaf litterfall */
        leaflitfallc = epv->day_leafc_litfall_increment;
        if (leaflitfallc > cs->leafc)
            leaflitfallc = cs->leafc;
        leaf_litfall (epc, leaflitfallc, cf, nf);

        /* fine root litterfall */
        frootlitfallc = epv->day_frootc_litfall_increment;
        if (frootlitfallc > cs->frootc)
            frootlitfallc = cs->frootc;
        froot_litfall (epc, frootlitfallc, cf, nf);

        /* turnover of live wood to dead wood also happens every day, at a
         * rate determined once each year, using the annual maximum livewoody
         * compartment masses and the specified livewood turnover rate */
        if (epc->woody)
        {
            /* turnover from live stem wood to dead stem wood */
            epv->day_livestemc_turnover_increment = cs->livestemc * epc->livewood_turnover / 365.;
            livestemtovrc = epv->day_livestemc_turnover_increment;
            livestemtovrn = livestemtovrc / epc->livewood_cn;
            if (livestemtovrc > cs->livestemc)
                livestemtovrc = cs->livestemc;
            if (livestemtovrn > ns->livestemn)
                livestemtovrn = ns->livestemn;
            if (livestemtovrc && livestemtovrn)
            {
                cf->livestemc_to_deadstemc = livestemtovrc;
                nf->livestemn_to_deadstemn = livestemtovrc / epc->deadwood_cn;
                nf->livestemn_to_retransn = livestemtovrn - nf->livestemn_to_deadstemn;
            }

            /* turnover from live coarse root wood to dead coarse root wood */
            epv->day_livecrootc_turnover_increment = cs->livecrootc * epc->livewood_turnover / 365.;
            livecroottovrc = epv->day_livecrootc_turnover_increment;
            livecroottovrn = livecroottovrc / epc->livewood_cn;
            if (livecroottovrc > cs->livecrootc)
                livecroottovrc = cs->livecrootc;
            if (livecroottovrn > ns->livecrootn)
                livecroottovrn = ns->livecrootn;
            if (livecroottovrc && livecroottovrn)
            {
                cf->livecrootc_to_deadcrootc = livecroottovrc;
                nf->livecrootn_to_deadcrootn = livecroottovrc / epc->deadwood_cn;
                nf->livecrootn_to_retransn = livecroottovrn - nf->livecrootn_to_deadcrootn;
            }
        }

    }
    else
    {
        /* Cases that have variable phenological signals between years */
        /* Use the phenology model described in White et al., 1997 */
        /* the two cases that make it to this block are:
         * model, deciduous, woody   and
         * model, deciduous, non-woody (grass), which are the two cases
         * handled by the White et al. paper */
        if (woody)
        {
            /* Use DECIDUOUS TREE PHENOLOGY MODEL */

            /* set flag for solstice period
             * (winter->summer = 1, summer->winter = 0) */
            if (metv->dayl >= metv->prev_dayl)
                ws_flag = 1;
            else
                ws_flag = 0;

            /* update offset_counter and test for the end of the offset
             * period */
            if (epv->offset_flag == 1.)
            {
                /* decrement counter for offset period */
                epv->offset_counter = epv->offset_counter - 1.;

                /* if this is the end of the offset_period, reset phenology
                 * flags and indices */
                if (epv->offset_counter == 0.)
                {
                    epv->offset_flag = 0.;
                    epv->offset_counter = 0.;
                    epv->dormant_flag = 1.;
                    epv->days_active = 0.;

                    /* reset the previous timestep litterfall flux memory */
                    epv->prev_leafc_to_litter = 0.;
                    epv->prev_frootc_to_litter = 0.;
                }
            }

            /* update onset_counter and test for the end of the onset
             * period */
            if (epv->onset_flag == 1.)
            {
                /* decrement counter for onset period */
                epv->onset_counter = epv->onset_counter - 1.;

                /* if this is the end of the onset period, reset phenology
                 * flags and indices */
                if (epv->onset_counter == 0.0)
                {
                    epv->onset_flag = 0.0;
                    epv->onset_counter = 0.0;
                }
            }

            /* test for switching from dormant period to growth period */
            if (epv->dormant_flag == 1.)
            {

                /* Test to turn on growing degree-day sum, if off.
                 * switch on the growing degree day sum on the winter
                 * solstice */

                if (epv->onset_gddflag == 0. && ws_flag == 1)
                {
                    epv->onset_gddflag = 1.;
                    epv->onset_gdd = 0.;
                }

                /* Test to turn off growing degree-day sum, if on.
                 * This test resets the growing degree day sum if it gets past
                 * the summer solstice without reaching the threshold value.
                 * In that case, it will take until the next winter solstice
                 * before the growing degree-day summation starts again. */

                if (epv->onset_gddflag == 1. && ws_flag == 0)
                {
                    epv->onset_gddflag = 0.;
                    epv->onset_gdd = 0.;
                }

                /* if the gdd flag is set, and if the soil is above freezing
                 * then accumulate growing degree days for onset trigger */

                if (epv->onset_gddflag == 1. && metv->tsoil > 0.)
                    epv->onset_gdd = epv->onset_gdd + metv->tsoil * 1.;

                /* set onset_flag if critical growing degree-day sum is
                 * exceeded */
                if (epv->onset_gdd > onset_critsum)
                {
                    epv->onset_flag = 1.;
                    epv->dormant_flag = 0.;
                    epv->onset_gddflag = 0.;
                    epv->onset_gdd = 0.;
                    epv->onset_counter = epc->transfer_days;
                }
            }

            /* test for switching from growth period to offset period */
            else if (epv->offset_flag == 0.)
            {
                /* only begin to test for offset daylength once past the
                 * summer sol */
                if (ws_flag == 0 && metv->dayl < critdayl)
                {
                    epv->offset_flag = 1.;
                    epv->offset_counter = epc->litfall_days;
                    epv->prev_leafc_to_litter = 0.;
                    epv->prev_frootc_to_litter = 0.;
                }
            }
        }                       /* end if woody (tree phenology model) */
        //        else
        //        {
        //            /* non-woody, use the GRASS PHENOLOGY MODEL to calculate the
        //               array of onset and offset days */
        //            /* loop through the entire tavg timeseries to calculate long-term
        //               average tavg and long-term average annual total precip */
        //            mean_tavg = 0.0;
        //            ann_prcp = 0.0;
        //            for (i=0 ; i<ndays ; i++)
        //            {
        //                mean_tavg += metarr->tavg[i];
        //                ann_prcp += metarr->prcp[i];
        //            }
        //            mean_tavg /= (double)ndays;
        //            ann_prcp /= (double)ndays / 365.0;
        //
        //            /* grass onset equation from White et al., 1997, with parameter
        //               values specified by Mike White, Aug. 1997 */
        //            t1 = exp(grass_a * (mean_tavg - grass_tmid));
        //            grass_stsumcrit = ((grass_stsummax - grass_stsummin)* 0.5 *
        //                    ((t1-1)/(t1+1))) + grass_stsummid;
        //            grass_prcpcrit = ann_prcp * grass_k;
        //
        //            /* now go through the phenological years and generate onset
        //               and offset days */
        //
        //            /* calculate the long-term average annual high temperature
        //               for use in grass offset prediction */
        //            tmax_ann = 0.0;
        //            tmin_annavg = 0.0;
        //            for (py=0 ; py<phenyears ; py++)
        //            {
        //                new_tmax = -1000.0;
        //                for (pday=0 ; pday<365 ; pday++)
        //                {
        //                    if (south)
        //                    {
        //                        if (py==0 && pday<182)
        //                        {
        //                            /* use the end of the first year to fill the 
        //                               beginning of a southern hemisphere phenological
        //                               year */
        //                            tmax = metarr->tmax[183+pday];
        //                            tmin_annavg += metarr->tmin[183+pday];
        //                        }
        //                        else if (py==phenyears-1 && pday>181)
        //                        {
        //                            /* use the beginning of the last year to fill the
        //                               end of the last phenological year */
        //                            tmax = metarr->tmax[ndays-547+pday];
        //                            tmin_annavg += metarr->tmin[ndays-547+pday];
        //                        }
        //                        else
        //                        {
        //                            tmax = metarr->tmax[py*365-182+pday];
        //                            tmin_annavg += metarr->tmin[py*365-182+pday];
        //                        }
        //                    }
        //                    else /* north */
        //                    {
        //                        tmax = metarr->tmax[py*365+pday];
        //                        tmin_annavg += metarr->tmin[py*365+pday];
        //                    }
        //
        //                    if (tmax > new_tmax) new_tmax = tmax;
        //
        //                } /* end pday loop */
        //
        //                tmax_ann += new_tmax;
        //            } /* end py loop */
        //            tmax_ann /= (double) phenyears;
        //            /* 92% of tmax_ann is the threshold used in grass offset below */
        //            tmax_ann *= 0.92;
        //            tmin_annavg /= (double) phenyears * 365.0;
        //
        //            /* loop through phenyears again, fill onset and offset arrays */
        //            for (py=0 ; py<phenyears ; py++)
        //            {
        //                sum_soilt = 0.0;
        //                sum_prcp = 0.0;
        //                onset_day = offset_day = -1;
        //                for (pday=0 ; pday<365 ; pday++)
        //                {
        //                    if (south)
        //                    {
        //                        if (py==0 && pday<182)
        //                        {
        //                            /* use the end of the first year to fill the 
        //                               beginning of a southern hemisphere phenological
        //                               year */
        //                            phensoilt = metarr->tavg_ra[183+pday];
        //                            phenprcp = metarr->prcp[183+pday];
        //                            grass_prcpyear[pday] = phenprcp;
        //                            grass_tminyear[pday] = metarr->tmin[183+pday];
        //                            grass_tmaxyear[pday] = metarr->tmax[183+pday];
        //                        }
        //                        else if (py==phenyears-1 && pday>181)
        //                        {
        //                            /* use the beginning of the last year to fill the
        //                               end of the last phenological year */
        //                            phensoilt = metarr->tavg_ra[ndays-547+pday];
        //                            phenprcp = metarr->prcp[ndays-547+pday];
        //                            grass_prcpyear[pday] = phenprcp;
        //                            grass_tminyear[pday] = metarr->tmin[ndays-547+pday];
        //                            grass_tmaxyear[pday] = metarr->tmax[ndays-547+pday];
        //                        }
        //                        else
        //                        {
        //                            phensoilt = metarr->tavg_ra[py*365-182+pday];
        //                            phenprcp = metarr->prcp[py*365-182+pday];
        //                            grass_prcpyear[pday] = phenprcp;
        //                            grass_tminyear[pday] = metarr->tmin[py*365-182+pday];
        //                            grass_tmaxyear[pday] = metarr->tmax[py*365-182+pday];
        //                        }
        //                    }
        //                    else /* north */
        //                    {
        //                        phensoilt = metarr->tavg_ra[py*365+pday];
        //                        phenprcp = metarr->prcp[py*365+pday];
        //                        grass_prcpyear[pday] = phenprcp;
        //                        grass_tminyear[pday] = metarr->tmin[py*365+pday];
        //                        grass_tmaxyear[pday] = metarr->tmax[py*365+pday];
        //                    }
        //
        //                    /* grass onset test */
        //                    if (onset_day == -1)
        //                    {
        //                        if (phensoilt > 0.0) sum_soilt += phensoilt;
        //                        sum_prcp += phenprcp;
        //                        if (sum_soilt >= grass_stsumcrit &&
        //                                sum_prcp >= grass_prcpcrit) onset_day = pday;
        //                    }
        //
        //                } /* end pday loop */
        //
        //                /* do averaging operations on grass_prcpyear and grass_tminyear,
        //                   and do tests for offset day. Offset due to hot & dry can't
        //                   happen within one month after the onset day, and offset due
        //                   to cold can't happen before midyear (yearday 182) */
        //                if (onset_day != -1)
        //                {
        //                    /* calculate three-day boxcar average of tmin */
        //                    if (boxcar_smooth(grass_tminyear, grass_3daytmin, 365,3,0))
        //                    {
        //                        bgc_printf(BV_ERROR, "Error in prephenology() call to boxcar()\n");
        //                        ok=0;
        //                    }
        //
        //                    for (pday=onset_day+30 ; pday<365 ; pday++)
        //                    {
        //                        /* calculate the previous 31-day prcp total */
        //                        psum_startday = pday - 30;
        //                        grass_prcpprev = 0.0;
        //                        for (i=psum_startday ; i<=pday ; i++)
        //                        {
        //                            grass_prcpprev += grass_prcpyear[i];
        //                        }
        //
        //                        /* calculate the next 7-day prcp total */
        //                        if (pday > 358) psum_stopday = 364;
        //                        else psum_stopday = pday + 6;
        //                        grass_prcpnext = 0.0;
        //                        for (i=pday ; i<=psum_stopday ; i++)
        //                        {
        //                            grass_prcpnext += grass_prcpyear[i];
        //                        }
        //
        //                        /* test for hot and dry conditions */
        //                        if (offset_day == -1)
        //                        {
        //                            if (grass_prcpprev < grass_prcpprevcrit && 
        //                                    grass_prcpnext < grass_prcpnextcrit &&
        //                                    grass_tmaxyear[pday] > tmax_ann)
        //                                offset_day = pday;
        //                        }
        //
        //                        /* test for cold offset condition */
        //                        if (offset_day == -1)
        //                        {
        //                            if (pday > 182 &&
        //                                    grass_3daytmin[pday] <= tmin_annavg)
        //                                offset_day = pday;
        //                        }
        //
        //                    } /* end of pdays loop for grass offset testing */
        //                } /* end of if onset_day != -1 block */
        //
        //                /* now do some exception handling for this year's phenology */
        //                if (onset_day != -1)
        //                {
        //                    /* leaves are turned on sometime this year */
        //                    /* subtract 15 days from onset day to approximate the
        //                       start of the new growth period, instead of the middle of
        //                       the new growth period, as is used in the White et al. ms. */
        //                    if (onset_day >= 15)
        //                    {
        //                        onset_day -= 15;
        //                    }
        //                    else onset_day = 0;
        //
        //                    /* if leaves never got turned off, force off on last day */
        //                    if (offset_day == -1) offset_day = 364;
        //
        //                    /* force onset and offset to be at least one day apart */
        //                    if (onset_day == offset_day)
        //                    {
        //                        if (onset_day > 0) onset_day--;
        //                        else offset_day++;
        //                    }
        //                }
        //                else
        //                {
        //                    /* leaves never got turned on, this is a non-growth
        //                       year.  This probably indicates a problem with the
        //                       phenology model */
        //                    onset_day = -1;
        //                    offset_day = -1;
        //                }
        //
        //                /* save these onset and offset days and go to the next
        //                   phenological year */
        //                onday_arr[py] = onset_day;
        //                offday_arr[py] = offset_day;
        //
        //            } /* end phenyears loop for filling onset and offset arrays */
        //        } /* end else !woody (grass phenology model) */

        /* now the onset and offset days are established for each phenyear,
         * either by the deciduous tree or the grass model.  Next loop through
         * phenyears filling the phenology signal arrays and copying them to 
         * the permanent phen struct arrays */

        /* deciduous */
        /* transfer growth fluxes */
        /* check for days left in transfer growth period */
        /* AAN - yes, this is an assignment */

        phen->remdays_transfer = epv->onset_counter;
        phen->remdays_litfall = epv->offset_counter;

        if ((ndays = phen->remdays_transfer))
        {
            /* transfer rate is defined to be a linearly decreasing
             * function that reaches zero on the last day of the transfer
             * period */
            cf->leafc_transfer_to_leafc = 2.0 * cs->leafc_transfer / ndays;
            nf->leafn_transfer_to_leafn = 2.0 * ns->leafn_transfer / ndays;
            cf->frootc_transfer_to_frootc = 2.0 * cs->frootc_transfer / ndays;
            nf->frootn_transfer_to_frootn = 2.0 * ns->frootn_transfer / ndays;
            if (epc->woody)
            {
                cf->livestemc_transfer_to_livestemc = 2.0 * cs->livestemc_transfer / ndays;
                nf->livestemn_transfer_to_livestemn = 2.0 * ns->livestemn_transfer / ndays;
                cf->deadstemc_transfer_to_deadstemc = 2.0 * cs->deadstemc_transfer / ndays;
                nf->deadstemn_transfer_to_deadstemn = 2.0 * ns->deadstemn_transfer / ndays;
                cf->livecrootc_transfer_to_livecrootc = 2.0 * cs->livecrootc_transfer / ndays;
                nf->livecrootn_transfer_to_livecrootn = 2.0 * ns->livecrootn_transfer / ndays;
                cf->deadcrootc_transfer_to_deadcrootc = 2.0 * cs->deadcrootc_transfer / ndays;
                nf->deadcrootn_transfer_to_deadcrootn = 2.0 * ns->deadcrootn_transfer / ndays;
            }
        }

        /* litterfall */
        /* defined such that all live material is removed by the end of the
         * litterfall period, with a linearly ramping removal rate. assumes that
         * the initial rate on the first day of litterfall is 0.0. */
        /* AAN - yes, this is an assignment */

        if ((ndays = phen->remdays_litfall))
        {
            if (ndays == 1.0)
            {
                /* last day of litterfall, special case to gaurantee
                 * that pools go to 0.0 */
                leaflitfallc = cs->leafc;
                frootlitfallc = cs->frootc;
            }
            else
            {
                /* otherwise, assess litterfall 
                 * rates as described above */
                leaflitfallc = epv->day_leafc_litfall_increment;
                drate = 2.0 * (cs->leafc - leaflitfallc * ndays) / (ndays * ndays);
                //printf ("drate = 2.0 * (%lf - %lf * %lf) / (%lf * %lf)\n", cs->leafc, epv->day_leafc_litfall_increment, ndays, ndays, ndays);
                epv->day_leafc_litfall_increment += drate;
                leaflitfallc = epv->day_leafc_litfall_increment;
                //printf ("days = %lf, phen->remdays_litfall = %lf, drate = %lf, leafc = %lf, litfallc = %lf\n", ndays, phen->remdays_litfall, drate, cs->leafc, leaflitfallc);
                frootlitfallc = epv->day_frootc_litfall_increment;
                drate = 2.0 * (cs->frootc - frootlitfallc * ndays) / (ndays * ndays);
                epv->day_frootc_litfall_increment += drate;
                frootlitfallc = epv->day_frootc_litfall_increment;
            }
            /* leaf litterfall */
            if (leaflitfallc > cs->leafc)
                leaflitfallc = cs->leafc;
            //            if (leaflitfallc)
            //printf ("leafc = %lf, leaflitfallc = %lf\n", cs->leafc, leaflitfallc);
            leaf_litfall (epc, leaflitfallc, cf, nf);
            /* fine root litterfall */
            if (frootlitfallc > cs->frootc)
                frootlitfallc = cs->frootc;
            //            if (frootlitfallc)
            froot_litfall (epc, frootlitfallc, cf, nf);
        }                       /* end if deciduous litterfall day */

        /* turnover of livewood to deadwood happens each day, just as for
         * evergreen types, at a rate determined from the annual maximum
         * livewood mass and the specified turnover rate */
        if (epc->woody)
        {
            /* turnover from live stem wood to dead stem wood */
            epv->day_livestemc_turnover_increment = cs->livestemc * epc->livewood_turnover / 365.;
            livestemtovrc = epv->day_livestemc_turnover_increment;
            livestemtovrn = livestemtovrc / epc->livewood_cn;
            if (livestemtovrc > cs->livestemc)
                livestemtovrc = cs->livestemc;
            if (livestemtovrn > ns->livestemn)
                livestemtovrn = ns->livestemn;
            if (livestemtovrc && livestemtovrn)
            {
                cf->livestemc_to_deadstemc = livestemtovrc;
                nf->livestemn_to_deadstemn = livestemtovrc / epc->deadwood_cn;
                nf->livestemn_to_retransn = livestemtovrn - nf->livestemn_to_deadstemn;
            }

            /* turnover from live coarse root wood to dead coarse root wood */
            epv->day_livecrootc_turnover_increment = cs->livecrootc * epc->livewood_turnover / 365.;
            livecroottovrc = epv->day_livecrootc_turnover_increment;
            livecroottovrn = livecroottovrc / epc->livewood_cn;
            if (livecroottovrc > cs->livecrootc)
                livecroottovrc = cs->livecrootc;
            if (livecroottovrn > ns->livecrootn)
                livecroottovrn = ns->livecrootn;
            if (livecroottovrc && livecroottovrn)
            {
                cf->livecrootc_to_deadcrootc = livecroottovrc;
                nf->livecrootn_to_deadcrootn = livecroottovrc / epc->deadwood_cn;
                nf->livecrootn_to_retransn = livecroottovrn - nf->livecrootn_to_deadcrootn;
            }
        }

    }                           /* end else phenology model block */

        //printf ("evergreen = %d, woody = %d, w_s_flag = %1.1d, onset_flag = %1.0lf, dormant_flag = %1.0lf, tsoil = %lf, onset_gddflag = %1.0lf, onset_gdd = %lf, onset_counter = %lf, remdays_transfer = %lf, remdays_litfall = %lf\n", evergreen, woody, ws_flag, epv->onset_flag, epv->dormant_flag, metv->tsoil, epv->onset_gddflag, epv->onset_gdd, epv->onset_counter, phen->remdays_transfer, phen->remdays_litfall);

    /* for woody types, find annual maximum value for live stemc and live crootc
     * calculation of livewood turnover rates */
    if (epc->woody)
    {
        if (epv->annmax_livestemc < cs->livestemc)
            epv->annmax_livestemc = cs->livestemc;
        if (epv->annmax_livecrootc < cs->livecrootc)
            epv->annmax_livecrootc = cs->livecrootc;
    }

    /* for all types, find annual maximum leafc */
    if (epv->annmax_leafc < cs->leafc)
        epv->annmax_leafc = cs->leafc;
    if (epv->annmax_frootc < cs->frootc)
        epv->annmax_frootc = cs->frootc;
}

int free_phenmem (phenarray_struct * phen)
{
    int             ok = 1;

    /* free memory in phenology arrays */
    free (phen->remdays_curgrowth);
    free (phen->remdays_transfer);
    free (phen->remdays_litfall);
    free (phen->predays_transfer);
    free (phen->predays_litfall);

    return (!ok);
}

//int phenology(const epconst_struct* epc, const phenology_struct* phen, epvar_struct* epv, cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf)
//{
//  int ok=1;
//  

void leaf_litfall (const epconst_struct * epc, double litfallc, cflux_struct * cf, nflux_struct * nf)
{
    double          c1, c2, c3, c4;
    double          n1, n2, n3, n4;
    double          nretrans;
    double          avg_cn;
    double          litfalln;

    avg_cn = epc->leaf_cn;
    litfalln = litfallc / epc->leaflitr_cn;

    c1 = litfallc * epc->leaflitr_flab;
    n1 = litfalln * epc->leaflitr_flab;
    c2 = litfallc * epc->leaflitr_fucel;
    n2 = litfalln * epc->leaflitr_fucel;
    c3 = litfallc * epc->leaflitr_fscel;
    n3 = litfalln * epc->leaflitr_fscel;
    c4 = litfallc * epc->leaflitr_flig;
    n4 = litfalln * epc->leaflitr_flig;
    nretrans = (litfallc / avg_cn) - (litfalln);
    ////printf ("litfallc = %lf, avg_cn = %lf, litfalln = %lf, litr_cn = %lf, nretrans = %lf\n", litfallc, avg_cn, litfalln, epc->leaflitr_cn, nretrans);
    //printf ("litfallc = %lf, c1 = %lf, c2 = %lf, c3 = %lf, c4 = %lf, c1+c2+c3+c4 = %lf\n", litfallc, c1, c2, c3, c4, c1 + c2 + c3 + c4);

    /* set fluxes in daily flux structure */
    cf->leafc_to_litr1c = c1;
    cf->leafc_to_litr2c = c2;
    cf->leafc_to_litr3c = c3;
    cf->leafc_to_litr4c = c4;
    nf->leafn_to_litr1n = n1;
    nf->leafn_to_litr2n = n2;
    nf->leafn_to_litr3n = n3;
    nf->leafn_to_litr4n = n4;
    nf->leafn_to_retransn = nretrans;
}

void froot_litfall (const epconst_struct * epc, double litfallc, cflux_struct * cf, nflux_struct * nf)
{
    double          c1, c2, c3, c4;
    double          n1, n2, n3, n4;
    double          avg_cn;

    avg_cn = epc->froot_cn;

    c1 = litfallc * epc->frootlitr_flab;
    n1 = c1 / avg_cn;
    c2 = litfallc * epc->frootlitr_fucel;
    n2 = c2 / avg_cn;
    c3 = litfallc * epc->frootlitr_fscel;
    n3 = c3 / avg_cn;
    c4 = litfallc * epc->frootlitr_flig;
    n4 = c4 / avg_cn;

    /* set fluxes in daily flux structure */
    cf->frootc_to_litr1c = c1;
    cf->frootc_to_litr2c = c2;
    cf->frootc_to_litr3c = c3;
    cf->frootc_to_litr4c = c4;
    nf->frootn_to_litr1n = n1;
    nf->frootn_to_litr2n = n2;
    nf->frootn_to_litr3n = n3;
    nf->frootn_to_litr4n = n4;
}
