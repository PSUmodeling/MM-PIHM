/* 
phenology.c
daily phenology fluxes

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

void phenology(const epconst_struct* epc, const phenology_struct* phen, epvar_struct* epv, cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf)
{
    int woody,evergreen,south;
    double t1;
    char round[80];
    int i,pday,ndays,py;
    int nyears,phenyears;
    int ngrowthdays,ntransferdays,nlitfalldays;
    int onday,offday;
    int counter;
    /* phenology model variables */
    int *onday_arr, *offday_arr;
    int fall_tavg_count;
    int onset_day, offset_day;
    double mean_tavg,fall_tavg;
    double phensoilt,phendayl;
    double onset_critsum, sum_soilt;
    double critdayl = 39300.0; /* seconds */
    /* grass model parameters */
    double ann_prcp;
    double sum_prcp, phenprcp;
    double grass_stsumcrit;
    double grass_prcpcrit;
    double grass_stsummax = 1380.0;
    double grass_stsummid = 900.0;
    double grass_stsummin = 418.0;
    double grass_a = 32.9;
    double grass_k = 0.15;
    double grass_tmid = 9.0;
    double grass_prcpyear[365];
    double grass_prcpprevcrit = 1.14;
    double grass_prcpprev;
    double grass_prcpnextcrit = 0.97;
    double grass_prcpnext;
    double grass_tmaxyear[365];
    double grass_tminyear[365];
    double grass_3daytmin[365];
    int psum_startday, psum_stopday;
    double tmax_ann, tmax, new_tmax;
    double tmin_annavg;

    /* set some local flags to control the phenology model behavior */
    /* woody=1 --> woody veg type        woody=0 --> non-woody veg type */
    /* evergreen=1 --> evergreen type    evergreen=0 --> deciduous type */
    /* south=1 --> southern hemisphere   south=0 --> northern hemisphere */
    woody = epc->woody;
    evergreen = epc->evergreen;

    onset_critsum = exp(4.795 + 0.129 * epv->annavg_t2m);

    /* define the phenology signals for cases in which the phenology signals
       are constant between years */
    if (evergreen)
    {
        epv->bglfr = 1. / (epc->leaf_turnover * 365.);
        epv->bgtr = 0.;
        epv->lgsf = 0.;
    } 
    else
    {
        /* Cases that have variable phenological signals between years */
        /* Use the phenology model described in White et al., 1997 */
        /* the two cases that make it to this block are:
           model, deciduous, woody   and
           model, deciduous, non-woody (grass), which are the two cases
           handled by the White et al. paper */
        if (woody)
        {
            /* Use DECIDUOUS TREE PHENOLOGY MODEL */
            /* set background litterfall rate, background transfer rate,
             * and long growing season factor to 0 for seasonal deciduous
             * types */
            epv->bglfr = 0.;
            epv->bgtr = 0.;
            epv->lgsf = 0.;

            /* set flag for solstice period
             * (winter->summer = 1, summer->winter = 0) */
            if (metv->dayl >= epv->prev_dayl)
                ws_flag = 1.;
            else
                ws_flag = 0.;

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
//            if (onset_counter(p) == 0.0_r8) then
//               ! this code block was originally handled by call cn_onset_cleanup(p)
//               ! inlined during vectorization
//
//               onset_flag(p) = 0.0_r8
//               onset_counter(p) = 0.0_r8
//               ! set all transfer growth rates to 0.0
//               leafc_xfer_to_leafc(p)   = 0.0_r8
//               frootc_xfer_to_frootc(p) = 0.0_r8
//               leafn_xfer_to_leafn(p)   = 0.0_r8
//               frootn_xfer_to_frootn(p) = 0.0_r8
//               if (woody(ivt(p)) == 1.0_r8) then
//                  livestemc_xfer_to_livestemc(p)   = 0.0_r8
//                  deadstemc_xfer_to_deadstemc(p)   = 0.0_r8
//                  livecrootc_xfer_to_livecrootc(p) = 0.0_r8
//                  deadcrootc_xfer_to_deadcrootc(p) = 0.0_r8
//                  livestemn_xfer_to_livestemn(p)   = 0.0_r8
//                  deadstemn_xfer_to_deadstemn(p)   = 0.0_r8
//                  livecrootn_xfer_to_livecrootn(p) = 0.0_r8
//                  deadcrootn_xfer_to_deadcrootn(p) = 0.0_r8
//               end if
//               ! set transfer pools to 0.0
//               leafc_xfer(p) = 0.0_r8
//               leafn_xfer(p) = 0.0_r8
//               frootc_xfer(p) = 0.0_r8
//               frootn_xfer(p) = 0.0_r8
//               if (woody(ivt(p)) == 1.0_r8) then
//                  livestemc_xfer(p) = 0.0_r8
//                  livestemn_xfer(p) = 0.0_r8
//                  deadstemc_xfer(p) = 0.0_r8
//                  deadstemn_xfer(p) = 0.0_r8
//                  livecrootc_xfer(p) = 0.0_r8
//                  livecrootn_xfer(p) = 0.0_r8
//                  deadcrootc_xfer(p) = 0.0_r8
//                  deadcrootn_xfer(p) = 0.0_r8
//               end if
//            end if
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

                if (epv->onset_gddflag == 1. && ws_flag == 0.)
                {
                    epv->onset_gddflag = 0.;
                    epv->onset_gdd = 0.;
                }

                /* if the gdd flag is set, and if the soil is above freezing
                 * then accumulate growing degree days for onset trigger */

                if (epv->onset_gddflag == 1. && metv->tsoil > 0.)
                    epv->onset_gdd = onset_gdd + metv->tsoil * 1.;

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

//               ! move all the storage pools into transfer pools,
//               ! where they will be transfered to displayed growth over the onset period.
//               ! this code was originally handled with call cn_storage_to_xfer(p)
//               ! inlined during vectorization
//
//               ! set carbon fluxes for shifting storage pools to transfer pools
//               leafc_storage_to_xfer(p)  = fstor2tran * leafc_storage(p)/dt
//               frootc_storage_to_xfer(p) = fstor2tran * frootc_storage(p)/dt
//               if (woody(ivt(p)) == 1.0_r8) then
//                  livestemc_storage_to_xfer(p)  = fstor2tran * livestemc_storage(p)/dt
//                  deadstemc_storage_to_xfer(p)  = fstor2tran * deadstemc_storage(p)/dt
//                  livecrootc_storage_to_xfer(p) = fstor2tran * livecrootc_storage(p)/dt
//                  deadcrootc_storage_to_xfer(p) = fstor2tran * deadcrootc_storage(p)/dt
//                  gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt
//               end if
//
//               ! set nitrogen fluxes for shifting storage pools to transfer pools
//               leafn_storage_to_xfer(p)  = fstor2tran * leafn_storage(p)/dt
//               frootn_storage_to_xfer(p) = fstor2tran * frootn_storage(p)/dt
//               if (woody(ivt(p)) == 1.0_r8) then
//                  livestemn_storage_to_xfer(p)  = fstor2tran * livestemn_storage(p)/dt
//                  deadstemn_storage_to_xfer(p)  = fstor2tran * deadstemn_storage(p)/dt
//                  livecrootn_storage_to_xfer(p) = fstor2tran * livecrootn_storage(p)/dt
//                  deadcrootn_storage_to_xfer(p) = fstor2tran * deadcrootn_storage(p)/dt
//               end if
            }

            /* test for switching from growth period to offset period */
            else if (epv->offset_flag == 0.)
            {
                /* only begin to test for offset daylength once past the
                 * summer sol */
                if (ws_flag == 0. && metv->dayl < critdayl)
                {
                    epv->offset_flag = 1.;
                    epv->offset_counter = epc->litfall_days;
                    epv->prev_leafc_to_litter = 0.;
                    epv->prev_frootc_to_litter = 0.;
                }
            }
        } /* end if woody (tree phenology model) */
        else
        {
            /* non-woody, use the GRASS PHENOLOGY MODEL to calculate the
               array of onset and offset days */
            /* loop through the entire tavg timeseries to calculate long-term
               average tavg and long-term average annual total precip */
            mean_tavg = 0.0;
            ann_prcp = 0.0;
            for (i=0 ; i<ndays ; i++)
            {
                mean_tavg += metarr->tavg[i];
                ann_prcp += metarr->prcp[i];
            }
            mean_tavg /= (double)ndays;
            ann_prcp /= (double)ndays / 365.0;

            /* grass onset equation from White et al., 1997, with parameter
               values specified by Mike White, Aug. 1997 */
            t1 = exp(grass_a * (mean_tavg - grass_tmid));
            grass_stsumcrit = ((grass_stsummax - grass_stsummin)* 0.5 *
                    ((t1-1)/(t1+1))) + grass_stsummid;
            grass_prcpcrit = ann_prcp * grass_k;

            /* now go through the phenological years and generate onset
               and offset days */

            /* calculate the long-term average annual high temperature
               for use in grass offset prediction */
            tmax_ann = 0.0;
            tmin_annavg = 0.0;
            for (py=0 ; py<phenyears ; py++)
            {
                new_tmax = -1000.0;
                for (pday=0 ; pday<365 ; pday++)
                {
                    if (south)
                    {
                        if (py==0 && pday<182)
                        {
                            /* use the end of the first year to fill the 
                               beginning of a southern hemisphere phenological
                               year */
                            tmax = metarr->tmax[183+pday];
                            tmin_annavg += metarr->tmin[183+pday];
                        }
                        else if (py==phenyears-1 && pday>181)
                        {
                            /* use the beginning of the last year to fill the
                               end of the last phenological year */
                            tmax = metarr->tmax[ndays-547+pday];
                            tmin_annavg += metarr->tmin[ndays-547+pday];
                        }
                        else
                        {
                            tmax = metarr->tmax[py*365-182+pday];
                            tmin_annavg += metarr->tmin[py*365-182+pday];
                        }
                    }
                    else /* north */
                    {
                        tmax = metarr->tmax[py*365+pday];
                        tmin_annavg += metarr->tmin[py*365+pday];
                    }

                    if (tmax > new_tmax) new_tmax = tmax;

                } /* end pday loop */

                tmax_ann += new_tmax;
            } /* end py loop */
            tmax_ann /= (double) phenyears;
            /* 92% of tmax_ann is the threshold used in grass offset below */
            tmax_ann *= 0.92;
            tmin_annavg /= (double) phenyears * 365.0;

            /* loop through phenyears again, fill onset and offset arrays */
            for (py=0 ; py<phenyears ; py++)
            {
                sum_soilt = 0.0;
                sum_prcp = 0.0;
                onset_day = offset_day = -1;
                for (pday=0 ; pday<365 ; pday++)
                {
                    if (south)
                    {
                        if (py==0 && pday<182)
                        {
                            /* use the end of the first year to fill the 
                               beginning of a southern hemisphere phenological
                               year */
                            phensoilt = metarr->tavg_ra[183+pday];
                            phenprcp = metarr->prcp[183+pday];
                            grass_prcpyear[pday] = phenprcp;
                            grass_tminyear[pday] = metarr->tmin[183+pday];
                            grass_tmaxyear[pday] = metarr->tmax[183+pday];
                        }
                        else if (py==phenyears-1 && pday>181)
                        {
                            /* use the beginning of the last year to fill the
                               end of the last phenological year */
                            phensoilt = metarr->tavg_ra[ndays-547+pday];
                            phenprcp = metarr->prcp[ndays-547+pday];
                            grass_prcpyear[pday] = phenprcp;
                            grass_tminyear[pday] = metarr->tmin[ndays-547+pday];
                            grass_tmaxyear[pday] = metarr->tmax[ndays-547+pday];
                        }
                        else
                        {
                            phensoilt = metarr->tavg_ra[py*365-182+pday];
                            phenprcp = metarr->prcp[py*365-182+pday];
                            grass_prcpyear[pday] = phenprcp;
                            grass_tminyear[pday] = metarr->tmin[py*365-182+pday];
                            grass_tmaxyear[pday] = metarr->tmax[py*365-182+pday];
                        }
                    }
                    else /* north */
                    {
                        phensoilt = metarr->tavg_ra[py*365+pday];
                        phenprcp = metarr->prcp[py*365+pday];
                        grass_prcpyear[pday] = phenprcp;
                        grass_tminyear[pday] = metarr->tmin[py*365+pday];
                        grass_tmaxyear[pday] = metarr->tmax[py*365+pday];
                    }

                    /* grass onset test */
                    if (onset_day == -1)
                    {
                        if (phensoilt > 0.0) sum_soilt += phensoilt;
                        sum_prcp += phenprcp;
                        if (sum_soilt >= grass_stsumcrit &&
                                sum_prcp >= grass_prcpcrit) onset_day = pday;
                    }

                } /* end pday loop */

                /* do averaging operations on grass_prcpyear and grass_tminyear,
                   and do tests for offset day. Offset due to hot & dry can't
                   happen within one month after the onset day, and offset due
                   to cold can't happen before midyear (yearday 182) */
                if (onset_day != -1)
                {
                    /* calculate three-day boxcar average of tmin */
                    if (boxcar_smooth(grass_tminyear, grass_3daytmin, 365,3,0))
                    {
                        bgc_printf(BV_ERROR, "Error in prephenology() call to boxcar()\n");
                        ok=0;
                    }

                    for (pday=onset_day+30 ; pday<365 ; pday++)
                    {
                        /* calculate the previous 31-day prcp total */
                        psum_startday = pday - 30;
                        grass_prcpprev = 0.0;
                        for (i=psum_startday ; i<=pday ; i++)
                        {
                            grass_prcpprev += grass_prcpyear[i];
                        }

                        /* calculate the next 7-day prcp total */
                        if (pday > 358) psum_stopday = 364;
                        else psum_stopday = pday + 6;
                        grass_prcpnext = 0.0;
                        for (i=pday ; i<=psum_stopday ; i++)
                        {
                            grass_prcpnext += grass_prcpyear[i];
                        }

                        /* test for hot and dry conditions */
                        if (offset_day == -1)
                        {
                            if (grass_prcpprev < grass_prcpprevcrit && 
                                    grass_prcpnext < grass_prcpnextcrit &&
                                    grass_tmaxyear[pday] > tmax_ann)
                                offset_day = pday;
                        }

                        /* test for cold offset condition */
                        if (offset_day == -1)
                        {
                            if (pday > 182 &&
                                    grass_3daytmin[pday] <= tmin_annavg)
                                offset_day = pday;
                        }

                    } /* end of pdays loop for grass offset testing */
                } /* end of if onset_day != -1 block */

                /* now do some exception handling for this year's phenology */
                if (onset_day != -1)
                {
                    /* leaves are turned on sometime this year */
                    /* subtract 15 days from onset day to approximate the
                       start of the new growth period, instead of the middle of
                       the new growth period, as is used in the White et al. ms. */
                    if (onset_day >= 15)
                    {
                        onset_day -= 15;
                    }
                    else onset_day = 0;

                    /* if leaves never got turned off, force off on last day */
                    if (offset_day == -1) offset_day = 364;

                    /* force onset and offset to be at least one day apart */
                    if (onset_day == offset_day)
                    {
                        if (onset_day > 0) onset_day--;
                        else offset_day++;
                    }
                }
                else
                {
                    /* leaves never got turned on, this is a non-growth
                       year.  This probably indicates a problem with the
                       phenology model */
                    onset_day = -1;
                    offset_day = -1;
                }

                /* save these onset and offset days and go to the next
                   phenological year */
                onday_arr[py] = onset_day;
                offday_arr[py] = offset_day;

            } /* end phenyears loop for filling onset and offset arrays */
        } /* end else !woody (grass phenology model) */

        /* now the onset and offset days are established for each phenyear,
           either by the deciduous tree or the grass model.  Next loop through
           phenyears filling the phenology signal arrays and copying them to 
           the permanent phen struct arrays */
        for (py=0 ; py<phenyears ; py++)
        {
            /* zero the 365-day phen arrays */
            for (pday=0 ; pday<365 ; pday++)
            {
                remdays_curgrowth[pday] = 0;
                remdays_transfer[pday] = 0;
                predays_transfer[pday] = 0;
                remdays_litfall[pday] = 0;
                predays_litfall[pday] = 0;
            }

            onday = onday_arr[py];
            offday = offday_arr[py];

            if (onday == -1 && offday == -1)
            {
                /* this is the special signal to repress all vegetation
                   growth */
                for (pday=0 ; pday<365 ; pday++)
                {
                    remdays_curgrowth[pday] = 0;
                    remdays_transfer[pday] = 0;
                    predays_transfer[pday] = 0;
                    remdays_litfall[pday] = 0;
                    predays_litfall[pday] = 0;
                }
            } /* end if special no-growth signal */
            else
            {
                /* normal growth year */
                ngrowthdays = offday - onday;
                if (ngrowthdays < 1)
                {
                    bgc_printf(BV_ERROR, "FATAL ERROR: ngrowthdays < 1\n");
                    bgc_printf(BV_ERROR, "ngrowthdays = %d\n",ngrowthdays);
                    bgc_printf(BV_ERROR, "onday = %d\toffday = %d\tphenyear = %d\n",
                            onday,offday,py);
                    ok=0;
                }
                /* define the length of the transfer and litfall periods */
                /* calculate number of transfer days and number of litfall days
                   as proportions of number of growth days, as specified by user.
                   Round and truncate to force these values between 1 and 
                   ngrowthdays */
                t1 = epc->transfer_pdays * (double)ngrowthdays;
                sprintf(round,"%.0f",t1);
                ntransferdays = atoi(round);
                if (ntransferdays < 1) ntransferdays = 1;
                if (ntransferdays > ngrowthdays) ntransferdays = ngrowthdays;
                t1 = epc->litfall_pdays * (double)ngrowthdays;
                sprintf(round,"%.0f",t1);
                nlitfalldays = atoi(round);
                if (nlitfalldays < 1) nlitfalldays = 1;
                if (nlitfalldays > ngrowthdays) nlitfalldays = ngrowthdays;

                for (pday=0 ; pday<onday ; pday++)
                {
                    remdays_curgrowth[pday] = 0;
                    remdays_transfer[pday] = 0;
                    remdays_litfall[pday] = 0;
                    predays_transfer[pday] = 0;
                    predays_litfall[pday] = 0;
                }
                counter = ngrowthdays;
                for (pday=onday ; pday<offday ; pday++)
                {
                    remdays_curgrowth[pday] = counter;
                    counter--;
                }
                for (pday=offday ; pday<365 ; pday++)
                {
                    remdays_curgrowth[pday] = 0;
                }
                counter = ntransferdays;
                for (pday=onday ; pday<onday+ntransferdays ; pday++)
                {
                    remdays_transfer[pday] = counter;
                    predays_transfer[pday] = ntransferdays - counter;
                    counter--;
                }
                for (pday=onday+ntransferdays ; pday<=offday ; pday++)
                {
                    remdays_transfer[pday] = 0;
                    predays_transfer[pday] = ntransferdays;
                }
                for (pday=offday+1 ; pday<365 ; pday++)
                {
                    remdays_transfer[pday] = 0;
                    predays_transfer[pday] = 0;
                }
                for (pday=onday ; pday<offday-nlitfalldays+1 ; pday++)
                {
                    remdays_litfall[pday] = 0;
                    predays_litfall[pday] = 0;
                }
                counter = nlitfalldays;
                for (pday=offday-nlitfalldays+1 ; pday<=offday ; pday++)
                {
                    remdays_litfall[pday] = counter;
                    predays_litfall[pday] = nlitfalldays - counter;
                    counter--;
                }
                for (pday=offday+1 ; pday<365 ; pday++)
                {
                    remdays_litfall[pday] = 0;
                    predays_litfall[pday] = 0;
                }
            } /* end else normal growth year */

            /* now put the signals for this phenological year into the
               right place in the permanent phen struct arrays */ 
            if (south)
            {
                if (py==0)
                {
                    /* only copy the second half of this phenological
                       year to the permanent phenology array */
                    for (pday=182 ; pday<365 ; pday++)
                    {
                        phen->remdays_curgrowth[pday-182] = remdays_curgrowth[pday];
                        phen->remdays_transfer[pday-182] = remdays_transfer[pday];
                        phen->remdays_litfall[pday-182] = remdays_litfall[pday];
                        phen->predays_transfer[pday-182] = predays_transfer[pday];
                        phen->predays_litfall[pday-182] = predays_litfall[pday];
                    }
                }
                else if (py==phenyears-1)
                {
                    /* only copy the first half of this phenological
                       year to the permanent phenology array */
                    for (pday=0 ; pday<182 ; pday++)
                    {
                        phen->remdays_curgrowth[py*365-182+pday] = remdays_curgrowth[pday];
                        phen->remdays_transfer[py*365-182+pday] = remdays_transfer[pday];
                        phen->remdays_litfall[py*365-182+pday] = remdays_litfall[pday];
                        phen->predays_transfer[py*365-182+pday] = predays_transfer[pday];
                        phen->predays_litfall[py*365-182+pday] = predays_litfall[pday];
                    }
                }
                else
                {
                    for (pday=0 ; pday<365 ; pday++)
                    {
                        phen->remdays_curgrowth[py*365-182+pday] = remdays_curgrowth[pday];
                        phen->remdays_transfer[py*365-182+pday] = remdays_transfer[pday];
                        phen->remdays_litfall[py*365-182+pday] = remdays_litfall[pday];
                        phen->predays_transfer[py*365-182+pday] = predays_transfer[pday];
                        phen->predays_litfall[py*365-182+pday] = predays_litfall[pday];
                    }
                }
            } /* end if south */
            else
            {
                /* north */
                for (pday=0 ; pday<365 ; pday++)
                {
                    phen->remdays_curgrowth[py*365+pday] = remdays_curgrowth[pday];
                    phen->remdays_transfer[py*365+pday] = remdays_transfer[pday];
                    phen->remdays_litfall[py*365+pday] = remdays_litfall[pday];
                    phen->predays_transfer[py*365+pday] = predays_transfer[pday];
                    phen->predays_litfall[py*365+pday] = predays_litfall[pday];
                }
            } /* end if north */
        } /* end phenyears loop for filling permanent arrays */
    } /* end else phenology model block */

    /* free the local array memory */
    free(onday_arr);
    free(offday_arr);

    return (!ok);
}

int free_phenmem(phenarray_struct* phen)
{
	int ok=1;
	
	/* free memory in phenology arrays */
	free(phen->remdays_curgrowth);
	free(phen->remdays_transfer);
	free(phen->remdays_litfall);
	free(phen->predays_transfer);
	free(phen->predays_litfall);
	
	return (!ok);
}
int phenology(const epconst_struct* epc, const phenology_struct* phen, epvar_struct* epv, cstate_struct* cs, cflux_struct* cf, nstate_struct* ns, nflux_struct* nf)
{
	int ok=1;
	double ndays;
	double leaflitfallc, frootlitfallc;
	double livestemtovrc, livestemtovrn;
	double livecroottovrc, livecroottovrn;
	double drate;
	
	/* phenological control for EVERGREENS */
	if (epc->evergreen)
	{
		/* transfer growth fluxes */
		/* check for days left in transfer growth period */
		/* AAN - yes, this is an assignment */
		if ((ndays = phen->remdays_transfer))
		{
			/* calculate rates required to empty each transfer
			compartment by the end of transfer period, at approximately a
			constant rate of transfer */
			cf->leafc_transfer_to_leafc = cs->leafc_transfer / ndays;
			nf->leafn_transfer_to_leafn = ns->leafn_transfer / ndays;
			cf->frootc_transfer_to_frootc = cs->frootc_transfer / ndays;
			nf->frootn_transfer_to_frootn = ns->frootn_transfer / ndays;
			if (epc->woody)
			{
				cf->livestemc_transfer_to_livestemc = cs->livestemc_transfer / ndays;
				nf->livestemn_transfer_to_livestemn = ns->livestemn_transfer / ndays;
				cf->deadstemc_transfer_to_deadstemc = cs->deadstemc_transfer / ndays;
				nf->deadstemn_transfer_to_deadstemn = ns->deadstemn_transfer / ndays;
				cf->livecrootc_transfer_to_livecrootc = cs->livecrootc_transfer / ndays;
				nf->livecrootn_transfer_to_livecrootn = ns->livecrootn_transfer / ndays;
				cf->deadcrootc_transfer_to_deadcrootc = cs->deadcrootc_transfer / ndays;
				nf->deadcrootn_transfer_to_deadcrootn = ns->deadcrootn_transfer / ndays;
			}
		}

		/* litterfall happens everyday, at a rate determined each year
		on the annual allocation day.  To prevent litterfall from driving
		pools negative in the case of a very high mortality, fluxes are
		checked and set to zero when the pools get too small. */

		/* leaf litterfall */
		leaflitfallc = epv->day_leafc_litfall_increment;
		if (leaflitfallc > cs->leafc) leaflitfallc = cs->leafc;
		if (ok && leaf_litfall(epc,leaflitfallc,cf,nf))
		{
			bgc_printf(BV_ERROR, "Error in call to leaf_litfall() from phenology()\n");
			ok=0;
		}

		/* fine root litterfall */
		frootlitfallc = epv->day_frootc_litfall_increment;
		if (frootlitfallc > cs->frootc) frootlitfallc = cs->frootc;
		if (ok && froot_litfall(epc,frootlitfallc,cf,nf))
		{
			bgc_printf(BV_ERROR, "Error in call to froot_litfall() from phenology()\n");
			ok=0;
		}
		
		/* turnover of live wood to dead wood also happens every day, at a
		rate determined once each year, using the annual maximum livewoody
		compartment masses and the specified livewood turnover rate */
		if (epc->woody)
		{
			/* turnover from live stem wood to dead stem wood */
			livestemtovrc = epv->day_livestemc_turnover_increment;
			livestemtovrn = livestemtovrc / epc->livewood_cn;
			if (livestemtovrc > cs->livestemc) livestemtovrc = cs->livestemc;
			if (livestemtovrn > ns->livestemn) livestemtovrn = ns->livestemn;
			if (livestemtovrc && livestemtovrn)
			{
				cf->livestemc_to_deadstemc = livestemtovrc;
				nf->livestemn_to_deadstemn = livestemtovrc / epc->deadwood_cn;
				nf->livestemn_to_retransn = livestemtovrn - nf->livestemn_to_deadstemn;
			}

			/* turnover from live coarse root wood to dead coarse root wood */
			livecroottovrc = epv->day_livecrootc_turnover_increment;
			livecroottovrn = livecroottovrc / epc->livewood_cn;
			if (livecroottovrc > cs->livecrootc) livecroottovrc = cs->livecrootc;
			if (livecroottovrn > ns->livecrootn) livecroottovrn = ns->livecrootn;
			if (livecroottovrc && livecroottovrn)
			{
				cf->livecrootc_to_deadcrootc = livecroottovrc;
				nf->livecrootn_to_deadcrootn = livecroottovrc / epc->deadwood_cn;
				nf->livecrootn_to_retransn = livecroottovrn - nf->livecrootn_to_deadcrootn;
			}
		}
		
	} /* end if evergreen */
	else
	{
		/* deciduous */
		/* transfer growth fluxes */
		/* check for days left in transfer growth period */
		/* AAN - yes, this is an assignment */
		if ((ndays = phen->remdays_transfer))
		{
			/* transfer rate is defined to be a linearly decreasing
			function that reaches zero on the last day of the transfer
			period */
			cf->leafc_transfer_to_leafc = 2.0*cs->leafc_transfer / ndays;
			nf->leafn_transfer_to_leafn = 2.0*ns->leafn_transfer / ndays;
			cf->frootc_transfer_to_frootc = 2.0*cs->frootc_transfer / ndays;
			nf->frootn_transfer_to_frootn = 2.0*ns->frootn_transfer / ndays;
			if (epc->woody)
			{
				cf->livestemc_transfer_to_livestemc = 2.0*cs->livestemc_transfer / ndays;
				nf->livestemn_transfer_to_livestemn = 2.0*ns->livestemn_transfer / ndays;
				cf->deadstemc_transfer_to_deadstemc = 2.0*cs->deadstemc_transfer / ndays;
				nf->deadstemn_transfer_to_deadstemn = 2.0*ns->deadstemn_transfer / ndays;
				cf->livecrootc_transfer_to_livecrootc = 2.0*cs->livecrootc_transfer / ndays;
				nf->livecrootn_transfer_to_livecrootn = 2.0*ns->livecrootn_transfer / ndays;
				cf->deadcrootc_transfer_to_deadcrootc = 2.0*cs->deadcrootc_transfer / ndays;
				nf->deadcrootn_transfer_to_deadcrootn = 2.0*ns->deadcrootn_transfer / ndays;
			}
		}
		
		/* litterfall */
		/* defined such that all live material is removed by the end of the
		litterfall period, with a linearly ramping removal rate. assumes that
		the initial rate on the first day of litterfall is 0.0. */
		/* AAN - yes, this is an assignment */
		if ((ndays = phen->remdays_litfall))
		{
			if (ndays == 1.0)
			{
				/* last day of litterfall, special case to gaurantee
				that pools go to 0.0 */
				leaflitfallc = cs->leafc;
				frootlitfallc = cs->frootc;
			}
			else
			{
				/* otherwise, assess litterfall 
				rates as described above */
				leaflitfallc = epv->day_leafc_litfall_increment;
				drate = 2.0*(cs->leafc - leaflitfallc*ndays)/(ndays*ndays);
				epv->day_leafc_litfall_increment += drate;
				frootlitfallc = epv->day_frootc_litfall_increment;
				drate = 2.0*(cs->frootc - frootlitfallc*ndays)/(ndays*ndays);
				epv->day_frootc_litfall_increment += drate;
			}
			/* leaf litterfall */
			if (leaflitfallc > cs->leafc) leaflitfallc = cs->leafc;
			if (ok && leaflitfallc && leaf_litfall(epc,leaflitfallc,cf,nf))
			{
				bgc_printf(BV_ERROR, "Error in call to leaf_litfall() from phenology()\n");
				ok=0;
			}
			/* fine root litterfall */
			if (frootlitfallc > cs->frootc) frootlitfallc = cs->frootc;
			if (ok && frootlitfallc && froot_litfall(epc,frootlitfallc,cf,nf))
			{
				bgc_printf(BV_ERROR, "Error in call to froot_litfall() from phenology()\n");
				ok=0;
			}
		} /* end if deciduous litterfall day */
		
		/* turnover of livewood to deadwood happens each day, just as for
		evergreen types, at a rate determined from the annual maximum
		livewood mass and the specified turnover rate */
		if (epc->woody)
		{
			/* turnover from live stem wood to dead stem wood */
			livestemtovrc = epv->day_livestemc_turnover_increment;
			livestemtovrn = livestemtovrc / epc->livewood_cn;
			if (livestemtovrc > cs->livestemc) livestemtovrc = cs->livestemc;
			if (livestemtovrn > ns->livestemn) livestemtovrn = ns->livestemn;
			if (livestemtovrc && livestemtovrn)
			{
				cf->livestemc_to_deadstemc = livestemtovrc;
				nf->livestemn_to_deadstemn = livestemtovrc / epc->deadwood_cn;
				nf->livestemn_to_retransn = livestemtovrn - nf->livestemn_to_deadstemn;
			}

			/* turnover from live coarse root wood to dead coarse root wood */
			livecroottovrc = epv->day_livecrootc_turnover_increment;
			livecroottovrn = livecroottovrc / epc->livewood_cn;
			if (livecroottovrc > cs->livecrootc) livecroottovrc = cs->livecrootc;
			if (livecroottovrn > ns->livecrootn) livecroottovrn = ns->livecrootn;
			if (livecroottovrc && livecroottovrn)
			{
				cf->livecrootc_to_deadcrootc = livecroottovrc;
				nf->livecrootn_to_deadcrootn = livecroottovrc / epc->deadwood_cn;
				nf->livecrootn_to_retransn = livecroottovrn - nf->livecrootn_to_deadcrootn;
			}
		}
		
	} /* end if deciduous */
	
	/* for woody types, find annual maximum value for live stemc and live crootc
	calculation of livewood turnover rates */
	if (epc->woody)
	{
		if (epv->annmax_livestemc < cs->livestemc) epv->annmax_livestemc = cs->livestemc;
		if (epv->annmax_livecrootc < cs->livecrootc) epv->annmax_livecrootc = cs->livecrootc;
	}	
	
	/* for all types, find annual maximum leafc */
	if (epv->annmax_leafc < cs->leafc) epv->annmax_leafc = cs->leafc;
	if (epv->annmax_frootc < cs->frootc) epv->annmax_frootc = cs->frootc;
	
	return (!ok);
}

int leaf_litfall(const epconst_struct* epc, double litfallc,
cflux_struct* cf, nflux_struct* nf)
{
	int ok=1;
	double c1,c2,c3,c4;
	double n1,n2,n3,n4;
	double nretrans;
	double avg_cn;
	double litfalln;
	
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
	nretrans = (litfallc/avg_cn) - (litfalln);
	
	if (ok)
	{
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
	
	return (!ok);
}

int froot_litfall(const epconst_struct* epc, double litfallc, 
cflux_struct* cf, nflux_struct* nf)
{
	int ok=1;
	double c1,c2,c3,c4;
	double n1,n2,n3,n4;
	double avg_cn;
	
	avg_cn = epc->froot_cn;
	
	c1 = litfallc * epc->frootlitr_flab;
	n1 = c1 / avg_cn;
	c2 = litfallc * epc->frootlitr_fucel;
	n2 = c2 / avg_cn;
	c3 = litfallc * epc->frootlitr_fscel;
	n3 = c3 / avg_cn;
	c4 = litfallc * epc->frootlitr_flig;
	n4 = c4 / avg_cn;
	
	if (ok)
	{
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
	
	return (!ok);
}

