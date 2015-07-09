#include "bgc.h"

void metarr_init (bgc_struct bgc, pihm_struct pihm, lsm_struct noah, int start_time, int end_time)
{
    /*
     * Generate meteorological forcing array for spin-up
     */
    double          sfctmp;
    double          solar;
    double          dayl, prev_dayl;
    double          rh;
    metarr_struct  *metarr;
    double       ***met_forcing;
    double       ***radn_forcing;
    int             hour;
    double          t;
    int             i, j, k;
    int             length;
    int             daylight_counter;
    spa_data        spa;
    int             spa_result;
    time_t          rawtime;
    struct tm      *timestamp;
    double         *swc;
    double         *stc;
    double         *soilw;
    double         *fluxsub[3];
    double          es;
    double          pres;
    int             ind;

    printf ("Initialize meteorological forcing array for model spin-up ...\n");

    length = (int)((end_time - start_time) / 24. / 3600.);

    swc = (double *) malloc (pihm->numele * sizeof (double));
    stc = (double *) malloc (pihm->numele * sizeof (double));
    soilw = (double *) malloc (pihm->numele * sizeof (double));

    met_forcing = (double ***) malloc (pihm->forcing.nts[METEO_TS] * sizeof (double **));
    radn_forcing = (double ***) malloc (pihm->forcing.nts[METEO_TS] * sizeof (double **));

    for (j = 0; j < pihm->forcing.nts[METEO_TS]; j++)
    {
        met_forcing[j] = (double **) malloc (24 * sizeof (double *));
        radn_forcing[j] = (double **) malloc (24 * sizeof (double *));

        for (hour = 0; hour < 24; hour++)
        {
            met_forcing[j][hour] = (double *) malloc (NUM_METEO_TS * sizeof (double));
            radn_forcing[j][hour] = (double *) malloc (4 * sizeof (double));
        }
    }

    for (k = 0; k < 3; k++)
    {
        fluxsub[k] = (double *) malloc (pihm->numele * sizeof (double));
    }

    for (i = 0; i < pihm->numele; i++)
    {
        bgc->grid[i].metarr.tmax = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.tmin = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.prcp = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.vpd = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.q2d = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.swavgfd = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.par = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.dayl = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.prev_dayl = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.tavg = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.tday = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.tnight = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.tsoil = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.swc = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.pa = (double *) malloc (length * sizeof (double));
        bgc->grid[i].metarr.soilw = (double *) malloc (length * sizeof (double));
        for (j = 0; j < 3; j++)
        {
            bgc->grid[i].metarr.subflux[j] = (double *) malloc (length * sizeof (double));
        }
    }

    for (j = 0; j < length; j++)
    {
        t = start_time + j * 24 * 3600;
        rawtime = t;
        timestamp = gmtime (&rawtime);
        spa.year = timestamp->tm_year + 1900;
        spa.month = timestamp->tm_mon + 1;
        spa.day = timestamp->tm_mday;
        spa.hour = timestamp->tm_hour;
        spa.minute = timestamp->tm_min;
        spa.second = timestamp->tm_sec;

        spa.timezone = 0;
        spa.delta_t = 67;
        spa.delta_ut1 = 0;
        spa.atmos_refract = 0.5667;

        spa.longitude = bgc->grid[0].sitec.lon;
        spa.latitude = bgc->grid[0].sitec.lat;
        spa.elevation = 0.0;
        for (i = 0; i < pihm->numele; i++)
        {
            spa.elevation += pihm->elem[i].topo.zmax;
        }
        spa.elevation /= (double) pihm->numele;
        /*
         * Calculate surface pressure based on FAO 1998 method (Narasimhan 2002) 
         */
        spa.pressure = 1013.25 * pow ((293.0 - 0.0065 * spa.elevation) / 293.0, 5.26);
        spa.temperature = noah->genprmt.tbot_data;

        spa.function = SPA_ZA_RTS;
        spa_result = spa_calculate (&spa);

        /* daylength (s) */
        dayl = (spa.sunset - spa.sunrise) * 3600.0;
        dayl = (dayl < 0.0) ? (dayl + 24.0 * 3600.0) : dayl;

        if (j == 0)
        {
            rawtime = rawtime - 24 * 3600;
            timestamp = gmtime (&rawtime);
            spa.year = timestamp->tm_year + 1900;
            spa.month = timestamp->tm_mon + 1;
            spa.day = timestamp->tm_mday;
            spa.hour = timestamp->tm_hour;
            spa.minute = timestamp->tm_min;
            spa.second = timestamp->tm_sec;
            spa_result = spa_calculate (&spa);
            prev_dayl = (spa.sunset - spa.sunrise) * 3600.0;
            prev_dayl = (prev_dayl < 0.0) ? (prev_dayl + 12.0 * 3600.0) : prev_dayl;
        }

        for (hour = 0; hour < 24; hour++)
        {
            for (k = 0; k < pihm->forcing.nts[METEO_TS]; k++)
            {
                IntrplForcing (pihm->forcing.ts[METEO_TS][k], t + hour * 3600, NUM_METEO_TS, &met_forcing[k][hour][0]);
            }

            if (noah->rad_mode > 0)
            {
                rawtime = t + hour * 3600;
                timestamp = gmtime (&rawtime);
                spa.year = timestamp->tm_year + 1900;
                spa.month = timestamp->tm_mon + 1;
                spa.day = timestamp->tm_mday;
                spa.hour = timestamp->tm_hour;
                spa.minute = timestamp->tm_min;
                spa.second = timestamp->tm_sec;
                spa.function = SPA_ZA;
                spa_result = spa_calculate (&spa);

                if (spa_result != 0)
                {
                    printf ("SPA Error Code: %d\n", spa_result);
                    exit (0);
                }
                spa.azimuth180 = mod ((360.0 + spa.azimuth180), 360.0);

                for (k = 0; k < pihm->forcing.nts[METEO_TS]; k++)
                {
                    IntrplForcing (noah->forcing.ts[k], t + hour * 3600, 2, &radn_forcing[k][hour][0]);
                    radn_forcing[k][hour][3] = spa.azimuth180;
                    radn_forcing[k][hour][2] = spa.zenith;
                }
            }
        }

        IntrplForcing (bgc->forcing.ts[SWC_TS][0], t + 24 * 3600, pihm->numele, &swc[0]);
        IntrplForcing (bgc->forcing.ts[STC_TS][0], t + 24 * 3600, pihm->numele, &stc[0]);
        IntrplForcing (bgc->forcing.ts[SOILM_TS][0], t + 24 * 3600, pihm->numele, &soilw[0]);
        for (k = 0; k < 3; k++)
        {
            IntrplForcing (bgc->forcing.ts[SUBFLX_TS][k], t + 24 * 3600, pihm->numele, &fluxsub[k][0]);
        }

        for (i = 0; i < pihm->numele; i++)
        {
            ind = pihm->attrib_tbl.meteo[i] - 1;

            metarr = &(bgc->grid[i].metarr);
            metarr->dayl[j] = dayl;
            if (j == 0)
            {
                metarr->prev_dayl[j] = prev_dayl;
            }
            else
            {
                metarr->prev_dayl[j] = metarr->prev_dayl[j - 1];
            }

            metarr->prcp[j] = 0.0;
            metarr->tmax[j] = -999.0;
            metarr->tmin[j] = 999.0;
            metarr->tavg[j] = 0.0;
            metarr->tday[j] = 0.0;
            metarr->tnight[j] = 0.0;
            metarr->vpd[j] = 0.0;
            metarr->q2d[j] = 0.0;
            metarr->swavgfd[j] = 0.0;
            metarr->par[j] = 0.0;
            metarr->pa[j] = 0.0;
            rh = 0.0;
            daylight_counter = 0;

            for (hour = 0; hour < 24; hour++)
            {
                metarr->prcp[j] += met_forcing[ind][hour][PRCP_TS] * 3600.0;    /* Convert from kg m-2 s-1 to kg m-2 */
                sfctmp = met_forcing[ind][hour][SFCTMP_TS] - 273.15;
                metarr->tmax[j] = (sfctmp > metarr->tmax[j]) ? sfctmp : metarr->tmax[j];
                metarr->tmin[j] = (sfctmp < metarr->tmin[j]) ? sfctmp : metarr->tmin[j];
                metarr->tavg[j] += sfctmp;
                if (noah->rad_mode > 0)
                {
                    solar = TopoRadiation (radn_forcing[ind][hour][0], radn_forcing[ind][hour][1], radn_forcing[ind][hour][2], radn_forcing[ind][hour][3], noah->grid[i].slope, noah->grid[i].aspect, noah->grid[i].h_phi, noah->grid[i].svf);
                }
                else
                {
                    solar = met_forcing[ind][hour][SOLAR_TS];
                }

                if (solar > 1.0)
                {
                    metarr->tday[j] += sfctmp;
                    rh = met_forcing[ind][hour][RH_TS] / 100.;
                    metarr->vpd[j] += (1.0 - rh) * 611.2 * exp (17.67 * sfctmp / (sfctmp + 243.5));
                    es = 611.2 * exp (17.67 * sfctmp / (sfctmp + 243.5));
                    pres = met_forcing[ind][hour][PRES_TS];
                    metarr->q2d[j] += (0.622 * es) / (pres - (1.0 - 0.622) * es) - (0.622 * es * rh) / (pres - (1.0 - 0.622) * es * rh);
                    metarr->pa[j] += met_forcing[ind][hour][PRES_TS];
                    metarr->swavgfd[j] += solar;
                    daylight_counter++;
                }
                else
                {
                    metarr->tnight[j] += sfctmp;
                }
            }

            metarr->tavg[j] /= 24.0;
            metarr->tday[j] /= (double) daylight_counter;
            metarr->tnight[j] /= (24.0 - (double) daylight_counter);
            metarr->vpd[j] /= (double) daylight_counter;
            metarr->q2d[j] /= (double) daylight_counter;
            metarr->pa[j] /= (double) daylight_counter;
            metarr->swavgfd[j] /= (double) daylight_counter;

            metarr->par[j] = metarr->swavgfd[j] * RAD2PAR;

            metarr->tsoil[j] = stc[i] - 273.15;
            metarr->swc[j] = swc[i];
            metarr->soilw[j] = soilw[i] * 1000.0;
            for (k = 0; k < 3; k++)
            {
                /* Convert from m3/s to kg/m2/d */
                //if (pihm->elem[i].forc.bc_type[k] < 0.0 && fluxsub[k][i] < 0.0)
                //{
                //    metarr->subflux[k][j] = 0.0;
                //}
                //else
                //{
                    metarr->subflux[k][j] = 1000.0 * fluxsub[k][i] * 24.0 * 3600.0 / pihm->elem[i].topo.area;
                //}
            }
        }
    }

    free (swc);
    free (stc);
    free (soilw);

    for (j = 0; j < pihm->forcing.nts[METEO_TS]; j++)
    {
        for (hour = 0; hour < 24; hour++)
        {
            free (met_forcing[j][hour]);
            free (radn_forcing[j][hour]);
        }
        free (met_forcing[j]);
        free (radn_forcing[j]);
    }
    free (met_forcing);
    free (radn_forcing);
}
