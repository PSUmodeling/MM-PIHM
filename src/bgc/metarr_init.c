#include "bgc.h"
#include "../noah/noah.h"

void metarr_init (bgc_struct BGCM, Model_Data PIHM, LSM_STRUCT LSM, double start_time, double end_time)
{
    /*
     * Generate meteorological forcing array for spin-up
     */

    double          sfctmp;
    double          solar;
    double          dayl, prev_dayl;
    double          RH;
    metarr_struct  *metarr;
    double          PIHM_forcing[PIHM->NumTS][24][7];
    double          rad_forcing[PIHM->NumTS][24][4];
    int             hour;
    double          t;
    int             i, j, k;
    int             length;
    int             daylight_coutner;
    spa_data        spa;
    int             spa_result;
    time_t          rawtime;
    struct tm      *timestamp;
    double          swc[PIHM->NumEle];
    double          stc[PIHM->NumEle];
    double          soilw[PIHM->NumEle];
    double          subflux[3][PIHM->NumEle];

    printf ("Initialize meteorological forcing array for model spin-up ...\n");

    length = (int)((end_time - start_time) / 24. / 3600.);

    //printf ("length = %d\n", length);

    for (i = 0; i < PIHM->NumEle; i++)
    {
        BGCM->grid[i].metarr.tmax = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.tmin = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.prcp = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.vpd = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.swavgfd = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.par = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.dayl = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.prev_dayl = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.tavg = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.tday = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.tnight = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.tsoil = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.swc = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.pa = (double *)malloc (length * sizeof (double));
        BGCM->grid[i].metarr.soilw = (double *)malloc (length * sizeof (double));
        for (j = 0; j < 3; j++)
        {
            BGCM->grid[i].metarr.subflux[j] = (double *)malloc (length * sizeof (double));
        }
    }

    for (j = 0; j < length; j++)
    {
        //printf ("DAY %d\n", j);
        t = start_time + j * 24. * 3600.;
        rawtime = (int)t;
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

        spa.longitude = BGCM->grid[0].sitec.lon;
        spa.latitude = BGCM->grid[0].sitec.lat;
        spa.elevation = 0.;
        for (i = 0; i < PIHM->NumEle; i++)
            spa.elevation = spa.elevation + (double)PIHM->Ele[i].zmax;
        spa.elevation = spa.elevation / (double)PIHM->NumEle;
        /*
         * Calculate surface pressure based on FAO 1998 method (Narasimhan 2002) 
         */
        spa.pressure = 1013.25 * pow ((293. - 0.0065 * spa.elevation) / 293., 5.26);
        spa.temperature = LSM->GENPRMT.TBOT_DATA;

        spa.function = SPA_ZA_RTS;
        spa_result = spa_calculate (&spa);

        /* daylength (s) */
        dayl = (spa.sunset - spa.sunrise) * 3600.;
        dayl = dayl < 0. ? (dayl + 24. * 3600.) : dayl;

        if (j == 0)
        {
            rawtime = rawtime - 24. * 3600.;
            timestamp = gmtime (&rawtime);
            spa.year = timestamp->tm_year + 1900;
            spa.month = timestamp->tm_mon + 1;
            spa.day = timestamp->tm_mday;
            spa.hour = timestamp->tm_hour;
            spa.minute = timestamp->tm_min;
            spa.second = timestamp->tm_sec;
            spa_result = spa_calculate (&spa);
            prev_dayl = (spa.sunset - spa.sunrise) * 3600.;
            prev_dayl = prev_dayl < 0. ? (prev_dayl + 12. * 3600.) : prev_dayl;
        }

        for (hour = 0; hour < 24; hour++)
        {
            for (k = 0; k < PIHM->NumTS; k++)
                MultiInterpolation (&PIHM->TSD_meteo[k], t + hour * 3600., &PIHM_forcing[k][hour][0], 7);
            if (LSM->RAD_MODE > 0)
            {
                rawtime = (int) (t + hour * 3600.);
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
                spa.azimuth180 = mod ((360. + spa.azimuth180), 360.);

                for (k = 0; k < PIHM->NumTS; k++)
                {
                    MultiInterpolation (&LSM->TSD_rad[k], t + hour * 3600., &rad_forcing[k][hour][0], 2);
                    rad_forcing[k][hour][3] = spa.azimuth180;
                    //printf ("Azimuth %lf\n", rad_forcing[k][hour][3]);
                    rad_forcing[k][hour][2] = spa.zenith;
                }
            }
        }

        MultiInterpolation (&BGCM->Forcing[SWC_TS][0], t + 24. * 3600., &swc[0], PIHM->NumEle);
        MultiInterpolation (&BGCM->Forcing[STC_TS][0], t + 24. * 3600., &stc[0], PIHM->NumEle);
        MultiInterpolation (&BGCM->Forcing[SOILM_TS][0], t + 24. * 3600., &soilw[0], PIHM->NumEle);
        for (k = 0; k < 3; k++)
        {
            MultiInterpolation (&BGCM->Forcing[SUBFLX_TS][k], t + 24. * 3600., &subflux[k][0], PIHM->NumEle);
        }

        for (i = 0; i < PIHM->NumEle; i++)
        {
            metarr = &(BGCM->grid[i].metarr);
            metarr->dayl[j] = dayl;
            if (j == 0)
                metarr->prev_dayl[j] = prev_dayl;
            else
                metarr->prev_dayl[j] = metarr->prev_dayl[j - 1];

            metarr->prcp[j] = 0.;
            metarr->tmax[j] = -999.;
            metarr->tmin[j] = 999.;
            metarr->tavg[j] = 0.;
            metarr->tday[j] = 0.;
            metarr->tnight[j] = 0.;
            metarr->vpd[j] = 0.;
            metarr->swavgfd[j] = 0.;
            metarr->par[j] = 0.;
            metarr->pa[j] = 0.;
            RH = 0.;
            daylight_coutner = 0;

            for (hour = 0; hour < 24; hour++)
            {
                metarr->prcp[j] = metarr->prcp[j] + PIHM_forcing[PIHM->Ele[i].meteo - 1][hour][PRCP_TS] * 3600.;    /* Convert from kg m-2 s-1 to kg m-2 */
                sfctmp = PIHM_forcing[PIHM->Ele[i].meteo - 1][hour][SFCTMP_TS] - 273.15;
                metarr->tmax[j] = sfctmp > metarr->tmax[j] ? sfctmp : metarr->tmax[j];
                metarr->tmin[j] = sfctmp < metarr->tmin[j] ? sfctmp : metarr->tmin[j];
                metarr->tavg[j] = metarr->tavg[j] + sfctmp;
                if (LSM->RAD_MODE > 0)
                    solar = topo_radiation (rad_forcing[PIHM->Ele[i].meteo - 1][hour][0], rad_forcing[PIHM->Ele[i].meteo - 1][hour][1], rad_forcing[PIHM->Ele[i].meteo - 1][hour][2], rad_forcing[PIHM->Ele[i].meteo - 1][hour][3], LSM->GRID[i].SLOPE, LSM->GRID[i].ASPECT, LSM->GRID[i].H_PHI, LSM->GRID[i].SVF);
                else
                    solar = PIHM_forcing[PIHM->Ele[i].meteo - 1][hour][SOLAR_TS];

                if (PIHM_forcing[PIHM->Ele[i].meteo - 1][hour][SOLAR_TS] > 0)
                {
                    metarr->tday[j] = metarr->tday[j] + sfctmp;
                    RH = PIHM_forcing[PIHM->Ele[i].meteo - 1][hour][RH_TS] / 100.;
                    metarr->vpd[j] = metarr->vpd[j] + (1. - RH) * 611.2 * exp (17.67 * sfctmp / (sfctmp + 243.5));
                    metarr->pa[j] = metarr->pa[j] + PIHM_forcing[PIHM->Ele[i].meteo - 1][hour][PRES_TS];
                    metarr->swavgfd[j] = metarr->swavgfd[j] + solar;
                    daylight_coutner++;
                }
                else
                    metarr->tnight[j] = metarr->tnight[j] + sfctmp;
            }

            metarr->tavg[j] = metarr->tavg[j] / 24.;
            metarr->tday[j] = metarr->tday[j] / (double)daylight_coutner;
            metarr->tnight[j] = metarr->tnight[j] / (24. - (double)daylight_coutner);
            metarr->vpd[j] = metarr->vpd[j] / (double)daylight_coutner;
            metarr->pa[j] = metarr->pa[j] / (double)daylight_coutner;
            metarr->swavgfd[j] = metarr->swavgfd[j] / (double)daylight_coutner;

            metarr->par[j] = metarr->swavgfd[j] * RAD2PAR;


            metarr->tsoil[j] = stc[i] - 273.15;
            metarr->swc[j] = swc[i];
            metarr->soilw[j] = soilw[j] * 1000.0;
            for (k = 0; k < 3; k++)
            {
                metarr->subflux[k][j] = subflux[k][j];
            }
            //if (i==0) printf ("%lf %lf\t", metarr->tsoil[j], metarr->swc[j]);
            //metarr->tsoil[j] = Interpolation (&BGCM->Forcing[STC_TS][i], t) - 273.15;
            //metarr->swc[j] = Interpolation (&BGCM->Forcing[SWC_TS][i], t);
        }
        //printf ("\n");
    }
}
