#include "bgc.h"

void daymet (bgc_struct BGCM, Model_Data PIHM, LSM_STRUCT LSM, double t, int spinup)
{
    double          sfctmp;
    double          solar;
    double          dayl, prev_dayl;
    double          RH;
    metvar_struct  *metv;
    int             hour;
    int             i;
    int             daylight_coutner;
    spa_data        spa;
    int             spa_result;
    time_t         *rawtime;
    struct tm      *timestamp;

    rawtime = (time_t *) malloc (sizeof (time_t));

    *rawtime = (int)t;
    timestamp = gmtime (rawtime);
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

    *rawtime = *rawtime - 24. * 3600.;
    timestamp = gmtime (rawtime);
    spa.year = timestamp->tm_year + 1900;
    spa.month = timestamp->tm_mon + 1;
    spa.day = timestamp->tm_mday;
    spa.hour = timestamp->tm_hour;
    spa.minute = timestamp->tm_min;
    spa.second = timestamp->tm_sec;
    spa_result = spa_calculate (&spa);
    prev_dayl = (spa.sunset - spa.sunrise) * 3600.;
    prev_dayl = prev_dayl < 0. ? (prev_dayl + 12. * 3600.) : prev_dayl;

    free (rawtime);

//    for (i = 0; i < PIHM->NumEle; i++)
    for (i = 0; i < 10; i++)
    {
        metv = &(BGCM->grid[i].metv);
        if (spinup == 1)
        {
            metv->dayl = dayl;
            metv->prev_dayl = prev_dayl;

            metv->prcp = 0.;
            metv->tmax = -999.;
            metv->tmin = 999.;
            metv->tavg = 0.;
            metv->tday = 0.;
            metv->tnight = 0.;
            metv->vpd = 0.;
            metv->swavgfd = 0.;
            metv->par = 0.;
            metv->pa = 0.;
            RH = 0.;
            daylight_coutner = 0;

            for (hour = 0; hour < 24; hour++)
            {
                metv->prcp = metv->prcp + Interpolation (&PIHM->Forcing[PRCP_TS][PIHM->Ele[i].prep - 1], t + hour * 3600.) * 3600.; /* Convert from kg m-2 s-1 to kg m-2 */
                sfctmp = Interpolation (&PIHM->Forcing[SFCTMP_TS][PIHM->Ele[i].temp - 1], t + hour * 3600.) - 273.15;
                metv->tmax = sfctmp > metv->tmax ? sfctmp : metv->tmax;
                metv->tmin = sfctmp < metv->tmin ? sfctmp : metv->tmin;
                metv->tavg = metv->tavg + sfctmp;
                solar = Interpolation (&PIHM->Forcing[SOLAR_TS][PIHM->Ele[i].Sdown - 1], t + hour * 3600.);
                metv->swavgfd = metv->swavgfd + solar;
                if (solar > 0)
                {
                    metv->tday = metv->tday + sfctmp;
                    RH = Interpolation (&PIHM->Forcing[RH_TS][PIHM->Ele[i].humidity - 1], t + hour * 3600.) / 100.;
                    metv->vpd = metv->vpd + (1. - RH) * 611.2 * exp (17.67 * sfctmp / (sfctmp + 243.5));
                    metv->pa = metv->pa + Interpolation (&PIHM->Forcing[PRES_TS][PIHM->Ele[i].pressure - 1], t + hour * 3600.);
                    daylight_coutner++;
                }
                else
                    metv->tnight = metv->tnight + sfctmp;
            }

            metv->tavg = metv->tavg / 24.;
            metv->tday = metv->tday / (double)daylight_coutner;
            metv->tnight = metv->tnight / (24. - (double)daylight_coutner);

            metv->swavgfd = metv->swavgfd / 24.;
            metv->par = metv->swavgfd * RAD2PAR;

            metv->vpd = metv->vpd / (double)daylight_coutner;
            metv->pa = metv->pa / (double)daylight_coutner;

            metv->tsoil = Interpolation (&BGCM->Forcing[STC_TS][i], t) - 273.15;
            metv->swc = Interpolation (&BGCM->Forcing[SWC_TS][i], t);
        }
        else
        {
            metv->dayl = dayl;
            metv->prev_dayl = prev_dayl;
            //        metv->prcp = ;
            //        metv->tmax = ;
            //        metv->tmin = ;
            //        metv->tavg = ;
            //        metv->tday = ;
            //        metv->tnight = ;
            //        metv->tsoil = ;
            //          metv->swc = ;
            //        metv->swavgfd = ;
            //        metv->par = ;
            //        metv->vpd = ;
            //        metv->pa = ;
        }

        //    printf ("%d: dayl = %lf, prev_dayl = %f, prcp = %lf, tmax = %lf, tmin = %lf, tavg = %lf, tday = %lf, tnight = %lf, tsoil = %lf, swc = %lf, swavgfd = %lf, par = %lf, vpd = %lf, pa = %lf\n", i, metv->dayl, metv->prev_dayl, metv->prcp, metv->tmax, metv->tmin, metv->tavg, metv->tday, metv->tnight, metv->tsoil, metv->swc, metv->swavgfd, metv->par, metv->vpd, metv->pa);
    }
}
