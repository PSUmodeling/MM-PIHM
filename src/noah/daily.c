#include "pihm.h"

void DailyVar (int t, int start_time, pihm_struct pihm)
{
    daily_struct   *daily;
    elem_struct    *elem;
    river_struct   *riv;
    double          dayl;
    double          prev_dayl;
    spa_data        spa;
    int             spa_result;
    time_t          rawtime;
    struct tm      *timestamp;
    int             i, j, k;
    double          sfctmp;
    double          solar;

    /*
     * Cumulates daily variables
     */
    /* Triangular grids */
    for (i = 0; i < pihm->numele; i++)
    {
        daily = &pihm->elem[i].daily;
        elem = &pihm->elem[i];

        /* Air temperature */
        sfctmp = elem->es.sfctmp;
        daily->es.sfctmp += sfctmp;
        daily->es.tmax = (daily->es.tmax > sfctmp) ? daily->es.tmax : sfctmp;
        daily->es.tmin = (daily->es.tmin < sfctmp) ? daily->es.tmin : sfctmp;

        /* Solar radiation */
        solar = elem->ef.soldn;

        /* Wind speed */
        daily->ps.sfcspd += elem->ps.sfcspd;

        /* Soil moisture, temperature, and ET */
        for (k = 0; k < elem->ps.nsoil; k++)
        {
            daily->es.stc[k] += elem->es.stc[k];
            daily->ws.sh2o[k] += elem->ws.sh2o[k];
            daily->wf.smflxv[k] += elem->wf.smflxv[k];
        }

#ifdef _CYCLES_
        for (k = 0; k < elem->ps.nsoil; k++)
        {
            daily->wf.et[k] += elem->wf.et[k];
        }
        daily->ps.sncovr += elem->ps.sncovr;
#endif

        /* Water storage terms */
        daily->ws.surf += elem->ws.surf;
        daily->ws.unsat += elem->ws.unsat;
        daily->ws.gw += elem->ws.gw;


        /* Lateral flux */
        for (k = 0; k < 3; k++)
        {
            daily->wf.subsurf[k] += elem->wf.subsurf[k];
            daily->wf.surf[k] += elem->wf.surf[k];
        }
        daily->wf.subsurf[3] = 0.0;
        daily->wf.surf[3] = 0.0;

        if (solar > 1.0)
        {
            daily->es.tday += sfctmp;
            daily->ps.q2d += elem->ps.q2sat - elem->ps.q2;
            daily->ps.ch += elem->ps.ch;
            daily->ps.sfcprs += elem->ps.sfcprs;
            daily->ps.albedo += elem->ps.albedo;
            daily->ef.soldn += solar;
            daily->ef.solar_total += solar * pihm->ctrl.stepsize;
            (daily->daylight_counter)++;
        }
        else
        {
            daily->es.tnight += sfctmp;
        }

        (daily->counter)++;
    }

    /* River segments */
    for (i = 0; i < pihm->numriv; i++)
    {
        daily = &pihm->riv[i].daily;
        riv = &pihm->riv[i];

        /* Water storage terms */
        daily->ws.surf += riv->ws.stage;
        daily->ws.gw += riv->ws.gw;

        /* Lateral flux */
        for (j = 0; j < 11; j++)
        {
            daily->wf.river[j] += riv->wf.river[j];
        }
        //daily->wf.fluxsub[0] += riv->wf.fluxriv[0];
        //daily->wf.fluxsurf[0] += riv->wf.fluxriv[10];

        //daily->wf.fluxsub[1] += riv->wf.fluxriv[1];
        //daily->wf.fluxsurf[1] += riv->wf.fluxriv[9];

        //daily->wf.fluxsub[2] += riv->wf.fluxriv[2];
        //daily->wf.fluxsurf[2] += riv->wf.fluxriv[4] + riv->wf.fluxriv[7];

        //daily->wf.fluxsub[3] += riv->wf.fluxriv[3];
        //daily->wf.fluxsurf[3] += riv->wf.fluxriv[5] + riv->wf.fluxriv[8];

        (daily->counter)++;
    }

    /* Calculate daily variables */
    if ((t - start_time) % DAYINSEC == 0 && t > start_time)
    {
        rawtime = (int)(t - DAYINSEC);
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

        spa.longitude = pihm->longitude;
        spa.latitude = pihm->latitude;
        spa.elevation = pihm->elevation;

        /*
         * Calculate surface pressure based on FAO 1998 method (Narasimhan 2002) 
         */
        spa.pressure =
            1013.25 * pow ((293. - 0.0065 * spa.elevation) / 293.0, 5.26);
        spa.temperature = pihm->noahtbl.tbot;

        spa.function = SPA_ZA_RTS;
        spa_result = spa_calculate (&spa);

        /* daylength (s) */
        dayl = (spa.sunset - spa.sunrise) * 3600.0;
        dayl = (dayl < 0.0) ? (dayl + 24.0 * 3600.0) : dayl;

        rawtime = rawtime - 24 * 3600;
        timestamp = gmtime (&rawtime);
        spa.year = timestamp->tm_year + 1900;
        spa.month = timestamp->tm_mon + 1;
        spa.day = timestamp->tm_mday;
        spa.hour = timestamp->tm_hour;
        spa.minute = timestamp->tm_min;
        spa.second = timestamp->tm_sec;
        spa_result = spa_calculate (&spa);
        prev_dayl = (spa.sunset - spa.sunrise) * 3600.;
        prev_dayl =
            (prev_dayl < 0.0) ? (prev_dayl + 12.0 * 3600.0) : prev_dayl;

        for (i = 0; i < pihm->numele; i++)
        {
            daily = &(pihm->elem[i].daily);

            daily->es.sfctmp /= (double)daily->counter;
            daily->ps.dayl = dayl;
            daily->ps.prev_dayl = prev_dayl;

            daily->ps.sfcspd /= (double)daily->counter;

            for (k = 0; k < pihm->elem[i].ps.nsoil; k++)
            {
                daily->es.stc[k] /= (double)daily->counter;
                daily->ws.sh2o[k] /= (double)daily->counter;
                daily->wf.smflxv[k] /= (double)daily->counter;
            }

#ifdef _CYCLES_
            for (k = 0; k < pihm->elem[i].ps.nsoil; k++)
            {
                daily->wf.et[k] /= (double)daily->counter;
            }
            daily->ps.sncovr /= (double)daily->counter;
#endif

            daily->ws.surf /= (double)daily->counter;
            daily->ws.unsat /= (double)daily->counter;
            daily->ws.gw /= (double)daily->counter;

            for (k = 0; k < 4; k++)
            {
                daily->wf.subsurf[k] /= (double)daily->counter;
                daily->wf.surf[k] /= (double)daily->counter;
            }

            daily->es.tday /= (double)daily->daylight_counter;
            daily->ps.q2d /= (double)daily->daylight_counter;
            daily->ps.ch /= (double)daily->daylight_counter;
            daily->ps.sfcprs /= (double)daily->daylight_counter;
            daily->ps.albedo /= (double)daily->daylight_counter;
            daily->ef.soldn /= (double)daily->daylight_counter;

            daily->es.tnight /= (double)(daily->counter -
                daily->daylight_counter);
        }

        for (i = 0; i < pihm->numriv; i++)
        {
            daily = &(pihm->riv[i].daily);

            daily->ws.surf /= (double)daily->counter;
            daily->ws.gw /= (double)daily->counter;

            for (j = 0; j < 11; j++)
            {
                daily->wf.river[j] /= (double)daily->counter;
            }
        }
    }
}

void InitDailyStruct (pihm_struct pihm)
{
    int             i;
    daily_struct   *daily;

    for (i = 0; i < pihm->numele; i++)
    {
        daily = &pihm->elem[i].daily;

        daily->counter = 0;
        daily->daylight_counter = 0;
        InitWState (&daily->ws);
        InitWFlux (&daily->wf);
        InitPState (&daily->ps);
        InitEState (&daily->es);
        InitEFlux (&daily->ef);
    }

    for (i = 0; i < pihm->numriv; i++)
    {
        daily = &pihm->riv[i].daily;

        daily->counter = 0;
        daily->daylight_counter = 0;
        InitWState (&daily->ws);
        InitWFlux (&daily->wf);
        InitPState (&daily->ps);
        InitEState (&daily->es);
        InitEFlux (&daily->ef);
    }
}
