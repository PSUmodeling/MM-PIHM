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
    int             i, k;
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
        sfctmp = *elem->forc.meteo[SFCTMP_TS];
        daily->tmax = (daily->tmax > sfctmp) ? daily->tmax : sfctmp;
        daily->tmin = (daily->tmin < sfctmp) ? daily->tmin : sfctmp;
        daily->sfctmp += sfctmp;

        /* Solar radiation */
        solar = elem->ef.soldn;

        /* Wind speed */
        daily->sfcspd += elem->ps.sfcspd;

        /* Soil temperature */
        for (k = 0; k < elem->ps.nsoil; k++)
        {
            daily->stc[k] += elem->es.stc[k];
            daily->sh2o[k] += elem->ws.sh2o[k];
            daily->smflx[k] += elem->wf.smflx[k];
            //daily->smc[k] += elem->ws.smc[k];
        }

        /* Water storage terms */
        daily->surf += elem->ws.surf;
        daily->unsat += elem->ws.unsat;
        daily->gw += elem->ws.gw;

        /* Lateral flux */
        for (k = 0; k < 3; k++)
        {
            daily->fluxsub[k] += elem->wf.fluxsub[k];
            daily->fluxsurf[k] += elem->wf.fluxsurf[k];
        }
        daily->fluxsub[3] = 0.0;
        daily->fluxsurf[3] = 0.0;

        if (solar > 1.0)
        {
            daily->tday += sfctmp;
            daily->q2d += elem->ps.q2sat - elem->ps.q2;
            daily->ch += elem->ps.ch;
            daily->sfcprs += elem->ps.sfcprs;
            daily->albedo += elem->ps.albedo;
            daily->solar += solar;
            daily->solar_total += solar * pihm->ctrl.stepsize;
            (daily->daylight_counter)++;
        }
        else
        {
            daily->tnight += sfctmp;
        }

        (daily->counter)++;
    }

    /* River segments */
    for (i = 0; i < pihm->numriv; i++)
    {
        daily = &pihm->riv[i].daily;
        riv = &pihm->riv[i];

        /* Water storage terms */
        daily->surf += riv->ws.stage;
        daily->gw += riv->ws.gw;

        /* Lateral flux */
        daily->fluxsub[0] += riv->wf.fluxriv[0];
        daily->fluxsurf[0] += riv->wf.fluxriv[10];

        daily->fluxsub[1] += riv->wf.fluxriv[1];
        daily->fluxsurf[1] += riv->wf.fluxriv[9];

        daily->fluxsub[2] += riv->wf.fluxriv[2];
        daily->fluxsurf[2] += riv->wf.fluxriv[4] + riv->wf.fluxriv[7];

        daily->fluxsub[3] += riv->wf.fluxriv[3];
        daily->fluxsurf[3] += riv->wf.fluxriv[5] + riv->wf.fluxriv[8];

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
        spa.pressure = 1013.25 * pow ((293. - 0.0065 * spa.elevation) / 293.0, 5.26);
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
        prev_dayl = (prev_dayl < 0.0) ? (prev_dayl + 12.0 * 3600.0) : prev_dayl;

        for (i = 0; i < pihm->numele; i++)
        {
            daily = &(pihm->elem[i].daily);

            daily->sfctmp /= (double)daily->counter;
            daily->dayl = dayl;
            daily->prev_dayl = prev_dayl;

            daily->sfcspd /= (double) daily->counter;

            for (k = 0; k < pihm->elem[i].ps.nsoil; k++)
            {
                daily->stc[k] /= (double)daily->counter;
                daily->sh2o[k] /= (double)daily->counter;
                daily->smflx[k] /= (double)daily->counter;
            }
            daily->surf /= (double)daily->counter;
            daily->unsat /= (double)daily->counter;
            daily->gw /= (double)daily->counter;

            for (k = 0; k < 4; k++)
            {
                daily->fluxsub[k] /= (double)daily->counter;
                daily->fluxsurf[k] /= (double)daily->counter;
            }

            daily->tday /= (double)daily->daylight_counter;
            daily->q2d /= (double)daily->daylight_counter;
            daily->sfcprs /= (double)daily->daylight_counter;
            daily->ch /= (double)daily->daylight_counter;
            daily->albedo /= (double)daily->daylight_counter;
            daily->solar /= (double)daily->daylight_counter;

            daily->tnight /= (double)(daily->counter - daily->daylight_counter);
        }

        for (i = 0; i < pihm->numriv; i++)
        {
            daily = &(pihm->riv[i].daily);

            daily->surf /= (double)daily->counter;
            daily->gw /= (double)daily->counter;

            for (k = 0; k < 4; k++)
            {
                daily->fluxsub[k] /= (double)daily->counter;
                daily->fluxsurf[k] /= (double)daily->counter;
            }
        }
    }
}

void InitDailyStruct (pihm_struct pihm)
{
    int             i, k;
    daily_struct   *daily;

    for (i = 0; i < pihm->numele; i++)
    {
        daily = &pihm->elem[i].daily;

        daily->counter = 0;
        daily->daylight_counter = 0;

        daily->sfctmp = 0.0;
        daily->tday = 0.0;
        daily->tnight = 0.0;
        daily->tmax = -999.0;
        daily->tmin = 999.0;
        daily->sfcspd = 0.0;
        daily->solar = 0.0;
        daily->solar_total = 0.0;
        daily->sfcprs = 0.0;
        daily->surf = 0.0;
        daily->unsat = 0.0;
        daily->gw = 0.0;
        daily->infil = 0.0;
        daily->rechg = 0.0;
        daily->dayl = 0.0;
        daily->prev_dayl = 0.0;
        for (k = 0; k < MAXLYR; k++)
        {
            daily->stc[k] = 0.0;
            daily->sh2o[k] = 0.0;
        }
        daily->q2d = 0.0;
        daily->albedo = 0.0;
        daily->ch = 0.0;
        for (k = 0; k < 4; k++)
        {
            daily->fluxsurf[k] = 0.0;
            daily->fluxsub[k] = 0.0;
        }
    }

    for (i = 0; i < pihm->numriv; i++)
    {
        daily = &pihm->riv[i].daily;

        daily->counter = 0;
        daily->daylight_counter = 0;

        daily->sfctmp = 0.0;
        daily->tday = 0.0;
        daily->tnight = 0.0;
        daily->tmax = -999.0;
        daily->tmin = 999.0;
        daily->sfcspd = 0.0;
        daily->solar = 0.0;
        daily->solar_total = 0.0;
        daily->sfcprs = 0.0;
        daily->surf = 0.0;
        daily->unsat = 0.0;
        daily->gw = 0.0;
        daily->infil = 0.0;
        daily->rechg = 0.0;
        daily->dayl = 0.0;
        daily->prev_dayl = 0.0;
        for (k = 0; k < MAXLYR; k++)
        {
            daily->smflx[MAXLYR] = 0.0;
            daily->stc[k] = 0.0;
            daily->sh2o[k] = 0.0;
        }
        daily->q2d = 0.0;
        daily->albedo = 0.0;
        daily->ch = 0.0;
        for (k = 0; k < 4; k++)
        {
            daily->fluxsurf[k] = 0.0;
            daily->fluxsub[k] = 0.0;
        }
    }
}
