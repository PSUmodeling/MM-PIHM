#include "pihm.h"

void DailyVar (int t, int start_time, pihm_struct pihm)
{
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
        elem = &pihm->elem[i];

        /* Air temperature */
        sfctmp = elem->es.sfctmp;
        elem->daily.es.avg_sfctmp += sfctmp;
        elem->daily.es.tmax =
            (elem->daily.es.tmax > sfctmp) ? elem->daily.es.tmax : sfctmp;
        elem->daily.es.tmin =
            (elem->daily.es.tmin < sfctmp) ? elem->daily.es.tmin : sfctmp;

        /* Solar radiation */
        solar = elem->ef.soldn;

        /* Wind speed */
        elem->daily.ps.avg_sfcspd += elem->ps.sfcspd;

        /* Soil moisture, temperature, and ET */
        for (k = 0; k < elem->ps.nsoil; k++)
        {
            elem->daily.es.avg_stc[k] += elem->es.stc[k];
            elem->daily.ws.avg_sh2o[k] += elem->ws.sh2o[k];
            elem->daily.wf.avg_smflxv[k] += elem->wf.smflxv[k];
        }

#ifdef _CYCLES_
        for (k = 0; k < elem->ps.nsoil; k++)
        {
            elem->daily.wf.avg_et[k] += elem->wf.et[k];
        }
        elem->daily.ps.avg_sncovr += elem->ps.sncovr;
#endif

        /* Water storage terms */
        elem->daily.ws.avg_surf += elem->ws.surf;
        elem->daily.ws.avg_unsat += elem->ws.unsat;
        elem->daily.ws.avg_gw += elem->ws.gw;

        /* Lateral flux */
        for (k = 0; k < 3; k++)
        {
            elem->daily.wf.avg_subsurf[k] += elem->wf.subsurf[k];
            elem->daily.wf.avg_surf[k] += elem->wf.surf[k];
        }

        if (solar > 0.0)
        {
            elem->daily.es.tday += sfctmp;
            elem->daily.ps.avg_q2d += elem->ps.q2sat - elem->ps.q2;
            elem->daily.ps.avg_ch += elem->ps.ch;
            elem->daily.ps.avg_sfcprs += elem->ps.sfcprs;
            elem->daily.ps.avg_albedo += elem->ps.albedo;
            elem->daily.ef.avg_soldn += solar;
            elem->daily.ef.solar_total += solar * pihm->ctrl.stepsize;
            (elem->daily.daylight_counter)++;
        }
        else
        {
            elem->daily.es.tnight += sfctmp;
        }

        (elem->daily.counter)++;
    }

    /* River segments */
    for (i = 0; i < pihm->numriv; i++)
    {
        riv = &pihm->riv[i];

        /* Water storage terms */
        riv->daily.ws.avg_stage += riv->ws.stage;
        riv->daily.ws.avg_gw += riv->ws.gw;

        /* Lateral flux */
        for (j = 0; j < 11; j++)
        {
            riv->daily.wf.avg_river[j] += riv->wf.river[j];
        }
        //daily->wf.fluxsub[0] += riv->wf.fluxriv[0];
        //daily->wf.fluxsurf[0] += riv->wf.fluxriv[10];

        //daily->wf.fluxsub[1] += riv->wf.fluxriv[1];
        //daily->wf.fluxsurf[1] += riv->wf.fluxriv[9];

        //daily->wf.fluxsub[2] += riv->wf.fluxriv[2];
        //daily->wf.fluxsurf[2] += riv->wf.fluxriv[4] + riv->wf.fluxriv[7];

        //daily->wf.fluxsub[3] += riv->wf.fluxriv[3];
        //daily->wf.fluxsurf[3] += riv->wf.fluxriv[5] + riv->wf.fluxriv[8];

        (riv->daily.counter)++;
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
            elem = &pihm->elem[i];

            elem->daily.es.avg_sfctmp /= (double)elem->daily.counter;
            elem->daily.ps.dayl = dayl;
            elem->daily.ps.prev_dayl = prev_dayl;

            elem->daily.ps.avg_sfcspd /= (double)elem->daily.counter;

            for (k = 0; k < pihm->elem[i].ps.nsoil; k++)
            {
                elem->daily.es.avg_stc[k] /= (double)elem->daily.counter;
                elem->daily.ws.avg_sh2o[k] /= (double)elem->daily.counter;
                elem->daily.wf.avg_smflxv[k] /= (double)elem->daily.counter;
            }

#ifdef _CYCLES_
            for (k = 0; k < pihm->elem[i].ps.nsoil; k++)
            {
                elem->daily.wf.avg_et[k] /= (double)elem->daily.counter;
            }
            elem->daily.ps.avg_sncovr /= (double)elem->daily.counter;
#endif

            elem->daily.ws.avg_surf /= (double)elem->daily.counter;
            elem->daily.ws.avg_unsat /= (double)elem->daily.counter;
            elem->daily.ws.avg_gw /= (double)elem->daily.counter;

            for (k = 0; k < 3; k++)
            {
                elem->daily.wf.avg_subsurf[k] /= (double)elem->daily.counter;
                elem->daily.wf.avg_surf[k] /= (double)elem->daily.counter;
            }

            elem->daily.es.tday /= (double)elem->daily.daylight_counter;
            elem->daily.ps.avg_q2d /= (double)elem->daily.daylight_counter;
            elem->daily.ps.avg_ch /= (double)elem->daily.daylight_counter;
            elem->daily.ps.avg_sfcprs /= (double)elem->daily.daylight_counter;
            elem->daily.ps.avg_albedo /= (double)elem->daily.daylight_counter;
            elem->daily.ef.avg_soldn /= (double)elem->daily.daylight_counter;

            elem->daily.es.tnight /= (double)(elem->daily.counter -
                elem->daily.daylight_counter);
        }

        for (i = 0; i < pihm->numriv; i++)
        {
            riv = &pihm->riv[i];

            riv->daily.ws.avg_stage /= (double)riv->daily.counter;
            riv->daily.ws.avg_gw /= (double)riv->daily.counter;

            for (j = 0; j < 11; j++)
            {
                riv->daily.wf.avg_river[j] /= (double)riv->daily.counter;
            }
        }
    }
}

void InitDailyStruct (pihm_struct pihm)
{
    int             i;
    elem_struct    *elem;
    river_struct   *riv;

    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        elem->daily.counter = 0;
        elem->daily.daylight_counter = 0;
        InitDailyWState (&elem->daily.ws);
        InitDailyWFlux (&elem->daily.wf);
        InitDailyPState (&elem->daily.ps);
        InitDailyEState (&elem->daily.es);
        InitDailyEFlux (&elem->daily.ef);
    }

    for (i = 0; i < pihm->numriv; i++)
    {
        riv = &pihm->riv[i];

        riv->daily.counter = 0;
        InitRiverDailyWState (&riv->daily.ws);
        InitRiverDailyWFlux (&riv->daily.wf);
    }
}

void InitDailyWState (daily_wstate_struct *ws)
{
    int             k;

    ws->avg_surf = 0.0;
    ws->avg_unsat = 0.0;
    ws->avg_gw = 0.0;
    for (k = 0; k < MAXLYR; k++)
    {
        ws->avg_sh2o[k] = 0.0;
    }
}

void InitDailyWFlux (daily_wflux_struct *wf)
{
    int             k;

    for (k = 0; k < NUM_EDGE; k++)
    {
        wf->avg_surf[k] = 0.0;
        wf->avg_subsurf[k] = 0.0;
    }

    for (k = 0; k < MAXLYR; k++)
    {
        wf->avg_et[k] = 0.0;
        wf->avg_smflxv[k] = 0.0;
    }
}

void InitDailyPState (daily_pstate_struct *ps)
{
    ps->dayl = BADVAL;
    ps->prev_dayl = BADVAL;
    ps->avg_q2d = 0.0;
    ps->avg_sfcprs = 0.0;
    ps->avg_ch = 0.0;
    ps->avg_albedo = 0.0;
    ps->avg_sfcspd = 0.0;
}

void InitDailyEState (daily_estate_struct *es)
{
    int             k;

    es->tmax = -999.0;
    es->tmin = 999.0;
    es->avg_sfctmp = 0.0;
    es->tday = 0.0;
    es->tnight = 0.0;
    for (k = 0; k < MAXLYR; k++)
    {
        es->avg_stc[k] = 0.0;
    }
}

void InitDailyEFlux (daily_eflux_struct *ef)
{
    ef->avg_soldn = 0.0;
    ef->par = 0.0;
    ef->solar_total = 0.0;
}

void InitRiverDailyWState (river_daily_wstate_struct *ws)
{
    ws->avg_stage = 0.0;
    ws->avg_gw = 0.0;
}

void InitRiverDailyWFlux (river_daily_wflux_struct *wf)
{
    int             k;

    for (k = 0; k < NUM_RIVFLX; k++)
    {
        wf->avg_river[k] = 0.0;
    }
}
