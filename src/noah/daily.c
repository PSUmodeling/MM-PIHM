#include "pihm.h"

void DailyVar (int t, int start_time, pihm_struct pihm)
{
    double          dayl;
    double          prev_dayl;
    spa_data        spa;
    int             spa_result;
    time_t          rawtime;
    struct tm      *timestamp;
    int             i;

    /*
     * Cumulates daily variables
     */
    /* Triangular grids */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < pihm->numele; i++)
    {
        double      sfctmp;
        double      solar;
        int         k;
        elem_struct *elem;

        elem = &pihm->elem[i];

        /* Air temperature */
        sfctmp = elem->es.sfctmp;
        elem->daily.avg_sfctmp += sfctmp;
        elem->daily.tmax =
            (elem->daily.tmax > sfctmp) ? elem->daily.tmax : sfctmp;
        elem->daily.tmin =
            (elem->daily.tmin < sfctmp) ? elem->daily.tmin : sfctmp;

        /* Solar radiation */
        solar = elem->ef.soldn;

        /* Wind speed */
        elem->daily.avg_sfcspd += elem->ps.sfcspd;

        /* Soil moisture, temperature, and ET */
        for (k = 0; k < elem->ps.nsoil; k++)
        {
            elem->daily.avg_stc[k] += elem->es.stc[k];
            elem->daily.avg_sh2o[k] += elem->ws.sh2o[k];
            elem->daily.avg_smc[k] += elem->ws.smc[k];
            elem->daily.avg_smflxv[k] += elem->wf.smflxv[k];
#ifdef _CYCLES_
            elem->daily.avg_et[k] += elem->wf.et[k];
#endif
        }

#ifdef _CYCLES_
        elem->daily.avg_sncovr += elem->ps.sncovr;
#endif

        /* Water storage terms */
        elem->daily.avg_surf += elem->ws.surf;
        elem->daily.avg_unsat += elem->ws.unsat;
        elem->daily.avg_gw += elem->ws.gw;

        /* Lateral flux */
        for (k = 0; k < NUM_EDGE; k++)
        {
            elem->daily.avg_subsurf[k] += elem->wf.subsurf[k];
            elem->daily.avg_ovlflow[k] += elem->wf.ovlflow[k];
        }

        if (solar > 0.0)
        {
            elem->daily.tday += sfctmp;
            elem->daily.avg_q2d += elem->ps.q2sat - elem->ps.q2;
            elem->daily.avg_ch += elem->ps.ch;
            elem->daily.avg_sfcprs += elem->ps.sfcprs;
            elem->daily.avg_albedo += elem->ps.albedo;
            elem->daily.avg_soldn += solar;
            elem->daily.solar_total += solar * pihm->ctrl.stepsize;
            (elem->daily.daylight_counter)++;
        }
        else
        {
            elem->daily.tnight += sfctmp;
        }

        (elem->daily.counter)++;
    }

    /* River segments */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < pihm->numriv; i++)
    {
        river_struct *riv;
        int         j;

        riv = &pihm->riv[i];

        /* Water storage terms */
        riv->daily.avg_stage += riv->ws.stage;
        riv->daily.avg_gw += riv->ws.gw;

        /* Lateral flux */
        for (j = 0; j < NUM_RIVFLX; j++)
        {
            riv->daily.avg_rivflow[j] += riv->wf.rivflow[j];
        }

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

        /* Calculate surface pressure based on FAO 1998 method
         * (Narasimhan 2002) */
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

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < pihm->numele; i++)
        {
            int     k;
            elem_struct *elem;

            elem = &pihm->elem[i];

            elem->daily.avg_sfctmp /= (double)elem->daily.counter;
            elem->daily.dayl = dayl;
            elem->daily.prev_dayl = prev_dayl;

            elem->daily.avg_sfcspd /= (double)elem->daily.counter;

            for (k = 0; k < pihm->elem[i].ps.nsoil; k++)
            {
                elem->daily.avg_stc[k] /= (double)elem->daily.counter;
                elem->daily.avg_sh2o[k] /= (double)elem->daily.counter;
                elem->daily.avg_smc[k] /= (double)elem->daily.counter;
                elem->daily.avg_smflxv[k] /= (double)elem->daily.counter;
#ifdef _CYCLES_
                elem->daily.avg_et[k] /= (double)elem->daily.counter;
#endif
            }

#ifdef _CYCLES_
            elem->daily.avg_sncovr /= (double)elem->daily.counter;
#endif

            elem->daily.avg_surf /= (double)elem->daily.counter;
            elem->daily.avg_unsat /= (double)elem->daily.counter;
            elem->daily.avg_gw /= (double)elem->daily.counter;

            for (k = 0; k < NUM_EDGE; k++)
            {
                elem->daily.avg_subsurf[k] /= (double)elem->daily.counter;
                elem->daily.avg_ovlflow[k] /= (double)elem->daily.counter;
            }

            elem->daily.tday /= (double)elem->daily.daylight_counter;
            elem->daily.avg_q2d /= (double)elem->daily.daylight_counter;
            elem->daily.avg_ch /= (double)elem->daily.daylight_counter;
            elem->daily.avg_sfcprs /= (double)elem->daily.daylight_counter;
            elem->daily.avg_albedo /= (double)elem->daily.daylight_counter;
            elem->daily.avg_soldn /= (double)elem->daily.daylight_counter;

            elem->daily.tnight /= (double)(elem->daily.counter -
                elem->daily.daylight_counter);
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < pihm->numriv; i++)
        {
            int     j;
            river_struct *riv;

            riv = &pihm->riv[i];

            riv->daily.avg_stage /= (double)riv->daily.counter;
            riv->daily.avg_gw /= (double)riv->daily.counter;

            for (j = 0; j < NUM_RIVFLX; j++)
            {
                riv->daily.avg_rivflow[j] /= (double)riv->daily.counter;
            }
        }
    }
}

void InitDailyStruct (pihm_struct pihm)
{
    int             i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < pihm->numele; i++)
    {
        int         k;
        elem_struct *elem;

        elem = &pihm->elem[i];

        elem->daily.counter = 0;
        elem->daily.daylight_counter = 0;

        elem->daily.avg_surf = 0.0;
        elem->daily.avg_unsat = 0.0;
        elem->daily.avg_gw = 0.0;
        for (k = 0; k < MAXLYR; k++)
        {
            elem->daily.avg_sh2o[k] = 0.0;
            elem->daily.avg_smc[k] = 0.0;
            elem->daily.avg_et[k] = 0.0;
            elem->daily.avg_smflxv[k] = 0.0;
            elem->daily.avg_stc[k] = 0.0;
        }

        for (k = 0; k < NUM_EDGE; k++)
        {
            elem->daily.avg_ovlflow[k] = 0.0;
            elem->daily.avg_subsurf[k] = 0.0;
        }

        elem->daily.dayl = BADVAL;
        elem->daily.prev_dayl = BADVAL;
        elem->daily.avg_q2d = 0.0;
        elem->daily.avg_sfcprs = 0.0;
        elem->daily.avg_ch = 0.0;
        elem->daily.avg_albedo = 0.0;
        elem->daily.avg_sfcspd = 0.0;

        elem->daily.tmax = -999.0;
        elem->daily.tmin = 999.0;
        elem->daily.avg_sfctmp = 0.0;
        elem->daily.tday = 0.0;
        elem->daily.tnight = 0.0;

        elem->daily.avg_soldn = 0.0;
        elem->daily.solar_total = 0.0;
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < pihm->numriv; i++)
    {
        int         k;
        river_struct *riv;

        riv = &pihm->riv[i];

        riv->daily.counter = 0;

        riv->daily.avg_stage = 0.0;
        riv->daily.avg_gw = 0.0;
        for (k = 0; k < NUM_RIVFLX; k++)
        {
            riv->daily.avg_rivflow[k] = 0.0;
        }
    }
}
