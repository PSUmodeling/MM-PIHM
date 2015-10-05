#include "pihm.h"
#ifdef _NOAH_
#include "noah.h"
#include "spa.h"
#endif

#ifdef _NOAH_
void DailyVar (int t, int start_time, pihm_struct pihm, lsm_struct noah)
#else
void DailyVar (int t, int start_time, pihm_struct pihm)
#endif
{
    daily_struct   *daily;
    elem_struct    *elem;
#ifdef _NOAH_
    grid_struct    *grid;
#endif
    int             i, k;
    double          sfctmp;
    double          solar;
    double          dayl;
    double          prev_dayl;
    time_t          rawtime;
    struct tm      *timestamp;
    spa_data        spa;
    int             spa_result;

    for (i = 0; i < pihm->numele; i++)
    {
        daily = &pihm->elem[i].daily;
        elem = &pihm->elem[i];
#ifdef _NOAH_
        grid = &noah->grid[i];
#endif

        /* Initialize daily structure at the beginning of simulation */
        if (t == start_time)
        {
            InitDailyStruct (daily);
        }

        /* Calculate daily variables */
        /* Air temperature */
        sfctmp = *elem->forc.meteo[SFCTMP_TS];
        daily->tmax = (daily->tmax > sfctmp) ? daily->tmax : sfctmp;
        daily->tmin = (daily->tmin < sfctmp) ? daily->tmin : sfctmp;
        daily->sfctmp += sfctmp;

        /* Solar radiation */
#ifdef _NOAH_
        solar = grid->soldn;
#else
        solar = *elem->forc.meteo[SOLAR_TS];
#endif

#ifdef _NOAH_
        /* Soil temperature */
        daily->stc += grid->stc[0];
#endif

#ifdef _NOAH_
        /* Root zone soil moisture */
        daily->sh2o += grid->soilw;
#endif

        /* Water storage terms */
        daily->surf += elem->surf;
        daily->unsat += elem->unsat;
        daily->gw += elem->gw;

        /* Lateral flux */
        for (k = 0; k < 3; k++)
        {
            daily->fluxsub[k] += elem->fluxsub[k];
            daily->fluxsurf[k] += elem->fluxsurf[k];
        }
        daily->fluxsub[3] = 0.0;
        daily->fluxsurf[3] = 0.0;

        if (solar > 1.0)
        {
            daily->tday += sfctmp;
            daily->q2d += grid->q2sat - grid->q2;
            daily->sfcprs += grid->sfcprs;
            daily->solar += solar;
            (daily->daylight_counter)++;
        }
        else
        {
            daily->tnight += sfctmp;
        }

        (daily->counter)++;
    }

    /* Calculate daily variables */
    if ((t - start_time) % 86400 == 0 && t > start_time)
    {
#ifdef _NOAH_
        rawtime = (int)(t - 86400);
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

        spa.longitude = noah->longitude;
        spa.latitude = noah->latitude;
        spa.elevation = 0.;
        for (i = 0; i < pihm->numele; i++)
        {
            spa.elevation = spa.elevation + (double)pihm->elem[i].topo.zmax;
        }
        spa.elevation = spa.elevation / (double)pihm->numele;
        /*
         * Calculate surface pressure based on FAO 1998 method (Narasimhan 2002) 
         */
        spa.pressure = 1013.25 * pow ((293. - 0.0065 * spa.elevation) / 293.0, 5.26);
        spa.temperature = noah->genprmt.tbot_data;

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
#endif

        for (i = 0; i < pihm->numele; i++)
        {
            daily = &(pihm->elem[i].daily);

            daily->dayl = dayl;
            daily->prev_dayl = prev_dayl;

            daily->sfctmp /= (double)daily->counter;
            daily->stc /= (double)daily->counter;
            daily->sh2o /= (double)daily->counter;
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
            daily->solar /= (double)daily->daylight_counter;

            daily->tnight /= (double)(daily->counter - daily->daylight_counter);
        }
    }
}

void InitDailyStruct (daily_struct *daily)
{
    int             k;
    daily->counter = 0;
    daily->daylight_counter = 0;

#ifdef _NOAH_
    daily->tmax = -999.0;
    daily->tmin = 999.0;
    daily->sfctmp = 0.0;
    daily->stc = 0.0;
    daily->sh2o = 0.0;
    daily->tday = 0.0;
    daily->tnight = 0.0;
    daily->q2d = 0.0;
    daily->sfcprs = 0.0;
    daily->solar = 0.0;
#endif
    for (k = 0; k < 4; k++)
    {
        daily->fluxsurf[k] = 0.0;
        daily->fluxsub[k] = 0.0;
    }
}
