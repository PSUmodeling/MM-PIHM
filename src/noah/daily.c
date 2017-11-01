#include "pihm.h"

void DailyVar(int t, int start_time, elem_struct *elem, river_struct *river,
    double dt)
{
    int             i;

    /*
     * Sum daily variables
     */
#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             k;

        /* Air temperature */
        elem[i].daily.avg_sfctmp += elem[i].es.sfctmp;
        elem[i].daily.tmax = (elem[i].daily.tmax > elem[i].es.sfctmp) ?
            elem[i].daily.tmax : elem[i].es.sfctmp;
        elem[i].daily.tmin = (elem[i].daily.tmin < elem[i].es.sfctmp) ?
            elem[i].daily.tmin : elem[i].es.sfctmp;

        /* Wind speed */
        elem[i].daily.avg_sfcspd += elem[i].ps.sfcspd;

        /* Soil moisture, temperature, and ET */
        for (k = 0; k < elem[i].ps.nsoil; k++)
        {
            elem[i].daily.avg_stc[k] += elem[i].es.stc[k];
            elem[i].daily.avg_sh2o[k] += elem[i].ws.sh2o[k];
            elem[i].daily.avg_smc[k] += elem[i].ws.smc[k];
            elem[i].daily.avg_smflxv[k] += elem[i].wf.smflxv[k];
#ifdef _CYCLES_
            elem[i].daily.avg_et[k] += elem[i].wf.et[k];
#endif
        }

#ifdef _CYCLES_
        elem[i].daily.avg_sncovr += elem[i].ps.sncovr;
#endif

        /* Water storage terms */
        elem[i].daily.avg_surf += elem[i].ws.surf;
        elem[i].daily.avg_unsat += elem[i].ws.unsat;
        elem[i].daily.avg_gw += elem[i].ws.gw;

        if (elem[i].ef.soldn > 0.0)
        {
            elem[i].daily.tday += elem[i].es.sfctmp;
            elem[i].daily.avg_q2d += elem[i].ps.q2sat - elem[i].ps.q2;
            elem[i].daily.avg_ch += elem[i].ps.ch;
            elem[i].daily.avg_rc += elem[i].ps.rc;
            elem[i].daily.avg_sfcprs += elem[i].ps.sfcprs;
            elem[i].daily.avg_albedo += elem[i].ps.albedo;
            elem[i].daily.avg_soldn += elem[i].ef.soldn;
            elem[i].daily.solar_total += elem[i].ef.soldn * dt;
            (elem[i].daily.daylight_counter)++;
        }
        else
        {
            elem[i].daily.tnight += elem[i].es.sfctmp;
        }

        (elem[i].daily.counter)++;
    }

#ifdef _CYCLES_
    /* River segments */
# ifdef _OPENMP
#  pragma omp parallel for
# endif
    for (i = 0; i < nriver; i++)
    {
        int             j;

        /* Water storage terms */
        river[i].daily.avg_stage += river[i].ws.stage;
        river[i].daily.avg_gw += river[i].ws.gw;

        /* Lateral flux */
        for (j = 0; j < NUM_RIVFLX; j++)
        {
            river[i].daily.avg_rivflow[j] += river[i].wf.rivflow[j];
        }

        (river[i].daily.counter)++;
    }
#endif

    /* Calculate daily variables */
    if ((t - start_time) % DAYINSEC == 0 && t > start_time)
    {
#ifdef _OPENMP
# pragma omp parallel for
#endif
        for (i = 0; i < nelem; i++)
        {
            int             k;

            elem[i].daily.avg_sfctmp /= (double)elem[i].daily.counter;

            elem[i].daily.avg_sfcspd /= (double)elem[i].daily.counter;

            for (k = 0; k < elem[i].ps.nsoil; k++)
            {
                elem[i].daily.avg_stc[k] /= (double)elem[i].daily.counter;
                elem[i].daily.avg_sh2o[k] /= (double)elem[i].daily.counter;
                elem[i].daily.avg_smc[k] /= (double)elem[i].daily.counter;
                elem[i].daily.avg_smflxv[k] /= (double)elem[i].daily.counter;
#ifdef _CYCLES_
                elem[i].daily.avg_et[k] /= (double)elem[i].daily.counter;
#endif
            }

#ifdef _CYCLES_
            elem[i].daily.avg_sncovr /= (double)elem[i].daily.counter;
#endif

            elem[i].daily.avg_surf /= (double)elem[i].daily.counter;
            elem[i].daily.avg_unsat /= (double)elem[i].daily.counter;
            elem[i].daily.avg_gw /= (double)elem[i].daily.counter;

            elem[i].daily.tday /= (double)elem[i].daily.daylight_counter;
            elem[i].daily.avg_q2d /= (double)elem[i].daily.daylight_counter;
            elem[i].daily.avg_ch /= (double)elem[i].daily.daylight_counter;
            elem[i].daily.avg_rc /= (double)elem[i].daily.daylight_counter;
            elem[i].daily.avg_sfcprs /= (double)elem[i].daily.daylight_counter;
            elem[i].daily.avg_albedo /= (double)elem[i].daily.daylight_counter;
            elem[i].daily.avg_soldn /= (double)elem[i].daily.daylight_counter;

            elem[i].daily.tnight /= (double)(elem[i].daily.counter -
                elem[i].daily.daylight_counter);
        }

#ifdef _CYCLES_
# ifdef _OPENMP
#  pragma omp parallel for
# endif
        for (i = 0; i < nriver; i++)
        {
            int             j;

            river[i].daily.avg_stage /= (double)river[i].daily.counter;
            river[i].daily.avg_gw /= (double)river[i].daily.counter;

            for (j = 0; j < NUM_RIVFLX; j++)
            {
                river[i].daily.avg_rivflow[j] /=
                    (double)river[i].daily.counter;
            }
        }
#endif
    }
}

void InitDailyStruct(elem_struct *elem, river_struct *river)
{
    int             i;

#ifdef _OPENMP
# pragma omp parallel for
#endif
#ifdef _LUMPED_
    for (i = 0; i < nelem + 1; i++)
#else
    for (i = 0; i < nelem; i++)
#endif
    {
        int             k;

        elem[i].daily.counter = 0;
        elem[i].daily.daylight_counter = 0;

        elem[i].daily.avg_surf = 0.0;
        elem[i].daily.avg_unsat = 0.0;
        elem[i].daily.avg_gw = 0.0;
        for (k = 0; k < MAXLYR; k++)
        {
            elem[i].daily.avg_sh2o[k] = 0.0;
            elem[i].daily.avg_smc[k] = 0.0;
            elem[i].daily.avg_et[k] = 0.0;
            elem[i].daily.avg_smflxv[k] = 0.0;
            elem[i].daily.avg_stc[k] = 0.0;
        }

        elem[i].daily.avg_q2d = 0.0;
        elem[i].daily.avg_sfcprs = 0.0;
        elem[i].daily.avg_ch = 0.0;
        elem[i].daily.avg_rc = 0.0;
        elem[i].daily.avg_albedo = 0.0;
        elem[i].daily.avg_sfcspd = 0.0;

        elem[i].daily.tmax = -999.0;
        elem[i].daily.tmin = 999.0;
        elem[i].daily.avg_sfctmp = 0.0;
        elem[i].daily.tday = 0.0;
        elem[i].daily.tnight = 0.0;

        elem[i].daily.avg_soldn = 0.0;
        elem[i].daily.solar_total = 0.0;
    }

#ifdef _CYCLES_
# ifdef _OPENMP
#  pragma omp parallel for
# endif
    for (i = 0; i < nriver; i++)
    {
        int             k;

        river[i].daily.counter = 0;

        river[i].daily.avg_stage = 0.0;
        river[i].daily.avg_gw = 0.0;
        for (k = 0; k < NUM_RIVFLX; k++)
        {
            river[i].daily.avg_rivflow[k] = 0.0;
        }
    }
#endif
}
