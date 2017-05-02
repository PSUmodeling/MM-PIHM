#include "pihm.h"

void DailyVar (int t, int start_time, pihm_struct pihm)
{
    int             i;

    /*
     * Cumulates daily variables
     */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int         k;
        elem_struct *elem;

        elem = &pihm->elem[i];

        /* Air temperature */
        elem->daily.avg_sfctmp += elem->es.sfctmp;
        elem->daily.tmax = (elem->daily.tmax > elem->es.sfctmp) ?
            elem->daily.tmax : elem->es.sfctmp;
        elem->daily.tmin = (elem->daily.tmin < elem->es.sfctmp) ?
            elem->daily.tmin : elem->es.sfctmp;

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

        if (elem->ef.soldn > 0.0)
        {
            elem->daily.tday += elem->es.sfctmp;
            elem->daily.avg_q2d += elem->ps.q2sat - elem->ps.q2;
            elem->daily.avg_ra += 1.0 / elem->ps.ch;
#ifdef _DEBUG_
            if (i == 0)
            {
                printf ("soldn = %lf, ra = %lf, q2sat = %lf, q2 = %lf\n", elem->ef.soldn, elem->ps.rc, elem->ps.q2sat, elem->ps.q2);
            }
#endif
            elem->daily.avg_rc += elem->ps.rc;
            elem->daily.avg_sfcprs += elem->ps.sfcprs;
            elem->daily.avg_albedo += elem->ps.albedo;
            elem->daily.avg_soldn += elem->ef.soldn;
            elem->daily.solar_total += elem->ef.soldn * pihm->ctrl.stepsize;
            (elem->daily.daylight_counter)++;
        }
        else
        {
            elem->daily.tnight += elem->es.sfctmp;
        }

        (elem->daily.counter)++;
    }

#ifdef _CYCLES_
    /* River segments */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
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
#endif

    /* Calculate daily variables */
    if ((t - start_time) % DAYINSEC == 0 && t > start_time)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < nelem; i++)
        {
            int     k;
            elem_struct *elem;

            elem = &pihm->elem[i];

            elem->daily.avg_sfctmp /= (double)elem->daily.counter;

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

            elem->daily.tday /= (double)elem->daily.daylight_counter;
            elem->daily.avg_q2d /= (double)elem->daily.daylight_counter;
            elem->daily.avg_ra /= (double)elem->daily.daylight_counter;
            elem->daily.avg_rc /= (double)elem->daily.daylight_counter;
            elem->daily.avg_sfcprs /= (double)elem->daily.daylight_counter;
            elem->daily.avg_albedo /= (double)elem->daily.daylight_counter;
            elem->daily.avg_soldn /= (double)elem->daily.daylight_counter;

            elem->daily.tnight /= (double)(elem->daily.counter -
                elem->daily.daylight_counter);
        }

#ifdef _CYCLES_
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < nriver; i++)
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
#endif
    }
}

void InitDailyStruct (pihm_struct pihm)
{
    int             i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
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

        elem->daily.avg_q2d = 0.0;
        elem->daily.avg_sfcprs = 0.0;
        elem->daily.avg_ra = 0.0;
        elem->daily.avg_rc = 0.0;
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

#ifdef _CYCLES_
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
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
#endif
}
