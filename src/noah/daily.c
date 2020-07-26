#include "pihm.h"

void DailyVar(int t, int start_time, elem_struct *elem)
{
    int             i;

    /*
     * Sum daily variables
     */
#if defined(_OPENMP)
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

        /* Soil moisture, temperature, and ET */
        for (k = 0; k < elem[i].ps.nlayers; k++)
        {
            elem[i].daily.avg_stc[k] += elem[i].es.stc[k];
            elem[i].daily.avg_sh2o[k] += elem[i].ws.swc[k];
#if defined(_CYCLES_OBSOLETE_)
            elem[i].daily.avg_et[k] += elem[i].wf.et[k];
#endif
        }

#if defined(_CYCLES_OBSOLETE_)
        elem[i].daily.avg_sncovr += elem[i].ps.sncovr;
#endif

        if (elem[i].ef.soldn > 0.0)
        {
            elem[i].daily.tday += elem[i].es.sfctmp;
            elem[i].daily.avg_q2d += elem[i].ps.q2sat - elem[i].ps.q2;
            elem[i].daily.avg_ch += elem[i].ps.ch;
            elem[i].daily.avg_rc += elem[i].ps.rc;
            elem[i].daily.avg_sfcprs += elem[i].ps.sfcprs;
            elem[i].daily.avg_albedo += elem[i].ps.albedo;
            elem[i].daily.avg_soldn += elem[i].ef.soldn;
            (elem[i].daily.daylight_counter)++;
        }
        else
        {
            elem[i].daily.tnight += elem[i].es.sfctmp;
        }

        (elem[i].daily.counter)++;
    }

    /* Calculate daily variables */
    if ((t - start_time) % DAYINSEC == 0 && t > start_time)
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (i = 0; i < nelem; i++)
        {
            int             k;

            elem[i].daily.avg_sfctmp /= (double)elem[i].daily.counter;

            for (k = 0; k < elem[i].ps.nlayers; k++)
            {
                elem[i].daily.avg_stc[k] /= (double)elem[i].daily.counter;
                elem[i].daily.avg_sh2o[k] /= (double)elem[i].daily.counter;
#if defined(_CYCLES_OBSOLETE_)
                elem[i].daily.avg_et[k] /= (double)elem[i].daily.counter;
#endif
            }

#if defined(_CYCLES_OBSOLETE_)
            elem[i].daily.avg_sncovr /= (double)elem[i].daily.counter;
#endif

            elem[i].daily.tday /= (double)elem[i].daily.daylight_counter;
            elem[i].daily.avg_q2d /= (double)elem[i].daily.daylight_counter;
            elem[i].daily.avg_ch /= (double)elem[i].daily.daylight_counter;
            elem[i].daily.avg_rc /= (double)elem[i].daily.daylight_counter;
            elem[i].daily.avg_sfcprs /= (double)elem[i].daily.daylight_counter;
            elem[i].daily.avg_albedo /= (double)elem[i].daily.daylight_counter;
#if defined(_CYCLES_OBSOLETE_)
            elem[i].daily.avg_soldn /= (double)elem[i].daily.counter;
#else
            elem[i].daily.avg_soldn /= (double)elem[i].daily.daylight_counter;
#endif

            elem[i].daily.tnight /= (double)(elem[i].daily.counter -
                elem[i].daily.daylight_counter);
        }

#if defined(_LUMPEDBGC_)
        int             k;

        elem[LUMPEDBGC].daily.tmax = 0.0;
        elem[LUMPEDBGC].daily.tmin = 0.0;
        for (i = 0; i < nelem; i++)
        {
            elem[LUMPEDBGC].daily.tmax += elem[i].daily.tmax * elem[i].topo.area;
            elem[LUMPEDBGC].daily.tmin += elem[i].daily.tmin * elem[i].topo.area;
            elem[LUMPEDBGC].daily.avg_sfctmp +=
                elem[i].daily.avg_sfctmp * elem[i].topo.area;
            /* When running lumped model, only average root zone to avoid uneven
             * layers */
            for (k = 0; k < elem[i].ps.nlayers; k++)
            {
                elem[LUMPEDBGC].daily.avg_stc[k] +=
                    elem[i].daily.avg_stc[k] * elem[i].topo.area;
                elem[LUMPEDBGC].daily.avg_sh2o[k] +=
                    elem[i].daily.avg_sh2o[k] * elem[i].topo.area;
                elem[LUMPEDBGC].daily.avg_smc[k] +=
                    elem[i].daily.avg_smc[k] * elem[i].topo.area;
            }

            elem[LUMPEDBGC].daily.tday +=
                elem[i].daily.tday * elem[i].topo.area;
            elem[LUMPEDBGC].daily.avg_q2d +=
                elem[i].daily.avg_q2d * elem[i].topo.area;
            elem[LUMPEDBGC].daily.avg_ch +=
                elem[i].daily.avg_ch * elem[i].topo.area;
            elem[LUMPEDBGC].daily.avg_rc +=
                elem[i].daily.avg_rc * elem[i].topo.area;
            elem[LUMPEDBGC].daily.avg_sfcprs +=
                elem[i].daily.avg_sfcprs * elem[i].topo.area;
            elem[LUMPEDBGC].daily.avg_albedo +=
                elem[i].daily.avg_albedo * elem[i].topo.area;
            elem[LUMPEDBGC].daily.avg_soldn +=
                elem[i].daily.avg_soldn * elem[i].topo.area;
            elem[LUMPEDBGC].daily.tnight +=
                elem[i].daily.tnight * elem[i].topo.area;
        }

        elem[LUMPEDBGC].daily.tmax /= elem[LUMPEDBGC].topo.area;
        elem[LUMPEDBGC].daily.tmin /= elem[LUMPEDBGC].topo.area;
        elem[LUMPEDBGC].daily.avg_sfctmp /= elem[LUMPEDBGC].topo.area;
        /* When running lumped model, only average root zone to avoid uneven
         * layers */
        for (k = 0; k < elem[LUMPEDBGC].ps.nlayers; k++)
        {
            elem[LUMPEDBGC].daily.avg_stc[k] /= elem[LUMPEDBGC].topo.area;
            elem[LUMPEDBGC].daily.avg_sh2o[k] /= elem[LUMPEDBGC].topo.area;
            elem[LUMPEDBGC].daily.avg_smc[k] /= elem[LUMPEDBGC].topo.area;
        }

        elem[LUMPEDBGC].daily.tday /= elem[LUMPEDBGC].topo.area;
        elem[LUMPEDBGC].daily.avg_q2d /= elem[LUMPEDBGC].topo.area;
        elem[LUMPEDBGC].daily.avg_ch /= elem[LUMPEDBGC].topo.area;
        elem[LUMPEDBGC].daily.avg_rc /= elem[LUMPEDBGC].topo.area;
        elem[LUMPEDBGC].daily.avg_sfcprs /= elem[LUMPEDBGC].topo.area;
        elem[LUMPEDBGC].daily.avg_albedo /= elem[LUMPEDBGC].topo.area;
        elem[LUMPEDBGC].daily.avg_soldn /= elem[LUMPEDBGC].topo.area;
        elem[LUMPEDBGC].daily.tnight /= elem[LUMPEDBGC].topo.area;
#endif
    }
}

void InitDailyStruct(elem_struct *elem)
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
#if defined(_LUMPEDBGC_)
    for (i = 0; i < nelem + 1; i++)
#else
    for (i = 0; i < nelem; i++)
#endif
    {
        int             k;

        elem[i].daily.counter = 0;
        elem[i].daily.daylight_counter = 0;

        for (k = 0; k < MAXLYR; k++)
        {
            elem[i].daily.avg_sh2o[k] = 0.0;
            elem[i].daily.avg_stc[k] = 0.0;
#if defined(_CYCLES_OBSOLETE_)
            elem[i].daily.avg_et[k] = 0.0;
#endif
        }

        elem[i].daily.avg_q2d = 0.0;
        elem[i].daily.avg_sfcprs = 0.0;
        elem[i].daily.avg_ch = 0.0;
        elem[i].daily.avg_rc = 0.0;
        elem[i].daily.avg_albedo = 0.0;
#if defined(_CYCLES_OBSOLETE_)
        elem[i].daily.avg_sncovr = 0.0;
#endif

        elem[i].daily.tmax = -999.0;
        elem[i].daily.tmin = 999.0;
        elem[i].daily.avg_sfctmp = 0.0;
        elem[i].daily.tday = 0.0;
        elem[i].daily.tnight = 0.0;

        elem[i].daily.avg_soldn = 0.0;

#if defined(_CYCLES_OBSOLETE_)
        for (k = 0; k < MAXCROP && '\0' != elem[i].crop[k].epc->cropn[0]; k++)
        {
            elem[i].crop[k].cwf.transp = 0.0;
            elem[i].crop[k].cwf.transp_pot = 0.0;
        }
#endif
    }
}
