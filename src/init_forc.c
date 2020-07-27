#include "pihm.h"

#if defined(_RT_)
void InitForcing(const rttbl_struct *rttbl, const calib_struct *cal,
    forc_struct *forc, elem_struct elem[])
#else
void InitForcing(const calib_struct *cal, forc_struct *forc, elem_struct elem[])
#endif
{
    int             i, j;

    /* Apply climate scenarios */
    for (i = 0; i < forc->nmeteo; i++)
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (j = 0; j < forc->meteo[i].length; j++)
        {
            forc->meteo[i].data[j][PRCP_TS] *= cal->prcp;
            forc->meteo[i].data[j][SFCTMP_TS] += cal->sfctmp;
        }
    }

    if (forc->nbc > 0)
    {
        for (i = 0; i < forc->nbc; i++)
        {
#if defined(_RT_)
            forc->bc[i].value =
                (double *)malloc((rttbl->num_stc + 1) * sizeof(double));
#else
            forc->bc[i].value = (double *)malloc(sizeof(double));
#endif
        }
    }
    if (forc->nmeteo > 0)
    {
        for (i = 0; i < forc->nmeteo; i++)
        {
            forc->meteo[i].value =
                (double *)malloc(NUM_METEO_VAR * sizeof(double));
        }
    }
    if (forc->nlai > 0)
    {
        for (i = 0; i < forc->nlai; i++)
        {
            forc->lai[i].value = (double *)malloc(sizeof(double));
        }
    }
    if (forc->nriverbc > 0)
    {
        for (i = 0; i < forc->nriverbc; i++)
        {
            forc->riverbc[i].value = (double *)malloc(sizeof(double));
        }
    }
    if (forc->nsource > 0)
    {
        for (i = 0; i < forc->nsource; i++)
        {
            forc->source[i].value = (double *)malloc(sizeof(double));
        }
    }
#if defined(_NOAH_)
    if (forc->nrad > 0)
    {
        for (i = 0; i < forc->nrad; i++)
        {
            forc->rad[i].value = (double *)malloc(2 * sizeof(double));
        }
    }
#endif

#if defined(_BGC_)
    if (forc->nco2 > 0)
    {
        forc->co2[0].value = (double *)malloc(sizeof(double));
    }

    if (forc->nndep > 0)
    {
        forc->ndep[0].value = (double *)malloc(sizeof(double));
    }
#endif

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem[i].ps.zlvl_wind =
            forc->meteo[elem[i].attrib.meteo_type - 1].zlvl_wind;
    }
}
