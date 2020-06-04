#include "pihm.h"

void Cycles(int t, elem_struct elem[])
{
    int             i;
    int             year, month, day;
    int             doy;
    pihm_t_struct   pihm_t;

    pihm_t = PIHMTime(t);
    year = pihm_t.year;
    month = pihm_t.month;
    day = pihm_t.day;
    doy = Doy(year, month, day);

#if defined(_DEBUG_)
    pihm_printf(VL_BRIEF, "DOY %d\n", doy);
#endif

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        /*
         * Run daily cycles processes
         */
        DailyOper(year, doy, elem[i].mgmt.auto_n, &elem[i].weather,
            &elem[i].mgmt, elem[i].crop, &elem[i].soil, &elem[i].ws,
            &elem[i].wf, &elem[i].es, &elem[i].cs, &elem[i].cf, &elem[i].ns,
            &elem[i].nf, &elem[i].ps);

        /* Calculate daily sink/source terms for NO3 and NH4 */
        CalSnkSrc(elem[i].ps.nlayers, &elem[i].nf, elem[i].solute);
    }
}

void CalSnkSrc(int nlayers, const nflux_struct *nf, solute_struct solute[])
{
    int             k;

    for (k = 0; k < MAXLYR; k++)
    {
        solute[NO3].snksrc[k] = 0.0;
        solute[NH4].snksrc[k] = 0.0;
    }

    solute[NO3].snksrc[0] += nf->surplus;
    solute[NH4].snksrc[0] += nf->urine;
    for (k = 0; k < nlayers; k++)
    {
        solute[NO3].snksrc[k] += (nf->nitrif[k] - nf->n2o_from_nitrif[k]) +
            (-nf->denitrif[k]) + (-nf->no3_uptake[k]) + nf->no3_fert[k] +
            nf->no3_immobil[k];

        solute[NH4].snksrc[k] += (-nf->nitrif[k]) + (-nf->volatil[k]) +
            (-nf->nh4_uptake[k]) + nf->nh4_fert[k] +
            nf->nh4_immobil[k] + nf->mineral[k];

        solute[NO3].snksrc[k] /= DAYINSEC;
        solute[NH4].snksrc[k] /= DAYINSEC;
    }
}
