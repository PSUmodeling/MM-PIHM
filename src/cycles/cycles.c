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
    PIHMprintf(VL_BRIEF, "DOY %d\n", doy);
#endif

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             ind;

        /*
         * Run daily cycles processes
         */
        DailyOper(year, doy, elem[i].mgmt.auto_n, &elem[i].weather,
            &elem[i].mgmt, elem[i].crop, &elem[i].soil, &elem[i].ws,
            &elem[i].wf, &elem[i].es, &elem[i].cs, &elem[i].cf, &elem[i].ns,
            &elem[i].nf, &elem[i].ps);

#if TEMP_DISABLED
        /* Calculate daily sink/source terms for NO3 and NH4 */
        CalSnkSrc(&elem[i].nf, elem[i].ps.nlayers, &elem[i].no3sol, &elem[i].nh4sol);
#endif
    }
}

#if TEMP_DISABLED
void CalSnkSrc(const nflux_struct *nf, int nlayers, solute_struct *no3sol,
    solute_struct *nh4sol)
{
    int             k;

    no3sol->snksrc = 0.0;
    nh4sol->snksrc = 0.0;

    no3sol->snksrc += nf->surplusn;
    nh4sol->snksrc += nf->urine;
    for (k = 0; k < nlayers; k++)
    {
        no3sol->snksrc += nf->uptake_no3[k] + nf->fert_no3[k] +
            nf->immob_no3[k] + nf->nitrif_nh4_to_no3[k] + nf->denitn[k] +
            nf->till_no3[k];

        nh4sol->snksrc += nf->uptake_nh4[k] + nf->fert_nh4[k] +
            nf->till_nh4[k] + nf->immob_nh4[k] - nf->nitrif_nh4_to_no3[k] -
            nf->nitrif_nh4_to_n2o[k] - nf->nh4volat[k];
    }
}
#endif