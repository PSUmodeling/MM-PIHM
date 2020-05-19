#include "pihm.h"

void DailyCycles(int t, pihm_struct pihm)
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
        elem_struct    *elem;

        elem = &pihm->elem[i];
        ind = elem->attrib.op_type - 1;

        /*
         * Run daily cycles processes
         */
        DailyOperations(doy, &pihm->opertbl[ind], &elem->daily,
            &elem->soil, &elem->mgmt, elem->crop, &elem->ps, &elem->ws,
            &elem->wf, &elem->cs, &elem->cf, &elem->ns, &elem->nf);

        /* Calculate daily sink/source terms for NO3 and NH4 */
        CalSnkSrc(&elem->nf, elem->ps.nlayers, &elem->no3sol, &elem->nh4sol);
    }
}

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
