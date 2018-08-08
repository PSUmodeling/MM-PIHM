#include "pihm.h"

void DailyCycles(int t, pihm_struct pihm)
{
    int             i;
    int             year, month, day;
    int             doy;
    pihm_t_struct   pihm_t;

    pihm_t = PIHMTime(t - DAYINSEC);
    year = pihm_t.year;
    month = pihm_t.month;
    day = pihm_t.day;
    doy = Doy(year, month, day);

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             ind;
        elem_struct    *elem;

        elem = &pihm->elem[i];
        ind = elem->attrib.op_type - 1;

        DailyOperations(year, doy, &pihm->opertbl[ind], &elem->daily,
            &elem->soil, &elem->mgmt, elem->crop, &elem->ps, &elem->ws,
            &elem->wf, &elem->cs, &elem->cf, &elem->ns, &elem->nf);
    }
}
