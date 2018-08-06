#include "pihm.h"

void DailyCycles(int t, pihm_struct pihm)
{
    int             i;
    int             year, month, day;
    int             doy;
    elem_struct    *elem;
    pihm_t_struct   pihm_t;

    elem = pihm->elem;

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

        ind = elem[i].attrib.op_type - 1;
    }
}
