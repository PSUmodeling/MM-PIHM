#include "pihm.h"

void Cycles(int t, const co2control_struct *co2ctrl, forc_struct *forc, elem_struct elem[])
{
    int             i;
    int             year, month, day;
    int             doy;
    double          co2;
    pihm_t_struct   pihm_t;

    pihm_t = PIHMTime(t);
    year = pihm_t.year;
    month = pihm_t.month;
    day = pihm_t.day;
    doy = Doy(year, month, day);

#if defined(_DEBUG_)
    pihm_printf(VL_BRIEF, "DOY %d\n", doy);
#endif

    co2 = (co2ctrl->varco2 == 1) ? GetCO2(t, &forc->co2[0]) : co2ctrl->co2ppm;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem[i].weather.co2 = co2;

        // Run daily cycles processes
        DailyOper(year, doy, elem[i].mgmt.auto_n, &elem[i].weather, &elem[i].mgmt, elem[i].crop, &elem[i].soil,
            &elem[i].ws, &elem[i].wf, &elem[i].es, &elem[i].cs, &elem[i].cf, &elem[i].ns, &elem[i].nf, &elem[i].ps);

        // Calculate daily sink/source terms for NO3 and NH4
        CalSnkSrc(elem[i].ps.nlayers, &elem[i].nf, elem[i].solute);
    }
}

void CalSnkSrc(int nlayers, const nflux_struct *nf, solute_struct solute[])
{
    int             kz;

    for (kz = 0; kz < MAXLYR; kz++)
    {
        solute[NO3].snksrc[kz] = 0.0;
        solute[NH4].snksrc[kz] = 0.0;
    }

    solute[NO3].snksrc[0] += nf->surplus;
    solute[NH4].snksrc[0] += nf->urine;
    for (kz = 0; kz < nlayers; kz++)
    {
        solute[NO3].snksrc[kz] += (nf->nitrif[kz] - nf->n2o_from_nitrif[kz]) + (-nf->denitrif[kz]) +
            (-nf->no3_uptake[kz]) + nf->no3_fert[kz] + nf->no3_immobil[kz];

        solute[NH4].snksrc[kz] += (-nf->nitrif[kz]) + (-nf->volatil[kz]) + (-nf->nh4_uptake[kz]) + nf->nh4_fert[kz] +
            nf->nh4_immobil[kz] + nf->mineral[kz];

        solute[NO3].snksrc[kz] /= DAYINSEC;
        solute[NH4].snksrc[kz] /= DAYINSEC;
    }
}
