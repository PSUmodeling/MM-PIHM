#include "pihm.h"

void BgcSpinup(pihm_struct pihm, N_Vector CV_Y, void *cvode_mem)
{
    int             i;
    int             spinyears = 0;
    int             first_spin_cycle = 1;
    int             steady;
    int             metyears;

    metyears = (pihm->ctrl.endtime - pihm->ctrl.starttime) / DAYINSEC / 365;

    do
    {
        PIHMprintf(VL_NORMAL, "Spinup year: %6d\n", spinyears + 1);

        ResetSpinupStat(pihm->elem);

        for (i = 0; i < pihm->ctrl.nstep; i++)
        {
            PIHM(pihm, cvode_mem, CV_Y, pihm->ctrl.tout[i],
                pihm->ctrl.tout[i + 1], 0.0);
        }

        /* Reset solver parameters */
        SetCVodeParam(pihm, cvode_mem, CV_Y);

        spinyears += metyears;

        steady = CheckBgcSteadyState(pihm->elem, pihm->siteinfo.area,
            first_spin_cycle, pihm->ctrl.endtime - pihm->ctrl.starttime,
            spinyears);

        first_spin_cycle = 0;
    } while (spinyears < pihm->ctrl.maxspinyears);
}

void ResetSpinupStat(elem_struct *elem)
{
    int             i;

#if defined(_LUMPED_)
    i = LUMPED;
#else
# ifdef _OPENMP
#  pragma omp parallel for
# endif
    for (i = 0; i < nelem; i++)
#endif
    {
        elem[i].spinup.soilc = 0.0;
        elem[i].spinup.totalc = 0.0;
    }
}

int CheckBgcSteadyState(const elem_struct *elem, double total_area,
    int first_cycle, int totalt, int spinyears)
{
    int             i;
    double          t1;
    double          soilc = 0.0;
    double          totalc = 0.0;
    int             steady;
    static double   soilc_prev = 0.0;

    /* Convert soilc and totalc to average daily soilc */
#if defined(_LUMPED_)
    i = LUMPED;
#else
# ifdef _OPENMP
#  pragma omp parallel for
# endif
    for (i = 0; i < nelem; i++)
#endif
    {
        soilc += elem[i].spinup.soilc * elem[i].topo.area /
            (double)(totalt / DAYINSEC) / total_area;
        totalc += elem[i].spinup.totalc * elem[i].topo.area /
            (double)(totalt / DAYINSEC) / total_area;
    }

    if (!first_cycle)
    {
        t1 = (soilc - soilc_prev) / (double)(totalt / DAYINSEC / 365);

        /* Check if domain reaches steady state */
        steady = (fabs(t1) < SPINUP_TOLERANCE);

        PIHMprintf(VL_NORMAL,
            "spinyears = %d soilc_prev = %lg soilc = %lg pdif = %lg\n",
            spinyears, soilc_prev, soilc, t1);
    }
    else
    {
        steady = 0;
    }

    soilc_prev = soilc;

    if (steady)
    {
        PIHMprintf(VL_NORMAL, "Reaches steady state after %d year.\n",
            spinyears);
    }

    return steady;
}
