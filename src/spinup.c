#include "pihm.h"

void Spinup(pihm_struct pihm, N_Vector CV_Y, void *cvode_mem)
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

#if defined(_BGC_)
        ResetSpinupStat(pihm->elem);
#endif

        for (i = 0; i < pihm->ctrl.nstep; i++)
        {
            PIHM(pihm, cvode_mem, CV_Y, pihm->ctrl.tout[i],
                pihm->ctrl.tout[i + 1], 0.0);
        }

        /* Reset solver parameters */
        SetCVodeParam(pihm, cvode_mem, CV_Y);

        spinyears += metyears;

        steady = CheckSteadyState(pihm->elem, pihm->siteinfo.area,
            first_spin_cycle, pihm->ctrl.endtime - pihm->ctrl.starttime,
            spinyears);

        first_spin_cycle = 0;
    } while (spinyears < pihm->ctrl.maxspinyears && (!steady));
}

#if defined(_BGC_)
void ResetSpinupStat(elem_struct *elem)
{
    int             i;

#if defined(_LUMPED_)
    i = LUMPED;
#else
# if defined(_OPENMP)
#  pragma omp parallel for
# endif
    for (i = 0; i < nelem; i++)
#endif
    {
        elem[i].spinup.soilc = 0.0;
        elem[i].spinup.totalc = 0.0;
    }
}
#endif

int CheckSteadyState(const elem_struct *elem, double total_area,
    int first_cycle, int totalt, int spinyears)
{
    int             i;
    double          totalw = 0.0;
    static double   totalw_prev = 0.0;
    int             steady;
#if defined(_BGC_)
    double          t1;
    double          soilc = 0.0;
    double          totalc = 0.0;
    static double   soilc_prev = 0.0;
#endif

#if defined(_LUMPED_)
    i = LUMPED;
#else
    for (i = 0; i < nelem; i++)
#endif
    {
        totalw += (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity *
            elem[i].topo.area / total_area;

#if defined(_BGC_)
        /* Convert soilc and totalc to average daily soilc */
        soilc += elem[i].spinup.soilc * elem[i].topo.area /
            (double)(totalt / DAYINSEC) / total_area;
        totalc += elem[i].spinup.totalc * elem[i].topo.area /
            (double)(totalt / DAYINSEC) / total_area;
#endif
    }

    if (!first_cycle)
    {
#if defined(_BGC_)
        t1 = (soilc - soilc_prev) / (double)(totalt / DAYINSEC / 365);
#endif

        /* Check if domain reaches steady state */
#if defined(_BGC_)
        steady = ((fabs(t1) < SPINUP_C_TOLERANCE) &&
            (totalw - totalw_prev < SPINUP_W_TOLERANCE));
#else
        steady = (totalw - totalw_prev < SPINUP_W_TOLERANCE);
#endif

        PIHMprintf(VL_NORMAL, "spinyears = %d ", spinyears);
        PIHMprintf(VL_NORMAL, "totalw_prev = %lg totalw = %lg wdif = %lg\n",
            totalw_prev, totalw, totalw - totalw_prev);
#if defined(_BGC_)
        PIHMprintf(VL_NORMAL, "soilc_prev = %lg soilc = %lg pdif = %lg\n",
            soilc_prev, soilc, t1);
#endif
    }
    else
    {
        steady = 0;
    }

    totalw_prev = totalw;
#if defined(_BGC_)
    soilc_prev = soilc;
#endif

    if (steady)
    {
        PIHMprintf(VL_NORMAL, "Reaches steady state after %d year.\n",
            spinyears);
    }

    return steady;
}
