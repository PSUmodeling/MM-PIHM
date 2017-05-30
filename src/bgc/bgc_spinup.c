#include "pihm.h"

void BgcSpinup (pihm_struct pihm, N_Vector CV_Y, void *cvode_mem)
{
    int             i;
    int             spinyears = 0;
    int             first_spin_cycle = 1;
    int             ss_total;
    double          metyears;

    metyears =
        (pihm->ctrl.endtime - pihm->ctrl.starttime) / DAYINSEC / 365;

    do
    {
        PIHMprintf (VL_NORMAL, "Spinup year: %6d\n", spinyears + 1);

        ResetSpinupStat (pihm->elem);

        for (i = 0; i < pihm->ctrl.nstep; i++)
        {
            PIHM (pihm, cvode_mem, CV_Y, pihm->ctrl.tout[i],
                pihm->ctrl.tout[i + 1]);
        }

        /* Reset solver parameters */
        SetCVodeParam (pihm, cvode_mem, CV_Y);

        spinyears += metyears;

        ss_total = CheckBgcSS (pihm->elem, first_spin_cycle,
            pihm->ctrl.endtime - pihm->ctrl.starttime, spinyears);

        first_spin_cycle = 0;
    } while (spinyears < pihm->ctrl.maxspinyears && ss_total < nelem);
}

void ResetSpinupStat (elem_struct *elem)
{
    int             i;

    for (i = 0; i < nelem; i++)
    {
        elem[i].spinup.soilc = 0.0;
        elem[i].spinup.totalc = 0.0;
    }
}

int CheckBgcSS (elem_struct *elem, int first_cycle, int totalt,
    int spinyears)
{
    int             i;
    double          t1;
    int             total_complete = 0;

    for (i = 0; i < nelem; i++)
    {
        elem[i].spinup.soilc /= (double)(totalt / DAYINSEC);
        elem[i].spinup.totalc /= (double)(totalt / DAYINSEC);

        if (!first_cycle)
        {
            /* Convert soilc and totalc to average daily soilc */
            t1 = (elem[i].spinup.soilc - elem[i].spinup.soilc_prev) /
                (double)(totalt / DAYINSEC / 365);

            /* Check if element reaches steady state */
            elem[i].spinup.steady = (fabs (t1) < SPINUP_TOLERANCE);

            PIHMprintf (VL_NORMAL, "Elem %d: "
                "spinyears = %d soilc_prev = %lg soilc = %lg pdif = %lg\n",
                i + 1, spinyears, elem[i].spinup.soilc_prev,
                elem[i].spinup.soilc, t1);
        }
        else
        {
            elem[i].spinup.steady = 0;
        }

        elem[i].spinup.soilc_prev = elem[i].spinup.soilc;

        total_complete += elem[i].spinup.steady;
    }

    PIHMprintf (VL_NORMAL, "%d elements steady, %d elements to go.\n",
            total_complete, nelem - total_complete);

    if (total_complete == nelem)
    {
        PIHMprintf (VL_NORMAL, "All elements steady after %d year.\n",
                spinyears);
    }

    return (total_complete);
}
