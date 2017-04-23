#include "pihm.h"

void PIHM (pihm_struct pihm, void *cvode_mem, N_Vector CV_Y, int t, int next_t)
{
    /* Apply forcing */
    ApplyForcing (&pihm->forc, pihm->elem, pihm->riv, t
#ifdef _NOAH_
        , &pihm->ctrl, pihm->latitude, pihm->longitude, pihm->elevation,
        pihm->noahtbl.tbot
#endif
        );

    /* Determine if land surface simulation is needed */
    if ((t - pihm->ctrl.starttime) % pihm->ctrl.etstep == 0)
    {
#ifdef _NOAH_
        /* Calculate surface energy balance */
        Noah (pihm);
#else
        /* Calculate Interception storage and ET */
        IntcpSnowET (t, (double)pihm->ctrl.etstep, pihm);
#endif
    }

    /*
     * Solve PIHM hydrology ODE using CVODE
     */
    SolveCVode (&t, next_t, pihm->ctrl.stepsize,
            cvode_mem, CV_Y);

    /* Use mass balance to calculate model fluxes or variables */
    Summary (pihm, CV_Y, (double)pihm->ctrl.stepsize);

#ifdef _NOAH_
    NoahHydrol (pihm->elem, (double)pihm->ctrl.stepsize);
#endif

    /*
     * Daily timestep modules
     */
#ifdef _DAILY_
    DailyVar (t, pihm->ctrl.starttime, pihm);

    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
#ifdef _CYCLES_
        DailyCycles (t - DAYINSEC, pihm);
#endif
#ifdef _BGC_
        if (pihm->ctrl.bgc_spinup)
        {
            Save2Stor (pihm, t, pihm->ctrl.spinupstart);
        }
        else
        {
            DailyBgc (pihm, t, pihm->ctrl.starttime);
            first_day = 0;
        }
#endif
    }
#endif

#ifdef _CYCLES_
    SoluteTransport (pihm->elem, pihm->riv, (double)pihm->ctrl.stepsize);
#endif

#ifdef _DAILY_
    /*
     * Initialize daily structures
     */
    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
        InitDailyStruct (pihm);
    }
#endif
}
