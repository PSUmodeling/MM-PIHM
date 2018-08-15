#include "pihm.h"

void PIHM(pihm_struct pihm, void *cvode_mem, N_Vector CV_Y, int t,
    int next_t, double cputime)
{
    /* Apply boundary conditions */
    ApplyBc(t, &pihm->forc, pihm->elem, pihm->river);

    /*
     * Apply forcing and simulate land surface processes
     */
    if ((t - pihm->ctrl.starttime) % pihm->ctrl.etstep == 0)
    {
        /* Apply forcing */
#if defined(_NOAH_)
        ApplyForc(t, pihm->ctrl.rad_mode, &pihm->siteinfo, &pihm->forc,
            pihm->elem);
#else
        ApplyForc(t, &pihm->forc, pihm->elem);
#endif

#if defined(_NOAH_)
        /* Calculate surface energy balance */
        Noah(pihm->elem, (double)pihm->ctrl.etstep);
#else
        /* Calculate Interception storage and ET */
        IntcpSnowEt(t, (double)pihm->ctrl.etstep, &pihm->cal, pihm->elem);
#endif

        /* Update print variables for land surface step variables */
        UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, LS_STEP);
        UpdPrintVar(pihm->print.tp_varctrl, pihm->print.ntpprint, LS_STEP);
    }

    /*
     * Solve PIHM hydrology ODE using CVode
     */
    SolveCVode(pihm->ctrl.starttime, &t, next_t, cputime, cvode_mem, CV_Y);

    /* Use mass balance to calculate model fluxes or variables */
    Summary(pihm->elem, pihm->river, CV_Y, (double)pihm->ctrl.stepsize);

#if defined(_NOAH_)
    NoahHydrol(pihm->elem, (double)pihm->ctrl.stepsize);
#endif

#if defined(_CYCLES_)
//    SoluteTransport(pihm->elem, pihm->river, (double)pihm->ctrl.stepsize);
#endif

    /* Update print variables for hydrology step variables */
    UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, HYDROL_STEP);
    UpdPrintVar(pihm->print.tp_varctrl, pihm->print.ntpprint, HYDROL_STEP);

#if defined(_DAILY_)
    DailyVar(t, pihm->ctrl.starttime, pihm->elem);

    /*
     * Daily timestep modules
     */
    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
# if defined(_BGC_)
        /* Daily BGC processes */
        DailyBgc(pihm, t - DAYINSEC);
# endif

# if defined(_CYCLES_)
        DailyCycles(t - DAYINSEC, pihm);
# endif

        /* Update print variables for CN (daily) step variables */
        UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, CN_STEP);
        UpdPrintVar(pihm->print.tp_varctrl, pihm->print.ntpprint, CN_STEP);

        /* Initialize daily structures */
        InitDailyStruct(pihm->elem);
    }
#endif

    /*
     * Print outputs
     */
    /* Print water balance */
    if (pihm->ctrl.waterbal)
    {
        PrintWaterBal(pihm->print.watbal_file, t, pihm->ctrl.starttime,
            pihm->ctrl.stepsize, pihm->elem, pihm->river);
    }

    /* Print binary and txt output files */
    PrintData(pihm->print.varctrl, pihm->print.nprint, t,
        t - pihm->ctrl.starttime, pihm->ctrl.ascii);

    /* Print tecplot output files */
    if (tecplot)
    {
        PrintDataTecplot(pihm->print.tp_varctrl, pihm->print.ntpprint, t,
            t - pihm->ctrl.starttime);
    }
}
