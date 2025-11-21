#include "pihm.h"

void PIHM(double cputime, pihm_struct *pihm, cvode_struct *cvode)
{
    int             t;

    t = pihm->ctrl.tout[pihm->ctrl.cstep];

    // Apply boundary conditions
    ApplyBc(t, &pihm->forc, pihm->elem, pihm->river);

    // Apply forcing and simulate land surface processes
    if ((t - pihm->ctrl.starttime) % pihm->ctrl.etstep == 0)
    {
        // Apply forcing
#if defined(_NOAH_)
        ApplyForcing(t, pihm->ctrl.rad_mode, &pihm->siteinfo, &pihm->forc, pihm->elem);
#else
        ApplyForcing(t, &pihm->forc, pihm->elem);
#endif

#if defined(_NOAH_)
        // Calculate surface energy balance
        Noah((double)pihm->ctrl.etstep, &pihm->lctbl, &pihm->calib, pihm->elem);
#else
        // Calculate Interception storage and ET
        IntcpSnowEt(t, (double)pihm->ctrl.etstep, &pihm->calib, pihm->elem);
#endif

        // Update print variables for land surface step variables
        UpdatePrintVar(pihm->print.nprint, LS_STEP, pihm->print.varctrl);
    }

#if defined(_CYCLES_)
    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
        Cycles(t, &pihm->co2ctrl, &pihm->forc, pihm->elem);

        // Update print variables for CN (daily) step variables
        UpdatePrintVar(pihm->print.nprint, CN_STEP, pihm->print.varctrl);
    }
#endif

    // Solve PIHM hydrology ODE using CVode
    SolveCVode(cputime, &pihm->ctrl, &t, cvode);

    // Use mass balance to calculate model fluxes or variables
    UpdateVar((double)pihm->ctrl.stepsize, pihm->elem, pihm->river, cvode);

#if defined(_NOAH_)
    NoahHydrol((double)pihm->ctrl.stepsize, pihm->elem);
#endif

    // Update print variables for hydrology step variables
    UpdatePrintVar(pihm->print.nprint, HYDROL_STEP, pihm->print.varctrl);

#if defined(_DAILY_)
    DailyVar(t, pihm->ctrl.starttime, pihm->elem);

    // Daily timestep modules
    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
# if defined(_BGC_)
        // Daily BGC processes
        DailyBgc(t - DAYINSEC, pihm);

        // Update print variables for CN (daily) step variables
        UpdatePrintVar(pihm->print.nprint, CN_STEP, pihm->print.varctrl);
# endif

        // Initialize daily structures
        InitDailyStruct(pihm->elem);
    }
#endif

    // Print outputs
    // Print water balance
    if (pihm->ctrl.waterbal)
    {
        PrintWaterBalance(t, pihm->ctrl.starttime, pihm->ctrl.stepsize, pihm->elem, pihm->river, pihm->print.watbal_file);
    }

    // Print binary and txt output files
    PrintData(pihm->print.nprint, t, t - pihm->ctrl.starttime, pihm->ctrl.ascii, pihm->print.varctrl);
}
