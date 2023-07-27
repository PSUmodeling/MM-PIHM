#include "pihm.h"

void PIHM(double cputime, pihm_struct pihm, void *cvode_mem, N_Vector CV_Y)
{
    int             t;
#if defined(_RT_)
    const int       SPECIATION_STEP = 3600;
#endif

    t = pihm->ctrl.tout[pihm->ctrl.cstep];

    // Apply boundary conditions
#if defined(_RT_)
    ApplyBc(t, &pihm->rttbl, &pihm->forc, pihm->elem, pihm->river);
#else
    ApplyBc(t, &pihm->forc, pihm->elem, pihm->river);
#endif

    // Apply forcing and simulate land surface processes
    if ((t - pihm->ctrl.starttime) % pihm->ctrl.etstep == 0)
    {
        // Apply forcing
#if defined(_RT_)
        ApplyForcing(t, pihm->ctrl.rad_mode, &pihm->siteinfo, &pihm->rttbl, &pihm->forc, pihm->elem);
#elif defined(_NOAH_)
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

#if defined(_RT_)
    // Reaction
    if (pihm->rttbl.transpt_flag == KIN_REACTION)
    {
        if ((t - pihm->ctrl.starttime) % pihm->ctrl.AvgScl == 0)
        {
            Reaction((double)pihm->ctrl.AvgScl, pihm->chemtbl, pihm->kintbl, &pihm->rttbl, pihm->elem);
        }
    }
#endif

#if defined(_CYCLES_)
    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
        Cycles(t, &pihm->co2ctrl, &pihm->forc, pihm->elem);

        // Update print variables for CN (daily) step variables
        UpdatePrintVar(pihm->print.nprint, CN_STEP, pihm->print.varctrl);
    }
#endif

    // Solve PIHM hydrology ODE using CVode
    SolveCVode(cputime, &pihm->ctrl, &t, cvode_mem, CV_Y);

    // Use mass balance to calculate model fluxes or variables
    UpdateVar((double)pihm->ctrl.stepsize, pihm->elem, pihm->river, CV_Y);

#if defined(_NOAH_)
    NoahHydrol((double)pihm->ctrl.stepsize, pihm->elem);
#endif

    // Update print variables for hydrology step variables
    UpdatePrintVar(pihm->print.nprint, HYDROL_STEP, pihm->print.varctrl);

#if defined(_RT_)
    // Update chemical concentrations
    if (pihm->rttbl.transpt_flag == KIN_REACTION)
    {
        if ((t - pihm->ctrl.starttime) % SPECIATION_STEP == 0)
        {
            // Speciation
            RiverSpeciation(pihm->chemtbl, &pihm->rttbl, pihm->river);
        }
    }
    else
    {
        UpdatePrimConc(&pihm->rttbl, pihm->elem, pihm->river);
    }

    UpdatePrintVar(pihm->print.nprint, RT_STEP, pihm->print.varctrl);
#endif

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
        PrintWaterBalance(t, pihm->ctrl.starttime, pihm->ctrl.stepsize, pihm->elem, pihm->river,
            pihm->print.watbal_file);
    }

    // Print binary and txt output files
    PrintData(pihm->print.nprint, t, t - pihm->ctrl.starttime, pihm->ctrl.ascii, pihm->print.varctrl);
}
