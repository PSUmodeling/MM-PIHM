#include "pihm.h"

#if defined(_RT_)
void PIHM(pihm_struct pihm, Chem_Data rt, void *cvode_mem, N_Vector CV_Y,
    double cputime)
#else
void PIHM(pihm_struct pihm, void *cvode_mem, N_Vector CV_Y, double cputime)
#endif
{
    int             t;
#if defined(_RT_)
    time_t          t_start_rt, t_end_rt;
    time_t          t_end_hydro, t_start_hydro;
#endif

    t = pihm->ctrl.tout[pihm->ctrl.cstep];

#if defined(_RT_)
    t_start_hydro = time(NULL);                        // 12.30 timing hydro
#endif
    /* Apply boundary conditions */
    ApplyBc(&pihm->forc, pihm->elem, pihm->river, t);

    /*
     * Apply forcing and simulate land surface processes
     */
    if ((t - pihm->ctrl.starttime) % pihm->ctrl.etstep == 0)
    {
        /* Apply forcing */
#if defined(_NOAH_)
        ApplyForc(&pihm->forc, pihm->elem, t , pihm->ctrl.rad_mode,
            &pihm->siteinfo);
#else
        ApplyForc(&pihm->forc, pihm->elem, t);
#endif

#if defined(_NOAH_)
        /* Calculate surface energy balance */
        Noah(pihm->elem, (double)pihm->ctrl.etstep);
#else
        /* Calculate Interception storage and ET */
        IntcpSnowEt(t, (double)pihm->ctrl.etstep, pihm->elem, &pihm->cal);
#endif

        /* Update print variables for land surface step variables */
        UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, LS_STEP);
        UpdPrintVar(pihm->print.tp_varctrl, pihm->print.ntpprint, LS_STEP);
    }

    /*
     * Solve PIHM hydrology ODE using CVode
     */
    SolveCVode(pihm->ctrl.starttime, &t, pihm->ctrl.tout[pihm->ctrl.cstep + 1],
        cputime, cvode_mem, CV_Y);

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

#if defined(_RT_)
    t_end_hydro = time(NULL);                          // 12.30 timing hydro
    t_duration_hydro += t_end_hydro - t_start_hydro;   // 12.30 timing hydro
#endif

#if defined(_RT_)
    t_start_rt = time(NULL);
    fluxtrans(pihm->ctrl.tout[pihm->ctrl.cstep + 1]/60, pihm->ctrl.stepsize/60,
        pihm, rt, &t_duration_transp, &t_duration_react);
    UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, RT_STEP);
    UpdPrintVar(pihm->print.tp_varctrl, pihm->print.ntpprint, RT_STEP);
    t_end_rt = time(NULL);
    t_duration_rt += t_end_rt - t_start_rt;

    //PrintChem(outputdir, project, rt, pihm->ctrl.tout[pihm->ctrl.cstep + 1]/60);
#endif

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
