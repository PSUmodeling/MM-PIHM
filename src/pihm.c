#include "pihm.h"

void PIHM(pihm_struct pihm, void *cvode_mem, N_Vector CV_Y, int t,
    int next_t, double cputime)
{
    /* Apply boundary conditions */
    ApplyBc(&pihm->forc, pihm->elem, pihm->river, t);

    /* Determine if land surface simulation is needed */
    if ((t - pihm->ctrl.starttime) % pihm->ctrl.etstep == 0)
    {
        /* Apply forcing */
#ifdef _NOAH_
        ApplyForc(&pihm->forc, pihm->elem, t , pihm->ctrl.rad_mode,
            &pihm->siteinfo);
#else
        ApplyForc(&pihm->forc, pihm->elem, t);
#endif

#ifdef _NOAH_
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

    /* Solve PIHM hydrology ODE using CVode */
    SolveCVode(pihm->ctrl.starttime, &t, next_t, cputime, cvode_mem, CV_Y);

    /* Use mass balance to calculate model fluxes or variables */
    Summary(pihm->elem, pihm->river, CV_Y, (double)pihm->ctrl.stepsize);

#ifdef _NOAH_
    NoahHydrol(pihm->elem, (double)pihm->ctrl.stepsize);
#endif

#if defined(_BGC_) && defined(_LUMPED_)
    NLeaching(pihm->elem, pihm->river, (double)pihm->ctrl.stepsize);
#endif

#ifdef _CYCLES_
    SoluteTransport(pihm->elem, pihm->river, (double)pihm->ctrl.stepsize);
#endif

    /* Update print variables for hydrology step variables */
    UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, HYDROL_STEP);
    UpdPrintVar(pihm->print.tp_varctrl, pihm->print.ntpprint, HYDROL_STEP);

#if defined(_BGC_)
    int             i;

    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
# if defined(_LUMPED_)
        i = LUMPED;
# else
#  ifdef _OPENMP
#   pragma omp parallel for
#  endif
        for (i = 0; i < nelem; i++)
# endif
        {
            /* Test for nitrogen balance */
            CheckNitrogenBalance(&pihm->elem[i].ns,
                &pihm->elem[i].epv.old_n_balance);
        }
    }
#endif

    /* Daily timestep modules */
#ifdef _DAILY_
    DailyVar(t, pihm->ctrl.starttime, pihm->elem, pihm->ctrl.stepsize);

    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
# ifdef _BGC_
        DailyBgc(pihm, t - DAYINSEC);
        first_balance = 0;
# endif

# ifdef _CYCLES_
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
