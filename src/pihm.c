#include "pihm.h"

void PIHM(pihm_struct pihm, void *cvode_mem, N_Vector CV_Y, int t,
    int next_t, double cputime)
{
    /*
     * Apply boundary conditions
     */
    ApplyBC(&pihm->forc, pihm->elem, pihm->riv, t);

    /* Determine if land surface simulation is needed */
    if ((t - pihm->ctrl.starttime) % pihm->ctrl.etstep == 0)
    {
        /* Apply forcing */
        ApplyForcing(&pihm->forc, pihm->elem, t
#ifdef _NOAH_
            , &pihm->ctrl, &pihm->siteinfo
#endif
            );

#ifdef _NOAH_
        /* Calculate surface energy balance */
        Noah(pihm);
#else
        /* Calculate Interception storage and ET */
        IntcpSnowET(t, (double)pihm->ctrl.etstep, pihm);
#endif
        /*
         * Update print variables for land surface step variables
         */
        UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, LS_STEP);
    }

    /*
     * Solve PIHM hydrology ODE using CVODE
     */
    SolveCVode(pihm->ctrl.starttime, &t, next_t, cputime, cvode_mem, CV_Y);

    /* Use mass balance to calculate model fluxes or variables */
    Summary(pihm, CV_Y, (double)pihm->ctrl.stepsize);

    if (pihm->ctrl.waterB)
    {
        /* Print water balance */
        PrintWaterBalance(pihm->print.walbal_file, t, pihm->ctrl.starttime,
            pihm->ctrl.stepsize, pihm->elem, nelem, pihm->riv, nriver);
    }
#ifdef _NOAH_
    NoahHydrol(pihm->elem, (double)pihm->ctrl.stepsize);
#endif

#ifdef _CYCLES_
    SoluteTransport(pihm->elem, pihm->riv, (double)pihm->ctrl.stepsize);
#endif

    /*
     * Update print variables for hydrology step variables
     */
    UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, HYDROL_STEP);

#ifdef _BGC_
    int             i;

    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
        for (i = 0; i < nelem; i++)
        {
            /* Test for nitrogen balance */
            CheckNitrogenBalance(&pihm->elem[i].ns,
                &pihm->elem[i].epv.old_n_balance);
        }
    }
#endif

    /*
     * Daily timestep modules
     */
#ifdef _DAILY_
    DailyVar(t, pihm->ctrl.starttime, pihm);

    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
#ifdef _BGC_
        DailyBgc(pihm, t - DAYINSEC);

        first_balance = 0;
#endif

#ifdef _CYCLES_
        DailyCycles(t - DAYINSEC, pihm);
#endif

        /*
         * Update print variables for CN (daily) step variables
         */
        UpdPrintVar(pihm->print.varctrl, pihm->print.nprint, CN_STEP);
    }
#endif

#ifdef _DAILY_
    /*
     * Initialize daily structures
     */
    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
        InitDailyStruct(pihm);
    }
#endif
    /*
     * Print outputs
     */
    PrintData(pihm->print.varctrl, pihm->print.nprint, t,
        t - pihm->ctrl.starttime, pihm->ctrl.ascii);
    if (pihm->ctrl.tecplot)
    {
        UpdPrintVarT(pihm->print.tp_varctrl, pihm->print.ntpprint);
        PrintDataTecplot(pihm->print.tp_varctrl, pihm->print.ntpprint, t,
            t - pihm->ctrl.starttime);
    }
}
