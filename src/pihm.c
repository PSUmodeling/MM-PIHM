#include "pihm.h"

void PIHM (pihm_struct pihm, void *cvode_mem, N_Vector CV_Y, N_Vector abstol, int t,
    int next_t, char *outputdir, char *simulation, double cputime_dt, double cputime)

{
	long int		nst, nfe, nfeLS, nni, ncfn, netf; /*Variables for monitoring performance */
	int				flag;
	char			WBname[100];
	FILE            *WaterBalance; /* Water balance file */
    /*
     * Apply boundary conditions
     */
    ApplyBC (&pihm->forc, pihm->elem, pihm->riv, t);

	sprintf(WBname, "%s%s_WaterBalance.plt", outputdir, simulation);
	WaterBalance = fopen(WBname, "w");
	CheckFile(WaterBalance, WBname);

    /* Determine if land surface simulation is needed */
    if ((t - pihm->ctrl.starttime) % pihm->ctrl.etstep == 0)
    {
        /* Apply forcing */
        ApplyForcing (&pihm->forc, pihm->elem, t
    #ifdef _NOAH_
            , &pihm->ctrl, pihm->latitude, pihm->longitude, pihm->elevation,
            pihm->noahtbl.tbot
    #endif
            );

#ifdef _NOAH_
        /* Calculate surface energy balance */
        Noah (pihm);
#else
        /* Calculate Interception storage and ET */
        IntcpSnowET (t, (double)pihm->ctrl.etstep, pihm);
#endif
        /*
         * Update print variables for land surface step variables
         */
        UpdPrintVar (pihm->prtctrl, pihm->ctrl.nprint, LS_STEP);
    }

    /*
     * Solve PIHM hydrology ODE using CVODE
     */		
         SolveCVode (pihm->ctrl.starttime, &t, next_t, pihm->ctrl.stepsize, cputime_dt, cputime,
            cvode_mem, CV_Y, pihm->ctrl.cvode_perf, simulation, outputdir);

		 if (pihm->ctrl.waterB)
		 {

			 /* Print water balance */
			/* PrintWaterBalance(WaterBalance, (next_t-pihm->ctrl.starttime), pihm->elem, nelem, pihm->riv, nriver, t);*/
		 }

    /* Use mass balance to calculate model fluxes or variables */
    Summary (pihm, CV_Y, (double)pihm->ctrl.stepsize);

#ifdef _NOAH_
    NoahHydrol (pihm->elem, (double)pihm->ctrl.stepsize);
#endif

#ifdef _CYCLES_
    SoluteTransport (pihm->elem, pihm->riv, (double)pihm->ctrl.stepsize);
#endif

    /*
     * Update print variables for hydrology step variables
     */
    UpdPrintVar (pihm->prtctrl, pihm->ctrl.nprint, HYDROL_STEP);

#ifdef _BGC_
    int             i;

    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
        for (i = 0; i < nelem; i++)
        {
            /* Test for nitrogen balance */
            CheckNitrogenBalance (&pihm->elem[i].ns,
                &pihm->elem[i].epv.old_n_balance);
        }
    }
#endif

    /*
     * Daily timestep modules
     */
#ifdef _DAILY_
    DailyVar (t, pihm->ctrl.starttime, pihm);

    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
#ifdef _BGC_
        DailyBgc (pihm, t - DAYINSEC);

        first_balance = 0;
#endif

#ifdef _CYCLES_
        DailyCycles (t - DAYINSEC, pihm);
#endif

        /*
         * Update print variables for CN (daily) step variables
         */
        UpdPrintVar (pihm->prtctrl, pihm->ctrl.nprint, CN_STEP);
    }
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
        /*
         * Print outputs
         */

		if (pihm->ctrl.tecplot)
		{
			PrintDataTecplot(pihm->prtctrlT, pihm->ctrl.nprintT, (next_t - pihm->ctrl.starttime),
				t - pihm->ctrl.starttime, pihm->ctrl.stepsize, t);
		}

        PrintData (pihm->prtctrl, pihm->ctrl.nprint, t,
            t - pihm->ctrl.starttime, pihm->ctrl.stepsize, pihm->ctrl.ascii);

		if (pihm->ctrl.prtvrbl[IC_CTRL] > 0)
		{
			/*
			* Write init files
			*/
			if (pihm->ctrl.write_ic && (t - pihm->ctrl.starttime) % pihm->ctrl.prtvrbl[IC_CTRL] == 0)
			{
				PrtInit(pihm, simulation, t);
			}
		}


	flag = CVodeGetNumSteps(cvode_mem, &nst);
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);	
	flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
	flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
	flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

	printf("nst = %-6ld nfe  = %-6ld nfeLS = %-6ld\n",
		nst, nfe,  nfeLS);
	printf("nni = %-6ld ncfn = %-6ld netf = %-6ld\n \n",
		nni, ncfn, netf);


	fclose(WaterBalance);
}
