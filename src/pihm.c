#include "pihm.h"

void PIHM (char *simulation, char *outputdir, int first_cycle
#ifdef _ENKF_
    , int starttime, int endtime, int startmode, double *param
#endif
    )
{
    pihm_struct     pihm;
    N_Vector        CV_Y, abstol;       /* State Variables Vector and Absolute Tolerance*/
    void           *cvode_mem;			/* Model Data Pointer */
    int             nsv;				/* Problem size */
    int             i;					/* Loop index */
    int             t;					/* Simulation time */
    realtype        cputime, cputime_dt;/* Time cpu duration */
    clock_t         ptime, start, ctime;
#ifdef _OPENMP
	realtype		ptime_omp, start_omp, ctime_omp;
#endif
	int			    num_threads;       /* Number of threads for OMP */
	long int		nst, nfe, nfeLS, nni, ncfn, netf; /*Variables for monitoring performance */
	int				flag;
	char			WBname[50];
	FILE            *WaterBalance; /* Water balance file */


#ifdef _BGC_
    int             first_balance = 1;
#endif

	/* Set the number of threads to use */
	num_threads = 1;     /* default value */
#ifdef _OPENMP
	num_threads = omp_get_max_threads();  /* Overwrite with OMP_NUM_THREADS environment variable */
	printf("\nOpenMP with %d threads\n", num_threads);
#endif

    /* Allocate memory for model data structure */
    pihm = (pihm_struct)malloc (sizeof *pihm);

    /* Read PIHM input files */
    ReadAlloc (simulation, pihm);

	sprintf(WBname, "%s%s_WaterBalance.plt", outputdir, simulation);
	WaterBalance = fopen(WBname, "w");
	CheckFile(WaterBalance, WBname);

#ifdef _ENKF_
    /* When running in ensemble mode, use parameters and calibration
     * determined by EnKF module */
    pihm->ctrl.init_type = startmode;
    pihm->ctrl.starttime = starttime;
    pihm->ctrl.endtime = endtime;

    Mbr2Cal (&pihm->cal, param);
#endif

    /* Number of state variables */
    nsv = 3 * pihm->numele + 2 * pihm->numriv;

    /* Initialize CVode state variables */

#ifdef _OPENMP
	CV_Y = N_VNew_OpenMP(nsv, num_threads);
	abstol = N_VNew_OpenMP(nsv, num_threads);
#else
	CV_Y = N_VNew_Serial(nsv);
	abstol = N_VNew_Serial(nsv);
#endif

    /* Initialize PIHM structure */
    Initialize (pihm, CV_Y);

	/* Set the vector absolute tolerance */
	for (i = 0; i < nsv; i++) {
#ifdef _OPENMP
		NV_Ith_OMP(abstol, i) = pihm->ctrl.abstol;
#else
		NV_Ith_S(abstol, i) = pihm->ctrl.abstol;
#endif
	}

    /* Allocate memory for solver */
    cvode_mem = CVodeCreate (CV_BDF, CV_NEWTON);
    if (cvode_mem == NULL)
    {
        PIHMprintf (VL_ERROR, "Error in allocating memory for solver.\n");
        PIHMexit (EXIT_FAILURE);
    }

    /* Create output structures */
    MapOutput (simulation, pihm, outputdir);

    if (first_cycle)
    {
        /* Backup input files */

		#if !defined(_MSC_VER)
			BKInput(simulation, outputdir);
		#endif

        /* Initialize output files and structures */
        InitOutputFile (pihm->prtctrl, pihm->ctrl.nprint, pihm->ctrl.ascii);
    }

    PIHMprintf (VL_VERBOSE, "\n\nSolving ODE system ... \n\n");

    /* Set solver parameters */
    SetCVodeParam (pihm, cvode_mem, CV_Y, abstol);

#ifdef _OPENMP
	start_omp = omp_get_wtime();
	ptime_omp = start_omp;
#else
	start = clock();
	ptime = start;
#endif


    /* Start solver in loops */
    for (i = 0; i < pihm->ctrl.nstep; i++)
    {
        /* Determine current step and next step */
        t = pihm->ctrl.tout[i];

        /* Apply forcing */
        ApplyForcing (&pihm->forc, pihm->elem, pihm->numele, pihm->riv,
            pihm->numriv, t
#ifdef _BGC_
            , &pihm->ctrl
#endif
            );

        /* Determine if land surface simulation is needed */
        if ((t - pihm->ctrl.starttime) % pihm->ctrl.etstep == 0)
        {
#ifdef _NOAH_
            /* Calculate surface energy balance */
            Noah (t, pihm);
#else
            /* Calculate Interception storage and ET */
            IntcpSnowET (t, (double)pihm->ctrl.etstep, pihm);
#endif
        }

        /*
         * Solve PIHM hydrology ODE using CVODE
         */

         /* CPU time */
#ifdef _OPENMP
		ctime_omp = omp_get_wtime();
		cputime_dt = ((double)(ctime_omp - ptime_omp));
		cputime = ((double)(ctime_omp - start_omp));
		ptime_omp = ctime_omp;
#else
		ctime = clock();
		cputime_dt = ((double)(ctime - ptime)) / CLOCKS_PER_SEC;
		cputime = ((double)(ctime - start)) / CLOCKS_PER_SEC;
		ptime = ctime;
#endif
		
         SolveCVode (pihm->ctrl.starttime, &t, pihm->ctrl.tout[i + 1], pihm->ctrl.stepsize, cputime_dt, cputime,
            cvode_mem, CV_Y, pihm->ctrl.cvode_perf, simulation, outputdir);

		 if (pihm->ctrl.waterB)
		 {

			 /* Print water balance */
			 PrintWaterBalance(WaterBalance, (pihm->ctrl.tout[i + 1]-pihm->ctrl.starttime), pihm->elem, pihm->numele, pihm->riv, pihm->numriv, i);
		 }

        /*
         * Use mass balance to calculate model fluxes or variables
         */
        Summary (pihm, CV_Y, (double)pihm->ctrl.stepsize);

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
                Save2Stor (pihm, t, pihm->ctrl.spinupstart,
                    pihm->ctrl.spinupend);
            }
            else
            {
                DailyBgc (pihm, t, pihm->ctrl.starttime, first_balance);

                first_balance = 0;
            }
#endif
        }
#endif

#ifdef _CYCLES_
        SoluteTransport (pihm->elem, pihm->numele, pihm->riv, pihm->numriv,
            (double)pihm->ctrl.stepsize);
#endif

        /*
         * Print outputs
         */

		if (pihm->ctrl.tecplot)
		{
			PrintDataTecplot(pihm->prtctrlT, pihm->ctrl.nprintT, (pihm->ctrl.tout[i + 1] - pihm->ctrl.starttime),
				t - pihm->ctrl.starttime, pihm->ctrl.stepsize, i);
		}

        PrintData (pihm->prtctrl, pihm->ctrl.nprint, t,
            t - pihm->ctrl.starttime, pihm->ctrl.stepsize, pihm->ctrl.ascii);
		

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

#ifdef _BGC_
    if (pihm->ctrl.bgc_spinup)
    {
        BGCSpinup (simulation, pihm, outputdir);
    }
#endif

    /*
     * Write init files
     */
    if (pihm->ctrl.write_ic)
    {
        PrtInit (pihm, simulation);
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

    /* Free memory */
#ifdef _OPENMP
	N_VDestroy_OpenMP(CV_Y);
	N_VDestroy_OpenMP(abstol);
#else
	N_VDestroy_Serial(CV_Y);
	N_VDestroy_Serial(abstol);
#endif
	fclose(WaterBalance);
    /* Free integrator memory */
    CVodeFree (&cvode_mem);
    FreeData (pihm);
    free (pihm);
}
