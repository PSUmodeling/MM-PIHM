#include "pihm.h"

void PIHM (char *simulation, char *outputdir, int first_cycle
#ifdef _ENKF_
    , int starttime, int endtime, int startmode, double *param
#endif
    )
{
    pihm_struct     pihm;
    N_Vector        CV_Y;       /* State Variables Vector */
    void           *cvode_mem;  /* Model Data Pointer */
    int             nsv;        /* Problem size */
    int             i;          /* Loop index */
    int             t;          /* Simulation time */
#ifdef _BGC_
    int             first_balance = 1;
#endif

    /* Allocate memory for model data structure */
    pihm = (pihm_struct)malloc (sizeof *pihm);

    /* Read PIHM input files */
    ReadAlloc (simulation, pihm);

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
    CV_Y = N_VNew_Serial (nsv);

    /* Initialize PIHM structure */
    Initialize (pihm, CV_Y);

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
        BKInput (simulation, outputdir);

        /* Initialize output files and structures */
        InitOutputFile (pihm->prtctrl, pihm->ctrl.nprint, pihm->ctrl.ascii);
    }

    PIHMprintf (VL_VERBOSE, "\n\nSolving ODE system ... \n\n");

    /* Set solver parameters */
    SetCVodeParam (pihm, cvode_mem, CV_Y);

    /* Start solver in loops */
    for (i = 0; i < pihm->ctrl.nstep; i++)
    {
        /* Determine current step and next step */
        t = pihm->ctrl.tout[i];

        /* Apply forcing */
        ApplyForcing (&pihm->forc, pihm->elem, pihm->numele, pihm->riv,
            pihm->numriv, t
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
        SolveCVode (&t, pihm->ctrl.tout[i + 1], pihm->ctrl.stepsize,
            cvode_mem, CV_Y);

        /*
         * Use mass balance to calculate model fluxes or variables
         */
        Summary (pihm, CV_Y, (double)pihm->ctrl.stepsize);

#ifdef _NOAH_
        NoahHydrol (pihm, (double)pihm->ctrl.stepsize);
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
        BGCSpinup (pihm, outputdir);
    }
#endif

    /*
     * Write init files
     */
    if (pihm->ctrl.write_ic)
    {
        PrtInit (pihm, simulation);

#ifdef _BGC_
        if (!pihm->ctrl.bgc_spinup)
        {
            WriteBGCIC (pihm->filename.bgcic, pihm->elem, pihm->numele,
                pihm->riv, pihm->numriv);
        }
#endif
    }

    /* Free memory */
    N_VDestroy_Serial (CV_Y);

    /* Free integrator memory */
    CVodeFree (&cvode_mem);

    FreeData (pihm);
    free (pihm);
}
