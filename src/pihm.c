#include "pihm.h"

#ifdef _ENKF_
void PIHMRun (char *simulation, char *outputdir, int first_cycle,
    int starttime, int endtime, int startmode, double *param)
#else
void PIHMRun (char *simulation, char *outputdir, int first_cycle)
#endif
{
    pihm_struct     pihm;
    N_Vector        CV_Y;       /* State Variables Vector */
    void           *cvode_mem;  /* Model Data Pointer */
    int             nsv;        /* Problem size */
    int             i;          /* Loop index */
    int             t;          /* Simulation time */

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
        printf ("Fatal error: CVodeMalloc failed. \n");
        PihmExit (1);
    }

    /* Create output structures */
    MapOutput (simulation, pihm, outputdir);

    if (first_cycle == 1)
    {
        /* Backup input files */
        BKInput (simulation, outputdir);

        /* Initialize output files and structures */
        InitOutputFile (pihm->prtctrl, pihm->ctrl.nprint, pihm->ctrl.ascii);
    }

    if (verbose_mode)
    {
        printf ("\n\nSolving ODE system ... \n\n");
    }

    /* Set solver parameters */
    SetCVodeParam (pihm, cvode_mem, CV_Y);

    /* Start solver in loops */
    for (i = 0; i < pihm->ctrl.nstep; i++)
    {
        /* Determine current step and next step */
        t = pihm->ctrl.tout[i];

        /* Apply forcing */
        ApplyForcing (&pihm->forc, pihm->elem, pihm->numele, pihm->riv, pihm->numriv, t);

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
        SolveCVode (&t, pihm->ctrl.tout[i + 1], pihm->ctrl.stepsize,
            cvode_mem, CV_Y);

        /*
         * Use mass balance to calculate model fluxes or variables
         */
        Summary (pihm, CV_Y, (double)pihm->ctrl.stepsize);

#ifdef _DAILY_
        DailyVar (t, pihm->ctrl.starttime, pihm);
#endif

        /*
         * Daily timestep modules
         */
#ifdef _DAILY_
        if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
        {
    #ifdef _CYCLES_
            DailyCycles (t - DAYINSEC, pihm);
    #endif
    #ifdef _BGC_
            if (pihm->ctrl.spinup)
            {
                Save2Stor (pihm, t, pihm->ctrl.spinupstart, pihm->ctrl.spinupend);
            }
            else
            {
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
    if (pihm->ctrl.spinup)
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

    /* Free memory */
    N_VDestroy_Serial (CV_Y);

    /* Free integrator memory */
    CVodeFree (&cvode_mem);

    FreeData (pihm);
    free (pihm);
}

void SetCVodeParam (pihm_struct pihm, void *cvode_mem, N_Vector CV_Y)
{
    int             flag;

    flag = CVodeSetFdata (cvode_mem, pihm);
    flag = CVodeSetInitStep (cvode_mem, (realtype) pihm->ctrl.initstep);
    flag = CVodeSetStabLimDet (cvode_mem, TRUE);
    flag = CVodeSetMaxStep (cvode_mem, (realtype) pihm->ctrl.maxstep);
    flag = CVodeMalloc (cvode_mem, Hydrol, (realtype) pihm->ctrl.starttime,
        CV_Y, CV_SS, (realtype) pihm->ctrl.reltol, &pihm->ctrl.abstol);
    flag = CVSpgmr (cvode_mem, PREC_NONE, 0);
}

void SolveCVode (int *t, int nextptr, int stepsize, void *cvode_mem,
    N_Vector CV_Y)
{
    realtype        solvert;
    realtype        cvode_val;
    struct tm      *timestamp;
    time_t          rawtime;
    int             flag;

    solvert = (realtype) (*t);

    flag = CVodeSetMaxNumSteps (cvode_mem, (long int)(stepsize * 20));
    flag = CVodeSetStopTime (cvode_mem, (realtype) nextptr);
    flag = CVode (cvode_mem, (realtype) nextptr, CV_Y, &solvert,
        CV_NORMAL_TSTOP);
    flag = CVodeGetCurrentTime (cvode_mem, &cvode_val);

    *t = (int)solvert;
    rawtime = (time_t) (*t);
    timestamp = gmtime (&rawtime);

    if (verbose_mode)
    {
        printf (" Step = %4.4d-%2.2d-%2.2d %2.2d:%2.2d (%d)\n",
            timestamp->tm_year + 1900, timestamp->tm_mon + 1,
            timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min, *t);
    }
#ifndef _ENKF_
    else if (rawtime % 3600 == 0)
    {
        printf (" Step = %4.4d-%2.2d-%2.2d %2.2d:%2.2d\n",
            timestamp->tm_year + 1900, timestamp->tm_mon + 1,
            timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
    }
#endif
}
