#include "pihm.h"

/* Global variables */
int             verbose_mode;
int             debug_mode;
int             append_mode;
int             corr_mode;
int             spinup_mode;
int             fixed_length;
int             tecplot;
char            project[MAXSTRING];
int             nelem;
int             nriver;
#if defined(_OPENMP)
int             nthreads = 1;    /* Default value */
#endif
#if defined(_BGC_)
int             first_balance;
#endif
#if defined(_RT_)
int             NumSpc;
#endif

int main(int argc, char *argv[])
{
    char            outputdir[MAXSTRING];
    pihm_struct     pihm;
    ctrl_struct    *ctrl;
    N_Vector        CV_Y;
    void           *cvode_mem;
#if defined(_OPENMP)
    double          start_omp;
#else
    clock_t         start;
#endif
    double          cputime, cputime_dt;    /* Time cpu duration */

#if defined(_OPENMP)
    /* Set the number of threads to use */
    nthreads = omp_get_max_threads();
#endif

    memset(outputdir, 0, MAXSTRING);

    /* Read command line arguments */
    ParseCmdLineParam(argc, argv, outputdir);

    /* Print AscII art */
    StartupScreen();

    /* Allocate memory for model data structure */
    pihm = (pihm_struct)malloc(sizeof(*pihm));
#if defined(_RT_)
    pihm->rt = (Chem_Data)malloc(sizeof (*pihm->rt));     // 12.30 RT use
#endif

    /* Read PIHM input files */
    ReadAlloc(pihm, pihm->rt);

    /* Initialize CVode state variables */
    CV_Y = N_VNew(NumStateVar());
    if (CV_Y == NULL)
    {
        PIHMprintf(VL_ERROR, "Error creating CVODE state variable vector.\n");
        PIHMexit(EXIT_FAILURE);
    }

    /* Initialize PIHM structure */
#if defined(_RT_)
    Initialize(pihm, pihm->rt, CV_Y, &cvode_mem);
#else
    Initialize(pihm, CV_Y, &cvode_mem);
#endif

    /* Create output directory */
    CreateOutputDir(outputdir);

    /* Create output structures */
#if defined(_CYCLES_)
    MapOutput(pihm->ctrl.prtvrbl, pihm->ctrl.tpprtvrbl, pihm->epctbl,
        pihm->elem, pihm->river, &pihm->meshtbl, outputdir, &pihm->print);
#elif defined(_RT_)
    MapOutput(pihm->ctrl.prtvrbl, pihm->ctrl.tpprtvrbl, pihm->chemtbl,
        &pihm->rttbl, pihm->rt,
        pihm->elem, pihm->river, &pihm->meshtbl, outputdir, &pihm->print);
#else
    MapOutput(pihm->ctrl.prtvrbl, pihm->ctrl.tpprtvrbl, pihm->elem, pihm->river,
        &pihm->meshtbl, outputdir, &pihm->print);
#endif

    /* Backup input files */
#if !defined(_MSC_VER)
    if (!append_mode)
    {
        BackupInput(outputdir, &pihm->filename);
    }
#endif

    InitOutputFile(&pihm->print, outputdir, pihm->ctrl.waterbal,
        pihm->ctrl.ascii);

    PIHMprintf(VL_VERBOSE, "\n\nSolving ODE system ... \n\n");

    /* Set solver parameters */
    SetCVodeParam(pihm, cvode_mem, CV_Y);

#if defined(_BGC_)
    first_balance = 1;
#endif

    /*
     * Run PIHM
     */
#if defined(_OPENMP)
    start_omp = omp_get_wtime();
#else
    start = clock();
#endif

    ctrl = &pihm->ctrl;

    if (spinup_mode)
    {
#if defined(_RT_)
        Spinup(pihm, pihm->rt, CV_Y, cvode_mem);
#else
        Spinup(pihm, CV_Y, cvode_mem);
#endif

        /* In spin-up mode, initial conditions are always printed */
        PrintInit(pihm->elem, pihm->river, outputdir,
            ctrl->endtime, ctrl->starttime,
            ctrl->endtime, ctrl->prtvrbl[IC_CTRL]);
#if defined(_BGC_)
        WriteBgcIc(outputdir, pihm->elem, pihm->river);
#endif
    }
    else
    {
        for (ctrl->cstep = 0; ctrl->cstep < ctrl->nstep; ctrl->cstep++)
        {
#if defined(_OPENMP)
            RunTime(start_omp, &cputime, &cputime_dt);
#else
            RunTime(start, &cputime, &cputime_dt);
#endif

#if defined(_RT_)
            PIHM(pihm, pihm->rt, cvode_mem, CV_Y, cputime);
#else
            PIHM(pihm, cvode_mem, CV_Y, cputime);
#endif

            /* Adjust CVODE max step to reduce oscillation */
            AdjCVodeMaxStep(cvode_mem, &pihm->ctrl);

            /* Print CVODE performance and statistics */
            if (debug_mode)
            {
                PrintPerf(cvode_mem, ctrl->tout[ctrl->cstep + 1],
                    ctrl->starttime, cputime_dt, cputime,
                    ctrl->maxstep, pihm->print.cvodeperf_file);
            }

            /* Write init files */
            if (ctrl->write_ic)
            {
                PrintInit(pihm->elem, pihm->river, outputdir,
                    ctrl->tout[ctrl->cstep + 1], ctrl->starttime,
                    ctrl->endtime, ctrl->prtvrbl[IC_CTRL]);
            }

        }

#if defined(_BGC_)
        if (ctrl->write_bgc_restart)
        {
            WriteBgcIc(outputdir, pihm->elem, pihm->river);
        }
#endif

# if TEMP_DISABLED
#if defined(_CYCLES_)
        if (ctrl->write_cycles_restart)
        {
            WriteCyclesIC(pihm->filename.cyclesic, pihm->elem, pihm->river);
        }
# endif
#endif
    }

    if (debug_mode)
    {
        PrintCVodeFinalStats(cvode_mem);
    }

    /* Free memory */
    N_VDestroy(CV_Y);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);
    FreeMem(pihm);
#if defined(_RT_)
    FreeChem(pihm->rt);
    free(pihm->rt);
#endif
    free(pihm);

    PIHMprintf(VL_BRIEF, "\nSimulation completed.\n");

    return EXIT_SUCCESS;
}
