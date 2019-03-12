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
double          t_duration, t_duration_hydro, t_duration_rt;
double          t_duration_transp, t_duration_react;
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

#if defined(_RT_)
    Chem_Data       chData;     // 12.30 RT use
    time_t          t_start, t_end, t_start_hydro, t_end_hydro;

    t_start = time(NULL);  // 12.30 timing
#endif

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
    chData = (Chem_Data)malloc(sizeof (*chData));     // 12.30 RT use
#endif

    /* Read PIHM input files */
    ReadAlloc(pihm, chData);

    /* Initialize CVode state variables */
    CV_Y = N_VNew(NumStateVar());
    if (CV_Y == NULL)
    {
        PIHMprintf(VL_ERROR, "Error creating CVODE state variable vector.\n");
        PIHMexit(EXIT_FAILURE);
    }

    /* Initialize PIHM structure */
    Initialize(pihm, CV_Y, &cvode_mem);
#if defined(_RT_)
    // 12.30 RT use, must be placed after Initialize()
    chem_alloc(project, pihm->filename.chem, pihm, chData);
#endif

    /* Create output directory */
    CreateOutputDir(outputdir);

    /* Create output structures */
#if defined(_CYCLES_)
    MapOutput(pihm->ctrl.prtvrbl, pihm->ctrl.tpprtvrbl, pihm->epctbl,
        pihm->elem, pihm->river, &pihm->meshtbl, outputdir, &pihm->print);
#elif defined(_RT_)
    MapOutput(pihm->ctrl.prtvrbl, pihm->ctrl.tpprtvrbl, chData,
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
        Spinup(pihm, chData, CV_Y, cvode_mem);
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
            PIHM(pihm, chData, cvode_mem, CV_Y, cputime);
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
    free(pihm);
#if defined(_RT_)
    // 12.30, RT use
    FreeChem(chData);  // 12.30, RT use
    free(chData);
    //
    // 12.30, RT use
    t_end = time(NULL);
    t_duration = t_end - t_start;
    fprintf(stderr, "Wall time of simulation = %.3f [min] or %.3f [hr]. \n", t_duration/60, t_duration/3600);
    fprintf(stderr, "Wall time of hydro step = %.3f [min] or %.3f [hr]. \n", t_duration_hydro/60, t_duration_hydro/3600);
    fprintf(stderr, "Wall time of rt step = %.3f [min] or %.3f [hr]. \n", t_duration_rt/60, t_duration_rt/3600);
    fprintf(stderr, "Wall time of rt_transp = %.3f [min] or %.3f [hr]. \n", t_duration_transp/60, t_duration_transp/3600);
    fprintf(stderr, "Wall time of rt_react = %.3f [min] or %.3f [hr]. \n", t_duration_react/60, t_duration_react/3600);
#endif

    PIHMprintf(VL_BRIEF, "\nSimulation completed.\n");

    return EXIT_SUCCESS;
}
