#include "pihm.h"

// Global variables
int             verbose_mode;
int             debug_mode;
int             append_mode;
int             corr_mode;
int             spinup_mode;
int             fixed_length;
char            project[MAXSTRING];
int             nelem;
int             nriver;
#if defined(_OPENMP)
int             nthreads = 1;               // Default value
#endif
#if defined(_BGC_)
int             nsolute = 1;
#elif defined(_CYCLES_)
int             nsolute = 2;
#elif defined(_RT_)
int             nsolute;
#endif
#if defined(_BGC_)
int             first_balance;
#endif

int main(int argc, char *argv[])
{
    char            outputdir[MAXSTRING];
    pihm_struct     pihm;
    ctrl_struct    *ctrl;
    N_Vector        CV_Y;
    void           *cvode_mem;
    SUNLinearSolver sun_ls;
#if defined(_OPENMP)
    double          start_omp;
#else
    clock_t         start;
#endif
    double          cputime, cputime_dt;    // Time cpu duration

#if defined(unix) || defined(__unix__) || defined(__unix)
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

#if defined(_OPENMP)
    // Set the number of threads to use
    nthreads = omp_get_max_threads();
#endif

#if defined(_DEBUG_)
    // When in debug mode, print PID and host name to a text file for gdb to attach to
    char            hostname[MAXSTRING];
    FILE           *fp;

    gethostname(hostname, MAXSTRING);
    fp = fopen("debug.txt", "w");
    fprintf(fp, "pid %d at %s\n", getpid(), hostname);
    fflush(fp);
    fclose(fp);
#endif


    memset(outputdir, 0, MAXSTRING);

    // Read command line arguments
    ParseCmdLineParam(argc, argv, outputdir);

    // Print AscII art
    StartupScreen();

    // Allocate memory for model data structure
    pihm = (pihm_struct)malloc(sizeof(*pihm));

    // Read PIHM input files
    ReadAlloc(pihm);

    // Initialize CVODE state variables
    CV_Y = N_VNew(NumStateVar());
    if (CV_Y == NULL)
    {
        pihm_printf(VL_ERROR, "Error creating CVODE state variable vector.\n");
        pihm_exit(EXIT_FAILURE);
    }

    // Initialize PIHM structure
    Initialize(pihm, CV_Y, &cvode_mem);

    // Create output directory
    CreateOutputDir(outputdir);

    // Create output structures
#if defined(_CYCLES_)
    MapOutput(outputdir, pihm->ctrl.prtvrbl, pihm->croptbl, pihm->elem, pihm->river, &pihm->print);
#elif defined(_RT_)
    MapOutput(outputdir, pihm->ctrl.prtvrbl, pihm->chemtbl, &pihm->rttbl, pihm->elem, pihm->river, &pihm->print);
#else
    MapOutput(outputdir, pihm->ctrl.prtvrbl, pihm->elem, pihm->river, &pihm->print);
#endif

    // Backup input files
#if !defined(_MSC_VER)
    if (!append_mode)
    {
        BackupInput(outputdir, &pihm->filename);
    }
#endif

#if defined(_LUMPED_) && defined(_RT_)
    InitOutputFiles(outputdir, pihm->ctrl.waterbal, pihm->ctrl.ascii, pihm->chemtbl, &pihm->rttbl, &pihm->print);
#else
    InitOutputFiles(outputdir, pihm->ctrl.waterbal, pihm->ctrl.ascii, &pihm->print);
#endif

    pihm_printf(VL_VERBOSE, "\n\nSolving ODE system ... \n\n");

    // Set solver parameters
    SetCVodeParam(pihm, cvode_mem, &sun_ls, CV_Y);

#if defined(_BGC_)
    first_balance = 1;
#endif

    // Run PIHM
#if defined(_OPENMP)
    start_omp = omp_get_wtime();
#else
    start = clock();
#endif

    ctrl = &pihm->ctrl;

    if (spinup_mode)
    {
        Spinup(pihm, CV_Y, cvode_mem, &sun_ls);

        // In spin-up mode, initial conditions are always printed
        PrintInit(outputdir, ctrl->endtime, ctrl->starttime, ctrl->endtime, ctrl->prtvrbl[IC_CTRL], pihm->elem,
            pihm->river);

#if defined(_BGC_)
        WriteBgcIc(outputdir, pihm->elem, pihm->river);
#endif

#if defined(_CYCLES_)
        WriteCyclesIc(outputdir, pihm->elem);
#endif

#if defined(_RT_)
        WriteRtIc(outputdir, pihm->chemtbl, &pihm->rttbl, pihm->elem);
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

            // Run PIHM time step
            PIHM(cputime, pihm, cvode_mem, CV_Y);

            // Adjust CVODE max step to reduce oscillation
            AdjCVodeMaxStep(cvode_mem, &pihm->ctrl);

            // Print CVODE performance and statistics
            if (debug_mode)
            {
                PrintPerf(ctrl->tout[ctrl->cstep + 1], ctrl->starttime, cputime_dt, cputime, ctrl->maxstep,
                    pihm->print.cvodeperf_file, cvode_mem);
            }

            // Write init files
            if (ctrl->write_ic)
            {
                PrintInit(outputdir, ctrl->tout[ctrl->cstep + 1], ctrl->starttime, ctrl->endtime,
                    ctrl->prtvrbl[IC_CTRL], pihm->elem, pihm->river);
            }

        }

#if defined(_BGC_)
        if (ctrl->write_bgc_restart)
        {
            WriteBgcIc(outputdir, pihm->elem, pihm->river);
        }
#endif

#if defined(_CYCLES_)
        if (ctrl->write_cycles_restart)
        {
            WriteCyclesIc(outputdir, pihm->elem);
        }
#endif

#if defined(_RT_)
        if (ctrl->write_rt_restart)
        {
            WriteRtIc(outputdir, pihm->chemtbl, &pihm->rttbl, pihm->elem);
        }
#endif
    }

    if (debug_mode)
    {
        PrintCVodeFinalStats(cvode_mem);
    }

    // Free memory
    N_VDestroy(CV_Y);

    // Free integrator memory
    CVodeFree(&cvode_mem);
    SUNLinSolFree(sun_ls);
    FreeMem(pihm);
    free(pihm);

    pihm_printf(VL_BRIEF, "Simulation completed.\n");

    return EXIT_SUCCESS;
}
