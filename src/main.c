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
#endif
#if defined(_BGC_)
int             first_balance;
#endif

int main(int argc, char *argv[])
{
    char            outputdir[MAXSTRING];
#if defined(_OPENMP)
    double          start_omp;
#else
    clock_t         start;
#endif
    double          cputime, cputime_dt;    // Time cpu duration
    ctrl_struct    *ctrl;
    pihm_struct    *pihm = (pihm_struct *)malloc(sizeof(pihm_struct));
    cvode_struct   *cvode = (cvode_struct *)malloc(sizeof(cvode_struct));

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

    // Read PIHM input files
    ReadAlloc(pihm);

    // Initialize PIHM structure
    Initialize(pihm, cvode);

    // Create output directory
    CreateOutputDir(outputdir);

    WriteMetadata(outputdir);

    // Create output structures
#if defined(_CYCLES_)
    MapOutput(outputdir, pihm->ctrl.prtvrbl, pihm->croptbl, pihm->elem, pihm->river, &pihm->print);
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

    InitOutputFiles(outputdir, pihm->ctrl.waterbal, pihm->ctrl.ascii, &pihm->print);

    pihm_printf(VL_VERBOSE, "\n\nSolving ODE system ... \n\n");

    // Set solver parameters
    SetCVodeParam(pihm, cvode);

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
        Spinup(pihm, cvode);

        // In spin-up mode, initial conditions are always printed
        PrintInit(outputdir, ctrl->endtime, ctrl->starttime, ctrl->endtime, ctrl->prtvrbl[IC_CTRL], pihm->elem, pihm->river);

#if defined(_BGC_)
        WriteBgcIc(outputdir, pihm->elem, pihm->river);
#endif

#if defined(_CYCLES_)
        WriteCyclesIc(outputdir, pihm->elem);
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
            PIHM(cputime, pihm, cvode);

            // Adjust CVODE max step to reduce oscillation
            AdjCVodeMaxStep(cvode, &pihm->ctrl);

            // Print CVODE performance and statistics
            if (debug_mode)
            {
                PrintPerf(ctrl->tout[ctrl->cstep + 1], ctrl->starttime, cputime_dt, cputime, ctrl->maxstep, pihm->print.cvodeperf_file, cvode);
            }

            // Write init files
            if (ctrl->write_ic)
            {
                PrintInit(outputdir, ctrl->tout[ctrl->cstep + 1], ctrl->starttime, ctrl->endtime, ctrl->prtvrbl[IC_CTRL], pihm->elem, pihm->river);
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
    }

    if (debug_mode)
    {
        PrintCVodeFinalStats(cvode);
    }

    // Free memory
    N_VDestroy(cvode->CV_Y);

    // Free integrator memory
    CVodeFree(&cvode->memory_block);
    SUNContext_Free(&cvode->sunctx);
    FreeMem(pihm);
    free(pihm);

    pihm_printf(VL_BRIEF, "Simulation completed.\n");

    return EXIT_SUCCESS;
}
