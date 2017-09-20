#include "pihm.h"

/* Global variables */
int             verbose_mode;
int             debug_mode;
int             corr_mode;
int             spinup_mode;
int             tecplot;
char            project[MAXSTRING];
int             nelem;
int             nriver;
#ifdef _OPENMP
int             nthreads = 1;    /* Default value */
#endif
#if defined(_BGC_) || defined (_CYCLES_)
int             first_balance;
#endif

int main(int argc, char *argv[])
{
    char            outputdir[MAXSTRING];
    pihm_struct     pihm;
    N_Vector        CV_Y;
    void           *cvode_mem;
    int             i;
#ifdef _OPENMP
    double          start_omp;
#else
    clock_t         start;
#endif
    double          cputime, cputime_dt;    /* Time cpu duration */
    long int        nst;
    long int        nfe;
    long int        nni;
    long int        ncfn;
    long int        netf;
    int             flag;

#ifdef _OPENMP
    /* Set the number of threads to use */
    nthreads = omp_get_max_threads();
#endif

    memset(outputdir, 0, MAXSTRING);

    /* Print AscII art */
    AsciiArt();

    /* Read command line arguments */
    ParseCmdLineParam(argc, argv, outputdir);

    /* Allocate memory for model data structure */
    pihm = (pihm_struct)malloc(sizeof(*pihm));

    /* Read PIHM input files */
    ReadAlloc(pihm);

    /* Initialize CVode state variables */
    CV_Y = N_VNew(NSV);

    /* Initialize PIHM structure */
    Initialize(pihm, CV_Y, &cvode_mem);

    /* Create output directory */
    CreateOutputDir(outputdir);

    /* Create output structures */
    MapOutput(pihm, outputdir);

    /* Backup input files */
#if !defined(_MSC_VER)
    BKInput(outputdir);
#endif

    InitOutputFile(&pihm->print, outputdir, pihm->ctrl.waterB,
        pihm->ctrl.ascii);

    PIHMprintf(VL_VERBOSE, "\n\nSolving ODE system ... \n\n");

    /* Set solver parameters */
    SetCVodeParam(pihm, cvode_mem, CV_Y);

#if defined(_BGC_) || defined (_CYCLES_)
    first_balance = 1;
#endif

    /*
     * Run PIHM
     */
#ifdef _OPENMP
        start_omp = omp_get_wtime();
#else
        start = clock();
#endif

#ifdef _BGC_
    if (spinup_mode)
    {
        BgcSpinup(pihm, CV_Y, cvode_mem);
    }
    else
    {
#endif
        for (i = 0; i < pihm->ctrl.nstep; i++)
        {
#ifdef _OPENMP
            RunTime(start_omp, &cputime, &cputime_dt);
#else
            RunTime(start, &cputime, &cputime_dt);
#endif

            PIHM(pihm, cvode_mem, CV_Y, pihm->ctrl.tout[i],
                pihm->ctrl.tout[i + 1], cputime);

            /* Adjust CVODE max step to reduce oscillation */
            AdjCVodeMaxStep(cvode_mem, &pihm->ctrl);

            if (debug_mode)
            {
                if (pihm->ctrl.tout[i] % 3600 == 0)
                {
                    fprintf(pihm->print.cvodeperf_file, "%d %f %f %f\n",
                        pihm->ctrl.tout[i] - pihm->ctrl.starttime, cputime_dt,
                        cputime, pihm->ctrl.maxstep);
                    fflush(pihm->print.cvodeperf_file);
                }
                /* Print CVODE statistics */
                PrintStats(cvode_mem, pihm->print.cvodeconv_file);
            }

            /*
             * Write init files
             */
            if (pihm->ctrl.write_ic &&
                ((pihm->ctrl.tout[i] - pihm->ctrl.starttime) %
                pihm->ctrl.prtvrbl[IC_CTRL] == 0 ||
                i == pihm->ctrl.nstep - 1))
            {
                PrtInit(pihm->elem, pihm->riv, outputdir, pihm->ctrl.tout[i]);
            }
        }
#ifdef _BGC_
    }
#endif


#ifdef _BGC_
    if (pihm->ctrl.write_bgc_restart)
    {
        WriteBgcIC(outputdir, pihm->elem, pihm->riv);
    }
#endif

#ifdef _CYCLES_
    if (pihm->ctrl.write_cycles_restart)
    {
        WriteCyclesIC(pihm->filename.cyclesic, pihm->elem, pihm->riv);
    }
#endif

    if (debug_mode)
    {
        /* Print some CVODE statistics */
        flag = CVodeGetNumSteps(cvode_mem, &nst);
        flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
        flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
        flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
        flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);

        PIHMprintf(VL_NORMAL,
            "number of steps = %-6ld "
            "number of rhs evals = %-6ld\n",
            nst, nfe);
        PIHMprintf(VL_NORMAL,
            "number of nonlin solv iters = %-6ld "
            "number of nonlin solv conv fails = %-6ld "
            "number of err test fails = %-6ld\n",
            nni, ncfn, netf);
    }

    /* Free memory */
    N_VDestroy(CV_Y);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);
    FreeData(pihm);
    free(pihm);
    PIHMprintf(VL_NORMAL, "\nSimulation completed.\n");

    return EXIT_SUCCESS;
}
