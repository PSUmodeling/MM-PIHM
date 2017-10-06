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
    CV_Y = N_VNew(NumStateVar());

    /* Initialize PIHM structure */
    Initialize(pihm, CV_Y, &cvode_mem);

    /* Create output directory */
    CreateOutputDir(outputdir);

    /* Create output structures */
    MapOutput(pihm->ctrl.prtvrbl, pihm->ctrl.tpprtvrbl, pihm->elem, pihm->river,
        &pihm->meshtbl, outputdir, &pihm->print);

    /* Backup input files */
#if !defined(_MSC_VER)
    BackupInput(outputdir);
#endif

    InitOutputFile(&pihm->print, outputdir, pihm->ctrl.waterbal,
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

            /* Print CVODE performance and statistics */
            if (debug_mode)
            {
                PrintPerf(pihm->ctrl.tout[i + 1], pihm->ctrl.starttime,
                    cputime_dt, cputime, pihm->ctrl.maxstep,
                    pihm->print.cvodeperf_file);

                PrintStats(cvode_mem, pihm->print.cvodeconv_file);
            }

            /* Write init files */
            if (pihm->ctrl.write_ic)
            {
                PrintInit(pihm->elem, pihm->river, outputdir,
                    pihm->ctrl.tout[i + 1], pihm->ctrl.starttime,
                    pihm->ctrl.endtime, pihm->ctrl.prtvrbl[IC_CTRL]);
            }
        }
#ifdef _BGC_
    }
#endif

#ifdef _BGC_
    if (pihm->ctrl.write_bgc_restart)
    {
        WriteBgcIc(outputdir, pihm->elem, pihm->river);
    }
#endif

#ifdef _CYCLES_
    if (pihm->ctrl.write_cycles_restart)
    {
        WriteCyclesIC(pihm->filename.cyclesic, pihm->elem, pihm->river);
    }
#endif

    if (debug_mode)
    {
        PrintCVodeFinalStats(cvode_mem);
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
