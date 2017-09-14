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
    double          ptime_omp, start_omp, ct_omp;
#else
    clock_t         ptime, start, ct;
#endif
    double          cputime, cputime_dt;    /* Time cpu duration */
    static double   dtime = 0.0;
    long int        nst;
    long int        nfe;
    long int        nni;
    long int        ncfn;
    long int        netf;
    long int        ncfni = 0;
    long int        nnii = 10;
    int             flag;
    char            watbal_fn[MAXSTRING];
    char            perf_fn[MAXSTRING];
    char            conv_fn[MAXSTRING];
    double          maxstep;

#ifdef _OPENMP
    /* Set the number of threads to use */
    nthreads = omp_get_max_threads();
#endif

    memset(outputdir, 0, MAXSTRING);

    /* Print AscII art */
    AsciiArt();

    /* Read command line arguments */
    ParseCmdLineParam(argc, argv, outputdir);

    /* Create output directory */
    CreateOutputDir(outputdir);

    /* Allocate memory for model data structure */
    pihm = (pihm_struct)malloc(sizeof(*pihm));

    /* Read PIHM input files */
    ReadAlloc(pihm);

    /* Initialize CVode state variables */
    CV_Y = N_VNew(NSV);

    /* Initialize PIHM structure */
    Initialize(pihm, CV_Y);

    /* Allocate memory for solver */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (cvode_mem == NULL)
    {
        PIHMprintf(VL_ERROR, "Error in allocating memory for solver.\n");
        PIHMexit(EXIT_FAILURE);
    }

    /* Create output structures */
    MapOutput(pihm, outputdir);
    MapTecplotOutput(pihm, outputdir);

    /* Backup input files */
#if !defined(_MSC_VER)
    BKInput(outputdir);
#endif

    /* Initialize output files and structures */
    if (pihm->ctrl.waterB)
    {
        sprintf(watbal_fn, "%s%s_WaterBalance.plt", outputdir, project);
        pihm->print.walbal_file = fopen(watbal_fn, "w");
        CheckFile(pihm->print.walbal_file, watbal_fn);
    }

    if (debug_mode)
    {
        sprintf(perf_fn, "%s%s.perf.txt", outputdir, project);
        pihm->print.cvodeperf_file = fopen(perf_fn, "w");
        CheckFile(pihm->print.cvodeperf_file, perf_fn);

        fprintf(pihm->print.cvodeperf_file,
            " Time step, cpu_dt, cpu_time, solver_step\n");

        sprintf(conv_fn, "%s%s.CVODE.log", outputdir, project);
        pihm->print.cvodeconv_file = fopen(conv_fn, "w");
        CheckFile(pihm->print.cvodeconv_file, conv_fn);
    }

    InitOutputFile(&pihm->print, pihm->ctrl.ascii);

    if (pihm->ctrl.write_ic)
    {
        char            icdir[MAXSTRING];
        sprintf(icdir, "input/%s/ic/", project);
        pihm_mkdir(icdir);
    }

    PIHMprintf(VL_VERBOSE, "\n\nSolving ODE system ... \n\n");

    /* Set solver parameters */
    SetCVodeParam(pihm, cvode_mem, CV_Y);
    maxstep = pihm->ctrl.maxstep;

#if defined(_BGC_) || defined (_CYCLES_)
    first_balance = 1;
#endif

    /*
     * Run PIHM
     */
#ifdef _BGC_
    if (spinup_mode)
    {
        BgcSpinup(pihm, CV_Y, cvode_mem);
    }
    else
    {
#endif
#ifdef _OPENMP
        start_omp = omp_get_wtime();
        ptime_omp = start_omp;
#else
        start = clock();
        ptime = start;
#endif
        for (i = 0; i < pihm->ctrl.nstep; i++)
        {
#ifdef _OPENMP
            ct_omp = omp_get_wtime();
            cputime_dt = ((double)(ct_omp - ptime_omp));
            cputime = ((double)(ct_omp - start_omp));
            ptime_omp = ct_omp;
#else
            ct = clock();
            cputime_dt = ((double)(ct - ptime)) / CLOCKS_PER_SEC;
            cputime = ((double)(ct - start)) / CLOCKS_PER_SEC;
            ptime = ct;
#endif
            ncfni = ncfn;
            nnii = nni;
            PIHM(pihm, cvode_mem, CV_Y, pihm->ctrl.tout[i],
                pihm->ctrl.tout[i + 1], cputime);

            if (debug_mode)
            {
                dtime += cputime_dt;

                if (pihm->ctrl.tout[i] % 3600 == 0)
                {
                    fprintf(pihm->print.cvodeperf_file, "%d %f %f %f\n",
                        pihm->ctrl.tout[i] - pihm->ctrl.starttime, cputime_dt,
                        cputime, maxstep);
                    dtime = 0.0;
                }
                /* Print CVODE statistics */
                PrintStats(cvode_mem, pihm->print.cvodeconv_file);
            }

            flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
            flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);

            /* Variable CVODE max step (to reduce oscillations) */
            if ((ncfn - ncfni > pihm->ctrl.nncfn ||
                nni - nnii > pihm->ctrl.nnimax) &&
                maxstep > pihm->ctrl.stmin)
            {
                maxstep /= pihm->ctrl.decr;
                flag = CVodeSetMaxStep(cvode_mem, (realtype)maxstep);
            }
            if (ncfn == ncfni && maxstep < pihm->ctrl.maxstep &&
                nni - nnii < pihm->ctrl.nnimin)
            {
                maxstep *= pihm->ctrl.incr;
                flag = CVodeSetMaxStep(cvode_mem, (realtype)maxstep);
            }

            /*
             * Write init files
             */
            if (pihm->ctrl.write_ic &&
                (pihm->ctrl.tout[i] - pihm->ctrl.starttime) %
                pihm->ctrl.prtvrbl[IC_CTRL] == 0)
            {
                PrtInit(pihm->elem, pihm->riv, pihm->ctrl.tout[i]);
            }
        }
#ifdef _BGC_
    }
#endif


#ifdef _BGC_
    if (pihm->ctrl.write_bgc_restart)
    {
        WriteBgcIC(pihm->filename.bgcic, pihm->elem, pihm->riv);
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

        printf("nst = %-6ld nfe  = %-6ld \n", nst, nfe);
        printf("nni = %-6ld ncfn = %-6ld netf = %-6ld\n \n", nni, ncfn, netf);
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
