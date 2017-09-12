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
    N_Vector        CV_Y, abstol;
    void           *cvode_mem;
    int             i;
    clock_t         ptime, start, ct;
    realtype        cputime, cputime_dt;    /* Time cpu duration */
    static double   dtime = 0;
    long int        nst;
    long int        nfe;
    long int        nfeLS;
    long int        nni;
    long int        ncfn;
    long int        netf;
    long int        ncfni = 0;
    long int        nnii = 10;
    int             flag;
    char            WBname[MAXSTRING];
    char            Perfname[MAXSTRING];
    char            Convname[MAXSTRING];
    double          maxstep;
#ifdef _OPENMP
    realtype        ptime_omp, start_omp, ct_omp;
#endif

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
    ReadAlloc(project, pihm);

    /* Initialize CVode state variables */
    CV_Y = N_VNew(NSV);
    abstol = N_VNew(NSV);

    /* Initialize PIHM structure */
    Initialize(pihm, CV_Y);

    /* Set the vector absolute tolerance */
    for (i = 0; i < NSV; i++)
    {
        NV_Ith(abstol, i) = pihm->ctrl.abstol;
    }

    /* Allocate memory for solver */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (cvode_mem == NULL)
    {
        PIHMprintf(VL_ERROR, "Error in allocating memory for solver.\n");
        PIHMexit(EXIT_FAILURE);
    }

    /* Create output structures */
    MapOutput(project, pihm, outputdir);

    /* Backup input files */
#if !defined(_MSC_VER)
    BKInput(project, outputdir);
#endif

    /* Initialize output files and structures */
    if (pihm->ctrl.waterB)
    {
        sprintf(WBname, "%s%s_WaterBalance.plt", outputdir, project);
        pihm->print.walbal_file = fopen(WBname, "w");
        CheckFile(pihm->print.walbal_file, WBname);
    }
    if (pihm->ctrl.cvode_perf)
    {
        sprintf(Perfname, "%s%s_Performance.txt", outputdir, project);
        pihm->print.cvodeperf_file = fopen(Perfname, "w");
        CheckFile(pihm->print.cvodeperf_file, Perfname);
        dtime = dtime + cputime_dt;
        fprintf(pihm->print.cvodeperf_file, " Time step, cpu_dt, cpu_time, solver_step\n");
        sprintf(Convname, "%s%s_CVODE.log", outputdir, project);
        pihm->print.cvodeconv_file = fopen(Convname, "w");
        CheckFile(pihm->print.cvodeconv_file, Convname);
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

            if (pihm->ctrl.cvode_perf)
            {
                dtime = dtime + cputime_dt;
                if (pihm->ctrl.tout[i] % 3600 == 0)
                {
                    fprintf(pihm->print.cvodeperf_file, "%d %f %f %f\n",
                        (pihm->ctrl.tout[i] - pihm->ctrl.starttime), cputime_dt,
                        cputime, maxstep);
                    dtime = 0.0;
                }
                /* Print CVODE statistics */
                PrintStats(cvode_mem, pihm->print.cvodeconv_file);
            }
            flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
            flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);

            /* Variable CVode Max Step (to reduce oscillations) */
            if ((ncfn - ncfni > pihm->ctrl.nncfn ||
                nni - nnii > pihm->ctrl.nnimax) &&
                maxstep > pihm->ctrl.stmin)
            {
                maxstep = maxstep / pihm->ctrl.decr;
                flag = CVodeSetMaxStep(cvode_mem, (realtype)maxstep);
            }
            if (ncfn == ncfni && maxstep < pihm->ctrl.maxstep &&
                nni - nnii < pihm->ctrl.nnimin)
            {
                maxstep = maxstep * pihm->ctrl.incr;
                flag = CVodeSetMaxStep(cvode_mem, (realtype)maxstep);
            }

            /*
             * Write init files
             */
            if (pihm->ctrl.write_ic &&
                (pihm->ctrl.tout[i] - pihm->ctrl.starttime) %
                pihm->ctrl.prtvrbl[IC_CTRL] == 0)
            {
                PrtInit(pihm->elem, pihm->riv, project, pihm->ctrl.tout[i]);
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

    /* Print some CVODE statistics */
    flag = CVodeGetNumSteps(cvode_mem, &nst);
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    flag = CVodeGetNumErrTestFails(cvode_mem, &netf);

    printf("nst = %-6ld nfe  = %-6ld \n", nst, nfe);
    printf("nni = %-6ld ncfn = %-6ld netf = %-6ld\n \n", nni, ncfn, netf);

    /* Free memory */
    N_VDestroy(CV_Y);
    N_VDestroy(abstol);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);
    if (pihm->ctrl.waterB)
    {
        fclose(pihm->print.walbal_file);
    }
    if (pihm->ctrl.cvode_perf)
    {
        fclose(pihm->print.cvodeperf_file);
        fclose(pihm->print.cvodeconv_file);
    }
    FreeData(pihm);
    free(pihm);
    PIHMprintf(VL_NORMAL, "\nSimulation completed.\n");

    return EXIT_SUCCESS;
}
