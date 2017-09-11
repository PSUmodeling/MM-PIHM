#include "pihm.h"

int             verbose_mode;
int             debug_mode;
int             corr_mode;
int             spinup_mode;
char            project[MAXSTRING];
int             nelem;
int             nriver;
clock_t         ptime, start, ct;
realtype        cputime, cputime_dt;    /* Time cpu duration */
static double   dtime = 0;
long int        nst, nfe, nfeLS, nni, ncfn, netf, ncfni = 0, nnii = 10; /*Variables for monitoring performance */
int             flag;
char            WBname[100];
char            Perfname[50];
char            Convname[50];
double          maxstep;
FILE           *WaterBalance;   /* Water balance file */
FILE           *Perf;           /* Performance file */
FILE           *Conv;           /* CVODE convergence file */
#ifdef _OPENMP
realtype        ptime_omp, start_omp, ct_omp;
#endif
#ifdef _OPENMP
int             nthreads = 1;   /* default value */
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

    memset(outputdir, 0, MAXSTRING);

    /* Set the number of threads to use */
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif

    /* Read command line arguments */
    ParseCmdLineParam(argc, argv, outputdir);

    /* Print AscII art */
    AsciiArt();

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
        WaterBalance = fopen(WBname, "w");
        CheckFile(WaterBalance, WBname);
    }
    if (pihm->ctrl.cvode_perf)
    {
        sprintf(Perfname, "%s%s_Performance.txt", outputdir, project);
        Perf = fopen(Perfname, "w");
        CheckFile(Perf, Perfname);
        dtime = dtime + cputime_dt;
        fprintf(Perf, " Time step, cpu_dt, cpu_time, solver_step\n");
        sprintf(Convname, "%s%s_CVODE.log", outputdir, project);
        Conv = fopen(Convname, "w");
        CheckFile(Conv, Convname);
    }

    InitOutputFile(pihm->prtctrl, pihm->ctrl.nprint, pihm->ctrl.ascii,
        pihm->prtctrlT, pihm->ctrl.nprintT, pihm->ctrl.tecplot);

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
                pihm->ctrl.tout[i + 1], cputime, WaterBalance);

            if (pihm->ctrl.cvode_perf)
            {
                dtime = dtime + cputime_dt;
                if (pihm->ctrl.tout[i] % 3600 == 0)
                {
                    fprintf(Perf, "%d %f %f %f\n",
                        (pihm->ctrl.tout[i] - pihm->ctrl.starttime), cputime_dt,
                        cputime, maxstep);
                    dtime = 0.0;
                }
                /* Print CVODE statistics */
                PrintStats(cvode_mem, Conv);
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
        fclose(WaterBalance);
    }
    if (pihm->ctrl.cvode_perf)
    {
        fclose(Perf);
        fclose(Conv);
    }
    FreeData(pihm);
    free(pihm);
    PIHMprintf(VL_NORMAL, "\nSimulation completed.\n");

    return EXIT_SUCCESS;
}
