#include "pihm.h"

int             verbose_mode;
int             debug_mode;
int             corr_mode;
int             spinup_mode;
int             num_threads;
char            project[MAXSTRING];
int             nelem;
int             nriver;
clock_t         ptime, start, ct;
/* Time cpu duration */
#ifdef _OPENMP
	realtype	ptime_omp, start_omp, ct_omp;
#else
    realtype    cputime, cputime_dt;
#endif
#ifdef _OPENMP
int             nthreads;
#endif
#if defined(_BGC_) || defined (_CYCLES_)
int             first_balance;
#endif

int main(int argc, char *argv[])
{
	char            outputdir[MAXSTRING];
	int             spec_output_mode = 0;
	pihm_struct     pihm;
	N_Vector        CV_Y, abstol;       /* State Variables Vector */
	void           *cvode_mem;  /* Model Data Pointer */
	int             i;          /* Loop index */

	/* Set the number of threads to use */
	num_threads = 1;     /* default value */
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif

    /* Read command line arguments */
    ParseCmdLineParam (argc, argv, &spec_output_mode, outputdir);

    /* Print AscII art */
    AsciiArt ();

    /* Create output directory */
    CreateOutputDir (outputdir);

    /* Allocate memory for model data structure */
    pihm = (pihm_struct)malloc (sizeof *pihm);

    /* Read PIHM input files */
    ReadAlloc (project, pihm);

/* Initialize CVode state variables */
#ifdef _OPENMP
	CV_Y = N_VNew_OpenMP(NSV, num_threads);
	abstol = N_VNew_OpenMP(NSV, num_threads);
#else
	CV_Y = N_VNew_Serial(NSV);
	abstol = N_VNew_Serial(NSV);
#endif

    /* Initialize PIHM structure */
    Initialize (pihm, CV_Y);
	/* Set the vector absolute tolerance */
	for (i = 0; i < NSV; i++) {
#ifdef _OPENMP
		NV_Ith_OMP(abstol, i) = pihm->ctrl.abstol;
#else
		NV_Ith_S(abstol, i) = pihm->ctrl.abstol;
#endif
	}
    /* Allocate memory for solver */
    cvode_mem = CVodeCreate (CV_BDF, CV_NEWTON);
    if (cvode_mem == NULL)
    {
        PIHMprintf (VL_ERROR, "Error in allocating memory for solver.\n");
        PIHMexit (EXIT_FAILURE);
    }

    /* Create output structures */
    MapOutput (project, pihm, outputdir);

    /* Backup input files */
#if !defined(_MSC_VER)
    BKInput (project, outputdir);
#endif

    /* Initialize output files and structures */
    InitOutputFile (pihm->prtctrl, pihm->ctrl.nprint, pihm->ctrl.ascii);

    PIHMprintf (VL_VERBOSE, "\n\nSolving ODE system ... \n\n");

    /* Set solver parameters */
    SetCVodeParam (pihm, cvode_mem, CV_Y, abstol);

#if defined(_BGC_) || defined (_CYCLES_)
    first_balance = 1;
#endif

    /*
     * Run PIHM
     */
#ifdef _BGC_
    if (spinup_mode)
    {
        BgcSpinup (pihm, CV_Y, cvode_mem);
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
        PIHM (pihm, cvode_mem, CV_Y, abstol, pihm->ctrl.tout[i],
            pihm->ctrl.tout[i + 1], outputdir, project, cputime_dt, cputime);
    }
#ifdef _BGC_
    }
#endif

    /*
     * Write init files
     */
    if (pihm->ctrl.write_ic)
    {
        PrtInit (pihm->elem, pihm->riv, project);
    }

#ifdef _BGC_
    if (pihm->ctrl.write_bgc_restart)
    {
        WriteBgcIC (pihm->filename.bgcic, pihm->elem, pihm->riv);
    }
#endif

#ifdef _CYCLES_
    if (pihm->ctrl.write_cycles_restart)
    {
        WriteCyclesIC (pihm->filename.cyclesic, pihm->elem, pihm->riv);
    }
#endif

    /* Free memory */
#ifdef _OPENMP
	N_VDestroy_OpenMP(CV_Y);
	N_VDestroy_OpenMP(abstol);
#else
	N_VDestroy_Serial(CV_Y);
	N_VDestroy_Serial(abstol);
#endif
    /* Free integrator memory */
    CVodeFree (&cvode_mem);

    FreeData (pihm);
    free (pihm);
    PIHMprintf (VL_NORMAL, "\nSimulation completed.\n");

    return (EXIT_SUCCESS);
 }
