#include "pihm.h"

int             verbose_mode;
int             debug_mode;
int             corr_mode;
char            project[MAXSTRING];
int             nelem;
int             nriver;
#ifdef _OPENMP
int             nthreads;
#endif
#ifdef _DAILY_
int             first_day;
#endif

int main (int argc, char *argv[])
{
    char            outputdir[MAXSTRING];
    pihm_struct     pihm;
    N_Vector        CV_Y;       /* State Variables Vector */
    void           *cvode_mem;  /* Model Data Pointer */
    int             i;          /* Loop index */

#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif

    /* Read command line arguments */
    ParseCmdLineParam (argc, argv, outputdir);

    /* Print AscII art */
    AsciiArt ();

    /* Create output directory */
    CreateOutputDir (outputdir);

    /* Allocate memory for model data structure */
    pihm = (pihm_struct)malloc (sizeof *pihm);

    /* Read PIHM input files */
    ReadAlloc (project, pihm);

    /* Initialize CVode state variables */
    CV_Y = N_VNew (NSV);

    /* Initialize PIHM structure */
    Initialize (pihm, CV_Y);

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
    BKInput (project, outputdir);

    /* Initialize output files and structures */
    InitOutputFile (pihm->prtctrl, pihm->ctrl.nprint, pihm->ctrl.ascii);

    PIHMprintf (VL_VERBOSE, "\n\nSolving ODE system ... \n\n");

    /* Set solver parameters */
    SetCVodeParam (pihm, cvode_mem, CV_Y);

    /*
     * Run PIHM
     */
#ifdef _DAILY_
    first_day = 1;
#endif
    for (i = 0; i < pihm->ctrl.nstep; i++)
    {
        PIHM (pihm, cvode_mem, CV_Y,
            pihm->ctrl.tout[i], pihm->ctrl.tout[i + 1]);

        /* Print outputs */
        PrintData (pihm->prtctrl, pihm->ctrl.nprint, pihm->ctrl.tout[i + 1],
            pihm->ctrl.tout[i + 1] - pihm->ctrl.starttime,
            pihm->ctrl.stepsize, pihm->ctrl.ascii);
    }

#ifdef _BGC_
    if (pihm->ctrl.bgc_spinup)
    {
        //BGCSpinup (pihm, outputdir);
    }
#endif

    /*
     * Write init files
     */
    if (pihm->ctrl.write_ic)
    {
        PrtInit (pihm->elem, pihm->riv, project);

#ifdef _BGC_
        if (pihm->ctrl.bgc_spinup)
        {
            WriteBgcIC (pihm->filename.bgcic, pihm->elem, pihm->riv);
        }
#endif
    }

    /* Free memory */
    N_VDestroy (CV_Y);

    /* Free integrator memory */
    CVodeFree (&cvode_mem);

    FreeData (pihm);
    free (pihm);
    PIHMprintf (VL_NORMAL, "\nSimulation completed.\n");

    return (EXIT_SUCCESS);
}
