#include "pihm.h"

#ifdef _NOAH_
#include "noah.h"
#endif

#ifdef _RT_
#include "rt.h"
#endif

#ifdef _BGC_
#include "bgc.h"
#endif

#ifdef _CYCLES_
#include "cycles.h"
#endif

#ifdef _ENKF_
#include "enkf.h"
#endif

int                 verbose_mode;
int                 debug_mode;

int main (int argc, char *argv[])
{
    char            project[MAXSTRING];
    char            simulation[MAXSTRING];
    char            outputdir[MAXSTRING];
    int             overwrite_mode = 0;
    int             c;
    int             first_cycle;
#ifdef _ENKF_
    int             ii;
    int             id;
    int             ierr;
    int             p;
    int             job_per_node;
    int             startmode;
    int             starttime;
    int             endtime;
    double         *param;
    int             success;
    enkf_struct     ens;

    ierr = MPI_Init (&argc, &argv);
    ierr = MPI_Comm_rank (MPI_COMM_WORLD, &id);
    ierr = MPI_Comm_size (MPI_COMM_WORLD, &p);

    //{
    //    int ii = 0;
    //    char hostname[256];
    //    gethostname(hostname, sizeof(hostname));
    //    printf("PID %d (%d) on %s ready for attach\n", getpid(), id, hostname);
    //    fflush(stdout);
    //    while (0 == ii)
    //        sleep(5);
    //}

    if (id == 0)
    {
#endif
    printf ("\n");
    printf ("\t\t########  #### ##     ## ##     ##\n");
    printf ("\t\t##     ##  ##  ##     ## ###   ###\n");
    printf ("\t\t##     ##  ##  ##     ## #### ####\n");
    printf ("\t\t########   ##  ######### ## ### ##\n");
    printf ("\t\t##         ##  ##     ## ##     ##\n");
    printf ("\t\t##         ##  ##     ## ##     ##\n");
    printf ("\t\t##        #### ##     ## ##     ##\n");
    printf ("\n\t    The Penn State Integrated Hydrologic Model\n");

#ifdef _NOAH_
    printf ("\n\t    * Land surface module turned on.\n");
#endif
#ifdef _RT_
    printf ("\n\t    * Reactive transport land surface hydrological mode.\n");
#endif
#ifdef _BGC_
    printf ("\n\t    * Biogeochemistry module turned on.\n");
#endif
#ifdef _ENKF_
    printf ("\n\t    * Ensemble Kalman filter turned on.\n");
#endif
#ifdef _CYCLES_
    printf ("\n\t    * Crop module turned on.\n");
#endif
#ifdef _ENKF_
    }
#endif

    /*
     * Read command line arguments
     */
    while ((c = getopt (argc, argv, "odv")) != -1)
    {
        if (optind >= argc)
        {
            printf ("\nUsage: ./pihm [-o -d -v] <project name>\n");
            printf
                ("\t-o Output will be written in the \"output\" directory and overwrite existing output.\n");
            printf ("\t-v Verbose mode\n");
            printf ("\t-d Debug mode\n");
            PihmExit (1);
        }
        switch (c)
        {
            case 'o':
                overwrite_mode = 1;
                printf
                    ("Overwrite mode turned on. Output directory is \"./output\".\n");
                break;
            case 'd':
                debug_mode = 1;
                printf ("Debug mode turned on.\n");
                break;
            case 'v':
                verbose_mode = 1;
                printf ("Verbose mode turned on.\n");
                break;
            case '?':
                printf ("Option not recognisable %s\n", argv[optind]);
                break;
            default:
                break;
        }
    }

    if (optind >= argc)
    {
        printf ("\nERROR:You must specify the name of project!\n");
        printf ("Usage: ./pihm [-o -d -v] <project name>\n");
        printf ("\t-o Output will be written in the \"output\" directory"
            "and overwrite existing output.\n");
        printf ("\t-v Verbose mode\n");
        printf ("\t-d Debug mode\n");
        PihmExit (1);
    }
    else
    {
        strcpy (project, argv[optind]);
    }

    first_cycle = 1;

#ifdef _ENKF_
    if (id == 0)
    {
#endif
    /*
     * Create output directory
     */
    CreateOutputDir (project, outputdir, overwrite_mode);
#ifdef _ENKF_
    }
#endif

#ifdef _ENKF_
    if (id == 0)
    {
        /*
         * EnKF initialization
         */
        ens = (enkf_struct) malloc (sizeof *ens);

        /* Read EnKF input file */
        EnKFRead (project, ens);
        
        /* Check if node number is appropriate */
        if (ens->ne % (p - 1) != 0)
        {
            printf ("ERROR: Please specify a correct node number!\n");
            fflush (stdout);
            PihmExit (1);
        }
        else
        {
            job_per_node = ens->ne / (p - 1);
        }
        
        /* Initialize observation operator vector */
        InitOper (project, ens);

        /* Perturb model parameters */
        Perturb (project, ens, outputdir);
    }

    /* Broadcast jobs per node and output directory to all nodes */
    MPI_Bcast (&job_per_node, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (outputdir, MAXSTRING, MPI_CHAR, 0, MPI_COMM_WORLD);
    param = (double *) malloc (job_per_node * (p - 1) * MAXPARAM * sizeof (double));

    /*
     * EnKF cycles
     */
    while (1)
    {
        if (id == 0)
        {
            if (ens->cycle_start_time >= ens->end_time)
            {
                /* Special case when EnKF cycle ends */
                JobHandout (ens->cycle_start_time, ens->cycle_end_time,
                    BADVAL, ens->member, param, ens->ne, p - 1);
            
                break;
            }
            else
            {
                /* Send required parameters to different nodes for PIHM
                 * runs */
                JobHandout (ens->cycle_start_time, ens->cycle_end_time,
                    ens->mbr_start_mode, ens->member, param, ens->ne, p - 1);
            }

            /* Screen output */
            PrintEnKFStatus (ens->cycle_start_time, ens->cycle_end_time);

            /* Waiting for different nodes to send signals indicating PIHM
             * simulations done */
            JobHandIn (p - 1);

            ens->update_param = 1;
            ens->update_var = 1;

            /*
             * Read variables from PIHM output files
             */
            ReadVar (project, outputdir, ens, ens->cycle_end_time);

            /* EnKF data assimilation */
            EnKF (project, ens, ens->cycle_end_time, outputdir);

            /* Proceed to next cycle */
            ens->mbr_start_mode = 3;
            ens->cycle_start_time = ens->cycle_end_time;
            ens->cycle_end_time += ens->interval;
        }
        else
        {
            /* Receive required parameters from Node 0 */
            JobRecv (&starttime, &endtime, &startmode, param, job_per_node * (p - 1));

            if (startmode == BADVAL)
            {
                /* Special case indicating end of EnKF cycles */
                break;
            }

            for (ii = (id - 1) * job_per_node; ii < id * job_per_node; ii++)
            {
                /* Determine name of simulation */
                sprintf (simulation, "%s.%3.3d", project, ii + 1);

                /* Run PIHM */
                PIHMRun (simulation, outputdir, first_cycle,
                    starttime, endtime, startmode, param + ii * MAXPARAM);
            }

            first_cycle = 0;
            success = 1;

            /* Notify Node 0 PIHM run is completed */
            ierr = MPI_Send (&success, 1, MPI_INT, 0, SUCCESS_TAG, MPI_COMM_WORLD);
        }
    }

    if (id == 0)
    {
        printf ("\nSimulation completed.\n");

        FreeEns (ens);
        free (ens);
    }

    free (param);
    ierr = MPI_Finalize ();
#else
    /* The name of the simulation is the same as the project */
    strcpy (simulation, project);

    PIHMRun (simulation, outputdir, first_cycle);

    printf ("\nSimulation completed.\n");
#endif

    return (0);
}

#ifdef _ENKF_
void PIHMRun (char *simulation, char *outputdir, int first_cycle,
    int starttime, int endtime, int startmode, double *param)
#else
void PIHMRun (char *simulation, char *outputdir, int first_cycle)
#endif
{
    pihm_struct     pihm;
    N_Vector        CV_Y;       /* State Variables Vector */
    void           *cvode_mem;  /* Model Data Pointer */
    int             flag;       /* Flag to test return value */
    int             nsv;        /* Problem size */
    int             i;       /* Loop index */
    int             t;          /* Simulation time */
    realtype        solvert;
    struct tm      *timestamp;
    time_t          rawtime;
    int             nextptr;
    realtype        cvode_val;

#ifdef _NOAH_
    lsm_struct      noah;
#endif

#ifdef _RT_
    Chem_Data       chData;
#endif

#ifdef _BGC_
    bgc_struct      bgc;
#endif

#ifdef _CYCLES_
    CyclesStruct    cycles;
#endif

    /* Allocate memory for model data structure */
    pihm = (pihm_struct)malloc (sizeof *pihm);
#ifdef _NOAH_
    noah = (lsm_struct)malloc (sizeof *noah);
#endif
#ifdef _RT_
    chData = (Chem_Data) malloc (sizeof *chData);
#endif
#ifdef _BGC_
    bgc = (bgc_struct) malloc (sizeof *bgc);
#endif
#ifdef _CYCLES_
    cycles = (CyclesStruct) malloc (sizeof *cycles);
#endif

    /* Read PIHM input files */
    ReadAlloc (simulation, pihm);
#ifdef _ENKF_
    /* When running in ensemble mode, use parameters and calibration
     * determined by EnKF module */
    pihm->ctrl.init_type = startmode;
    pihm->ctrl.starttime = starttime;
    pihm->ctrl.endtime = endtime;

    Mbr2Cal (&pihm->cal, param);
#endif

    if (pihm->ctrl.unsat_mode == 2)
    {
        /* State variable size */
        nsv = 3 * pihm->numele + 2 * pihm->numriv;
    }

    /* Initialize CVode state variables */
    CV_Y = N_VNew_Serial (nsv);

    /* Initialize PIHM structure */
    Initialize (pihm, CV_Y, simulation);

    /* Allocate memory for solver */
    cvode_mem = CVodeCreate (CV_BDF, CV_NEWTON);
    if (cvode_mem == NULL)
    {
        printf ("Fatal error: CVodeMalloc failed. \n");
        PihmExit (1);
    }

#ifdef _NOAH_
    /* Read LSM input and initialize LSM structure */
    LsmRead (simulation, noah, pihm);
    LsmInitialize (simulation, pihm, noah);
#endif
#ifdef _BGC_
    /* Read BGC input and initialize BGC structure */
    BgcRead (simulation, bgc, pihm);
    BgcInit (simulation, pihm, noah, bgc);
#endif

#ifdef _CYCLES_
    /* Read Cycles input and initialize Cycles structure */
    CyclesRead (simulation, cycles, pihm);
    CyclesInit (cycles, pihm);
#endif
//#ifdef _RT_
//    chem_alloc(filename, mData, cData, chData, cData->StartTime/60);
//    /* Prepare chem output files */
//    InitialChemFile(filename, chData->NumBTC, chData->BTC_loc);
//#endif

    /* Create output structures */
    MapOutput (simulation, pihm, outputdir);
#ifdef _NOAH_
    MapLsmOutput (simulation, noah, pihm->numele, outputdir);
#endif
#ifdef _BGC_
    MapBgcOutput (simulation, bgc, pihm->numele, outputdir);
#endif

    if (first_cycle == 1)
    {
        //BKInput (simulation, outputdir);
        
        /* initialize output files and structures */
        InitOutputFile (pihm->prtctrl, pihm->ctrl.nprint, pihm->ctrl.ascii);
#ifdef _NOAH_
        InitOutputFile (noah->prtctrl, noah->nprint, pihm->ctrl.ascii);
#endif
#ifdef _BGC_
        InitOutputFile (bgc->prtctrl, bgc->ctrl.nprint, pihm->ctrl.ascii);
#endif
    }

    if (verbose_mode)
    {
        printf ("\n\nSolving ODE system ... \n\n");
    }

    /* Set solver parameters */
    flag = CVodeSetFdata (cvode_mem, pihm);
    flag = CVodeSetInitStep (cvode_mem, (realtype)pihm->ctrl.initstep);
    flag = CVodeSetStabLimDet (cvode_mem, TRUE);
    flag = CVodeSetMaxStep (cvode_mem, (realtype)pihm->ctrl.maxstep);
    flag = CVodeMalloc (cvode_mem, f, (realtype)pihm->ctrl.starttime,
        CV_Y, CV_SS, (realtype)pihm->ctrl.reltol, &pihm->ctrl.abstol);
    flag = CVSpgmr (cvode_mem, PREC_NONE, 0);

    /* Start solver in loops */
    for (i = 0; i < pihm->ctrl.nstep; i++)
    {
        /* Determine current step and next step */
        t = pihm->ctrl.tout[i];
        nextptr = pihm->ctrl.tout[i + 1];

        /* Apply forcing */
        ApplyForcing (&pihm->forcing, t);

        /* Determine if land surface simulation is needed */
        if ((t - pihm->ctrl.starttime) % pihm->ctrl.etstep == 0)
        {
#ifdef _NOAH_
            /* Calculate surface energy balance */
            PIHMxNoah (t, (double)pihm->ctrl.etstep, pihm, noah);

//#ifdef _BGC_
//            BgcCoupling ((int)t, (int)pihm->ctrl.starttime, pihm, noah, bgc);
//#endif
#else
            /* Calculate Interception storage and ET */
            IntcpSnowET (t, (double)pihm->ctrl.etstep, pihm);
#endif
        }

        solvert = (realtype)t;
        flag = CVodeSetMaxNumSteps (cvode_mem,
            (long int)(pihm->ctrl.stepsize * 20));
        flag = CVodeSetStopTime (cvode_mem, (realtype)nextptr);
        flag = CVode (cvode_mem, (realtype)nextptr, CV_Y, &solvert,
            CV_NORMAL_TSTOP);
        flag = CVodeGetCurrentTime (cvode_mem, &cvode_val);

        t = (int)solvert;
        rawtime = (time_t)t;
        timestamp = gmtime (&rawtime);

        if (verbose_mode)
        {
            printf (" Step = %4.4d-%2.2d-%2.2d %2.2d:%2.2d (%d)\n",
                timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min,
                t);
        }
#ifndef _ENKF_
        else if (rawtime % 3600 == 0)
        {
            printf (" Step = %4.4d-%2.2d-%2.2d %2.2d:%2.2d\n",
                timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                timestamp->tm_mday, timestamp->tm_hour,
                timestamp->tm_min);
        }
#endif

        /*
         * Use mass balance to calculate model fluxes or variables
         */
        Summary (pihm, CV_Y, (double) pihm->ctrl.stepsize);

#ifdef _NOAH_
        AvgFlux (noah, pihm);
#endif

#ifdef _NOAH_
        DailyVar (t, pihm->ctrl.starttime, pihm, noah);
#else
        DailyVar (t, pihm->ctrl.starttime, pihm);
#endif

#ifdef _RT_
        /* PIHM-rt control file */
        fluxtrans (t / 60.0, StepSize / 60.0, mData, chData, CV_Y);
#endif

        /*
         * Print outputs
         */
        PrintData (pihm->prtctrl, pihm->ctrl.nprint, t,
            t - pihm->ctrl.starttime, pihm->ctrl.stepsize,
            pihm->ctrl.ascii);
#ifdef _NOAH_
        PrintData (noah->prtctrl, noah->nprint, t,
            t - pihm->ctrl.starttime, pihm->ctrl.stepsize,
            pihm->ctrl.ascii);
#endif
#ifdef _RT_
        /* PIHM-rt output routine */
        PrintChem (filename, chData, t / 60);
#endif

        if ((t - pihm->ctrl.starttime) % 86400 == 0)
        {
#ifdef _BGC_
            if (bgc->ctrl.spinup)
            {
                pihm2metarr (bgc, pihm, t, bgc->ctrl.spinupstart, bgc->ctrl.spinupend);
            }
            else
            {
            }
#endif
            InitDailyStruct (pihm);
        }
    }

#ifdef _BGC_
    if (bgc->ctrl.spinup == 1)
    {
        BgcSpinup (simulation, bgc, pihm, noah);
    }
#endif

    /* Free memory */
    N_VDestroy_Serial (CV_Y);

    /* Free integrator memory */
    CVodeFree (&cvode_mem);

    /*
     * Write init files
     */
    if (pihm->ctrl.write_ic)
    {
        PrtInit (pihm, simulation);
#ifdef _NOAH_
        LsmPrtInit (pihm, noah, simulation);
#endif
    }

#ifdef _NOAH_
    LsmFreeData (pihm, noah);
    free (noah);
#endif
#ifdef _BGC_
    free (bgc);
#endif
    FreeData (pihm);
    free (pihm);
}
