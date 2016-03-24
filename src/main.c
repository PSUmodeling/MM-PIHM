#include "pihm.h"

int             verbose_mode;
int             debug_mode;
char            project[MAXSTRING];

int main (int argc, char *argv[])
{
    char            outputdir[MAXSTRING];
    char            simulation[MAXSTRING];
    int             spec_output_mode = 0;
    int             c;
    int             first_cycle = 1;
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
        printf
            ("\n\t    * Reactive transport land surface hydrological mode.\n");
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
    while ((c = getopt (argc, argv, "o:dvl")) != -1)
    {
        switch (c)
        {
            case 'o':
                /* Specify output directory */
                sprintf (outputdir, "output/%s/", optarg);
                spec_output_mode = 1;
                printf
                    ("Output directory is specified as \"%s\".\n", outputdir);
                break;
            case 'd':
                /* Debug mode */
                debug_mode = 1;
                printf ("Debug mode turned on.\n");
                break;
            case 'v':
                /* Verbose mode */
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
        printf ("Usage: ./pihm [-o output_dir] [-d] [-v] <project name>\n");
        printf ("\t-o Specify output directory.\n");
        printf ("\t-v Verbose mode\n");
        printf ("\t-d Debug mode\n");
        PihmExit (1);
    }
    else
    {
        strcpy (project, argv[optind]);
    }

#ifdef _ENKF_
    if (id == 0)
    {
#endif
        /*
         * Create output directory
         */
        CreateOutputDir (outputdir, spec_output_mode);

#ifdef _ENKF_
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
        InitEns (ens);

        /* Perturb model parameters */
        Perturb (project, ens, outputdir);
    }

    /* Broadcast jobs per node and output directory to all nodes */
    MPI_Bcast (&job_per_node, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (outputdir, MAXSTRING, MPI_CHAR, 0, MPI_COMM_WORLD);
    param =
        (double *)malloc (job_per_node * (p - 1) * MAXPARAM *
        sizeof (double));

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
            EnKF (ens, ens->cycle_end_time, outputdir);

            /* Proceed to next cycle */
            ens->mbr_start_mode = 3;
            ens->cycle_start_time = ens->cycle_end_time;
            ens->cycle_end_time += ens->interval;
        }
        else
        {
            /* Receive required parameters from Node 0 */
            JobRecv (&starttime, &endtime, &startmode, param,
                job_per_node * (p - 1));

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
            ierr =
                MPI_Send (&success, 1, MPI_INT, 0, SUCCESS_TAG,
                MPI_COMM_WORLD);
        }
    }
#else
    /* The name of the simulation is the same as the project */
    strcpy (simulation, project);

    PIHMRun (simulation, outputdir, first_cycle);
#endif

#ifdef _ENKF_
    if (id == 0)
    {
#endif
        printf ("\nSimulation completed.\n");

#ifdef _ENKF_
        FreeEns (ens);
        free (ens);
    }

    free (param);
    ierr = MPI_Finalize ();
#endif

    return (0);
}
