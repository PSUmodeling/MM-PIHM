#include "pihm.h"

void PIHMParal (int id, int p, char *outputdir)
{
    char            simulation[MAXSTRING];
    int             job_per_node;
    int             startmode;
    int             starttime;
    int             endtime;
    double         *param;
    int             first_cycle = 1;
    int             success;
    int             ierr;
    int             ii;
    enkf_struct     ens;

    MPI_Bcast (project, MAXSTRING, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast (outputdir, MAXSTRING, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast (&debug_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&verbose_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (0 == id)
    {
        /*
         * EnKF initialization
         */
        ens = (enkf_struct)malloc (sizeof *ens);

        /* Read EnKF input file */
        ReadEnKF (ens);

        /* Check if node number is appropriate */
        if (ens->ne % (p - 1) != 0)
        {
            fprintf (stderr,
                "Error: Please specify a correct node number!\n");
            PIHMexit (EXIT_FAILURE);
        }
        else
        {
            job_per_node = ens->ne / (p - 1);
        }

        /* Initialize observation operator vector */
        InitEns (ens);

        /* Perturb model parameters */
        Perturb (ens, outputdir);
    }

    /* Broadcast jobs per node and output directory to all nodes */
    MPI_Bcast (&job_per_node, 1, MPI_INT, 0, MPI_COMM_WORLD);
    param =
        (double *)malloc (job_per_node * (p - 1) * MAXPARAM *
        sizeof (double));

    /*
     * EnKF cycles
     */
    while (1)
    {
        if (0 == id)
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
            ReadVar (outputdir, ens, ens->cycle_end_time);

            /* EnKF data assimilation */
            EnKFDA (ens, ens->cycle_end_time, outputdir);

            /* Proceed to next cycle */
            ens->mbr_start_mode = RST_FILE;
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
                PIHM (simulation, outputdir, first_cycle,
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

    if (0 == id)
    {
        FreeEns (ens);
        free (ens);
    }

    free (param);
}
