#include "pihm.h"
#ifdef _NOAH_
#include "noah.h"
#endif
#include "enkf.h"

void JobHandout (int starttime, int endtime, int startmode, ens_mbr_struct *member,
    double *param, int ne, int total_jobs)
{
    int             dest;
    int             msg[3];
    int             i, j;

    msg[0] = starttime;
    msg[1] = endtime;
    msg[2] = startmode;

    for (i = 0; i < ne; i++)
    {
        for (j = 0; j < MAXPARAM; j++)
        {
            param[i * MAXPARAM + j] = member[i].param[j];
        }
    }

    for (dest = 1; dest < total_jobs + 1; dest++)
    {
        MPI_Send (msg, 3, MPI_INT, dest, CYCLE_TAG, MPI_COMM_WORLD);
        MPI_Send (param, ne * MAXPARAM, MPI_DOUBLE, dest, PARAM_TAG, MPI_COMM_WORLD);
    }
}

void JobRecv (int *starttime, int *endtime, int *startmode, double *param, int ne)
{
    int             msg[3];
    MPI_Status      status;

    MPI_Recv (msg, 3, MPI_INT, 0, CYCLE_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv (param, ne * MAXPARAM, MPI_DOUBLE, 0, PARAM_TAG, MPI_COMM_WORLD, &status);

    *starttime = msg[0];
    *endtime = msg[1];
    *startmode = msg[2];
}

void PrintEnKFStatus (int starttime, int endtime)
{
    time_t          rawtime;
    struct tm      *timestamp;

    printf ("\nRunning ensemble members from ");

    rawtime = (int) starttime;
    timestamp = gmtime (&rawtime);
    printf ("%4.4d-%2.2d-%2.2d %2.2d:%2.2d to ",
        timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday,
        timestamp->tm_hour, timestamp->tm_min);

    rawtime = (int) endtime;
    timestamp = gmtime (&rawtime);
    printf ("%4.4d-%2.2d-%2.2d %2.2d:%2.2d\n",
        timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday,
        timestamp->tm_hour, timestamp->tm_min);
}

void JobHandIn (int total_jobs)
{
    int             ierr;
    int             success;
    int             received = 0;
    int             source;

    MPI_Status      status;

    while (received < total_jobs)
    {
        ierr = MPI_Recv (&success, 1, MPI_INT, MPI_ANY_SOURCE, SUCCESS_TAG, MPI_COMM_WORLD, &status);
        received++;
        source = status.MPI_SOURCE;
        //printf("PIHM job handed in from Node: %d\n", source);
    }
}

void WritePara (char *project, int start_mode, int start_time, int end_time)
{
    ctrl_struct     ctrl;
    char            fn[MAXSTRING];
    time_t          rawtime;
    struct tm      *timestamp;
    FILE           *fid;

    ReadPara (project, &ctrl);

    sprintf (fn, "input/%s/%s.para", project, project);
    fid = fopen (fn, "w");
    CheckFile (fid, fn);

    fprintf (fid, "%-20s%d\n", "INIT_MODE", start_mode);
    fprintf (fid, "%-20s%d\n", "ASCII_OUTPUT", ctrl.ascii);
    fprintf (fid, "%-20s%d\n", "WRITE_IC", ctrl.write_ic);
    fprintf (fid, "%-20s%d\n", "UNSAT_MODE", ctrl.unsat_mode);
    fprintf (fid, "%-20s%d\n", "SURF_MODE", ctrl.surf_mode);
    fprintf (fid, "%-20s%d\n", "RIV_MODE", ctrl.riv_mode);
    fprintf (fid, "%-20s%d\n", "SOLVER", ctrl.solver);
    fprintf (fid, "%-20s%lg\n", "ABSTOL", ctrl.abstol);
    fprintf (fid, "%-20s%lg\n", "RELTOL", ctrl.reltol);
    fprintf (fid, "%-20s%lg\n", "INIT_SOLVER_STEP", ctrl.initstep);
    fprintf (fid, "%-20s%lf\n", "MAX_SOLVER_STEP", ctrl.maxstep);
    fprintf (fid, "%-20s%d\n", "LSM_STEP", ctrl.etstep);

    rawtime = (time_t) start_time;
    timestamp = gmtime (&rawtime);
    fprintf (fid, "%-20s%4.4d-%2.2d-%2.2d %2.2d:%2.2d\n", "START",
                timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
    rawtime = (time_t) end_time;
    timestamp = gmtime (&rawtime);
    fprintf (fid, "%-20s%4.4d-%2.2d-%2.2d %2.2d:%2.2d\n", "END",
                timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);

    fprintf (fid, "%-20s%d\n", "MODEL_STEPSIZE", ctrl.stepsize);

    fprintf (fid, "%-20s%d\n", "GW", ctrl.prtvrbl[GW_CTRL]);
    fprintf (fid, "%-20s%d\n", "SURF", ctrl.prtvrbl[SURF_CTRL]);
    fprintf (fid, "%-20s%d\n", "SNOW", ctrl.prtvrbl[SNOW_CTRL]);
    fprintf (fid, "%-20s%d\n", "RIVSTG", ctrl.prtvrbl[RIVSTG_CTRL]);
    fprintf (fid, "%-20s%d\n", "INFIL", ctrl.prtvrbl[INFIL_CTRL]);
    fprintf (fid, "%-20s%d\n", "RECHARGE", ctrl.prtvrbl[RECHARGE_CTRL]);
    fprintf (fid, "%-20s%d\n", "CMC", ctrl.prtvrbl[CMC_CTRL]);
    fprintf (fid, "%-20s%d\n", "UNSAT", ctrl.prtvrbl[UNSAT_CTRL]);
    fprintf (fid, "%-20s%d\n", "EC", ctrl.prtvrbl[EC_CTRL]);
    fprintf (fid, "%-20s%d\n", "ETT", ctrl.prtvrbl[ETT_CTRL]);
    fprintf (fid, "%-20s%d\n", "EDIR", ctrl.prtvrbl[EDIR_CTRL]);
    fprintf (fid, "%-20s%d\n", "RIVFLX0", ctrl.prtvrbl[RIVFLX0_CTRL]);
    fprintf (fid, "%-20s%d\n", "RIVFLX1", ctrl.prtvrbl[RIVFLX1_CTRL]);
    fprintf (fid, "%-20s%d\n", "RIVFLX2", ctrl.prtvrbl[RIVFLX2_CTRL]);
    fprintf (fid, "%-20s%d\n", "RIVFLX3", ctrl.prtvrbl[RIVFLX3_CTRL]);
    fprintf (fid, "%-20s%d\n", "RIVFLX4", ctrl.prtvrbl[RIVFLX4_CTRL]);
    fprintf (fid, "%-20s%d\n", "RIVFLX5", ctrl.prtvrbl[RIVFLX5_CTRL]);
    fprintf (fid, "%-20s%d\n", "RIVFLX6", ctrl.prtvrbl[RIVFLX6_CTRL]);
    fprintf (fid, "%-20s%d\n", "RIVFLX7", ctrl.prtvrbl[RIVFLX7_CTRL]);
    fprintf (fid, "%-20s%d\n", "RIVFLX8", ctrl.prtvrbl[RIVFLX8_CTRL]);
    fprintf (fid, "%-20s%d\n", "RIVFLX9", ctrl.prtvrbl[RIVFLX9_CTRL]);
    fprintf (fid, "%-20s%d\n", "RIVFLX10", ctrl.prtvrbl[RIVFLX10_CTRL]);
    fprintf (fid, "%-20s%d\n", "SUBFLX", ctrl.prtvrbl[SUBFLX_CTRL]);
    fprintf (fid, "%-20s%d\n", "TOTALFLX", ctrl.prtvrbl[TOTALFLX_CTRL]);

    fclose (fid);
}

void InitOper (char *project, enkf_struct ens)
{
    pihm_struct     pihm;
    N_Vector        CV_Y;       /* State Variables Vector */
#ifdef _NOAH_
    lsm_struct      noah;
#endif
    int             nsv;
    int             i;

    char            outputdir[MAXSTRING];

    outputdir[0] = '\0';

    pihm = (pihm_struct)malloc (sizeof *pihm);
#ifdef _NOAH_
    noah = (lsm_struct)malloc (sizeof *noah);
#endif
    ReadAlloc (project, pihm);

    if (pihm->ctrl.unsat_mode == 2)
    {
        /* problem size */
        nsv = 3 * pihm->numele + 2 * pihm->numriv;
    }

    CV_Y = N_VNew_Serial (nsv);

    Initialize (pihm, CV_Y, project);

#ifdef _NOAH_
    /* Read LSM input and initialize LSM structure */
    LsmRead (project, noah, pihm);
    LsmInitialize (project, pihm, noah);
#endif
    MapOutput (project, pihm, outputdir);
#ifdef _NOAH_
    MapLsmOutput (project, noah, pihm->numele, outputdir);
#endif

    ens->numele = pihm->numele;
    ens->numriv = pihm->numriv;
    ens->ascii = pihm->ctrl.ascii;


    for (i = 0; i < ens->nobs; i++)
    {
        if (strcasecmp (ens->obs[i].name, "discharge") == 0)
        {
            ens->obs[i].type = RUNOFF_OBS;
            DisOper (&ens->obs[i], pihm);
        }
        else if (strcasecmp (ens->obs[i].name, "skin_temperature") == 0)
        {
            ens->obs[i].type = TSKIN_OBS;
            LandSfcTmpOper (&ens->obs[i], pihm);
        }
        else
        {
            printf ("ERROR: Cannot find the operator for %s!\n",
                ens->obs[i].name);
            PihmExit (1);
        }
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
