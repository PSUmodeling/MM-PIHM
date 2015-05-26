/*****************************************************************************
 * File		: pihm.c
 * Function	: Main program file
 * Developer of MM-PIHM     :   Yuning Shi  (yshi@psu.edu)
 *                              Chen Bao    (cub200@psu.edu)
 *                              Yu Zhang    (yzz120@psu.edu)
 * Developer of PIHM 2.2    :	Xuan Yu	    (xxy113@psu.edu)
 * Developer of PIHM 2.0    :	Mukesh Kumar	(muk139@psu.edu)
 * Developer of PIHM 1.0    :	Yizhong Qu	(quyizhong@gmail.com)
 *----------------------------------------------------------------------------
 * This code is free for research purpose only.
 * Please provide relevant references if you use this code in your research
 * 	work
 *----------------------------------------------------------------------------
 * DEVELOPMENT RELATED REFERENCES:
 * Flux-PIHM:
 *  a) Shi, Y. et al., 2013, "Development of a coupled land surface
 *      hydrologic model and evaluation at a critical zone observatory",
 *      Journal of Hydrometeorology, 14, 1401--1420.
 *  b) Shi, Y. et al., 2014, "Evaluation of the parameter sensitivity of 
 *	a coupled land surface hydrologic model. Journal of Hydrometeorology,
 *	15, 279--299.
 * PIHM V2:
 *  a) Kumar, M., 2008, "Development and Implementation of a Multiscale, 
 *	Multiprocess Hydrologic Model". PhD Thesis, Penn State University
 * PIHM V1:
 *  a) Qu, Y., 2005, "An Integrated hydrologic model for multiprocec
 *	simulation using semi-discrete finite volume approach".PhD Thesis, PSU
 *  b) Qu, Y. & C. Duffy, 2007, "A semidiscrete finite volume formulation
 *	for multiprocess watershed simulation". Water Resources Research
 ****************************************************************************/

#include "pihm.h"               /* Data Model and Variable Declarations */

#ifdef _FLUX_PIHM_
#include "noah.h"
#endif

#ifdef _RT_
#include "rt.h"
#endif

#ifdef _BGC_
#include "bgc.h"
#endif

#ifdef _ENKF_
#include "EnKF.h"
#endif

int main (int argc, char *argv[])
{
    struct tm      *timestamp;
    time_t         *rawtime;
    char            project[50];
    char           *simulation;
    char           *outputdir;
    char            str[11];
    char            system_cmd[1024];
    char            source_file[1024];
    int             overwrite_mode = 0;
    int             c;
    int             verbose = 0;
    int             debug = 0;
    int             first_cycle = 1;

#ifdef _ENKF_
    MPI_Status      status;
    ensemble_struct ensemble;
#endif

#ifdef _ENKF_
    ierr = MPI_Init ( &argc, &argv );
    ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );
    ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );

    if (id == 0)
    {
#endif
        rawtime = (time_t *) malloc (sizeof (time_t));

        printf ("\n");
        printf ("\t\t########  #### ##     ## ##     ##\n");
        printf ("\t\t##     ##  ##  ##     ## ###   ###\n"); 
        printf ("\t\t##     ##  ##  ##     ## #### ####\n");
        printf ("\t\t########   ##  ######### ## ### ##\n");
        printf ("\t\t##         ##  ##     ## ##     ##\n");
        printf ("\t\t##         ##  ##     ## ##     ##\n"); 
        printf ("\t\t##        #### ##     ## ##     ##\n");
        printf ("\n\t    The Penn State Integrated Hydrologic Model\n");

#ifdef _FLUX_PIHM_
        printf ("\n\t    * Land surface module turned on.\n");
#endif
#ifdef _RT_
        printf ("\n\t    * Reactive transport land surface hydrological mode.\n");
#endif
#ifdef _BGC_
        printf ("\n\t    * Biogeochemistry module turned on.\n");
#endif

        if (0 == (mkdir ("output", 0755)))
            printf (" Output directory was created.\n\n");

        while((c = getopt(argc, argv, "odv")) != -1)
        {
            if (optind >= argc)
            {
                printf ("\nUsage: ./pihm [-o -d -v] <project name>\n");
                printf ("\t-o Output will be written in the \"output\" directory and overwrite existing output.\n");
                printf ("\t-v Verbose mode\n");
                printf ("\t-d Debug mode\n");
                exit (1);
            }
            switch(c)
            {
                case 'o':
                    overwrite_mode = 1;
                    printf ("Overwrite mode turned on. Output directory is \"./output\".\n");
                    break;
                case 'd':
                    debug = 1;
                    printf ("Debug mode turned on.\n");
                    break;
                case 'v':
                    verbose = 1;
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
            printf ("\t-o Output will be written in the \"output\" directory and overwrite existing output.\n");
            printf ("\t-v Verbose mode\n");
            printf ("\t-d Debug mode\n");
            exit (1);
        }
        else
            strcpy (project, argv[optind]);

        time (rawtime);
        timestamp = localtime (rawtime);

        if (overwrite_mode == 1)
        {
            outputdir = (char *)malloc (8 * sizeof (char));
            strcpy (outputdir, "output/");
        }
        else
        {
            /* Create output directory based on projectname and time */
            sprintf (str, "%2.2d%2.2d%2.2d%2.2d%2.2d", timestamp->tm_year + 1900 - 2000, timestamp->tm_mon + 1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
            outputdir = (char *)malloc ((strlen (project) + 20) * sizeof (char));
            sprintf (outputdir, "output/%s.%s/", project, str);
            printf ("\nOutput directory: %s\n", outputdir);
            mkdir (outputdir, 0755);
        }

#ifdef _ENKF_
        mData = (Model_Data) malloc (sizeof *mData);

        ReadPara ();

        En = (ensemble) malloc (sizeof *En);

        read_EnKF (project, En);

        if (En->ne % (p - 1) != 0)
        {
            printf("\n ERROR: Please specify a correct node number!");
            ierr = MPI_Finalize();
            exit(1);
        }

        /* Perturb model parameters */
        perturb(project, mData, En);
    }

    while (1)
    {
        if (id == 0)
        {
            start_time = (int) cData->StartTime;
            end_time = (int) cData->EndTime;
            init_type = cData->init_type;
            ne = ensemble->ne;

            if (cData->StartTime >= ensemble->EndTime)
            {
                init_type = -999;

                for (j = 1; j < p; j++)
                {
                    MPI_Send (&StartTime, 1, MPI_INT, j, START_TIME, MPI_COMM_WORLD);
                    MPI_Send (&EndTime, 1, MPI_INT, j, END_TIME, MPI_COMM_WORLD);
                    MPI_Send (&init_type, 1, MPI_INT, j, INIT_TYPE, MPI_COMM_WORLD);
                    MPI_Send (&ne, 1, MPI_INT, j, NE, MPI_COMM_WORLD);
                    MPI_Send (&project, 50, MPI_CHAR, j, PROJECT, MPI_COMM_WORLD);
                }
                break;
            }
            else
            {
                for (j = 1; j < p; j++)
                {
                    MPI_Send (&StartTime, 1, MPI_INT, j, START_TIME, MPI_COMM_WORLD);
                    MPI_Send (&EndTime, 1, MPI_INT, j, END_TIME, MPI_COMM_WORLD);
                    MPI_Send (&init_type, 1, MPI_INT, j, INIT_TYPE, MPI_COMM_WORLD);
                    MPI_Send (&ne, 1, MPI_INT, j, NE, MPI_COMM_WORLD);
                    MPI_Send (&project, 50, MPI_CHAR, j, PROJECT, MPI_COMM_WORLD);
                }
            }

            job_hand_in (p - 1);

            obs_time = (int)cData->EndTime;

            read_initial (project, mData, En, obs_time);

            EnKF (projname, mData, En, obs_time);

            cData->init_type = 3;
            cData->StartTime = cData->EndTime;

            interval = En->interval;

            cData->EndTime = cData->EndTime + interval * 60.0 * 60.0;
        }
        else
        {
            MPI_Recv(&StartTime, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&EndTime, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&init_type, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(&ne, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);
            MPI_Recv(&project, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);

            if (init_type == -999)
                break;

            for (i = (id - 1) * ne / (p - 1); i < id * ne / (p - 1); i++)
            {
                simulation = (char *) malloc ((strlen (project) + 5) * sizeof(char));
                sprintf (indchar, "%3.3d", i+1);
                strcpy (simulation, project);
                strcat (simulation, ".");
                strcat (simulation, indchar);
#else
                simulation = (char *) malloc ((strlen (project) + 1) * sizeof (char));
                strcpy (simulation, project);
#endif

                if (first_cycle)
                {
                    /* Save input files into output directory */
                    sprintf (source_file, "input/%s/%s.para", project, project);
                    if (access (source_file, F_OK) != -1)
                    {
                        sprintf (system_cmd, "cp %s %s/%s.para.bak", source_file, outputdir, project);
                        system (system_cmd);
                    }
                    sprintf (source_file, "input/%s/%s.calib", project, simulation);
                    if (access (source_file, F_OK) != -1)
                    {
                        sprintf (system_cmd, "cp %s %s/%s.calib.bak", source_file, outputdir, simulation);
                        system (system_cmd);
                    }
                    sprintf (source_file, "input/%s/%s.init", project, project);
                    if (access (source_file, F_OK) != -1)
                    {
                        sprintf (system_cmd, "cp %s %s/%s.init.bak", source_file, outputdir, simulation);
                        system (system_cmd);
                    }
                }

                pihm (project, verbose, debug, outputdir, first_cycle);

                first_cycle = 0;
#ifdef _ENKF_
            }

            ierr = MPI_Send (&EndTime, 1, MPI_INT, 0, 5, MPI_COMM_WORLD);
        }
    }

    if (id == 0)
    {
        freeEnsemble(mData, En);
        free(En);
        free(rawtime);
        FreeData(mData, &cData);
        free(mData);
    }

    ierr = MPI_Finalize();
#endif

    return (0);
}

void pihm (char *project, int verbose, int debug, char *output_dir, int first_cycle)
{
    Model_Data      mData;      /* Model Data */
    Control_Data    cData;      /* Solver Control Data */
    N_Vector        CV_Y;       /* State Variables Vector */
    void           *cvode_mem;  /* Model Data Pointer */
    int             flag;       /* flag to test return value */
    int             N;          /* Problem size */
    int             i, j;       /* loop index */
    realtype        t;          /* simulation time */
    struct tm      *timestamp;
    time_t         *rawtime;
    realtype        NextPtr;
    realtype        StepSize;   /* stress period & step size */
    realtype        cvode_val;
    char           *filename;

#ifdef _FLUX_PIHM_
    LSM_STRUCT      LSM;
#endif

#ifdef _RT_
    Chem_Data       chData;
#endif 

#ifdef _BGC_
    bgc_struct      BGCM;
#endif


    rawtime = (time_t *) malloc (sizeof (time_t));
    filename = (char *) malloc ((strlen (project) + 1) * sizeof (char));
    strcpy (filename, project);


    /* Allocate memory for model data structure */
    mData = (Model_Data) malloc (sizeof *mData);

    cData = (Control_Data) malloc (sizeof *cData);

    cData->Verbose = verbose;
    cData->Debug = debug;

#ifdef _FLUX_PIHM_
    LSM = (LSM_STRUCT) malloc (sizeof *LSM);
#endif
#ifdef _RT_
    chData = (Chem_Data) malloc (sizeof * chData);
#endif
#ifdef _BGC_
    BGCM = (bgc_struct) malloc (sizeof *BGCM);
#endif
    read_alloc (filename, mData, cData);
#ifdef _FLUX_PIHM_
    LSM_read (filename, LSM, cData);
#endif
#ifdef _BGC_
    BGC_read (filename, BGCM, mData);
#endif

    if (mData->UnsatMode == 2)
    {
        /* problem size */
        N = 3 * mData->NumEle + 2 * mData->NumRiv;
        mData->DummyY = (realtype *) malloc ((3 * mData->NumEle + 2 * mData->NumRiv) * sizeof (realtype));
    }
    /* initial state variable depending on machine */
    CV_Y = N_VNew_Serial (N);

    /* initialize mode data structure */
    initialize (filename, mData, cData, CV_Y);
#ifdef _FLUX_PIHM_
    LSM_initialize (filename, mData, cData, LSM);
#endif
#ifdef _RT_
    chem_alloc(filename, mData, cData, chData, cData->StartTime/60);
    /* Prepare chem output files */
    InitialChemFile(filename, chData->NumBTC, chData->BTC_loc);
#endif
#ifdef _BGC_
    BGC_init (filename, mData, LSM, BGCM);
#endif

    if (first_cycle == 1)
    {
        /* initialize output files and structures */
        initialize_output (filename, mData, cData, output_dir);
#ifdef _FLUX_PIHM_
        LSM_initialize_output (filename, mData, cData, LSM, output_dir);
#endif
    }

#ifdef _BGC_
    if (BGCM->ctrl.spinup == 1)
    {
        bgc_spinup (filename, BGCM, mData, LSM);
    }
    else
    {
#endif
        if (cData->Verbose == 1)
            printf ("\n\nSolving ODE system ... \n\n");

        /* allocate memory for solver */
        cvode_mem = CVodeCreate (CV_BDF, CV_NEWTON);
        if (cvode_mem == NULL)
        {
            printf ("Fatal error: CVodeMalloc failed. \n");
            exit (1);
        }

        flag = CVodeSetFdata (cvode_mem, mData);
        flag = CVodeSetInitStep (cvode_mem, cData->InitStep);
        flag = CVodeSetStabLimDet (cvode_mem, TRUE);
        flag = CVodeSetMaxStep (cvode_mem, cData->MaxStep);
        flag = CVodeMalloc (cvode_mem, f, cData->StartTime, CV_Y, CV_SS, cData->reltol, &cData->abstol);
        flag = CVSpgmr (cvode_mem, PREC_NONE, 0);
        //flag = CVSpgmrSetGSType(cvode_mem, MODIFIED_GS);

        /* set start time */
        t = cData->StartTime;

        /* start solver in loops */
        for (i = 0; i < cData->NumSteps; i++)
        {
            /* inner loops to next output points with ET step size control */
            while (t < cData->Tout[i + 1])
            {
                if (t + cData->ETStep >= cData->Tout[i + 1])
                    NextPtr = cData->Tout[i + 1];
                else
                    NextPtr = t + cData->ETStep;
                StepSize = NextPtr - t;

                mData->dt = StepSize;

                if ((int)t % (int)cData->ETStep == 0)
                {
#ifdef _FLUX_PIHM_
                    for (j = 0; j < mData->NumEle; j++)
                    {
                        mData->avg_inf[j] = (mData->avg_inf[j] + mData->EleViR[j]) / (cData->ETStep / StepSize);
                        mData->avg_subflux[j][0] = (mData->avg_subflux[j][0] + mData->FluxSub[j][0]) / (cData->ETStep / StepSize);
                        mData->avg_subflux[j][1] = (mData->avg_subflux[j][1] + mData->FluxSub[j][1]) / (cData->ETStep / StepSize);
                        mData->avg_subflux[j][2] = (mData->avg_subflux[j][2] + mData->FluxSub[j][2]) / (cData->ETStep / StepSize);
                    }

                    /* calculate surface energy balance */
                    PIHM2Noah (t, cData->ETStep, mData, LSM);
                    Noah2PIHM (mData, LSM);

                    for (j = 0; j < mData->NumEle; j++)
                    {
                        mData->avg_inf[j] = 0.0; 
                        mData->avg_subflux[j][0] = 0.0;
                        mData->avg_subflux[j][1] = 0.0;
                        mData->avg_subflux[j][2] = 0.0;
                    }
                }
                else
                {
                    for (j = 0; j < mData->NumEle; j++)
                    {
                        mData->avg_inf[j] = mData->avg_inf[j] + mData->EleViR[j];
                        mData->avg_subflux[j][0] = mData->avg_subflux[j][0] + mData->FluxSub[j][0];
                        mData->avg_subflux[j][1] = mData->avg_subflux[j][1] + mData->FluxSub[j][1];
                        mData->avg_subflux[j][2] = mData->avg_subflux[j][2] + mData->FluxSub[j][2];
                    }
#else
                    /* calculate Interception Storage and ET */
                    is_sm_et (t, cData->ETStep, mData, CV_Y);
#endif
                }


#ifdef COUPLE_I
                t = NextPtr;
#else
                /* Added to adatpt to larger time step. YS */
                flag = CVodeSetMaxNumSteps(cvode_mem, (long int)(StepSize* 20));
                flag = CVode (cvode_mem, NextPtr, CV_Y, &t, CV_NORMAL);
                flag = CVodeGetCurrentTime(cvode_mem, &cvode_val);
#endif
                *rawtime = (int)t;
                timestamp = gmtime (rawtime);

                if (cData->Verbose)
                    printf (" Time = %4.4d-%2.2d-%2.2d %2.2d:%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
                else if ((int)*rawtime % 3600 == 0)
                    printf (" Time = %4.4d-%2.2d-%2.2d %2.2d:%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);

                summary (mData, CV_Y, t - StepSize, StepSize);
#ifdef _RT_
                /* PIHM-rt control file */
                fluxtrans(t / 60.0, StepSize / 60.0, mData, chData, CV_Y);
#endif
                update (t, mData);
            }
#ifdef _BGC_
        }
#endif
        /* Print outputs */
        for (j = 0; j < cData->NumPrint; j++)
            PrintData (cData->PCtrl[j], t, StepSize, cData->Ascii);
#ifdef _FLUX_PIHM_
        for (j = 0; j < LSM->NPRINT; j++)
            PrintData (LSM->PCtrl[j], t, StepSize, cData->Ascii);
#endif
#ifdef _RT_
        /* PIHM-rt output routine */
        PrintChem(filename, chData,t/60);
#endif
    }

    if (cData->Spinup)
    {
        PrintInit (mData, filename);
#ifdef _FLUX_PIHM_
        LSM_PrintInit (mData, LSM, filename);
#endif
    }

    printf ("\n Done. \n");
    /* Free memory */
    N_VDestroy_Serial (CV_Y);

    /* Free integrator memory */
    CVodeFree (&cvode_mem);

    free (rawtime);
    free (filename);
    FreeData (mData, cData);
#ifdef _FLUX_PIHM_
    LSM_FreeData (mData, LSM);
    free (LSM);
#endif
    free (mData);
    free (cData);
}
