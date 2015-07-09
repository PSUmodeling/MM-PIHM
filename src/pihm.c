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

#ifdef _ENKF_
#include "enkf.h"
#endif

int             verbose_mode;
int             debug_mode;
int             first_cycle;

int main (int argc, char *argv[])
{
    char            project[MAXSTRING];
    char            simulation[MAXSTRING];
    char            outputdir[MAXSTRING];
    char            system_cmd[MAXSTRING];
    char            source_file[MAXSTRING];
    int             overwrite_mode = 0;
    int             c;

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

    while ((c = getopt (argc, argv, "odv")) != -1)
    {
        if (optind >= argc)
        {
            printf ("\nUsage: ./pihm [-o -d -v] <project name>\n");
            printf
                ("\t-o Output will be written in the \"output\" directory and overwrite existing output.\n");
            printf ("\t-v Verbose mode\n");
            printf ("\t-d Debug mode\n");
            exit (1);
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
        exit (1);
    }
    else
    {
        strcpy (project, argv[optind]);
    }

    first_cycle = 1;

    /* Create output directory */
    CreateOutputDir (project, outputdir, overwrite_mode);

    /* The name of the simulation is the same as the project */
    strcpy (simulation, project);

    /* Save input files into output directory */
    sprintf (source_file, "input/%s/%s.para", project, project);
    if (access (source_file, F_OK) != -1)
    {
        sprintf (system_cmd, "cp %s %s/%s.para.bak", source_file, outputdir,
            project);
        system (system_cmd);
    }
    sprintf (source_file, "input/%s/%s.calib", project, simulation);
    if (access (source_file, F_OK) != -1)
    {
        sprintf (system_cmd, "cp %s %s/%s.calib.bak", source_file, outputdir,
            simulation);
        system (system_cmd);
    }
    sprintf (source_file, "input/%s/%s.init", project, project);
    if (access (source_file, F_OK) != -1)
    {
        sprintf (system_cmd, "cp %s %s/%s.init.bak", source_file, outputdir,
            simulation);
        system (system_cmd);
    }

    PIHMRun (simulation, outputdir, first_cycle);

    return (0);
}

void PIHMRun (char *simulation, char *outputdir, int first_cycle)
{
    pihm_struct     pihm;
    N_Vector        CV_Y;       /* State Variables Vector */
    void           *cvode_mem;  /* Model Data Pointer */
    int             flag;       /* Flag to test return value */
    int             nsv;        /* Problem size */
    int             i, j;       /* Loop index */
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

    /* Read PIHM input files */
    ReadAlloc (simulation, pihm);

    if (pihm->ctrl.unsat_mode == 2)
    {
        /* problem size */
        nsv = 3 * pihm->numele + 2 * pihm->numriv;
    }

    /* Initial state variable depending on machine */
    CV_Y = N_VNew_Serial (nsv);

    /* Initialize PIHM model structure */
    Initialize (pihm, CV_Y);

    /* Allocate memory for solver */
    cvode_mem = CVodeCreate (CV_BDF, CV_NEWTON);
    if (cvode_mem == NULL)
    {
        printf ("Fatal error: CVodeMalloc failed. \n");
        exit (1);
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
        /* initialize output files and structures */
        InitOutputFile (pihm->prtctrl, pihm->ctrl.nprint, pihm->ctrl.ascii);
#ifdef _NOAH_
        InitOutputFile (noah->prtctrl, noah->nprint, pihm->ctrl.ascii);
#endif
#ifdef _BGC_
        InitOutputFile (bgc->prtctrl, bgc->ctrl.nprint, pihm->ctrl.ascii);
#endif
        first_cycle = 0;
    }

#ifdef _BGC_
    if (bgc->ctrl.spinup == 1)
    {
        BgcSpinup (simulation, bgc, pihm, noah);
    }
    else
    {
#endif
        if (verbose_mode)
            printf ("\n\nSolving ODE system ... \n\n");

        /* Set solver parameters */
        flag = CVodeSetFdata (cvode_mem, pihm);
        flag = CVodeSetInitStep (cvode_mem, (realtype)pihm->ctrl.initstep);
        flag = CVodeSetStabLimDet (cvode_mem, TRUE);
        flag = CVodeSetMaxStep (cvode_mem, (realtype)pihm->ctrl.maxstep);
        flag = CVodeMalloc (cvode_mem, f, (realtype)pihm->ctrl.starttime,
            CV_Y, CV_SS, (realtype)pihm->ctrl.reltol, &pihm->ctrl.abstol);
        flag = CVSpgmr (cvode_mem, PREC_NONE, 0);

        /* set start time */
        t = pihm->ctrl.starttime;

        /* start solver in loops */
        for (i = 0; i < pihm->ctrl.nstep; i++)
        {
            nextptr = pihm->ctrl.tout[i + 1];

            ApplyForcing (&pihm->forcing, t);

            if ((t - pihm->ctrl.starttime) % pihm->ctrl.etstep == 0)
            {
#ifdef _NOAH_
                /* Calculate surface energy balance */
                PIHMxNoah (t, (double) pihm->ctrl.etstep, pihm, noah);

                for (j = 0; j < noah->nprint; j++)
                {
                    PrintData (&noah->prtctrl[j], t, pihm->ctrl.etstep,
                        pihm->ctrl.ascii);
                }

//#ifdef _BGC_
//                    BgcCoupling ((int) t, (int) cData->StartTime, mData, LSM, BGCM);
//#endif
#else
                /* Calculate Interception Storage and ET */
                IntcpSnowET (t, (double) pihm->ctrl.etstep, pihm);
#endif
            }

            /* Added to adatpt to larger time step. */
            flag = CVodeSetMaxNumSteps (cvode_mem, (long int) (pihm->ctrl.stepsize * 20));
            solvert = (realtype) t;
            flag = CVode (cvode_mem, (realtype) nextptr, CV_Y, &solvert, CV_NORMAL);
            flag = CVodeGetCurrentTime (cvode_mem, &cvode_val);

            t = (int) solvert;
            rawtime = t;
            timestamp = gmtime (&rawtime);

            if (verbose_mode)
            {
                printf (" Time = %4.4d-%2.2d-%2.2d %2.2d:%2.2d (%d)\n",
                    timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                    timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min,
                    t);
            }
            else if (rawtime % 3600 == 0)
            {
                printf (" Time = %4.4d-%2.2d-%2.2d %2.2d:%2.2d\n",
                    timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                    timestamp->tm_mday, timestamp->tm_hour,
                    timestamp->tm_min);
            }

            Summary (pihm, CV_Y, (double) pihm->ctrl.stepsize);
#ifdef _NOAH_
            AvgFlux (noah, pihm);
#endif
#ifdef _RT_
            /* PIHM-rt control file */
            fluxtrans (t / 60.0, StepSize / 60.0, mData, chData, CV_Y);
#endif

            /* Print outputs */
            for (j = 0; j < pihm->ctrl.nprint; j++)
            {
                PrintData (&pihm->prtctrl[j], t, (int) pihm->ctrl.stepsize,
                    pihm->ctrl.ascii);
            }
#ifdef _RT_
            /* PIHM-rt output routine */
            PrintChem (filename, chData, t / 60);
#endif
        }
#ifdef _BGC_
    }
#endif
//    }
//
//    if (cData->Spinup)
//    {
//        PrintInit (mData, filename);
//#ifdef _FLUX_PIHM_
//        LSM_PrintInit (mData, LSM, filename);
//#endif
//    }
//
//    printf ("\n Done. \n");
    /* Free memory */
    N_VDestroy_Serial (CV_Y);

    /* Free integrator memory */
    CVodeFree (&cvode_mem);

    FreeData (pihm);
//#ifdef _FLUX_PIHM_
//    LSM_FreeData (mData, LSM);
//    free (LSM);
//#endif
    free (pihm);
}

void CreateOutputDir (char *project, char *outputdir, int overwrite_mode)
{
    time_t          rawtime;
    struct tm      *timestamp;
    char            str[11];

    if (0 == (mkdir ("output", 0755)))
        printf (" Output directory was created.\n\n");

    time (&rawtime);
    timestamp = localtime (&rawtime);

    if (overwrite_mode)
    {
        strcpy (outputdir, "output/");
    }
    else
    {
        /* Create output directory based on projectname and time */
        sprintf (str, "%2.2d%2.2d%2.2d%2.2d%2.2d",
            timestamp->tm_year + 1900 - 2000, timestamp->tm_mon + 1,
            timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
        sprintf (outputdir, "output/%s.%s/", project, str);
        printf ("\nOutput directory: %s\n", outputdir);
        mkdir (outputdir, 0755);
    }
}
