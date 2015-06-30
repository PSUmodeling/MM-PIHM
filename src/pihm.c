#include "pihm.h"

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

#ifdef _FLUX_PIHM_
    printf ("\n\t    * Land surface module turned on.\n");
#endif
#ifdef _RT_
    printf ("\n\t    * Reactive transport land surface hydrological mode.\n");
#endif
#ifdef _BGC_
    printf ("\n\t    * Biogeochemistry module turned on.\n");
#endif

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
        printf ("\t-o Output will be written in the \"output\" directory and overwrite existing output.\n");
        printf ("\t-v Verbose mode\n");
        printf ("\t-d Debug mode\n");
        exit (1);
    }
    else
        strcpy (project, argv[optind]);

    first_cycle = 1;

    /* Create output directory */
    CreateOutputDir (project, outputdir, overwrite_mode);

    strcpy (simulation, project);

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

    PIHMRun (project, outputdir, first_cycle);

    first_cycle = 0;

    return (0);
}

void PIHMRun (char *simulation, char *outputdir, int first_cycle)
{
    pihm_struct     pihm;
    N_Vector        CV_Y;       /* State Variables Vector */
    void           *cvode_mem;  /* Model Data Pointer */
    int             flag;       /* flag to test return value */
    int             nsv;          /* Problem size */
    int             i, j;       /* loop index */
    int             t;          /* simulation time */
    realtype        solvert;
    struct tm      *timestamp;
    time_t          rawtime;
    realtype        nextptr;
    realtype        stepsize;   /* stress period & step size */
    realtype        cvode_val;

#ifdef _FLUX_PIHM_
    LSM_STRUCT      LSM;
#endif

#ifdef _RT_
    Chem_Data       chData;
#endif 

#ifdef _BGC_
    bgc_struct      BGCM;
#endif

    /* Allocate memory for model data structure */
    pihm = (pihm_struct) malloc (sizeof *pihm);

#ifdef _FLUX_PIHM_
    LSM = (LSM_STRUCT) malloc (sizeof *LSM);
#endif
#ifdef _RT_
    chData = (Chem_Data) malloc (sizeof * chData);
#endif
#ifdef _BGC_
    BGCM = (bgc_struct) malloc (sizeof *BGCM);
#endif

    /* Read PIHM input files */
    ReadAlloc (simulation, pihm);
//#ifdef _FLUX_PIHM_
//    LSM_read (filename, LSM, cData);
//#endif
//#ifdef _BGC_
//    BGC_read (filename, BGCM, mData);
//#endif
//
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
//#ifdef _FLUX_PIHM_
//    LSM_initialize (filename, mData, cData, LSM);
//#endif
//#ifdef _RT_
//    chem_alloc(filename, mData, cData, chData, cData->StartTime/60);
//    /* Prepare chem output files */
//    InitialChemFile(filename, chData->NumBTC, chData->BTC_loc);
//#endif
//#ifdef _BGC_
//    BGC_init (filename, mData, LSM, BGCM);
//#endif
//
    MapOutput (simulation, pihm, outputdir);

    if (first_cycle == 1)
    {
        /* initialize output files and structures */
        InitOutputFile (pihm->prtctrl, pihm->ctrl.nprint, pihm->ctrl.ascii);
//#ifdef _FLUX_PIHM_
//        LSM_initialize_output (filename, mData, cData, LSM, output_dir);
//#endif
//#ifdef _BGC_
//        bgc_initialize_output (filename, mData, cData, BGCM, output_dir);
//#endif
        first_cycle = 0;
    }
//
//#ifdef _BGC_
//    if (BGCM->ctrl.spinup == 1)
//    {
//        bgc_spinup (filename, BGCM, mData, LSM);
//    }
//    else
//    {
//#endif
        if (verbose_mode)
            printf ("\n\nSolving ODE system ... \n\n");

        /* Set solver parameters */
        flag = CVodeSetFdata (cvode_mem, pihm);
        flag = CVodeSetInitStep (cvode_mem, (realtype) pihm->ctrl.initstep);
        flag = CVodeSetStabLimDet (cvode_mem, TRUE);
        flag = CVodeSetMaxStep (cvode_mem, (realtype) pihm->ctrl.maxstep);
        flag = CVodeMalloc (cvode_mem, f, (realtype ) pihm->ctrl.starttime,
            CV_Y, CV_SS, (realtype) pihm->ctrl.reltol, &pihm->ctrl.abstol);
        flag = CVSpgmr (cvode_mem, PREC_NONE, 0);

        /* set start time */
        t = pihm->ctrl.starttime;

        /* start solver in loops */
        for (i = 0; i < pihm->ctrl.nstep; i++)
        {
            /* inner loops to next output points with ET step size control */
            while (t < pihm->ctrl.tout[i + 1])
            {
                if (t + pihm->ctrl.etstep >= pihm->ctrl.tout[i + 1])
                {
                    nextptr = pihm->ctrl.tout[i + 1];
                }
                else
                {
                    nextptr = t + pihm->ctrl.etstep;
                }

                stepsize = nextptr - (realtype) t;

                pihm->dt = (double) stepsize;

                ApplyForcing (&pihm->forcing, t);

                if (t % pihm->ctrl.etstep == 0)
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

#ifdef _BGC_
                    BgcCoupling ((int) t, (int) cData->StartTime, mData, LSM, BGCM);
#endif
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
                    IntcpSnowET (t, (double) pihm->ctrl.etstep, pihm);
#endif
                }

                /* Added to adatpt to larger time step. */
                flag = CVodeSetMaxNumSteps(cvode_mem, (long int)(stepsize* 20));
                solvert = (realtype) t;
                flag = CVode (cvode_mem, nextptr, CV_Y, &solvert, CV_NORMAL);
                flag = CVodeGetCurrentTime(cvode_mem, &cvode_val);

                t = (int) solvert;
                rawtime = t;
                timestamp = gmtime (&rawtime);

                if (verbose_mode)
                    printf (" Time = %4.4d-%2.2d-%2.2d %2.2d:%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
                else if (rawtime % 3600 == 0)
                    printf (" Time = %4.4d-%2.2d-%2.2d %2.2d:%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);

                //summary (mData, CV_Y, t - StepSize, StepSize);
#ifdef _RT_
                /* PIHM-rt control file */
                fluxtrans(t / 60.0, StepSize / 60.0, mData, chData, CV_Y);
#endif
                //update (t, mData);
            }
//
//            /* Print outputs */
//            for (j = 0; j < cData->NumPrint; j++)
//                PrintData (cData->PCtrl[j], t, StepSize, cData->Ascii);
//#ifdef _FLUX_PIHM_
//            for (j = 0; j < LSM->NPRINT; j++)
//                PrintData (LSM->PCtrl[j], t, StepSize, cData->Ascii);
//#endif
//#ifdef _BGC_
//#endif
//#ifdef _RT_
//            /* PIHM-rt output routine */
//            PrintChem(filename, chData,t/60);
//#endif
//#ifdef _BGC_
        }
//#endif
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
//    /* Free memory */
//    N_VDestroy_Serial (CV_Y);
//
//    /* Free integrator memory */
//    CVodeFree (&cvode_mem);
//
//    free (rawtime);
//    free (filename);
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
        sprintf (str, "%2.2d%2.2d%2.2d%2.2d%2.2d", timestamp->tm_year + 1900 - 2000, timestamp->tm_mon + 1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
        sprintf (outputdir, "output/%s.%s/", project, str);
        printf ("\nOutput directory: %s\n", outputdir);
        mkdir (outputdir, 0755);
    }
}
