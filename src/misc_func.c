#include "pihm.h"

void ParseCmdLineParam (int argc, char *argv[], int *spec_output_mode,
    char *outputdir)
{
    int             c;

    while ((c = getopt (argc, argv, "o:dvl")) != -1)
    {
        switch (c)
        {
            case 'o':
                /* Specify output directory */
                sprintf (outputdir, "output/%s/", optarg);
                *spec_output_mode = 1;
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

        fflush (stdout);
    }

    if (optind >= argc)
    {
        fprintf (stderr, "Error:You must specify the name of project!\n");
        fprintf (stderr,
            "Usage: ./pihm [-o output_dir] [-d] [-v] <project name>\n");
        fprintf (stderr, "\t-o Specify output directory.\n");
        fprintf (stderr, "\t-v Verbose mode\n");
        fprintf (stderr, "\t-d Debug mode\n");
        PIHMexit (EXIT_FAILURE);
    }
    else
    {
        strcpy (project, argv[optind]);
    }
}

void CreateOutputDir (char *outputdir, int spec_output_mode)
{
    time_t          rawtime;
    struct tm      *timestamp;
    char            str[11];

    if (0 == (mkdir ("output", 0755)))
    {
        PIHMprintf (VL_NORMAL, " Output directory was created.\n\n");
    }

    if (!spec_output_mode)
    {
        time (&rawtime);
        timestamp = localtime (&rawtime);

        /* Create output directory based on projectname and time */
        strftime (str, 11, "%y%m%d%H%M", timestamp);
        sprintf (outputdir, "output/%s.%s/", project, str);

        PIHMprintf (VL_NORMAL, "\nOutput directory: %s\n", outputdir);
    }
    else
    {
        
        PIHMprintf
            (VL_NORMAL, "Output directory is specified as \"%s\".\n",
            outputdir);
    }

    mkdir (outputdir, 0755);
}

void BKInput (char *simulation, char *outputdir)
{
    char            project[MAXSTRING];
    char           *token;
    char            tempname[MAXSTRING];
    char            system_cmd[MAXSTRING];
    char            source_file[MAXSTRING];

    strcpy (tempname, simulation);
    if (strstr (tempname, ".") != 0)
    {
        token = strtok (tempname, ".");
        strcpy (project, token);
    }
    else
    {
        strcpy (project, simulation);
    }

    /* Save input files into output directory */
    sprintf (source_file, "input/%s/%s.para", project, project);
    if (access (source_file, F_OK) != -1)
    {
        sprintf (system_cmd, "cp %s ./%s/%s.para.bak", source_file, outputdir,
            project);
        system (system_cmd);
    }
    sprintf (source_file, "input/%s/%s.calib", project, simulation);
    if (access (source_file, F_OK) != -1)
    {
        sprintf (system_cmd, "cp %s ./%s/%s.calib.bak", source_file,
            outputdir, simulation);
        system (system_cmd);
    }
    sprintf (source_file, "input/%s/%s.init", project, simulation);
    if (access (source_file, F_OK) != -1)
    {
        sprintf (system_cmd, "cp %s ./%s/%s.init.bak", source_file, outputdir,
            simulation);
        system (system_cmd);
    }
}

void SetCVodeParam (pihm_struct pihm, void *cvode_mem, N_Vector CV_Y)
{
    int             flag;

    flag = CVodeInit (cvode_mem, Hydrol, (realtype)pihm->ctrl.starttime,
        CV_Y);
    flag = CVodeSStolerances (cvode_mem,(realtype) pihm->ctrl.reltol,
        pihm->ctrl.abstol);
    flag = CVodeSetUserData (cvode_mem, pihm);
    flag = CVodeSetInitStep (cvode_mem, (realtype) pihm->ctrl.initstep);
    flag = CVodeSetStabLimDet (cvode_mem, TRUE);
    flag = CVodeSetMaxStep (cvode_mem, (realtype) pihm->ctrl.maxstep);
    //flag = CVodeMalloc (cvode_mem, Hydrol, (realtype) pihm->ctrl.starttime,
    //    CV_Y, CV_SS, ;
    flag = CVSpgmr (cvode_mem, PREC_NONE, 0);
}

void SolveCVode (int *t, int nextptr, int stepsize, void *cvode_mem,
    N_Vector CV_Y)
{
    realtype        solvert;
    realtype        cvode_val;
    char            timestr[MAXSTRING];
    struct tm      *timestamp;
    time_t          rawtime;
    int             flag;

    solvert = (realtype) (*t);

    flag = CVodeSetMaxNumSteps (cvode_mem, (long int)(stepsize * 20));
    flag = CVodeSetStopTime (cvode_mem, (realtype) nextptr);
    flag = CVode (cvode_mem, (realtype) nextptr, CV_Y, &solvert,
        CV_NORMAL);
    flag = CVodeGetCurrentTime (cvode_mem, &cvode_val);

    *t = (int)solvert;
    rawtime = (time_t) (*t);
    timestamp = gmtime (&rawtime);
    strftime (timestr, 17, "%Y-%m-%d %H:%M", timestamp);

    if (debug_mode)
    {
        PIHMprintf (VL_VERBOSE, " Step = %s (%d)\n", timestr, *t);
    }
#ifndef _ENKF_
    else if (rawtime % 3600 == 0)
    {
        PIHMprintf (VL_NORMAL, " Step = %s\n", timestr);
    }
#endif
}

void _PIHMexit (const char *fn, int lineno, const char *func, int error)
{
#ifdef _ENKF_
    int             id;
    int             ierr;

    ierr = MPI_Comm_rank (MPI_COMM_WORLD, &id);
    PIHMprintf (VL_ERROR, "Error in %s\n", func);
    PIHMprintf (VL_ERROR, "Exiting from Node %d\n", id);
    MPI_Abort (MPI_COMM_WORLD, error);
#else
    PIHMprintf (VL_ERROR, "\n");
    PIHMprintf (VL_ERROR, "Exiting from %s", func);
    if (debug_mode)
    {
        PIHMprintf (VL_ERROR, " (%s, Line %d)", fn, lineno);
    }
    PIHMprintf (VL_ERROR, "...\n\n");

    exit (error);
#endif
}
