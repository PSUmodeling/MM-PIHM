#include "pihm.h"
#include "optparse.h"


#if defined(_WIN32)
/* Do windows stuff here */
#include <direct.h>
#define pihm_mkdir(path) _mkdir((path))
#define pihm_access(path, amode) _access((path), (amode))
#define F_OK    0
#else 
#define pihm_mkdir(path) mkdir(path, 0755)
#define pihm_access(path, amode) access((path), (amode))
#endif


extern char            project[MAXSTRING];


void ParseCmdLineParam (int argc, char *argv[], int *spec_output_mode,
    char *outputdir)
{
	struct optparse options;
	optparse_init(&options, argv);
	struct optparse_long longopts[] = {
		{ "output", 'o', OPTPARSE_REQUIRED },
		{ "debug", 'd', OPTPARSE_NONE },
		{ "verbose", 'v', OPTPARSE_NONE },
		{ 0 }
	};

    int option;
	while ((option = optparse_long (&options, longopts, NULL)) != -1)
	{
	    switch (option)
	    {
	        case 'o':
	            /* Specify output directory */
				sprintf(outputdir, "output/%s/", options.optarg);
				*spec_output_mode = 1;
				printf ("Output directory is specified as \"%s\".\n", outputdir);
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
	             printf ("Option not recognisable %s %s \n", argv[0], options.errmsg);
	             exit(EXIT_FAILURE);
	            break;
	    }

        fflush (stdout);
    }

    if (options.optind >= argc)
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
        //Parse remaining arguments

        strcpy(project, optparse_arg(&options));
        if (project == NULL)
            PIHMexit(EXIT_FAILURE);
    }
}

void CreateOutputDir (char *outputdir, int spec_output_mode)
{
    time_t          rawtime;
    struct tm      *timestamp;
    char            str[11];

    if (0 == (pihm_mkdir(outputdir)))
	{
		PIHMprintf(VL_NORMAL, " Output directory was created.\n\n");
	}

    if (!spec_output_mode)
    {
        time (&rawtime);
        timestamp = localtime (&rawtime);

        /* Create output directory based on projectname and time */
        sprintf (str, "%2.2d%2.2d%2.2d%2.2d%2.2d",
            timestamp->tm_year + 1900 - 2000, timestamp->tm_mon + 1,
            timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
        sprintf (outputdir, "output/%s.%s/", project, str);

        PIHMprintf (VL_NORMAL, "\nOutput directory: %s\n", outputdir);
    }
    pihm_mkdir(outputdir);
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
    if (pihm_access(source_file, F_OK) != -1)
    {
		sprintf (system_cmd, "cp %s ./%s/%s.para.bak", source_file, outputdir, project);
		system(system_cmd);
	}
    sprintf (source_file, "input/%s/%s.calib", project, simulation);
    if (pihm_access(source_file, F_OK) != -1)
    {
        sprintf (system_cmd, "cp %s ./%s/%s.calib.bak", source_file,
            outputdir, simulation);
    }
    sprintf (source_file, "input/%s/%s.init", project, simulation);
    if (pihm_access(source_file, F_OK) != -1)
    {

        sprintf (system_cmd, "cp %s ./%s/%s.init.bak", source_file, outputdir,
            simulation);
        system (system_cmd);
    }
}

void SetCVodeParam (pihm_struct pihm, void *cvode_mem, N_Vector CV_Y, N_Vector abstol)
{
    int             flag;

	flag = CVodeSetUserData(cvode_mem, pihm);
    flag = CVodeSetInitStep (cvode_mem, (realtype) pihm->ctrl.initstep);
    flag = CVodeSetStabLimDet (cvode_mem, TRUE);
    flag = CVodeSetMaxStep (cvode_mem, (realtype) pihm->ctrl.maxstep);
	//flag = CVodeInit(cvode_mem, Hydrol, (realtype)pihm->ctrl.starttime, CV_Y);
    flag = CVodeInit(cvode_mem, Hydrol, (realtype)0.0, CV_Y);
	flag = CVodeSVtolerances(cvode_mem, (realtype)pihm->ctrl.reltol, abstol);
    flag = CVSpgmr (cvode_mem, PREC_NONE, 0);
}

void SolveCVode(int starttime, int *t, int nextptr, int stepsize, double cputime_dt, double cputime, void *cvode_mem,
	N_Vector CV_Y, int cvode_perf, char *simulation, char *outputdir)
{
	realtype        solvert;
	realtype        cvode_val;
    realtype        tout = (realtype)(nextptr - starttime);
	struct tm      *timestamp;
	time_t          rawtime;
	int             flag;
	static double	dtime = 0;
	char			Perfname[50];
	FILE            *Perf; /* Performance file */
	
	flag = CVodeSetMaxNumSteps(cvode_mem, (long int)(stepsize * 20));
	flag = CVodeSetStopTime(cvode_mem, tout);
	flag = CVode(cvode_mem, (realtype)nextptr, CV_Y, &solvert,
		CV_NORMAL);

	flag = CVodeGetCurrentTime(cvode_mem, &cvode_val);

    *t = (int)(solvert + starttime);
	rawtime = (time_t)(*t);
	timestamp = gmtime(&rawtime);

    if (flag != CV_SUCCESS && flag != CV_TSTOP_RETURN && flag != CV_ROOT_RETURN)
    {
        PIHMprintf(VL_NORMAL, " FAILED CVode at %4.4d-%2.2d-%2.2d %2.2d:%2.2d (%d) returned [%d]\n",
            timestamp->tm_year + 1900, timestamp->tm_mon + 1,
            timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min, *t, nextptr, flag);
    }

    dtime = dtime + cputime_dt;
	if (verbose_mode)
	{
		PIHMprintf(VL_VERBOSE, " Step = %4.4d-%2.2d-%2.2d %2.2d:%2.2d (%d) Time %d \n",
			timestamp->tm_year + 1900, timestamp->tm_mon + 1,
			timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min, *t, nextptr);
	}
#ifndef _ENKF_
	else if (rawtime % 3600 == 0)
	{
		PIHMprintf(VL_NORMAL, " Step = %4.4d-%2.2d-%2.2d %2.2d:%2.2d \n",
			timestamp->tm_year + 1900, timestamp->tm_mon + 1,
			timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
	}
#endif
	sprintf(Perfname, "%s%s_Performance.txt", outputdir, simulation);
	Perf = fopen(Perfname, "a");
	CheckFile(Perf, Perfname);
	if (rawtime % 3600 == 0)
	{
		fprintf(Perf, " Step = %4.4d-%2.2d-%2.2d %2.2d:%2.2d CPU time =  %f %f \n",
			timestamp->tm_year + 1900, timestamp->tm_mon + 1,
			timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min, dtime, cputime);
		dtime = 0.;
	}
	fclose(Perf);
	if (cvode_perf)
	{
		/* Print some CVODE statistics */
		PrintStats(cvode_mem);
	}

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
