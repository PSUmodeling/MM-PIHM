#include "pihm.h"

void ParseCmdLineParam (int argc, char *argv[], char *outputdir)
{
    int             c;

    outputdir[0] = '\0';

    while ((c = getopt (argc, argv, "o:cdvV")) != -1)
    {
        switch (c)
        {
            case 'o':
                /* Specify output directory */
                sprintf (outputdir, "output/%s/", optarg);
                break;
            case 'c':
                /* Surface elevatoin correction mode */
                corr_mode = 1;
                printf ("Surface elevation correction mode turned on.\n");
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
            case 'V':
                /* Print version number */
                printf ("\nMM-PIHM Version %s.\n", VERSION);
                PIHMexit (EXIT_SUCCESS);
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
            "Usage: ./pihm [-o output_dir] [-c] [-d] [-v] <project name>\n");
        fprintf (stderr, "\t-o Specify output directory\n");
        fprintf (stderr, "\t-c Correct surface elevation\n");
        fprintf (stderr, "\t-v Verbose mode\n");
        fprintf (stderr, "\t-d Debug mode\n");
        PIHMexit (EXIT_FAILURE);
    }
    else
    {
        strcpy (project, argv[optind]);
    }
}

void CreateOutputDir (char *outputdir)
{
    time_t          rawtime;
    struct tm      *timestamp;
    char            str[11];

    if (0 == (mkdir ("output", 0755)))
    {
        PIHMprintf (VL_NORMAL, " Output directory was created.\n\n");
    }

    if (outputdir[0] == '\0')
    {
        /* Create default output directory name based on project and time */
        time (&rawtime);
        timestamp = localtime (&rawtime);
        strftime (str, 11, "%y%m%d%H%M", timestamp);
        sprintf (outputdir, "output/%s.%s/", project, str);
    }

    PIHMprintf (VL_NORMAL, "\nOutput directory: %s\n", outputdir);

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

void _PIHMexit (const char *fn, int lineno, const char *func, int error)
{
    PIHMprintf (VL_ERROR, "\n");
    PIHMprintf (VL_ERROR, "Exiting from %s", func);
    if (debug_mode)
    {
        PIHMprintf (VL_ERROR, " (%s, Line %d)", fn, lineno);
    }
    PIHMprintf (VL_ERROR, "...\n\n");

    exit (error);
}

