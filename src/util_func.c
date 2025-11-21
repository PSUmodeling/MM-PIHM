#include "pihm.h"
#include "optparse.h"

void ParseCmdLineParam(int argc, char *argv[], char outputdir[])
{
    int             option;
    struct optparse options;
    struct optparse_long longopts[] = {
        {"append",     'a', OPTPARSE_NONE},
        {"brief",      'b', OPTPARSE_NONE},
        {"correction", 'c', OPTPARSE_NONE},
        {"debug",      'd', OPTPARSE_NONE},
        {"fixed",      'f', OPTPARSE_NONE},
        {"output",     'o', OPTPARSE_REQUIRED},
        {"silent",     's', OPTPARSE_NONE},
        {"version",    'V', OPTPARSE_NONE},
        {"verbose",    'v', OPTPARSE_NONE},
        {0, 0, 0}
    };

    optparse_init(&options, argv);

    while ((option = optparse_long(&options, longopts, NULL)) != -1)
    {
        switch (option)
        {
            case 'o':
                // Specify output directory
                sprintf(outputdir, "output/%s/", options.optarg);
                break;
            case 'c':
                // Surface elevation correction mode
                corr_mode = 1;
                break;
            case 'd':
                // Debug mode
                debug_mode = 1;
                break;
            case 'f':
                // Fixed length spin-up
                fixed_length = 1;
                break;
            case 'v':
                // Verbose mode
                verbose_mode = VL_VERBOSE;
                break;
            case 'b':
                // Brief mode
                verbose_mode = VL_BRIEF;
                break;
            case 's':
                // Silent mode
                verbose_mode = VL_SILENT;
                break;
            case 'a':
                // Append mode
                append_mode = 1;
                break;
            case 'V':
                // Print version number
                printf("MM-PIHM Version %s\n", VERSION);
#if defined(_CYCLES_)
                printf("Cycles Version %s coupled\n", CYCLES_VERSION);
#endif
#if defined(_OPENMP)
                printf("Paralleled with OpenMP\n");
#endif
                pihm_exit(EXIT_SUCCESS);
                break;
            case '?':
                pihm_printf(VL_ERROR, "Option not recognizable %s\n", options.errmsg);
                pihm_exit(EXIT_FAILURE);
                break;
            default:
                break;
        }

        fflush(stdout);
    }

    if (options.optind >= argc)
    {
        pihm_printf(VL_ERROR, "Error:You must specify the name of project!\n"
            "Usage: ./pihm [-o output_dir] [-c] [-d] [-t] [-v] [-V]"
            " <project name>\n"
            "    -o Specify output directory\n"
            "    -b Brief mode\n"
            "    -c Correct surface elevation\n"
            "    -d Debug mode\n"
            "    -V Version number\n"
            "    -v Verbose mode\n");
        pihm_exit(EXIT_FAILURE);
    }
    else
    {
        // Parse remaining arguments
        strcpy(project, optparse_arg(&options));
    }
}

void CreateOutputDir(char outputdir[])
{
    time_t          rawtime;
    struct tm      *timestamp;
    char            str[11];
    char            icdir[MAXSTRING];
    char            proj[MAXSTRING];
    char           *token;

    strcpy(proj, project);
    if (strstr(proj, ".") != 0)
    {
        token = strtok(proj, ".");
        strcpy(proj, token);
    }
    else
    {
        strcpy(proj, project);
    }

    if (0 == (pihm_mkdir("output")))
    {
        pihm_printf(VL_NORMAL, "Output directory was created.\n\n");
    }

    if (outputdir[0] == '\0')
    {
        // Create default output directory name based on project and time
        time(&rawtime);
        timestamp = localtime(&rawtime);
        strftime(str, 11, "%y%m%d%H%M", timestamp);
        sprintf(outputdir, "output/%s.%s/", proj, str);
    }

    if (pihm_mkdir(outputdir) != 0)
    {
        if (errno != EEXIST)
        {
            pihm_printf(VL_ERROR, "Error creating output directory %s\n", outputdir);
            pihm_exit(EXIT_FAILURE);
        }
        else
        {
            pihm_printf(VL_BRIEF, "Output directory %s already exists. Overwriting.\n", outputdir);
        }
    }
    else
    {
        pihm_printf(VL_BRIEF, "Output directory %s was created.\n", outputdir);
    }

    sprintf(icdir, "%srestart/", outputdir);
    if (pihm_mkdir(icdir) != 0 && errno != EEXIST)
    {
        pihm_printf(VL_ERROR, "Error creating restart directory %s\n", outputdir);
        pihm_exit(EXIT_FAILURE);
    }
}

void BackupInput(const char outputdir[], const filename_struct *filename)
{
    char            system_cmd[MAXSTRING];

    // Save input files into output directory
    if (pihm_access(filename->para, F_OK) != -1)
    {
        sprintf(system_cmd, "cp %s ./%s/%s.para.bak", filename->para, outputdir, project);
        system(system_cmd);
    }
    if (pihm_access(filename->calib, F_OK) != -1)
    {
        sprintf(system_cmd, "cp %s ./%s/%s.calib.bak", filename->calib, outputdir, project);
        system(system_cmd);
    }
    if (pihm_access(filename->ic, F_OK) != -1)
    {
        sprintf(system_cmd, "cp %s ./%s/%s.ic.bak", filename->ic, outputdir, project);
        system(system_cmd);
    }
#if defined(_BGC_)
    if (pihm_access(filename->bgcic, F_OK) != -1)
    {
        sprintf(system_cmd, "cp %s ./%s/%s.bgcic.bak", filename->bgcic, outputdir, project);
        system(system_cmd);
    }
#endif
#if defined(_CYCLES_)
    if (pihm_access(filename->cyclesic, F_OK) != -1)
    {
        sprintf(system_cmd, "cp %s ./%s/%s.cyclesic.bak", filename->cyclesic, outputdir, project);
        system(system_cmd);
    }
#endif
}

void CheckCVodeFlag(int cv_flag)
{
    if (cv_flag < 0)
    {
        pihm_printf(VL_ERROR, "CVODE error %s\n", CVodeGetReturnFlagName(cv_flag));
        pihm_exit(EXIT_FAILURE);
    }
}

// This function is the same as in Cycles to ensure successful coupling. If this function is modified, the corresponding
// function in Cycles should be modified, too.
int roundi(double x)
{
    return (int)((x < 0.0) ? x - 0.5 : x + 0.5);
}
