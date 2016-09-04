#include "pihm.h"

void EnKFRead (enkf_struct ens)
{
    char            fn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    FILE           *enkf_file;
    int             i;
    int             n;
    ctrl_struct     ctrl;
    int             lno = 0;

    /*
     * Open .enkf file
     */
    sprintf (fn, "input/%s/%s.enkf", project, project);
    enkf_file = fopen (fn, "r");
    CheckFile (enkf_file, fn);
    if (verbose_mode)
    {
        printf (" Reading %s\n", fn);
    }

    /*
     * Read .enkf file
     */
    FindLine (enkf_file, "BOF", &lno, fn);

    NextLine (enkf_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "NUM_ENSEMBLE_MEMBER", &ens->ne, 'i', fn, lno);

    NextLine (enkf_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "ASSIMILATION_INTERVAL", &ens->interval, 'i', fn,
        lno);
    ens->interval *= 3600;

    NextLine (enkf_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "START_MODE", &ens->start_mode, 'i', fn, lno);

    NextLine (enkf_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "INFLATION_WEIGHT", &ens->weight, 'd', fn, lno);

    NextLine (enkf_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "ASSIMILATION_END_TIME", &ens->end_time, 't', fn,
        lno);

    FindLine (enkf_file, "PARAMETER", &lno, fn);
    n = CountLine (enkf_file, cmdstr, 1, "NUM_OBS");

    /* Rewind to read */
    FindLine (enkf_file, "PARAMETER", &lno, fn);

    /* Start reading EnKF file */
    for (i = 0; i < n; i++)
    {
        NextLine (enkf_file, cmdstr, &lno);
        sscanf (cmdstr, "%s %d %d %lf %lf %lf %lf %d",
            ens->param[i].name,
            &ens->param[i].perturb, &ens->param[i].update,
            &ens->param[i].perturb_min, &ens->param[i].perturb_max,
            &ens->param[i].min, &ens->param[i].max, &ens->param[i].type);
    }

    NextLine (enkf_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "NUM_OBS", &ens->nobs, 'i', fn, lno);

    if (ens->nobs > 0)
    {
        ens->obs = (obs_struct *)malloc (ens->nobs * sizeof (obs_struct));

        for (i = 0; i < ens->nobs; i++)
        {
            NextLine (enkf_file, cmdstr, &lno);
            ReadKeyword (cmdstr, "OBS_TYPE", ens->obs[i].name, 's', fn, lno);

            NextLine (enkf_file, cmdstr, &lno);
            ReadKeyword (cmdstr, "OBS_FILE", ens->obs[i].fn, 's', fn, lno);


            NextLine (enkf_file, cmdstr, &lno);
            ReadKeyword (cmdstr, "OBS_LOCATION_X", &ens->obs[i].x, 'd', fn,
                lno);

            NextLine (enkf_file, cmdstr, &lno);
            ReadKeyword (cmdstr, "OBS_LOCATION_Y", &ens->obs[i].y, 'd', fn,
                lno);

            NextLine (enkf_file, cmdstr, &lno);
            ReadKeyword (cmdstr, "OBS_RADIUS", &ens->obs[i].rad, 'd', fn,
                lno);

            NextLine (enkf_file, cmdstr, &lno);
            ReadKeyword (cmdstr, "OBS_DEPTH", &ens->obs[i].depth, 'd', fn,
                lno);
        }
    }

    fclose (enkf_file);

    /*
     * Read .para to obtain the first cycle's start and end times
     */
    sprintf (fn, "input/%s/%s.para", project, project);
    ReadPara (fn, &ctrl);
    ens->mbr_start_mode = ctrl.init_type;
    ens->cycle_start_time = ctrl.starttime;
    ens->cycle_end_time = ctrl.endtime;
}
