#include "pihm.h"

void ReadEnKF (enkf_struct ens)
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
    PIHMprintf (VL_VERBOSE, " Reading %s\n", fn);

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
    FindLine (enkf_file, "BOF", &lno, fn);
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

void ReadObs (int obs_time, char *fn, double *obs, double *obs_error)
{
    time_t          rawtime;
    struct tm      *timeinfo;
    double          temp1, temp2;
    FILE           *fid;
    char            cmdstr[MAXSTRING];
    int             match;
    int             lno = 0;

    timeinfo = (struct tm *)malloc (sizeof (struct tm));

    fid = fopen (fn, "r");
    CheckFile (fid, fn);
    PIHMprintf (VL_VERBOSE, " Reading %s\n", fn);

    FindLine (fid, "BOF", &lno, fn);

    while (1)
    {
        NextLine (fid, cmdstr, &lno);
        match = sscanf (cmdstr, "%d-%d-%d %d:%d %lf %lf",
            &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday,
            &timeinfo->tm_hour, &timeinfo->tm_min, &temp1, &temp2);
        timeinfo->tm_year = timeinfo->tm_year - 1900;
        timeinfo->tm_mon = timeinfo->tm_mon - 1;
        timeinfo->tm_sec = 0;
        rawtime = timegm (timeinfo);

        if (rawtime == obs_time)
        {
            *obs = temp1;
            *obs_error = temp2;
            break;
        }
        else if (strcasecmp (cmdstr, "EOF") == 0)
        {
            PIHMprintf (VL_ERROR,
                "\nError finding observation in %s.\n", fn);
            PIHMexit (EXIT_FAILURE);
        }
        else if (match != 7)
        {
            PIHMprintf (VL_ERROR,
                "Error reading observation in %s near Line %d.\n", fn, lno);
            PIHMexit (EXIT_FAILURE);
        }
    }

    fclose (fid);
    free (timeinfo);
}

void ReadFcst (enkf_struct ens, obs_struct obs, double *xf)
{
    int             i, j, k;
    int             ne;
    int             var_ind;
    double          xj;

    ne = ens->ne;

    for (i = 0; i < ne + 1; i++)
    {
        xf[i] = 0.0;
    }

    for (i = 0; i < ne; i++)
    {
        for (j = 0; j < obs.nlyr; j++)
        {
            var_ind = obs.var_ind[j];
            xj = 0.0;

            for (k = 0; k < ens->var[var_ind].dim; k++)
            {
                xj += obs.weight[k] * (ens->member[i].var[var_ind][k] *
                    obs.k[k][j] + obs.b[k][j]);
            }

            xf[i] += xj;
        }

        if (obs.type == RUNOFF_OBS)
        {
            xf[i] = log (xf[i] + 1.0);
        }

        xf[ne] += xf[i];
    }

    xf[ne] /= (double)ne;
}

void ReadVar (char *outputdir, enkf_struct ens, int obs_time)
{
    int             i, j, k;
    int             ii;
    int             ne;
    int             success = 0;
    int             length;
    double         *buffer;

    char            fn[MAXSTRING];
    FILE           *fid;

    ne = ens->ne;

    buffer =
        (double *)malloc ((ens->numele + ens->numriv + 1) * sizeof (double));

    for (i = 0; i < ne; i++)
    {
        for (k = 0; k < MAXVAR; k++)
        {
            if (ens->var[k].dim > 0)
            {
                sprintf (fn, "%s%s.%3.3d.%s.dat",
                    outputdir, project, i + 1, ens->var[k].name);
                fid = fopen (fn, "rb");
                CheckFile (fid, fn);

                fseek (fid, 0L, SEEK_END);

                length = (int)(ftell (fid) / (ens->var[k].dim + 1) / 8);

                rewind (fid);

                for (ii = 0; ii < length; ii++)
                {
                    fread (buffer, sizeof (double), ens->var[k].dim + 1, fid);

                    if ((int)buffer[0] == obs_time)
                    {
                        success = 1;

                        for (j = 0; j < ens->var[k].dim; j++)
                        {
                            ens->member[i].var[k][j] = buffer[j + 1];
                        }
                        break;
                    }
                }

                fclose (fid);

                if (success == 0)
                {
                    PIHMprintf (VL_ERROR,
                        "Error finding %s output for member %d.",
                        ens->var[k].name, i + 1);
                    PIHMexit (EXIT_FAILURE);
                }
            }
        }
    }

    free (buffer);
}
