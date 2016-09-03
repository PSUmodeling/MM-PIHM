#include "pihm.h"

void EnKFCore (double *xa, double obs, double obs_err, double *xf, int ne)
{
    int             i;
    double          y_hxm;
    double          km;
    double         *hxa;
    double         *xp;
    double          var;
    double          fac;
    double          s;
    double          beta;
    double          fac2;

    hxa = (double *)malloc (ne * sizeof (double));
    xp = (double *)malloc (ne * sizeof (double));

    for (i = 0; i < ne; i++)
    {
        xa[i] -= xa[ne];
    }

    y_hxm = obs - xf[ne];

    var = 0;

    for (i = 0; i < ne; i++)
    {
        hxa[i] = xf[i] - xf[ne];
        var += hxa[i] * hxa[i];
    }

    fac = 1.0 / ((double)ne - 1.0);
    s = fac * var + obs_err * obs_err;
    beta = 1.0 / (1.0 + sqrt ((obs_err * obs_err) / s));
    fac2 = fac / s;

    km = 0.0;

    for (i = 0; i < ne; i++)
    {
        xp[i] = xa[i];
        km += xp[i] * hxa[i];
    }

    km *= fac2;

    xa[ne] += km * y_hxm;

    km *= beta;

    for (i = 0; i < ne; i++)
    {
        xa[i] += km * (0.0 - hxa[i]);
        xa[i] += xa[ne];
    }

    free (hxa);
    free (xp);
}

void EnKF (enkf_struct ens, int obs_time, char *outputdir)
{
    double         *xf;
    int             i, j, k;
    int             ne;
    char            obsfn[MAXSTRING];
    char            obsin_fn[MAXSTRING];
    FILE           *obsfile;
    double          obs;
    double          obs_error;
    struct tm      *timestamp;
    time_t          rawtime;

    enkf_struct     ens0;

    ne = ens->ne;

    ens0 = (enkf_struct)malloc (sizeof *ens0);
    ens0->ne = ne;
    ens0->member = (ensmbr_struct *)malloc (ne * sizeof (ensmbr_struct));

    for (i = 0; i < ne; i++)
    {
        for (j = 0; j < MAXVAR; j++)
        {
            if (ens->var[j].dim > 0)
            {
                ens0->member[i].var[j] =
                    (double *)malloc (ens->var[j].dim * sizeof (double));
            }
        }
    }

    xf = (double *)malloc (sizeof (double) * (ne + 1));
    rawtime = obs_time;
    timestamp = gmtime (&rawtime);

    printf ("\nStarting EnKF ... \n");

    /*
     * Copy prior from ens to ens0
     */
    for (i = 0; i < ne; i++)
    {
        for (j = 0; j < MAXPARAM; j++)
        {
            ens0->member[i].param[j] = ens->member[i].param[j];
        }
        for (j = 0; j < MAXVAR; j++)
        {
            if (ens->var[j].dim > 0)
            {
                for (k = 0; k < ens->var[j].dim; k++)
                {
                    ens0->member[i].var[j][k] = ens->member[i].var[j][k];
                }
            }
        }
    }

    //for (j=0; j<ne; j++)
    //{
    //    ens->TotalWaterVolume[j] = 0;
    //    for (k = 0; k<DS->NumEle; k++)
    //    {
    //        ens->TotalWaterVolume[j] = ens->TotalWaterVolume[j] + DS->Ele[k].area*(ens->member[j].variable[0][k] + ens->member[j].variable[2][k]*DS->Ele[k].Porosity + ens->member[j].variable[3][k] + ens->member[j].variable[8][k] + ens->member[j].variable[7][k] + ens->member[j].variable[22][k]/Lv/1000*3600);
    //    }
    //    for (k = 0; k<DS->NumRiv; k++)
    //    {
    //        ens->TotalWaterVolume[j] = ens->TotalWaterVolume[j] + DS->Riv[k].coeff*DS->Riv[k].Length*(ens->member[j].variable[1][k] + ens->member[j].variable[8][k+DS->NumEle] + ens->member[j].variable[3][k+DS->NumEle]);
    //    }
    //    ens->TotalWaterVolume[j] = ens->TotalWaterVolume[j] + ens->member[j].variable[10][DS->NumRiv-1]/24;
    //    En0->TotalWaterVolume[j] = ens->TotalWaterVolume[j];
    //}

    if (ens->nobs > 0)
    {
        sprintf (obsfn, "%sobs.dat", outputdir);
        obsfile = fopen (obsfn, "a");
        fprintf (obsfile, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
            timestamp->tm_year + 1900, timestamp->tm_mon + 1,
            timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);

        for (i = 0; i < ens->nobs; i++)
        {
            printf ("\n*****%s******\n", ens->obs[i].name);

            /* 
             * Read observations
             */
            sprintf (obsin_fn, "input/%s/%s", project, ens->obs[i].fn);
            ReadObs (obs_time, obsin_fn, &obs, &obs_error);

            printf ("observation = %lf\n", obs);
            printf ("error = %lf\n", obs_error);

            /*
             * Read ensemble forecasts
             */
            ReadFcst (ens, ens->obs[i], xf);

            /* 
             * Prepare forecast vectors
             */
            printf ("prediction = ");
            for (j = 0; j < ne; j++)
            {
                printf ("%f\t", xf[j]);
            }
            printf ("mean: %f\n", xf[ne]);

            /* 
             * Write observations to files
             */
            fprintf (obsfile, "\t%lf", obs);

            /* 
             * Update analysis
             */
            UpdAnlys (ens, obs, obs_error, xf);
        }

        fprintf (obsfile, "\n");
        fflush (obsfile);
        fclose (obsfile);

        /* 
         * Covariance inflation
         */
        CovInflt (ens, ens0);
    }

    WriteEnKFOut (project, ens, outputdir, obs_time);

    for (i = 0; i < ne; i++)
    {
        for (j = 0; j < MAXVAR; j++)
        {
            if (ens->var[j].dim > 0)
            {
                free (ens0->member[i].var[j]);
            }
        }
    }

    free (ens0->member);
    free (ens0);
    free (xf);
}

void UpdAnlys (enkf_struct ens, double obs, double obs_error, double *xf)
{
    int             i, j, k;
    int             ne;
    double         *xa;
    double        **x;

    ne = ens->ne;

    xa = (double *)malloc (sizeof (double) * (ne + 1));

    x = (double **)malloc (ne * sizeof (double));

    /*
     * Optimize parameters using EnKF
     */
    for (i = 0; i < MAXPARAM; i++)
    {
        if (ens->param[i].update == 1 && ens->update_param == 1)
        {
            /* Make pointers point to the parameter that needs to be
             * updated */
            for (j = 0; j < ne; j++)
            {
                x[j] = &ens->member[j].param[i];
            }

            /* Take log if the parameter is conductivity or Czil */
            if (ens->param[i].type == LOG_TYPE)
            {
                for (j = 0; j < ne; j++)
                {
                    *x[j] = log10 (*x[j]);
                }
            }

            xa[ne] = 0.0;

            for (j = 0; j < ne; j++)
            {
                xa[j] = *x[j];
                xa[ne] += xa[j];
            }

            xa[ne] /= (double)ne;

            EnKFCore (xa, obs, obs_error, xf, ne);

            for (j = 0; j < ne; j++)
            {
                *x[j] = xa[j];
                if (ens->param[i].type == LOG_TYPE)
                {
                    *x[j] = pow (10.0, *x[j]);
                }
            }
        }
    }

    /*
     * Optimize model variables
     */
    if (ens->update_var == 1)
    {
        for (i = 0; i < MAXVAR; i++)
        {
            if (ens->var[i].dim > 0)
            {
                for (k = 0; k < ens->var[i].dim; k++)
                {
                    for (j = 0; j < ne; j++)
                    {
                        x[j] = &ens->member[j].var[i][k];
                    }

                    xa[ne] = 0.0;

                    for (j = 0; j < ne; j++)
                    {
                        xa[j] = *x[j];
                        xa[ne] += xa[j];
                    }
                    xa[ne] /= (double)ne;

                    EnKFCore (xa, obs, obs_error, xf, ne);

                    for (j = 0; j < ne; j++)
                    {
                        *x[j] = xa[j];
                    }
                }
            }
        }
    }


    free (xa);
    free (x);
}

void CovInflt (enkf_struct ens, enkf_struct ens0)
{
    int             i, j, k;
    int             ne;

    double         *xp0, *xp;

    double        **x, *x0;
    double          average0, average;
    double          param_min, param_max;
    double          param_std;
    double          c1, c2, c;

    //double TotalWaterVolume[ens->ne];
    //double MassCoeff[ens->ne];

    ne = ens->ne;

    xp0 = (double *)malloc (sizeof (double) * ne);
    xp = (double *)malloc (sizeof (double) * ne);

    x = (double **)malloc (ne * sizeof (double));
    x0 = (double *)malloc (ne * sizeof (double));

    printf ("\n*****Parameters********\n");

    for (i = 0; i < MAXPARAM; i++)
    {
        if (ens->param[i].update == 1 && ens->update_param == 1)
        {
            printf ("%s\n", ens->param[i].name);
            /* Make pointers point to the parameter that needs to be updated */
            for (j = 0; j < ne; j++)
            {
                x[j] = &ens->member[j].param[i];
                x0[j] = ens0->member[j].param[i];
            }

            /* Take log if the parameter is conductivity or Czil */
            if (ens->param[i].type == LOG_TYPE)
            {
                for (j = 0; j < ne; j++)
                {
                    *x[j] = log10 (*x[j]);
                    x0[j] = log10 (x0[j]);
                }
            }

            average0 = 0.0;
            average = 0.0;

            for (j = 0; j < ne; j++)
            {
                average0 += x0[j];
                average += *x[j];
            }

            average /= (double)ne;
            average0 /= (double)ne;

            for (j = 0; j < ne; j++)
            {
                xp0[j] = x0[j] - average0;
                xp[j] = *x[j] - average;
            }

            if (average < ens->param[i].max - 0.25 * ens->param[i].init_std &&
                average > ens->param[i].min + 0.25 * ens->param[i].init_std)
            {
                /* Calculate new perturbations and standard deviation */
                param_std = 0.0;

                for (j = 0; j < ne; j++)
                {
                    xp[j] =
                        (1.0 - ens->weight) * xp[j] + ens->weight * xp0[j];
                    param_std += xp[j] * xp[j];
                }
                param_std = sqrt (param_std / ((double)ne - 1.0));

                /* Coavariance inflation */
                param_min = 999.0;
                param_max = -999.0;

                for (j = 0; j < ne; j++)
                {
                    if (param_std < 0.25 * ens->param[i].init_std)
                    {
                        xp[j] =
                            0.25 * ens->param[i].init_std / param_std * xp[j];
                    }
                    param_min =
                        (param_min <
                        average + xp[j]) ? param_min : (average + xp[j]);
                    param_max =
                        (param_max >
                        average + xp[j]) ? param_max : (average + xp[j]);
                }

                c1 = (ens->param[i].max - 1.0E-6 - average) / (param_max -
                    average);
                c2 = (average - ens->param[i].min - 1.0E-6) / (average -
                    param_min);
                c = (c1 < c2) ? c1 : c2;
                c = (c < 1.0) ? c : 1.0;

                for (j = 0; j < ne; j++)
                {
                    *x[j] = average + c * xp[j];
                    if (ens->param[i].type == LOG_TYPE)
                    {
                        *x[j] = pow (10.0, *x[j]);
                    }
                    printf ("%lf\t", *x[j]);
                }
                if (ens->param[i].type == LOG_TYPE)
                {
                    average = pow (10.0, average);
                }
                printf ("mean: %lf\n", average);
            }
            else
            {
                printf
                    ("EnKF analysis %lf is out of range. Parameter is not updated\n",
                    average);
                average = average0;
                for (j = 0; j < ne; j++)
                {
                    *x[j] = x0[j];
                    if (ens->param[i].type == LOG_TYPE)
                    {
                        *x[j] = pow (10.0, x0[j]);
                    }
                }

                if (ens->param[i].type == LOG_TYPE)
                {
                    average = pow (10.0, average);
                }
            }
        }
    }

    if (ens->update_var == 1)
    {
        for (i = 0; i < MAXVAR; i++)
        {
            if (ens->var[i].dim > 0)
            {
                for (k = 0; k < ens->var[i].dim; k++)
                {
                    for (j = 0; j < ne; j++)
                    {
                        x[j] = &ens->member[j].var[i][k];
                        x0[j] = ens0->member[j].var[i][k];
                    }

                    average = 0.0;
                    average0 = 0.0;

                    for (j = 0; j < ne; j++)
                    {
                        average += *x[j];
                        average0 += x0[j];
                    }
                    average /= (double)ne;
                    average0 /= (double)ne;

                    for (j = 0; j < ne; j++)
                    {
                        xp0[j] = x0[j] - average0;
                        xp[j] = *x[j] - average;
                        *x[j] =
                            average + (1.0 - ens->weight) * xp[j] +
                            ens->weight * xp0[j];
                        if ((int)ens->var[i].max != BADVAL)
                        {
                            *x[j] =
                                (*x[j] <
                                ens->var[i].max) ? *x[j] : ens->var[i].max;
                        }

                        if ((int)ens->var[i].min != BADVAL)
                        {
                            *x[j] =
                                (*x[j] >
                                ens->var[i].min) ? *x[j] : ens->var[i].min;
                        }
                    }
                }
            }
        }
    }

    free (xp0);
    free (xp);
    free (x);
    free (x0);

    printf ("Inflation done!\n");
}

void ReadObs (int obs_time, char *fn, double *obs, double *obs_error)
{
    time_t          rawtime;
    struct tm      *timeinfo;
    double          temp1, temp2;
    FILE           *fid;
    char            cmdstr[MAXSTRING];
    int             match;

    timeinfo = (struct tm *)malloc (sizeof (struct tm));

    fid = fopen (fn, "r");
    CheckFile (fid, fn);

    FindLine (fid, "BOF");

    while (1)
    {
        NextLine (fid, cmdstr);
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
            printf ("\nFATAL ERROR: No observation availablein %s!\n", fn);
            PIHMExit (EXIT_FAILURE);
        }
        else if (match != 7)
        {
            printf ("ERROR: Observation file %s format error!\n", fn);
            PIHMExit (EXIT_FAILURE);
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
                    printf
                        ("Fatal Error: No %s output available for member %d at %d (%d)!",
                        ens->var[k].name, i + 1, obs_time, (int)buffer[0]);
                    PIHMExit (EXIT_FAILURE);
                }
            }
        }
    }

    free (buffer);
}

void WriteEnKFOut (char *project, enkf_struct ens, char *outputdir, int t)
{
    int             i, j, k;
    int             id;
    time_t          rawtime;
    struct tm      *timestamp;
    FILE           *fid;
    FILE           *fid1;
    char            varn[MAXSTRING];
    char            fn[MAXSTRING];
    double         *x;
    int             ne;
    double          outtime;

    rawtime = (time_t) t;
    timestamp = gmtime (&rawtime);

    ne = ens->ne;

    x = (double *)malloc (ne * sizeof (double));

    /*
     * Write parameter output
     */
    for (i = 0; i < MAXPARAM; i++)
    {
        if (ens->param[i].update == 1 && ens->update_param == 1)
        {
            for (j = 0; j < ne; j++)
            {
                x[j] = ens->member[j].param[i];
            }

            sprintf (fn, "%s%s.txt", outputdir, ens->param[i].name);
            fid = fopen (fn, "a");
            fprintf (fid, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
                timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);

            for (j = 0; j < ne; j++)
            {
                fprintf (fid, "\t%lf", x[j]);
            }
            fprintf (fid, "\n");
            fflush (fid);
            fclose (fid);
        }
    }

    for (i = 0; i < ne; i++)
    {
        /*
         * Write variable/flux output
         */
        for (k = 0; k < MAXVAR; k++)
        {
            if (ens->var[k].dim > 0)
            {
                sprintf (fn, "%s%s.%3.3d.%s.dat",
                    outputdir, project, i + 1, ens->var[k].name);
                fid = fopen (fn, "ab");
                outtime = (double)t;
                fwrite (&outtime, sizeof (double), 1, fid);
                for (j = 0; j < ens->var[k].dim; j++)
                {
                    fwrite (&ens->member[i].var[k][j], sizeof (double), 1,
                        fid);
                }
                fflush (fid);
                fclose (fid);

                if (ens->ascii)
                {
                    sprintf (fn, "%s%s.%3.3d.%s.txt",
                        outputdir, project, i + 1, ens->var[k].name);
                    fid = fopen (fn, "a");
                    fprintf (fid, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
                        timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                        timestamp->tm_mday, timestamp->tm_hour,
                        timestamp->tm_min);
                    for (j = 0; j < ens->var[k].dim; j++)
                    {
                        fprintf (fid, "\t%lf", ens->member[i].var[k][j]);
                    }
                    fprintf (fid, "\n");
                    fflush (fid);
                    fclose (fid);
                }
            }
        }

        /*
         * Write PIHM initial condition
         */
        sprintf (fn, "input/%s/%s.%3.3d.ic", project, project, i + 1);
        fid = fopen (fn, "wb");
        CheckFile (fid, fn);

        sprintf (fn, "input/%s/%s.%3.3d.ic.txt", project, project, i + 1);
        fid1 = fopen (fn, "w");
        CheckFile (fid1, fn);

        for (j = 0; j < ens->numele; j++)
        {
            id = FindVar (ens->var, "is");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "snow");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "surf");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "unsat");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "gw");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "t1");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "snowh");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            for (k = 0; k < MAXLYR; k++)
            {
                sprintf (varn, "stc%d", k);
                id = FindVar (ens->var, varn);
                fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
                fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);
            }
            for (k = 0; k < MAXLYR; k++)
            {
                sprintf (varn, "smc%d", k);
                id = FindVar (ens->var, varn);
                fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
                fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);
            }
            for (k = 0; k < MAXLYR; k++)
            {
                sprintf (varn, "swc%d", k);
                id = FindVar (ens->var, varn);
                fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
                fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);
            }

            fprintf (fid1, "\n");
        }

        for (j = 0; j < ens->numriv; j++)
        {
            id = FindVar (ens->var, "stage");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "rivgw");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\n", ens->member[i].var[id][j]);
        }
        fflush (fid);
        fclose (fid);
        fflush (fid1);
        fclose (fid1);
    }

    free (x);
}
