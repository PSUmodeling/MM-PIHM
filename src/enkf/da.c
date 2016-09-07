#include "pihm.h"

void EnKFDA (enkf_struct ens, int obs_time, char *outputdir)
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

    PIHMprintf (VL_NORMAL, "\nStarting EnKF ... \n");

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
            PIHMprintf (VL_NORMAL, "\n*****%s******\n", ens->obs[i].name);

            /* 
             * Read observations
             */
            sprintf (obsin_fn, "input/%s/%s", project, ens->obs[i].fn);
            ReadObs (obs_time, obsin_fn, &obs, &obs_error);

            PIHMprintf (VL_NORMAL, "observation = %lf\n", obs);
            PIHMprintf (VL_NORMAL, "error = %lf\n", obs_error);

            /*
             * Read ensemble forecasts
             */
            ReadFcst (ens, ens->obs[i], xf);

            /* 
             * Prepare forecast vectors
             */
            PIHMprintf (VL_NORMAL, "prediction = ");
            for (j = 0; j < ne; j++)
            {
                PIHMprintf (VL_NORMAL, "%f\t", xf[j]);
            }
            PIHMprintf (VL_NORMAL, "mean: %f\n", xf[ne]);

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

            EnKF (xa, obs, obs_error, xf, ne);

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

                    EnKF (xa, obs, obs_error, xf, ne);

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
