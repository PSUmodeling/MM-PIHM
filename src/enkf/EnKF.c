#include "pihm.h"
#include "noah.h"
#include "enkf.h"

//#define NUMVRBL 25+2*DS->NumSoilLayer
//#define NUMPRMT 38
//
//int datenum(int year, int month, int day, int hour, int min, int sec)
//{
//    struct tm *timeinfo;
//    time_t rawtime;
//
//    timeinfo = (struct tm *)malloc(sizeof(struct tm));
//
//    timeinfo->tm_year = year-1900;
//    timeinfo->tm_mon = month - 1;
//    timeinfo->tm_mday = day;
//    timeinfo->tm_hour = hour;
//    timeinfo->tm_min = min;
//    timeinfo->tm_sec = sec;
//
//    rawtime = timegm(timeinfo);
//
//    free (timeinfo):
//
//        return rawtime;
//}

void EnKFCore(double *xa, double obs, double obs_err, double *xf, int ne)
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

    hxa = (double *) malloc (ne * sizeof (double));
    xp = (double *) malloc (ne * sizeof (double));

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

void EnKF (char *project, enkf_struct ens, int obs_time, char *outputdir)
{
    double         *xf;
    int             i, j;
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
    ens0->member = (ens_mbr_struct *)malloc (ne * sizeof(ens_mbr_struct));

    for (i = 0; i < ne; i++)
    {
        for (j = 0; j < MAXVAR; j++)
        {
            if (ens->var[j].dim > 0)
            {
                ens0->member[i].var[j] =
                    (double *)malloc(ens->var[j].dim * sizeof(double));
            }
        }
    }

    xf = (double *)malloc (sizeof (double) * (ne + 1));
    rawtime = obs_time;
    timestamp = gmtime(&rawtime);

    //srand(time(NULL));

    printf("\nStarting EnKF ... \n");

    /* Copy prior from ens to En0 */

    for (i = 0; i < ne; i++)
    {
        ens0->member[i] = ens->member[i];
    }
    //    for (j=0; j<NUMPRMT; j++)
    //    {
    //        En0->member[i].parameter[j] = ens->member[i].parameter[j];
    //    }
    //    for (j=0; j<NUMVRBL; j++)
    //    {
    //        for(k=0; k<ens->vrbl[j].vrbl_dim; k++)
    //        {
    //            En0->member[i].variable[j][k] = ens->member[i].variable[j][k];
    //        }
    //    }
    //}


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
            printf("\n*****%s******\n", ens->obs[i].name);

            /* Read observations */
            sprintf (obsin_fn, "input/%s/%s", project, ens->obs[i].fn);
            ReadObs (obs_time, obsin_fn, &obs, &obs_error);

            //if (ens->obs[i].type == RUNOFF_OBS)
            //{
            //    obs = log (obs + 1.0);
            //}
            printf("observation = %lf\n", obs);
            printf("error = %lf\n", obs_error);

            /* Read ensemble forecasts */

            ReadFcst (ens, ens->obs[i], xf);

            /* Prepare forecast vectors */

            printf("prediction = ");
            for (j = 0; j < ne; j++)
            {
                printf("%f\t", xf[j]);
            }
            printf("mean: %f\n", xf[ne]);

            /* Write observations to files */

            fprintf (obsfile, "\t%lf", obs);

            UpdAnlys (ens, obs, obs_error, xf);
        }

        fprintf(obsfile, "\n");
        fflush (obsfile);
        fclose (obsfile);

        /* Covariance inflation */
        CovInflt(ens, ens0);
    }

    WriteEnKFOut (project, ens, outputdir, obs_time);
    //printf("\n\nEnKF done.");
    //free(obs);
    //free(obs_error); 
    //free(En0->TotalWaterVolume);
    //for (i=0; i<ens->ne; i++)
    //{
    //    for (j=0; j<NUMVRBL; j++)
    //    {
    //        free(En0->member[i].variable[j]);
    //    }
    //    free(En0->member[i].variable);
    //    free(En0->member[i].parameter);
    //}
    //free(En0->member);
    //free(En0);
    free(xf);
}

void UpdAnlys (enkf_struct ens, double obs, double obs_error, double *xf)
{
    int             i, j, k;
    int             ne;

    double         *xa;

    double        **x;

    ne = ens->ne;

    xa = (double *)malloc (sizeof (double) * (ne + 1));

    x = (double **)malloc (ne * sizeof(double));

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

            if (ens->param[i].type == 1)
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
                if (ens->param[i].type == 1)
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

    x = (double **)malloc (ne * sizeof(double));
    x0 = (double *)malloc (ne * sizeof(double));

    printf("\n*****Parameters********\n");

    for (i = 0; i < MAXPARAM; i++)
    {
        if (ens->param[i].update == 1 && ens->update_param == 1)
        {
            printf("%s\n", ens->param[i].name);
            /* Make pointers point to the parameter that needs to be updated */
            for (j = 0; j < ne; j++)
            {
                x[j] = &ens->member[j].param[i];
                x0[j] = ens0->member[j].param[i];
            }

            /* Take log if the parameter is conductivity or Czil */

            if (ens->param[i].type == 1)
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
                    xp[j] = (1.0 - ens->weight) * xp[j] + ens->weight * xp0[j];
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
                        xp[j] = 0.25 * ens->param[i].init_std / param_std * xp[j];
                    }
                    param_min = (param_min < average + xp[j]) ? param_min: (average + xp[j]);
                    param_max = (param_max > average + xp[j]) ? param_max: (average + xp[j]);
                }

                c1 = (ens->param[i].max - 1.0E-6 - average) / (param_max - average);
                c2 = (average - ens->param[i].min - 1.0E-6) / (average - param_min);
                c = (c1<c2) ? c1 : c2;
                c = (c < 1.0) ? c : 1.0;

                for (j = 0; j < ne; j++)
                {
                    *x[j] = average + c * xp[j];
                    if (ens->param[i].type == 1)
                    {
                        *x[j] = pow (10.0, *x[j]);
                    }
                    printf("%lf\t",*x[j]);
                }
                if (ens->param[i].type == 1)
                {
                    average = pow (10.0, average);
                }
                printf("mean: %lf\n", average);
            }
            else
            {
                printf("EnKF analysis %lf is out of range. Parameter is not updated\n", average);
                average = average0;
                for (j = 0; j < ne; j++)
                {
                    *x[j] = x0[j];
                    if (ens->param[i].type == 1)
                    {
                        *x[j] = pow (10.0, x0[j]);
                    }
                }

                if (ens->param[i].type == 1)
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
                        *x[j] = average + (1.0 - ens->weight) * xp[j] + ens->weight * xp0[j];
                        if ((int)ens->var[i].max != BADVAL)
                        {
                            *x[j] = (*x[j] < ens->var[i].max) ? *x[j] : ens->var[i].max;
                        }
                        
                        if ((int)ens->var[i].min != BADVAL)
                        {
                            *x[j] = (*x[j] > ens->var[i].min) ? *x[j] : ens->var[i].min;
                        }
                    }
                }
            }
        }
    }

    free(xp0);
    free(xp);
    free(x);
    free(x0);

    printf("Inflation done!\n");
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
            &timeinfo->tm_hour,&timeinfo->tm_min, &temp1, &temp2);
        timeinfo->tm_year = timeinfo->tm_year - 1900;
        timeinfo->tm_mon = timeinfo->tm_mon - 1;
        timeinfo->tm_sec = 0;
        rawtime = timegm(timeinfo);

        if (rawtime == obs_time)
        {
            *obs = temp1;
            *obs_error = temp2;
            break;
        }
        else if (strcasecmp (cmdstr, "EOF") == 0)
        {
            printf("\nFATAL ERROR: No observation availablein %s!\n",
                fn);
            PihmExit(1);
        }
        else if (match != 7)
        {
            printf ("ERROR: Observation file %s format error!\n", fn);
            PihmExit(1);
        }
    }

    fclose (fid);
    free(timeinfo);
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

    xf[ne] /= (double) ne;
}

void ReadVar (char *project, char *outputdir, enkf_struct ens, int obs_time)
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

    buffer = (double *) malloc ((ens->numele + ens->numriv + 1) * sizeof (double));

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

                length = (int) (ftell (fid) / (ens->var[k].dim + 1) / 8);

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

                if (success == 0)
                {
                    printf("Fatal Error: No %s output available for member %d at %d (%d)!", ens->var[k].name, i + 1, obs_time, (int)buffer[0]);
                    PihmExit(1);
                }

                fclose (fid);
            }
        }
    }

    free (buffer);
}

void WriteEnKFOut (char *project, enkf_struct ens, char *outputdir, int t)
{
    //char indchar[3], *prmtfn, *ofn, *initfn, *calibfn;
    //FILE *Ofile, *initfile, *prmtfile, *calibfile;
    //double *x, prmtmean;
    int             i, j, k;
    time_t          rawtime;
    struct tm      *timestamp;
    FILE           *fid;
    char            fn[MAXSTRING];
    double         *x;
    int             ne;
    double          outtime;

    rawtime = (time_t)t;
    timestamp = gmtime(&rawtime);

    ne = ens->ne;

    x = (double *)malloc (ne * sizeof(double));

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
                timestamp->tm_year+1900, timestamp->tm_mon + 1,
                timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);

            for (j = 0; j < ne; j++)
            {
                fprintf (fid, "\t%lf", x[j]);
            }
            fprintf (fid, "\n");
            fflush(fid);
            fclose(fid);
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
                outtime = (double) t;
                fwrite (&outtime, sizeof (double), 1, fid);
                for (j = 0; j < ens->var[k].dim; j++)
                {
                    fwrite (&ens->member[i].var[k][j], sizeof (double), 1, fid);
                }
                fflush(fid);
                fclose(fid);

                if (ens->ascii)
                {
                    sprintf (fn, "%s%s.%3.3d.%s.txt",
                        outputdir, project, i + 1, ens->var[k].name);
                    fid = fopen (fn, "a");
                    fprintf (fid, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
                        timestamp->tm_year+1900, timestamp->tm_mon + 1,
                        timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
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
        sprintf (fn, "input/%s/%s.%3.3d.init", project, project, i + 1);
        fid = fopen (fn, "wb");
        CheckFile (fid, fn);

        for (j = 0; j < ens->numele; j++)
        {
            fwrite (&ens->member[i].var[7][j], sizeof (double), 1, fid);
            fwrite (&ens->member[i].var[8][j], sizeof (double), 1, fid);
            fwrite (&ens->member[i].var[0][j], sizeof (double), 1, fid);
            fwrite (&ens->member[i].var[2][j], sizeof (double), 1, fid);
            fwrite (&ens->member[i].var[3][j], sizeof (double), 1, fid);
        }

        for (j = 0; j < ens->numriv; j++)
        {
            fwrite (&ens->member[i].var[1][j], sizeof (double), 1, fid);
            fwrite (&ens->member[i].var[3][j + ens->numele], sizeof (double), 1, fid);
        } 
        fflush(fid);
        fclose(fid);

        /*
         * Write Noah initial condition
         */
        sprintf (fn, "input/%s/%s.%3.3d.lsminit", project, project, i + 1);
        fid = fopen (fn, "wb");
        CheckFile (fid, fn);

        for (j = 0; j < ens->numele; j++)
        {
            fwrite (&ens->member[i].var[24][j], sizeof (double), 1, fid);
            fwrite (&ens->member[i].var[25][j], sizeof (double), 1, fid);
            for (k = 0; k < MAXLYR; k++)
            {
                fwrite (&ens->member[i].var[26 + k][j], sizeof (double), 1, fid);
            }
            for (k = 0; k < MAXLYR; k++)
            {
                fwrite (&ens->member[i].var[26 + MAXLYR + k][j], sizeof (double), 1, fid);
            }
            for (k = 0; k < MAXLYR; k++)
            {
                fwrite (&ens->member[i].var[26 + 2 * MAXLYR + k][j], sizeof (double), 1, fid);
            }
        }
        fflush (fid);
        fclose (fid);
    }

    free (x);
}
//void freeEnsemble(Model_Data DS, ensemble ens)
//{
//    int i, j;
//
//    free(ens->TotalWaterVolume);
//    for (i=0; i<ens->ne; i++)
//    {
//        for (j=0; j<NUMVRBL; j++)
//        {
//            free(ens->member[i].variable[j]);
//        }
//        free(ens->member[i].variable);
//        free(ens->member[i].parameter);
//    }
//    for (i=0; i<ens->no_obs; i++)
//    {
//        free(ens->observation[i].vrbl_ind);
//        free(ens->observation[i].grid_ind);
//        free(ens->observation[i].weight);
//    }
//    free(ens->observation);
//    free(ens->member);
//    free(ens->vrbl);
//}
