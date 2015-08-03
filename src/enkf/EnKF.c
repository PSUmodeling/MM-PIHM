#include "pihm.h"
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

    printf("\n\nStarting EnKF ... \n");

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

    /* Define observation errors and create "observations" from truth */

    if (ens->nobs > 0)
    {
        sprintf (obsfn, "%s/obs.dat", outputdir);
        obsfile = fopen (obsfn, "a");
        fprintf (obsfile, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
            timestamp->tm_year + 1900, timestamp->tm_mon + 1,
            timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);

        for (i = 0; i < ens->nobs; i++)
        {
            sprintf (obsin_fn, "input/%s/%s", project, ens->obs[i].fn);
            ReadObs (obs_time, obsin_fn, &obs, &obs_error);

            if (ens->obs[i].type == RUNOFF_OBS)
            {
                obs = log (obs + 1.0);
            }
            printf("%s observation = %lf\n", ens->obs[i].fn, obs);
            printf("%s error = %lf\n", ens->obs[i].fn, obs_error);

            /* Read ensemble forecasts */

            ReadFcst (ens, ens->obs[i], xf);

            /* Prepare forecast vectors */

            printf("\n*****Predictions******\n");

            printf("Predictions of %s\n", ens->obs[i].name);
            for (j = 0; j < ne; j++)
            {
                printf("%f\t", xf[j]);
            }
            printf("%f\n", xf[ne]);

            /* Write observations to files */

            fprintf (obsfile, "\t%lf", obs);

    //        update_analysis(DS, ens, *obs, *obs_error, xf);
        }
        fprintf(obsfile, "\n");
        fflush (obsfile);
        fclose (obsfile);
    //    inflate(DS, ens, En0);

    }
    //write_EnKF_output(filename, DS, ens, timestamp);
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

//void update_analysis(ensemble ens, double obs, double obs_error, double *xf)
//{
//    int i, j, k;
//    int ne;
//
//    double *xa;
//
//    double **x;
//
//    ne = ens->ne;
//
//    xa = (double *)malloc(sizeof(double)*(ne+1));
//
//    x = (double **)malloc(ne*sizeof(double));
//
//    /* Optimize Parameters from observations */
//
//    for (i=0; i<NUMPRMT; i++)
//    {
//        if (ens->prmt[i].update == 1 & ens->update_prmt == 1)
//        {
//            /* Make pointers point to the parameter that needs to be updated */
//            for (j=0; j<ne; j++)
//            {
//                x[j] = &ens->member[j].parameter[i];
//            }
//
//            /* Take log if the parameter is conductivity or Czil */
//
//            if (ens->prmt[i].type == 1)
//            {
//                for (j=0; j<ne; j++)
//                {
//                    *x[j] = log10(*x[j]);
//                }
//            }
//
//            /* Convert parameters to -inf:inf space and calculate initial pertubations */
//
//            xa[ne] = 0;
//
//            for (j = 0; j<ne; j++)
//            {
//                xa[j] = *x[j];
//                xa[ne] = xa[ne]+xa[j];
//            }
//
//            xa[ne] = xa[ne]/(double)ne;
//
//            EnKF_core(xa, obs, obs_error, xf, ne);
//
//            for (j = 0; j<ne; j++)
//            {
//                *x[j] = xa[j];
//                if (ens->prmt[i].type == 1)
//                {
//                    *x[j] = pow(10,*x[j]);
//                }
//            }
//        }
//
//    }
//
//    if (ens->update_stvrbl == 1)
//    {
//        for (i=0; i<NUMVRBL; i++)
//        {
//            for (k=0; k<ens->vrbl[i].vrbl_dim; k++)
//            {
//                for (j=0; j<ne; j++)
//                {
//                    x[j] = &ens->member[j].variable[i][k];
//                }
//
//                xa[ne] = 0;
//
//                for (j = 0; j<ne; j++)
//                {
//                    xa[j] = *x[j];
//                    xa[ne] = xa[ne]+xa[j];
//                }
//                xa[ne] = xa[ne]/(double)ne;
//
//                EnKF_core(xa, obs, obs_error, xf, ne);
//                for (j = 0; j<ne; j++)
//                {
//                    *x[j] = xa[j];
//                }
//            }
//        }
//    }
//
//
//    free(xa);
//    free(x);
//}
//
//void inflate(Model_Data DS, ensemble ens, ensemble En0)
//{
//    char *prmtn[NUMPRMT];		/* Modified by Y. Shi */
//    char *prmtfn;
//    FILE *prmtfile;
//
//    int i, j, k;
//    int ne;
//
//    double *xp0, *xp;
//
//    double vrbl_min, vrbl_max;
//    double vrbl_mean, vrbl_mean0;
//
//    double **x, *x0;
//    double prmt_mean, prmt_mean0, prmt_min, prmt_max;
//    double prmt_std;
//    double vrbl_temp;
//    double c1, c2, c;
//
//    double TotalWaterVolume[ens->ne];
//    double MassCoeff[ens->ne];
//
//
//    ne = ens->ne;
//
//    xp0 = (double *)malloc(sizeof(double)*ne);
//    xp = (double *)malloc(sizeof(double)*ne);
//
//    x = (double **)malloc(ne*sizeof(double));
//    x0 = (double *)malloc(ne*sizeof(double));
//
//
//    /* Print out parameter names */
//
//    printf("\n\n*****Parameters********");
//
//    for (i=0; i<NUMPRMT; i++)
//    {
//        if (ens->prmt[i].update == 1 & ens->update_prmt == 1)
//        {
//            printf("\n %s\n", ens->prmt[i].prmtn);
//            /* Make pointers point to the parameter that needs to be updated */
//            for (j=0; j<ne; j++)
//            {
//                x[j] = &ens->member[j].parameter[i];
//                x0[j] = En0->member[j].parameter[i];
//            }
//
//
//            /* Take log if the parameter is conductivity or Czil */
//
//            if (ens->prmt[i].type == 1)
//            {
//                for (j=0; j<ne; j++)
//                {
//                    *x[j] = log10(*x[j]);
//                    x0[j] = log10(x0[j]);
//                }
//            }
//
//            prmt_mean0 = 0;
//            prmt_mean = 0;
//
//            for (j = 0; j<ne; j++)
//            {
//                prmt_mean0 = prmt_mean0+x0[j];
//                prmt_mean = prmt_mean+*x[j];
//            }
//
//            prmt_mean = prmt_mean/(double)ne;
//            prmt_mean0 = prmt_mean0/(double)ne;
//
//            for (j=0; j<ne; j++)
//            {
//                xp0[j] = x0[j] - prmt_mean0;
//                xp[j] = *x[j] - prmt_mean;
//            }
//
//            if (prmt_mean<ens->prmt[i].max - 0.25*ens->prmt[i].init_std && prmt_mean>ens->prmt[i].min + 0.25*ens->prmt[i].init_std)
//            {
//                // Calculate new perturbations and standard deviation
//
//                prmt_std = 0;
//
//                for (j = 0; j<ne; j++)
//                {
//                    xp[j] = (1-ens->weight)*xp[j]+ens->weight*xp0[j];
//                    prmt_std = prmt_std + xp[j]*xp[j];
//                }
//                prmt_std = sqrt(prmt_std/((double)ne-1));
//
//                // Coavariance inflation
//
//                prmt_min = 999;
//                prmt_max = -999;
//
//                for (j = 0; j<ne; j++)
//                {
//                    if (prmt_std<0.25*ens->prmt[i].init_std)
//                    {
//                        xp[j] = 0.25*ens->prmt[i].init_std/prmt_std*xp[j];
//                    }
//                    prmt_min = prmt_min<(prmt_mean+xp[j])?prmt_min:(prmt_mean+xp[j]);
//                    prmt_max = prmt_max>(prmt_mean+xp[j])?prmt_max:(prmt_mean+xp[j]);
//                }
//                c1 = (ens->prmt[i].max-EPSILON-prmt_mean)/(prmt_max-prmt_mean);
//                c2 = (prmt_mean-ens->prmt[i].min-EPSILON)/(prmt_mean-prmt_min);
//                c = c1<c2?c1:c2;
//                c = c<1?c:1.0;
//
//                for (j = 0; j<ne; j++)
//                {
//                    *x[j] = prmt_mean + c*xp[j];
//                    if (ens->prmt[i].type == 1)
//                    {
//                        *x[j] = pow(10,*x[j]);
//                    }
//                    printf("%lf\t",*x[j]);
//                }
//                if (ens->prmt[i].type == 1)
//                {
//                    prmt_mean = pow(10,prmt_mean);
//                }
//                printf("%lf",prmt_mean);
//            }
//            else
//            {
//                printf("\tParameter not updated");
//                prmt_mean = prmt_mean0;
//                for (j = 0; j<ne; j++)
//                {
//                    *x[j] = x0[j];
//                    if (ens->prmt[i].type == 1)
//                    {
//                        *x[j] = pow(10, x0[j]);
//                    }
//                }
//
//                if (ens->prmt[i].type == 1)
//                {
//                    prmt_mean = pow(10,prmt_mean);
//                }
//            }
//        }
//    }
//
//    if (ens->update_stvrbl == 1)
//    {
//        for (i=0; i<NUMVRBL; i++)
//        {
//            for (k=0; k<ens->vrbl[i].vrbl_dim; k++)
//            {
//                for (j=0; j<ne; j++)
//                {
//                    x[j] = &ens->member[j].variable[i][k];
//                    x0[j] = En0->member[j].variable[i][k];
//                }
//                vrbl_mean = 0;
//                vrbl_mean0 = 0;
//
//                for (j = 0; j<ne; j++)
//                {
//                    vrbl_mean = vrbl_mean+*x[j];
//                    vrbl_mean0 = vrbl_mean0+x0[j];
//                }
//                vrbl_mean = vrbl_mean/(double)ne;
//                vrbl_mean0 = vrbl_mean0/(double)ne;
//
//                for (j = 0; j<ne; j++)
//                {
//                    xp0[j] = x0[j]-vrbl_mean0;
//                    xp[j] = *x[j]-vrbl_mean;
//                    *x[j] = vrbl_mean + (1-ens->weight)*xp[j]+ens->weight*xp0[j];
//                    /*
//                       if (i < 2 | i > 11)
//                       {
//                     *x[j] = *x[j]<vrbl_min?vrbl_min:*x[j];
//                     }
//                     else
//                     {
//                     *x[j] = *x[j]>vrbl_min?(*x[j]<vrbl_max?*x[j]:x0[j]):x0[j];
//                     }
//                     */
//                }
//            }
//        }
//
//        printf("\nMass conservation coefficient\n");
//        for (j=0; j<ne; j++)
//        {
//            TotalWaterVolume[j] = En0->TotalWaterVolume[j];
//            for (k = 0; k<DS->NumEle; k++)
//            {
//                ens->TotalWaterVolume[j] = ens->TotalWaterVolume[j] + DS->Ele[k].area*(ens->member[j].variable[0][k] + ens->member[j].variable[2][k]*DS->Ele[k].Porosity + ens->member[j].variable[3][k] + ens->member[j].variable[8][k] + ens->member[j].variable[7][k] + ens->member[j].variable[22][k]/Lv/1000*3600);
//            }
//            for (k = 0; k<DS->NumRiv; k++)
//            {
//                ens->TotalWaterVolume[j] = ens->TotalWaterVolume[j] + DS->Riv[k].coeff*DS->Riv[k].Length*(ens->member[j].variable[1][k] + ens->member[j].variable[8][k+DS->NumEle] + ens->member[j].variable[3][k+DS->NumEle]);
//            }
//            ens->TotalWaterVolume[j] = ens->TotalWaterVolume[j] + ens->member[j].variable[10][DS->NumRiv-1]/24;
//            MassCoeff[j] = ens->TotalWaterVolume[j]/TotalWaterVolume[j];
//            printf("%f\t", MassCoeff[j]);
//        }
//        /*
//           for (j=0; j<ne; j++)
//           {
//           for (k = 0; k<DS->NumEle; k++)
//           {
//           ens->IC[j].Overland[k] = MassCoeff[j]*ens->IC[j].Overland[k];
//           ens->IC[j].Unsat[k] = MassCoeff[j]*ens->IC[j].Unsat[k];
//           ens->IC[j].IS[k] = MassCoeff[j]*ens->IC[j].IS[k];
//           ens->IC[j].LE[k] = MassCoeff[j]*ens->IC[j].LE[k];
//           }
//           for (k = 0; k<DS->NumEle+DS->NumRiv; k++)
//           {
//           ens->IC[j].Sat[k] = MassCoeff[j]*ens->IC[j].Sat[k];
//           ens->IC[j].Snow[k] = MassCoeff[j]*ens->IC[j].Snow[k];
//           }
//           ens->IC[j].rivFlx[1] = MassCoeff[j]*ens->IC[j].rivFlx[1];
//           }
//           */
//    }
//    free(xp0);
//    free(xp);
//    free(x);
//    free(x0);
//
//    printf("\nInflate done!");
//
//}
//
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
            exit(1);
        }
        else if (match != 7)
        {
            printf ("ERROR: Observation file %s format error!\n", fn);
            exit(1);
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
    double          x;

    ne = ens->ne;

    for (i = 0; i < ne + 1; i++)
    {
        xf[i] = 0.0;
    }

    for (i = 0; i < ne; i++)
    {
        for (j = 0; j < obs.nctrl; j++)
        {
            var_ind = obs.var_ind[j];

            for (k = 0; k < ens->var[var_ind].dim; k++)
            {
                x += ens->member[i].var[var_ind][k] *
                    obs.k[j] + obs.b[j];
                xf[i] += x * obs.weight[k];
            }
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
                sprintf (fn, "%s/%s.%3.3d.%s.dat",
                    outputdir, project, i + 1, ens->var[k].name);
                fid = fopen (fn, "rb");
                CheckFile (fid, fn);

                fseek (fid, 0L, SEEK_END);

                length = (int) (ftell (fid) / (ens->var[k].dim + 1) / 8);

                rewind (fid);

                for (ii = 0; ii < length; ii++)
                {
                    fread (buffer, sizeof (double), ens->var[k].dim + 1, fid);

                    if ((int) buffer[0] == obs_time)
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
                    printf("Fatal Error: No %s output available!", ens->var[k].name);
                    exit(1);
                }
            }
        }
    }

    free (buffer);
}
//void write_EnKF_output(char *projectname, Model_Data DS, ensemble ens, struct tm *timestamp)
//{
//    char indchar[3], *prmtfn, *ofn, *initfn, *calibfn;
//    FILE *Ofile, *initfile, *prmtfile, *calibfile;
//    double *x, prmtmean;
//    int i, j, k;
//    int ne;
//
//    ne = ens->ne;
//
//    x = (double *)malloc(ne*sizeof(double));
//
//    for (i=0; i<NUMPRMT; i++)
//    {
//        if (ens->prmt[i].update == 1 & ens->update_prmt == 1)
//        {
//            for (j=0; j<ne; j++)
//            {
//                x[j] = ens->member[j].parameter[i];
//            }
//
//            prmtfn = (char *)malloc((strlen(ens->prmt[i].prmtn)+12)*sizeof(char));
//            strcpy(prmtfn,"output/");
//            strcat(prmtfn, ens->prmt[i].prmtn);
//            prmtfile = fopen(strcat(prmtfn,".dat"),"a");
//            free(prmtfn);
//            fprintf(prmtfile,"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",timestamp->tm_year+1900,timestamp->tm_mon+1,timestamp->tm_mday,timestamp->tm_hour,timestamp->tm_min);
//
//            prmtmean = 0;
//            for (j=0; j<ne; j++)
//            {
//                fprintf(prmtfile,"\t%lf",x[j]);
//                if (ens->prmt[i].type == 1)
//                {
//                    prmtmean = prmtmean + log10(x[j]);
//                }
//                else
//                {
//                    prmtmean = prmtmean + x[j];
//                }
//            }
//            prmtmean = prmtmean/(double)(ens->ne);
//            if (ens->prmt[i].type == 1)
//            {
//                prmtmean = pow(10, prmtmean);
//            }
//            fprintf(prmtfile, "\t%lf\n",prmtmean);
//            fflush(prmtfile);
//            fclose(prmtfile);
//        }
//    }
//
//    for (i=0; i<ne; i++)
//    {
//        sprintf(indchar, "%3.3d", i+1);
//
//        for (k=0; k<NUMVRBL; k++)
//        {
//            ofn = (char *)malloc((strlen(projectname)+strlen(ens->vrbl[k].vrbln)+13)*sizeof(char));
//            strcpy(ofn,"output/");
//            strcat(ofn, projectname);
//            strcat(ofn,".");
//            strcat(ofn, indchar);
//            strcat(ofn,".");
//            Ofile=fopen(strcat(ofn, ens->vrbl[k].vrbln),"a");
//            free(ofn);
//            fprintf(Ofile,"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",timestamp->tm_year+1900,timestamp->tm_mon+1,timestamp->tm_mday,timestamp->tm_hour,timestamp->tm_min);
//            for (j=0; j<ens->vrbl[k].vrbl_dim; j++)
//            {
//                fprintf(Ofile, "\t%lf", ens->member[i].variable[k][j]);
//            }
//            fprintf(Ofile, "\n");
//            fflush(Ofile);
//            fclose(Ofile);
//        }
//
//
//        initfn = (char *)malloc((strlen(projectname)+16)*sizeof(char));
//        strcpy(initfn, "input/");
//        strcat(initfn, projectname);
//        strcat(initfn, ".");
//        strcat(initfn, indchar);
//        initfile = fopen(strcat(initfn, ".init"),"w");
//        free(initfn);
//
//        for(j=0; j<DS->NumEle; j++)
//        {
//            fprintf(initfile, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf", ens->member[i].variable[7][j],ens->member[i].variable[8][j],ens->member[i].variable[0][j],ens->member[i].variable[2][j],ens->member[i].variable[3][j],ens->member[i].variable[23][j]);
//            for (k=0; k<DS->NumSoilLayer; k++)
//            {
//                fprintf(initfile, "\t%lf", ens->member[i].variable[24+k][j]);
//            }
//            for (k=0; k<DS->NumSoilLayer+1;k++)
//            {
//                fprintf(initfile, "\t%lf", ens->member[i].variable[24+DS->NumSoilLayer+k][j]);
//            }
//            fprintf(initfile, "\n");
//        }
//        for(j=0; j<DS->NumRiv; j++)
//        {
//            fprintf(initfile, "%lf\t%lf\t%lf\n",ens->member[i].variable[1][j],ens->member[i].variable[3][j+DS->NumEle],ens->member[i].variable[8][j+DS->NumEle]);
//        } 
//        fflush(initfile);
//        fclose(initfile);
//
//        calibfn = (char *)malloc((strlen(projectname)+17)*sizeof(char));
//        strcpy(calibfn,"input/");
//        strcat(calibfn,projectname);
//        strcat(calibfn,".");
//        strcat(calibfn,indchar);
//        calibfile = fopen(strcat(calibfn, ".calib"), "w");
//        free(calibfn);
//
//        fprintf(calibfile,"%lf\t%lf\t%lf\t%lf\t%lf",ens->member[i].parameter[0],ens->member[i].parameter[1],ens->member[i].parameter[2],ens->member[i].parameter[3],ens->member[i].parameter[4]);
//        fprintf(calibfile,"\n%lf\t%lf\t%lf",ens->member[i].parameter[5],ens->member[i].parameter[6],ens->member[i].parameter[7]);
//        fprintf(calibfile,"\n%lf\t%lf\t%lf",ens->member[i].parameter[8],ens->member[i].parameter[9],ens->member[i].parameter[10]);
//        fprintf(calibfile,"\n%lf\t%lf",ens->member[i].parameter[11],ens->member[i].parameter[12]);
//        fprintf(calibfile,"\n%lf\t%lf\t%lf",ens->member[i].parameter[13],ens->member[i].parameter[14],ens->member[i].parameter[15]);
//        fprintf(calibfile,"\n%lf\t%lf",ens->member[i].parameter[16],ens->member[i].parameter[17]);
//        fprintf(calibfile,"\n%lf\t%lf\t%lf",ens->member[i].parameter[18],ens->member[i].parameter[19],  ens->member[i].parameter[20]);
//        fprintf(calibfile,"\n%lf\t%lf\t%lf\t%lf",ens->member[i].parameter[21],ens->member[i].parameter[22],ens->member[i].parameter[23],ens->member[i].parameter[24]);
//        fprintf(calibfile,"\n%lf\t%lf",ens->member[i].parameter[25],ens->member[i].parameter[26]);
//        fprintf(calibfile,"\n%lf\t%lf",ens->member[i].parameter[27],ens->member[i].parameter[28]);
//        fprintf(calibfile,"\n%lf\t%lf\t%lf\t%lf",ens->member[i].parameter[29],ens->member[i].parameter[30],ens->member[i].parameter[31],ens->member[i].parameter[32]);
//        fprintf(calibfile,"\n%lf\t%lf\t%lf\t%lf\t%lf",ens->member[i].parameter[33],ens->member[i].parameter[34],ens->member[i].parameter[35],ens->member[i].parameter[36],ens->member[i].parameter[37]);
//        fflush(calibfile);
//        fclose(calibfile);
//
//    }
//    free(x);
//}
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
