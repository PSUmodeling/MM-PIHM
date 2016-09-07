#include "pihm.h"

void Perturb (enkf_struct ens, char *outputdir)
{
    int             ne;
    int             i, j;
    int             n = 0;
    char            fn[MAXSTRING];

    calib_struct    cal;

    double        **x;
    double        **randnum;
    int             ind[MAXPARAM];

    double          prior;
    double          prior_std;

    ne = ens->ne;

    x = (double **)malloc (ne * sizeof (double));

    /*
     * Generate initial parameter values 
     */
    for (i = 0; i < MAXPARAM; i++)
    {
        if (ens->param[i].type == LOG_TYPE)
        {
            ens->param[i].min = log10 (ens->param[i].min);
            ens->param[i].max = log10 (ens->param[i].max);
            ens->param[i].perturb_min = log10 (ens->param[i].perturb_min);
            ens->param[i].perturb_max = log10 (ens->param[i].perturb_max);
        }
    }

    if (ens->start_mode == 3)
    {
        /* Using existing .calib files to initialize model parameters */
        for (i = 0; i < ne; i++)
        {
            sprintf (fn, "input/%s/%s.%3.3d.calib", project, project, i + 1);
            ReadCalib (fn, &cal);
            Calib2Mbr (cal, ens->member[i].param);
        }

        /* Calculate initial standard deviation */
        for (i = 0; i < MAXPARAM; i++)
        {
            if (ens->param[i].perturb == 1)
            {
                prior = 0.0;
                prior_std = 0.0;

                for (j = 0; j < ne; j++)
                {
                    if (ens->param[i].type == 1)
                    {
                        prior += log10 (ens->member[j].param[i]);
                    }
                    else
                    {
                        prior += ens->member[j].param[i];
                    }
                }

                prior /= (double)ne;

                for (j = 0; j < ne; j++)
                {
                    if (ens->param[i].type == LOG_TYPE)
                    {
                        prior_std +=
                            (log10 (ens->member[j].param[i]) -
                            prior) * (log10 (ens->member[j].param[i]) -
                            prior);
                    }
                    else
                    {
                        prior_std +=
                            (ens->member[j].param[i] -
                            prior) * (ens->member[j].param[i] - prior);
                    }
                }
                prior_std = sqrt (prior_std / ((double)ne - 1.0));

                ens->param[i].init_std = prior_std;
            }
        }
    }
    else if (ens->start_mode == 2)
    {
        /* Perturb around .calib values */
//
//      ReadCalib (project, project, &cal);
//
//        srand(time(NULL));
//
//        printf("\n\nInitial parameters");
//
//        for (i=0; i<NUMPRMT; i++)
//        {
//            if (En->prmt[i].perturb == 1)
//            {
//                printf("\n perturbed prmt %s", En->prmt[i].prmtn);
//                prmtfn = (char *)malloc((strlen(En->prmt[i].prmtn)+12)*sizeof(char));
//                strcpy(prmtfn,"output/");
//                strcat(prmtfn, En->prmt[i].prmtn);
//                prmtfile = fopen(strcat(prmtfn,".dat"),"w");
//                free(prmtfn);
//                fprintf(prmtfile,"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_hour,timeinfo->tm_min);
//                printf("\n%s\n", En->prmt[i].prmtn);
//                for (j=0; j<ne; j++)
//                {
//                    x[j] = &En->member[j].parameter[i];
//
//                }
//
//                if (En->prmt[i].type == 1)
//                {
//                    En->prmt[i].min = log10(En->prmt[i].min);
//                    En->prmt[i].max = log10(En->prmt[i].max);
//                    En->prmt[i].perturb_min = log10(En->prmt[i].perturb_min);
//                    En->prmt[i].perturb_max = log10(En->prmt[i].perturb_max);
//                }
//
//                prmt_mean = 0;
//
//                for (j=0; j<ne; j++)
//                {
//                    *x[j] = (j+1)*(En->prmt[i].max - En->prmt[i].min)/(ne+1)+En->prmt[i].min;
//                    prmt_mean = prmt_mean + *x[j];
//                }
//                prmt_mean = prmt_mean/(double)ne;
//
//
//                /* Calculate initial standard deviation */
//                En->prmt[i].init_std = 0;
//                for (j=0; j<ne; j++)
//                {
//                    En->prmt[i].init_std = (*x[j]-prmt_mean)*(*x[j]-prmt_mean)+En->prmt[i].init_std;
//                    if (En->prmt[i].type == 1)
//                    {
//                        *x[j] = pow(10, *x[j]);
//                    }
//                    fprintf(prmtfile,"\t%lf",*x[j]);
//                    printf("%lf\t",*x[j]);
//                }
//                En->prmt[i].init_std = sqrt(En->prmt[i].init_std/((double)ne-1));
//
//                if (En->prmt[i].type == 1)
//                {
//                    prmt_mean = pow(10, prmt_mean);
//                }
//                fprintf(prmtfile,"\t%lf\n",prmt_mean);
//                printf("%lf",prmt_mean);
//                printf("\nInitial std %f", En->prmt[i].init_std);
//                fflush(prmtfile);
//                fclose(prmtfile);
//            }
//        }
//        for (i=0; i<ne; i++)
//        {
//            caliboutfn = (char *)malloc((strlen(projectname)+17)*sizeof(char));
//            sprintf(tmpLName,"%3.3d",i+1);
//            strcpy(caliboutfn,"input/");
//            strcat(caliboutfn,projectname);
//            strcat(caliboutfn,".");
//            strcat(caliboutfn,tmpLName);
//            caliboutfile = fopen(strcat(caliboutfn, ".calib"), "w");
//            free(caliboutfn);
//
//            fprintf(caliboutfile,"%lf\t%lf\t%lf\t%lf\t%lf",En->member[i].parameter[0],En->member[i].parameter[1],En->member[i].parameter[2],En->member[i].parameter[3],En->member[i].parameter[4]);
//            fprintf(caliboutfile,"\n%lf\t%lf\t%lf",En->member[i].parameter[5],En->member[i].parameter[6],En->member[i].parameter[7]);
//            fprintf(caliboutfile,"\n%lf\t%lf\t%lf",En->member[i].parameter[8],En->member[i].parameter[9],En->member[i].parameter[10]);
//            fprintf(caliboutfile,"\n%lf\t%lf",En->member[i].parameter[11],En->member[i].parameter[12]);
//            fprintf(caliboutfile,"\n%lf\t%lf\t%lf",En->member[i].parameter[13],En->member[i].parameter[14],En->member[i].parameter[15]);
//            fprintf(caliboutfile,"\n%lf\t%lf",En->member[i].parameter[16],En->member[i].parameter[17]);
//            fprintf(caliboutfile,"\n%lf\t%lf\t%lf",En->member[i].parameter[18],En->member[i].parameter[19],  En->member[i].parameter[20]);
//            fprintf(caliboutfile,"\n%lf\t%lf\t%lf\t%lf",En->member[i].parameter[21],En->member[i].parameter[22],En->member[i].parameter[23],En->member[i].parameter[24]);
//            fprintf(caliboutfile,"\n%lf\t%lf",En->member[i].parameter[25],En->member[i].parameter[26]);
//            fprintf(caliboutfile,"\n%lf\t%lf",En->member[i].parameter[27],En->member[i].parameter[28]);
//            fprintf(caliboutfile,"\n%lf\t%lf\t%lf\t%lf",En->member[i].parameter[29],En->member[i].parameter[30],En->member[i].parameter[31],En->member[i].parameter[32]);
//            fprintf(caliboutfile,"\n%lf\t%lf\t%lf\t%lf\t%lf",En->member[i].parameter[33],En->member[i].parameter[34],En->member[i].parameter[35],En->member[i].parameter[36],En->member[i].parameter[37]);
//            fclose(caliboutfile);
//        }
    }
    else if (ens->start_mode < 2)
    {
        /* Starting from initial mean determined by start_mode */
        sprintf (fn, "input/%s/%s.calib", project, project);
        ReadCalib (fn, &cal);

        for (i = 0; i < ne; i++)
        {
            Calib2Mbr (cal, ens->member[i].param);
        }

        n = 0;

        for (i = 0; i < MAXPARAM; i++)
        {
            if (ens->param[i].perturb == 1)
            {
                ind[n] = i;
                n++;
            }
        }

        randnum = (double **)malloc (n * sizeof (double *));
        for (i = 0; i < n; i++)
        {
            randnum[i] = (double *)malloc (ne * sizeof (double));
        }

        /* The initial standard deviation of perturbed parameter is 1/5 of
         * the parameter range. The allowed range for generated random number
         * is therefore [-2.5 - mode, 2.5 - mode] */

        GenRandNum (ne, n, randnum,
            -2.5 + (double)ens->start_mode, 2.5 - (double)ens->start_mode);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < ne; j++)
            {
                x[j] = &ens->member[j].param[ind[i]];
            }

            prior_std =
                0.2 * (ens->param[ind[i]].max - ens->param[ind[i]].min);
            prior =
                (ens->param[ind[i]].min + ens->param[ind[i]].max) / 2.0 +
                (double)ens->start_mode * prior_std;

            for (j = 0; j < ne; j++)
            {
                *x[j] = prior + prior_std * randnum[i][j];

                if (ens->param[ind[i]].type == LOG_TYPE)
                {
                    *x[j] = pow (10, *x[j]);
                }
            }

            if (ens->param[ind[i]].type == LOG_TYPE)
            {
                prior = pow (10.0, prior);
            }

            ens->param[ind[i]].init_std = prior_std;

        }

        WriteCalFile (ens, project);

        for (i = 0; i < n; i++)
        {
            free (randnum[i]);
        }
        free (randnum);
    }

    PIHMprintf (VL_NORMAL, "\nInitial parameters\n");

    for (i = 0; i < MAXPARAM; i++)
    {
        if (ens->param[i].perturb == 1)
        {
            PIHMprintf (VL_NORMAL, "%s:\n", ens->param[i].name);

            prior = 0.0;

            for (j = 0; j < ne; j++)
            {
                PIHMprintf (VL_NORMAL, "%lf\t", ens->member[j].param[i]);

                prior += ens->member[j].param[i];
            }

            prior /= (double)ne;

            PIHMprintf (VL_NORMAL, "mean: %lf\n", prior);

            PIHMprintf (VL_NORMAL, "Initial std %lf\n", ens->param[i].init_std);

            WriteParamOutput (ens->cycle_start_time, ens, i, outputdir);
        }
    }

    free (x);
}
