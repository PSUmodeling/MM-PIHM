#include "pihm.h"

void JobHandout (int starttime, int endtime, int startmode,
    ensmbr_struct *member, double *param, int ne, int total_jobs)
{
    int             dest;
    int             msg[3];
    int             i, j;

    msg[0] = starttime;
    msg[1] = endtime;
    msg[2] = startmode;

    for (i = 0; i < ne; i++)
    {
        for (j = 0; j < MAXPARAM; j++)
        {
            param[i * MAXPARAM + j] = member[i].param[j];
        }
    }

    for (dest = 1; dest < total_jobs + 1; dest++)
    {
        MPI_Send (msg, 3, MPI_INT, dest, CYCLE_TAG, MPI_COMM_WORLD);
        MPI_Send (param, ne * MAXPARAM, MPI_DOUBLE, dest, PARAM_TAG,
            MPI_COMM_WORLD);
    }
}

void JobRecv (int *starttime, int *endtime, int *startmode, double *param,
    int ne)
{
    int             msg[3];
    MPI_Status      status;

    MPI_Recv (msg, 3, MPI_INT, 0, CYCLE_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv (param, ne * MAXPARAM, MPI_DOUBLE, 0, PARAM_TAG, MPI_COMM_WORLD,
        &status);

    *starttime = msg[0];
    *endtime = msg[1];
    *startmode = msg[2];
}

void PrintEnKFStatus (int starttime, int endtime)
{
    time_t          rawtime;
    struct tm      *timestamp;

    printf ("\nRunning ensemble members from ");

    rawtime = (int)starttime;
    timestamp = gmtime (&rawtime);
    printf ("%4.4d-%2.2d-%2.2d %2.2d:%2.2d to ",
        timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday,
        timestamp->tm_hour, timestamp->tm_min);

    rawtime = (int)endtime;
    timestamp = gmtime (&rawtime);
    printf ("%4.4d-%2.2d-%2.2d %2.2d:%2.2d\n",
        timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday,
        timestamp->tm_hour, timestamp->tm_min);
}

void JobHandIn (int total_jobs)
{
    int             ierr;
    int             success;
    int             received = 0;
    int             source;

    MPI_Status      status;

    while (received < total_jobs)
    {
        ierr =
            MPI_Recv (&success, 1, MPI_INT, MPI_ANY_SOURCE, SUCCESS_TAG,
            MPI_COMM_WORLD, &status);
        received++;
        source = status.MPI_SOURCE;
        //printf("PIHM job handed in from Node: %d\n", source);
    }
}

void InitOper (pihm_struct pihm, enkf_struct ens)
{
    int             i;

    for (i = 0; i < ens->nobs; i++)
    {
        if (strcasecmp (ens->obs[i].name, "discharge") == 0)
        {
            ens->obs[i].type = RUNOFF_OBS;
            DisOper (&ens->obs[i], ens->var, pihm);
        }
        else if (strcasecmp (ens->obs[i].name, "skin_temperature") == 0)
        {
            ens->obs[i].type = TSKIN_OBS;
            LandSfcTmpOper (&ens->obs[i], ens->var, pihm);
        }
        else
        {
            printf ("ERROR: Cannot find the operator for %s!\n",
                ens->obs[i].name);
            PihmExit (1);
        }
    }
}

void FreeEns (enkf_struct ens)
{
    int             i, j;

    for (i = 0; i < ens->nobs; i++)
    {
        free (ens->obs[i].var_ind);
        free (ens->obs[i].weight);
        for (j = 0; j < ens->obs[i].dim; j++)
        {
            free (ens->obs[i].k[j]);
            free (ens->obs[i].b[j]);
        }
        free (ens->obs[i].k);
        free (ens->obs[i].b);
    }

    free (ens->obs);

    for (i = 0; i < ens->ne; i++)
    {
        for (j = 0; j < MAXVAR; j++)
        {
            if (ens->var[j].dim > 0)
            {
                free (ens->member[i].var[j]);
            }
        }
    }

    free (ens->member);
}

void InitEns (enkf_struct ens)
{
    int             ne;
    pihm_struct     pihm;
    N_Vector        CV_Y;       /* State Variables Vector */
    int             nsv;
    int             i, j;
    char            outputdir[MAXSTRING];

    outputdir[0] = '\0';

    pihm = (pihm_struct)malloc (sizeof *pihm);

    ReadAlloc (project, pihm);

    /* problem size */
    nsv = 3 * pihm->numele + 2 * pihm->numriv;

    CV_Y = N_VNew_Serial (nsv);

    Initialize (pihm, CV_Y);

    MapOutput (project, pihm, outputdir);

    ens->numele = pihm->numele;
    ens->numriv = pihm->numriv;
    ens->ascii = pihm->ctrl.ascii;

    /* Initialize ensemble members */
    ne = ens->ne;

    ens->member = (ensmbr_struct *)malloc (ne * sizeof (ensmbr_struct));

    /*
     * Define variable controls: vairable names, variable dimension, etc.
     */
    MapVar (ens->var, ens->numele, ens->numriv);

    InitOper (pihm, ens);

    printf ("Ensemble members: %d\n", ne);
    printf ("Default observation cycle: %-d hour(s)\n", ens->interval / 3600);
    printf ("Observations:");
    if (ens->nobs == 0)
    {
        printf (" none");
    }
    else
    {
        for (i = 0; i < ens->nobs - 1; i++)
        {
            printf (" %s,", ens->obs[i].name);
        }
        printf (" %s\n", ens->obs[ens->nobs - 1].name);
    }

    for (i = 0; i < ne; i++)
    {
        for (j = 0; j < MAXVAR; j++)
        {
            if (ens->var[j].dim > 0)
            {
                ens->member[i].var[j] =
                    (double *)malloc (ens->var[j].dim * sizeof (double));
            }
        }

    }

    N_VDestroy_Serial (CV_Y);

    FreeData (pihm);
    free (pihm);
}

void Perturb (char *project, enkf_struct ens, char *outputdir)
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

    printf ("\nInitial parameters\n");

    for (i = 0; i < MAXPARAM; i++)
    {
        if (ens->param[i].perturb == 1)
        {
            printf ("%s:\n", ens->param[i].name);

            prior = 0.0;

            for (j = 0; j < ne; j++)
            {
                printf ("%lf\t", ens->member[j].param[i]);

                prior += ens->member[j].param[i];
            }

            prior /= (double)ne;

            printf ("mean: %lf\n", prior);

            printf ("Initial std %lf\n", ens->param[i].init_std);

            WriteParamOutput (ens->cycle_start_time, ens, i, outputdir);
        }
    }

    free (x);
}

void MapVar (var_struct *var, int numele, int numriv)
{
    int             i, k;
    int             n = 0;

    for (i = 0; i < MAXVAR; i++)
    {
        var[i].dim = 0;
    }

    for (i = 0; i < MAXVAR; i++)
    {
        switch (i)
        {
            case SURF_CTRL:
                strcpy (var[n].name, "surf");
                var[n].dim = numele;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case UNSAT_CTRL:
                strcpy (var[n].name, "unsat");
                var[n].dim = numele;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case GW_CTRL:
                strcpy (var[n].name, "gw");
                var[n].dim = numele;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case RIVSTG_CTRL:
                strcpy (var[n].name, "stage");
                var[n].dim = numriv;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case RIVGW_CTRL:
                strcpy (var[n].name, "rivgw");
                var[n].dim = numriv;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case SNOW_CTRL:
                strcpy (var[n].name, "snow");
                var[n].dim = numele;
                var[n].min = 0.0;
                var[n].max = BADVAL;
                n++;
                break;
            case CMC_CTRL:
                strcpy (var[n].name, "is");
                var[n].dim = numele;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case INFIL_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RECHARGE_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case EC_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case ETT_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case EDIR_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX0_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX1_CTRL:
                strcpy (var[n].name, "rivflx1");
                var[n].dim = numriv;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case RIVFLX2_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX3_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX4_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX5_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX6_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX7_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX8_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX9_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case RIVFLX10_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case SUBFLX_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case SURFFLX_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case T1_CTRL:
                strcpy (var[n].name, "t1");
                var[n].dim = numele;
                var[n].min = BADVAL;
                var[n].max = BADVAL;
                n++;
                break;
            case STC_CTRL:
                for (k = 0; k < MAXLYR; k++)
                {
                    sprintf (var[n].name, "stc%d", k);
                    var[n].dim = numele;
                    var[n].min = BADVAL;
                    var[n].max = BADVAL;
                    n++;
                }
                break;
            case SMC_CTRL:
                for (k = 0; k < MAXLYR; k++)
                {
                    sprintf (var[n].name, "smc%d", k);
                    var[n].dim = numele;
                    var[n].min = BADVAL;
                    var[n].max = BADVAL;
                    n++;
                }
                break;
            case SH2O_CTRL:
                for (k = 0; k < MAXLYR; k++)
                {
                    sprintf (var[n].name, "swc%d", k);
                    var[n].dim = numele;
                    var[n].min = BADVAL;
                    var[n].max = BADVAL;
                    n++;
                }
                break;
            case SNOWH_CTRL:
                strcpy (var[n].name, "snowh");
                var[n].dim = numele;
                var[n].min = 0.0;
                var[n].max = BADVAL;
                n++;
                break;
            case ALBEDO_CTRL:
                strcpy (var[n].name, "albedo");
                var[n].dim = numele;
                var[n].min = 0.0;
                var[n].max = 1.0;
                n++;
                break;
            case LE_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case SH_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case G_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case ETP_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case ESNOW_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case ROOTW_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case SOILM_CTRL:
                var[n].dim = 0;
                n++;
                break;
            case SOLAR_CTRL:
                var[n].dim = 0;
                n++;
                break;
            default:
                break;
        }

    }
}

void Calib2Mbr (calib_struct cal, double *param)
{
    param[KSATH] = cal.ksath;
    param[KSATV] = cal.ksatv;
    param[KINF] = cal.kinfv;
    param[KMACH] = cal.kmach;
    param[KMACV] = cal.kmacv;
    param[DINF] = cal.dinf;
    param[RZD] = cal.rzd;
    param[DMAC] = cal.dmac;
    param[POROSITY] = cal.porosity;
    param[ALPHA] = cal.alpha;
    param[BETA] = cal.beta;
    param[AREAFV] = cal.areafv;
    param[AREAFH] = cal.areafh;
    param[VEGFRAC] = cal.vegfrac;
    param[ALBEDO] = cal.albedo;
    param[ROUGH] = cal.rough;
    param[PRCP] = cal.prcp;
    param[SFCTMP] = cal.sfctmp;
    param[EC] = cal.ec;
    param[ETT] = cal.ett;
    param[EDIR] = cal.edir;
    param[RIVROUGH] = cal.rivrough;
    param[RIVKSATH] = cal.rivksath;
    param[RIVKSATV] = cal.rivksatv;
    param[RIVBEDTHICK] = cal.rivbedthick;
    param[RIVDEPTH] = cal.rivdepth;
    param[RIVSHPCOEFF] = cal.rivshpcoeff;
#ifdef _NOAH_
    param[THETAREF] = cal.thetaref;
    param[THETAW] = cal.thetaw;
    param[RSMIN] = cal.rsmin;
    param[DRIP] = cal.drip;
    param[INTCP] = cal.intcp;
    param[CZIL] = cal.czil;
    param[FXEXP] = cal.fxexp;
    param[CFACTR] = cal.cfactr;
    param[RGL] = cal.rgl;
    param[HS] = cal.hs;
#endif
}

void Mbr2Cal (calib_struct *cal, const double *param)
{
    cal->ksath = param[KSATH];
    cal->ksatv = param[KSATV];
    cal->kinfv = param[KINF];
    cal->kmach = param[KMACH];
    cal->kmacv = param[KMACV];
    cal->dinf = param[DINF];
    cal->rzd = param[RZD];
    cal->dmac = param[DMAC];
    cal->porosity = param[POROSITY];
    cal->alpha = param[ALPHA];
    cal->beta = param[BETA];
    cal->areafv = param[AREAFV];
    cal->areafh = param[AREAFH];
    cal->vegfrac = param[VEGFRAC];
    cal->albedo = param[ALBEDO];
    cal->rough = param[ROUGH];
    cal->prcp = param[PRCP];
    cal->sfctmp = param[SFCTMP];
    cal->ec = param[EC];
    cal->ett = param[ETT];
    cal->edir = param[EDIR];
    cal->rivrough = param[RIVROUGH];
    cal->rivksath = param[RIVKSATH];
    cal->rivksatv = param[RIVKSATV];
    cal->rivbedthick = param[RIVBEDTHICK];
    cal->rivdepth = param[RIVDEPTH];
    cal->rivshpcoeff = param[RIVSHPCOEFF];
#ifdef _NOAH_
    cal->thetaref = param[THETAREF];
    cal->thetaw = param[THETAW];
    cal->rsmin = param[RSMIN];
    cal->drip = param[DRIP];
    cal->intcp = param[INTCP];
    cal->czil = param[CZIL];
    cal->fxexp = param[FXEXP];
    cal->cfactr = param[CFACTR];
    cal->rgl = param[RGL];
    cal->hs = param[HS];
#endif
}

void WriteParamOutput (int rawtime, enkf_struct ens, int ind, char *outputdir)
{
    char            fn[MAXSTRING];
    time_t          timevar;
    struct tm      *timeinfo;
    FILE           *fid;
    int             i;

    timevar = (time_t) rawtime;
    timeinfo = gmtime (&timevar);

    sprintf (fn, "%s/%s.txt", outputdir, ens->param[ind].name);
    fid = fopen (fn, "w");

    fprintf (fid, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
        timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday,
        timeinfo->tm_hour, timeinfo->tm_min);
    for (i = 0; i < ens->ne; i++)
    {
        fprintf (fid, "\t%lf", ens->member[i].param[ind]);
    }
    fprintf (fid, "\n");
    fflush (fid);
    fclose (fid);
}

double Randn ()
{
    double          temp1, temp2;
    double          x;

    temp1 = (double)(rand () % MAXINT + 1) / MAXINT;
    temp2 = (double)(rand () % MAXINT + 1) / MAXINT;

    x = sqrt (-2.0 * log (temp1)) * cos (2.0 * PI * temp2);

    return (x);
}

void GenRandNum (int ne, int nparam, double **randnum, double lower,
    double upper)
{
    int             i, j, k;
    double          corr[MAXPARAM][MAXPARAM];
    double          mean[MAXPARAM];
    double          std;
    int             corr_flag;
    int             std_flag = 1;
    double          s1, s2;
    double          max = -999.0;

    if (nparam > 1)
    {
        corr_flag = 1;
    }
    else
    {
        corr_flag = 0;
    }

    srand (time (NULL));

    do
    {
        max = -999.0;

        for (j = 0; j < nparam; j++)
        {
            std_flag = 1;

            while (std_flag)
            {
                mean[j] = 0.0;
                std = 0.0;

                for (i = 0; i < ne; i++)
                {
                    randnum[j][i] = Randn ();
                    mean[j] += randnum[j][i];
                }

                mean[j] /= (double)ne;

                for (i = 0; i < ne; i++)
                {
                    std +=
                        (randnum[j][i] - mean[j]) * (randnum[j][i] - mean[j]);
                }
                std = sqrt (std / ((double)ne - 1.0));

                std_flag = 0;

                for (i = 0; i < ne; i++)
                {
                    randnum[j][i] = (randnum[j][i] - mean[j]) / std * 1.0;

                    if (randnum[j][i] <= lower || randnum[j][i] >= upper)
                    {
                        std_flag = 1;
                        break;
                    }
                }
            }
        }

        for (i = 0; i < nparam; i++)
        {
            for (j = 0; j < nparam; j++)
            {
                corr[i][j] = 0.0;
                s1 = 0.0;
                s2 = 0.0;

                for (k = 0; k < ne; k++)
                {
                    corr[i][j] =
                        (randnum[i][k] - mean[i]) * (randnum[j][k] -
                        mean[j]) + corr[i][j];
                    s1 = s1 + (randnum[i][k] - mean[i]) * (randnum[i][k] -
                        mean[i]);
                    s2 = s2 + (randnum[j][k] - mean[j]) * (randnum[j][k] -
                        mean[j]);
                }

                corr[i][j] = corr[i][j] / sqrt (s1 * s2);

                if (fabs (corr[i][j]) > max && i != j)
                {
                    max = corr[i][j];
                }
            }
        }

        if (max < CORRMAX)
        {
            corr_flag = 0;
        }
    } while (corr_flag);
}

void WriteCalFile (enkf_struct ens, char *project)
{
    char            fn[MAXSTRING];
    FILE           *fid;
    int             i, j;

    for (i = 0; i < ens->ne; i++)
    {
        sprintf (fn, "input/%s/%s.%3.3d.calib", project, project, i + 1);
        fid = fopen (fn, "w");

        for (j = 0; j < MAXPARAM; j++)
        {
            if (ens->param[j].name[0] != '\0')
            {
                fprintf (fid, "%-16s%-6lf\n", ens->param[j].name,
                    ens->member[i].param[j]);
            }
        }

        fflush (fid);
        fclose (fid);
    }

}
