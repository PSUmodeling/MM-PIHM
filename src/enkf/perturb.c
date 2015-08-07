#include "pihm.h" 
#include "noah.h"
#include "enkf.h"

void Perturb(char *project, enkf_struct ens, char *outputdir)
{
    int             ne;
    int             i, j;
    int             n = 0;
    time_t          rawtime;
    char            simulation[MAXSTRING];

    calib_struct    cal;

    double        **x;
    double        **randnum;
    int             ind[MAXPARAM];

    double          prior;
    double          prior_std;

    /* Initialize ensemble members */
    ne = ens->ne;

    ens->member = (ens_mbr_struct *) malloc (ne * sizeof (ens_mbr_struct));

    /*
     * Define variable controls: vairable names, variable dimension, etc.
     */
    MapVar (ens->var, ens->numele, ens->numriv);

//    En->TotalWaterVolume = (double *)malloc(ne*sizeof(double));
//
    x = (double **) malloc (ne * sizeof (double));
    
    rawtime = (int) ens->cycle_start_time;

    printf("Ensemble members: %d\n", ne);
    printf("Default observation cycle: %-d hour(s)\n", ens->interval /  3600);
    printf("Observations:");
    if (ens->nobs == 0)
    {
        printf(" none");
    }
    else
    {
        for (i = 0; i < ens->nobs - 1; i++)
        {
            printf(" %s,", ens->obs[i].name);
        }
        printf(" %s\n", ens->obs[ens->nobs - 1].name);
    }

    for (i = 0; i < ne; i++)
    {
        for (j = 0; j < MAXVAR; j++)
        {
            if (ens->var[j].dim > 0)
            {
                ens->member[i].var[j] =
                    (double *) malloc (ens->var[j].dim * sizeof(double));
            }
        }

//        if (En->memStartMode == 3)
//        {
//            initinfn = (char *)malloc((strlen(projectname)+17)*sizeof(char));
//            sprintf(tmpLName,"%3.3d",i+1);
//            strcpy(initinfn, "input/");
//            strcat(initinfn,projectname);
//            strcat(initinfn,".");
//            strcat(initinfn,tmpLName);
//            initinfile = fopen(strcat(initinfn, ".init0"), "r");
//            free(initinfn);
//
//            if(initinfile == NULL)
//            {
//                printf("\n  Fatal Error: %s.init0 is in use or does not exist!\n", projectname);
//                fclose(initinfile);
//                PihmExit(1);
//            }
//            else
//            {
//                initoutfn = (char *)malloc((strlen(projectname)+16)*sizeof(char));
//                sprintf(tmpLName,"%3.3d",i+1);
//                strcpy(initoutfn, "input/");
//                strcat(initoutfn,projectname);
//                strcat(initoutfn,".");
//                strcat(initoutfn,tmpLName);
//                initoutfile = fopen(strcat(initoutfn, ".init"), "w");
//                free(initoutfn);
//                for(j=0; j<DS->NumEle; j++)
//                {
//                    /* Read from init0 file */
//                    fscanf(initinfile, "%lf %lf %lf %lf %lf %lf", &En->member[i].variable[7][j],&En->member[i].variable[8][j],&En->member[i].variable[0][j],&En->member[i].variable[2][j],&En->member[i].variable[3][j],&En->member[i].variable[23][j]);
//                    for (k=0; k<DS->NumSoilLayer; k++)
//                    {
//                        fscanf(initinfile, "%lf", &En->member[i].variable[24+k][j]);
//                    }
//                    for (k=0; k<DS->NumSoilLayer+1;k++)
//                    {
//                        fscanf(initinfile, "%lf", &En->member[i].variable[24+DS->NumSoilLayer+k][j]);
//                    }
//                    /* Write to init file */
//                    fprintf(initoutfile, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf", En->member[i].variable[7][j],En->member[i].variable[8][j],En->member[i].variable[0][j],En->member[i].variable[2][j],En->member[i].variable[3][j],En->member[i].variable[23][j]);
//                    //					fprintf(initoutfile, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", En->IC[i].IS[j],En->IC[i].Snow[j],En->IC[i].Overland[j],En->IC[i].Unsat[j],En->IC[i].Sat[j],En->IC[i].Tsfc[j],En->IC[i].Tsoil[j][0],En->IC[i].Tsoil[j][1],En->IC[i].Tsoil[j][2],En->IC[i].Tsoil[j][3],En->IC[i].SM[j][0],En->IC[i].SM[j][1],En->IC[i].SM[j][2],En->IC[i].SM[j][3],En->IC[i].SM[j][4]);
//                    for (k=0; k<DS->NumSoilLayer; k++)
//                    {
//                        fprintf(initoutfile, "\t%lf", En->member[i].variable[24+k][j]);
//                    }
//                    for (k=0; k<DS->NumSoilLayer+1;k++)
//                    {
//                        fprintf(initoutfile, "\t%lf", En->member[i].variable[24+DS->NumSoilLayer+k][j]);
//                    }
//                    fprintf(initoutfile, "\n");
//                }
//                for(j=0; j<DS->NumRiv; j++)
//                {
//                    fscanf(initinfile, "%lf %lf %lf",&En->member[i].variable[1][j],&En->member[i].variable[3][j+DS->NumEle],&En->member[i].variable[8][j+DS->NumEle]);
//                    fprintf(initoutfile, "%lf\t%lf\t%lf\n",En->member[i].variable[1][j],En->member[i].variable[3][j+DS->NumEle],En->member[i].variable[8][j+DS->NumEle]);
//                    //					fprintf(initoutfile, "%lf\t%lf\t%lf\n",En->IC[i].RiverState[j],En->IC[i].Sat[j+DS->NumEle],En->IC[i].Snow[j+DS->NumEle]);
//                } 
//            }
//            fflush(initoutfile);
//            fclose(initinfile); 
//            fclose(initoutfile);
//
//        }
    }
//
//
    /*
     * Generate initial parameter values 
     */

    for (i = 0; i < MAXPARAM; i++)
    {
        if (ens->param[i].type == 1)
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
            sprintf (simulation, "%s.%3.3d", project, i + 1);
            ReadCalib (project, simulation, &cal);
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

                prior /= (double) ne;

                for (j = 0; j < ne; j++)
                {
                    if (ens->param[i].type == 1)
                    {
                        prior_std += (log10 (ens->member[j].param[i]) - prior) * (log10 (ens->member[j].param[i]) - prior);
                    }
                    else
                    {
                        prior_std += (ens->member[j].param[i] - prior) * (ens->member[j].param[i] - prior);
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
        ReadCalib (project, project, &cal);

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
            -2.5 + (double) ens->start_mode, 2.5 - (double) ens->start_mode);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < ne; j++)
            {
                x[j] = &ens->member[j].param[ind[i]];
            }

            prior_std = 0.2 * (ens->param[ind[i]].max - ens->param[ind[i]].min);
            prior = (ens->param[ind[i]].min + ens->param[ind[i]].max) / 2.0 + (double)ens->start_mode * prior_std;

            for (j = 0; j < ne; j++)
            {
                *x[j] = prior + prior_std * randnum[i][j];

                if (ens->param[ind[i]].type == 1)
                {
                    *x[j] = pow (10, *x[j]);
                }
            }

            if (ens->param[ind[i]].type == 1)
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
    
    printf("\nInitial parameters\n");

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

            printf("Initial std %lf\n", ens->param[i].init_std);

            WriteParamOutput (ens->cycle_start_time, ens, i, outputdir);
        }
    }

    free(x);
}

void MapVar (var_struct *var, int numele, int numriv)
{
    int             i, j;

    for (i = 0; i < MAXVAR; i++)
    {
        var[i].dim = 0;
    }

    for (i = 0; i < MAXVAR; i++)
    {
        switch (i)
        {
            case 0:
                strcpy(var[i].name, "surf");
                var[i].dim = numele;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 1:
                strcpy(var[i].name, "stage");
                var[i].dim = numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 2:
                strcpy(var[i].name, "unsat");
                var[i].dim = numele;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 3:
                strcpy(var[i].name, "gw");
                var[i].dim = numele + numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 4:
                strcpy(var[i].name, "et0");
                var[i].dim = numele;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 5:
                strcpy(var[i].name, "et1");
                var[i].dim = numele;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 6:
                strcpy(var[i].name, "et2");
                var[i].dim = numele;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 7:
                strcpy(var[i].name, "is");
                var[i].dim = numele;
                var[i].min = 0.0;
                var[i].max = BADVAL;
                break;
            case 8:
                strcpy(var[i].name,  "snow");
                var[i].dim = numele;
                var[i].min = 0.0;
                var[i].max = BADVAL;
                break;
            case 9:
                strcpy(var[i].name, "rivflx0");
                var[i].dim = numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 10:
                strcpy(var[i].name,  "rivflx1");
                var[i].dim = numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 11:
                strcpy(var[i].name,  "rivflx2");
                var[i].dim = numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 12:
                strcpy(var[i].name, "rivflx3");
                var[i].dim = numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 13:
                strcpy(var[i].name,  "rivflx4");
                var[i].dim = numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 14:
                strcpy(var[i].name, "rivflx5");
                var[i].dim = numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 15:
                strcpy(var[i].name,  "rivflx6");
                var[i].dim = numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 16:
                strcpy(var[i].name,  "rivflx7");
                var[i].dim = numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 17:
                strcpy(var[i].name,  "rivflx8");
                var[i].dim = numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 18:
                strcpy(var[i].name, "rivflx9");
                var[i].dim = numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 19:
                strcpy(var[i].name, "rivflx10");
                var[i].dim = numriv;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 20:
                strcpy(var[i].name,  "recharge");
                var[i].dim = numele;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 21:
                strcpy(var[i].name,  "g");
                var[i].dim = numele;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 22:
                strcpy(var[i].name,  "sh");
                var[i].dim = numele;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 23:
                strcpy(var[i].name,  "le");
                var[i].dim = numele;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 24:
                strcpy(var[i].name, "t1");
                var[i].dim = numele;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
                break;
            case 25:
                strcpy(var[i].name, "snowh");
                var[i].dim = numele;
                var[i].min = 0.0;
                var[i].max = BADVAL;
                break;
            default:
                break;
        }

        for (j = 0; j < MAXLYR; j++)
        {
            if (i == 26 + j)
            {
                sprintf(var[i].name, "stc%-d", j);
                var[i].dim = numele;
                var[i].min = BADVAL;
                var[i].max = BADVAL;
            }

            if (i == 26 + MAXLYR + j)
            {
                sprintf(var[i].name, "smc%-d", j);
                var[i].dim = numele;
                var[i].min = 0.0;
                var[i].max = BADVAL;
            }

            if (i == 26 + 2 * MAXLYR + j)
            {
                sprintf (var[i].name, "swc%-d", j);
                var[i].dim = numele;
                var[i].min = 0.0;
                var[i].max = BADVAL;
            }
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
    param[ET0] = cal.et[0];
    param[ET1] = cal.et[1];
    param[ET2] = cal.et[2];
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
    cal->et[0] = param[ET0];
    cal->et[1] = param[ET1];
    cal->et[2] = param[ET2];
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
    struct tm  *timeinfo;
    FILE           *fid;
    int             i;

    timevar = (time_t) rawtime;
    timeinfo = gmtime(&timevar);

    sprintf (fn, "%s/%s.txt", outputdir, ens->param[ind].name);
    fid = fopen (fn, "w");

    fprintf (fid, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
        timeinfo->tm_year+1900, timeinfo->tm_mon + 1, timeinfo->tm_mday,
        timeinfo->tm_hour, timeinfo->tm_min);
    for (i = 0; i < ens->ne; i++)
    {
        fprintf (fid, "\t%lf", ens->member[i].param[ind]);
    }
    fprintf (fid, "\n");
    fflush(fid);
    fclose(fid); 
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
                fprintf (fid, "%-16s%-6lf\n", ens->param[j].name, ens->member[i].param[j]);
            }
        }

        fflush (fid);
        fclose (fid);
    }

}

double Randn()
{
    double temp1, temp2;
    double x;

    temp1 = (double) (rand() % MAXINT + 1) / MAXINT;
    temp2 = (double) (rand() % MAXINT + 1) / MAXINT;

    x = sqrt (-2.0 * log (temp1)) * cos (2.0 * PI * temp2);

    return (x);
}

void GenRandNum (int ne, int nparam, double **randnum, double lower, double upper)
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

    srand(time(NULL));

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
                    randnum[j][i] = Randn();
                    mean[j] += randnum[j][i];
                }

                mean[j] /= (double) ne;

                for (i = 0; i < ne; i++)
                {
                    std += (randnum[j][i] - mean[j]) * (randnum[j][i] - mean[j]);
                }
                std = sqrt (std / ((double) ne - 1.0));

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
                    corr[i][j] = (randnum[i][k] - mean[i]) * (randnum[j][k] - mean[j]) + corr[i][j];
                    s1 = s1 + (randnum[i][k] - mean[i]) * (randnum[i][k] - mean[i]);
                    s2 = s2 + (randnum[j][k] - mean[j]) * (randnum[j][k] - mean[j]);
                }

                corr[i][j] = corr[i][j] / sqrt (s1 * s2);

                if (fabs(corr[i][j]) > max && i != j)
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
