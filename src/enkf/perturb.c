#include "pihm.h" 
#include "noah.h"
#include "enkf.h"

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
    double          temp;

    if (nparam > 1)
    {
        corr_flag = 1;
    }
    else
    {
        corr_flag = 0;
    }

    do
    {
        srand(time(NULL));
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
                    printf ("corr = %lf, max = %lf\n", corr[i][j], max);
                }
            }
        }

        if (max < CORRMAX)
        {
            corr_flag = 0;
        }
    } while (corr_flag);

    printf ("lower = %lf, upper = %lf\n", lower, upper);
    for (j = 0; j < ne; j++)
    {
        printf ("%lf, ", randnum[0][j]);
    }
    printf ("\n");
}
    
void Perturb(char *project, enkf_struct ens, char *outputdir)
{
    int ne;
    int i, j, k;
    int n = 0;
    time_t      rawtime;

    int         numele;
    int         numriv;
    riv_att_tbl_struct  riv_att_tbl;
    riv_shp_tbl_struct  riv_shp_tbl;
    riv_matl_tbl_struct riv_matl_tbl;
    riv_ic_tbl_struct   riv_ic_tbl;
    ic_struct           ic;
    forcing_ts_struct   forcing;
    mesh_tbl_struct     mesh_tbl;
    calib_struct        cal;

    double    **x;
    double    **randnum;
    int         ind[MAXPARAM];
    //double prmt_mean;
    //double temp1, temp2, temp3;
    //int NumEle, NumNode, NumRiv;
    //char *fn;
    //char str[3];

    //FILE *prmtfile;
    //FILE *calibfile;
    //FILE *caliboutfile;
    //FILE *calibinfile;
    //FILE *obsfile;

    //char *obsfn;
    //char *calibinfn;
    //char *caliboutfn;
    //FILE *initoutfile;
    //FILE *initinfile;
    //char *initinfn;
    //char *initoutfn;
    //char *prmtfn;
    //char tmpLName[11];

    double          prior;
    double          prior_std;
    //int prmt_flag;

    /* Initialize ensemble members' initial conditions using DS */
    ne = ens->ne;

    ens->member = (ens_mbr_struct *) malloc (ne * sizeof (ens_mbr_struct));

    /*
     * Define variable controls: vairable names, variable dimension, etc.
     */
    ReadRiv (project, &riv_att_tbl, &riv_shp_tbl, &riv_matl_tbl, &riv_ic_tbl,
        &ic, &forcing);
    numriv = riv_att_tbl.number;

    ReadMesh (project, &mesh_tbl);
    numele = mesh_tbl.numele;

    MapVar (ens->var, numele, numriv);

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

//    En->TotalWaterVolume = (double *)malloc(ne*sizeof(double));
//
    x = (double **) malloc (ne * sizeof (double));
    
    rawtime = (int) ens->cycle_start_time;

    printf("Ensemble members: %d\n", ne);
    printf("Default observation cycle: %-d hour(s)\n", ens->interval);
//    printf("\nObservations:");
//    if (En->no_obs == 0)
//    {
//        printf(" none");
//    }
//    else
//    {
//        for (i=0; i<En->no_obs-1; i++)
//        {
//            printf(" %s,", En->observation[i].obsn);
//        }
//        printf(" %s.", En->observation[En->no_obs-1].obsn);
//    }
//    obsfn = "output/obs.dat";
//    obsfile = fopen(obsfn, "w");
//    fclose(obsfile);

    for (i = 0; i < ne; i++)
    {
        ens->member[i].param = (double *) malloc (MAXPARAM * sizeof(double));
        ens->member[i].var = (double **) malloc (MAXVAR * sizeof(double));

        for (j = 0; j < MAXVAR; j++)
        {
            ens->member[i].var[j] =
                (double *) malloc (ens->var[j].dim * sizeof(double));
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
//                exit(1);
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
    if (ens->start_mode == 3)
    {
        /* Using .calib0 files to initialize model parameters */
//        for (i=0; i<ne; i++)
//        {
//            calibinfn = (char *)malloc((strlen(projectname)+18)*sizeof(char));
//            sprintf(tmpLName,"%3.3d",i+1);
//            strcpy(calibinfn,"input/");
//            strcat(calibinfn,projectname);
//            strcat(calibinfn,".");
//            strcat(calibinfn,tmpLName);
//            calibinfile = fopen(strcat(calibinfn, ".calib0"), "r");
//            free(calibinfn);
//
//            if(calibinfile == NULL)
//            {
//                printf("\n  Fatal Error: %s.calib0 is in use or does not exist!\n", projectname);
//                exit(1);
//            }
//
//            /* start reading calib_file */
//            for (j=0; j<NUMPRMT; j++)
//            {
//                fscanf(calibinfile,"%lf", &En->member[i].parameter[j]);
//            }
//            //			 %lf %lf %lf %lf",&En->member[i].KsatH,&En->member[i].KsatV,&En->member[i].infKsatV,&En->member[i].macKsatH,&En->member[i].macKsatV);
//            //			fscanf(calibinfile,"%lf %lf %lf",&En->member[i].infD,&En->member[i].RzD,&En->member[i].macD);
//            //			fscanf(calibinfile,"%lf %lf %lf",&En->member[i].Porosity,&En->member[i].Alpha,&En->member[i].Beta);
//            //			fscanf(calibinfile,"%lf %lf",&En->member[i].vAreaF,&En->member[i].hAreaF);
//            //			fscanf(calibinfile,"%lf %lf %lf",&En->member[i].VegFrac,&En->member[i].Albedo,&En->member[i].Rough);
//            //			fscanf(calibinfile,"%lf %lf",&En->member[i].Prep,&En->member[i].Temp);
//            //			fscanf(calibinfile,"%lf %lf %lf",&temp1,&temp2,&temp3);
//            //			fscanf(calibinfile,"%lf %lf %lf %lf",&En->member[i].rivRough,&En->member[i].rivKsatH,&En->member[i].rivKsatV,&En->member[i].rivbedThick);
//            //			fscanf(calibinfile,"%lf %lf",&En->member[i].rivDepth,&En->member[i].rivShapeCoeff);
//            //			fscanf(calibinfile,"%lf %lf",&En->member[i].TF,&En->member[i].IS);
//            //			fscanf(calibinfile,"%lf %lf %lf %lf",&En->member[i].Rmin,&En->member[i].Czil,&En->member[i].fx_soil,&En->member[i].fx_canopy);
//            //			fscanf(calibinfile,"%lf %lf %lf %lf %lf",&En->member[i].Rs_ref,&En->member[i].h_s,&En->member[i].Tref,&En->member[i].ThetaRef,&En->member[i].ThetaW);
//
//            /* finish reading calib file */  
//            fclose(calibinfile);
//
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
//
//            /* finish reading calib file */  
//            fclose(caliboutfile);
//        }
//
//        for (i=0; i<NUMPRMT; i++)
//        {
//            if (En->prmt[i].perturb == 1)
//            {
//                prmtfn = (char *)malloc((strlen(En->prmt[i].prmtn)+12)*sizeof(char));
//                strcpy(prmtfn,"output/");
//                strcat(prmtfn, En->prmt[i].prmtn);
//                prmtfile = fopen(strcat(prmtfn,".dat"),"w");
//                free(prmtfn);
//                fprintf(prmtfile,"\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_hour,timeinfo->tm_min);
//                printf("\n%s\n", En->prmt[i].prmtn);
//
//                for (j=0; j<ne; j++)
//                {
//                    x[j] = &En->member[j].parameter[i];
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
//                for (j=0; j<ne; j++)
//                {
//                    if (En->prmt[i].type == 1)
//                    {
//                        *x[j] = log10(*x[j]);
//                    }
//                    prmt_mean = prmt_mean + *x[j];
//                }
//                prmt_mean = prmt_mean/(double)ne;
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
//
    }
    else if (ens->start_mode == 2)
    {
        /* Perturb around .calib values */
//
//        for (i=0; i<ne; i++)
//        {
//            fn = (char *)malloc((strlen(projectname)+13)*sizeof(char));
//            strcpy(fn,"input/");
//            strcat(fn,projectname);
//            calibfile = fopen(strcat(fn, ".calib"), "r");
//            free(fn);
//
//            if(calibfile == NULL)
//            {
//                printf("\n  Fatal Error: %s.calib is in use or does not exist!\n", projectname);
//                exit(1);
//            }
//
//            /* start reading calib_file */
//            for (j=0; j<NUMPRMT; j++)
//            {
//                fscanf(calibfile,"%lf", &En->member[i].parameter[j]);
//            }
//            fclose(calibfile);
//        }
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
        /* Initialize ensemble members' calibration coefficient using .calib */
        ReadCalib (project, project, &cal);

        for (i = 0; i < ne; i++)
        {
            Calib2Mbr (cal, ens->member[i].param);
        }

        printf("Initial parameters\n");

        n = 0;

        for (i = 0; i < MAXPARAM; i++)
        {
            if (ens->param[i].perturb == 1)
            {
                ind[n] = i;
                n++;
            }
        }

        randnum = (double **) malloc (n * sizeof (double *));
        for (i = 0; i < n; i++)
        {
            randnum[i] = (double *) malloc (ne * sizeof (double));
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
            prior = (ens->param[ind[i]].min + ens->param[ind[i]].max) / 2.0 + (double) ens->start_mode * prior_std;

            for (j = 0; j < ne; j++)
            {
                *x[j] = prior + prior_std * randnum[i][j];

                if (ens->param[ind[i]].type == 1)
                {
                    *x[j] = pow (10, *x[j]);
                }
                printf("%lf\t", *x[j]);
            }

            ens->param[ind[i]].init_std = prior_std;

            //fprintf(paramfile,"\t%lf",*x[j]);

            printf("\nInitial std %f", ens->param[ind[i]].init_std);
            //fflush(paramfile);
            //fclose(paramfile);
            WriteParamFile (ens->cycle_start_time, ens, ind[i], outputdir);
        }

        WriteCalFile (ens, project);

    }

    for (i = 0; i < n; i++)
    {
        free (randnum[i]);
    }
    free (randnum);

    free(x);
}

void MapVar (var_struct *var, int numele, int numriv)
{
    int             i, j;

    for (i = 0; i < MAXVAR; i++)
    {
        switch (i)
        {
            case 0:
                strcpy(var[i].name, "surf");
                var[i].dim = numele;
                break;
            case 1:
                strcpy(var[i].name, "stage");
                var[i].dim = numriv;
                break;
            case 2:
                strcpy(var[i].name, "unsat");
                var[i].dim = numele;
                break;
            case 3:
                strcpy(var[i].name, "gw");
                var[i].dim = numele+numriv;
                break;
            case 4:
                strcpy(var[i].name, "et0");
                var[i].dim = numele;
                break;
            case 5:
                strcpy(var[i].name, "et1");
                var[i].dim = numele;
                break;
            case 6:
                strcpy(var[i].name, "et2");
                var[i].dim = numele;
                break;
            case 7:
                strcpy(var[i].name, "is");
                var[i].dim = numele;
                break;
            case 8:
                strcpy(var[i].name,  "snow");
                var[i].dim = numele+numriv;
                break;
            case 9:
                strcpy(var[i].name, "rivflx0");
                var[i].dim = numriv;
                break;
            case 10:
                strcpy(var[i].name,  "rivflx1");
                var[i].dim = numriv;
                break;
            case 11:
                strcpy(var[i].name,  "rivflx2");
                var[i].dim = numriv;
                break;
            case 12:
                strcpy(var[i].name, "rivflx3");
                var[i].dim = numriv;
                break;
            case 13:
                strcpy(var[i].name,  "rivflx4");
                var[i].dim = numriv;
                break;
            case 14:
                strcpy(var[i].name, "rivflx5");
                var[i].dim = numriv;
                break;
            case 15:
                strcpy(var[i].name,  "rivflx6");
                var[i].dim = numriv;
                break;
            case 16:
                strcpy(var[i].name,  "rivflx7");
                var[i].dim = numriv;
                break;
            case 17:
                strcpy(var[i].name,  "rivflx8");
                var[i].dim = numriv;
                break;
            case 18:
                strcpy(var[i].name, "rivflx9");
                var[i].dim = numriv;
                break;
            case 19:
                strcpy(var[i].name, "rivflx10");
                var[i].dim = numriv;
                break;
            case 20:
                strcpy(var[i].name,  "recharge");
                var[i].dim = numele;
                break;
            case 21:
                strcpy(var[i].name,  "g");
                var[i].dim = numele;
                break;
            case 22:
                strcpy(var[i].name,  "sh");
                var[i].dim = numele;
                break;
            case 23:
                strcpy(var[i].name,  "le");
                var[i].dim = numele;
                break;
            case 24:
                strcpy(var[i].name, "t1");
                var[i].dim = numele;
                break;
            case 25:
                strcpy(var[i].name, "snowh");
                var[i].dim = numele;
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
            }

            if (i == 26 + MAXLYR + j)
            {
                sprintf(var[i].name, "smc%-d", j);
                var[i].dim = numele;
            }

            if (i == 26 + 2 * MAXLYR + j)
            {
                sprintf (var[i].name, "swc%-d", j);
                var[i].dim = numele;
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
}

void WriteParamFile (int rawtime, enkf_struct ens, int ind, char *outputdir)
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
        fprintf(fid, "\t%lf", ens->member[i].param[ind]);
    }

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
        sprintf (fn, "input/%s/%s.%3.3d.calib", project, project, i);
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
