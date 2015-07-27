#include "pihm.h"  
#include "enkf.h"

void EnKFRead(char *project, enkf_struct ens)
{
    char            fn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    FILE           *enkf_file;
    int             i;
    int             n;
    ctrl_struct     ctrl;

    //int i, j;
    //int tempindex;

    //time_t rawtime;
    //struct tm *timeinfo;

    //int NumTout;
    //char tempchar[10];


    //int Verbose, Debug, init_type, gwD, surfD, snowD, rivStg, Rech, IsD, usD, et0, et1, et2, lsv, rivFlx[10];
    //int gwDInt, surfDInt, snowDInt, rivStgInt, RechInt, IsDInt, usDInt, etInt, rivFlxInt, UnsatMode, SurfMode, RivMode, Solver, GSType, MaxK;
    //realtype delt, abstol, reltol, InitStep, MaxStep, ETStep, a, b;
    //int outtype;

    //timeinfo = (struct tm *)malloc(sizeof(struct tm));

    /*========== open *.enkf file ==========*/ 
    sprintf (fn, "input/%s/%s.enkf", project, project);
    enkf_file = fopen (fn, "r");
    CheckFile (enkf_file, fn);

    FindLine (enkf_file, "BOF");

    NextLine (enkf_file, cmdstr);
    ReadKeywordInt (cmdstr, "NUM_ENSEMBLE_MEMBER", &ens->ne);

    NextLine (enkf_file, cmdstr);
    ReadKeywordInt (cmdstr, "ASSIMILATION_INTERVAL", &ens->interval);

    NextLine (enkf_file, cmdstr);
    ReadKeywordInt (cmdstr, "START_MODE", &ens->mbr_start_mode);

    NextLine (enkf_file, cmdstr);
    ReadKeywordDouble (cmdstr, "INFLATION_WEIGHT", &ens->weight);

    NextLine (enkf_file, cmdstr);
    ReadKeywordTime (cmdstr, "ASSIMILATION_END_TIME", &ens->end_time);

    FindLine (enkf_file, "PARAMETER");
    n = CountLine (enkf_file, 1, "NUM_OBS");

    /* Rewind to read */
    FindLine (enkf_file, "PARAMETER");

    /* Start reading EnKF file */ 
    for (i = 0; i < n; i++)
    {
        NextLine (enkf_file, cmdstr);
        sscanf (cmdstr, "%s %d %d %lf %lf %lf %lf %d", 
            ens->param[i].name,
            &ens->param[i].perturb, &ens->param[i].update,
            &ens->param[i].perturb_min, &ens->param[i].perturb_max,
            &ens->param[i].min, &ens->param[i].max,
            &ens->param[i].type);
    }

    //fscanf(EnKF_file,"%d", &En->no_obs);

    //if (En->no_obs>0)
    //{
    //    En->observation = (ObsCtrl *)malloc(En->no_obs*sizeof(ObsCtrl));
    //    for (i=0; i<En->no_obs; i++)
    //    {
    //        fscanf(EnKF_file, "%s %s %d %d", &En->observation[i].obsn, &En->observation[i].obsfn, &En->observation[i].obs_type, &En->observation[i].no_ctrl);
    //        En->observation[i].vrbl_ind = (int *)malloc(En->observation[i].no_ctrl*sizeof(int));
    //        En->observation[i].grid_ind = (int *)malloc(En->observation[i].no_ctrl*sizeof(int));
    //        En->observation[i].weight = (realtype *)malloc(En->observation[i].no_ctrl*sizeof(realtype));
    //        for (j=0; j<En->observation[i].no_ctrl; j++)
    //        {
    //            fscanf(EnKF_file, "%d %d %lf", &En->observation[i].vrbl_ind[j], &En->observation[i].grid_ind[j], &En->observation[i].weight[j]);
    //        }
    //    }
    //}

    fclose (enkf_file);

    /*
     * Read .para to obtain the first cycle's start and end times
     */
    ReadPara (project, &ctrl);
    ens->mbr_start_mode = ctrl.init_type;
    ens->cycle_start_time = ctrl.starttime;
    ens->cycle_end_time = ctrl.endtime;
}

//void read_para_en(char *project, Model_Data DS, Control_Data *CS,realtype StartTime, realtype EndTime, int init_type)
//{
//    int i, j;
//    int tempindex;
//
//    struct tm *timeinfo;
//
//    int NumTout;
//    char *fn;
//    char tempchar[10];
//
//    FILE *para_file;
//
//    timeinfo = (struct tm *)malloc(sizeof(struct tm));
//
//    fn = (char *)malloc((strlen(project)+12)*sizeof(char));
//    strcpy(fn,"input/");
//    strcat(fn, project);
//    para_file = fopen(strcat(fn, ".para"), "r");  
//
//    if(para_file == NULL)
//    {
//        printf("\n  Fatal Error: %s.para is in use or does not exist!\n", project);
//        exit(1);
//    }
//    fscanf(para_file, "%d %d", &(CS->Verbose), &(CS->Debug));
//    fscanf(para_file, "%d", &tempindex);
//    CS->init_type = init_type;
//    fscanf(para_file, "%d %d %d %d", &(CS->gwD), &(CS->surfD), &(CS->snowD), &(CS->rivStg));
//    fscanf(para_file, "%d %d %d", &(CS->Rech), &(CS->IsD), &(CS->usD));
//    fscanf(para_file, "%d %d %d %d", &(CS->et[0]), &(CS->et[1]), &(CS->et[2]), &(CS->lsv));
//    for(i=0;i<10;i++)
//    {
//        fscanf(para_file, "%d", &(CS->rivFlx[i]));
//    }
//    fscanf(para_file, "%d %d %d %d", &(CS->gwDInt), &(CS->surfDInt), &(CS->snowDInt), &(CS->rivStgInt));
//    fscanf(para_file, "%d %d %d %d", &(CS->RechInt), &(CS->IsDInt), &(CS->usDInt), &(CS->etInt));
//    fscanf(para_file, "%d",&(CS->rivFlxInt));
//
//    fscanf(para_file, "%d %d %d", &DS->UnsatMode, &DS->SurfMode, &DS->RivMode);
//
//    fscanf(para_file, "%d", &(CS->Solver));
//    if(CS->Solver == 2)
//    {
//        fscanf(para_file, "%d %d %lf", &CS->GSType, &CS->MaxK, &CS->delt);
//    }
//    fscanf(para_file, "%lf %lf", &(CS->abstol), &(CS->reltol));
//    fscanf(para_file, "%lf %lf %lf", &(CS->InitStep), &(CS->MaxStep), &(CS->ETStep));
//    fscanf(para_file, "%d-%d-%d %d:%d", &timeinfo->tm_year,&timeinfo->tm_mon,&timeinfo->tm_mday,&timeinfo->tm_hour,&timeinfo->tm_min);
//    CS->StartTime = StartTime;
//    fscanf(para_file, "%d-%d-%d %d:%d", &timeinfo->tm_year,&timeinfo->tm_mon,&timeinfo->tm_mday,&timeinfo->tm_hour,&timeinfo->tm_min);
//    CS->EndTime = EndTime;
//
//    //	printf("start = %f, end = %f",CS->StartTime,CS->EndTime);
//    fscanf(para_file, "%d", &(CS->outtype));
//    if(CS->outtype == 0)
//    {
//        fscanf(para_file, "%lf %lf", &CS->a, &CS->b);
//    }
//
//    if(CS->a != 1.0)
//    {
//        NumTout = (int)(log(1 - (CS->EndTime - CS->StartTime)*(1 -  CS->a)/CS->b)/log(CS->a));
//    }
//    else
//    {
//        if((CS->EndTime - CS->StartTime)/CS->b - ((int) (CS->EndTime - CS->StartTime)/CS->b) > 0)
//        {
//            NumTout = (int) ((CS->EndTime - CS->StartTime)/CS->b);
//        }
//        else
//        {
//            NumTout = (int) ((CS->EndTime - CS->StartTime)/CS->b - 1);
//        }  
//    }
//
//    CS->NumSteps = NumTout + 1;
//
//    CS->Tout = (realtype *)malloc((CS->NumSteps + 1)*sizeof(realtype));
//
//    for(i=0; i<CS->NumSteps+1; i++)
//    {
//        if(i == 0)
//        {
//            CS->Tout[i] = CS->StartTime;
//        }
//        else
//        {
//            CS->Tout[i] = CS->Tout[i-1] + pow(CS->a, i)*CS->b;
//        }  
//    }
//
//    if(CS->Tout[CS->NumSteps] < CS->EndTime)
//    {
//        CS->Tout[CS->NumSteps] = CS->EndTime;
//    }
//
//    fclose(para_file);
//    free(fn);
//    free(timeinfo);
//    //	printf("done.\n");
//
//}

