#include "pihm.h"

void ReadCyclesCtrl(const char filen[], agtbl_struct *agtbl, ctrl_struct *ctrl)
{
    FILE           *fp;
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];
    int             i, n;
    int             match;
    int             index;
    int             lno = 0;

    /* Open simulation control file */
    fp = fopen(filen, "r");
    CheckFile(fp, filen);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filen);

    agtbl->oper = (int *)malloc(nelem * sizeof(int));

    /* Read simulation control file */
    FindLine(fp, "BOF", &lno, filen);
    /* Skip header line */
    NextLine(fp, cmdstr, &lno);

    for (i = 0; i < nelem; i++)
    {
        NextLine(fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%d %d", &index, &agtbl->oper[i]) != 2)
        {
            PIHMprintf(VL_ERROR,
                "Error reading information of the %dth element for Cycles.\n",
                i + 1);
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filen, lno);
            PIHMexit(EXIT_FAILURE);
        }
    }

    FindLine(fp, "OPERATION_FILE", &lno, filen);

    n = 0;
    while (1)
    {
        NextLine(fp, cmdstr, &lno);
        sscanf(cmdstr, "%s", optstr);

        if (strcasecmp(cmdstr, "EOF") == 0 ||
            strcasecmp(optstr, "PRINT_CTRL") == 0)
        {
            break;
        }

        match = sscanf(cmdstr, "%d %s", &index, agtbl->oper_filen[n]);
        if (match != 2 || n != index - 1)
        {
            PIHMprintf(VL_ERROR, "Error reading operation description.\n");
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filen, lno);
            PIHMexit(EXIT_FAILURE);
        }
        n++;
    }

    agtbl->noper = n;

    /* Output control */
    FindLine(fp, "BOF", &lno, filen);

#if TEMP_DISABLED
    FindLine(fp, "RESTART_CTRL", &lno, filen);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "READ_IC", &ctrl->read_cycles_restart, 'i', filen,
        lno);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "WRITE_IC", &ctrl->write_cycles_restart, 'i',
        filen, lno);
#endif

    FindLine(fp, "PRINT_CTRL", &lno, filen);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[BIOMASS_CTRL] = ReadPrtCtrl(cmdstr, "BIOMASS", filen, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[LAI_CTRL] = ReadPrtCtrl(cmdstr, "LAI", filen, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[RADNINTCP_CTRL] = ReadPrtCtrl(cmdstr, "RADN_INTCP",
        filen, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[WATER_STS_CTRL] = ReadPrtCtrl(cmdstr, "WATER_STRESS",
        filen, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[N_STS_CTRL] = ReadPrtCtrl(cmdstr, "N_STRESS", filen, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[CROP_TR_CTRL] = ReadPrtCtrl(cmdstr, "CROP_TRANSP", filen, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[CROP_POTTR_CTRL] = ReadPrtCtrl(cmdstr, "CROP_POT_TRANSP",
        filen, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[RES_EVAP_CTRL] = ReadPrtCtrl(cmdstr, "RES_EVAP", filen,
        lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[N_PROFILE_CTRL] = ReadPrtCtrl(cmdstr, "N_PROFILE", filen,
        lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[N_RIVER_CTRL] = ReadPrtCtrl(cmdstr, "N_RIVER",
        filen, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[DENITRIF_CTRL] = ReadPrtCtrl(cmdstr, "DENITRIF",
        filen, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[LEACHING_CTRL] = ReadPrtCtrl(cmdstr, "LEACHING",
        filen, lno);

    fclose(fp);
}

#if _CYCLES_OBSOLETE_
void ReadSoilInit(const char filename[], soiltbl_struct *soiltbl)
{
    FILE           *soil_fp;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             layer;
    int             i, k;
    int             lno = 0;

    /*
     * Open soil initialization file
     */
    soil_fp = fopen(filename, "r");
    CheckFile(soil_fp, filename);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filename);

    soiltbl->totalLayers = (int *)malloc(soiltbl->number * sizeof(int));
    soiltbl->clay_lyr = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->sand_lyr = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->iom_lyr = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->bd_lyr = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->no3_lyr = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->nh4_lyr = (double **)malloc(soiltbl->number * sizeof(double *));

    /* Read soil file */
    FindLine(soil_fp, "BOF", &lno, filename);

    for (i = 0; i < soiltbl->number; i++)
    {
        NextLine(soil_fp, cmdstr, &lno);
        ReadKeyword(cmdstr, "SOIL_TYPE", &index, 'i', filename, lno);

        if (i != index - 1)
        {
            PIHMprintf(VL_ERROR,
                "Error reading soil description of the %dth soil type.\n",
                i + 1);
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            PIHMexit(EXIT_FAILURE);
        }

        NextLine(soil_fp, cmdstr, &lno);
        ReadKeyword(cmdstr, "TOTAL_LAYERS", &soiltbl->totalLayers[i], 'i',
            filename, lno);

        soiltbl->clay_lyr[i] = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->sand_lyr[i] = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->iom_lyr[i] = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->bd_lyr[i] = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->no3_lyr[i] = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->nh4_lyr[i] = (double *)malloc(MAXLYR * sizeof(double));

        /* Skip header */
        NextLine(soil_fp, cmdstr, &lno);

        for (k = 0; k < soiltbl->totalLayers[i]; k++)
        {
            NextLine(soil_fp, cmdstr, &lno);
            match = sscanf(cmdstr, "%d %lf %lf %lf %lf %lf %lf", &layer,
                &soiltbl->clay_lyr[i][k], &soiltbl->sand_lyr[i][k],
                &soiltbl->iom_lyr[i][k], &soiltbl->bd_lyr[i][k],
                &soiltbl->no3_lyr[i][k], &soiltbl->nh4_lyr[i][k]);

            if (match != 7 || k != layer - 1)
            {
                PIHMprintf(VL_ERROR,
                    "Error reading description of the %dth layer of the %dth"
                    "soil type.\n", k + 1, i + 1);
                PIHMprintf(VL_ERROR,
                    "Error in %s near Line %d.\n", filename, lno);
                PIHMexit(EXIT_FAILURE);
            }

            soiltbl->clay_lyr[i][k] = (soiltbl->clay_lyr[i][k] < 0.0) ?
                soiltbl->clay[i] : soiltbl->clay_lyr[i][k];
            soiltbl->sand_lyr[i][k] = (soiltbl->sand_lyr[i][k] < 0.0) ?
                100.0 - soiltbl->clay[i] - soiltbl->silt[i] :
                soiltbl->sand_lyr[i][k];
            soiltbl->iom_lyr[i][k] = (soiltbl->iom_lyr[i][k] < 0.0) ?
                soiltbl->om[i] : soiltbl->iom_lyr[i][k];
            soiltbl->bd_lyr[i][k] = (soiltbl->bd_lyr[i][k] < 0.0) ?
                soiltbl->bd[i] : soiltbl->bd_lyr[i][k];
        }

        for (k = soiltbl->totalLayers[i] - 1; k < MAXLYR; k++)
        {
            soiltbl->clay_lyr[i][k] = soiltbl->clay_lyr[i][k - 1];
            soiltbl->sand_lyr[i][k] = soiltbl->sand_lyr[i][k - 1];
            soiltbl->iom_lyr[i][k] = soiltbl->iom_lyr[i][k - 1];
            soiltbl->bd_lyr[i][k] = soiltbl->bd_lyr[i][k - 1];
            soiltbl->no3_lyr[i][k] = soiltbl->no3_lyr[i][k - 1];
            soiltbl->nh4_lyr[i][k] = soiltbl->nh4_lyr[i][k - 1];
        }
    }

    fclose(soil_fp);
}

void ReadMultOper(const agtbl_struct *agtbl, const epconst_struct epctbl[],
    opertbl_struct opertbl[])
{
    int             i;
    char            filename[MAXSTRING];

    for (i = 0; i < agtbl->nopfile; i++)
    {
        sprintf(filename, "input/%s/%s", project, agtbl->opfilen[i]);

        ReadOperation(filename, 999, epctbl, &opertbl[i]);
    }
}
#endif

//int CropExist(char *cropName, const croptbl_struct *croptbl)
//{
//    int             i;
//    int             exist = 0;
//
//    for (i = 0; i < croptbl->number; i++)
//    {
//        if (strcmp(cropName, croptbl->cropName[i]) == 0)
//        {
//            exist = 1;
//            break;
//        }
//    }
//
//    return (exist);
//}
//
//void ReadCyclesIC(char *fn, elem_struct *elem, river_struct *riv)
//{
//    FILE           *init_file;
//    int             i;
//
//    init_file = fopen(fn, "rb");
//    CheckFile(init_file, fn);
//    PIHMprintf(VL_VERBOSE, " Reading %s\n", fn);
//
//    for (i = 0; i < nelem; i++)
//    {
//        fread(&elem[i].cycles_restart, sizeof(cyclesic_struct), 1, init_file);
//    }
//
//    for (i = 0; i < nriver; i++)
//    {
//        fread(&riv[i].cycles_restart, sizeof(river_cyclesic_struct), 1,
//            init_file);
//    }
//
//    fclose(init_file);
//}
