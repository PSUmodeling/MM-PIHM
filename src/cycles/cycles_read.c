#include "pihm.h"

void ReadCyclesCtrl(const char filename[], agtbl_struct *agtbl,
    ctrl_struct *ctrl)
{
    FILE           *cycles_fp;
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];
    int             i;
    int             match;
    int             index;
    int             lno = 0;

    /* Open simulation control file */
    cycles_fp = fopen(filename, "r");
    CheckFile(cycles_fp, filename);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filename);

    agtbl->op = (int *)malloc(nelem * sizeof(int));
    agtbl->rotsz = (int *)malloc(nelem * sizeof(int));
    agtbl->auto_N = (int *)malloc(nelem * sizeof(int));
    agtbl->auto_P = (int *)malloc(nelem * sizeof(int));
    agtbl->auto_S = (int *)malloc(nelem * sizeof(int));

    /* Read simulation control file */
    FindLine(cycles_fp, "BOF", &lno, filename);

    NextLine(cycles_fp, cmdstr, &lno);
    for (i = 0; i < nelem; i++)
    {
        NextLine(cycles_fp, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %d %d %d %d %d", &index, &agtbl->op[i],
            &agtbl->rotsz[i], &agtbl->auto_N[i], &agtbl->auto_P[i],
            &agtbl->auto_S[i]);
        if (match != 6)
        {
            PIHMprintf(VL_ERROR,
                "Error reading information of the %dth element for Cycles.\n",
                i + 1);
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            PIHMexit(EXIT_FAILURE);
        }
    }

    FindLine(cycles_fp, "OPERATION_FILE", &lno, filename);

    i = 0;
    while (1)
    {
        NextLine(cycles_fp, cmdstr, &lno);
        sscanf(cmdstr, "%s", optstr);

        if (strcasecmp(cmdstr, "EOF") == 0 ||
            strcasecmp(optstr, "PRINT_CTRL") == 0)
        {
            break;
        }

        match = sscanf(cmdstr, "%d %s", &index, agtbl->opfilen[i]);
        if (match != 2 || i != index - 1)
        {
            PIHMprintf(VL_ERROR, "Error reading operation description.\n");
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            PIHMexit(EXIT_FAILURE);
        }
        i++;
    }

    agtbl->nopfile = i;

    /* Output control */
    FindLine(cycles_fp, "BOF", &lno, filename);

    FindLine(cycles_fp, "RESTART_CTRL", &lno, filename);

    NextLine(cycles_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "READ_IC", &ctrl->read_cycles_restart, 'i', filename,
        lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "WRITE_IC", &ctrl->write_cycles_restart, 'i',
        filename, lno);

    FindLine(cycles_fp, "PRINT_CTRL", &lno, filename);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[BIOMASS_CTRL] = ReadPrtCtrl(cmdstr, "BIOMASS", filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[LAI_CTRL] = ReadPrtCtrl(cmdstr, "LAI", filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[RADNINTCP_CTRL] = ReadPrtCtrl(cmdstr, "RADN_INTCP",
        filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[WATER_STS_CTRL] = ReadPrtCtrl(cmdstr, "WATER_STRESS",
        filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[N_STS_CTRL] = ReadPrtCtrl(cmdstr, "N_STRESS", filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[CROP_TR_CTRL] = ReadPrtCtrl(cmdstr, "CROP_TR", filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[CROP_POTTR_CTRL] = ReadPrtCtrl(cmdstr, "CROP_POT_TR",
        filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[RES_EVAP_CTRL] = ReadPrtCtrl(cmdstr, "RES_EVAP", filename,
        lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[NO3_PROF_CTRL] = ReadPrtCtrl(cmdstr, "NO3_PROF", filename,
        lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[NO3_RIVER_CTRL] = ReadPrtCtrl(cmdstr, "NO3_RIVER",
        filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[NH4_PROF_CTRL] = ReadPrtCtrl(cmdstr, "NH4_PROF", filename,
        lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[NH4_RIVER_CTRL] = ReadPrtCtrl(cmdstr, "NH4_RIVER",
        filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[NO3_DENIT_CTRL] = ReadPrtCtrl(cmdstr, "NO3_DENITRIF",
        filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[NO3_LEACH_CTRL] = ReadPrtCtrl(cmdstr, "NO3_LEACH",
        filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[NH4_LEACH_CTRL] = ReadPrtCtrl(cmdstr, "NH4_LEACH",
        filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[NO3_LEACH_RIVER_CTRL] = ReadPrtCtrl(cmdstr,
        "NO3_LEACH_RIVER", filename, lno);

    NextLine(cycles_fp, cmdstr, &lno);
    ctrl->prtvrbl[NH4_LEACH_RIVER_CTRL] = ReadPrtCtrl(cmdstr,
        "NH4_LEACH_RIVER", filename, lno);

    fclose(cycles_fp);
}

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
