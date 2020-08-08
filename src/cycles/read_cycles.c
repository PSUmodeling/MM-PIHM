#include "pihm.h"

void ReadCyclesCtrl(const char fn[], agtbl_struct *agtbl, ctrl_struct *ctrl)
{
    FILE           *fp;
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];
    int             i, n;
    int             match;
    int             index;
    int             lno = 0;

    /* Open simulation control file */
    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    agtbl->oper = (int *)malloc(nelem * sizeof(int));

    /* Read simulation control file */
    FindLine(fp, "BOF", &lno, fn);
    /* Skip header line */
    NextLine(fp, cmdstr, &lno);

    for (i = 0; i < nelem; i++)
    {
        NextLine(fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%d %d", &index, &agtbl->oper[i]) != 2)
        {
            pihm_printf(VL_ERROR,
                "Error reading information of the %dth element for Cycles.\n",
                i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", fn, lno);
            pihm_exit(EXIT_FAILURE);
        }
    }

    FindLine(fp, "OPERATION_FILE", &lno, fn);

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
            pihm_printf(VL_ERROR, "Error reading operation description.\n");
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", fn, lno);
            pihm_exit(EXIT_FAILURE);
        }
        n++;
    }

    agtbl->noper = n;

    /* Output control */
    FindLine(fp, "BOF", &lno, fn);

#if TEMP_DISABLED
    FindLine(fp, "RESTART_CTRL", &lno, fn);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "READ_IC", 'i', fn, lno, &ctrl->read_cycles_restart);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "WRITE_IC", 'i', fn, lno, &ctrl->write_cycles_restart);
#endif

    FindLine(fp, "PRINT_CTRL", &lno, fn);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[YIELD_CTRL] = ReadPrintCtrl(cmdstr, "YIELD", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[BIOMASS_CTRL] = ReadPrintCtrl(cmdstr, "BIOMASS", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[LAI_CTRL] = ReadPrintCtrl(cmdstr, "LAI", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[RADNINTCP_CTRL] = ReadPrintCtrl(cmdstr, "RADN_INTCP",
        fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[WATER_STS_CTRL] = ReadPrintCtrl(cmdstr, "WATER_STRESS",
        fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[N_STS_CTRL] = ReadPrintCtrl(cmdstr, "N_STRESS", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[CROP_TR_CTRL] = ReadPrintCtrl(cmdstr, "CROP_TRANSP", fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[CROP_POTTR_CTRL] = ReadPrintCtrl(cmdstr, "CROP_POT_TRANSP",
        fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[N_PROFILE_CTRL] = ReadPrintCtrl(cmdstr, "N_PROFILE", fn,
        lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[N_RIVER_CTRL] = ReadPrintCtrl(cmdstr, "N_RIVER",
        fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[DENITRIF_CTRL] = ReadPrintCtrl(cmdstr, "DENITRIF",
        fn, lno);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[LEACHING_CTRL] = ReadPrintCtrl(cmdstr, "LEACHING",
        fn, lno);

    fclose(fp);
}

void ReadSoilInit(const char fn[], soiltbl_struct *soiltbl)
{
    FILE           *fp;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             i, kz;
    int             lno = 0;

    /*
     * Open soil initialization file
     */
    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    soiltbl->nlayers    = (int *)malloc(soiltbl->number * sizeof(int));
    soiltbl->clay_layer = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->sand_layer = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->om_layer   = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->bd_layer   = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->no3        = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->nh4        = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->fc         = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->pwp        = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->air_entry_pot
                        = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->b          = (double **)malloc(soiltbl->number * sizeof(double *));

    FindLine(fp, "BOF", &lno, fn);

    /* Read soil information for the ith soil type */
    for (i = 0; i < soiltbl->number; i++)
    {
        NextLine(fp, cmdstr, &lno);
        ReadKeyword(cmdstr, "SOIL_TYPE", 'i', fn, lno, &index);

        if (i != index - 1)
        {
            pihm_printf(VL_ERROR,
                "Error reading soil description of the %dth soil type.\n",
                i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", fn, lno);
            pihm_exit(EXIT_FAILURE);
        }

        NextLine(fp, cmdstr, &lno);
        ReadKeyword(cmdstr, "TOTAL_LAYERS", 'i', fn, lno,
            &soiltbl->nlayers[i]);

        soiltbl->clay_layer[i]    = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->sand_layer[i]    = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->om_layer[i]      = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->bd_layer[i]      = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->no3[i]           = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->nh4[i]           = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->fc[i]            = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->pwp[i]           = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->air_entry_pot[i] = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->b[i]             = (double *)malloc(MAXLYR * sizeof(double));

        /* Skip header line */
        NextLine(fp, cmdstr, &lno);

        for (kz = 0; kz < soiltbl->nlayers[i]; kz++)
        {
            int             layer;
            double          bd;             /* Saxton's bulk density */
            double          wc33;           /* Saxton's volumetric WC at 33 J/kg
                                             */
            double          wc1500;         /* Saxton's volumetric WC at
                                             * 1500 J/kg */
            double          fc_water_pot;

            NextLine(fp, cmdstr, &lno);
            match = sscanf(cmdstr, "%d %lf %lf %lf %lf %lf %lf", &layer,
                &soiltbl->clay_layer[i][kz], &soiltbl->sand_layer[i][kz],
                &soiltbl->om_layer[i][kz], &soiltbl->bd_layer[i][kz],
                &soiltbl->no3[i][kz], &soiltbl->nh4[i][kz]);

            if (match != 7 || kz != layer - 1)
            {
                pihm_printf(VL_ERROR,
                    "Error reading description of the %dth layer of the %dth"
                    "soil type.\n", kz + 1, i + 1);
                pihm_printf(VL_ERROR,
                    "Error in %s near Line %d.\n", fn, lno);
                pihm_exit(EXIT_FAILURE);
            }

            soiltbl->clay_layer[i][kz] = (soiltbl->clay_layer[i][kz] < 0.0) ?
                soiltbl->clay[i] : soiltbl->clay_layer[i][kz];
            soiltbl->clay_layer[i][kz] /= 100.0;
            soiltbl->sand_layer[i][kz] = (soiltbl->sand_layer[i][kz] < 0.0) ?
                100.0 - soiltbl->clay[i] - soiltbl->silt[i] :
                soiltbl->sand_layer[i][kz];
            soiltbl->sand_layer[i][kz] /= 100.0;
            soiltbl->om_layer[i][kz] = (soiltbl->om_layer[i][kz] < 0.0) ?
                soiltbl->om[i] : soiltbl->om_layer[i][kz];

            bd = BulkDensity(soiltbl->clay_layer[i][kz],
                soiltbl->sand_layer[i][kz], soiltbl->om_layer[i][kz]);
            wc33 = VolWCAt33Jkg(soiltbl->clay_layer[i][kz],
                soiltbl->sand_layer[i][kz], soiltbl->om_layer[i][kz]);
            wc1500 = VolWCAt1500Jkg(soiltbl->clay_layer[i][kz],
                soiltbl->sand_layer[i][kz], soiltbl->om_layer[i][kz]);

            /* Saxton's b */
            soiltbl->b[i][kz] = (log(1500.0) - log(33.0)) /
                (log(wc33) - log(wc1500));

            /* Air entry potential */
            soiltbl->air_entry_pot[i][kz] = -33.0 *
                pow(wc33 / (1.0 - bd / 2.65), soiltbl->b[i][kz]) *
                ((roundi(soiltbl->bd_layer[i][kz]) == BADVAL) ?
                1.0 :
                pow(soiltbl->bd_layer[i][kz] / bd, 0.67 * soiltbl->b[i][kz]));

            /* Bulk density */
            soiltbl->bd_layer[i][kz] =
                (roundi(soiltbl->bd_layer[i][kz]) == BADVAL) ?
                bd : soiltbl->bd_layer[i][kz];

            /* Field capacity */
            fc_water_pot = -0.35088 * soiltbl->clay_layer[i][kz] * 100.0 -
                28.947;
            soiltbl->fc[i][kz] = SoilWaterContent(soiltbl->smcmax[i],
                soiltbl->air_entry_pot[i][kz], soiltbl->b[i][kz], fc_water_pot);

            /* Permanent wilting point */
            soiltbl->pwp[i][kz] = SoilWaterContent(soiltbl->smcmax[i],
                soiltbl->air_entry_pot[i][kz], soiltbl->b[i][kz], -1500.0);
        }

        /* When the number of soil layers in Flux-PIHM is larger than in
         * soilinit file, use the last described layer to fill the rest of
         * layers */
        if (soiltbl->nlayers[i] < MAXLYR)
        {
            for (kz = soiltbl->nlayers[i]; kz < MAXLYR; kz++)
            {
                soiltbl->clay_layer[i][kz]    = soiltbl->clay_layer[i][kz - 1];
                soiltbl->sand_layer[i][kz]    = soiltbl->sand_layer[i][kz - 1];
                soiltbl->om_layer[i][kz]      = soiltbl->om_layer[i][kz - 1];
                soiltbl->bd_layer[i][kz]      = soiltbl->bd_layer[i][kz - 1];
                soiltbl->no3[i][kz]           = soiltbl->no3[i][kz - 1];
                soiltbl->nh4[i][kz]           = soiltbl->nh4[i][kz - 1];
                soiltbl->b[i][kz]             = soiltbl->b[i][kz - 1];
                soiltbl->air_entry_pot[i][kz] =
                    soiltbl->air_entry_pot[i][kz - 1];
                soiltbl->fc[i][kz]            = soiltbl->fc[i][kz - 1];
                soiltbl->pwp[i][kz]           = soiltbl->pwp[i][kz - 1];
            }
        }
    }

    fclose(fp);
}

void ReadMultOper(const agtbl_struct *agtbl, mgmt_struct mgmttbl[],
    crop_struct croptbl[])
{
    int             i;
    char            fn[MAXSTRING];

    for (i = 0; i < agtbl->noper; i++)
    {
        sprintf(fn, "input/%s/%s", project, agtbl->oper_filen[i]);

        ReadOper(fn, BADVAL, BADVAL, &mgmttbl[i], croptbl);
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
//    init_file = pihm_fopen(fn, "rb");
//    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);
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
