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
    fp = pihm_fopen(filen, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filen);

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
            pihm_printf(VL_ERROR,
                "Error reading information of the %dth element for Cycles.\n",
                i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filen, lno);
            pihm_exit(EXIT_FAILURE);
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
            pihm_printf(VL_ERROR, "Error reading operation description.\n");
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filen, lno);
            pihm_exit(EXIT_FAILURE);
        }
        n++;
    }

    agtbl->noper = n;

    /* Output control */
    FindLine(fp, "BOF", &lno, filen);

#if TEMP_DISABLED
    FindLine(fp, "RESTART_CTRL", &lno, filen);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "READ_IC", 'i', filen, lno, &ctrl->read_cycles_restart);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "WRITE_IC", 'i', filen, lno, &ctrl->write_cycles_restart);
#endif

    FindLine(fp, "PRINT_CTRL", &lno, filen);

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[YIELD_CTRL] = ReadPrtCtrl(cmdstr, "YIELD", filen, lno);

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

void ReadSoilInit(const char filen[], soiltbl_struct *soiltbl)
{
    FILE           *fp;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             i, k;
    int             lno = 0;

    /*
     * Open soil initialization file
     */
    fp = pihm_fopen(filen, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filen);

    soiltbl->nlayers    = (int *)malloc(soiltbl->number * sizeof(int));
    soiltbl->clay_layer = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->sand_layer = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->om_layer  = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->bd_layer   = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->no3        = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->nh4        = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->fc         = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->pwp        = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->air_entry_pot
                        = (double **)malloc(soiltbl->number * sizeof(double *));
    soiltbl->b          = (double **)malloc(soiltbl->number * sizeof(double *));

    FindLine(fp, "BOF", &lno, filen);

    /* Read soil information for the ith soil type */
    for (i = 0; i < soiltbl->number; i++)
    {
        NextLine(fp, cmdstr, &lno);
        ReadKeyword(cmdstr, "SOIL_TYPE", 'i', filen, lno, &index);

        if (i != index - 1)
        {
            pihm_printf(VL_ERROR,
                "Error reading soil description of the %dth soil type.\n",
                i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filen, lno);
            pihm_exit(EXIT_FAILURE);
        }

        NextLine(fp, cmdstr, &lno);
        ReadKeyword(cmdstr, "TOTAL_LAYERS", 'i', filen, lno, &soiltbl->nlayers[i]);

        soiltbl->clay_layer[i]    = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->sand_layer[i]    = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->om_layer[i]     = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->bd_layer[i]      = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->no3[i]           = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->nh4[i]           = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->fc[i]            = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->pwp[i]           = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->air_entry_pot[i] = (double *)malloc(MAXLYR * sizeof(double));
        soiltbl->b[i]             = (double *)malloc(MAXLYR * sizeof(double));

        /* Skip header line */
        NextLine(fp, cmdstr, &lno);

        for (k = 0; k < soiltbl->nlayers[i]; k++)
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
                &soiltbl->clay_layer[i][k], &soiltbl->sand_layer[i][k],
                &soiltbl->om_layer[i][k], &soiltbl->bd_layer[i][k],
                &soiltbl->no3[i][k], &soiltbl->nh4[i][k]);

            if (match != 7 || k != layer - 1)
            {
                pihm_printf(VL_ERROR,
                    "Error reading description of the %dth layer of the %dth"
                    "soil type.\n", k + 1, i + 1);
                pihm_printf(VL_ERROR,
                    "Error in %s near Line %d.\n", filen, lno);
                pihm_exit(EXIT_FAILURE);
            }

            soiltbl->clay_layer[i][k] = (soiltbl->clay_layer[i][k] < 0.0) ?
                soiltbl->clay[i] : soiltbl->clay_layer[i][k];
            soiltbl->clay_layer[i][k] /= 100.0;
            soiltbl->sand_layer[i][k] = (soiltbl->sand_layer[i][k] < 0.0) ?
                100.0 - soiltbl->clay[i] - soiltbl->silt[i] :
                soiltbl->sand_layer[i][k];
            soiltbl->sand_layer[i][k] /= 100.0;
            soiltbl->om_layer[i][k] = (soiltbl->om_layer[i][k] < 0.0) ?
                soiltbl->om[i] : soiltbl->om_layer[i][k];

            bd = BulkDensity(soiltbl->clay_layer[i][k],
                soiltbl->sand_layer[i][k], soiltbl->om_layer[i][k]);
            wc33 = VolWCAt33Jkg(soiltbl->clay_layer[i][k],
                soiltbl->sand_layer[i][k], soiltbl->om_layer[i][k]);
            wc1500 = VolWCAt1500Jkg(soiltbl->clay_layer[i][k],
                soiltbl->sand_layer[i][k], soiltbl->om_layer[i][k]);

            /* Saxton's b */
            soiltbl->b[i][k] = (log(1500.0) - log(33.0)) /
                (log(wc33) - log(wc1500));

            /* Air entry potential */
            soiltbl->air_entry_pot[i][k] = -33.0 *
                pow(wc33 / (1.0 - bd / 2.65), soiltbl->b[i][k]) *
                ((roundi(soiltbl->bd_layer[i][k]) == BADVAL) ?
                1.0 :
                pow(soiltbl->bd_layer[i][k] / bd, 0.67 * soiltbl->b[i][k]));

            /* Bulk density */
            soiltbl->bd_layer[i][k] =
                (roundi(soiltbl->bd_layer[i][k]) == BADVAL) ?
                bd : soiltbl->bd_layer[i][k];

            /* Field capacity */
            fc_water_pot = -0.35088 * soiltbl->clay_layer[i][k] * 100.0 -
                28.947;
            soiltbl->fc[i][k] = SoilWaterContent(soiltbl->smcmax[i],
                soiltbl->air_entry_pot[i][k], soiltbl->b[i][k], fc_water_pot);

            /* Permanent wilting point */
            soiltbl->pwp[i][k] = SoilWaterContent(soiltbl->smcmax[i],
                soiltbl->air_entry_pot[i][k], soiltbl->b[i][k], -1500.0);
        }

        /* When the number of soil layers in Flux-PIHM is larger than in
         * soilinit file, use the last described layer to fill the rest of
         * layers */
        if (soiltbl->nlayers[i] < MAXLYR)
        {
            for (k = soiltbl->nlayers[i]; k < MAXLYR; k++)
            {
                soiltbl->clay_layer[i][k]    = soiltbl->clay_layer[i][k - 1];
                soiltbl->sand_layer[i][k]    = soiltbl->sand_layer[i][k - 1];
                soiltbl->om_layer[i][k]     = soiltbl->om_layer[i][k - 1];
                soiltbl->bd_layer[i][k]      = soiltbl->bd_layer[i][k - 1];
                soiltbl->no3[i][k]           = soiltbl->no3[i][k - 1];
                soiltbl->nh4[i][k]           = soiltbl->nh4[i][k - 1];
                soiltbl->b[i][k]             = soiltbl->b[i][k - 1];
                soiltbl->air_entry_pot[i][k] = soiltbl->air_entry_pot[i][k - 1];
                soiltbl->fc[i][k]            = soiltbl->fc[i][k - 1];
                soiltbl->pwp[i][k]           = soiltbl->pwp[i][k - 1];
            }
        }
    }

    fclose(fp);
}

void ReadMultOper(const agtbl_struct *agtbl, mgmt_struct mgmttbl[],
    crop_struct croptbl[])
{
    int             i;
    char            filen[MAXSTRING];

    for (i = 0; i < agtbl->noper; i++)
    {
        sprintf(filen, "input/%s/%s", project, agtbl->oper_filen[i]);

        ReadOper(filen, BADVAL, BADVAL, &mgmttbl[i], croptbl);
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
