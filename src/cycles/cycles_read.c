#include "pihm.h"

void ReadCyclesCtrl (char *filename, agtbl_struct *agtbl, ctrl_struct *ctrl, int numele)
{
    FILE           *simctrl_file;
    //time_t          rawtime;
    //struct tm      *timestamp;
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];
    int             i;
    int             match;
    int             index;

    /* Open simulation control file */
    simctrl_file = fopen (filename, "r");

    if (NULL == simctrl_file)
    {
        fprintf (stderr, "Error opening %s.\n", filename);
        PIHMError (1);
    }

    if (verbose_mode)
    {
        printf ("Reading %s...\n", filename);
    }

    agtbl->op = (int *)malloc (numele * sizeof (int));
    agtbl->rotsz = (int *)malloc (numele * sizeof (int));
    agtbl->auto_N = (int *)malloc (numele * sizeof (int));
    agtbl->auto_P = (int *)malloc (numele * sizeof (int));
    agtbl->auto_S = (int *)malloc (numele * sizeof (int));

    /* Read simulation control file */
    FindLine (simctrl_file, "BOF");

    NextLine (simctrl_file, cmdstr);
    for (i = 0; i < numele; i++)
    {
        NextLine (simctrl_file, cmdstr);
        match =
            sscanf (cmdstr, "%d %d %d %d %d %d", &index, &agtbl->op[i],
            &agtbl->rotsz[i], &agtbl->auto_N[i], &agtbl->auto_P[i],
            &agtbl->auto_S[i]);
        if (match != 6)
        {
            printf ("Cannot read information of the %dth element!\n", i + 1);
            printf (".cycles file format error!\n");
            exit (1);
        }
    }

    FindLine (simctrl_file, "OPERATION_FILE");

    i = 0;
    while (1)
    {
        NextLine (simctrl_file, cmdstr);
        sscanf (cmdstr, "%s", optstr);

        if (strcasecmp (cmdstr, "EOF") == 0 || strcasecmp (optstr, "PRINT_CTRL") == 0)
        {
            break;
        }

        match = sscanf (cmdstr, "%d %s", &index, agtbl->opfilen[i]);
        if (match != 2 || i != index - 1)
        {
            fprintf (stderr, "Error reading %s.\n", filename),
            fprintf (stderr, "Please check file format.\n");
            PIHMError (1);
        }
        i++;
    }

    agtbl->nopfile = i;

    /* Output control */
    FindLine (simctrl_file, "PRINT_CTRL");

    NextLine (simctrl_file, cmdstr);
    ReadKeyword (cmdstr, "BIOMASS", &ctrl->prtvrbl[BIOMASS_CTRL], 'i');

    NextLine (simctrl_file, cmdstr);
    ReadKeyword (cmdstr, "RADN_INTCP", &ctrl->prtvrbl[RADNINTCP_CTRL], 'i');

    NextLine (simctrl_file, cmdstr);
    ReadKeyword (cmdstr, "WATER_STRESS", &ctrl->prtvrbl[WATER_STS_CTRL], 'i');

    NextLine (simctrl_file, cmdstr);
    ReadKeyword (cmdstr, "N_STRESS", &ctrl->prtvrbl[N_STS_CTRL], 'i');

    NextLine (simctrl_file, cmdstr);
    ReadKeyword (cmdstr, "CROP_TR", &ctrl->prtvrbl[CROP_TR_CTRL], 'i');

    NextLine (simctrl_file, cmdstr);
    ReadKeyword (cmdstr, "CROP_POT_TR", &ctrl->prtvrbl[CROP_POTTR_CTRL], 'i');

    NextLine (simctrl_file, cmdstr);
    ReadKeyword (cmdstr, "RES_EVAP", &ctrl->prtvrbl[RES_EVAP_CTRL], 'i');

    NextLine (simctrl_file, cmdstr);
    ReadKeyword (cmdstr, "NO3_PROF", &ctrl->prtvrbl[NO3_PROF_CTRL], 'i');

    NextLine (simctrl_file, cmdstr);
    ReadKeyword (cmdstr, "NO3_RIVER", &ctrl->prtvrbl[NO3_RIVER_CTRL], 'i');

    NextLine (simctrl_file, cmdstr);
    ReadKeyword (cmdstr, "NH4_PROF", &ctrl->prtvrbl[NH4_PROF_CTRL], 'i');

    NextLine (simctrl_file, cmdstr);
    ReadKeyword (cmdstr, "NH4_RIVER", &ctrl->prtvrbl[NH4_RIVER_CTRL], 'i');


    fclose (simctrl_file);

    //rawtime = (int)ctrl->starttime;
    //timestamp = gmtime (&rawtime);
    //ctrl->simStartYear = timestamp->tm_year + 1900;

    //rawtime = (int)ctrl->endtime;
    //timestamp = gmtime (&rawtime);
    //ctrl->simEndYear = timestamp->tm_year + 1900 - 1;

    //ctrl->totalYears = ctrl->simEndYear - ctrl->simStartYear + 1;
}

void ReadSoilInit (char *filename, soiltbl_struct *soiltbl)
{
    FILE           *soil_file;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             layer;
    int             i, j;

    /* 
     * Open soil initialization file
     */
    soil_file = fopen (filename, "r");

    if (NULL == soil_file)
    {
        fprintf (stderr, "Error opening %s.\n", filename);
        PIHMError (1);
    }

    if (verbose_mode)
    {
        printf ("Reading %s...\n", filename);
    }

    soiltbl->totalLayers = (int *)malloc (soiltbl->number * sizeof (int));
    soiltbl->clay_lyr =
        (double **)malloc (soiltbl->number * sizeof (double *));
    soiltbl->sand_lyr =
        (double **)malloc (soiltbl->number * sizeof (double *));
    soiltbl->iom_lyr =
        (double **)malloc (soiltbl->number * sizeof (double *));
    soiltbl->bd_lyr = (double **)malloc (soiltbl->number * sizeof (double *));
    soiltbl->NO3_lyr =
        (double **)malloc (soiltbl->number * sizeof (double *));
    soiltbl->NH4_lyr =
        (double **)malloc (soiltbl->number * sizeof (double *));

    /* Read soil file */
    FindLine (soil_file, "BOF");

    for (i = 0; i < soiltbl->number; i++)
    {
        NextLine (soil_file, cmdstr);
        ReadKeyword (cmdstr, "SOIL_TYPE", &index, 'i');

        if (i != index - 1)
        {
            printf ("Cannot read information of the %dth soil type!\n",
                i + 1);
            printf (".soilinit file format error!\n");
            exit (1);
        }

        NextLine (soil_file, cmdstr);
        ReadKeyword (cmdstr, "TOTAL_LAYERS", &soiltbl->totalLayers[i], 'i');

        soiltbl->clay_lyr[i] =
            (double *)malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->sand_lyr[i] =
            (double *)malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->iom_lyr[i] =
            (double *)malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->bd_lyr[i] =
            (double *)malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->NO3_lyr[i] =
            (double *)malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->NH4_lyr[i] =
            (double *)malloc (soiltbl->totalLayers[i] * sizeof (double));

        /* Skip header */
        NextLine (soil_file, cmdstr);

        for (j = 0; j < soiltbl->totalLayers[i]; j++)
        {
            NextLine (soil_file, cmdstr);
            match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf %lf", &layer,
                &soiltbl->clay_lyr[i][j], &soiltbl->sand_lyr[i][j],
                &soiltbl->iom_lyr[i][j], &soiltbl->bd_lyr[i][j],
                &soiltbl->NO3_lyr[i][j], &soiltbl->NH4_lyr[i][j]);

            if (match != 7 || j != layer - 1)
            {
                printf
                    ("Cannot read information of the %dth layer of the %dth"
                    "soil type!\n", j + 1, i + 1);
                printf (".soilinit file format error!\n");
                exit (1);
            }
        }
    }

    fclose (soil_file);
}

void ReadCrop (char *filename, croptbl_struct *croptbl)
{
    FILE           *crop_file;
    char            cmdstr[MAXSTRING];
    char            temp[MAXSTRING];
    int             j;

    crop_file = fopen (filename, "r");

    if (NULL == crop_file)
    {
        fprintf (stderr, "Error opening %s.\n", filename);
        PIHMError (1);
    }

    if (verbose_mode)
    {
        printf ("Reading %s...\n", filename);
    }

    /* Read crop description file */
    /* First count how many crop types are there in the description file */
    croptbl->number = CountOccurance (crop_file, "NAME");

    croptbl->cropName = (char **)malloc (croptbl->number * sizeof (char *));
    croptbl->userFloweringTT =
        (double *)malloc (croptbl->number * sizeof (int));
    croptbl->userMaturityTT =
        (double *)malloc (croptbl->number * sizeof (int));
    croptbl->userMaximumSoilCoverage =
        (double *)malloc (croptbl->number * sizeof (double));

    croptbl->userMaximumRootingDepth =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userExpectedYieldAvg =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userExpectedYieldMax =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userExpectedYieldMin =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userPercentMoistureInYield =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userFractionResidueStanding =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userFractionResidueRemoved =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userClippingBiomassThresholdUpper =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userClippingBiomassThresholdLower =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userClippingTiming =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userClippingDestiny =
        (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userTranspirationMinTemperature =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userTranspirationThresholdTemperature =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userColdDamageMinTemperature =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userColdDamageThresholdTemperature =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userTemperatureBase =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userTemperatureOptimum =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userTemperatureMaximum =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userShootPartitionInitial =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userShootPartitionFinal =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userRadiationUseEfficiency =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userTranspirationUseEfficiency =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userHIx = (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userHIo = (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userHIk = (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userEmergenceTT =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userNMaxConcentration =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userNDilutionSlope =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userKc = (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userAnnual = (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userLegume = (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userC3orC4 = (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userExtinctionCoefficient =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userPlantingDensity =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userClippingStart =
        (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userClippingEnd = (int *)malloc (croptbl->number * sizeof (int));
    croptbl->LWP_StressOnset =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->LWP_WiltingPoint =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->transpirationMax =
        (double *)malloc (croptbl->number * sizeof (double));

    /* Rewind to the beginning of file */
    FindLine (crop_file, "BOF");

    /* Read crop properties */
    for (j = 0; j < croptbl->number; j++)
    {
        croptbl->cropName[j] = (char *)malloc (MAXSTRING * sizeof (char));
        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "NAME", croptbl->cropName[j], 's');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "FLOWERING_TT",
            &croptbl->userFloweringTT[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "MATURITY_TT",
            &croptbl->userMaturityTT[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "MAXIMUM_SOIL_COVERAGE",
            &croptbl->userMaximumSoilCoverage[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "MAXIMUM_ROOTING_DEPTH",
            &croptbl->userMaximumRootingDepth[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "AVERAGE_EXPECTED_YIELD",
            &croptbl->userExpectedYieldAvg[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "MAXIMUM_EXPECTED_YIELD",
            &croptbl->userExpectedYieldMax[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "MINIMUM_EXPECTED_YIELD",
            &croptbl->userExpectedYieldMin[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "COMMERCIAL_YIELD_MOISTURE",
            &croptbl->userPercentMoistureInYield[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "STANDING_RESIDUE_AT_HARVEST",
            &croptbl->userFractionResidueStanding[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "RESIDUE_REMOVED",
            &croptbl->userFractionResidueRemoved[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "CLIPPING_BIOMASS_THRESHOLD_UPPER",
            &croptbl->userClippingBiomassThresholdUpper[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "CLIPPING_BIOMASS_THRESHOLD_LOWER",
            &croptbl->userClippingBiomassThresholdLower[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "HARVEST_TIMING",
            &croptbl->userClippingTiming[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "CLIPPING_BIOMASS_DESTINY", temp, 's');
        if (strcasecmp ("REMOVE", temp) == 0)
        {
            croptbl->userClippingDestiny[j] = REMOVE_CLIPPING;
        }
        else if (strcasecmp ("RETURN", temp) == 0)
        {
            croptbl->userClippingDestiny[j] = RETURN_CLIPPING;
        }
        else if (strcasecmp ("GRAZING", temp) == 0)
        {
            croptbl->userClippingDestiny[j] = GRAZING_CLIPPING;
        }
        else
        {
            printf ("Option %s not recoganized!\n", temp);
            exit (1);
        }

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "MIN_TEMPERATURE_FOR_TRANSPIRATION",
            &croptbl->userTranspirationMinTemperature[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "THRESHOLD_TEMPERATURE_FOR_TRANPIRATION",
            &croptbl->userTranspirationThresholdTemperature[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "MIN_TEMPERATURE_FOR_COLD_DAMAGE",
            &croptbl->userColdDamageMinTemperature[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "THRESHOLD_TEMPERATURE_FOR_COLD_DAMAGE",
            &croptbl->userColdDamageThresholdTemperature[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "BASE_TEMPERATURE_FOR_DEVELOPMENT",
            &croptbl->userTemperatureBase[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "OPTIMUM_TEMPERATURE_FOR_DEVELOPEMENT",
            &croptbl->userTemperatureOptimum[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "MAX_TEMPERATURE_FOR_DEVELOPMENT",
            &croptbl->userTemperatureMaximum[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "INITIAL_PARTITIONING_TO_SHOOT",
            &croptbl->userShootPartitionInitial[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "FINAL_PARTITIONING_TO_SHOOT",
            &croptbl->userShootPartitionFinal[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "RADIATION_USE_EFFICIENCY",
            &croptbl->userRadiationUseEfficiency[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "TRANSPIRATION_USE_EFFICIENCY",
            &croptbl->userTranspirationUseEfficiency[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "MAXIMUM_HARVEST_INDEX",
            &croptbl->userHIx[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "MINIMUM_HARVEST_INDEX",
            &croptbl->userHIo[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "HARVEST_INDEX", &croptbl->userHIk[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "THERMAL_TIME_TO_EMERGENCE",
            &croptbl->userEmergenceTT[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "N_MAX_CONCENTRATION",
            &croptbl->userNMaxConcentration[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "N_DILUTION_SLOPE",
            &croptbl->userNDilutionSlope[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "KC", &croptbl->userKc[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "ANNUAL", &croptbl->userAnnual[j], 'i');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "LEGUME", &croptbl->userLegume[j], 'i');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "C3", &croptbl->userC3orC4[j], 'i');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "LWP_STRESS_ONSET",
            &croptbl->LWP_StressOnset[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "LWP_WILTING_POINT",
            &croptbl->LWP_WiltingPoint[j], 'd');

        NextLine (crop_file, cmdstr);
        ReadKeyword (cmdstr, "TRANSPIRATION_MAX",
            &croptbl->transpirationMax[j], 'd');
    }

    fclose (crop_file);
}

void ReadOperation (const agtbl_struct *agtbl, mgmttbl_struct *mgmttbl,
    const croptbl_struct *croptbl)
{
    FILE           *op_file;
    char            cmdstr[MAXSTRING];
    char            filename[MAXSTRING];
    int             ntill;
    int             nplnt;
    int             nirrg;
    int             nfert;
    int             nautoirrg;
    int             i, j, k;
    cropmgmt_struct *cropmgmt;
    op_struct      *q;

    mgmttbl->number = agtbl->nopfile;

    mgmttbl->cropmgmt =
        (cropmgmt_struct *)malloc (mgmttbl->number *
        sizeof (cropmgmt_struct));

    for (i = 0; i < mgmttbl->number; i++)
    {
        cropmgmt = &mgmttbl->cropmgmt[i];

        cropmgmt->usingAutoIrr = 0;

        sprintf (filename, "input/%s/%s", project, agtbl->opfilen[i]);
        op_file = fopen (filename, "r");

        if (NULL == op_file)
        {
            fprintf (stderr, "Error opening %s.\n", filename);
            PIHMError (1);
        }

        if (verbose_mode)
        {
            printf ("Reading %s...\n", filename);
        }

        FindLine (op_file, "BOF");
        nplnt = CountOccurance (op_file, "PLANTING");

        FindLine (op_file, "BOF");
        ntill = CountOccurance (op_file, "TILLAGE");

        FindLine (op_file, "BOF");
        nirrg = CountOccurance (op_file, "FIXED_IRRIGATION");

        FindLine (op_file, "BOF");
        nfert = CountOccurance (op_file, "FIXED_FERTILIZATION");

        FindLine (op_file, "BOF");
        nautoirrg = CountOccurance (op_file, "AUTO_IRRIGATION");

        cropmgmt->totalCropsPerRotation = nplnt;
        if (nplnt > 0)
        {
            cropmgmt->plantingOrder =
                (op_struct *)malloc (nplnt * sizeof (op_struct));

            /* Rewind to the beginning of file and read all planting operations */
            FindLine (op_file, "BOF");

            for (j = 0; j < nplnt; j++)
            {
                q = &(cropmgmt->plantingOrder[j]);

                FindLine (op_file, "PLANTING");

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "YEAR", &q->opYear, 'i');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "DOY", &q->opDay, 'i');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "CROP", q->cropName, 's');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "USE_AUTO_IRR",
                    &q->usesAutoIrrigation, 'i');
                if (q->usesAutoIrrigation == 0)
                {
                    q->usesAutoIrrigation = -1;
                }
                else
                {
                    cropmgmt->usingAutoIrr = 1;
                }

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "USE_AUTO_FERT",
                    &q->usesAutoFertilization, 'i');
                if (q->usesAutoFertilization == 0)
                {
                    q->usesAutoFertilization = -1;
                }

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "FRACTION", &q->plantingDensity, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "CLIPPING_START", &q->clippingStart, 'i');
                if (q->clippingStart > 366 || q->clippingStart < 1)
                {
                    printf
                        ("ERROR: Please specify valid DOY for clipping start date!\n");
                    exit (1);
                }

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "CLIPPING_END", &q->clippingEnd, 'i');
                if (q->clippingEnd > 366 || q->clippingEnd < 1)
                {
                    printf
                        ("ERROR: Please specify valid DOY for clipping end date!\n");
                    exit (1);
                }

                q->status = 0;

                /* Link planting order and crop description */
                for (k = 0; k < croptbl->number; k++)
                {
                    if (strcmp (q->cropName, croptbl->cropName[k]) == 0)
                    {
                        q->plantID = k;
                        break;
                    }
                }
                if (k >= croptbl->number)
                {
                    printf
                        ("ERROR: Cannot find the plant description of %s, please check your input file\n",
                        q->cropName);
                    exit (1);
                }
            }
        }

        cropmgmt->numTillage = ntill;
        if (ntill > 0)
        {
            cropmgmt->Tillage =
                (op_struct *)malloc (ntill * sizeof (op_struct));

            /* Rewind to the beginning of file and read all tillage operations */
            FindLine (op_file, "BOF");

            for (j = 0; j < ntill; j++)
            {
                q = &(cropmgmt->Tillage[j]);

                FindLine (op_file, "TILLAGE");

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "YEAR", &q->opYear, 'i');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "DOY", &q->opDay, 'i');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "TOOL", q->opToolName, 's');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "DEPTH", &q->opDepth, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "SOIL_DISTURB_RATIO", &q->opSDR, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "MIXING_EFFICIENCY",
                    &q->opMixingEfficiency, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "CROP_NAME", q->cropNameT, 's');

                /* Check if the specified crop exists */
                if (strcasecmp (q->cropNameT, "N/A") != 0 &&
                    strcasecmp (q->cropNameT, "All") != 0 &&
                    !CropExist (q->cropNameT, croptbl))
                {
                    printf ("ERROR: Crop name %s not recognized!\n",
                        q->cropNameT);
                    exit (1);
                }

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "FRAC_THERMAL_TIME",
                    &q->fractionThermalTime, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "KILL_EFFICIENCY",
                    &q->killEfficiency, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "GRAIN_HARVEST", &q->grainHarvest, 'i');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "FORAGE_HARVEST",
                    &q->forageHarvest, 'd');

                q->status = 0;
            }
        }

        cropmgmt->numFertilization = nfert;
        if (nfert > 0)
        {
            cropmgmt->FixedFertilization =
                (op_struct *)malloc (nfert * sizeof (op_struct));

            /* Rewind to the beginning of file and read all fertilization
             * operations */
            FindLine (op_file, "BOF");

            for (j = 0; j < nfert; j++)
            {
                q = &(cropmgmt->FixedFertilization[j]);

                FindLine (op_file, "FIXED_FERTILIZATION");

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "YEAR", &q->opYear, 'i');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "DOY", &q->opDay, 'i');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "SOURCE", q->opSource, 's');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "MASS", &q->opMass, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "FORM", q->opForm, 's');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "METHOD", q->opMethod, 's');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "LAYER", &q->opLayer, 'i');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "C_ORGANIC", &q->opC_Organic, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "C_CHARCOAL", &q->opC_Charcoal, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "N_ORGANIC", &q->opN_Organic, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "N_CHARCOAL", &q->opN_Charcoal, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "N_NH4", &q->opN_NH4, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "N_NO3", &q->opN_NO3, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "P_ORGANIC", &q->opP_Organic, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "P_CHARCOAL", &q->opP_Charcoal, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "P_INORGANIC", &q->opP_Inorganic, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "K", &q->opK, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "S", &q->opS, 'd');

                q->status = 0;

                if (q->opC_Organic + q->opC_Charcoal + q->opN_Organic +
                    q->opN_Charcoal + q->opN_NH4 + q->opN_NO3 +
                    q->opP_Organic + q->opP_Charcoal + q->opP_Inorganic +
                    q->opK + q->opS <= 1.0)
                {
                    q->opMass /= 1000.0;
                }
                else
                {
                    printf
                        ("ERROR: Added fertilization fractions must be <= 1\n");
                    exit (1);
                }
            }
        }

        cropmgmt->numIrrigation = nirrg;
        if (nirrg > 0)
        {
            cropmgmt->FixedIrrigation =
                (op_struct *)malloc (nirrg * sizeof (op_struct));

            /* Rewind to the beginning of file and read all irrigation
             * operations */
            FindLine (op_file, "BOF");

            for (j = 0; j < nirrg; j++)
            {
                q = &(cropmgmt->FixedIrrigation[j]);

                FindLine (op_file, "FIXED_IRRIGATION");

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "YEAR", &q->opYear, 'i');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "DOY", &q->opDay, 'i');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "VOLUME", &q->opVolume, 'd');

                q->status = 0;
            }
        }

        cropmgmt->numAutoIrrigation = nautoirrg;
        if (nautoirrg > 0)
        {
            cropmgmt->autoIrrigation =
                (autoirr_struct *)malloc (nautoirrg *
                sizeof (autoirr_struct));
            /* Rewind to the beginning of file and read all planting operations */
            FindLine (op_file, "BOF");

            for (j = 0; j < nautoirrg; j++)
            {
                FindLine (op_file, "AUTO_IRRIGATION");

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "CROP",
                    cropmgmt->autoIrrigation[j].cropName, 's');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "START_DAY",
                    &cropmgmt->autoIrrigation[j].startDay, 'i');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "STOP_DAY",
                    &cropmgmt->autoIrrigation[j].stopDay, 'i');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "WATER_DEPLETION",
                    &cropmgmt->autoIrrigation[j].waterDepletion, 'd');

                NextLine (op_file, cmdstr);
                ReadKeyword (cmdstr, "LAST_SOIL_LAYER",
                    &cropmgmt->autoIrrigation[j].lastSoilLayer, 'i');
            }
        }

        /* Link plating order and auto irrigation */
        for (j = 0; j < cropmgmt->totalCropsPerRotation; j++)
        {
            if (cropmgmt->usingAutoIrr &&
                cropmgmt->plantingOrder[j].usesAutoIrrigation == 1)
            {
                for (k = 0; k < nautoirrg; k++)
                {
                    if (strcmp (cropmgmt->plantingOrder[j].cropName,
                            cropmgmt->autoIrrigation[k].cropName) == 0)
                    {
                        cropmgmt->plantingOrder[j].usesAutoIrrigation = k;
                        break;
                    }
                }
                if (k >= nautoirrg)
                {
                    printf
                        ("ERROR: Cannot find the description of auto irrigation for %s!\n",
                        cropmgmt->plantingOrder[j].cropName);
                    exit (1);
                }
            }
            else
                cropmgmt->plantingOrder[j].usesAutoIrrigation = -1;
        }

        fclose (op_file);
    }
}

//
//void DailyCycles (CyclesStruct cycles, pihm_struct pihm, int t, char *project)
//{
//    int             y, d;
//    time_t          rawtime;
//    struct tm      *timestamp;
//
//    CropManagementStruct *cropmgmt;
//    CommunityStruct *Community;
//    ResidueStruct  *Residue;
//    SimControlStruct *SimControl;
//    SnowStruct     *Snow;
//    SoilStruct     *Soil;
//    SoilCarbonStruct *SoilCarbon;
//    WeatherStruct  *Weather;
//    SummaryStruct  *Summary;
//
//    FieldOperationStruct *plantingOrder;
//    FieldOperationStruct *FixedFertilization;
//    FieldOperationStruct *Tillage;
//    FieldOperationStruct *FixedIrrigation;
//    int             operation_index;
//    int             i, c;
//    int             kill_all = 0;
//
//    rawtime = (time_t)(t - 24 * 3600);
//    timestamp = gmtime (&rawtime);
//
//    y = timestamp->tm_year + 1900 - cycles->SimControl.simStartYear;
//
//    d = doy (timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday, 1);
//
//    printf ("Year %d doy %d\n", y, d);
//
//    for (i = 0; i < pihm->numele; i++)
//    {
//        kill_all = 0;
//
//        cropmgmt = &cycles->grid[i].cropmgmt;
//        Community = &cycles->grid[i].Community;
//        Residue = &cycles->grid[i].Residue;
//        SimControl = &cycles->SimControl;
//        Snow = &cycles->grid[i].Snow;
//        Soil = &cycles->grid[i].Soil;
//        SoilCarbon = &cycles->grid[i].SoilCarbon;
//        Weather = &cycles->grid[i].Weather;
//        Summary = &cycles->grid[i].Summary;
//
//        if (d == 1)
//        { 
//            FirstDOY (&cropmgmt->rotationYear, SimControl->yearsInRotation, Soil->totalLayers, SoilCarbon, Residue, Soil);
//        }
//
//        /* If any crop in the community is growing, run the growing crop subroutine */
//        if (Community->NumActiveCrop > 0)
//            GrowingCrop (cropmgmt->rotationYear, y, d, cropmgmt->ForcedHarvest, cropmgmt->numHarvest, Community, Residue, SimControl, Soil, SoilCarbon, Weather, Snow, project);
//
//        while (IsOperationToday (cropmgmt->rotationYear, d, cropmgmt->plantingOrder, cropmgmt->totalCropsPerRotation, &operation_index))
//        {
//            plantingOrder = &cropmgmt->plantingOrder[operation_index];
//            PlantingCrop (Community, cropmgmt, operation_index);
//            if (verbose_mode)
//                printf ("DOY %3.3d %-20s %s\n", d, "Planting", plantingOrder->cropName);
//        }
//        UpdateOperationStatus (cropmgmt->plantingOrder, cropmgmt->totalCropsPerRotation);
//
//        while (IsOperationToday (cropmgmt->rotationYear, d, cropmgmt->FixedFertilization, cropmgmt->numFertilization, &operation_index))
//        {
//            FixedFertilization = &cropmgmt->FixedFertilization[operation_index];
//            if (verbose_mode)
//                printf ("DOY %3.3d %-20s %s\n", d, "Fixed Fertilization", FixedFertilization->opSource);
//
//            ApplyFertilizer (FixedFertilization, Soil, Residue);
//        }
//        UpdateOperationStatus (cropmgmt->FixedFertilization, cropmgmt->numFertilization);
//
//        while (IsOperationToday (cropmgmt->rotationYear, d, cropmgmt->Tillage, cropmgmt->numTillage, &operation_index))
//        {
//            Tillage = &(cropmgmt->Tillage[operation_index]);
//            if (verbose_mode)
//                printf ("DOY %3.3d %-20s %s\n", d, "Tillage", Tillage->opToolName);
//
//            if (strcasecmp (Tillage->opToolName, "Kill_Crop") != 0)
//                ExecuteTillage (SoilCarbon->abgdBiomassInput, Tillage, cropmgmt->tillageFactor, Soil, Residue);
//            else if (Community->NumActiveCrop > 0)
//            {
//                if (strcasecmp (Tillage->cropNameT, "N/A") == 0 ||
//                    strcasecmp (Tillage->cropNameT, "All") == 0)
//                {
//                    kill_all = 1;
//                }
//
//                for (c = 0; c < Community->NumCrop; c++)
//                {
//                    if (Community->Crop[c].stageGrowth > NO_CROP)
//                    {
//                        if (kill_all || strcasecmp (Tillage->cropNameT, Community->Crop[c].cropName) == 0)
//                        {
//                            HarvestCrop (y, d, SimControl->simStartYear, &Community->Crop[c], Residue, Soil, SoilCarbon, Weather, project);
//                            Community->NumActiveCrop--;
//                        }
//                    }
//                }
//            }
//        }
//        UpdateOperationStatus (cropmgmt->Tillage, cropmgmt->numTillage);
//
//        UpdateCommunity (Community);
//
//        Soil->irrigationVol = 0.0;
//        while (IsOperationToday (cropmgmt->rotationYear, d, cropmgmt->FixedIrrigation, cropmgmt->numIrrigation, &operation_index))
//        {
//            FixedIrrigation = &(cropmgmt->FixedIrrigation[operation_index]);
//            if (verbose_mode)
//                printf ("DOY %3.3d %-20s %lf\n", d, "Irrigation", FixedIrrigation->opVolume);
//
//            Soil->irrigationVol += FixedIrrigation->opVolume;
//        }
//        UpdateOperationStatus (cropmgmt->FixedIrrigation, cropmgmt->numIrrigation);
//
//        ComputeResidueCover (Residue);
//
//        TillageFactorSettling (cropmgmt->tillageFactor, Soil->totalLayers, Soil->waterContent, Soil->Porosity);
//
//        SnowProcesses (Snow, y, d, Weather, Residue->stanResidueTau, Community->svRadiationInterception);
//
//        Redistribution (y, d, Weather->precipitation[y][d - 1], Snow->snowFall, Snow->snowMelt, SimControl->hourlyInfiltration, Community, Soil, Residue);
//
//        ResidueEvaporation (Residue, Soil, Community, Weather->ETref[y][d - 1], Snow->snowCover);
//
//        Evaporation (Soil, Community, Residue, Weather->ETref[y][d - 1], Snow->snowCover);
//
//        Temperature (y, d, Snow->snowCover, Community->svRadiationInterception, Soil, Weather, Residue);
//
//        ComputeFactorComposite (SoilCarbon, d, y, Weather->lastDoy[y], Soil);
//
//        ComputeSoilCarbonBalanceMB (SoilCarbon, y, Residue, Soil, cropmgmt->tillageFactor);
//
//        NitrogenTransformation (y, d, Soil, Community, Residue, Weather, SoilCarbon);
//
//        if (d == Weather->lastDoy[y])
//        {
//            LastDOY (y, SimControl->simStartYear, Soil->totalLayers, Soil, SoilCarbon, Residue, Summary, project);
//        }
//    }
//}

int CropExist (char *cropName, const croptbl_struct *croptbl)
{
    int             i;
    int             exist = 0;

    for (i = 0; i < croptbl->number; i++)
    {
        if (strcmp (cropName, croptbl->cropName[i]) == 0)
        {
            exist = 1;
            break;
        }
    }

    return (exist);
}
