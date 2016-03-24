#include "pihm.h"

//void ReadCycles (pihm_struct pihm)
//{
//    /*
//     * Read simulation control file
//     */
//    ReadSimControl (pihm->filename.cycles, &pihm->agtbl, pihm->numele);
//
//    /*
//     * Read soil initialization file
//     */
//    ReadSoilInit (pihm->filename.soilinit, &pihm->soiltbl);
//
//    /*
//     * Read crop description file
//     */
//    ReadCrop (pihm->filename.crop, &pihm->croptbl);
//
//    /*
//     * Read operation file
//     */
//    ReadOperation (pihm->filename.op, &pihm->mgmttbl, &pihm->croptbl);
//
//    ///* Copy operation to all model grids */
//    //for (i = 0; i < pihm->numele; i++)
//    //{
//    //    cycles->grid[i].CropManagement.plantingOrder = (FieldOperationStruct *)malloc (CropManagement.totalCropsPerRotation * sizeof (FieldOperationStruct));
//    //    cycles->grid[i].CropManagement.ForcedHarvest = (FieldOperationStruct *)malloc (CropManagement.numHarvest * sizeof (FieldOperationStruct));
//    //    cycles->grid[i].CropManagement.FixedFertilization = (FieldOperationStruct *)malloc (CropManagement.numFertilization * sizeof (FieldOperationStruct));
//    //    cycles->grid[i].CropManagement.FixedIrrigation = (FieldOperationStruct *)malloc (CropManagement.numIrrigation * sizeof (FieldOperationStruct));
//    //    cycles->grid[i].CropManagement.Tillage = (FieldOperationStruct *)malloc (CropManagement.numTillage * sizeof (FieldOperationStruct));
//    //    cycles->grid[i].CropManagement.autoIrrigation = (autoIrrigationStruct *)malloc (CropManagement.numAutoIrrigation * sizeof (autoIrrigationStruct));
//
//    //    cycles->grid[i].CropManagement = CropManagement;
//    //}
//
//    //ReadWeather (simulation, &weather, cycles->SimControl.simStartYear, cycles->SimControl.totalYears);
//
//    //for (i = 0; i < pihm->numele; i++)
//    //{
//    //    cycles->grid[i].Weather.wind = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
//    //    cycles->grid[i].Weather.ETref = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
//    //    cycles->grid[i].Weather.precipitation = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
//    //    cycles->grid[i].Weather.RHmax = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
//    //    cycles->grid[i].Weather.RHmin = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
//    //    cycles->grid[i].Weather.solarRadiation = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
//    //    cycles->grid[i].Weather.tMax = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
//    //    cycles->grid[i].Weather.tMin = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
//    //    cycles->grid[i].Weather.yearlyAmplitude = (double *)malloc (cycles->SimControl.totalYears * sizeof (double));
//    //    cycles->grid[i].Weather.annualAverageTemperature = (double *)malloc (cycles->SimControl.totalYears * sizeof (double));
//    //    cycles->grid[i].Weather.lastDoy = (int *)malloc (cycles->SimControl.totalYears * sizeof (int));
//    //    for (y = 0; y < cycles->SimControl.totalYears; y++)
//    //    {
//    //        cycles->grid[i].Weather.wind[y] = (double *)malloc (366 * sizeof (double));
//    //        cycles->grid[i].Weather.ETref[y] = (double *)malloc (366 * sizeof (double));
//    //        cycles->grid[i].Weather.precipitation[y] = (double *)malloc (366 * sizeof (double));
//    //        cycles->grid[i].Weather.RHmax[y] = (double *)malloc (366 * sizeof (double));
//    //        cycles->grid[i].Weather.RHmin[y] = (double *)malloc (366 * sizeof (double));
//    //        cycles->grid[i].Weather.solarRadiation[y] = (double *)malloc (366 * sizeof (double));
//    //        cycles->grid[i].Weather.tMax[y] = (double *)malloc (366 * sizeof (double));
//    //        cycles->grid[i].Weather.tMin[y] = (double *)malloc (366 * sizeof (double));
//    //    }
//
//    //    cycles->grid[i].Weather = weather;
//    //}
//
//}

//void CyclesInit (CyclesStruct cycles, pihm_struct pihm)
//{
//    int         i, c;
//
//    for (i = 0; i < pihm->numele; i++)
//    {
//        CalculateDerivedWeather (&cycles->grid[i].Weather, cycles->SimControl.totalYears);
//
//        InitializeSoil (&cycles->grid[i].Soil, &cycles->grid[i].Weather, &cycles->SimControl, &pihm->elem[i].soil);
//
//        InitializeResidue (&cycles->grid[i].Residue, cycles->SimControl.totalYears, cycles->grid[i].Soil.totalLayers);
//
//        InitializeSoilCarbon (&cycles->grid[i].SoilCarbon, cycles->grid[i].Soil.totalLayers);
//
//        cycles->grid[i].Community.NumActiveCrop = 0;
//
//        for (c = 0; c < cycles->grid[i].Community.NumCrop; c++)
//        {
//            cycles->grid[i].Community.Crop[c].stageGrowth = NO_CROP;
//        }
//
//        cycles->grid[i].CropManagement.tillageFactor = (double *)malloc (cycles->grid[i].Soil.totalLayers * sizeof (double));
//
//        cycles->grid[i].Snow.Snow = 0.0;
//
//        cycles->grid[i].Summary = (SummaryStruct) {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//
//        ComputeThermalTime (cycles->SimControl.totalYears, &cycles->grid[i].Community, &cycles->grid[i].Weather);
//    }
//
//}

void ReadCyclesCtrl (char *filename, agtbl_struct *agtbl, int numele)
{
    FILE           *simctrl_file;
    //time_t          rawtime;
    //struct tm      *timestamp;
    char            cmdstr[MAXSTRING];
    int             i;
    int             match;
    int             index;

    /* Open simulation control file */
    simctrl_file = fopen (filename, "r");
    CheckFile (simctrl_file, filename);

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
        match = sscanf (cmdstr, "%d %d %d %d %d %d", &index, &agtbl->op[i], &agtbl->rotsz[i], &agtbl->auto_N[i], &agtbl->auto_P[i], &agtbl->auto_S[i]);
        if (match != 6)
        {
            printf ("Cannot read information of the %dth element!\n", i + 1);
            printf (".att file format error!\n");
            exit (1);
        }
    }

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
    CheckFile (soil_file, filename);

    soiltbl->totalLayers = (int *) malloc (soiltbl->number * sizeof (int));
    soiltbl->clay_lyr = (double **) malloc (soiltbl->number * sizeof (double *));
    soiltbl->sand_lyr = (double **) malloc (soiltbl->number * sizeof (double *));
    soiltbl->iom_lyr = (double **) malloc (soiltbl->number * sizeof (double *));
    soiltbl->bd_lyr = (double **) malloc (soiltbl->number * sizeof (double *));
    soiltbl->NO3_lyr = (double **) malloc (soiltbl->number * sizeof (double *));
    soiltbl->NH4_lyr = (double **) malloc (soiltbl->number * sizeof (double *));

    /* Read soil file */
    FindLine (soil_file, "BOF");

    for (i = 0; i < soiltbl->number; i++)
    {
        NextLine (soil_file, cmdstr);
        ReadKeywordInt (cmdstr, "SOIL_TYPE", &index);

        if (i != index - 1)
        {
            printf ("Cannot read information of the %dth soil type!\n",
                i + 1);
            printf (".soilinit file format error!\n");
            exit (1);
        }

        NextLine (soil_file, cmdstr);
        ReadKeywordInt (cmdstr, "TOTAL_LAYERS", &soiltbl->totalLayers[i]);

        soiltbl->clay_lyr[i] =
            (double *) malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->sand_lyr[i] =
            (double *) malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->iom_lyr[i] =
            (double *) malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->bd_lyr[i] =
            (double *) malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->NO3_lyr[i] =
            (double *) malloc (soiltbl->totalLayers[i] * sizeof (double));
        soiltbl->NH4_lyr[i] =
            (double *) malloc (soiltbl->totalLayers[i] * sizeof (double));
        
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
                printf ("Cannot read information of the %dth layer of the %dth"
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
    CheckFile (crop_file, filename);

    /* Read crop description file */
    /* First count how many crop types are there in the description file */
    croptbl->number = CountOccurance (crop_file, "NAME");

    croptbl->cropName = (char **)malloc (croptbl->number * sizeof (char *));
    croptbl->userSeedingDate =
        (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userFloweringDate =
        (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userMaturityDate =
        (int *)malloc (croptbl->number * sizeof (int));
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
    croptbl->userHIx =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userHIo =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userHIk =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userEmergenceTT =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userNMaxConcentration =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userNDilutionSlope =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userKc =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userAnnual =
        (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userLegume =
        (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userC3orC4 =
        (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userExtinctionCoefficient =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userPlantingDensity =
        (double *)malloc (croptbl->number * sizeof (double));
    croptbl->userClippingStart =
        (int *)malloc (croptbl->number * sizeof (int));
    croptbl->userClippingEnd =
        (int *)malloc (croptbl->number * sizeof (int));
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
        ReadKeywordStr (cmdstr, "NAME", croptbl->cropName[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "AVERAGE_SEEDING_DATE", &croptbl->userSeedingDate[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "AVERAGE_50%_FLOWERING_DATE", &croptbl->userFloweringDate[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "AVERAGE_MATURITY_DATE", &croptbl->userMaturityDate[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAXIMUM_SOIL_COVERAGE", &croptbl->userMaximumSoilCoverage[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAXIMUM_ROOTING_DEPTH", &croptbl->userMaximumRootingDepth[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "AVERAGE_EXPECTED_YIELD", &croptbl->userExpectedYieldAvg[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAXIMUM_EXPECTED_YIELD", &croptbl->userExpectedYieldMax[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MINIMUM_EXPECTED_YIELD", &croptbl->userExpectedYieldMin[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "COMMERCIAL_YIELD_MOISTURE", &croptbl->userPercentMoistureInYield[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "STANDING_RESIDUE_AT_HARVEST", &croptbl->userFractionResidueStanding[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "RESIDUE_REMOVED", &croptbl->userFractionResidueRemoved[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "CLIPPING_BIOMASS_THRESHOLD_UPPER", &croptbl->userClippingBiomassThresholdUpper[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "CLIPPING_BIOMASS_THRESHOLD_LOWER", &croptbl->userClippingBiomassThresholdLower[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "HARVEST_TIMING", &croptbl->userClippingTiming[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordStr (cmdstr, "CLIPPING_BIOMASS_DESTINY", temp);
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
        ReadKeywordDouble (cmdstr, "MIN_TEMPERATURE_FOR_TRANSPIRATION", &croptbl->userTranspirationMinTemperature[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "THRESHOLD_TEMPERATURE_FOR_TRANPIRATION", &croptbl->userTranspirationThresholdTemperature[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MIN_TEMPERATURE_FOR_COLD_DAMAGE", &croptbl->userColdDamageMinTemperature[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "THRESHOLD_TEMPERATURE_FOR_COLD_DAMAGE", &croptbl->userColdDamageThresholdTemperature[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "BASE_TEMPERATURE_FOR_DEVELOPMENT", &croptbl->userTemperatureBase[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "OPTIMUM_TEMPERATURE_FOR_DEVELOPEMENT", &croptbl->userTemperatureOptimum[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAX_TEMPERATURE_FOR_DEVELOPMENT", &croptbl->userTemperatureMaximum[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "INITIAL_PARTITIONING_TO_SHOOT", &croptbl->userShootPartitionInitial[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "FINAL_PARTITIONING_TO_SHOOT", &croptbl->userShootPartitionFinal[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "RADIATION_USE_EFFICIENCY", &croptbl->userRadiationUseEfficiency[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "TRANSPIRATION_USE_EFFICIENCY", &croptbl->userTranspirationUseEfficiency[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAXIMUM_HARVEST_INDEX", &croptbl->userHIx[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MINIMUM_HARVEST_INDEX", &croptbl->userHIo[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "HARVEST_INDEX", &croptbl->userHIk[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "THERMAL_TIME_TO_EMERGENCE", &croptbl->userEmergenceTT[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "N_MAX_CONCENTRATION", &croptbl->userNMaxConcentration[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "N_DILUTION_SLOPE", &croptbl->userNDilutionSlope[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "KC", &croptbl->userKc[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "ANNUAL", &croptbl->userAnnual[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "LEGUME", &croptbl->userLegume[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "C3", &croptbl->userC3orC4[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "LWP_STRESS_ONSET", &croptbl->LWP_StressOnset[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "LWP_WILTING_POINT", &croptbl->LWP_WiltingPoint[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "TRANSPIRATION_MAX", &croptbl->transpirationMax[j]);
    }

    fclose (crop_file);
}

void ReadOperation (char *filename, mgmttbl_struct *mgmttbl, const croptbl_struct *croptbl)
{
    FILE           *op_file;
    char            cmdstr[MAXSTRING];
    int             ntill;
    int             nplnt;
    int             nirrg;
    int             nfert;
    int             nautoirrg;
    int             i, j;
    op_struct      *q;

    op_file = fopen (filename, "r");
    CheckFile (op_file, filename);

    FindLine (op_file, "BOF");

    mgmttbl->number = CountOccurance (op_file, "OPERATION_TS");

    FindLine (op_file, "BOF");
    if (mgmttbl->number > 0)
    {
        mgmttbl->cropmgmt =
            (cropmgmt_struct *)malloc (mgmttbl->number * sizeof (cropmgmt_struct));

        for (i = 0; i < mgmttbl->number; i++)
        {
            mgmttbl->cropmgmt[i].totalCropsPerRotation =
                CountOccuranceBefore (op_file, "PLANTING", "OPERATION_TS");
            mgmttbl->cropmgmt[i].numTillage =
                CountOccuranceBefore (op_file, "TILLAGE", "OPERATION_TS");
            mgmttbl->cropmgmt[i].numIrrigation =
                CountOccuranceBefore (op_file, "FIXED_IRRIGATION", "OPERATION_TS");
            mgmttbl->cropmgmt[i].numFertilization =
                CountOccuranceBefore (op_file, "FIXED_FERTILIZATION", "OPERATION_TS");
            mgmttbl->cropmgmt[i].numAutoIrrigation =
                CountOccuranceBefore (op_file, "AUTO_IRRIGATION", "OPERATION_TS");
        }
    }

//
//    /* Read field operation file and count numbers of operations */
//    FindLine (op_file, "BOF");
//    planting_counter = CountOccurance (op_file, "PLANTING");
//
//    FindLine (op_file, "BOF");
//    harvest_counter = CountOccurance (op_file, "FORCED_HARVEST");
//
//    FindLine (op_file, "BOF");
//    tillage_counter = CountOccurance (op_file, "TILLAGE");
//    
//    FindLine (op_file, "BOF");
//    irrigation_counter = CountOccurance (op_file, "FIXED_IRRIGATION");
//
//    FindLine (op_file, "BOF");
//    fertilization_counter = CountOccurance (op_file, "FIXED_FERTILIZATION");
//
//    FindLine (op_file, "BOF");
//    auto_irrigation_counter = CountOccurance (op_file, "AUTO_IRRIGATION");
//
//    /* Allocate memories for field operation classes */
//    CropManagement->totalCropsPerRotation = planting_counter;
//    CropManagement->plantingOrder = (FieldOperationStruct *)malloc (planting_counter * sizeof (FieldOperationStruct));
//
//    CropManagement->numHarvest = harvest_counter;
//    CropManagement->ForcedHarvest = (FieldOperationStruct *)malloc (harvest_counter * sizeof (FieldOperationStruct));
//
//    CropManagement->numFertilization = fertilization_counter;
//    CropManagement->FixedFertilization = (FieldOperationStruct *)malloc (fertilization_counter * sizeof (FieldOperationStruct));
//
//    CropManagement->numIrrigation = irrigation_counter;
//    CropManagement->FixedIrrigation = (FieldOperationStruct *)malloc (irrigation_counter * sizeof (FieldOperationStruct));
//
//    CropManagement->numTillage = tillage_counter;
//    CropManagement->Tillage = (FieldOperationStruct *)malloc (tillage_counter * sizeof (FieldOperationStruct));
//
//    CropManagement->numAutoIrrigation = auto_irrigation_counter;
//    CropManagement->autoIrrigation = (autoIrrigationStruct *)malloc (auto_irrigation_counter * sizeof (autoIrrigationStruct));
//    CropManagement->usingAutoIrr = 0;
//
//    if (planting_counter)
//    {
//        /* Rewind to the beginning of file and read all planting operations */
//        FindLine (op_file, "BOF");
//
//        for (i = 0; i < planting_counter; i++)
//        {
//            q = &(CropManagement->plantingOrder[i]);
//
//            FindLine (op_file, "PLANTING");
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "YEAR", &q->opYear);
//            if (q->opYear > yearsInRotation)
//            {
//                printf ("ERROR: Operation year is larger than years in rotation!\n");
//                printf ("Please remove this operation and retry.\n");
//            }
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "DOY", &q->opDay);
//            
//            NextLine (op_file, cmdstr);
//            ReadKeywordStr (cmdstr, "CROP", q->cropName);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "AUTO_IRRIGATION", &q->usesAutoIrrigation);
//            if (CropManagement->plantingOrder[i].usesAutoIrrigation == 0)
//                CropManagement->plantingOrder[i].usesAutoIrrigation = -1;
//            else
//                CropManagement->usingAutoIrr = 1;
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "AUTO_FERTILIZATION", &q->usesAutoFertilization);
//            if (CropManagement->plantingOrder[i].usesAutoFertilization == 0)
//                CropManagement->plantingOrder[i].usesAutoFertilization = -1;
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "FRACTION", &q->plantingDensity);
//
//            q->status = 0;
//
//            /* Link planting order and crop description */
//            for (j = 0; j < Community->NumCrop; j++)
//            {
//                if (strcmp (CropManagement->plantingOrder[i].cropName, Community->Crop[j].cropName) == 0)
//                {
//                    CropManagement->plantingOrder[i].plantID = j;
//                    break;
//                }
//            }
//            if (j >= Community->NumCrop)
//            {
//                printf ("ERROR: Cannot find the plant description of %s, please check your input file\n", CropManagement->plantingOrder[i].cropName);
//                exit (1);
//            }
//        }
//    }
//
//    if (harvest_counter)
//    {
//        /* Rewind to the beginning of file and read all forced harvest operations */
//        FindLine (op_file, "BOF");
//
//        for (i = 0; i < harvest_counter; i++)
//        {
//            q = &(CropManagement->ForcedHarvest[i]);
//
//            FindLine (op_file, "FORCED_HARVEST");
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "YEAR", &q->opYear);
//            if (q->opYear > yearsInRotation)
//            {
//                printf ("ERROR: Operation year is larger than years in rotation!\n");
//                printf ("Please remove this operation and retry.\n");
//            }
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "DOY", &q->opDay);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordStr (cmdstr, "CROP", q->cropName);
//
//            q->status = 0;
//
//            /* Link forced harvest and crop description */
//            for (j = 0; j < Community->NumCrop; j++)
//            {
//                if (strcmp (CropManagement->ForcedHarvest[i].cropName, Community->Crop[j].cropName) == 0)
//                {
//                    CropManagement->ForcedHarvest[i].plantID = j;
//                    break;
//                }
//            }
//            if (j >= Community->NumCrop)
//            {
//                printf ("ERROR: Cannot find the plant description of %s, please check your input file\n", CropManagement->ForcedHarvest[i].cropName);
//                exit (1);
//            }
//        }
//    }
//
//    if (tillage_counter)
//    {
//        /* Rewind to the beginning of file and read all tillage operations */
//        FindLine (op_file, "BOF");
//
//        for (i = 0; i < tillage_counter; i++)
//        {
//            q = &(CropManagement->Tillage[i]);
//
//            FindLine (op_file, "TILLAGE");
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "YEAR", &q->opYear);
//            if (q->opYear > yearsInRotation)
//            {
//                printf ("ERROR: Operation year is larger than years in rotation!\n");
//                printf ("Please remove this operation and retry.\n");
//            }
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "DOY", &q->opDay);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordStr (cmdstr, "TOOL", q->opToolName);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "DEPTH", &q->opDepth);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "SOIL_DISTURB_RATIO", &q->opSDR);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "MIXING_EFFICIENCY", &q->opMixingEfficiency);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordStr (cmdstr, "CROP_NAME", q->cropNameT);
//
//            /* Check if the specified crop exists */
//            if (strcasecmp (q->cropNameT, "N/A") != 0 &&
//                strcasecmp (q->cropNameT, "All") != 0 &&
//                !CropExist (q->cropNameT, Community))
//            {
//                printf ("ERROR: Crop name %s not recognized!\n", q->cropNameT);
//                exit (1);
//            }
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "FRAC_THERMAL_TIME", &q->fractionThermalTime);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "KILL_EFFICIENCY", &q->killEfficiency);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "GRAIN_HARVEST", &q->grainHarvest);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "FORAGE_HARVEST", &q->forageHarvest);
//
//            q->status = 0;
//        }
//    }
//
//    if (irrigation_counter)
//    {
//        /* Rewind to the beginning of file and read all irrigation
//         * operations */
//        FindLine (op_file, "BOF");
//
//        for (i = 0; i < irrigation_counter; i++)
//        {
//            q = &(CropManagement->FixedIrrigation[i]);
//
//            FindLine (op_file, "FIXED_IRRIGATION");
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "YEAR", &q->opYear);
//            if (q->opYear > yearsInRotation)
//            {
//                printf ("ERROR: Operation year is larger than years in rotation!\n");
//                printf ("Please remove this operation and retry.\n");
//            }
//            
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "DOY", &q->opDay);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "VOLUME", &q->opVolume);
//
//            q->status = 0;
//        }
//    }
//
//    if (fertilization_counter)
//    {
//        /* Rewind to the beginning of file and read all fertilization
//         * operations */
//        FindLine (op_file, "BOF");
//
//        for (i = 0; i < fertilization_counter; i++)
//        {
//            q = &(CropManagement->FixedFertilization[i]);
//
//            FindLine (op_file, "FIXED_FERTILIZATION");
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "YEAR", &q->opYear);
//            if (q->opYear > yearsInRotation)
//            {
//                printf ("ERROR: Operation year is larger than years in rotation!\n");
//                printf ("Please remove this operation and retry.\n");
//            }
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "DOY", &q->opDay);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordStr (cmdstr, "SOURCE", q->opSource);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "MASS", &q->opMass);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordStr (cmdstr, "FORM", q->opForm);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordStr (cmdstr, "METHOD", q->opMethod);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "LAYER", &q->opLayer);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "C_ORGANIC", &q->opC_Organic);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "C_CHARCOAL", &q->opC_Charcoal);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "N_ORGANIC", &q->opN_Organic);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "N_CHARCOAL", &q->opN_Charcoal);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "N_NH4", &q->opN_NH4);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "N_NO3", &q->opN_NO3);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "P_ORGANIC", &q->opP_Organic);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "P_CHARCOAL", &q->opP_Charcoal);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "P_INORGANIC", &q->opP_Inorganic);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "K", &q->opK);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "S", &q->opS);
//
//            q->status = 0;
//
//            if (q->opC_Organic + q->opC_Charcoal + q->opN_Organic + q->opN_Charcoal + q->opN_NH4 + q->opN_NO3 + q->opP_Organic + q->opP_Charcoal + q->opP_Inorganic + q->opK + q->opS <= 1.0)
//            {
//                q->opMass /= 1000.0;
//            }
//            else
//            {
//                printf ("ERROR: Added fertilization fractions must be <= 1\n");
//                exit (1);
//            }
//        }
//    }
//
//    if (CropManagement->usingAutoIrr)
//    {
//        /* Rewind to the beginning of file and read all planting operations */
//        FindLine (op_file, "BOF");
//
//        for (i = 0; i < auto_irrigation_counter; i++)
//        {
//            FindLine (op_file, "AUTO_IRRIGATION");
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordStr (cmdstr, "CROP", CropManagement->autoIrrigation[i].cropName);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "START_DAY", &CropManagement->autoIrrigation[i].startDay);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "STOP_DAY", &CropManagement->autoIrrigation[i].stopDay);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordDouble (cmdstr, "WATER_DEPLETION", &CropManagement->autoIrrigation[i].waterDepletion);
//
//            NextLine (op_file, cmdstr);
//            ReadKeywordInt (cmdstr, "LAST_SOIL_LAYER", &CropManagement->autoIrrigation[i].lastSoilLayer);
//        }
//    }
//
//    /* Link plating order and auto irrigation */
//    for (i = 0; i < CropManagement->totalCropsPerRotation; i++)
//    {
//        if (CropManagement->usingAutoIrr && CropManagement->plantingOrder[i].usesAutoIrrigation == 1)
//        {
//            for (j = 0; j < auto_irrigation_counter; j++)
//            {
//                if (strcmp (CropManagement->plantingOrder[i].cropName, CropManagement->autoIrrigation[j].cropName) == 0)
//                {
//                    CropManagement->plantingOrder[i].usesAutoIrrigation = j;
//                    break;
//                }
//            }
//            if (j >= auto_irrigation_counter)
//            {
//                printf ("ERROR: Cannot find the description of auto irrigation for %s!\n", CropManagement->plantingOrder[i].cropName);
//                exit (1);
//            }
//        }
//        else
//            CropManagement->plantingOrder[i].usesAutoIrrigation = -1;
//    }
//
    fclose (op_file);
}
//
//int CropExist (char *cropName, const CommunityStruct *Community)
//{
//    int             i;
//    int             exist = 0;
//
//
//    for (i = 0; i < Community->NumCrop; i++)
//    {
//        if (strcmp (cropName, Community->Crop[i].cropName) == 0)
//        {
//            exist = 1;
//            break;
//        }
//    }
//
//    return (exist);
//}
//
//
//void DailyCycles (CyclesStruct cycles, pihm_struct pihm, int t, char *project)
//{
//    int             y, d;
//    time_t          rawtime;
//    struct tm      *timestamp;
//
//    CropManagementStruct *CropManagement;
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
//        CropManagement = &cycles->grid[i].CropManagement;
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
//            FirstDOY (&CropManagement->rotationYear, SimControl->yearsInRotation, Soil->totalLayers, SoilCarbon, Residue, Soil);
//        }
//
//        /* If any crop in the community is growing, run the growing crop subroutine */
//        if (Community->NumActiveCrop > 0)
//            GrowingCrop (CropManagement->rotationYear, y, d, CropManagement->ForcedHarvest, CropManagement->numHarvest, Community, Residue, SimControl, Soil, SoilCarbon, Weather, Snow, project);
//
//        while (IsOperationToday (CropManagement->rotationYear, d, CropManagement->plantingOrder, CropManagement->totalCropsPerRotation, &operation_index))
//        {
//            plantingOrder = &CropManagement->plantingOrder[operation_index];
//            PlantingCrop (Community, CropManagement, operation_index);
//            if (verbose_mode)
//                printf ("DOY %3.3d %-20s %s\n", d, "Planting", plantingOrder->cropName);
//        }
//        UpdateOperationStatus (CropManagement->plantingOrder, CropManagement->totalCropsPerRotation);
//
//        while (IsOperationToday (CropManagement->rotationYear, d, CropManagement->FixedFertilization, CropManagement->numFertilization, &operation_index))
//        {
//            FixedFertilization = &CropManagement->FixedFertilization[operation_index];
//            if (verbose_mode)
//                printf ("DOY %3.3d %-20s %s\n", d, "Fixed Fertilization", FixedFertilization->opSource);
//
//            ApplyFertilizer (FixedFertilization, Soil, Residue);
//        }
//        UpdateOperationStatus (CropManagement->FixedFertilization, CropManagement->numFertilization);
//
//        while (IsOperationToday (CropManagement->rotationYear, d, CropManagement->Tillage, CropManagement->numTillage, &operation_index))
//        {
//            Tillage = &(CropManagement->Tillage[operation_index]);
//            if (verbose_mode)
//                printf ("DOY %3.3d %-20s %s\n", d, "Tillage", Tillage->opToolName);
//
//            if (strcasecmp (Tillage->opToolName, "Kill_Crop") != 0)
//                ExecuteTillage (SoilCarbon->abgdBiomassInput, Tillage, CropManagement->tillageFactor, Soil, Residue);
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
//        UpdateOperationStatus (CropManagement->Tillage, CropManagement->numTillage);
//
//        UpdateCommunity (Community);
//
//        Soil->irrigationVol = 0.0;
//        while (IsOperationToday (CropManagement->rotationYear, d, CropManagement->FixedIrrigation, CropManagement->numIrrigation, &operation_index))
//        {
//            FixedIrrigation = &(CropManagement->FixedIrrigation[operation_index]);
//            if (verbose_mode)
//                printf ("DOY %3.3d %-20s %lf\n", d, "Irrigation", FixedIrrigation->opVolume);
//
//            Soil->irrigationVol += FixedIrrigation->opVolume;
//        }
//        UpdateOperationStatus (CropManagement->FixedIrrigation, CropManagement->numIrrigation);
//
//        ComputeResidueCover (Residue);
//
//        TillageFactorSettling (CropManagement->tillageFactor, Soil->totalLayers, Soil->waterContent, Soil->Porosity);
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
//        ComputeSoilCarbonBalanceMB (SoilCarbon, y, Residue, Soil, CropManagement->tillageFactor);
//
//        NitrogenTransformation (y, d, Soil, Community, Residue, Weather, SoilCarbon);
//
//        if (d == Weather->lastDoy[y])
//        {
//            LastDOY (y, SimControl->simStartYear, Soil->totalLayers, Soil, SoilCarbon, Residue, Summary, project);
//        }
//    }
//}
