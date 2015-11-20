#include "Cycles.h"

void CyclesRead (char *simulation, CyclesStruct cycles, pihm_struct pihm)
{
    int             i, j;
    int             y;
    SoilStruct      soil;
    CommunityStruct community;
    CropManagementStruct CropManagement;
    WeatherStruct   weather;

    cycles->grid = (grid_struct *)malloc (pihm->numele * sizeof (grid_struct));

    /*
     * Read simulation control file
     */
    ReadSimControl (simulation, &cycles->SimControl, pihm);

    /*
     * Read soil initialization file
     */
    ReadSoilInit (simulation, &soil);

    /* Copy soil initialization to all model grids */
    for (i = 0; i < pihm->numele; i++)
    {
        cycles->grid[i].Soil.totalLayers = soil.totalLayers;

        cycles->grid[i].Soil.layerThickness = (double *)malloc (soil.totalLayers * sizeof (double));
        cycles->grid[i].Soil.Clay = (double *)malloc (soil.totalLayers * sizeof (double));
        cycles->grid[i].Soil.Sand = (double *)malloc (soil.totalLayers * sizeof (double));
        cycles->grid[i].Soil.IOM = (double *)malloc (soil.totalLayers * sizeof (double));
        cycles->grid[i].Soil.BD = (double *)malloc (soil.totalLayers * sizeof (double));
        cycles->grid[i].Soil.FC = (double *)malloc (soil.totalLayers * sizeof (double));
        cycles->grid[i].Soil.PWP = (double *)malloc (soil.totalLayers * sizeof (double));
        cycles->grid[i].Soil.NO3 = (double *)malloc (soil.totalLayers * sizeof (double));
        cycles->grid[i].Soil.NH4 = (double *)malloc (soil.totalLayers * sizeof (double));

        for (j = 0; j < soil.totalLayers; j++)
        {
            cycles->grid[i].Soil.layerThickness[j] = soil.layerThickness[j];
            cycles->grid[i].Soil.Clay[j] = soil.Clay[j];
            cycles->grid[i].Soil.Sand[j] = soil.Sand[j];
            cycles->grid[i].Soil.IOM[j] = soil.IOM[j];
            cycles->grid[i].Soil.BD[j] = soil.BD[j];
            cycles->grid[i].Soil.FC[j] = soil.FC[j];
            cycles->grid[i].Soil.PWP[j] = soil.PWP[j];
            cycles->grid[i].Soil.NO3[j] = soil.NO3[j];
            cycles->grid[i].Soil.NH4[j] = soil.NH4[j];
        }
    }

    /*
     * Read crop description file
     */

    ReadCrop (simulation, &community);

    /* Copy crop description to all model grids */
    for (i = 0; i < pihm->numele; i++)
    {
        cycles->grid[i].Community.Crop = (CropStruct *)malloc (community.NumCrop * sizeof (CropStruct));

        cycles->grid[i].Community = community;
    }

    /*
     * Read operation file
     */

    ReadOperation (simulation, &CropManagement, &community, cycles->SimControl.yearsInRotation);

    /* Copy operation to all model grids */
    for (i = 0; i < pihm->numele; i++)
    {
        cycles->grid[i].CropManagement.plantingOrder = (FieldOperationStruct *)malloc (CropManagement.totalCropsPerRotation * sizeof (FieldOperationStruct));
        cycles->grid[i].CropManagement.ForcedHarvest = (FieldOperationStruct *)malloc (CropManagement.numHarvest * sizeof (FieldOperationStruct));
        cycles->grid[i].CropManagement.FixedFertilization = (FieldOperationStruct *)malloc (CropManagement.numFertilization * sizeof (FieldOperationStruct));
        cycles->grid[i].CropManagement.FixedIrrigation = (FieldOperationStruct *)malloc (CropManagement.numIrrigation * sizeof (FieldOperationStruct));
        cycles->grid[i].CropManagement.Tillage = (FieldOperationStruct *)malloc (CropManagement.numTillage * sizeof (FieldOperationStruct));
        cycles->grid[i].CropManagement.autoIrrigation = (autoIrrigationStruct *)malloc (CropManagement.numAutoIrrigation * sizeof (autoIrrigationStruct));

        cycles->grid[i].CropManagement = CropManagement;
    }

    ReadWeather (simulation, &weather, cycles->SimControl.simStartYear, cycles->SimControl.totalYears);

    for (i = 0; i < pihm->numele; i++)
    {
        cycles->grid[i].Weather.wind = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
        cycles->grid[i].Weather.ETref = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
        cycles->grid[i].Weather.precipitation = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
        cycles->grid[i].Weather.RHmax = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
        cycles->grid[i].Weather.RHmin = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
        cycles->grid[i].Weather.solarRadiation = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
        cycles->grid[i].Weather.tMax = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
        cycles->grid[i].Weather.tMin = (double **)malloc (cycles->SimControl.totalYears * sizeof (double *));
        cycles->grid[i].Weather.yearlyAmplitude = (double *)malloc (cycles->SimControl.totalYears * sizeof (double));
        cycles->grid[i].Weather.annualAverageTemperature = (double *)malloc (cycles->SimControl.totalYears * sizeof (double));
        cycles->grid[i].Weather.lastDoy = (int *)malloc (cycles->SimControl.totalYears * sizeof (int));
        for (y = 0; y < cycles->SimControl.totalYears; y++)
        {
            cycles->grid[i].Weather.wind[y] = (double *)malloc (366 * sizeof (double));
            cycles->grid[i].Weather.ETref[y] = (double *)malloc (366 * sizeof (double));
            cycles->grid[i].Weather.precipitation[y] = (double *)malloc (366 * sizeof (double));
            cycles->grid[i].Weather.RHmax[y] = (double *)malloc (366 * sizeof (double));
            cycles->grid[i].Weather.RHmin[y] = (double *)malloc (366 * sizeof (double));
            cycles->grid[i].Weather.solarRadiation[y] = (double *)malloc (366 * sizeof (double));
            cycles->grid[i].Weather.tMax[y] = (double *)malloc (366 * sizeof (double));
            cycles->grid[i].Weather.tMin[y] = (double *)malloc (366 * sizeof (double));
        }

        cycles->grid[i].Weather = weather;
    }

}

void CyclesInit (CyclesStruct cycles, pihm_struct pihm)
{
    int         i, c;

    for (i = 0; i < pihm->numele; i++)
    {
        CalculateDerivedWeather (&cycles->grid[i].Weather, cycles->SimControl.totalYears);

        InitializeSoil (&cycles->grid[i].Soil, &cycles->grid[i].Weather, &cycles->SimControl, &pihm->elem[i].soil);

        InitializeResidue (&cycles->grid[i].Residue, cycles->SimControl.totalYears, cycles->grid[i].Soil.totalLayers);

        InitializeSoilCarbon (&cycles->grid[i].SoilCarbon, cycles->grid[i].Soil.totalLayers);

        cycles->grid[i].Community.NumActiveCrop = 0;

        for (c = 0; c < cycles->grid[i].Community.NumCrop; c++)
        {
            cycles->grid[i].Community.Crop[c].stageGrowth = NO_CROP;
        }

        cycles->grid[i].CropManagement.tillageFactor = (double *)malloc (cycles->grid[i].Soil.totalLayers * sizeof (double));

        cycles->grid[i].Snow.Snow = 0.0;

        cycles->grid[i].Summary = (SummaryStruct) {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        ComputeThermalTime (cycles->SimControl.totalYears, &cycles->grid[i].Community, &cycles->grid[i].Weather);
    }

}

void ReadSimControl (char *simulation, SimControlStruct *SimControl, const pihm_struct pihm)
{
    /*
     * Read simulation control file
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * simctrl_file	    FILE*	File pointer of simulation control
     *					  file
     * filename		    char*	Simulation control file name
     */
    FILE           *simctrl_file;
    char           *token;
    char            tempname[MAXSTRING];
    char            project[MAXSTRING];
    char            filename[MAXSTRING];
    char            cmdstr[MAXSTRING];
    time_t          rawtime;
    struct tm      *timestamp;

    /* Detect if model is running in ensemble mode */
    strcpy (tempname, simulation);
    if (strstr (tempname, ".") != 0)
    {
        token = strtok (tempname, ".");
        strcpy (project, token);
    }
    else
    {
        strcpy (project, simulation);
    }

    /* Open simulation control file */
    sprintf (filename, "input/%s/%s.cycles", project, project);
    simctrl_file = fopen (filename, "r");
    printf ("%-30s %s.\n", "Read simulation control file:", filename);

    if (simctrl_file == NULL)
    {
        printf ("ERROR: Cannot find the simulation control file %s!\n", filename);
        exit (1);
    }

    /* Read simulation control file */
    FindLine (simctrl_file, "BOF");

    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "ROTATION_SIZE", &SimControl->yearsInRotation);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "ADJUSTED_YIELDS", &SimControl->adjustedYields);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "HOURLY_INFILTRATION", &SimControl->hourlyInfiltration);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "AUTOMATIC_NITROGEN", &SimControl->automaticNitrogen);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "AUTOMATIC_PHOSPHORUS", &SimControl->automaticPhosphorus);
    
    NextLine (simctrl_file, cmdstr);
    ReadKeywordInt (cmdstr, "AUTOMATIC_SULFUR", &SimControl->automaticSulfur);

    fclose (simctrl_file);

    rawtime = (int)pihm->ctrl.starttime;
    timestamp = gmtime (&rawtime);
    SimControl->simStartYear = timestamp->tm_year + 1900;

    rawtime = (int)pihm->ctrl.endtime;
    timestamp = gmtime (&rawtime);
    SimControl->simEndYear = timestamp->tm_year + 1900 - 1;

    SimControl->totalYears = SimControl->simEndYear - SimControl->simStartYear + 1;
}

void ReadSoilInit (char *simulation, SoilStruct *soil)
{
    /*
     * Read soil description file
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * soil_file	    FILE*	File pointer of soil description file
     * fullname		    char*	Full file name of the soil description
     *					  file
     * cmdstr		    char[MAXSTRING]
     *					Command string
     * i		    int		Loop counter
     */
    char           *token;
    char            tempname[MAXSTRING];
    char            project[MAXSTRING];
    char            filename[MAXSTRING];
    FILE           *soil_file;
    char            cmdstr[MAXSTRING];
    int             j;

    /* Detect if model is running in ensemble mode */
    strcpy (tempname, simulation);
    if (strstr (tempname, ".") != 0)
    {
        token = strtok (tempname, ".");
        strcpy (project, token);
    }
    else
    {
        strcpy (project, simulation);
    }

    /* Open simulation control file */
    sprintf (filename, "input/%s/%s.soilinit", project, simulation);
    soil_file = fopen (filename, "r");

    if (soil_file == NULL)
    {
        printf ("\nERROR: Cannot find the soil file %s!\n", filename);
        exit (1);
    }
    else
        printf ("%-30s input/%s.\n", "Read soil initialization file:", filename);

    /* Read soil file */
    FindLine (soil_file, "BOF");

    NextLine (soil_file, cmdstr);
    ReadKeywordDouble (cmdstr, "CURVE_NUMBER", &soil->Curve_Number);

    NextLine (soil_file, cmdstr);
    ReadKeywordDouble (cmdstr, "SLOPE", &soil->Percent_Slope);
    soil->Percent_Slope /= 100.0;

    NextLine (soil_file, cmdstr);
    ReadKeywordInt (cmdstr, "TOTAL_LAYERS", &soil->totalLayers);

    /* Allocate memories for soil class */
    soil->layerThickness = (double *)malloc (soil->totalLayers * sizeof (double));
    soil->Clay = (double *)malloc (soil->totalLayers * sizeof (double));
    soil->Sand = (double *)malloc (soil->totalLayers * sizeof (double));
    soil->IOM = (double *)malloc (soil->totalLayers * sizeof (double));
    soil->BD = (double *)malloc (soil->totalLayers * sizeof (double));
    soil->FC = (double *)malloc (soil->totalLayers * sizeof (double));
    soil->PWP = (double *)malloc (soil->totalLayers * sizeof (double));
    soil->NO3 = (double *)malloc (soil->totalLayers * sizeof (double));
    soil->NH4 = (double *)malloc (soil->totalLayers * sizeof (double));

    /* Skip header */
    NextLine (soil_file, cmdstr);

    for (j = 0; j < soil->totalLayers; j++)
    {
        NextLine (soil_file, cmdstr);
        sscanf (cmdstr, "%*d %lg %lf %lf %lf %lf %lf %lf %lf %lf", &soil->layerThickness[j], &soil->Clay[j], &soil->Sand[j], &soil->IOM[j], &soil->BD[j], &soil->FC[j], &soil->PWP[j], &soil->NO3[j], &soil->NH4[j]);
    }

    fclose (soil_file);
}

void ReadCrop (char *simulation, CommunityStruct *community)
{
    /*
     * Read crop description file
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * crop_file	    FILE*	File pointer of crop description file
     * fullname		    char*	Full file name of the crop description
     *					  file
     * cmdstr		    char[MAXSTRING]
     *					Command string
     * optstr		    char[MAXSTRING]
     *					Option argument string
     * crop_counter	    int		Crop counter
     * i		    int		Loop counter
     */
    char           *token;
    char            tempname[MAXSTRING];
    char            project[MAXSTRING];
    char            filename[MAXSTRING];
    FILE           *crop_file;
    char            cmdstr[MAXSTRING];
    char            temp[MAXSTRING];
    int             crop_counter = 0;
    int             j;
    CropStruct     *crop;

    /* Detect if model is running in ensemble mode */
    strcpy (tempname, simulation);
    if (strstr (tempname, ".") != 0)
    {
        token = strtok (tempname, ".");
        strcpy (project, token);
    }
    else
    {
        strcpy (project, simulation);
    }

    /* Open simulation control file */
    sprintf (filename, "input/%s/%s.crop", project, simulation);
    crop_file = fopen (filename, "r");

    if (crop_file == NULL)
    {
        printf ("\nError: Cannot find the crop description file %s!\n", filename);
        exit (1);
    }
    else
        printf ("%-30s input/%s.\n", "Read crop description file:", filename);

    /* Read crop description file */
    /* First count how many crop types are there in the description file */
    crop_counter = CountOccurance (crop_file, "NAME");

    /* Allocate memories for Crop classes */
    community->NumCrop = crop_counter;
    community->Crop = (CropStruct *)malloc (crop_counter * sizeof (CropStruct));

    /* Rewind to the beginning of file */
    FindLine (crop_file, "BOF");

    /* Read crop properties */
    for (j = 0; j < community->NumCrop; j++)
    {
        crop = &(community->Crop[j]);

        NextLine (crop_file, cmdstr);
        ReadKeywordStr (cmdstr, "NAME", crop->cropName);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "AVERAGE_SEEDING_DATE", &crop->userSeedingDate);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "AVERAGE_50%_FLOWERING_DATE", &crop->userFloweringDate);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "AVERAGE_MATURITY_DATE", &crop->userMaturityDate);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAXIMUM_SOIL_COVERAGE", &crop->userMaximumSoilCoverage);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAXIMUM_ROOTING_DEPTH", &crop->userMaximumRootingDepth);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "AVERAGE_EXPECTED_YIELD", &crop->userExpectedYieldAvg);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAXIMUM_EXPECTED_YIELD", &crop->userExpectedYieldMax);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MINIMUM_EXPECTED_YIELD", &crop->userExpectedYieldMin);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "COMMERCIAL_YIELD_MOISTURE", &crop->userPercentMoistureInYield);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "STANDING_RESIDUE_AT_HARVEST", &crop->userFractionResidueStanding);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "RESIDUE_REMOVED", &crop->userFractionResidueRemoved);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "CLIPPING_BIOMASS_THRESHOLD_UPPER", &crop->userClippingBiomassThresholdUpper);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "CLIPPING_BIOMASS_THRESHOLD_LOWER", &crop->userClippingBiomassThresholdLower);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "HARVEST_TIMING", &crop->userClippingTiming);

        NextLine (crop_file, cmdstr);
        ReadKeywordStr (cmdstr, "CLIPPING_BIOMASS_DESTINY", temp);
        if (strcasecmp ("REMOVE", temp) == 0)
            crop->userClippingDestiny = REMOVE_CLIPPING;
        else if (strcasecmp ("RETURN", temp) == 0)
            crop->userClippingDestiny = RETURN_CLIPPING;
        else if (strcasecmp ("GRAZING", temp) == 0)
            crop->userClippingDestiny = GRAZING_CLIPPING;
        else
        {
            printf ("Option %s not recoganized!\n", temp);
            exit (1);
        }

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MIN_TEMPERATURE_FOR_TRANSPIRATION", &crop->userTranspirationMinTemperature);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "THRESHOLD_TEMPERATURE_FOR_TRANPIRATION", &crop->userTranspirationThresholdTemperature);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MIN_TEMPERATURE_FOR_COLD_DAMAGE", &crop->userColdDamageMinTemperature);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "THRESHOLD_TEMPERATURE_FOR_COLD_DAMAGE", &crop->userColdDamageThresholdTemperature);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "BASE_TEMPERATURE_FOR_DEVELOPMENT", &crop->userTemperatureBase);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "OPTIMUM_TEMPERATURE_FOR_DEVELOPEMENT", &crop->userTemperatureOptimum);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAX_TEMPERATURE_FOR_DEVELOPMENT", &crop->userTemperatureMaximum);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "INITIAL_PARTITIONING_TO_SHOOT", &crop->userShootPartitionInitial);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "FINAL_PARTITIONING_TO_SHOOT", &crop->userShootPartitionFinal);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "RADIATION_USE_EFFICIENCY", &crop->userRadiationUseEfficiency);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "TRANSPIRATION_USE_EFFICIENCY", &crop->userTranspirationUseEfficiency);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAXIMUM_HARVEST_INDEX", &crop->userHIx);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MINIMUM_HARVEST_INDEX", &crop->userHIo);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "HARVEST_INDEX", &crop->userHIk);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "THERMAL_TIME_TO_EMERGENCE", &crop->userEmergenceTT);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "N_MAX_CONCENTRATION", &crop->userNMaxConcentration);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "N_DILUTION_SLOPE", &crop->userNDilutionSlope);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "KC", &crop->userKc);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "ANNUAL", &crop->userAnnual);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "LEGUME", &crop->userLegume);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "C3", &crop->userC3orC4);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "LWP_STRESS_ONSET", &crop->LWP_StressOnset);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "LWP_WILTING_POINT", &crop->LWP_WiltingPoint);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "TRANSPIRATION_MAX", &crop->transpirationMax);

        /* Convert units */
        crop->userMaximumSoilCoverage = crop->userMaximumSoilCoverage * 0.94 / 100.0;
        crop->userPercentMoistureInYield = crop->userPercentMoistureInYield / 100.0;
        crop->userFractionResidueStanding = crop->userFractionResidueStanding / 100.0;
        crop->userFractionResidueRemoved = crop->userFractionResidueRemoved / 100.0;
        if (crop->userClippingTiming != BADVAL)
            crop->userClippingTiming = crop->userClippingTiming / 100.0;
        else
            crop->userClippingTiming = 0.0;

        crop->calculatedFloweringTT = 0.0;
        crop->calculatedMaturityTT = 0.0;
        crop->calculatedSimAvgYield = 0.0;
        crop->calculatedSimMaxYield = 0.0;
        crop->calculatedSimMinYield = 0.0;
    }

    //for (i = 0; i < numele; i++)
    //{
    //    cycles->grid[i].Community.Crop = (CropStruct *)malloc (crop_counter * sizeof (CropStruct));

    //    cycles->grid[i].Community.NumCrop = crop_counter;
    //    for (j = 0; j < crop_counter; j++)
    //    {
    //        cycles->grid[i].Community.Crop[j] = crop[j];
    //    }
    //}

    fclose (crop_file);
}

void ReadOperation (char *simulation, CropManagementStruct *CropManagement, const CommunityStruct *Community, int yearsInRotation)
{
    /*
     * Read field operation description file
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * operation_file	    FILE*	File pointer of fiedl operation file
     * fullname		    char*	Full file name of the field operation
     *					  file
     * cmdstr		    char[MAXSTRING]
     *					Command string
     * optstr		    char[MAXSTRING]
     *					Option argument string
     * tillage_counter	    int		Tillage operation counter
     * planting_counter	    int		Planting operatin counter
     * irrigation_counter   int		Fixed irrigation operation counter
     * fertilization_counter
     *			    int		Fixed fertilization operation counter
     * auto_irrigation_counter
     *			    int		Automatic irrigation operation counter
     * harvest_counter	    int		Forced harvest operation counter
     * tempyear		    int
     * i		    int		Loop counter
     * j		    int		Loop counter
     * q		    FieldOperationStruct*
     */
    FILE           *operation_file;
    char           *token;
    char            tempname[MAXSTRING];
    char            project[MAXSTRING];
    char            filename[MAXSTRING];
    char            cmdstr[MAXSTRING];
    int             tillage_counter = 0;
    int             planting_counter = 0;
    int             irrigation_counter = 0;
    int             fertilization_counter = 0;
    int             auto_irrigation_counter = 0;
    int             harvest_counter = 0;
    int             i, j;
    FieldOperationStruct *q;

    /* Detect if model is running in ensemble mode */
    strcpy (tempname, simulation);
    if (strstr (tempname, ".") != 0)
    {
        token = strtok (tempname, ".");
        strcpy (project, token);
    }
    else
    {
        strcpy (project, simulation);
    }

    /* Open simulation control file */
    sprintf (filename, "input/%s/%s.operation", project, project);
    operation_file = fopen (filename, "r");
    printf ("%-30s %s.\n", "Read simulation control file:", filename);

    if (operation_file == NULL)
    {
        printf ("ERROR: Cannot find the field operation file %s!\n", filename);
        exit (1);
    }
    else
        printf ("%-30s input/%s.\n", "Read field operation file:", filename);

    /* Read field operation file and count numbers of operations */
    FindLine (operation_file, "BOF");
    planting_counter = CountOccurance (operation_file, "PLANTING");

    FindLine (operation_file, "BOF");
    harvest_counter = CountOccurance (operation_file, "FORCED_HARVEST");

    FindLine (operation_file, "BOF");
    tillage_counter = CountOccurance (operation_file, "TILLAGE");
    
    FindLine (operation_file, "BOF");
    irrigation_counter = CountOccurance (operation_file, "FIXED_IRRIGATION");

    FindLine (operation_file, "BOF");
    fertilization_counter = CountOccurance (operation_file, "FIXED_FERTILIZATION");

    FindLine (operation_file, "BOF");
    auto_irrigation_counter = CountOccurance (operation_file, "AUTO_IRRIGATION");

    /* Allocate memories for field operation classes */
    CropManagement->totalCropsPerRotation = planting_counter;
    CropManagement->plantingOrder = (FieldOperationStruct *)malloc (planting_counter * sizeof (FieldOperationStruct));

    CropManagement->numHarvest = harvest_counter;
    CropManagement->ForcedHarvest = (FieldOperationStruct *)malloc (harvest_counter * sizeof (FieldOperationStruct));

    CropManagement->numFertilization = fertilization_counter;
    CropManagement->FixedFertilization = (FieldOperationStruct *)malloc (fertilization_counter * sizeof (FieldOperationStruct));

    CropManagement->numIrrigation = irrigation_counter;
    CropManagement->FixedIrrigation = (FieldOperationStruct *)malloc (irrigation_counter * sizeof (FieldOperationStruct));

    CropManagement->numTillage = tillage_counter;
    CropManagement->Tillage = (FieldOperationStruct *)malloc (tillage_counter * sizeof (FieldOperationStruct));

    CropManagement->numAutoIrrigation = auto_irrigation_counter;
    CropManagement->autoIrrigation = (autoIrrigationStruct *)malloc (auto_irrigation_counter * sizeof (autoIrrigationStruct));
    CropManagement->usingAutoIrr = 0;

    if (planting_counter)
    {
        /* Rewind to the beginning of file and read all planting operations */
        FindLine (operation_file, "BOF");

        for (i = 0; i < planting_counter; i++)
        {
            q = &(CropManagement->plantingOrder[i]);

            FindLine (operation_file, "PLANTING");

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "YEAR", &q->opYear);
            if (q->opYear > yearsInRotation)
            {
                printf ("ERROR: Operation year is larger than years in rotation!\n");
                printf ("Please remove this operation and retry.\n");
            }

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "DOY", &q->opDay);
            
            NextLine (operation_file, cmdstr);
            ReadKeywordStr (cmdstr, "CROP", q->cropName);

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "AUTO_IRRIGATION", &q->usesAutoIrrigation);
            if (CropManagement->plantingOrder[i].usesAutoIrrigation == 0)
                CropManagement->plantingOrder[i].usesAutoIrrigation = -1;
            else
                CropManagement->usingAutoIrr = 1;

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "AUTO_FERTILIZATION", &q->usesAutoFertilization);
            if (CropManagement->plantingOrder[i].usesAutoFertilization == 0)
                CropManagement->plantingOrder[i].usesAutoFertilization = -1;

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "FRACTION", &q->plantingDensity);

            q->status = 0;

            /* Link planting order and crop description */
            for (j = 0; j < Community->NumCrop; j++)
            {
                if (strcmp (CropManagement->plantingOrder[i].cropName, Community->Crop[j].cropName) == 0)
                {
                    CropManagement->plantingOrder[i].plantID = j;
                    break;
                }
            }
            if (j >= Community->NumCrop)
            {
                printf ("ERROR: Cannot find the plant description of %s, please check your input file\n", CropManagement->plantingOrder[i].cropName);
                exit (1);
            }
        }
    }

    if (harvest_counter)
    {
        /* Rewind to the beginning of file and read all forced harvest operations */
        FindLine (operation_file, "BOF");

        for (i = 0; i < harvest_counter; i++)
        {
            q = &(CropManagement->ForcedHarvest[i]);

            FindLine (operation_file, "FORCED_HARVEST");

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "YEAR", &q->opYear);
            if (q->opYear > yearsInRotation)
            {
                printf ("ERROR: Operation year is larger than years in rotation!\n");
                printf ("Please remove this operation and retry.\n");
            }

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "DOY", &q->opDay);

            NextLine (operation_file, cmdstr);
            ReadKeywordStr (cmdstr, "CROP", q->cropName);

            q->status = 0;

            /* Link forced harvest and crop description */
            for (j = 0; j < Community->NumCrop; j++)
            {
                if (strcmp (CropManagement->ForcedHarvest[i].cropName, Community->Crop[j].cropName) == 0)
                {
                    CropManagement->ForcedHarvest[i].plantID = j;
                    break;
                }
            }
            if (j >= Community->NumCrop)
            {
                printf ("ERROR: Cannot find the plant description of %s, please check your input file\n", CropManagement->ForcedHarvest[i].cropName);
                exit (1);
            }
        }
    }

    if (tillage_counter)
    {
        /* Rewind to the beginning of file and read all tillage operations */
        FindLine (operation_file, "BOF");

        for (i = 0; i < tillage_counter; i++)
        {
            q = &(CropManagement->Tillage[i]);

            FindLine (operation_file, "TILLAGE");

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "YEAR", &q->opYear);
            if (q->opYear > yearsInRotation)
            {
                printf ("ERROR: Operation year is larger than years in rotation!\n");
                printf ("Please remove this operation and retry.\n");
            }

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "DOY", &q->opDay);

            NextLine (operation_file, cmdstr);
            ReadKeywordStr (cmdstr, "TOOL", q->opToolName);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "DEPTH", &q->opDepth);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "SOIL_DISTURB_RATIO", &q->opSDR);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "MIXING_EFFICIENCY", &q->opMixingEfficiency);

            NextLine (operation_file, cmdstr);
            ReadKeywordStr (cmdstr, "CROP_NAME", q->cropNameT);

            /* Check if the specified crop exists */
            if (strcasecmp (q->cropNameT, "N/A") != 0 &&
                strcasecmp (q->cropNameT, "All") != 0 &&
                !CropExist (q->cropNameT, Community))
            {
                printf ("ERROR: Crop name %s not recognized!\n", q->cropNameT);
                exit (1);
            }

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "FRAC_THERMAL_TIME", &q->fractionThermalTime);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "KILL_EFFICIENCY", &q->killEfficiency);

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "GRAIN_HARVEST", &q->grainHarvest);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "FORAGE_HARVEST", &q->forageHarvest);

            q->status = 0;
        }
    }

    if (irrigation_counter)
    {
        /* Rewind to the beginning of file and read all irrigation
         * operations */
        FindLine (operation_file, "BOF");

        for (i = 0; i < irrigation_counter; i++)
        {
            q = &(CropManagement->FixedIrrigation[i]);

            FindLine (operation_file, "FIXED_IRRIGATION");

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "YEAR", &q->opYear);
            if (q->opYear > yearsInRotation)
            {
                printf ("ERROR: Operation year is larger than years in rotation!\n");
                printf ("Please remove this operation and retry.\n");
            }
            
            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "DOY", &q->opDay);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "VOLUME", &q->opVolume);

            q->status = 0;
        }
    }

    if (fertilization_counter)
    {
        /* Rewind to the beginning of file and read all fertilization
         * operations */
        FindLine (operation_file, "BOF");

        for (i = 0; i < fertilization_counter; i++)
        {
            q = &(CropManagement->FixedFertilization[i]);

            FindLine (operation_file, "FIXED_FERTILIZATION");

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "YEAR", &q->opYear);
            if (q->opYear > yearsInRotation)
            {
                printf ("ERROR: Operation year is larger than years in rotation!\n");
                printf ("Please remove this operation and retry.\n");
            }

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "DOY", &q->opDay);

            NextLine (operation_file, cmdstr);
            ReadKeywordStr (cmdstr, "SOURCE", q->opSource);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "MASS", &q->opMass);

            NextLine (operation_file, cmdstr);
            ReadKeywordStr (cmdstr, "FORM", q->opForm);

            NextLine (operation_file, cmdstr);
            ReadKeywordStr (cmdstr, "METHOD", q->opMethod);

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "LAYER", &q->opLayer);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "C_ORGANIC", &q->opC_Organic);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "C_CHARCOAL", &q->opC_Charcoal);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "N_ORGANIC", &q->opN_Organic);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "N_CHARCOAL", &q->opN_Charcoal);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "N_NH4", &q->opN_NH4);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "N_NO3", &q->opN_NO3);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "P_ORGANIC", &q->opP_Organic);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "P_CHARCOAL", &q->opP_Charcoal);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "P_INORGANIC", &q->opP_Inorganic);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "K", &q->opK);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "S", &q->opS);

            q->status = 0;

            if (q->opC_Organic + q->opC_Charcoal + q->opN_Organic + q->opN_Charcoal + q->opN_NH4 + q->opN_NO3 + q->opP_Organic + q->opP_Charcoal + q->opP_Inorganic + q->opK + q->opS <= 1.0)
            {
                q->opMass /= 1000.0;
            }
            else
            {
                printf ("ERROR: Added fertilization fractions must be <= 1\n");
                exit (1);
            }
        }
    }

    if (CropManagement->usingAutoIrr)
    {
        /* Rewind to the beginning of file and read all planting operations */
        FindLine (operation_file, "BOF");

        for (i = 0; i < auto_irrigation_counter; i++)
        {
            FindLine (operation_file, "AUTO_IRRIGATION");

            NextLine (operation_file, cmdstr);
            ReadKeywordStr (cmdstr, "CROP", CropManagement->autoIrrigation[i].cropName);

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "START_DAY", &CropManagement->autoIrrigation[i].startDay);

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "STOP_DAY", &CropManagement->autoIrrigation[i].stopDay);

            NextLine (operation_file, cmdstr);
            ReadKeywordDouble (cmdstr, "WATER_DEPLETION", &CropManagement->autoIrrigation[i].waterDepletion);

            NextLine (operation_file, cmdstr);
            ReadKeywordInt (cmdstr, "LAST_SOIL_LAYER", &CropManagement->autoIrrigation[i].lastSoilLayer);
        }
    }

    /* Link plating order and auto irrigation */
    for (i = 0; i < CropManagement->totalCropsPerRotation; i++)
    {
        if (CropManagement->usingAutoIrr && CropManagement->plantingOrder[i].usesAutoIrrigation == 1)
        {
            for (j = 0; j < auto_irrigation_counter; j++)
            {
                if (strcmp (CropManagement->plantingOrder[i].cropName, CropManagement->autoIrrigation[j].cropName) == 0)
                {
                    CropManagement->plantingOrder[i].usesAutoIrrigation = j;
                    break;
                }
            }
            if (j >= auto_irrigation_counter)
            {
                printf ("ERROR: Cannot find the description of auto irrigation for %s!\n", CropManagement->plantingOrder[i].cropName);
                exit (1);
            }
        }
        else
            CropManagement->plantingOrder[i].usesAutoIrrigation = -1;
    }

    fclose (operation_file);
}

int CropExist (char *cropName, const CommunityStruct *Community)
{
    int             i;
    int             exist = 0;


    for (i = 0; i < Community->NumCrop; i++)
    {
        if (strcmp (cropName, Community->Crop[i].cropName) == 0)
        {
            exist = 1;
            break;
        }
    }

    return (exist);
}

void ReadWeather (char *simulation, WeatherStruct *Weather, int start_year, int total_years)
{
    /*
     * Read weather file
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * weather_file	    FILE*	File pointer of weather file
     * fullname		    char*	Full file name of the weather file
     * cmdstr		    char[MAXSTRING]
     *					Command string
     * optstr		    char[MAXSTRING]
     *					Option argument string
     * y		    int		year counter
     * doy		    int		Day of year counter
     * temp_year	    int
     * temp_doy		    int
     * start_year_str	    char[5]
     */
    FILE           *weather_file;
    int             y, doy;
    int             temp_year, temp_doy;
    char            start_year_str[5];
    char           *token;
    char            tempname[MAXSTRING];
    char            project[MAXSTRING];
    char            filename[MAXSTRING];
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];

    /* Detect if model is running in ensemble mode */
    strcpy (tempname, simulation);
    if (strstr (tempname, ".") != 0)
    {
        token = strtok (tempname, ".");
        strcpy (project, token);
    }
    else
    {
        strcpy (project, simulation);
    }

    /* Open simulation control file */
    sprintf (filename, "input/%s/%s.weather", project, project);
    weather_file = fopen (filename, "r");
    printf ("%-30s %s.\n", "Read weather file:", filename);

    if (weather_file == NULL)
    {
        printf ("ERROR: Cannot find the weather file %s!\n", filename);
        exit (1);
    }

    /* Allocate */
    Weather->wind = (double **)malloc (total_years * sizeof (double *));
    Weather->ETref = (double **)malloc (total_years * sizeof (double *));
    Weather->precipitation = (double **)malloc (total_years * sizeof (double *));
    Weather->RHmax = (double **)malloc (total_years * sizeof (double *));
    Weather->RHmin = (double **)malloc (total_years * sizeof (double *));
    Weather->solarRadiation = (double **)malloc (total_years * sizeof (double *));
    Weather->tMax = (double **)malloc (total_years * sizeof (double *));
    Weather->tMin = (double **)malloc (total_years * sizeof (double *));
    Weather->yearlyAmplitude = (double *)malloc (total_years * sizeof (double));
    Weather->annualAverageTemperature = (double *)malloc (total_years * sizeof (double));
    Weather->lastDoy = (int *)malloc (total_years * sizeof (int));
    for (y = 0; y < total_years; y++)
    {
        Weather->wind[y] = (double *)malloc (366 * sizeof (double));
        Weather->ETref[y] = (double *)malloc (366 * sizeof (double));
        Weather->precipitation[y] = (double *)malloc (366 * sizeof (double));
        Weather->RHmax[y] = (double *)malloc (366 * sizeof (double));
        Weather->RHmin[y] = (double *)malloc (366 * sizeof (double));
        Weather->solarRadiation[y] = (double *)malloc (366 * sizeof (double));
        Weather->tMax[y] = (double *)malloc (366 * sizeof (double));
        Weather->tMin[y] = (double *)malloc (366 * sizeof (double));
    }

    sprintf (start_year_str, "%4.4d", start_year);

    /* Read weather file */
    FindLine (weather_file, "BOF");

    /* Read in site information and count number of weather records */
    NextLine (weather_file, cmdstr);

    ReadKeywordDouble (cmdstr, "LATITUDE", &Weather->siteLatitude);

    NextLine (weather_file, cmdstr);
    ReadKeywordDouble (cmdstr, "ALTITUDE", &Weather->siteAltitude);

    NextLine (weather_file, cmdstr);
    ReadKeywordDouble (cmdstr, "SCREENING_HEIGHT", &Weather->screeningHeight);

    /* Skip header */
    NextLine (weather_file, cmdstr);

    strcpy (cmdstr, "\0");

    NextLine (weather_file, cmdstr);
    while (strcasecmp (cmdstr, "EOF") != 0)
    {
        sscanf (cmdstr, "%s", optstr);
        if (strcasecmp (start_year_str, optstr) == 0)
        {
            for (y = 0; y < total_years; y++)
            {
                for (doy = 1; doy < 367; doy++)
                {
                    sscanf (cmdstr, "%d %d %*f %*f %*f %*f %*f %*f %*f", &temp_year, &temp_doy);
                    if (temp_year == y + start_year && temp_doy == doy)
                    {
                        sscanf (cmdstr, "%*d %*d %lf %lf %lf %lf %lf %lf %lf", &Weather->precipitation[y][doy - 1], &Weather->tMax[y][doy - 1], &Weather->tMin[y][doy - 1], &Weather->solarRadiation[y][doy - 1], &Weather->RHmax[y][doy - 1], &Weather->RHmin[y][doy - 1], &Weather->wind[y][doy - 1]);
                        if (doy == 366)
                            Weather->lastDoy[y] = 366;
                        NextLine (weather_file, cmdstr);
                    }
                    else if (doy == 366 && temp_year == y + start_year + 1 && temp_doy == 1)
                        Weather->lastDoy[y] = 365;
                    else if (doy == 366 && feof (weather_file))
                        Weather->lastDoy[y] = 365;
                    else
                    {
                        printf ("ERROR: Please check your weather input file near YEAR: %4.4d, DOY: %-d, expecting %4.4d-%-d, eof status %d\n", temp_year, temp_doy, y, doy, feof (weather_file));
                        exit (1);
                    }
                }
            }
            break;
        }

        NextLine (weather_file, cmdstr);
    }

    fclose (weather_file);
}

void DailyCycles (CyclesStruct cycles, pihm_struct pihm, int t, char *project)
{
    int             y, d;
    time_t          rawtime;
    struct tm      *timestamp;

    CropManagementStruct *CropManagement;
    CommunityStruct *Community;
    ResidueStruct  *Residue;
    SimControlStruct *SimControl;
    SnowStruct     *Snow;
    SoilStruct     *Soil;
    SoilCarbonStruct *SoilCarbon;
    WeatherStruct  *Weather;
    SummaryStruct  *Summary;

    FieldOperationStruct *plantingOrder;
    FieldOperationStruct *FixedFertilization;
    FieldOperationStruct *Tillage;
    FieldOperationStruct *FixedIrrigation;
    int             operation_index;
    int             i, c;
    int             kill_all = 0;

    rawtime = (time_t)(t - 24 * 3600);
    timestamp = gmtime (&rawtime);

    y = timestamp->tm_year + 1900 - cycles->SimControl.simStartYear;

    d = doy (timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday, 1);

    printf ("Year %d doy %d\n", y, d);

    for (i = 0; i < pihm->numele; i++)
    {
        kill_all = 0;

        CropManagement = &cycles->grid[i].CropManagement;
        Community = &cycles->grid[i].Community;
        Residue = &cycles->grid[i].Residue;
        SimControl = &cycles->SimControl;
        Snow = &cycles->grid[i].Snow;
        Soil = &cycles->grid[i].Soil;
        SoilCarbon = &cycles->grid[i].SoilCarbon;
        Weather = &cycles->grid[i].Weather;
        Summary = &cycles->grid[i].Summary;

        if (d == 1)
        { 
            FirstDOY (&CropManagement->rotationYear, SimControl->yearsInRotation, Soil->totalLayers, SoilCarbon, Residue, Soil);
        }

        /* If any crop in the community is growing, run the growing crop subroutine */
        if (Community->NumActiveCrop > 0)
            GrowingCrop (CropManagement->rotationYear, y, d, CropManagement->ForcedHarvest, CropManagement->numHarvest, Community, Residue, SimControl, Soil, SoilCarbon, Weather, Snow, project);

        while (IsOperationToday (CropManagement->rotationYear, d, CropManagement->plantingOrder, CropManagement->totalCropsPerRotation, &operation_index))
        {
            plantingOrder = &CropManagement->plantingOrder[operation_index];
            PlantingCrop (Community, CropManagement, operation_index);
            if (verbose_mode)
                printf ("DOY %3.3d %-20s %s\n", d, "Planting", plantingOrder->cropName);
        }
        UpdateOperationStatus (CropManagement->plantingOrder, CropManagement->totalCropsPerRotation);

        while (IsOperationToday (CropManagement->rotationYear, d, CropManagement->FixedFertilization, CropManagement->numFertilization, &operation_index))
        {
            FixedFertilization = &CropManagement->FixedFertilization[operation_index];
            if (verbose_mode)
                printf ("DOY %3.3d %-20s %s\n", d, "Fixed Fertilization", FixedFertilization->opSource);

            ApplyFertilizer (FixedFertilization, Soil, Residue);
        }
        UpdateOperationStatus (CropManagement->FixedFertilization, CropManagement->numFertilization);

        while (IsOperationToday (CropManagement->rotationYear, d, CropManagement->Tillage, CropManagement->numTillage, &operation_index))
        {
            Tillage = &(CropManagement->Tillage[operation_index]);
            if (verbose_mode)
                printf ("DOY %3.3d %-20s %s\n", d, "Tillage", Tillage->opToolName);

            if (strcasecmp (Tillage->opToolName, "Kill_Crop") != 0)
                ExecuteTillage (SoilCarbon->abgdBiomassInput, Tillage, CropManagement->tillageFactor, Soil, Residue);
            else if (Community->NumActiveCrop > 0)
            {
                if (strcasecmp (Tillage->cropNameT, "N/A") == 0 ||
                    strcasecmp (Tillage->cropNameT, "All") == 0)
                {
                    kill_all = 1;
                }

                for (c = 0; c < Community->NumCrop; c++)
                {
                    if (Community->Crop[c].stageGrowth > NO_CROP)
                    {
                        if (kill_all || strcasecmp (Tillage->cropNameT, Community->Crop[c].cropName) == 0)
                        {
                            HarvestCrop (y, d, SimControl->simStartYear, &Community->Crop[c], Residue, Soil, SoilCarbon, Weather, project);
                            Community->NumActiveCrop--;
                        }
                    }
                }
            }
        }
        UpdateOperationStatus (CropManagement->Tillage, CropManagement->numTillage);

        UpdateCommunity (Community);

        Soil->irrigationVol = 0.0;
        while (IsOperationToday (CropManagement->rotationYear, d, CropManagement->FixedIrrigation, CropManagement->numIrrigation, &operation_index))
        {
            FixedIrrigation = &(CropManagement->FixedIrrigation[operation_index]);
            if (verbose_mode)
                printf ("DOY %3.3d %-20s %lf\n", d, "Irrigation", FixedIrrigation->opVolume);

            Soil->irrigationVol += FixedIrrigation->opVolume;
        }
        UpdateOperationStatus (CropManagement->FixedIrrigation, CropManagement->numIrrigation);

        ComputeResidueCover (Residue);

        TillageFactorSettling (CropManagement->tillageFactor, Soil->totalLayers, Soil->waterContent, Soil->Porosity);

        SnowProcesses (Snow, y, d, Weather, Residue->stanResidueTau, Community->svRadiationInterception);

        Redistribution (y, d, Weather->precipitation[y][d - 1], Snow->snowFall, Snow->snowMelt, SimControl->hourlyInfiltration, Community, Soil, Residue);

        ResidueEvaporation (Residue, Soil, Community, Weather->ETref[y][d - 1], Snow->snowCover);

        Evaporation (Soil, Community, Residue, Weather->ETref[y][d - 1], Snow->snowCover);

        Temperature (y, d, Snow->snowCover, Community->svRadiationInterception, Soil, Weather, Residue);

        ComputeFactorComposite (SoilCarbon, d, y, Weather->lastDoy[y], Soil);

        ComputeSoilCarbonBalanceMB (SoilCarbon, y, Residue, Soil, CropManagement->tillageFactor);

        NitrogenTransformation (y, d, Soil, Community, Residue, Weather, SoilCarbon);

        if (d == Weather->lastDoy[y])
        {
            LastDOY (y, SimControl->simStartYear, Soil->totalLayers, Soil, SoilCarbon, Residue, Summary, project);
        }
    }
}
