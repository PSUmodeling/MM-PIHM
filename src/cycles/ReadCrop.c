#include "Cycles.h"

void ReadCrop (char *filename, CommunityStruct *Community)
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
    FILE           *crop_file;
    char            fullname[MAXSTRING];
    char            cmdstr[MAXSTRING];
    char            temp[MAXSTRING];
    int             crop_counter = 0;
    int             i;
    CropStruct     *Crop;

    /* Open simulation control file */
    sprintf (fullname, "input/%s", filename);
    crop_file = fopen (fullname, "r");

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
    Community->NumCrop = crop_counter;
    Community->Crop = (CropStruct *)malloc (crop_counter * sizeof (CropStruct));

    /* Rewind to the beginning of file */
    FindLine (crop_file, "BOF");

    /* Read crop properties */
    for (i = 0; i < Community->NumCrop; i++)
    {
        Crop = &(Community->Crop[i]);

        NextLine (crop_file, cmdstr);
        ReadKeywordStr (cmdstr, "NAME", Crop->cropName);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "AVERAGE_SEEDING_DATE", &Crop->userSeedingDate);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "AVERAGE_50%_FLOWERING_DATE", &Crop->userFloweringDate);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "AVERAGE_MATURITY_DATE", &Crop->userMaturityDate);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAXIMUM_SOIL_COVERAGE", &Crop->userMaximumSoilCoverage);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAXIMUM_ROOTING_DEPTH", &Crop->userMaximumRootingDepth);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "AVERAGE_EXPECTED_YIELD", &Crop->userExpectedYieldAvg);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAXIMUM_EXPECTED_YIELD", &Crop->userExpectedYieldMax);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MINIMUM_EXPECTED_YIELD", &Crop->userExpectedYieldMin);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "COMMERCIAL_YIELD_MOISTURE", &Crop->userPercentMoistureInYield);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "STANDING_RESIDUE_AT_HARVEST", &Crop->userFractionResidueStanding);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "RESIDUE_REMOVED", &Crop->userFractionResidueRemoved);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "CLIPPING_BIOMASS_THRESHOLD_UPPER", &Crop->userClippingBiomassThresholdUpper);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "CLIPPING_BIOMASS_THRESHOLD_LOWER", &Crop->userClippingBiomassThresholdLower);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "HARVEST_TIMING", &Crop->userClippingTiming);

        NextLine (crop_file, cmdstr);
        ReadKeywordStr (cmdstr, "CLIPPING_BIOMASS_DESTINY", temp);
        if (strcasecmp ("REMOVE", temp) == 0)
            Crop->userClippingDestiny = REMOVE_CLIPPING;
        else if (strcasecmp ("RETURN", temp) == 0)
            Crop->userClippingDestiny = RETURN_CLIPPING;
        else if (strcasecmp ("GRAZING", temp) == 0)
            Crop->userClippingDestiny = GRAZING_CLIPPING;
        else
        {
            printf ("Option %s not recoganized!\n", temp);
            exit (1);
        }

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MIN_TEMPERATURE_FOR_TRANSPIRATION", &Crop->userTranspirationMinTemperature);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "THRESHOLD_TEMPERATURE_FOR_TRANPIRATION", &Crop->userTranspirationThresholdTemperature);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MIN_TEMPERATURE_FOR_COLD_DAMAGE", &Crop->userColdDamageMinTemperature);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "THRESHOLD_TEMPERATURE_FOR_COLD_DAMAGE", &Crop->userColdDamageThresholdTemperature);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "BASE_TEMPERATURE_FOR_DEVELOPMENT", &Crop->userTemperatureBase);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "OPTIMUM_TEMPERATURE_FOR_DEVELOPEMENT", &Crop->userTemperatureOptimum);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAX_TEMPERATURE_FOR_DEVELOPMENT", &Crop->userTemperatureMaximum);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "INITIAL_PARTITIONING_TO_SHOOT", &Crop->userShootPartitionInitial);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "FINAL_PARTITIONING_TO_SHOOT", &Crop->userShootPartitionFinal);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "RADIATION_USE_EFFICIENCY", &Crop->userRadiationUseEfficiency);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "TRANSPIRATION_USE_EFFICIENCY", &Crop->userTranspirationUseEfficiency);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MAXIMUM_HARVEST_INDEX", &Crop->userHIx);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "MINIMUM_HARVEST_INDEX", &Crop->userHIo);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "HARVEST_INDEX", &Crop->userHIk);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "THERMAL_TIME_TO_EMERGENCE", &Crop->userEmergenceTT);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "N_MAX_CONCENTRATION", &Crop->userNMaxConcentration);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "N_DILUTION_SLOPE", &Crop->userNDilutionSlope);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "KC", &Crop->userKc);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "ANNUAL", &Crop->userAnnual);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "LEGUME", &Crop->userLegume);

        NextLine (crop_file, cmdstr);
        ReadKeywordInt (cmdstr, "C3", &Crop->userC3orC4);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "LWP_STRESS_ONSET", &Crop->LWP_StressOnset);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "LWP_WILTING_POINT", &Crop->LWP_WiltingPoint);

        NextLine (crop_file, cmdstr);
        ReadKeywordDouble (cmdstr, "TRANSPIRATION_MAX", &Crop->transpirationMax);

        /* Convert units */
        Crop->userMaximumSoilCoverage = Crop->userMaximumSoilCoverage * 0.94 / 100.0;
        Crop->userPercentMoistureInYield = Crop->userPercentMoistureInYield / 100.0;
        Crop->userFractionResidueStanding = Crop->userFractionResidueStanding / 100.0;
        Crop->userFractionResidueRemoved = Crop->userFractionResidueRemoved / 100.0;
        if (Crop->userClippingTiming != BADVAL)
            Crop->userClippingTiming = Crop->userClippingTiming / 100.0;
        else
            Crop->userClippingTiming = 0.0;

        Crop->calculatedFloweringTT = 0.0;
        Crop->calculatedMaturityTT = 0.0;
        Crop->calculatedSimAvgYield = 0.0;
        Crop->calculatedSimMaxYield = 0.0;
        Crop->calculatedSimMinYield = 0.0;
    }

    fclose (crop_file);
}
