#include "Cycles.h"

void DailyOperations (int rotationYear, int y, int doy, CropManagementStruct *CropManagement, CommunityStruct *Community, ResidueStruct *Residue, SimControlStruct *SimControl, SnowStruct *Snow, SoilStruct *Soil, SoilCarbonStruct *SoilCarbon, WeatherStruct *Weather, const char *project)
{
    /*
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * FixedFertilization   FieldOperationStruct*
     * Tillage		    FieldOperationStruct*
     * FixedIrrigation	    FieldOperationStruct*
     */
    FieldOperationStruct *plantingOrder;
    FieldOperationStruct *FixedFertilization;
    FieldOperationStruct *Tillage;
    FieldOperationStruct *FixedIrrigation;
    int             operation_index;
    int             i;

        /* If any crop in the community is growing, run the growing crop subroutine */
    if (Community->NumActiveCrop > 0)
        GrowingCrop (rotationYear, y, doy, CropManagement->ForcedHarvest, CropManagement->numHarvest, Community, Residue, SimControl, Soil, SoilCarbon, Weather, Snow, project);

    while (IsOperationToday (rotationYear, doy, CropManagement->plantingOrder, CropManagement->totalCropsPerRotation, &operation_index))
    {
        plantingOrder = &CropManagement->plantingOrder[operation_index];
        PlantingCrop (Community, CropManagement, operation_index);
        if (verbose_mode)
            printf ("DOY %3.3d %-20s %s\n", doy, "Planting", plantingOrder->cropName);
    }
    UpdateOperationStatus (CropManagement->plantingOrder, CropManagement->totalCropsPerRotation);

    while (IsOperationToday (rotationYear, doy, CropManagement->FixedFertilization, CropManagement->numFertilization, &operation_index))
    {
        FixedFertilization = &CropManagement->FixedFertilization[operation_index];
        if (verbose_mode)
            printf ("DOY %3.3d %-20s %s\n", doy, "Fixed Fertilization", FixedFertilization->opSource);

        ApplyFertilizer (FixedFertilization, Soil, Residue);
    }
    UpdateOperationStatus (CropManagement->FixedFertilization, CropManagement->numFertilization);

    while (IsOperationToday (rotationYear, doy, CropManagement->Tillage, CropManagement->numTillage, &operation_index))
    {
        Tillage = &(CropManagement->Tillage[operation_index]);
        if (verbose_mode)
            printf ("DOY %3.3d %-20s %s\n", doy, "Tillage", Tillage->opToolName);

        if (strcasecmp (Tillage->opToolName, "Kill_Crop") != 0)
            ExecuteTillage (SoilCarbon->abgdBiomassInput, Tillage, CropManagement->tillageFactor, Soil, Residue);
        else if (Community->NumActiveCrop > 0)
        {
            for (i = 0; i < Community->NumCrop; i++)
            {
                if (Community->Crop[i].stageGrowth > NO_CROP)
                {
                    HarvestCrop (y, doy, SimControl->simStartYear, &Community->Crop[i], Residue, Soil, SoilCarbon, Weather, project);
                    Community->NumActiveCrop--;
                }
            }
        }
    }
    UpdateOperationStatus (CropManagement->Tillage, CropManagement->numTillage);

    UpdateCommunity (Community);

    Soil->irrigationVol = 0.0;
    while (IsOperationToday (rotationYear, doy, CropManagement->FixedIrrigation, CropManagement->numIrrigation, &operation_index))
    {
        FixedIrrigation = &(CropManagement->FixedIrrigation[operation_index]);
        if (verbose_mode)
            printf ("DOY %3.3d %-20s %lf\n", doy, "Irrigation", FixedIrrigation->opVolume);

        Soil->irrigationVol += FixedIrrigation->opVolume;
    }
    UpdateOperationStatus (CropManagement->FixedIrrigation, CropManagement->numIrrigation);

    ComputeResidueCover (Residue);

    TillageFactorSettling (CropManagement->tillageFactor, Soil->totalLayers, Soil->waterContent, Soil->Porosity);

    SnowProcesses (Snow, y, doy, Weather, Residue->stanResidueTau, Community->svRadiationInterception);

    Redistribution (y, doy, Weather->precipitation[y][doy - 1], Snow->snowFall, Snow->snowMelt, SimControl->hourlyInfiltration, Community, Soil, Residue);

    ResidueEvaporation (Residue, Soil, Community, Weather->ETref[y][doy - 1], Snow->snowCover);

    Evaporation (Soil, Community, Residue, Weather->ETref[y][doy - 1], Snow->snowCover);

    Temperature (y, doy, Snow->snowCover, Community->svRadiationInterception, Soil, Weather, Residue);

    ComputeFactorComposite (SoilCarbon, doy, y, Weather->lastDoy[y], Soil);

    ComputeSoilCarbonBalanceMB (SoilCarbon, y, Residue, Soil, CropManagement->tillageFactor);

    NitrogenTransformation (y, doy, Soil, Community, Residue, Weather, SoilCarbon);
}

void GrowingCrop (int rotationYear, int y, int d, FieldOperationStruct *ForcedHarvest, int numHarvest, CommunityStruct *Community, ResidueStruct *Residue, const SimControlStruct *SimControl, SoilStruct *Soil, SoilCarbonStruct *SoilCarbon, const WeatherStruct *Weather, const SnowStruct *Snow, const char *project)
{
    /*
     * Processes that only occur while a crop is growing are performed
     * Any needed harvests/forages have occured
     * Final Harvest date set once crop maturity achieved
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * forcedHarvest	    int
     * i		    int
     */
    int             forcedHarvest = 0;
    int             i;
    int             clippingFlag = 0;
    double          totalBiomass = 0.0;
    double          clippingBiomassThresholdUpper = -999.0;

    //forcedHarvest = ForcedMaturity (rotationYear, d, Weather->lastDoy[y], *nextSeedingYear, *nextSeedingDate, SimControl->yearsInRotation);

    //for (i = 0; i < numHarvest; i++)
    //{
    //    if (ForcedHarvest[i].opYear == rotationYear && ForcedHarvest[i].opDay == d && ForcedHarvest[i].plantID == Crop->cropUniqueIdentifier)
    //    {
    //        forcedHarvest = 1;
    //        break;
    //    }
    //}

    for (i = 0; i < Community->NumCrop; i++)
    {
        if (Community->Crop[i].stageGrowth > NO_CROP)
        {
            if (Community->Crop[i].svTT_Cumulative < Community->Crop[i].userEmergenceTT)
                Community->Crop[i].stageGrowth = PRE_EMERGENCE;
            else if (Community->Crop[i].svTT_Cumulative < Community->Crop[i].calculatedFloweringTT)
            {
                if (Community->Crop[i].userAnnual)
                    Community->Crop[i].stageGrowth = VEGETATIVE_GROWTH;
                else
                    Community->Crop[i].stageGrowth = PERENNIAL;
            }
            else if (Community->Crop[i].svTT_Cumulative < Community->Crop[i].calculatedMaturityTT)
            {
                if (Community->Crop[i].userAnnual)
                    Community->Crop[i].stageGrowth = REPRODUCTIVE_GROWTH;
                else
                    Community->Crop[i].stageGrowth = PERENNIAL;
            }
            else if (Community->Crop[i].svTT_Cumulative >= Community->Crop[i].calculatedMaturityTT)
            {
                Community->Crop[i].stageGrowth = MATURITY;

                if (Community->Crop[i].harvestDateFinal < 1)
                {
                    /* SetCropStatusToMature */
                    Community->Crop[i].cropMature = 1;
                    Community->Crop[i].harvestDateFinal = FinalHarvestDate (Weather->lastDoy[y], d);
                }
            }
        }
    }

    Phenology (y, d, Weather, Community);

    RadiationInterception (y, d, Community);

    WaterUptake (y, d, Community, Soil, Weather);

    Processes (y, d, SimControl->automaticNitrogen, Community, Residue, Weather, Soil, SoilCarbon);

    /* Calculate total aboveground biomass of all of the community members
     * and find the smallest clipping threshold */
    for (i = 0; i < Community->NumCrop; i++)
    {
        clippingBiomassThresholdUpper = (clippingBiomassThresholdUpper < Community->Crop[i].userClippingBiomassThresholdUpper) ? clippingBiomassThresholdUpper : Community->Crop[i].userClippingBiomassThresholdUpper;
        totalBiomass += Community->Crop[i].svShoot;
    }

    for (i = 0; i < Community->NumCrop; i++)
    {
        if (Community->Crop[i].userClippingTiming > 0.0)
        {
            if (Community->Crop[i].userClippingTiming <= Community->Crop[i].svTT_Cumulative / Community->Crop[i].calculatedMaturityTT ||
                totalBiomass >= clippingBiomassThresholdUpper)
            {
                if ((Community->Crop[i].harvestCount < 3 && Community->Crop[i].userAnnual) || (!Community->Crop[i].userAnnual))
                {
                    clippingFlag = 1;
                    break;
                }
            }
        }
    }

    for (i = 0; i < Community->NumCrop; i++)
    {
        if (Community->Crop[i].stageGrowth > NO_CROP)
        {
            if (Weather->tMin[y][d - 1] < Community->Crop[i].userColdDamageThresholdTemperature)
            {
                if (Community->Crop[i].userAnnual && Community->Crop[i].svTT_Cumulative > Community->Crop[i].calculatedFloweringTT)
                {
                    GrainHarvest (y, d, SimControl->simStartYear, &Community->Crop[i], Residue, Soil, SoilCarbon, Weather, project);
                    Community->NumActiveCrop--;
                }
                else
                    ComputeColdDamage (y, d, &Community->Crop[i], Weather, Snow, Residue);
            }

            if (d == Community->Crop[i].harvestDateFinal || forcedHarvest)
            {
                if (Community->Crop[i].userAnnual)
                {
                    GrainHarvest (y, d, SimControl->simStartYear, &Community->Crop[i], Residue, Soil, SoilCarbon, Weather, project);
                    Community->NumActiveCrop--;
                }
                else
                {
                    ForageHarvest (y, d, SimControl->simStartYear, &Community->Crop[i], Residue, Soil, SoilCarbon, Weather, project);
                    HarvestCrop (y, d, SimControl->simStartYear, &Community->Crop[i], Residue, Soil, SoilCarbon, Weather, project);
                    Community->NumActiveCrop--;
                }
            }
            else if (clippingFlag == 1)
            {
                if (Community->Crop[i].svShoot >= Community->Crop[i].userClippingBiomassThresholdLower * (1.0 - exp (-Community->Crop[i].userPlantingDensity)))
                {
                    ForageHarvest (y, d, SimControl->simStartYear, &Community->Crop[i], Residue, Soil, SoilCarbon, Weather, project);
                    AddCrop (&Community->Crop[i]);
                    Community->Crop[i].stageGrowth = CLIPPING;
                    Community->Crop[i].harvestCount += 1;
                }
            }

            Community->Crop[i].rcCropTranspirationPotential += Community->Crop[i].svTranspirationPotential;
            Community->Crop[i].rcCropTranspiration += Community->Crop[i].svTranspiration;
            Community->Crop[i].rcSoilWaterEvaporation += Soil->evaporationVol;
        }
    }

    UpdateCommunity (Community);
}

//void PlantingCrop (int doy, int *nextSeedingYear, int *nextSeedingDate, CropManagementStruct *CropManagement, CommunityStruct *Community)
//{
//    /*
//     * New realized crop is created next crop in the rotation selected status
//     * set to growing. HarvestDate is reset to an unreachable date.
//     */
//
//    SelectNextCrop (CropManagement);
//    if (verbose_mode)
//        printf ("DOY %3.3d %-20s %s\n", doy, "Planting", CropManagement->plantingOrder[CropManagement->plantingIndex].cropName);
//
//    NewCrop (Community, CropManagement);
//
//    Community->NumActiveCrop++;
//
//    //AddCrop (Crop);
//
//    *nextSeedingYear = CropManagement->nextCropSeedingYear;
//    *nextSeedingDate = CropManagement->nextCropSeedingDate;
//    
//}

double FinalHarvestDate (int lastDoy, int d)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * harvestDate	    int		[return value]
     */
    int             harvestDate;

    harvestDate = d + 10;

    if (harvestDate > lastDoy)
        harvestDate -= lastDoy;

    return (harvestDate);
}

int ForcedMaturity (int rotationYear, int d, int lastDoy, int nextSeedingYear, int nextSeedingDate, int rotationSize)
{
    /*
     * Returns true if planted crop is within saftey margin of days of the
     * next crop to be planted
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * forced_maturity	    int
     * margin		    int
     * nextRotationYear	    int
     */
    int             forced_maturity = 0;
    int             margin = 9;
    int             nextRotationYear;

    nextRotationYear = rotationYear;

    if (nextRotationYear < rotationSize)
        nextRotationYear++;
    else
        nextRotationYear = 1;

    if (rotationYear == nextSeedingYear)
    {
        if (d < nextSeedingDate && (d + margin) >= nextSeedingDate)
            forced_maturity = 1;
    }
    else if (nextRotationYear == nextSeedingYear)
    {
        if (d >= (lastDoy + nextSeedingDate - margin))
            forced_maturity = 1;
    }

    return (forced_maturity);
}
