#include "Cycles.h"

//void SelectCropInitialPosition (CropManagementStruct *CropManagement)
//{
//    CropManagement->plantingIndex = -1;
//    //CropManagement->describedIndex = -1;
//    PeekNextCrop (CropManagement);
//}
//
//void SelectNextCrop (CropManagementStruct *CropManagement)
//{
//    /*
//     * Planting index set to next crop (if any) to be planted in the rotation
//     * Described index set to previously searched value inside plantingOrder
//     * Next crop in the rotation (if any)has limited values stored
//     */
//    if (CropManagement->plantingIndex < CropManagement->totalCropsPerRotation - 1)
//        CropManagement->plantingIndex++;
//    else
//        CropManagement->plantingIndex = 0;
//
//    //CropManagement->describedIndex = CropManagement->plantingOrder[CropManagement->plantingIndex].plantID;
//
//    PeekNextCrop (CropManagement);
//}
//
//void PeekNextCrop (CropManagementStruct *CropManagement)
//{
//    /*
//     * Store the name, year, day, and type values for the next crop
//     */
//    int             tempIndex;
//
//    tempIndex = CropManagement->plantingIndex;
//
//    if (CropManagement->totalCropsPerRotation > 0)
//    {
//        /* Select next crop */
//        if (tempIndex < CropManagement->totalCropsPerRotation - 1)
//            tempIndex++;
//        else
//            tempIndex = 0;
//
//        CropManagement->nextCropSeedingYear = CropManagement->plantingOrder[tempIndex].opYear;
//        CropManagement->nextCropSeedingDate = CropManagement->plantingOrder[tempIndex].opDay;
//        strcpy (CropManagement->nextCropName, CropManagement->plantingOrder[tempIndex].cropName);
//    }
//}

void PlantingCrop (CommunityStruct *Community, const CropManagementStruct *CropManagement, int plantingIndex)
{
    FieldOperationStruct *plantingOrder;
    autoIrrigationStruct *autoIrrigation;
    CropStruct           *Crop;
    //describedCropStruct *describedCrop;

    plantingOrder = &(CropManagement->plantingOrder[plantingIndex]);
    //describedCrop = &(CropManagement->describedCrop[CropManagement->describedIndex]);

    //Crop->cropUniqueIdentifier = plantingOrder->plantID;
    //strcpy (Crop->cropName, plantingOrder->cropName);

    Crop = &Community->Crop[plantingOrder->plantID];

    Crop->stageGrowth = PLANTING;

    Community->NumActiveCrop++;

    if (plantingOrder->usesAutoIrrigation > 0)
    {
        autoIrrigation = &(CropManagement->autoIrrigation[plantingOrder->usesAutoIrrigation]);
        Crop->autoIrrigationUsed = 1;
        Crop->autoIrrigationStartDay = autoIrrigation->startDay;
        Crop->autoIrrigationStopDay = autoIrrigation->stopDay;
        Crop->autoIrrigationWaterDepletion = autoIrrigation->waterDepletion;
        Crop->autoIrrigationLastSoilLayer = autoIrrigation->lastSoilLayer;
    }
    else
        Crop->autoIrrigationUsed = 0;

    Crop->userPlantingDensity = plantingOrder->plantingDensity;

    Crop->svTT_Daily = 0.0;
    Crop->svTT_Cumulative = 0.0;
    Crop->svRadiationInterception = 0.0;
    Crop->svBiomass = 0.0;
    Crop->svShoot = 0.0;
    Crop->svRoot = 0.0;
    Crop->svRizho = 0.0;
    Crop->svShootDailyGrowth = 0.0;
    Crop->svRootDailyGrowth = 0.0;
    Crop->svRizhoDailyDeposition = 0.0;
    Crop->svUnstressedShootDailyGrowth = 0.0;
    Crop->svUnstressedRootDailyGrowth = 0.0;
    Crop->svPostFloweringShootBiomass = 0.0;
    Crop->svRootingDepth = 0.0;
    Crop->svTranspiration = 0.0;
    Crop->svTranspirationPotential = 0.0;
    Crop->svN_Shoot = 0.0;
    Crop->svN_Root = 0.0;
    Crop->svN_Rhizo = 0.0;
    Crop->svN_RizhoDailyDeposition = 0.0;
    Crop->svN_AutoAdded = 0.0;
    Crop->svN_Fixation = 0.0;
    Crop->svWaterStressFactor = 0.0;
    Crop->svN_StressFactor = 0.0;
    Crop->svShootUnstressed = 0.0;
    Crop->svN_StressCumulative = 0.0;

    Crop->svRadiationInterception_nc = 0.0;

    Crop->rcForageYield = 0.0;
    Crop->rcGrainYield = 0.0;
    Crop->rcBiomass = 0.0;
    Crop->rcRoot = 0.0;
    Crop->rcResidueBiomass = 0.0;
    Crop->rcCropTranspiration = 0.0;
    Crop->rcCropTranspirationPotential = 0.0;
    Crop->rcSoilWaterEvaporation = 0.0;
    Crop->rcHarvestIndex = 0;
    Crop->rcTotalNitrogen = 0.0;
    Crop->rcRootNitrogen = 0.0;
    Crop->rcGrainNitrogenYield = 0.0;
    Crop->rcForageNitrogenYield = 0.0;
    Crop->rcNitrogenCumulative = 0.0;
    Crop->rcNitrogenInHarvest = 0.0;
    Crop->rcNitrogenInResidue = 0.0;
    Crop->rcNitrogenForageConc = 0.0;
    //Crop->userSeedingDate = describedCrop->userSeedingDate;
    //Crop->userFloweringDate = describedCrop->userFloweringDate;
    //Crop->userMaturityDate = describedCrop->userMaturityDate;
    //Crop->userMaximumSoilCoverage = describedCrop->userMaximumSoilCoverage;
    //Crop->userMaximumRootingDepth = describedCrop->userMaximumRootingDepth;
    //Crop->userExpectedYieldAvg = describedCrop->userExpectedYieldAvg;
    //Crop->userExpectedYieldMax = describedCrop->userExpectedYieldMax;
    //Crop->userExpectedYieldMin = describedCrop->userExpectedYieldMin;
    //Crop->userPercentMoistureInYield = describedCrop->userPercentMoistureInYield;
    //Crop->userFractionResidueStanding = describedCrop->userFractionResidueStanding;
    //Crop->userFractionResidueRemoved = describedCrop->userFractionResidueRemoved;
    //Crop->userClippingTiming = describedCrop->userClippingTiming;
    //Crop->userTranspirationMinTemperature = describedCrop->userTranspirationMinTemperature;
    //Crop->userTranspirationThresholdTemperature = describedCrop->userTranspirationThresholdTemperature;
    //Crop->userColdDamageMinTemperature = describedCrop->userColdDamageMinTemperature;
    //Crop->userColdDamageThresholdTemperature = describedCrop->userColdDamageThresholdTemperature;
    //Crop->userTemperatureBase = describedCrop->userTemperatureBase;
    //Crop->userTemperatureOptimum = describedCrop->userTemperatureOptimum;
    //Crop->userTemperatureMaximum = describedCrop->userTemperatureMaximum;
    //Crop->userShootPartitionInitial = describedCrop->userShootPartitionInitial;
    //Crop->userShootPartitionFinal = describedCrop->userShootPartitionFinal;
    //Crop->userRadiationUseEfficiency = describedCrop->userRadiationUseEfficiency;
    //Crop->userTranspirationUseEfficiency = describedCrop->userTranspirationUseEfficiency;
    //Crop->userHIx = describedCrop->userHIx;
    //Crop->userHIo = describedCrop->userHIo;
    //Crop->userHIk = describedCrop->userHIk;
    //Crop->userEmergenceTT = describedCrop->userEmergenceTT;
    //Crop->userNMaxConcentration = describedCrop->userNMaxConcentration;
    //Crop->userNDilutionSlope = describedCrop->userNDilutionSlope;
    //Crop->userKc = describedCrop->userKc;
    //Crop->userAnnual = describedCrop->userAnnual;
    //Crop->userLegume = describedCrop->userLegume;
    //Crop->userC3orC4 = describedCrop->userC3orC4;

    //Crop->calculatedFloweringTT = describedCrop->calculatedFloweringTT;
    //Crop->calculatedMaturityTT = describedCrop->calculatedMaturityTT;
    //Crop->calculatedSimAvgYield = describedCrop->calculatedSimAvgYield;
    //Crop->calculatedSimMaxYield = describedCrop->calculatedSimMaxYield;
    //Crop->calculatedSimMinYield = describedCrop->calculatedSimMinYield;

    Crop->harvestDateFinal = -1;
    Crop->harvestCount = 0;
    //Crop->stageGrowth = -1;
}

void AddCrop (CropStruct *Crop)
{
    Crop->rcForageYield = 0.0;
    Crop->rcGrainYield = 0.0;
    Crop->rcBiomass = 0.0;
    Crop->rcRoot = 0.0;
    Crop->rcResidueBiomass = 0.0;
    Crop->rcCropTranspiration = 0.0;
    Crop->rcCropTranspirationPotential = 0.0;
    Crop->rcSoilWaterEvaporation = 0.0;
    Crop->rcHarvestIndex = 0;
    Crop->rcTotalNitrogen = 0.0;
    Crop->rcRootNitrogen = 0.0;
    Crop->rcGrainNitrogenYield = 0.0;
    Crop->rcForageNitrogenYield = 0.0;
    Crop->rcNitrogenCumulative = 0.0;
    Crop->rcNitrogenInHarvest = 0.0;
    Crop->rcNitrogenInResidue = 0.0;
    Crop->rcNitrogenForageConc = 0.0;

}

void KillCrop (CropStruct *Crop)
{
    Crop->stageGrowth = NO_CROP;
    //strcpy (Crop->cropName, "\0");
    Crop->autoIrrigationUsed = -1;
    Crop->autoIrrigationStartDay = 0;
    Crop->autoIrrigationStopDay = 0;
    Crop->autoIrrigationWaterDepletion = 0.0;
    Crop->autoIrrigationLastSoilLayer = 0;

    Crop->userPlantingDensity = 0.0;
    //Crop->userSeedingDate = 0;
    //Crop->userFloweringDate = 0;
    //Crop->userMaturityDate = 0;
    //Crop->userMaximumSoilCoverage = 0.0;
    //Crop->userMaximumRootingDepth = 0.0;
    //Crop->userExpectedYieldAvg = 0.0;
    //Crop->userExpectedYieldMax = 0.0;
    //Crop->userExpectedYieldMin = 0.0;
    //Crop->userPercentMoistureInYield = 0.0;
    //Crop->userFractionResidueStanding = 0.0;
    //Crop->userFractionResidueRemoved = 0.0;
    //Crop->userClippingTiming = 0.0;
    //Crop->userTranspirationMinTemperature = 0.0;
    //Crop->userTranspirationThresholdTemperature = 0.0;
    //Crop->userColdDamageMinTemperature = 0.0;
    //Crop->userColdDamageThresholdTemperature = 0.0;
    //Crop->userTemperatureBase = 0.0;
    //Crop->userTemperatureOptimum = 0.0;
    //Crop->userTemperatureMaximum = 0.0;
    //Crop->userShootPartitionInitial = 0.0;
    //Crop->userShootPartitionFinal = 0.0;
    //Crop->userRadiationUseEfficiency = 0.0;
    //Crop->userTranspirationUseEfficiency = 0.0;
    //Crop->userHIx = 0.0;
    //Crop->userHIo = 0.0;
    //Crop->userHIk = 0.0;
    //Crop->userEmergenceTT = 0.0;
    //Crop->userNMaxConcentration = 0.0;
    //Crop->userNDilutionSlope = 0.0;
    //Crop->userKc = 0.0;
    //Crop->userAnnual = 0;
    //Crop->userLegume = 0;
    //Crop->userC3orC4 = 0;
    //Crop->calculatedFloweringTT = 0.0;
    //Crop->calculatedMaturityTT = 0.0;
    //Crop->calculatedSimAvgYield = 0.0;
    //Crop->calculatedSimMaxYield = 0.0;
    //Crop->calculatedSimMinYield = 0.0;
    Crop->harvestDateFinal = -1;
    Crop->harvestCount = -1;
    //Crop->stageGrowth = -1;

    Crop->svTT_Daily = 0.0;
    Crop->svTT_Cumulative = 0.0;
    Crop->svRadiationInterception = 0.0;
    Crop->svRadiationInterception_nc = 0.0;
    Crop->svBiomass = 0.0;
    Crop->svShoot = 0.0;
    Crop->svRoot = 0.0;
    Crop->svRizho = 0.0;
    Crop->svShootDailyGrowth = 0.0;
    Crop->svRootDailyGrowth = 0.0;
    Crop->svRizhoDailyDeposition = 0.0;
    Crop->svUnstressedShootDailyGrowth = 0.0;
    Crop->svUnstressedRootDailyGrowth = 0.0;
    Crop->svPostFloweringShootBiomass = 0.0;
    Crop->svRootingDepth = 0.0;
    Crop->svTranspiration = 0.0;
    Crop->svTranspirationPotential = 0.0;
    Crop->svN_Shoot = 0.0;
    Crop->svN_Root = 0.0;
    Crop->svN_Rhizo = 0.0;
    Crop->svN_RizhoDailyDeposition = 0.0;
    Crop->svN_AutoAdded = 0.0;
    Crop->svN_Fixation = 0.0;
    Crop->svWaterStressFactor = 0.0;
    Crop->svN_StressFactor = 0.0;

    Crop->svShootUnstressed = 0.0;
    Crop->svN_StressCumulative = 0.0;

    Crop->rcForageYield = 0.0;
    Crop->rcGrainYield = 0.0;
    Crop->rcBiomass = 0.0;
    Crop->rcRoot = 0.0;
    Crop->rcResidueBiomass = 0.0;
    Crop->rcCropTranspiration = 0.0;
    Crop->rcCropTranspirationPotential = 0.0;
    Crop->rcSoilWaterEvaporation = 0.0;
    Crop->rcHarvestIndex = 0.0;
    Crop->rcTotalNitrogen = 0.0;
    Crop->rcRootNitrogen = 0.0;
    Crop->rcGrainNitrogenYield = 0.0;
    Crop->rcForageNitrogenYield = 0.0;
    Crop->rcNitrogenCumulative = 0.0;
    Crop->rcNitrogenInHarvest = 0.0;
    Crop->rcNitrogenInResidue = 0.0;
    Crop->rcNitrogenForageConc = 0.0;
}

void UpdateCommunity (CommunityStruct *Community)
{
    int         i;

    Community->svRadiationInterception = 0.0;
    Community->svBiomass = 0.0;
    Community->svShoot = 0.0;
    Community->svRoot = 0.0;
    Community->svRizho = 0.0;
    Community->svShootDailyGrowth = 0.0;
    Community->svRootDailyGrowth = 0.0;
    Community->svRizhoDailyDeposition = 0.0;
    Community->svRootingDepth = 0.0;
    Community->svTranspiration = 0.0;
    Community->svTranspirationPotential = 0.0;
    Community->svN_Shoot = 0.0;
    Community->svN_Root = 0.0;
    Community->svN_Rhizo = 0.0;
    Community->svN_RizhoDailyDeposition = 0.0;
    Community->svN_AutoAdded = 0.0;
    Community->svN_Fixation = 0.0;
    Community->svWaterStressFactor = 0.0;
    Community->svN_StressFactor = 0.0;
    
    for (i = 0; i < Community->NumCrop; i++)
    {
        if (Community->Crop[i].stageGrowth > NO_CROP)
        {
            Community->svRadiationInterception += Community->Crop[i].svRadiationInterception;
            Community->svBiomass += Community->Crop[i].svBiomass;
            Community->svShoot += Community->Crop[i].svShoot;
            Community->svRoot += Community->Crop[i].svRoot;
            Community->svRizho += Community->Crop[i].svRizho;
            Community->svShootDailyGrowth += Community->Crop[i].svShootDailyGrowth;
            Community->svRootDailyGrowth += Community->Crop[i].svRootDailyGrowth;
            Community->svRizhoDailyDeposition += Community->Crop[i].svRizhoDailyDeposition;
            Community->svRootingDepth = Community->Crop[i].svRootingDepth > Community->svRootingDepth ? Community->Crop[i].svRootingDepth : Community->svRootingDepth;
            Community->svTranspiration += Community->Crop[i].svTranspiration;
            Community->svTranspirationPotential += Community->Crop[i].svTranspirationPotential;
            Community->svN_Shoot += Community->Crop[i].svN_Shoot;
            Community->svN_Root += Community->Crop[i].svN_Root;
            Community->svN_Rhizo += Community->Crop[i].svN_Rhizo;
            Community->svN_RizhoDailyDeposition += Community->Crop[i].svN_RizhoDailyDeposition;
            Community->svN_AutoAdded += Community->Crop[i].svN_AutoAdded;
            Community->svN_Fixation += Community->Crop[i].svN_Fixation;
            Community->svWaterStressFactor += Community->Crop[i].svWaterStressFactor;
            Community->svN_StressFactor += Community->Crop[i].svN_StressFactor;
        }
    }
} 
