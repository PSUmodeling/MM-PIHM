#include "pihm.h"

void InitCycles (elem_struct *elem, river_struct *riv,
    const ctrl_struct *ctrl, const mgmttbl_struct *mgmttbl,
    const agtbl_struct *agtbl, const croptbl_struct *croptbl,
    const soiltbl_struct *soiltbl)
{
    int             opind;
    int             i, j, k;
    int             start_year;
    int             end_year;
    int             total_years;
    int             y;
    struct tm      *timestamp;
    time_t          rawtime;
    comm_struct    *comm;
    weather_struct *weather;
    cropmgmt_struct *cropmgmt;

    rawtime = (time_t) ctrl->starttime;
    timestamp = gmtime (&rawtime);
    start_year = timestamp->tm_year + 1900;

    rawtime = (time_t) (ctrl->endtime - DAYINSEC);
    timestamp = gmtime (&rawtime);
    end_year = timestamp->tm_year + 1900;

    total_years = end_year - start_year + 1;

    for (i = 0; i < nelem; i++)
    {
        /*
         * Initialize weather structure
         */
        weather = &elem[i].weather;

        weather->wind = (double **)malloc (total_years * sizeof (double *));
        //weather->ETref = (double **)malloc (total_years * sizeof (double *));
        //weather->precipitation = (double **)malloc (total_years * sizeof (double *));
        //weather->RHmax = (double **)malloc (total_years * sizeof (double *));
        //weather->RHmin = (double **)malloc (total_years * sizeof (double *));
        weather->solarRadiation =
            (double **)malloc (total_years * sizeof (double *));
        weather->tMax = (double **)malloc (total_years * sizeof (double *));
        weather->tMin = (double **)malloc (total_years * sizeof (double *));
        weather->vpd = (double **)malloc (total_years * sizeof (double *));
        //weather->yearlyAmplitude = (double *)malloc (total_years * sizeof (double));
        //weather->annualAverageTemperature = (double *)malloc (total_years * sizeof (double));
        weather->lastDoy = (int *)malloc (total_years * sizeof (int));
        for (y = 0; y < total_years; y++)
        {
            weather->wind[y] = (double *)malloc (366 * sizeof (double));
            //weather->ETref[y] = (double *)malloc (366 * sizeof (double));
            //weather->precipitation[y] = (double *)malloc (366 * sizeof (double));
            //weather->RHmax[y] = (double *)malloc (366 * sizeof (double));
            //weather->RHmin[y] = (double *)malloc (366 * sizeof (double));
            weather->solarRadiation[y] =
                (double *)malloc (366 * sizeof (double));
            weather->tMax[y] = (double *)malloc (366 * sizeof (double));
            weather->tMin[y] = (double *)malloc (366 * sizeof (double));
            weather->vpd[y] = (double *)malloc (366 * sizeof (double));
        }

        /*
         * Copy crop management to each element
         */
        cropmgmt = &elem[i].cropmgmt;

        if (agtbl->op[i] > agtbl->nopfile)
        {
            Cycles_printf (VL_ERROR,
                "ERROR: Operation file for operation index %d "
                "is not provided.\n", agtbl->op[i]);
            Cycles_exit (EXIT_FAILURE);
        }

        opind = agtbl->op[i] - 1;

        cropmgmt->yearsInRotation = agtbl->rotsz[i];
        cropmgmt->automaticNitrogen = agtbl->auto_N[i];
        cropmgmt->automaticPhosphorus = agtbl->auto_P[i];
        cropmgmt->automaticSulfur = agtbl->auto_S[i];

        cropmgmt->rotationYear = 0;

        cropmgmt->totalCropsPerRotation = mgmttbl->cropmgmt[opind].totalCropsPerRotation;
        cropmgmt->plantingOrder = mgmttbl->cropmgmt[opind].plantingOrder;

        cropmgmt->numTillage = mgmttbl->cropmgmt[opind].numTillage;
        cropmgmt->Tillage = mgmttbl->cropmgmt[opind].Tillage;

        cropmgmt->numFertilization = mgmttbl->cropmgmt[opind].numFertilization;
        cropmgmt->FixedFertilization = mgmttbl->cropmgmt[opind].FixedFertilization;

        cropmgmt->numIrrigation = mgmttbl->cropmgmt[opind].numIrrigation;
        cropmgmt->FixedIrrigation = mgmttbl->cropmgmt[opind].FixedIrrigation;

        cropmgmt->numAutoIrrigation = mgmttbl->cropmgmt[opind].numAutoIrrigation;
        cropmgmt->autoIrrigation = mgmttbl->cropmgmt[opind].autoIrrigation;

        cropmgmt->op_status[PLANT_OP] =
            (int *)calloc (cropmgmt->totalCropsPerRotation, sizeof (int));
        cropmgmt->op_status[TILLAGE_OP] =
            (int *)calloc (cropmgmt->numTillage, sizeof (int));
        cropmgmt->op_status[FIXIRR_OP] =
            (int *)calloc (cropmgmt->numIrrigation, sizeof (int));
        cropmgmt->op_status[FIXFERT_OP] =
            (int *)calloc (cropmgmt->numFertilization, sizeof (int));

        for (k = 0; k < elem[i].ps.nsoil; k++)
        {
            cropmgmt->tillageFactor[k] = 0.0;
        }

        cropmgmt->usingAutoIrr = mgmttbl->cropmgmt[opind].usingAutoIrr;
        cropmgmt->usingAutoFert = mgmttbl->cropmgmt[opind].usingAutoFert;

        /*
         * Copy crop description to each element
         */
        comm = &elem[i].comm;

        comm->NumCrop = croptbl->number;
        comm->Crop =
            (crop_struct *)malloc (comm->NumCrop * sizeof (crop_struct));

        for (j = 0; j < comm->NumCrop; j++)
        {
            strcpy (comm->Crop[j].cropName, croptbl->cropName[j]);

            comm->Crop[j].userFloweringTT = croptbl->userFloweringTT[j];
            comm->Crop[j].userMaturityTT = croptbl->userMaturityTT[j];
            comm->Crop[j].userMaximumSoilCoverage =
                croptbl->userMaximumSoilCoverage[j];
            comm->Crop[j].userMaximumRootingDepth =
                croptbl->userMaximumRootingDepth[j];
            comm->Crop[j].userExpectedYieldAvg =
                croptbl->userExpectedYieldAvg[j];
            comm->Crop[j].userExpectedYieldMax =
                croptbl->userExpectedYieldMax[j];
            comm->Crop[j].userExpectedYieldMin =
                croptbl->userExpectedYieldMin[j];
            comm->Crop[j].userPercentMoistureInYield =
                croptbl->userPercentMoistureInYield[j];
            comm->Crop[j].userFractionResidueStanding =
                croptbl->userFractionResidueStanding[j];
            comm->Crop[j].userFractionResidueRemoved =
                croptbl->userFractionResidueRemoved[j];
            comm->Crop[j].userClippingBiomassThresholdUpper =
                croptbl->userClippingBiomassThresholdUpper[j];
            comm->Crop[j].userClippingBiomassThresholdLower =
                croptbl->userClippingBiomassThresholdLower[j];
            comm->Crop[j].userClippingTiming = croptbl->userClippingTiming[j];
            comm->Crop[j].userClippingDestiny =
                croptbl->userClippingDestiny[j];
            comm->Crop[j].userTranspirationMinTemperature =
                croptbl->userTranspirationMinTemperature[j];
            comm->Crop[j].userTranspirationThresholdTemperature =
                croptbl->userTranspirationThresholdTemperature[j];
            comm->Crop[j].userColdDamageMinTemperature =
                croptbl->userColdDamageMinTemperature[j];
            comm->Crop[j].userColdDamageThresholdTemperature =
                croptbl->userColdDamageThresholdTemperature[j];
            comm->Crop[j].userTemperatureBase =
                croptbl->userTemperatureBase[j];
            comm->Crop[j].userTemperatureOptimum =
                croptbl->userTemperatureOptimum[j];
            comm->Crop[j].userTemperatureMaximum =
                croptbl->userTemperatureMaximum[j];
            comm->Crop[j].userShootPartitionInitial =
                croptbl->userShootPartitionInitial[j];
            comm->Crop[j].userShootPartitionFinal =
                croptbl->userShootPartitionFinal[j];
            comm->Crop[j].userRadiationUseEfficiency =
                croptbl->userRadiationUseEfficiency[j];
            comm->Crop[j].userTranspirationUseEfficiency =
                croptbl->userTranspirationUseEfficiency[j];
            comm->Crop[j].userHIx = croptbl->userHIx[j];
            comm->Crop[j].userHIo = croptbl->userHIo[j];
            comm->Crop[j].userHIk = croptbl->userHIk[j];
            comm->Crop[j].userEmergenceTT = croptbl->userEmergenceTT[j];
            comm->Crop[j].userNMaxConcentration =
                croptbl->userNMaxConcentration[j];
            comm->Crop[j].userNDilutionSlope = croptbl->userNDilutionSlope[j];
            comm->Crop[j].userKc = croptbl->userKc[j];
            comm->Crop[j].userAnnual = croptbl->userAnnual[j];
            comm->Crop[j].userLegume = croptbl->userLegume[j];
            comm->Crop[j].userC3orC4 = croptbl->userC3orC4[j];
            comm->Crop[j].userExtinctionCoefficient =
                croptbl->userExtinctionCoefficient[j];
            comm->Crop[j].userPlantingDensity =
                croptbl->userPlantingDensity[j];
            comm->Crop[j].userClippingStart = croptbl->userClippingStart[j];
            comm->Crop[j].userClippingEnd = croptbl->userClippingEnd[j];
            comm->Crop[j].LWP_StressOnset = croptbl->LWP_StressOnset[j];
            comm->Crop[j].LWP_WiltingPoint = croptbl->LWP_WiltingPoint[j];
            comm->Crop[j].transpirationMax = croptbl->transpirationMax[j];

            comm->Crop[j].userMaximumSoilCoverage *= 0.94 / 100.0;
            comm->Crop[j].userPercentMoistureInYield /= 100.0;
            comm->Crop[j].userFractionResidueStanding /= 100.0;
            comm->Crop[j].userFractionResidueRemoved /= 100.0;
            if (comm->Crop[j].userClippingTiming != BADVAL)
            {
                comm->Crop[j].userClippingTiming /= 100.0;
            }
            else
            {
                comm->Crop[j].userClippingTiming = 0.0;
            }

            comm->Crop[j].calculatedSimAvgYield = 0.0;
            comm->Crop[j].calculatedSimMaxYield = 0.0;
            comm->Crop[j].calculatedSimMinYield = 0.0;

            comm->Crop[j].stageGrowth = NO_CROP;

            InitCropSV (&comm->Crop[j]);
        }

        UpdateCommunity (comm);

        InitializeSoil (&elem[i].soil, soiltbl, &elem[i].ps,
            elem[i].attrib.soil_type);

        InitializeResidue (&elem[i].residue, elem[i].ps.nsoil);

        InitializeSoilCarbon (&elem[i].soilc, elem[i].ps.nsoil);

        elem[i].comm.NumActiveCrop = 0;

    }

    for (i = 0; i < nriver; i++)
    {
        riv[i].NO3sol.soluteMass = 0.0;
        riv[i].NH4sol.soluteMass = 0.0;

        for (k = 0; k < elem[riv[i].rightele - 1].ps.nsoil; k++)
        {
            riv[i].NO3sol.soluteMass +=
                elem[riv[i].rightele - 1].soil.NO3[k];
            riv[i].NH4sol.soluteMass +=
                elem[riv[i].rightele - 1].soil.NH4[k];
        }

        for (k = 0; k < elem[riv[i].leftele - 1].ps.nsoil; k++)
        {
            riv[i].NO3sol.soluteMass +=
                elem[riv[i].leftele - 1].soil.NO3[k];
            riv[i].NH4sol.soluteMass +=
                elem[riv[i].leftele - 1].soil.NH4[k];
        }
    }
}
