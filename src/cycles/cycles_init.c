#include "pihm.h"

void InitCycles(const agtbl_struct *agtbl, const soiltbl_struct *soiltbl,
    epconst_struct epctbl[], elem_struct elem[], river_struct river[])
{
    int             soil_ind;
    int             i, j, k;

    for (i = 0; i < nelem; i++)
    {
        /*
         * Initialize initial soil variables
         */
        soil_ind = elem[i].attrib.soil_type - 1;

        for (k = 0; k < MAXLYR; k++)
        {
            if (k < elem[i].ps.nsoil)
            {
                elem[i].soil.clay[k] = soiltbl->clay_lyr[soil_ind][k];
                elem[i].soil.clay[k] *= 0.01;
                elem[i].soil.sand[k] = soiltbl->sand_lyr[soil_ind][k];
                elem[i].soil.sand[k] *= 0.01;
                elem[i].soil.iom[k] = soiltbl->iom_lyr[soil_ind][k];
                elem[i].soil.bd[k] = soiltbl->bd_lyr[soil_ind][k];
            }
            else
            {
                elem[i].soil.clay[k] = BADVAL;
                elem[i].soil.sand[k] = BADVAL;
                elem[i].soil.iom[k] = BADVAL;
                elem[i].soil.bd[k] = BADVAL;
            }
        }

        /*
         * Initialize crop structures
         */
        for (j = 0; j < MAXCROP; j++)
        {
            elem[i].crop[j].epc = &epctbl[j];

            InitCropSV(&elem[i].crop[j]);
        }

        /*
         * Initialize management structure
         */
        elem[i].mgmt.rot_size = agtbl->rotsz[i];
        elem[i].mgmt.auto_n = agtbl->auto_N[i];

        elem[i].mgmt.rot_year = 0;
        for (j = 0; j < 4; j++)
        {
            elem[i].mgmt.op_ptr[j] = 0;
        }
    }
}

void FirstDay(const soiltbl_struct *soiltbl, elem_struct elem[],
    river_struct river[])
{
    int             i, k;
    int             soil_ind;

    for (i = 0; i < nelem; i++)
    {
        soil_ind = elem[i].attrib.soil_type - 1;

        elem[i].restart_input.resw_stan = 0.0;
        elem[i].restart_input.resw_flat = 0.0;
        elem[i].restart_input.resm_stan = 0.0;
        elem[i].restart_input.resm_flat = 0.0;
        elem[i].restart_input.manuc_surf = 0.0;
        elem[i].restart_input.resn_stan = 0.0;
        elem[i].restart_input.resn_flat = 0.0;
        elem[i].restart_input.manun_surf = 0.0;

        for (k = 0; k < MAXLYR; k++)
        {
            elem[i].restart_input.res_abgd[k] = BADVAL;
            elem[i].restart_input.res_root[k] = BADVAL;
            elem[i].restart_input.res_rhizo[k] = BADVAL;
            elem[i].restart_input.manuc[k] = BADVAL;
            elem[i].restart_input.resn_abgd[k] = BADVAL;
            elem[i].restart_input.resn_root[k] = BADVAL;
            elem[i].restart_input.resn_rhizo[k] = BADVAL;
            elem[i].restart_input.manun[k] = BADVAL;
            elem[i].restart_input.soc[k] = BADVAL;
            elem[i].restart_input.son[k] = BADVAL;
            elem[i].restart_input.mbc[k] = BADVAL;
            elem[i].restart_input.mbn[k] = BADVAL;
            elem[i].restart_input.no3[k] = BADVAL;
            elem[i].restart_input.nh4[k] = BADVAL;
        }

        for (k = 0; k < elem[i].ps.nsoil; k++)
        {
            elem[i].restart_input.soc[k] = elem[i].soil.iom[k] / 100.0 * 0.58 *
                elem[i].ps.sldpth[k] * elem[i].soil.bd[k] * 1000.0;
            /* Initializes as 3% of SOC_Mass but "added" C */
            elem[i].restart_input.mbc[k] = 0.03 * elem[i].restart_input.soc[k];
            elem[i].restart_input.res_abgd[k] = 0.0;
            elem[i].restart_input.res_root[k] = 0.0;
            elem[i].restart_input.res_rhizo[k] = 0.0;
            elem[i].restart_input.manuc[k] = 0.0;

            elem[i].restart_input.no3[k] =
                soiltbl->no3_lyr[soil_ind][k] * 1.0E-4;
            elem[i].restart_input.nh4[k] =
                soiltbl->nh4_lyr[soil_ind][k] * 1.0E-4;
            /* Initializes with CN ratio = 10 */
            elem[i].restart_input.son[k] = elem[i].restart_input.soc[k] * 0.1;
            /* Initializes with CN ratio = 10 */
            elem[i].restart_input.mbn[k] = elem[i].restart_input.mbc[k] * 0.1;
            elem[i].restart_input.resn_abgd[k] = 0.0;
            elem[i].restart_input.resn_root[k] = 0.0;
            elem[i].restart_input.resn_rhizo[k] = 0.0;
            elem[i].restart_input.manun[k] = 0.0;
        }
    }

    for (i = 0; i < nriver; i++)
    {
        river[i].restart_input.streamno3 = 0.0;
        river[i].restart_input.streamnh4 = 0.0;
        river[i].restart_input.bedno3 = 0.0;
        river[i].restart_input.bednh4 = 0.0;
    }
}

void InitCyclesVar(elem_struct elem[], river_struct river[], N_Vector CV_Y)
{
    int             i, k;

    for (i = 0; i < nelem; i++)
    {
        RestartInput(&elem[i].restart_input, &elem[i].ps, &elem[i].ws,
            &elem[i].cs, &elem[i].ns);

        elem[i].np.no3 = 0.0;
        elem[i].np.nh4 = 0.0;
        for (k = 0; k < elem[i].ps.nsoil; k++)
        {
            elem[i].np.no3 += elem[i].ns.no3[k];
            elem[i].np.nh4 += elem[i].ns.nh4[k];
        }

        NV_Ith(CV_Y, NO3(i)) = elem[i].np.no3;
        NV_Ith(CV_Y, NH4(i)) = elem[i].np.nh4;

        elem[i].np0 = elem[i].np;

        elem[i].no3sol.snksrc = 0.0;
        elem[i].nh4sol.snksrc = 0.0;
    }

    for (i = 0; i < nriver; i++)
    {
        river[i].ns.streamno3 = river[i].restart_input.streamno3;
        river[i].ns.streamnh4 = river[i].restart_input.streamnh4;
        river[i].ns.bedno3 = river[i].restart_input.bedno3;
        river[i].ns.bednh4 = river[i].restart_input.bednh4;

        NV_Ith(CV_Y, STREAMNO3(i)) = river[i].ns.streamno3;
        NV_Ith(CV_Y, STREAMNH4(i)) = river[i].ns.streamnh4;
        NV_Ith(CV_Y, RIVBEDNO3(i)) = river[i].ns.bedno3;
        NV_Ith(CV_Y, RIVBEDNH4(i)) = river[i].ns.bednh4;
    }
}

//void InitCycles(elem_struct *elem, river_struct *riv,
//    const ctrl_struct *ctrl, const mgmttbl_struct *mgmttbl,
//    const agtbl_struct *agtbl, const croptbl_struct *croptbl,
//    const soiltbl_struct *soiltbl)
//{
//    int             opind;
//    int             i, j, k;
//    int             start_year;
//    int             end_year;
//    int             total_years;
//    int             y;
//    struct tm      *timestamp;
//    time_t          rawtime;
//    comm_struct    *comm;
//    weather_struct *weather;
//    cropmgmt_struct *cropmgmt;
//
//    rawtime = (time_t)ctrl->starttime;
//    timestamp = gmtime(&rawtime);
//    start_year = timestamp->tm_year + 1900;
//
//    rawtime = (time_t)(ctrl->endtime - DAYINSEC);
//    timestamp = gmtime(&rawtime);
//    end_year = timestamp->tm_year + 1900;
//
//    total_years = end_year - start_year + 1;
//
//    for (i = 0; i < nelem; i++)
//    {
//        /*
//         * Initialize weather structure
//         */
//        weather = &elem[i].weather;
//
//        weather->wind = (double **)malloc(total_years * sizeof(double *));
//        //weather->ETref = (double **)malloc (total_years * sizeof(double *));
//        //weather->precipitation = (double **)malloc (total_years * sizeof(double *));
//        //weather->RHmax = (double **)malloc (total_years * sizeof(double *));
//        //weather->RHmin = (double **)malloc (total_years * sizeof(double *));
//        weather->solarRadiation =
//            (double **)malloc(total_years * sizeof(double *));
//        weather->tMax = (double **)malloc(total_years * sizeof(double *));
//        weather->tMin = (double **)malloc(total_years * sizeof(double *));
//        weather->vpd = (double **)malloc(total_years * sizeof(double *));
//        //weather->yearlyAmplitude = (double *)malloc (total_years * sizeof(double));
//        //weather->annualAverageTemperature = (double *)malloc (total_years * sizeof(double));
//        weather->lastDoy = (int *)malloc(total_years * sizeof(int));
//        for (y = 0; y < total_years; y++)
//        {
//            weather->wind[y] = (double *)malloc(366 * sizeof(double));
//            //weather->ETref[y] = (double *)malloc (366 * sizeof(double));
//            //weather->precipitation[y] = (double *)malloc (366 * sizeof(double));
//            //weather->RHmax[y] = (double *)malloc (366 * sizeof(double));
//            //weather->RHmin[y] = (double *)malloc (366 * sizeof(double));
//            weather->solarRadiation[y] =
//                (double *)malloc(366 * sizeof(double));
//            weather->tMax[y] = (double *)malloc(366 * sizeof(double));
//            weather->tMin[y] = (double *)malloc(366 * sizeof(double));
//            weather->vpd[y] = (double *)malloc(366 * sizeof(double));
//        }
//
//        /*
//         * Copy crop management to each element
//         */
//        cropmgmt = &elem[i].cropmgmt;
//
//        if (agtbl->op[i] > agtbl->nopfile)
//        {
//            Cycles_printf(VL_ERROR,
//                "ERROR: Operation file for operation index %d "
//                "is not provided.\n", agtbl->op[i]);
//            Cycles_exit(EXIT_FAILURE);
//        }
//
//        opind = agtbl->op[i] - 1;
//
//        cropmgmt->yearsInRotation = agtbl->rotsz[i];
//        cropmgmt->automaticNitrogen = agtbl->auto_N[i];
//        cropmgmt->automaticPhosphorus = agtbl->auto_P[i];
//        cropmgmt->automaticSulfur = agtbl->auto_S[i];
//
//        cropmgmt->rotationYear = 0;
//
//        cropmgmt->totalCropsPerRotation =
//            mgmttbl->cropmgmt[opind].totalCropsPerRotation;
//        cropmgmt->plantingOrder = mgmttbl->cropmgmt[opind].plantingOrder;
//
//        cropmgmt->numTillage = mgmttbl->cropmgmt[opind].numTillage;
//        cropmgmt->Tillage = mgmttbl->cropmgmt[opind].Tillage;
//
//        cropmgmt->numFertilization = mgmttbl->cropmgmt[opind].numFertilization;
//        cropmgmt->FixedFertilization =
//            mgmttbl->cropmgmt[opind].FixedFertilization;
//
//        cropmgmt->numIrrigation = mgmttbl->cropmgmt[opind].numIrrigation;
//        cropmgmt->FixedIrrigation = mgmttbl->cropmgmt[opind].FixedIrrigation;
//
//        cropmgmt->numAutoIrrigation =
//            mgmttbl->cropmgmt[opind].numAutoIrrigation;
//        cropmgmt->autoIrrigation = mgmttbl->cropmgmt[opind].autoIrrigation;
//
//        cropmgmt->op_status[PLANT_OP] =
//            (int *)calloc(cropmgmt->totalCropsPerRotation, sizeof(int));
//        cropmgmt->op_status[TILLAGE_OP] =
//            (int *)calloc(cropmgmt->numTillage, sizeof(int));
//        cropmgmt->op_status[FIXIRR_OP] =
//            (int *)calloc(cropmgmt->numIrrigation, sizeof(int));
//        cropmgmt->op_status[FIXFERT_OP] =
//            (int *)calloc(cropmgmt->numFertilization, sizeof(int));
//
//        for (k = 0; k < elem[i].ps.nsoil; k++)
//        {
//            cropmgmt->tillageFactor[k] = 0.0;
//        }
//
//        cropmgmt->usingAutoIrr = mgmttbl->cropmgmt[opind].usingAutoIrr;
//        cropmgmt->usingAutoFert = mgmttbl->cropmgmt[opind].usingAutoFert;
//
//        /*
//         * Copy crop description to each element
//         */
//        comm = &elem[i].comm;
//
//        comm->NumCrop = croptbl->number;
//        comm->Crop =
//            (crop_struct *)malloc(comm->NumCrop * sizeof(crop_struct));
//
//        for (j = 0; j < comm->NumCrop; j++)
//        {
//            strcpy(comm->Crop[j].cropName, croptbl->cropName[j]);
//
//            comm->Crop[j].userFloweringTT = croptbl->userFloweringTT[j];
//            comm->Crop[j].userMaturityTT = croptbl->userMaturityTT[j];
//            comm->Crop[j].userMaximumSoilCoverage =
//                croptbl->userMaximumSoilCoverage[j];
//            comm->Crop[j].userMaximumRootingDepth =
//                croptbl->userMaximumRootingDepth[j];
//            comm->Crop[j].userExpectedYieldAvg =
//                croptbl->userExpectedYieldAvg[j];
//            comm->Crop[j].userExpectedYieldMax =
//                croptbl->userExpectedYieldMax[j];
//            comm->Crop[j].userExpectedYieldMin =
//                croptbl->userExpectedYieldMin[j];
//            comm->Crop[j].userPercentMoistureInYield =
//                croptbl->userPercentMoistureInYield[j];
//            comm->Crop[j].userFractionResidueStanding =
//                croptbl->userFractionResidueStanding[j];
//            comm->Crop[j].userFractionResidueRemoved =
//                croptbl->userFractionResidueRemoved[j];
//            comm->Crop[j].userClippingBiomassThresholdUpper =
//                croptbl->userClippingBiomassThresholdUpper[j];
//            comm->Crop[j].userClippingBiomassThresholdLower =
//                croptbl->userClippingBiomassThresholdLower[j];
//            comm->Crop[j].userClippingTiming = croptbl->userClippingTiming[j];
//            comm->Crop[j].userClippingDestiny = croptbl->userClippingDestiny[j];
//            comm->Crop[j].userTranspirationMinTemperature =
//                croptbl->userTranspirationMinTemperature[j];
//            comm->Crop[j].userTranspirationThresholdTemperature =
//                croptbl->userTranspirationThresholdTemperature[j];
//            comm->Crop[j].userColdDamageMinTemperature =
//                croptbl->userColdDamageMinTemperature[j];
//            comm->Crop[j].userColdDamageThresholdTemperature =
//                croptbl->userColdDamageThresholdTemperature[j];
//            comm->Crop[j].userTemperatureBase = croptbl->userTemperatureBase[j];
//            comm->Crop[j].userTemperatureOptimum =
//                croptbl->userTemperatureOptimum[j];
//            comm->Crop[j].userTemperatureMaximum =
//                croptbl->userTemperatureMaximum[j];
//            comm->Crop[j].userShootPartitionInitial =
//                croptbl->userShootPartitionInitial[j];
//            comm->Crop[j].userShootPartitionFinal =
//                croptbl->userShootPartitionFinal[j];
//            comm->Crop[j].userRadiationUseEfficiency =
//                croptbl->userRadiationUseEfficiency[j];
//            comm->Crop[j].userTranspirationUseEfficiency =
//                croptbl->userTranspirationUseEfficiency[j];
//            comm->Crop[j].userHIx = croptbl->userHIx[j];
//            comm->Crop[j].userHIo = croptbl->userHIo[j];
//            comm->Crop[j].userHIk = croptbl->userHIk[j];
//            comm->Crop[j].userEmergenceTT = croptbl->userEmergenceTT[j];
//            comm->Crop[j].userNMaxConcentration =
//                croptbl->userNMaxConcentration[j];
//            comm->Crop[j].userNDilutionSlope = croptbl->userNDilutionSlope[j];
//            comm->Crop[j].userKc = croptbl->userKc[j];
//            comm->Crop[j].userAnnual = croptbl->userAnnual[j];
//            comm->Crop[j].userLegume = croptbl->userLegume[j];
//            comm->Crop[j].userC3orC4 = croptbl->userC3orC4[j];
//            comm->Crop[j].userExtinctionCoefficient =
//                croptbl->userExtinctionCoefficient[j];
//            comm->Crop[j].userPlantingDensity = croptbl->userPlantingDensity[j];
//            comm->Crop[j].userClippingStart = croptbl->userClippingStart[j];
//            comm->Crop[j].userClippingEnd = croptbl->userClippingEnd[j];
//            comm->Crop[j].LWP_StressOnset = croptbl->LWP_StressOnset[j];
//            comm->Crop[j].LWP_WiltingPoint = croptbl->LWP_WiltingPoint[j];
//            comm->Crop[j].transpirationMax = croptbl->transpirationMax[j];
//
//            comm->Crop[j].userMaximumSoilCoverage *= 0.94 / 100.0;
//            comm->Crop[j].userPercentMoistureInYield /= 100.0;
//            comm->Crop[j].userFractionResidueStanding /= 100.0;
//            comm->Crop[j].userFractionResidueRemoved /= 100.0;
//            if (comm->Crop[j].userClippingTiming != BADVAL)
//            {
//                comm->Crop[j].userClippingTiming /= 100.0;
//            }
//            else
//            {
//                comm->Crop[j].userClippingTiming = 0.0;
//            }
//
//            comm->Crop[j].calculatedSimAvgYield = 0.0;
//            comm->Crop[j].calculatedSimMaxYield = 0.0;
//            comm->Crop[j].calculatedSimMinYield = 0.0;
//
//            comm->Crop[j].stageGrowth = NO_CROP;
//
//            InitCropSV(&comm->Crop[j]);
//        }
//
//        UpdateCommunity(comm);
//
//        InitializeSoil(&elem[i].soil, soiltbl, &elem[i].ps,
//            &elem[i].cycles_restart, elem[i].attrib.soil_type,
//            ctrl->read_cycles_restart);
//
//        InitializeResidue(&elem[i].residue, elem[i].ps.nsoil);
//
//        InitializeSoilCarbon(&elem[i].soilc, elem[i].ps.nsoil);
//
//        elem[i].comm.NumActiveCrop = 0;
//
//    }
//
//    if (ctrl->read_cycles_restart)
//    {
//        for (i = 0; i < nriver; i++)
//        {
//            riv[i].NO3sol.soluteMass = riv[i].cycles_restart.NO3_Mass;
//            riv[i].NH4sol.soluteMass = riv[i].cycles_restart.NH4_Mass;
//        }
//    }
//    else
//    {
//        for (i = 0; i < nriver; i++)
//        {
//            riv[i].NO3sol.soluteMass = 0.0;
//            riv[i].NH4sol.soluteMass = 0.0;
//
//            for (k = 0; k < elem[riv[i].rightele - 1].ps.nsoil; k++)
//            {
//                riv[i].NO3sol.soluteMass +=
//                    elem[riv[i].rightele - 1].soil.NO3[k];
//                riv[i].NH4sol.soluteMass +=
//                    elem[riv[i].rightele - 1].soil.NH4[k];
//            }
//
//            for (k = 0; k < elem[riv[i].leftele - 1].ps.nsoil; k++)
//            {
//                riv[i].NO3sol.soluteMass +=
//                    elem[riv[i].leftele - 1].soil.NO3[k];
//                riv[i].NH4sol.soluteMass +=
//                    elem[riv[i].leftele - 1].soil.NH4[k];
//            }
//        }
//    }
//}
//
//void WriteCyclesIC(char *restart_fn, elem_struct *elem, river_struct *riv)
//{
//    int             i, j;
//    FILE           *restart_file;
//
//    restart_file = fopen(restart_fn, "wb");
//    CheckFile(restart_file, restart_fn);
//    PIHMprintf(VL_VERBOSE, "Writing Cycles initial conditions.\n");
//
//    for (i = 0; i < nelem; i++)
//    {
//        for (j = 0; j < MAXLYR; j++)
//        {
//            fwrite(&elem[i].soil.SOC_Mass[j], sizeof(double), 1, restart_file);
//        }
//        for (j = 0; j < MAXLYR; j++)
//        {
//            fwrite(&elem[i].soil.SON_Mass[j], sizeof(double), 1, restart_file);
//        }
//        for (j = 0; j < MAXLYR; j++)
//        {
//            fwrite(&elem[i].soil.MBC_Mass[j], sizeof(double), 1, restart_file);
//        }
//        for (j = 0; j < MAXLYR; j++)
//        {
//            fwrite(&elem[i].soil.MBN_Mass[j], sizeof(double), 1, restart_file);
//        }
//        for (j = 0; j < MAXLYR; j++)
//        {
//            fwrite(&elem[i].soil.NO3[j], sizeof(double), 1, restart_file);
//        }
//        for (j = 0; j < MAXLYR; j++)
//        {
//            fwrite(&elem[i].soil.NH4[j], sizeof(double), 1, restart_file);
//        }
//    }
//
//    for (i = 0; i < nriver; i++)
//    {
//        fwrite(&riv[i].NO3sol.soluteMass, sizeof(double), 1, restart_file);
//        fwrite(&riv[i].NH4sol.soluteMass, sizeof(double), 1, restart_file);
//    }
//
//    fclose(restart_file);
//}
