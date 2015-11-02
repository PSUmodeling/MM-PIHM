#ifndef CYCLES_FUNC_HEADER
#define CYCLES_FUNC_HEADER

/* Crop.c */
//void            SelectCropInitialPosition (CropManagementStruct *CropManagement);
//void            SelectNextCrop (CropManagementStruct *CropManagement);
//void            PeekNextCrop (CropManagementStruct *CropManagement);
//void            NewCrop (CommunityStruct *Community, const CropManagementStruct *CropManagement);
void PlantingCrop (CommunityStruct *Community, const CropManagementStruct *CropManagement, int plantingIndex);
void            AddCrop (CropStruct *Crop);
void            KillCrop (CropStruct *Crop);
void            UpdateCommunity (CommunityStruct *Community);

/* CropHarvest.c */
void            GrainHarvest (int y, int doy, int startYear, CropStruct *Crop, ResidueStruct *Residue, const SoilStruct *Soil, SoilCarbonStruct *SoilCarbon, const WeatherStruct *Weather, const char *project);
void            ForageHarvest (int y, int doy, int startYear, CropStruct *Crop, ResidueStruct *Residue, const SoilStruct *Soil, SoilCarbonStruct *SoilCarbon, const WeatherStruct *Weather, const char *project);
void            HarvestCrop (int y, int doy, int startYear, CropStruct *Crop, ResidueStruct *Residue, const SoilStruct *Soil, SoilCarbonStruct *SoilCarbon, const WeatherStruct *Weather, const char *project);
void            DistributeRootDetritus (int y, double rootMass, double rhizoMass, double rootN, double rhizoN, const SoilStruct *Soil, const CropStruct *Crop, ResidueStruct *Residue, SoilCarbonStruct *SoilCarbon);
double          ComputeHarvestIndex (double HIx, double HIo, double HIk, double cumulativeShoot, double cumulativePostFloweringShootBiomass);

/* CropProcess.c */
void            Processes (const int y, const int doy, const int autoNitrogen, CommunityStruct *Community, ResidueStruct *Residue, const WeatherStruct *Weather, SoilStruct *Soil, SoilCarbonStruct *SoilCarbon);
void            CropGrowth (int y, int doy, double *DailyGrowth, const double Stage, CropStruct *Crop, ResidueStruct *Residue, const WeatherStruct *Weather);
void            CropNitrogenDemand (double N_AbgdConcReq, double N_RootConcReq, double *N_ReqAbgdGrowth, double *N_ReqRootGrowth, double *N_ReqRhizodeposition, double *N_CropDemand, CropStruct *Crop);
void            CropNitrogenConcentration (double *N_AbgdConcReq, double *N_RootConcReq, double *NaAbgd, double *NxAbgd, double *NcAbgd, double *NnAbgd, double *NxRoot, const double Stage, const CropStruct *Crop);
void            CropNitrogenStress (const double NaAbgd, const double NcAbgd, const double NnAbgd, CropStruct *Crop);
void CropNitrogenUptake (double *N_ReqAbgdGrowth, double *N_ReqRootGrowth, double *N_ReqRhizodeposition, double *NxAbgd, double *NxRoot, int autoNitrogen, double NO3supply, double NH4supply, double *NO3Uptake, double *NH4Uptake, double *N_CropDemand, CommunityStruct *Community, SoilStruct *Soil);
void            PotentialSoluteUptakeOption2 (double *SoluteSupply, double *SoluteUptake, double Kd, int totalLayers, const double *BD, const double *dz, const double *WaterUptake, const double *Solute, const double *WC);
double          ShootBiomassPartitioning (const double Stage, const double Po, const double Pf);
void            RadiationInterception (int y, int doy, CommunityStruct *Community);
void            Phenology (int y, int doy, const WeatherStruct *Weather, CommunityStruct *Community);
void            ComputeColdDamage (int y, int doy, CropStruct *Crop, const WeatherStruct *Weather, const SnowStruct *Snow, ResidueStruct *Residue);
double          ColdDamage (const double T, const double Crop_Tn, const double Crop_Tth);
double TemperatureFunctionGrowth (double tMax, double tOpt, double tMin, double T);

/* CropThermalTime.c */
void            ComputeThermalTime (int total_years, CommunityStruct *Community, WeatherStruct *Weather);
double          ThermalTime (double T_base, double T_op, double T_Max, double Temperature);

/* CropTranspiration.c */
void            WaterUptake (int y, int doy, CommunityStruct *Community, SoilStruct *Soil, const WeatherStruct *Weather);
void            CalcRootFraction (double *fractionRootsByLayer, SoilStruct *Soil, CropStruct *Crop);
double          TemperatureLimitation (double T, double T_Min, double T_Threshold);

/* DailyOperation.c */
void DailyOperations (int rotationYear, int y, int doy, CropManagementStruct *CropManagement, CommunityStruct *Community, ResidueStruct *Residue, SimControlStruct *SimControl, SnowStruct *Snow, SoilStruct *Soil, SoilCarbonStruct *SoilCarbon, WeatherStruct *Weather, const char *project);
void            GrowingCrop (int rotationYear, int y, int d, FieldOperationStruct *ForcedHarvest, int numHarvest, CommunityStruct *Community, ResidueStruct *Residue, const SimControlStruct *SimControl, SoilStruct *Soil, SoilCarbonStruct *SoilCarbon, const WeatherStruct *Weather, const SnowStruct *Snow, const char *project);
//void            PlantingCrop (int doy, int *nextSeedingYear, int *nextSeedingDate, CropManagementStruct *CropManagement, CommunityStruct *Community);
double          FinalHarvestDate (int lastDoy, int d);
int             ForcedMaturity (int rotationYear, int d, int lastDoy, int nextSeedingYear, int nextSeedingDate, int rotationSize);

/* Fertilization.c */
void            ApplyFertilizer (FieldOperationStruct *fixedFertilization, SoilStruct *Soil, ResidueStruct *Residue);

/* FieldOperations.c */
//void            SelectNextOperation (int NumOperation, int *operationIndex);
//void            SelectOperationYear (int rotationYear, const FieldOperationStruct *FieldOperation, int NumOperation, int *operationIndex);
//int             IsOperationToday (int rotationYear, int doy, const FieldOperationStruct *FieldOperation, int operationIndex);
int IsOperationToday (int rotationYear, int doy, FieldOperationStruct *FieldOperation, int numOperation, int *operationIndex);
void UpdateOperationStatus (FieldOperationStruct *FieldOperation, int numOperation);

/* Initialize.c */
void            Initialize (SimControlStruct *SimControl, WeatherStruct *Weather, SoilStruct *Soil, ResidueStruct *Residue, SoilCarbonStruct *SoilCarbon, CommunityStruct *Community, CropManagementStruct *CropManagement, SnowStruct *Snow);
void            FreeCyclesStruct (CyclesStruct Cycles, int total_years);

/* Irrigation.c */
double          FindIrrigationVolume (int opLayer, double opWaterDepletion, const SoilStruct *Soil);

/* MiscFunc.c */
int             IsLeapYear (int year);
int             doy (int year, int month, int mday, int leap_year_mode);
int             t2doy (time_t * rawtime);
void            doy2date (int year, int jday, int *month, int *mday, int leap_year_mode);
int             LT (double x, double y);
int             LE (double x, double y);
int             EQ (double x, double y);
int             GE (double x, double y);
int             GT (double x, double y);

/* Print.c */
void            InitializeOutput (char *project, const CommunityStruct *Community, int layers);
void            PrintDailyOutput (int y, int doy, int start_year, const WeatherStruct *Weather, const CommunityStruct *Community, const SoilStruct *Soil, const SnowStruct *Snow, const ResidueStruct *Residue, const char *project);
void            PrintSeasonOutput (int y, int doy, int start_year, const WeatherStruct *Weather, const CropStruct *Crop, const char *project);
void            PrintAnnualOutput (int y, int start_year, const SoilStruct *Soil, const SoilCarbonStruct *SoilCarbon, const char *project);
void            PrintCarbonEvolution (int y, int start_year, int total_layers, const SoilStruct *Soil, const SoilCarbonStruct *SoilCarbon, const ResidueStruct *Residue, char *project);
void            StoreSummary (SummaryStruct *Summary, const SoilCarbonStruct *SoilCarbon, const ResidueStruct *Residue, int totalLayers, int y);
void            PrintSummary (const SummaryStruct *Summary, int totalYears, const char *project);

/* ReadCrop.c */
void            ReadCrop (char *project, CommunityStruct *Community);

/* ReadFunc.c */
char FirstNonWhite (char *cmdstr);
int Readable (char *cmdstr);
int FindLine (FILE * fid, char *token);
void NextLine (FILE * fid, char *cmdstr);
int CountLine (FILE * fid, int num_arg, ...);
void CheckFile (FILE * fid, char *fn);
void ReadTS (char *cmdstr, int *ftime, double *data, int nvrbl);
int CountOccurance (FILE *fid, char *token);
void ReadKeywordDouble (char *buffer, char *keyword, double *value);
void ReadKeywordInt (char *buffer, char *keyword, int *value);
void ReadKeywordTime (char *buffer, char *keyword, int *value);
void ReadKeywordStr (char *buffer, char *keyword, char *value);

/* ReadOperation.c */
void            ReadOperation (char *project, CropManagementStruct *CropManagement, const CommunityStruct *Community, int yearsInRotation);

/* ReadSimControl.c */
void            ReadSimControl (char *project, SimControlStruct *SimControl);

/* ReadSoil.c */
void            ReadSoil (char *project, SoilStruct *Soil);

/* ReadWeather.c */
void            ReadWeather (char *project, WeatherStruct *Weather, int start_year, int total_years);

/* ReferenceET.c */
double          CalculatePMET (double lat, double pAtm, double screeningHeight, double Tmax, double Tmin, double sRad, double rhMax, double rhMin, double wind, double doy);
double          SaturatedVaporPressure (double T);
double          PotentialRadiation (double Lat, int doy);
double          NetRadiation (double Pot_Rad, double Solar_Rad, double Actual_VP, double TMax, double TMin);
double          Aero_Res (double uz, double z);

/* Residue.c */
void            InitializeResidue (ResidueStruct *Residue, int totalYears, int totalLayers);
void            ComputeResidueCover (ResidueStruct *Residue);
void            ResidueWetting (ResidueStruct *Residue, SoilStruct *Soil);
void            ResidueEvaporation (ResidueStruct *Residue, SoilStruct *Soil, const CommunityStruct *Community, double ETo, double snowCover);

/* Snow.c */
void            SnowProcesses (SnowStruct *Snow, int y, int doy, WeatherStruct *Weather, double TauStandingRes, double CropInterception);
void            CalculateSnowFall (SnowStruct *Snow, double Tavg, double PP);
void            CalculateSnowMelt (SnowStruct *Snow, double Tavg, double Tx, double Tn);
void            CalculateSnowEvaporation (SnowStruct *Snow, double TauStandingRes, double CropInterception, double ETo);
double          CalculateSnowCover (SnowStruct *Snow);

/* Soil.c */
void            InitializeSoil (SoilStruct *Soil, WeatherStruct *Weather, SimControlStruct *SimControl);
double          SoilWaterPotential (double SaturationWC, double AirEntryPot, double Campbell_b, double WC);
double          VolumetricWCAt33Jkg (double Clay, double Sand, double OM);
double          VolumetricWCAt1500Jkg (double Clay, double Sand, double OM);
double          SoilWaterContent (double SaturationWC, double AirEntryPot, double Campbell_b, double Water_Potential);
double          BulkDensity (double Clay, double Sand, double OM);

/* SoilCarbon.c */
void            InitializeSoilCarbon (SoilCarbonStruct *SoilCarbon, int totalLayers);
void            ComputeFactorComposite (SoilCarbonStruct *SoilCarbon, int doy, int y, int last_doy, SoilStruct *Soil);
void            ComputeSoilCarbonBalanceMB (SoilCarbonStruct *SoilCarbon, int y, ResidueStruct *Residue, SoilStruct *Soil, double *tillageFactor);
void            ComputeSoilCarbonBalance (SoilCarbonStruct *SoilCarbon, int y, ResidueStruct *Residue, SoilStruct *Soil, double *tillageFactor);
void            StoreOutput (SoilCarbonStruct *SoilCarbon, int y, int totalLayers, double *SOCMass);
double          Aeration (double AC);
double          Moisture (double wp);
double          TemperatureFunction (double T);
double          MaximumAbgdHumificationFactor (double clayFraction);
double          MaximumRootHumificationFactor (double clayFraction);
double          MaximumRhizHumificationFactor (double clayFraction);
double          MaximumManuHumificationFactor (double clayFraction);
double          NitrogenMineralization (double CNDecomposing, double CNnew, double humRate, double decomposedMass);
double          CNdestiny (double NmineralConc, double CNdecomposing);
double          PoolNitrogenMineralization (double NmineralConc, double CNRatioDecomposing, double humRate, double decomposedMass, double carbonConc);
double          Function_CNnew (double NmineralConc, double CNDecomposingPool);

/* SoilEvaporation.c */
void            Evaporation (SoilStruct *Soil, const CommunityStruct *Community, ResidueStruct *Residue, double ETo, double SnowCover);
double          Depth_Limitation_To_Evaporation (double Depth);
double          Water_Content_Limitation_To_Evaporation (double FC, double WC_AirDry, double WC);

/* SoilInfiltration.c */
void            Redistribution (int y, int doy, double precipitation, double snowFall, double snowMelt, int hourlyInfiltration, const CommunityStruct *Community, SoilStruct *Soil, ResidueStruct *Residue);
void            CascadeRedistribution (SoilStruct *Soil);
void            SubDailyRedistribution (SoilStruct *Soil);
void            CurveNumber (SoilStruct *Soil);
double          K_Sat (double WCATSAT, double WCAT33, double b);
double          AverageHydraulicConductance (double WCSAT, double WCFC, double aep, double b, double SWP1, double SWP2, double ks);
double          Numerator (double aep, double b, double SWP);
double          Denominator (double aep, double b, double SWP);

/* SoilNitrogen.c */
void            NitrogenTransformation (int y, int doy, SoilStruct *Soil, const CommunityStruct *Community, const ResidueStruct *Residue, const WeatherStruct *Weather, const SoilCarbonStruct *SoilCarbon);
void            Nitrification (double *Profile_N_Nitrified, double *Profile_N2O_Nitrified, SoilStruct *Soil, const SoilCarbonStruct *SoilCarbon);
void            Denitrification (double *Profile_N_Denitrified, double *Profile_N2O_Denitrified, SoilStruct *Soil, const SoilCarbonStruct *SoilCarbon);
void            Volatilization (int y, int doy, double *Profile_NH4_Volatilization, SoilStruct *Soil, const CommunityStruct *Community, const ResidueStruct *Residue, const WeatherStruct *Weather);
double          N2OFractionNitrification (double air);
double          pHFunction (double pH);
double          VolatilizationDepthFunction (double depth);
double          AirMolarDensity (double T, double P);
double          BoundaryLayerConductance (double RI, double RM, double WS, double AMD);

/* SoilSolute.c */
void            SoluteTransport (int totalLayers, double Sol_Kd, double WInConc, double *leachate, const double *WFlux, double *soluteMass, const double *BD, const double *thickness, const double *porosity, const double *WCinitial);
void            SoluteTransportEvaporation (int totalLayers, double Sol_Kd, const double *WFlux, double *soluteMass, const double *BD, const double *thickness, const double *porosity, double *WCinitial);
double          LinearEquilibriumConcentration (double Kd, double bulkDensity, double layerThickness, double waterContent, double soluteMass);
double          LinearEquilibriumSoluteMass (double Kd, double bulkDensity, double layerThickness, double waterContent, double concentration);

/* SoilTemperature.c */
void            Temperature (int y, int doy, double snowCover, double cropInterception, SoilStruct *Soil, WeatherStruct *Weather, ResidueStruct *Residue);
double          HeatCapacity (double bulkDensity, double volumetricWC);
double          HeatConductivity (double bulkDensity, double volumetricWC, double fractionClay);
double          EstimatedSoilTemperature (double nodeDepth, int doy, double annualAvgTemperature, double yearlyAmplitude, int phase, double dampingDepth);

/* Tillage.c */
void            ExecuteTillage (double *abgdBiomassInput, const FieldOperationStruct *Tillage, double *tillageFactor, SoilStruct *Soil, ResidueStruct *Residue);
void            TillageFactorSettling (double *tillageFactor, int totalLayers, const double *waterContent, const double *Porosity);
double          Fraction (double a, double b, double c, double d, double f);
void            ComputeTillageFactor (const FieldOperationStruct *Tillage, double *tillageFactor, const SoilStruct *Soil, const double *soilLayerBottom, double toolDepth);
double          ComputeTextureFactor (double Clay);

/* Weather.c */
double          SatVP (double T);
void            CalculateDerivedWeather (WeatherStruct *Weather, int total_years);

void            PrintSimContrl (SimControlStruct SimControl);
void            PrintSoil (SoilStruct Soil);
void            PrintCrop (CommunityStruct Community);
void            PrintOperation (FieldOperationStruct *plantedCrops, int NumPlanting, FieldOperationStruct *ForcedHarvest, int NumHarvest, FieldOperationStruct *Tillage, int NumTillage, FieldOperationStruct *FixedIrrigation, int NumIrrigation, FieldOperationStruct *FixedFertilization, int NumFertilization);
void            PrintWeather (WeatherStruct Weather);

#endif
