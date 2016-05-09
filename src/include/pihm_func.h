#ifndef PIHM_FUNC_HEADER
#define PIHM_FUNC_HEADER

#define _ARITH_

#define SURF(i)         i
#define UNSAT(i)        i + pihm->numele
#define GW(i)           i + 2 * pihm->numele
#define RIVSTG(i)       i + 3 * pihm->numele
#define RIVGW(i)        i + 3 * pihm->numele + pihm->numriv

/*
 * Function Declarations
 */
#ifdef _ENKF_
void            PIHMRun (char *, char *, int, int, int, int, double *);
#else
void            PIHMRun (char *, char *, int);
#endif
void            CreateOutputDir (char *, int);
void            ReadAlloc (char *, pihm_struct);
void            ReadRiv (char *, rivtbl_struct *, shptbl_struct *,
    matltbl_struct *, forc_struct *);
void            ReadMesh (char *, meshtbl_struct *);
void            ReadAtt (char *, atttbl_struct *, int);
void            ReadSoil (char *, soiltbl_struct *);
void            ReadGeol (char *, geoltbl_struct *);
void            ReadLC (char *, lctbl_struct *);
void            ReadForc (char *, forc_struct *);
void            ReadLAI (char *, forc_struct *, int numele,
    const atttbl_struct *);
void            ReadBC (char *, forc_struct *);
void            ReadPara (char *, ctrl_struct *);
void            ReadCalib (char *, calib_struct *);
void            ReadIC (char *, elem_struct *, int, river_struct *, int);
void            FreeData (pihm_struct pihm);
#ifdef _NOAH_
void            ReadLsm (char *, double *, double *, ctrl_struct *,
    noahtbl_struct *);
void            ReadRad (char *, forc_struct *);
#endif

int             Readable (char *token);
int             FindLine (FILE * fid, char *token);
void            NextLine (FILE * fid, char *cmdstr);
int             CountLine (FILE *, char *, int, ...);
void            ReadTS (char *cmdstr, int *ftime, double *data, int nvrbl);
void            CheckFile (FILE * fid, char *fn);
void            ReadKeywordDouble (char *buffer, char *keyword,
    double *value);
void            ReadKeywordInt (char *buffer, char *keyword, int *value);
void            ReadKeywordTime (char *buffer, char *keyword, int *value);
void            ReadKeywordStr (char *buffer, char *keyword, char *value);
int             CountOccurance (FILE * fid, char *token);

void            Initialize (pihm_struct, N_Vector);
void            InitMeshStruct (elem_struct *, int, meshtbl_struct);
void            InitTopo (elem_struct *, int, meshtbl_struct);
void            InitSoil (elem_struct *, int, atttbl_struct, soiltbl_struct,
#ifdef _NOAH_
    noahtbl_struct,
#endif
    calib_struct);
void            ZeroWaterFlux (wf_struct *wf);
double          FieldCapacity (double, double, double, double, double);
double          WiltingPoint (double, double, double, double);
void            InitLC (elem_struct *, int, atttbl_struct, lctbl_struct,
    calib_struct);
void            InitRiver (river_struct *, int, elem_struct *, rivtbl_struct,
    shptbl_struct, matltbl_struct, meshtbl_struct, calib_struct);
void            InitForcing (elem_struct *, int, river_struct *, int,
    atttbl_struct, rivtbl_struct, forc_struct *forc, calib_struct);
void            CorrectElevation (elem_struct *, int, river_struct *, int);
void            InitSurfL (elem_struct *, int, river_struct *,
    meshtbl_struct);
void            SaturationIC (elem_struct *, int, river_struct *, int);
void            InitVar (elem_struct *, int, river_struct *, int, N_Vector);
void            CalcModelStep (ctrl_struct *);

void            MapOutput (char *, pihm_struct, char *);
void            InitOutputFile (prtctrl_struct *, int, int);
void            ApplyForcing (forc_struct *, int);
void            IntcpSnowET (int, double, pihm_struct);
void            IntrplForcing (tsdata_struct, int, int);
double          MonthlyLAI (int, int);
double          MonthlyRL (int, int);
double          MonthlyMF (int);
int             Hydrol (realtype, N_Vector, N_Vector, void *);
void            LateralFlow (pihm_struct);
void            VerticalFlow (pihm_struct);
void            RiverFlow (pihm_struct);
void            RiverToEle (river_struct *, elem_struct *, elem_struct *,
    int, double *, double *, double *, double);
double          DhByDl (double *, double *, double *);
double          RivArea (int, double, double);
double          RivPerim (int, double, double);
double          EqWid (int, double, double);
double          OLFEleToRiv (double, double, double, double, double, double);
double          OverlandFlow (double, double, double, double, double);
double          AvgY (double, double, double);
double          AvgYsfc (double, double, double);
double          EffKinf (double, double, int, double, double, double);
double          EffKV (double, double, int, double, double, double);
double          AvgKV (double dmac, double deficit, double gw,
    double macp_status, double satn, double satkfunc, double kmacv,
    double ksatv, double areafh);
double          EffKH (double, double, double, double, double, double);

void            PrtInit (pihm_struct pihm, char *simulation);
void            PrintData (prtctrl_struct *, int, int, int, int, int);
int             MacroporeStatus (double, double, double, double, double,
    double, double);
double          KrFunc (double, double, double);
double          Psi (double, double, double);

void            Summary (pihm_struct pihm, N_Vector CV_Y, double stepsize);
void            SetCVodeParam (pihm_struct, void *, N_Vector);
void            SolveCVode (int *, int, int, void *, N_Vector);
int             SoilTex (double, double);
double          Qtz (int);
double          PtfKV (double, double, double, double, int);
double          PtfThetaS (double, double, double, double, int);
double          PtfThetaR (double, double, double, double, int);
double          PtfAlpha (double, double, double, double, int);
double          PtfBeta (double, double, double, double, int);
//
//#ifdef _DAILY_
//void InitDailyStruct (pihm_struct pihm);
//#endif
//
void            BKInput (char *, char *);
//void PihmFree (void **ptr);
void            PihmExit (int error);

#ifdef _NOAH_
void            InitLsm (elem_struct *, int, ctrl_struct, noahtbl_struct,
    calib_struct);
void            CalcLatFlx (const ws_struct *ws, const ps_struct *ps,
    wf_struct *wf);
int             FindWT (const double *sldpth, int nsoil, double gw,
    double *satdpth);
void            DefSldpth (double *, int *, double, const double *, int);
void            RootDist (const double *, int, int, double *);
void            CalcSlopeAspect (elem_struct *, int, meshtbl_struct);
int             FindLayer (const double *, int, double);
double          AvgElev (elem_struct *, int);
double          GWTransp (double, double *, int, int);
void            SunPos (int, double, double, double, double, double *,
    double *);
double          TopoRadn (double, double, double, double, double, double,
    const double *, double);
void            CalHum (ps_struct *, es_struct *);
void            Noah (int, pihm_struct);
double          FrozRain (double, double);
void            AvgFlux (elem_struct *, int, int);
void            SfcDifOff (ps_struct *, const lc_struct *, double, double,
    int);
void            SFlx (ws_struct *, wf_struct *, const wf_struct *,
                    es_struct *, ef_struct *, ps_struct *, lc_struct *,
                    soil_struct *,
#ifdef _CYCLES_
                    comm_struct *, residue_struct *,
#endif
                    int);
double          CSnow (double);
void            SnowNew (const es_struct *, double, ps_struct *);
double          SnFrac (double, double, double, double);
void            AlCalc (ps_struct *, double, int);
double          TDfCnd (double, double, double, double, double);
double          Snowz0 (double, double, double);
void            Penman (wf_struct *, es_struct *, ef_struct *, ps_struct *,
    double *, double, int, int);
void            CanRes (ws_struct *, es_struct *, ef_struct *, ps_struct *,
    const double *, const soil_struct *, const lc_struct *);
void            DEvap (const ws_struct *, wf_struct *, const ps_struct *,
    const lc_struct *, const soil_struct *);
void            Evapo (ws_struct *, wf_struct *, ps_struct *,
    const lc_struct *, soil_struct *,
#ifdef _CYCLES_
    comm_struct *, residue_struct *, const es_struct *es,
#endif
    const double *, double);
void            Transp (const ws_struct *, wf_struct *, const ps_struct *,
    const lc_struct *, const soil_struct *, const double *);
void            NoPac (ws_struct *, wf_struct *, const wf_struct *,
    es_struct *, ef_struct *, ps_struct *, lc_struct *,
    soil_struct *,
#ifdef _CYCLES_
    comm_struct *, residue_struct *,
#endif
    const double *, double, double);
double          TBnd (double, double, const double *, double, int, int);
double          TmpAvg (double, double, double, const double *, int, int);
void            SnkSrc (double *, double, double, double *,
    const soil_struct *, const double *, int, double, int, double);
void            Rosr12 (double *, double *, double *, double *, double *,
    double *, int);
void            ShFlx (ws_struct *, es_struct *, ef_struct *, ps_struct *,
    const lc_struct *, const soil_struct *, double, double,
    double, const double *, double);
void            SmFlx (ws_struct *, wf_struct *, const wf_struct *,
    ps_struct *, const lc_struct *, const soil_struct *,
#ifdef _CYCLES_
    residue_struct *,
#endif
    const double *, double, double);
void            HRT (ws_struct *, es_struct *, ef_struct *, ps_struct *,
    const lc_struct *, const soil_struct *, double *,
    const double *, double, double, double, double, double *,
    double *, double *);
void            SRT (ws_struct *, wf_struct *, const wf_struct *, ps_struct *,
    const soil_struct *,
#ifdef _CYCLES_
    residue_struct *,
#endif
    double *, double *, double *, double *, double *, const double *, double);
void            SStep (ws_struct *, wf_struct *, const wf_struct *,
                    ps_struct *, const soil_struct *, double *, double,
                    const double *, double *, double *, double *, double *,
                    double);
void            WDfCnd (double *, double *, double, double, double, int,
    const soil_struct *, const ps_struct *);
void            SnoPac (ws_struct *, wf_struct *, const wf_struct *,
    es_struct *, ef_struct *, ps_struct *, lc_struct *,
    soil_struct *,
#ifdef _CYCLES_
    comm_struct *, residue_struct *,
#endif
    int, const double *, double, double, double, double);
void            SnowPack (double, double, double *, double *, double, double);
double          EFFKV (double, double, int, double, double, double);
double          Pslmu (double);
double          Pslms (double);
double          Pslhu (double);
double          Pslhs (double);
double          Pspmu (double);
double          Pspms (double);
double          Psphu (double);
double          Psphs (double);
#ifdef _DAILY_
void            DailyVar (int, int, pihm_struct);
void            InitDailyStruct (pihm_struct);
#endif
#endif

#ifdef _ENKF_
void            EnKFRead (char *project, enkf_struct ens);
double          Randn ();
void            MapVar (var_struct * var, int numele, int numriv);
void            Perturb (char *project, enkf_struct ens, char *outputdir);
void            MapVar (var_struct * var, int numele, int numriv);
void            Calib2Mbr (calib_struct cal, double *param);
void            Mbr2Cal (calib_struct *cal, const double *param);
void            WriteParamOutput (int rawtime, enkf_struct ens, int ind,
    char *outputdir);
void            WriteCalFile (enkf_struct ens, char *project);
void            JobHandout (int starttime, int endtime, int startmode,
    ensmbr_struct * member, double *param, int ne, int total_jobs);
void            JobRecv (int *starttime, int *endtime, int *startmode,
    double *param, int ne);
void            PrintEnKFStatus (int starttime, int endtime);
void            JobHandIn (int total_jobs);
void            WritePara (char *project, int start_mode, int start_time,
    int end_time);

void            EnKFCore (double *xa, double obs, double obs_err, double *xf,
    int ne);
void            EnKF (enkf_struct ens, int obs_time, char *outputdir);
void            ReadObs (int obs_time, char *fn, double *obs,
    double *obs_error);
void            InitOper (pihm_struct pihm, enkf_struct ens);
void            DisOper (obs_struct * obs, var_struct *var, pihm_struct pihm);
void            ReadFcst (enkf_struct enkf, obs_struct obs, double *xf);
void            ReadVar (char *project, char *outputdir, enkf_struct ens,
    int obs_time);
void            UpdAnlys (enkf_struct ens, double obs, double obs_error,
    double *xf);
void            CovInflt (enkf_struct ens, enkf_struct ens0);
void            WriteEnKFOut (char *project, enkf_struct ens, char *outputdir,
    int t);
void            GenRandNum (int ne, int nparam, double **randnum,
    double lower, double upper);
double          Randn ();
void            LandSfcTmpOper (obs_struct * obs, var_struct *var, pihm_struct pihm);
void		COSMOSOper (obs_struct *obs, var_struct *var, pihm_struct pihm);
void            FreeEns (enkf_struct ens);
int             FindVar (var_struct *var, char *varname);
void InitEns (enkf_struct ens);
#endif

#ifdef _CYCLES_
void DailyCycles (int t, pihm_struct pihm);
void FirstDOY (int *rotationYear, int yearsInRotation, int totalLayers, soilc_struct *SoilCarbon, residue_struct *Residue, const soil_struct *Soil);
void            ReadCyclesCtrl (char *, agtbl_struct *, int);
void            ReadSoilInit (char *, soiltbl_struct *);
void            ReadCrop (char *, croptbl_struct *);
void            ReadOperation (const agtbl_struct *, mgmttbl_struct *, const croptbl_struct *);
int             CropExist (char *, const croptbl_struct *);
void            InitCycles (elem_struct *, int, const ctrl_struct *, const mgmttbl_struct *, const agtbl_struct *, const croptbl_struct *, const soiltbl_struct *);
void            InitializeSoil (soil_struct *, const soiltbl_struct *, const ps_struct *);
double          BulkDensity (double, double, double);
void            InitializeResidue (residue_struct *, int);
void            InitializeSoilCarbon (soilc_struct *, int);
void            ComputeFactorComposite (soilc_struct *SoilCarbon, int doy, int y, int last_doy, soil_struct *Soil);
void            ComputeSoilCarbonBalanceMB (soilc_struct *SoilCarbon, int y, residue_struct *Residue, soil_struct *Soil, double *tillageFactor);
void            ComputeSoilCarbonBalance (soilc_struct *SoilCarbon, int y, residue_struct *Residue, soil_struct *Soil, double *tillageFactor);
void            StoreOutput (soilc_struct *SoilCarbon, int y, int totalLayers, double *SOCMass);
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
void WaterUptake (comm_struct *Community, soil_struct *Soil, double sfctmp, wf_struct *wf, double pc, double dt);
double TemperatureLimitation (double T, double T_Min, double T_Threshold);
void CalcRootFraction (double *fractionRootsByLayer, soil_struct *Soil, crop_struct *Crop);
int DOY (int);
int IsLeapYear (int year);
void DailyOperations (int y, int d, cropmgmt_struct *CropManagement, comm_struct *Community, residue_struct *Residue, ctrl_struct *SimControl, snow_struct *snow, soil_struct *Soil, soilc_struct *SoilCarbon, weather_struct *Weather);
double Depth_Limitation_To_Evaporation (double Depth);
double Water_Content_Limitation_To_Evaporation (double FC, double WC_AirDry, double WC);
void Evaporation (soil_struct *Soil, const comm_struct *Community, residue_struct *Residue, double ETo, double SnowCover);
void LastDOY (int y, int totalLayers, soil_struct *Soil, soilc_struct *SoilCarbon, residue_struct *Residue);
void GrowingCrop (int y, int doy, comm_struct *Community, residue_struct *Residue, const ctrl_struct *SimControl, soil_struct *Soil, soilc_struct *SoilCarbon, cropmgmt_struct *, const weather_struct *Weather, const snow_struct *Snow);
void CropStage (int d, comm_struct *Community, int last_doy);
double FinalHarvestDate (int lastDoy, int d);
void Phenology (int y, int doy, const weather_struct *Weather, comm_struct *Community);
double ThermalTime (double T_base, double T_op, double T_Max, double Temperature);
void RadiationInterception (int y, int doy, comm_struct *Community);
void Processes (int y, int doy, int autoNitrogen, comm_struct *Community, residue_struct *Residue, const weather_struct *Weather, soil_struct *Soil, soilc_struct *SoilCarbon);
void CropNitrogenConcentration (double *N_AbgdConcReq, double *N_RootConcReq, double *NaAbgd, double *NxAbgd, double *NcAbgd, double *NnAbgd, double *NxRoot, double Stage, const crop_struct *Crop);
void CropNitrogenStress (double NaAbgd, double NcAbgd, double NnAbgd, crop_struct *Crop);
void CropGrowth (int y, int doy, double *DailyGrowth, double Stage, crop_struct *Crop, residue_struct *Residue, const weather_struct *Weather);
void CropNitrogenDemand (double N_AbgdConcReq, double N_RootConcReq, double *N_ReqAbgdGrowth, double *N_ReqRootGrowth, double *N_ReqRhizodeposition, double *N_CropDemand, crop_struct *Crop);
void PotentialSoluteUptakeOption2 (double *SoluteSupply, double *SoluteUptake, double Kd, int totalLayers, const double *BD, const double *dz, const double *WaterUptake, const double *Solute, const double *WC);
void CropNitrogenUptake (double *N_ReqAbgdGrowth, double *N_ReqRootGrowth, double *N_ReqRhizodeposition, double *NxAbgd, double *NxRoot, int autoNitrogen, double NO3supply, double NH4supply, double *NO3Uptake, double *NH4Uptake, double *N_CropDemand, comm_struct *Community, soil_struct *Soil);
void DistributeRootDetritus (double rootMass, double rhizoMass, double rootN, double rhizoN, const soil_struct *Soil, const crop_struct *Crop, residue_struct *Residue, soilc_struct *SoilCarbon);
double ShootBiomassPartitioning (double Stage, double Po, double Pf);
double TemperatureFunctionGrowth (double tMax, double tOpt, double tMin, double T);
int ForcedClipping (int d, comm_struct *Community);
void GrainHarvest (int y, int doy, crop_struct *Crop, residue_struct *Residue, soil_struct *Soil, soilc_struct *SoilCarbon);
void ComputeColdDamage (int y, int doy, crop_struct *Crop, const weather_struct *Weather, const snow_struct *Snow, residue_struct *Residue);
double ColdDamage (double T, double Crop_Tn, double Crop_Tth);
void ForageAndSeedHarvest (int y, int doy, crop_struct *Crop, residue_struct *Residue, soil_struct *Soil, soilc_struct *SoilCarbon);
void HarvestCrop (int y, int doy, crop_struct *Crop, residue_struct *Residue, soil_struct *Soil, soilc_struct *SoilCarbon);
void PlantingCrop (comm_struct *Community, const cropmgmt_struct *CropManagement, int plantingIndex);
void AddCrop (crop_struct *Crop);
void KillCrop (crop_struct *Crop);
void UpdateCommunity (comm_struct *Community);
double ComputeHarvestIndex (double HIx, double HIo, double HIk, double cumulativeShoot, double cumulativePostFloweringShootBiomass);
int IsOperationToday (int rotationYear, int doy, op_struct *FieldOperation, int numOperation, int *operationIndex);
void ApplyFertilizer (op_struct *fixedFertilization, soil_struct *Soil, residue_struct *Residue);
void UpdateOperationStatus (op_struct *FieldOperation, int numOperation);
void FieldOperation (int rotationYear, int y, int doy, cropmgmt_struct *CropManagement, comm_struct *Community, soil_struct *Soil, residue_struct *Residue, ctrl_struct *SimControl, soilc_struct *SoilCarbon, weather_struct *Weather);
void ExecuteTillage (double *abgdBiomassInput, const op_struct *Tillage, double *tillageFactor, soil_struct *Soil, residue_struct *Residue);
void TillageFactorSettling (double *tillageFactor, int totalLayers, const double *waterContent, const double *Porosity);
double Fraction (double a, double b, double c, double d, double f);
void ComputeTillageFactor (const op_struct *Tillage, double *tillageFactor, const soil_struct *Soil, const double *soilLayerBottom, double toolDepth);
double ComputeTextureFactor (double Clay);
void ComputeResidueCover (residue_struct *Residue);
void ResidueEvaporation (residue_struct *Residue, soil_struct *Soil, const comm_struct *Community, double ETo, double snowCover);
void NitrogenTransformation (int y, int doy, soil_struct *Soil, const comm_struct *Community, const residue_struct *Residue, const weather_struct *Weather, const soilc_struct *SoilCarbon);
void Nitrification (double *Profile_N_Nitrified, double *Profile_N2O_Nitrified, soil_struct *Soil, const soilc_struct *SoilCarbon);
void Denitrification (double *Profile_N_Denitrified, double *Profile_N2O_Denitrified, soil_struct *Soil, const soilc_struct *SoilCarbon);
void Volatilization (int y, int doy, double *Profile_NH4_Volatilization, soil_struct *Soil, const comm_struct *Community, const residue_struct *Residue, const weather_struct *Weather);
double N2OFractionNitrification (double air);
double pHFunction (double pH);
double VolatilizationDepthFunction (double depth);
double AirMolarDensity (double T, double P);
double BoundaryLayerConductance (double RI, double RM, double WS, double AMD);
void ResidueWetting (residue_struct *Residue, double *infil_vol);
double FindIrrigationVolume (int opLayer, double opWaterDepletion, const soil_struct *Soil);
void SoluteTransport (elem_struct *elem, int numele, double dt);
void Adsorption (const double *sldpth, const double *sh2o, const double *bd, int nsoil, double Sol_Kd, solute_struct *solute);
double LinearEquilibriumConcentration (double Kd, double bulkDensity, double layerThickness, double waterContent, double soluteMass);
double LinearEquilibriumSoluteMass (double Kd, double bulkDensity, double layerThickness, double waterContent, double concentration);
void Elem2ElemSolTrnsp (const elem_struct *src, const elem_struct *snk, double fluxsub, const double *conc, double dt, double *flux_sol_src, double *flux_sol_snk);
#endif

#endif
