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
void            PIHMRun (char *, char *, int
#ifdef _ENKF_
    , int, int, int, double *
#endif
    );
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
void            ReadLAI (char *, forc_struct *, int, const atttbl_struct *);
void            ReadBC (char *, forc_struct *);
void            ReadPara (char *, ctrl_struct *);
void            ReadCalib (char *, calib_struct *);
void            ReadIC (char *, elem_struct *, int, river_struct *, int);
void            FreeData (pihm_struct);
#ifdef _NOAH_
void            ReadLsm (char *, double *, double *, ctrl_struct *,
    noahtbl_struct *);
void            ReadRad (char *, forc_struct *);
#endif

int             Readable (char *);
int             FindLine (FILE *, char *);
void            NextLine (FILE *, char *);
int             CountLine (FILE *, char *, int, ...);
int             ReadTS (char *, int *, double *, int);
int             ReadKeyword (char *, char *, void *, char);
int             CountOccurance (FILE *, char *);

void            Initialize (pihm_struct, N_Vector);
void            InitMeshStruct (elem_struct *, int, meshtbl_struct);
void            InitTopo (elem_struct *, int, meshtbl_struct);
void            InitSoil (elem_struct *, int, soiltbl_struct,
#ifdef _NOAH_
    noahtbl_struct,
#endif
    calib_struct);
void            InitWState (wstate_struct *);
void            InitRiverWState (river_wstate_struct *);
void            InitPState (pstate_struct *);
void            InitWFlux (wflux_struct *);
void            InitRiverWFlux (river_wflux_struct *);
void            InitEState (estate_struct *);
void            InitEFlux (eflux_struct *);
double          FieldCapacity (double, double, double, double, double);
double          WiltingPoint (double, double, double, double);
void            InitLC (elem_struct *, int, lctbl_struct, calib_struct);
void            InitRiver (river_struct *, int, elem_struct *, rivtbl_struct,
    shptbl_struct, matltbl_struct, meshtbl_struct, calib_struct);
void            InitForcing (elem_struct *, int, river_struct *, int,
    atttbl_struct, rivtbl_struct, forc_struct *, calib_struct);
void            CorrectElevation (elem_struct *, int, river_struct *, int);
void            InitSurfL (elem_struct *, int, river_struct *,
    meshtbl_struct);
void            SaturationIC (elem_struct *, int, river_struct *, int);
void            InitVar (elem_struct *, int, river_struct *, int, N_Vector);
void            CalcModelStep (ctrl_struct *);

void            MapOutput (char *, pihm_struct, char *);
void            InitOutputFile (prtctrl_struct *, int, int);
void            ApplyForcing (forc_struct *, elem_struct *, int, river_struct *, int, int
#ifdef _BGC_
    , ctrl_struct *
#endif
    );
void            ApplyBC (forc_struct *, elem_struct *, int, int);
void            ApplyMeteoForc (forc_struct *, elem_struct *, int, int);
void            ApplyLAI (forc_struct *, elem_struct *, int, int);
void            ApplyRiverBC (forc_struct *, river_struct *, int, int);


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
double          AvgKV (double, double, double, double, double, double, double,
    double, double);
double          EffKH (double, double, double, double, double, double);

void            PrtInit (pihm_struct, char *);
void            PrintData (prtctrl_struct *, int, int, int, int, int);
int             MacroporeStatus (double, double, double, double, double,
    double, double);
double          KrFunc (double, double, double);
double          Psi (double, double, double);

void            Summary (pihm_struct, N_Vector, double);
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
#define PIHMError(i)  _PIHMError(__FILE__, __LINE__, __FUNCTION__, i)
void            _PIHMError (const char *, int, const char *, int);

#ifdef _NOAH_
void            InitLsm (elem_struct *, int, ctrl_struct, noahtbl_struct,
    calib_struct);
void            CalcLatFlx (const wstate_struct *, const pstate_struct *,
    wflux_struct *);
int             FindWT (const double *, int, double, double *);
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
void            CalHum (pstate_struct *, estate_struct *);
void            Noah (int, pihm_struct);
double          FrozRain (double, double);
void            AvgFlux (elem_struct *, int, int);
void            SfcDifOff (pstate_struct *, const lc_struct *, double, double,
    int);
void            SFlx (wstate_struct *, wflux_struct *, const wflux_struct *,
    estate_struct *, eflux_struct *, pstate_struct *, lc_struct *, epconst_struct *,
    soil_struct *,
#ifdef _CYCLES_
    comm_struct *, residue_struct *,
#endif
    int);
double          CSnow (double);
void            SnowNew (const estate_struct *, double, pstate_struct *);
double          SnFrac (double, double, double, double);
void            AlCalc (pstate_struct *, double, int);
double          TDfCnd (double, double, double, double, double);
double          Snowz0 (double, double, double);
void            Penman (wflux_struct *, estate_struct *, eflux_struct *, pstate_struct *,
    double *, double, int, int);
void            CanRes (wstate_struct *, estate_struct *, eflux_struct *, pstate_struct *,
    const double *, const soil_struct *, const lc_struct *, const epconst_struct *);
void            DEvap (const wstate_struct *, wflux_struct *, const pstate_struct *,
    const lc_struct *, const soil_struct *);
void            Evapo (wstate_struct *, wflux_struct *, pstate_struct *,
    const lc_struct *, soil_struct *,
#ifdef _CYCLES_
    comm_struct *, residue_struct *, const estate_struct *,
#endif
    const double *, double);
void            Transp (const wstate_struct *, wflux_struct *, const pstate_struct *,
    const lc_struct *, const soil_struct *, const double *);
void            NoPac (wstate_struct *, wflux_struct *, const wflux_struct *,
    estate_struct *, eflux_struct *, pstate_struct *, lc_struct *, soil_struct *,
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
void            ShFlx (wstate_struct *, estate_struct *, eflux_struct *, pstate_struct *,
    const lc_struct *, const soil_struct *, double, double,
    double, const double *, double);
void            SmFlx (wstate_struct *, wflux_struct *, const wflux_struct *,
    pstate_struct *, const lc_struct *, const soil_struct *,
#ifdef _CYCLES_
    residue_struct *,
#endif
    const double *, double, double);
void            HRT (wstate_struct *, estate_struct *, eflux_struct *, pstate_struct *,
    const lc_struct *, const soil_struct *, double *,
    const double *, double, double, double, double, double *,
    double *, double *);
void            SRT (wstate_struct *, wflux_struct *, const wflux_struct *, pstate_struct *,
    const soil_struct *,
#ifdef _CYCLES_
    residue_struct *,
#endif
    double *, double *, double *, double *, double *, const double *, double);
void            SStep (wstate_struct *, wflux_struct *, const wflux_struct *,
    pstate_struct *, const soil_struct *, double *, double,
    const double *, double *, double *, double *, double *, double);
void            WDfCnd (double *, double *, double, double, double, int,
    const soil_struct *, const pstate_struct *);
void            SnoPac (wstate_struct *, wflux_struct *, const wflux_struct *,
    estate_struct *, eflux_struct *, pstate_struct *, lc_struct *, soil_struct *,
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
void            EnKFRead (char *, enkf_struct);
double          Randn ();
void            MapVar (var_struct *, int, int);
void            Perturb (char *, enkf_struct, char *);
void            MapVar (var_struct *, int, int);
void            Calib2Mbr (calib_struct, double *);
void            Mbr2Cal (calib_struct *, const double *);
void            WriteParamOutput (int, enkf_struct, int, char *);
void            WriteCalFile (enkf_struct, char *);
void            JobHandout (int, int, int, ensmbr_struct *, double *, int,
    int);
void            JobRecv (int *, int *, int *, double *, int);
void            PrintEnKFStatus (int, int);
void            JobHandIn (int);
void            WritePara (char *, int, int, int);

void            EnKFCore (double *, double, double, double *, int);
void            EnKF (enkf_struct, int, char *);
void            ReadObs (int, char *, double *, double *);
void            InitOper (pihm_struct, enkf_struct);
void            DisOper (obs_struct *, var_struct *, pihm_struct);
void            ReadFcst (enkf_struct, obs_struct, double *);
void            ReadVar (char *, char *, enkf_struct, int);
void            UpdAnlys (enkf_struct, double, double, double *);
void            CovInflt (enkf_struct, enkf_struct);
void            WriteEnKFOut (char *, enkf_struct, char *, int);
void            GenRandNum (int, int, double **, double, double);
void            LandSfcTmpOper (obs_struct *, var_struct *, pihm_struct);
void            COSMOSOper (obs_struct *, var_struct *, pihm_struct);
void            FreeEns (enkf_struct);
int             FindVar (var_struct *, char *);
void            InitEns (enkf_struct);
void            Parallel (int, int, char *);
#endif

#ifdef _CYCLES_
void            DailyCycles (int, pihm_struct);
void            FirstDOY (int *, int, int, soilc_struct *, residue_struct *,
    const soil_struct *);
void            ReadCyclesCtrl (char *, agtbl_struct *, ctrl_struct *, int);
void            ReadSoilInit (char *, soiltbl_struct *);
void            ReadCrop (char *, croptbl_struct *);
void            ReadOperation (const agtbl_struct *, mgmttbl_struct *,
    const croptbl_struct *);
int             CropExist (char *, const croptbl_struct *);
void            InitCycles (elem_struct *, int, river_struct *, int,
    const ctrl_struct *, const mgmttbl_struct *, const agtbl_struct *,
    const croptbl_struct *, const soiltbl_struct *);
void            InitializeSoil (soil_struct *, const soiltbl_struct *,
    const pstate_struct *, int);
double          BulkDensity (double, double, double);
void            InitializeResidue (residue_struct *, int);
void            InitializeSoilCarbon (soilc_struct *, int);
void            ComputeFactorComposite (soilc_struct *, int, int, int,
    soil_struct *);
void            ComputeSoilCarbonBalanceMB (soilc_struct *, int,
    residue_struct *, soil_struct *, double *);
void            ComputeSoilCarbonBalance (soilc_struct *, int,
    residue_struct *, soil_struct *, double *);
void            StoreOutput (soilc_struct *, int, int, double *);
double          Aeration (double);
double          Moisture (double);
double          TemperatureFunction (double);
double          MaximumAbgdHumificationFactor (double);
double          MaximumRootHumificationFactor (double);
double          MaximumRhizHumificationFactor (double);
double          MaximumManuHumificationFactor (double);
double          NitrogenMineralization (double, double, double, double);
double          CNdestiny (double, double);
double          PoolNitrogenMineralization (double, double, double, double,
    double);
double          Function_CNnew (double, double);
void            WaterUptake (comm_struct *, soil_struct *, double,
    wflux_struct *, double, double);
double          TemperatureLimitation (double, double, double);
void            CalcRootFraction (double *, soil_struct *, crop_struct *);
int             DOY (int);
int             IsLeapYear (int);
void            DailyOperations (int, int, cropmgmt_struct *, comm_struct *,
    residue_struct *, ctrl_struct *, snow_struct *, soil_struct *,
    soilc_struct *, weather_struct *);
double          Depth_Limitation_To_Evaporation (double);
double          Water_Content_Limitation_To_Evaporation (double, double,
    double);
void            Evaporation (soil_struct *, const comm_struct *,
    residue_struct *, double, double);
void            LastDOY (int, int, soil_struct *, soilc_struct *,
    residue_struct *);
void            GrowingCrop (int, int, comm_struct *, residue_struct *,
    const ctrl_struct *, soil_struct *, soilc_struct *, cropmgmt_struct *,
    const weather_struct *, const snow_struct *);
void            CropStage (int, comm_struct *, int);
double          FinalHarvestDate (int, int, double, double, double);
void            Phenology (int, int, const weather_struct *, comm_struct *);
double          ThermalTime (double, double, double, double);
void            RadiationInterception (int, int, comm_struct *);
void            Processes (int, int, int, comm_struct *, residue_struct *,
    const weather_struct *, soil_struct *, soilc_struct *);
void            CropNitrogenConcentration (double *, double *, double *,
    double *, double *, double *, double *, double, const crop_struct *);
void            CropNitrogenStress (double, double, double, crop_struct *);
void            CropGrowth (int, int, double *, double, crop_struct *,
    residue_struct *, const weather_struct *);
void            CropNitrogenDemand (double, double, double *, double *,
    double *, double *, crop_struct *);
void            PotentialSoluteUptakeOption2 (double *, double *, double, int,
    const double *, const double *, const double *, const double *,
    const double *);
void            CropNitrogenUptake (double *, double *, double *, double *,
    double *, int, double, double, double *, double *, double *,
    comm_struct *, soil_struct *);
void            DistributeRootDetritus (double, double, double, double,
    const soil_struct *, const crop_struct *, residue_struct *,
    soilc_struct *);
double          ShootBiomassPartitioning (double, double, double);
double          TemperatureFunctionGrowth (double, double, double, double);
int             ForcedClipping (int, comm_struct *);
void            GrainHarvest (int, int, crop_struct *, residue_struct *,
    soil_struct *, soilc_struct *);
void            ComputeColdDamage (int, int, crop_struct *,
    const weather_struct *, const snow_struct *, residue_struct *);
double          ColdDamage (double, double, double);
void            ForageAndSeedHarvest (int, int, crop_struct *,
    residue_struct *, soil_struct *, soilc_struct *);
void            HarvestCrop (int, int, crop_struct *, residue_struct *,
    soil_struct *, soilc_struct *);
void            PlantingCrop (comm_struct *, const cropmgmt_struct *, int);
void            AddCrop (crop_struct *);
void            KillCrop (crop_struct *);
void            UpdateCommunity (comm_struct *);
double          ComputeHarvestIndex (double, double, double, double, double);
int             IsOperationToday (int, int, op_struct *, int, int *);
void            ApplyFertilizer (op_struct *, soil_struct *,
    residue_struct *);
void            UpdateOperationStatus (op_struct *, int);
void            FieldOperation (int, int, int, cropmgmt_struct *,
    comm_struct *, soil_struct *, residue_struct *, ctrl_struct *,
    soilc_struct *, weather_struct *);
void            ExecuteTillage (double *, const op_struct *, double *,
    soil_struct *, residue_struct *);
void            TillageFactorSettling (double *, int, const double *,
    const double *);
double          Fraction (double, double, double, double, double);
void            ComputeTillageFactor (const op_struct *, double *,
    const soil_struct *, const double *, double);
double          ComputeTextureFactor (double);
void            ComputeResidueCover (residue_struct *);
void            ResidueEvaporation (residue_struct *, soil_struct *,
    const comm_struct *, double, double);
void            NitrogenTransformation (int, int, soil_struct *,
    const comm_struct *, const residue_struct *, const weather_struct *,
    const soilc_struct *);
void            Nitrification (double *, double *, soil_struct *,
    const soilc_struct *);
void            Denitrification (double *, double *, soil_struct *,
    const soilc_struct *);
void            Volatilization (int, int, double *, soil_struct *,
    const comm_struct *, const residue_struct *, const weather_struct *);
double          N2OFractionNitrification (double);
double          pHFunction (double);
double          VolatilizationDepthFunction (double);
double          AirMolarDensity (double, double);
double          BoundaryLayerConductance (double, double, double, double);
void            ResidueWetting (residue_struct *, double *);
double          FindIrrigationVolume (int, double, const soil_struct *);
void            SoluteTransport (elem_struct *, int, river_struct *, int,
    double);
void            Adsorption (const double *sldpth, const double *,
    const double *, int, double, solute_struct *);
double          LinearEquilibriumConcentration (double, double, double,
    double, double);
double          LinearEquilibriumSoluteMass (double, double, double, double,
    double);
void            Elem2ElemSolTrnsp (const elem_struct *, const elem_struct *,
    double *, const double *, double, double *, double *);
void            Elem2RiverSolTrnsp (const elem_struct *, const river_struct *,
    double, double *, const double *, const double *, double, double *,
    double *);
void            River2RiverSolTrnsp (river_struct *, const river_struct *,
    double *, const double *, const double *, double, double *, double *);
void            InitCropSV (crop_struct *);
#endif

#ifdef _BGC_
void            ReadEPC (epctbl_struct *);
void            ReadBGC (char *, ctrl_struct *, co2control_struct *,
    ndepcontrol_struct *, char *, char *);
void            InitElemStor (stor_struct *, int, int);
void            InitRiverStor (river_stor_struct *, int, int);
void            ReadAnnFile (tsdata_struct *, char *);
void            InitBGC (elem_struct *, int, river_struct *, int,
    const epctbl_struct *, const ctrl_struct *);
void            restart_input (cstate_struct *, nstate_struct *, epvar_struct *, bgc_ic_struct *);
void            InitBGCVar (elem_struct *, int, river_struct *, int,
    cinit_struct, cstate_struct, nstate_struct, char *, int);
void            firstday (const epconst_struct *, const cinit_struct *, epvar_struct *, cstate_struct *, nstate_struct *);
void            zero_srcsnk (cstate_struct *, nstate_struct *, summary_struct *);
void            Save2Stor (pihm_struct, int, int, int);
void            BGCSpinup (char *, pihm_struct, char *);
void            restart_output (cstate_struct *, nstate_struct *, epvar_struct *, bgc_ic_struct *);
double          GetCO2 (tsdata_struct, int);
double          GetNdep (tsdata_struct, int);
void            DayMet (const stor_struct *, daily_struct *, int);
void            RiverDayMet (const river_stor_struct *, river_daily_struct *, int);
void            PrecisionControl (cstate_struct *cs, nstate_struct *ns);
void            MakeZeroFluxStruct (cflux_struct *, nflux_struct *);
void            Phenology (const epconst_struct *, const daily_struct *, phenology_struct *, epvar_struct *, cstate_struct *, cflux_struct *, nstate_struct *, nflux_struct *nf);
void            LeafLitFall (const epconst_struct *, double, cflux_struct *, nflux_struct *);
void            FRootLitFall (const epconst_struct *, double, cflux_struct *, nflux_struct *);
void            RadTrans (const cstate_struct *, const daily_struct *, eflux_struct *, pstate_struct *, const epconst_struct *, epvar_struct *);
void            SoilPsi (const soil_struct *, double, double *);
void            MaintResp (const cstate_struct *, const nstate_struct *, const epconst_struct *, const daily_struct *, cflux_struct *, epvar_struct *);
void            CanopyCond (const epconst_struct *, const daily_struct *, pstate_struct *, const soil_struct *, epvar_struct *);
void            Photosynthesis (psn_struct *);
void            TotalPhotosynthesis (const epconst_struct *, const daily_struct *, const pstate_struct *, epvar_struct *, cflux_struct *, psn_struct *, psn_struct *);
void            Decomp (double, const epconst_struct *, epvar_struct *, cstate_struct *, cflux_struct *, nstate_struct *, nflux_struct *, ntemp_struct *);
void            DailyAllocation (cflux_struct *, cstate_struct *, nflux_struct *, nstate_struct *, epconst_struct *, epvar_struct *, ntemp_struct *, const int);
void            AnnualRates (const epconst_struct *, epvar_struct *);
void            GrowthResp (epconst_struct *, cflux_struct *);
void            DailyCarbonStateUpdate (cflux_struct *, cstate_struct *, int, int, int);
void            DailyNitrogenStateUpdate (nflux_struct *, nstate_struct *, int alloc, int woody, int evergreen);
void            Mortality (const epconst_struct *, cstate_struct *, cflux_struct *, nstate_struct *, nflux_struct *);
void            CheckCarbonBalance (cstate_struct *, double *, int);
void            CheckNitrogenBalance (nstate_struct *, double *, int);
void            CSummary (cflux_struct *, cstate_struct *, summary_struct *);
void            NLeaching (elem_struct *, int, river_struct *, int);
void            DailyBgc (pihm_struct, int, int, int);
#endif

#endif
