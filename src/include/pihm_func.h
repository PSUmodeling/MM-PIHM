#ifndef PIHM_FUNC_HEADER
#define PIHM_FUNC_HEADER

#define _ARITH_

#ifdef _BGC_
#define NSV             4 * nelem + 4 * nriver
#else
#define NSV             3 * nelem + 2 * nriver
#endif

#define SURF(i)         i
#define UNSAT(i)        i + nelem
#define GW(i)           i + 2 * nelem
#define RIVSTG(i)       i + 3 * nelem
#define RIVGW(i)        i + 3 * nelem + nriver

#ifdef _BGC_
#define SMINN(i)        i + 3 * nelem + 2 * nriver
#define STREAMN(i)      i + 4 * nelem + 2 * nriver
#define RIVBEDN(i)      i + 4 * nelem + 3 * nriver
#endif

/*
 * Function Declarations
 */
void            ApplyBC (forc_struct *, elem_struct *, river_struct *, int);
void            ApplyElemBC (forc_struct *, elem_struct *, int);
void            ApplyForcing (forc_struct *, elem_struct *, int
#ifdef _NOAH_
    , ctrl_struct *, double, double, double, double
#endif
    );
void            ApplyLAI (forc_struct *, elem_struct *, int);
void            ApplyMeteoForc (forc_struct *, elem_struct *, int
#ifdef _NOAH_
    , int, double, double, double, double
#endif
    );
void            ApplyRiverBC (forc_struct *, river_struct *, int);
void            AsciiArt ();
double          AvgKV (double, double, double, double, double, double, double,
    double);
double          AvgYsfc (double, double, double);
double          AvgY (double, double, double);
void            BKInput (char *, char *);
void            CalcModelStep (ctrl_struct *);
void            CheckFile (FILE *, char *);
void            CorrectElevation (elem_struct *, river_struct *);
int             CountLine (FILE *, char *, int, ...);
int             CountOccurance (FILE *, char *);
void            CreateOutputDir (char *);
double          DhByDl (double *, double *, double *);
double          EffKH (double, double, double, double, double, double);
double          EffKinf (double, double, int, double, double, double);
double          EffKV (double, int, double, double, double);
double          FieldCapacity (double, double, double, double, double);
void            FindLine (FILE *, char *, int *, const char *);
void            FreeData (pihm_struct);
void            FrictSlope (elem_struct *, river_struct *, int, double *,
    double *);
void            Hydrol (pihm_struct);
void            Initialize (pihm_struct, N_Vector);
void            InitEFlux (eflux_struct *);
void            InitEState (estate_struct *);
void            InitForcing (elem_struct *, forc_struct *,
    const calib_struct *
#ifdef _BGC_
    , int, int
#endif
    );
void            InitLC (elem_struct *, const lctbl_struct *,
    const calib_struct *);
void            InitMeshStruct (elem_struct *, const meshtbl_struct *);
void            InitOutputFile (prtctrl_struct *, int, int);
void            InitRiver (river_struct *, elem_struct *, const rivtbl_struct *,
    const shptbl_struct *, const matltbl_struct *, const meshtbl_struct *,
    const calib_struct *);
void            InitRiverWFlux (river_wflux_struct *);
void            InitRiverWState (river_wstate_struct *);
void            InitSoil (elem_struct *, const soiltbl_struct *,
#ifdef _NOAH_
    const noahtbl_struct *,
#endif
    const calib_struct *);
void            InitSurfL (elem_struct *, river_struct *, const meshtbl_struct *);
void            InitTopo (elem_struct *, const meshtbl_struct *);
void            InitVar (elem_struct *, river_struct *, N_Vector);
void            InitWFlux (wflux_struct *);
void            InitWState (wstate_struct *);
void            IntcpSnowET (int, double, pihm_struct);
void            IntrplForcing (tsdata_struct *, int, int);
double          KrFunc (double, double, double);
void            LateralFlow (pihm_struct);
int             MacroporeStatus (double, double, double, double, double,
    double);
void            MapOutput (char *, pihm_struct, char *);
void            MassBalance (wstate_struct *, wstate_struct *, wflux_struct *,
    double *, const soil_struct *, double, double);
double          MonthlyLAI (int, int);
double          MonthlyMF (int);
double          MonthlyRL (int, int);
#ifdef _CVODE_OMP
#define N_VNew(N)       N_VNew_OpenMP(N, nthreads)
#else
#define N_VNew(N)       N_VNew_Serial(N);
#endif
void            NextLine (FILE *, char *, int *);
#ifdef _CVODE_OMP
#define NV_DATA         NV_DATA_OMP
#define NV_Ith          NV_Ith_OMP
#else
#define NV_DATA         NV_DATA_S
#define NV_Ith          NV_Ith_S
#endif
int             ODE (realtype, N_Vector, N_Vector, void *);
double          OverlandFlow (double, double, double, double, double);
double          OLFEleToRiv (double, double, double, double, double, double);
void            ParseCmdLineParam (int, char *[], char *);
#define PIHMexit(...)  _PIHMexit(__FILE__, __LINE__, __FUNCTION__, __VA_ARGS__)
void            _PIHMexit (const char *, int, const char *, int);
#define PIHMprintf(...)   _PIHMprintf(__FILE__, __LINE__, __FUNCTION__, __VA_ARGS__)
void            _PIHMprintf (const char *, int, const char *, int,
    const char *, ...);
void            PIHM (pihm_struct, void *, N_Vector, int, int);
pihm_t_struct   PIHMTime (int);
void            PrintData (prtctrl_struct *, int, int, int, int);
void            PrtInit (elem_struct *, river_struct *, char *);
double          Psi (double, double, double);
double          PtfAlpha (double, double, double, double, int);
double          PtfBeta (double, double, double, double, int);
double          PtfKV (double, double, double, double, int);
double          PtfThetaR (double, double, double, double, int);
double          PtfThetaS (double, double, double, double, int);
double          Qtz (int);
void            ReadAlloc (char *, pihm_struct);
void            ReadAtt (char *, atttbl_struct *);
void            ReadBC (char *, forc_struct *);
void            ReadCalib (char *, calib_struct *);
void            ReadForc (char *, forc_struct *);
void            ReadGeol (char *, geoltbl_struct *);
void            ReadIC (char *, elem_struct *, river_struct *);
int             ReadKeyword (char *, char *, void *, char, char *, int);
void            ReadLAI (char *, forc_struct *, const atttbl_struct *);
void            ReadLC (char *, lctbl_struct *);
void            ReadMesh (char *, meshtbl_struct *);
void            ReadPara (char *, ctrl_struct *);
int             ReadPrtCtrl (char *, char *, char *, int);
void            ReadRiv (char *, rivtbl_struct *, shptbl_struct *,
    matltbl_struct *, forc_struct *);
void            ReadSoil (char *, soiltbl_struct *);
int             ReadTS (char *, int *, double *, int);
int             Readable (char *);
void            RiverFlow (pihm_struct);
void            RiverToEle (river_struct *, elem_struct *, elem_struct *,
    int, double, double *, double *, double *);
double          _RivWdthAreaPerim (int, int, double, double);
#define RivArea(...)    _RivWdthAreaPerim(RIVER_AREA, __VA_ARGS__)
#define RivEqWid(...)   _RivWdthAreaPerim(RIVER_WDTH, __VA_ARGS__)
#define RivPerim(...)   _RivWdthAreaPerim(RIVER_PERIM, __VA_ARGS__)
void            SaturationIC (elem_struct *, river_struct *);
void            SetCVodeParam (pihm_struct, void *, N_Vector);
int             SoilTex (double, double);
void            SolveCVode (int *, int, int, void *, N_Vector);
int             StrTime (const char *);
void            Summary (pihm_struct, N_Vector, double);
void            UpdPrintVar (prtctrl_struct *, int, int);
void            VerticalFlow (pihm_struct);
double          WiltingPoint (double, double, double, double);

/*
 * Noah functions
 */
#ifdef _NOAH_
void            AlCalc (pstate_struct *, double, int);
double          AvgElev (elem_struct *);
double          CSnow (double);
void            CalHum (pstate_struct *, estate_struct *);
void            CalcLatFlx (const pstate_struct *, wflux_struct *, double);
void            CalcSlopeAspect (elem_struct *, const meshtbl_struct *);
void            CanRes (wstate_struct *, estate_struct *, eflux_struct *,
    pstate_struct *, const soil_struct *,
    const epconst_struct *);
void            DEvap (const wstate_struct *, wflux_struct *,
    const pstate_struct *, const lc_struct *, const soil_struct *);
void            DefSldpth (double *, int *, double *, double, const double *, int);
void            Evapo (wstate_struct *, wflux_struct *, pstate_struct *,
    const lc_struct *, soil_struct *,
#ifdef _CYCLES_
    comm_struct *, residue_struct *, const estate_struct *,
#endif
    double);
int             FindLayer (const double *, int, double);
int             FindWT (const double *, int, double, double *);
double          FrozRain (double, double);
double          GWTransp (double, double *, int, int);
void            HRT (wstate_struct *, estate_struct *, eflux_struct *,
    pstate_struct *, const lc_struct *, const soil_struct *, double *,
    double, double, double, double, double *, double *, double *);
void            InitLsm (elem_struct *, const ctrl_struct *,
    const noahtbl_struct *, const calib_struct *);
void            NoPac (wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, pstate_struct *, lc_struct *, soil_struct *,
#ifdef _CYCLES_
    comm_struct *, residue_struct *,
#endif
    double, double);
void            Noah (pihm_struct);
void            NoahHydrol (elem_struct *, double);
void            PcpDrp (wstate_struct *, wflux_struct *, const lc_struct *,
    double, double);
void            Penman (wflux_struct *, estate_struct *, eflux_struct *,
    pstate_struct *, double *, double, int, int);
double          Pslhs (double);
double          Pslhu (double);
double          Pslmu (double);
double          Pslms (double);
double          Psphs (double);
double          Psphu (double);
double          Pspms (double);
double          Pspmu (double);
void            ReadLsm (char *, double *, double *, ctrl_struct *,
    noahtbl_struct *);
void            ReadRad (char *, forc_struct *);
void            RootDist (const double *, int, int, double *);
void            Rosr12 (double *, double *, double *, double *, double *,
    double *, int);
void            SfcDifOff (pstate_struct *, const lc_struct *, double, double,
    int);
void            SFlx (wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, pstate_struct *, lc_struct *, epconst_struct *,
    soil_struct *,
#ifdef _CYCLES_
    comm_struct *, residue_struct *,
#endif
    int);
void            SRT (wstate_struct *, wflux_struct *, pstate_struct *,
    const soil_struct *,
#ifdef _CYCLES_
    residue_struct *, double,
#endif
    double *, double *, double *, double *, double *);
void            ShFlx (wstate_struct *, estate_struct *, eflux_struct *,
    pstate_struct *, const lc_struct *, const soil_struct *, double, double,
    double, double);
void            SmFlx (wstate_struct *, wflux_struct *, pstate_struct *,
    const soil_struct *,
#ifdef _CYCLES_
    residue_struct *,
#endif
    double );
double          SnFrac (double, double, double, double);
void            SnkSrc (double *, double, double, double *,
    const soil_struct *, const double *, double, int, double);
void            SnoPac (wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, pstate_struct *, lc_struct *, soil_struct *,
#ifdef _CYCLES_
    comm_struct *, residue_struct *,
#endif
    int, double, double, double, double);
void            SnowNew (const estate_struct *, double, pstate_struct *);
void            SnowPack (double, double, double *, double *, double, double);
double          Snowz0 (double, double, double);
void            SStep (wstate_struct *, wflux_struct *, pstate_struct *,
    const soil_struct *, double *, double *, double *, double *, double *,
    double);
void            SunPos (int, double, double, double, double, spa_data *);
double          TBnd (double, double, const double *, double, int, int);
double          TDfCnd (double, double, double, double, double);
double          TmpAvg (double, double, double, const double *, int);
double          TopoRadn (double, double, double, double, double, double,
    const double *, double);
void            Transp (const wstate_struct *, wflux_struct *,
    const pstate_struct *, const lc_struct *, const soil_struct *);
void            WDfCnd (double *, double *, double, double, int,
    const soil_struct *, const pstate_struct *);
#endif

#ifdef _DAILY_
void            DailyVar (int, int, pihm_struct);
void            InitDailyStruct (pihm_struct);
#endif

#ifdef _CYCLES_
#define Cycles_printf   PIHMprintf
#define Cycles_exit     PIHMexit
void            DailyCycles (int, pihm_struct);
void            FirstDOY (int *, int, int, soilc_struct *, residue_struct *,
    const soil_struct *);
void            ReadCyclesCtrl (char *, agtbl_struct *, ctrl_struct *);
void            ReadCyclesIC (char *, elem_struct *, river_struct *);
void            ReadSoilInit (char *, soiltbl_struct *);
void            ReadCrop (char *, croptbl_struct *);
void            ReadOperation (const agtbl_struct *, mgmttbl_struct *,
    const croptbl_struct *);
int             CropExist (char *, const croptbl_struct *);
void            InitCycles (elem_struct *, river_struct *,
    const ctrl_struct *, const mgmttbl_struct *, const agtbl_struct *,
    const croptbl_struct *, const soiltbl_struct *);
void            InitializeSoil (soil_struct *, const soiltbl_struct *,
    const pstate_struct *, const cyclesic_struct *, int, int);
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
double          ShootBiomassPartitioning (double, double, double, int);
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
int             IsOperationToday (int, int, const void *, int, int *, int *,
    int);
void            ApplyFertilizer (const fixfert_struct *, soil_struct *,
    residue_struct *);
void            UpdateOperationStatus (int *, int);
void            FieldOperation (int, int, int, cropmgmt_struct *,
    comm_struct *, soil_struct *, residue_struct *, ctrl_struct *,
    soilc_struct *, weather_struct *);
void            ExecuteTillage (double *, const tillage_struct *, double *,
    soil_struct *, residue_struct *);
void            TillageFactorSettling (double *, int, const double *,
    const double *);
double          Fraction (double, double, double, double, double);
void            ComputeTillageFactor (const tillage_struct *, double *,
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
void            SoluteTransport (elem_struct *, river_struct *, double);
void            Adsorption (const double *sldpth, const double *,
    const double *, int, double, solute_struct *);
double          LinearEquilibriumConcentration (double, double, double,
    double, double);
double          LinearEquilibriumSoluteMass (double, double, double, double,
    double);
void            Elem2ElemSolTrnsp (const elem_struct *, const elem_struct *,
    double *, const double *, double, double *, double *);
void            Elem2RiverSolTrnsp (const elem_struct *, const river_struct *,
    double, double *, const double *, double, double, double *, double *);
void            River2RiverSolTrnsp (river_struct *, const river_struct *,
    double *, double, double, double, double *, double *);
void            InitCropSV (crop_struct *);
void            WriteCyclesIC (char *, elem_struct *, river_struct *);
#endif

#ifdef _BGC_
void            BackgroundLitterfall (const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, nflux_struct *);
void            BgcSpinup (pihm_struct, N_Vector, void *);
void            CanopyCond (const epconst_struct *, epvar_struct *,
    const eflux_struct *, const pstate_struct *, const soil_struct *,
    const daily_struct *);
void            CheckBgcSS (elem_struct *, int, int, int);
void            CheckCarbonBalance (cstate_struct *, double *);
void            CheckNitrogenBalance (nstate_struct *, double *);
void            CSummary (cflux_struct *, cstate_struct *, summary_struct *);
void            DailyAllocation (cflux_struct *, const cstate_struct *,
    nflux_struct *, const nstate_struct *, const epconst_struct *,
    epvar_struct *, ntemp_struct *);
void            DailyBgc (pihm_struct, int);
void            DailyCarbonStateUpdate (cflux_struct *, cstate_struct *, int,
    int, int);
void            DailyNitrogenStateUpdate (nflux_struct *, nstate_struct *,
    solute_struct *, int, int, int);
void            Decomp (double, const epconst_struct *, epvar_struct *,
    cstate_struct *, cflux_struct *, nstate_struct *, nflux_struct *,
    ntemp_struct *);
void            EvergreenPhenology (const epconst_struct *, epvar_struct *,
    cstate_struct *);
void            FRootLitFall (const epconst_struct *, double, cflux_struct *,
    nflux_struct *);
void            FirstDay (elem_struct *, river_struct *,
    const cninit_struct *);
double          GetCO2 (tsdata_struct, int);
double          GetNdep (tsdata_struct, int);
void            GrowthResp (epconst_struct *, cflux_struct *);
void            InitBgc (elem_struct *, const epctbl_struct *);
void            InitBgcVar (elem_struct *, river_struct *, N_Vector);
void            LeafLitFall (const epconst_struct *, double, cflux_struct *,
    nflux_struct *);
void            LivewoodTurnover (const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, const nstate_struct *,
    nflux_struct *);
void            MaintResp (const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, const nstate_struct *,
    const daily_struct *);
void            MakeZeroFluxStruct (cflux_struct *, nflux_struct *);
void            Mortality (const epconst_struct *, cstate_struct *,
    cflux_struct *, nstate_struct *, nflux_struct *);
void            NTransport (pihm_struct);
void            OffsetLitterfall (const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, nflux_struct *);
void            OnsetGrowth (const epconst_struct *, const epvar_struct *,
    const cstate_struct *, cflux_struct *, const nstate_struct *,
    nflux_struct *);
void            Phenology (const epconst_struct *, epvar_struct *,
    cstate_struct *, cflux_struct *, nstate_struct *, nflux_struct *,
    const daily_struct *);
void            Photosynthesis (psn_struct *);
void            PrecisionControl (cstate_struct *cs, nstate_struct *ns);
void            RadTrans (const cstate_struct *, eflux_struct *,
    pstate_struct *, const epconst_struct *, epvar_struct *,
    const daily_struct *);
void            ReadAnnFile (tsdata_struct *, char *);
void            ReadBgc (char *, ctrl_struct *, co2control_struct *,
    ndepcontrol_struct *, cninit_struct *, char *, char *);
void            ReadBgcIC (char *, elem_struct *, river_struct *);
void            ReadEPC (epctbl_struct *);
void            ResetSpinupStat (elem_struct *);
void            RestartInput (cstate_struct *, nstate_struct *,
    epvar_struct *, bgcic_struct *);
void            RestartOutput (cstate_struct *, nstate_struct *,
    epvar_struct *, bgcic_struct *);
void            SeasonDecidPhenology (const epconst_struct *, epvar_struct *,
    const daily_struct *);
void            SoilPsi (const soil_struct *, double, double *);
void            TotalPhotosynthesis (const epconst_struct *, epvar_struct *,
    const pstate_struct *, cflux_struct *, psn_struct *, psn_struct *,
    daily_struct *);
void            WriteBgcIC (char *, elem_struct *, river_struct *);
void            ZeroSrcSnk (cstate_struct *, nstate_struct *,
    summary_struct *);
#endif
#endif
