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
void            ApplyBC (forc_struct *, elem_struct *, int, int);
void            ApplyForcing (forc_struct *, elem_struct *, int,
    river_struct *, int, int
#ifdef _NOAH_
    , ctrl_struct *, double, double, double, double
#endif
    );
void            ApplyLAI (forc_struct *, elem_struct *, int, int
#ifdef _BGC_
    ,int
#endif
    );
void            ApplyMeteoForc (forc_struct *, elem_struct *, int, int
#ifdef _NOAH_
    , int, double, double, double, double
#endif
    );
void            ApplyRiverBC (forc_struct *, river_struct *, int, int);
void            AsciiArt ();
double          AvgKV (double, double, double, double, double, double, double,
    double);
double          AvgYsfc (double, double, double);
double          AvgY (double, double, double);
void            BKInput (char *, char *);
void            CalcModelStep (ctrl_struct *);
void            CheckFile (FILE *, char *);
void            CorrectElevation (elem_struct *, int, river_struct *, int);
int             CountLine (FILE *, char *, int, ...);
int             CountOccurance (FILE *, char *);
void            CreateOutputDir (char *, int);
double          DhByDl (double *, double *, double *);
double          EffKH (double, double, double, double, double, double);
double          EffKinf (double, double, int, double, double, double);
double          EffKV (double, int, double, double, double);
double          FieldCapacity (double, double, double, double, double);
void            FindLine (FILE *, char *, int *, const char *);
void            FreeData (pihm_struct);
int             Hydrol (realtype, N_Vector, N_Vector, void *);
void            Initialize (pihm_struct, N_Vector);
void            InitEFlux (eflux_struct *);
void            InitEState (estate_struct *);
void            InitForcing (elem_struct *, int, forc_struct *, calib_struct
#ifdef _BGC_
    , int, int, int
#endif
    );
void            InitLC (elem_struct *, int, lctbl_struct, calib_struct);
void            InitMeshStruct (elem_struct *, int, meshtbl_struct);
void            InitOutputFile (prtctrl_struct *, int, int);
void            InitRiver (river_struct *, int, elem_struct *, rivtbl_struct,
    shptbl_struct, matltbl_struct, meshtbl_struct, calib_struct);
void            InitRiverWFlux (river_wflux_struct *);
void            InitRiverWState (river_wstate_struct *);
void            InitSoil (elem_struct *, int, soiltbl_struct,
#ifdef _NOAH_
    noahtbl_struct,
#endif
    calib_struct);
void            InitSurfL (elem_struct *, int, river_struct *,
    meshtbl_struct);
void            InitTopo (elem_struct *, int, meshtbl_struct);
void            InitVar (elem_struct *, int, river_struct *, int, N_Vector);
void            InitWFlux (wflux_struct *);
void            InitWState (wstate_struct *);
void            IntcpSnowET (int, double, pihm_struct);
void            IntrplForcing (tsdata_struct, int, int);
double          KrFunc (double, double, double);
void            LateralFlow (pihm_struct);
int             MacroporeStatus (double, double, double, double, double,
    double);
void            MapOutput (char *, pihm_struct, char *);
double          MonthlyLAI (int, int);
double          MonthlyMF (int);
double          MonthlyRL (int, int);
void            NextLine (FILE *, char *, int *);
double          OverlandFlow (double, double, double, double, double);
double          OLFEleToRiv (double, double, double, double, double, double);
void            ParseCmdLineParam (int, char *[], int *, char *);
#define PIHMexit(...)  _PIHMexit(__FILE__, __LINE__, __FUNCTION__, __VA_ARGS__)
void            _PIHMexit (const char *, int, const char *, int);
#define PIHMprintf(...)   _PIHMprintf(__FILE__, __LINE__, __FUNCTION__, __VA_ARGS__)
void            _PIHMprintf (const char *, int, const char *, int,
    const char *, ...);
void            PIHM (char *, char *, int
#ifdef _ENKF_
    , int, int, int, double *
#endif
    );
void            PrintData (prtctrl_struct *, int, int, int, int, int);
void            PrtInit (pihm_struct, char *);
double          Psi (double, double, double);
double          PtfAlpha (double, double, double, double, int);
double          PtfBeta (double, double, double, double, int);
double          PtfKV (double, double, double, double, int);
double          PtfThetaR (double, double, double, double, int);
double          PtfThetaS (double, double, double, double, int);
double          Qtz (int);
void            ReadAlloc (char *, pihm_struct);
void            ReadAtt (char *, atttbl_struct *, int);
void            ReadBC (char *, forc_struct *);
void            ReadCalib (char *, calib_struct *);
void            ReadForc (char *, forc_struct *);
void            ReadGeol (char *, geoltbl_struct *);
void            ReadIC (char *, elem_struct *, int, river_struct *, int);
int             ReadKeyword (char *, char *, void *, char, char *, int);
void            ReadLAI (char *, forc_struct *, int, const atttbl_struct *);
void            ReadLC (char *, lctbl_struct *);
void            ReadMesh (char *, meshtbl_struct *);
void            ReadPara (char *, ctrl_struct *);
void            ReadRiv (char *, rivtbl_struct *, shptbl_struct *,
    matltbl_struct *, forc_struct *);
void            ReadSoil (char *, soiltbl_struct *);
int             ReadTS (char *, int *, double *, int);
int             Readable (char *);
void            RiverFlow (pihm_struct);
void            RiverToEle (river_struct *, elem_struct *, elem_struct *,
    int, double *, double *, double *);
double          _RivWdthAreaPerim (int, int, double, double);
#define RivArea(...)    _RivWdthAreaPerim(RIVER_AREA, __VA_ARGS__)
#define RivEqWid(...)   _RivWdthAreaPerim(RIVER_WDTH, __VA_ARGS__)
#define RivPerim(...)   _RivWdthAreaPerim(RIVER_PERIM, __VA_ARGS__)
void            SaturationIC (elem_struct *, int, river_struct *, int);
void            SetCVodeParam (pihm_struct, void *, N_Vector);
int             SoilTex (double, double);
void            SolveCVode (int *, int, int, void *, N_Vector);
void            Summary (pihm_struct, N_Vector, double);
void            VerticalFlow (pihm_struct);
double          WiltingPoint (double, double, double, double);

/*
 * Noah functions
 */
#ifdef _NOAH_
void            AlCalc (pstate_struct *, double, int);
double          AvgElev (elem_struct *, int);
double          CSnow (double);
void            CalHum (pstate_struct *, estate_struct *);
void            CalcLatFlx (const pstate_struct *, wflux_struct *, double);
void            CalcSlopeAspect (elem_struct *, int, meshtbl_struct);
void            CanRes (wstate_struct *, estate_struct *, eflux_struct *,
    pstate_struct *, const double *, const soil_struct *, const lc_struct *,
    const epconst_struct *);
void            DEvap (const wstate_struct *, wflux_struct *,
    const pstate_struct *, const lc_struct *, const soil_struct *);
void            DefSldpth (double *, int *, double, const double *, int);
void            Evapo (wstate_struct *, wflux_struct *, pstate_struct *,
    const lc_struct *, soil_struct *,
#ifdef _CYCLES_
    comm_struct *, residue_struct *, const estate_struct *,
#endif
    const double *, double);
int             FindLayer (const double *, int, double);
int             FindWT (const double *, int, double, double *);
double          FrozRain (double, double);
double          GWTransp (double, double *, int, int);
void            HRT (wstate_struct *, estate_struct *, eflux_struct *,
    pstate_struct *, const lc_struct *, const soil_struct *, double *,
    const double *, double, double, double, double, double *, double *,
    double *);
void            InitLsm (elem_struct *, int, ctrl_struct, noahtbl_struct,
    calib_struct);
void            NoPac (wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, pstate_struct *, lc_struct *, soil_struct *,
#ifdef _CYCLES_
    comm_struct *, residue_struct *,
#endif
    const double *, double, double);
void            Noah (pihm_struct);
void            NoahHydrol (pihm_struct, double);
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
    residue_struct *,
#endif
    double *, double *, double *, double *, double *, const double *, double);
void            ShFlx (wstate_struct *, estate_struct *, eflux_struct *,
    pstate_struct *, const lc_struct *, const soil_struct *, double, double,
    double, const double *, double);
void            SmFlx (wstate_struct *, wflux_struct *, pstate_struct *,
    const soil_struct *,
#ifdef _CYCLES_
    residue_struct *,
#endif
    const double *, double );
double          SnFrac (double, double, double, double);
void            SnkSrc (double *, double, double, double *,
    const soil_struct *, const double *, int, double, int, double);
void            SnoPac (wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, pstate_struct *, lc_struct *, soil_struct *,
#ifdef _CYCLES_
    comm_struct *, residue_struct *,
#endif
    int, const double *, double, double, double, double);
void            SnowNew (const estate_struct *, double, pstate_struct *);
void            SnowPack (double, double, double *, double *, double, double);
double          Snowz0 (double, double, double);
void            SStep (wstate_struct *, wflux_struct *, pstate_struct *,
    const soil_struct *, double *, const double *, double *, double *,
    double *, double *, double);
void            SunPos (int, double, double, double, double, double *,
    double *);
double          TBnd (double, double, const double *, double, int, int);
double          TDfCnd (double, double, double, double, double);
double          TmpAvg (double, double, double, const double *, int, int);
double          TopoRadn (double, double, double, double, double, double,
    const double *, double);
void            Transp (const wstate_struct *, wflux_struct *,
    const pstate_struct *, const lc_struct *, const soil_struct *,
    const double *);
void            WDfCnd (double *, double *, double, double, double, int,
    const soil_struct *, const pstate_struct *);
#endif

#ifdef _DAILY_
void            DailyVar (int, int, pihm_struct);
void            InitDailyStruct (pihm_struct);
#endif

#ifdef _ENKF_
void            Calib2Mbr (calib_struct, double *);
void            COSMOSOper (obs_struct *, var_struct *, pihm_struct);
void            CovInflt (enkf_struct, enkf_struct);
void            DisOper (obs_struct *, var_struct *, pihm_struct);
void            EnKF (double *, double, double, double *, int);
void            EnKFDA (enkf_struct, int, char *);
int             FindVar (var_struct *, char *);
void            FreeEns (enkf_struct);
void            GenRandNum (int, int, double **, double, double);
void            InitEns (enkf_struct);
void            InitOper (pihm_struct, enkf_struct);
void            JobHandIn (int);
void            JobHandout (int, int, int, ensmbr_struct *, double *, int,
    int);
void            JobRecv (int *, int *, int *, double *, int);
void            LandSfcTmpOper (obs_struct *, var_struct *, pihm_struct);
void            MapVar (var_struct *, int, int);
void            Mbr2Cal (calib_struct *, const double *);
void            PauseParal (int);
void            ReadEnKF (enkf_struct);
void            ReadFcst (enkf_struct, obs_struct, double *);
void            ReadObs (int, char *, double *, double *);
void            ReadVar (char *, enkf_struct, int);
void            Perturb (enkf_struct, char *);
double          Randn ();
void            PIHMParal (int, int, char *);
void            PrintEnKFStatus (int, int);
void            UpdAnlys (enkf_struct, double, double, double *);
void            WriteCalFile (enkf_struct, char *);
void            WriteEnKFOut (char *, enkf_struct, char *, int);
void            WritePara (char *, int, int, int);
void            WriteParamOutput (int, enkf_struct, int, char *);
#endif

#ifdef _CYCLES_
#define Cycles_printf   PIHMprintf
#define Cycles_exit     PIHMexit
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
int             IsOperationToday (int, int, void *, int, int *, int);
void            ApplyFertilizer (fixfert_struct *, soil_struct *,
    residue_struct *);
void            UpdateOperationStatus (void *, int, int);
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
    double, double *, const double *, double, double, double *, double *);
void            River2RiverSolTrnsp (river_struct *, const river_struct *,
    double *, double, double, double, double *, double *);
void            InitCropSV (crop_struct *);
#endif

#ifdef _BGC_
void            AnnualRates (const epconst_struct *, epvar_struct *);
void            BGCSpinup (char *, pihm_struct, char *);
void            CanopyCond (const epconst_struct *, const daily_struct *,
    pstate_struct *, const soil_struct *, epvar_struct *);
void            CheckCarbonBalance (cstate_struct *, double *, int);
void            CheckNitrogenBalance (nstate_struct *, double *, int);
void            CSummary (cflux_struct *, cstate_struct *, summary_struct *);
void            DailyAllocation (cflux_struct *, cstate_struct *,
    nflux_struct *, nstate_struct *, epconst_struct *, epvar_struct *,
    ntemp_struct *, const int);
void            DailyBgc (pihm_struct, int, int, int);
void            DailyCarbonStateUpdate (cflux_struct *, cstate_struct *, int,
    int, int);
void            DailyNitrogenStateUpdate (nflux_struct *, nstate_struct *,
    int alloc, int woody, int evergreen);
void            DayMet (const stor_struct *, daily_struct *, int);
void            Decomp (double, const epconst_struct *, epvar_struct *,
    cstate_struct *, cflux_struct *, nstate_struct *, nflux_struct *,
    ntemp_struct *);
void            FRootLitFall (const epconst_struct *, double, cflux_struct *,
    nflux_struct *);
void            FirstDay (const epconst_struct *, const cinit_struct *,
    epvar_struct *, cstate_struct *, nstate_struct *);
double          GetCO2 (tsdata_struct, int);
double          GetNdep (tsdata_struct, int);
void            GrowthResp (epconst_struct *, cflux_struct *);
void            InitBGC (elem_struct *, int, river_struct *, int,
    const epctbl_struct *, const ctrl_struct *);
void            InitBGCVar (elem_struct *, int, river_struct *, int,
    cinit_struct, cstate_struct, nstate_struct, char *, int);
void            InitElemStor (stor_struct *, int, int);
void            InitRiverStor (river_stor_struct *, int, int);
void            LeafLitFall (const epconst_struct *, double, cflux_struct *,
    nflux_struct *);
void            MaintResp (const cstate_struct *, const nstate_struct *,
    const epconst_struct *, const daily_struct *, cflux_struct *,
    epvar_struct *);
void            MakeZeroFluxStruct (cflux_struct *, nflux_struct *);
void            Mortality (const epconst_struct *, cstate_struct *,
    cflux_struct *, nstate_struct *, nflux_struct *);
void            NLeaching (elem_struct *, int, river_struct *, int);
void            Phenology (const epconst_struct *, const daily_struct *,
    phenology_struct *, epvar_struct *, cstate_struct *, cflux_struct *,
    nstate_struct *, nflux_struct *nf);
void            Photosynthesis (psn_struct *);
void            PrecisionControl (cstate_struct *cs, nstate_struct *ns);
void            RadTrans (const cstate_struct *, const daily_struct *,
    eflux_struct *, pstate_struct *, const epconst_struct *, epvar_struct *);
void            ReadAnnFile (tsdata_struct *, char *);
void            ReadBGC (char *, ctrl_struct *, co2control_struct *,
    ndepcontrol_struct *, char *, char *);
void            ReadEPC (epctbl_struct *);
void            RestartInput (cstate_struct *, nstate_struct *,
    epvar_struct *, bgcic_struct *);
void            RestartOutput (cstate_struct *, nstate_struct *,
    epvar_struct *, bgcic_struct *);
void            RiverDayMet (const river_stor_struct *, river_daily_struct *,
    int);
void            Save2Stor (pihm_struct, int, int, int);
void            SoilPsi (const soil_struct *, double, double *);
void            TotalPhotosynthesis (const epconst_struct *,
    const daily_struct *, const pstate_struct *, epvar_struct *,
    cflux_struct *, psn_struct *, psn_struct *);
void            WriteBGCIC (char *, elem_struct *, int, river_struct *, int);
void            ZeroSrcSnk (cstate_struct *, nstate_struct *,
    summary_struct *);
#endif
#endif
