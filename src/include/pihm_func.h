#ifndef PIHM_FUNC_HEADER
#define PIHM_FUNC_HEADER

#define _ARITH_

// State variables
#define SURF(i)                 (i)
#define UNSAT(i)                (i + nelem)
#define GW(i)                   (i + 2 * nelem)
#define RIVER(i)                (i + 3 * nelem)

#if defined(_DGW_)
# define UNSAT_GEOL(i)          (i + 3 * nelem + nriver)
# define GW_GEOL(i)             (i + 4 * nelem + nriver)
#endif

#if defined(_BGC_) || defined(_CYCLES_)
# if defined(_DGW_)
#  define SOLUTE_SOIL(i, j)     ((i) * nsolute + j + 5 * nelem + nriver)
#  define SOLUTE_RIVER(i, j)    ((i) * nsolute + j + (5 + nsolute) * nelem + nriver)
#  define SOLUTE_GEOL(i, j)     ((i) * nsolute + j + (5 + nsolute) * nelem + (1 + nsolute) * nriver)
# else
#  define SOLUTE_SOIL(i, j)     ((i) * nsolute + j + 3 * nelem + nriver)
#  define SOLUTE_RIVER(i, j)    ((i) * nsolute + j + (3 + nsolute) * nelem + nriver)
# endif
#endif

#define AvgElev(...)            _WsAreaElev(WS_ZMAX, __VA_ARGS__)
#define AvgZmin(...)            _WsAreaElev(WS_ZMIN, __VA_ARGS__)
#define TotalArea(...)          _WsAreaElev(WS_AREA, __VA_ARGS__)

// CVode functions
#if defined(_CVODE_OMP)
# define N_VNew(N, sunctx)      N_VNew_OpenMP(N, nthreads, sunctx)
# define NV_DATA                NV_DATA_OMP
# define NV_Ith                 NV_Ith_OMP
#else
# define N_VNew(N, sunctx)      N_VNew_Serial(N, sunctx)
# define NV_DATA                NV_DATA_S
# define NV_Ith                 NV_Ith_S
#endif

// PIHM system function
#define pihm_error(...)         _error(__FILE__, __LINE__, __FUNCTION__, debug_mode, __VA_ARGS__)
#define pihm_exit(...)          _custom_exit(__FILE__, __LINE__, __FUNCTION__, debug_mode,  __VA_ARGS__)
#define pihm_printf(...)        _custom_printf(verbose_mode, __VA_ARGS__)
#define pihm_fopen              _custom_fopen
#if defined(_WIN32) || defined(_WIN64)
# define pihm_mkdir(path)       _mkdir((path))
# define pihm_access(path, amode) _access((path), (amode))
#else
# define pihm_mkdir(path)       mkdir(path, 0755)
# define pihm_access(path, amode) access((path), (amode))
#endif
#if defined(_MSC_VER)
# define timegm                 _mkgmtime
# define strcasecmp             _stricmp
# define strncasecmp            _strnicmp
#endif

#if defined(_CYCLES_)
#define cycles_error            pihm_error
#define cycles_exit             pihm_exit
#define cycles_fopen            pihm_fopen
#define cycles_printf           pihm_printf
#endif

#define MIN(x, y)               (((x) < (y)) ? (x) : (y))
#define MAX(x, y)               (((x) > (y)) ? (x) : (y))

int             feenableexcept(int);

// Function Declarations
void            _InitLc(const lctbl_struct *, const calib_struct *, elem_struct *);
double          _WsAreaElev(int, const elem_struct *);
void            AdjCVodeMaxStep(cvode_struct *, ctrl_struct *);
void            ApplyBc(int, forc_struct *, elem_struct [], river_struct []);
void            ApplyElemBc(int, forc_struct *, elem_struct []);
#if defined(_NOAH_)
void            ApplyForcing(int, int, const siteinfo_struct *, forc_struct *, elem_struct []);
#else
void            ApplyForcing(int, forc_struct *, elem_struct []);
#endif
#if defined(_BGC_) || defined(_CYCLES_)
void            ApplyLai(elem_struct []);
#else
void            ApplyLai(int, forc_struct *, elem_struct []);
#endif
#if defined(_NOAH_)
void            ApplyMeteoForcing(int, int, const siteinfo_struct *, forc_struct *, elem_struct []);
#else
void            ApplyMeteoForcing(int, forc_struct *, elem_struct []);
#endif
void            ApplyRiverBc(int, forc_struct *, river_struct []);
double          AvgKv(double, double, const soil_struct *);
double          AvgH(double, double, double);
double          AvgHsurf(double, double, double);
void            BackupInput(const char [], const filename_struct *);
void            BoundFluxElem(int, int, const topo_struct *, const soil_struct *, const bc_struct *,
    const wstate_struct *, wflux_struct *);
double          BoundFluxRiver(int, const river_topo_struct *, const shp_struct *, const matl_struct *,
    const river_bc_struct *, const river_wstate_struct *);
void            CalcModelSteps(ctrl_struct *);
double          ChannelFlowElemToRiver(double, double, const river_struct *, elem_struct *);
double          ChannelFlowRiverToRiver(const river_struct *, const river_struct *);
void            CheckCVodeFlag(int);
int             CheckHeader(const char [], int , ...);
#if defined(_BGC_)
int             CheckSteadyState(int, int, int, double, const elem_struct []);
#else
int             CheckSteadyState(int, int, double, const elem_struct []);
#endif
void            CorrectElev(const river_struct [], elem_struct []);
void            CreateOutputDir(char []);
double          DhByDl(const double [], const double [], const double []);
double          EffKh(double, const soil_struct *);
double          EffKinf(double, double, double, double, double, const soil_struct *);
double          EffKv(const soil_struct *, double, int);
void            EtUptake(elem_struct []);
double          FieldCapacity(double, double, double, double);
void            FreeAtttbl(atttbl_struct *);
void            FreeCtrl(ctrl_struct *);
void            FreeForc(forc_struct *);
void            FreeLctbl(lctbl_struct *);
void            FreeMatltbl(matltbl_struct *);
void            FreeMeshtbl(meshtbl_struct *);
void            FreeMem(pihm_struct *);
void            FreeRivtbl(rivtbl_struct *);
void            FreeShptbl(shptbl_struct *);
void            FreeSoiltbl(soiltbl_struct *);
void            FrictionSlope(const elem_struct [], const river_struct [],
    double [], double []);
void            Hydrol(const ctrl_struct *, elem_struct [], river_struct []);
double          Infil(double, const topo_struct *, const soil_struct *, const wstate_struct *, const wstate_struct *,
    const wflux_struct *);
void            InitEFlux(eflux_struct *);
void            InitEState(estate_struct *);
void            InitForcing(const calib_struct *, forc_struct *, elem_struct []);
void            Initialize(pihm_struct *, cvode_struct *);
void            InitLc(const lctbl_struct *, const calib_struct *, elem_struct []);
void            InitMesh(const meshtbl_struct *, elem_struct []);
void            InitOutputFiles(const char [], int, int, print_struct *);
void            InitPrintCtrl(const char [], const char [], int, int, int, varctrl_struct *);
void            InitRiver(const meshtbl_struct *, const rivtbl_struct *, const shptbl_struct *, const matltbl_struct *,
    const calib_struct *, elem_struct [], river_struct []);
void            InitRiverWFlux(river_wflux_struct *);
void            InitRiverWState(river_wstate_struct *);
#if defined(_NOAH_)
void            InitSoil(const soiltbl_struct *, const noahtbl_struct *, const calib_struct *, elem_struct []);
#else
void            InitSoil(const soiltbl_struct *, const calib_struct *, elem_struct []);
#endif
void            InitSurfL(const meshtbl_struct *, elem_struct []);
void            InitTopo(const meshtbl_struct *, elem_struct []);
void            InitVar(elem_struct [], river_struct [], cvode_struct *);
void            InitWbFile(char *, char *, FILE *);
void            InitWFlux(wflux_struct *);
void            InitWState(wstate_struct *);
void            IntcpSnowEt(int, double, const calib_struct *, elem_struct []);
void            IntrplForcing(int, int, int, int, tsdata_struct *);
double          KrFunc(double, double);
void            LateralFlow(const river_struct [], elem_struct []);
#if defined(_CYCLES_)
void            MapOutput(const char [], const int [], const crop_struct [], const elem_struct [],
    const river_struct [], print_struct *);
#else
void            MapOutput(const char [], const int [], const elem_struct [], const river_struct [],  print_struct *);
#endif
#if defined(_DGW_)
void            AdjustFluxes(double, double, const soil_struct *, const soil_struct *, const wstate_struct *,
    const wstate_struct *, double *, wflux_struct *);
#else
void            AdjustFluxes(double, double, const soil_struct *, const wstate_struct *, const wstate_struct *,
    double *, wflux_struct *);
#endif
double          MonthlyLai(int, int);
double          MonthlyMf(int);
double          MonthlyRl(int, int);
int             NumStateVar(void);
int             Ode(sunrealtype, N_Vector, N_Vector, void *);
double          OutletFlux(int, const river_topo_struct *, const shp_struct *, const matl_struct *,
    const river_bc_struct *, const river_wstate_struct *);
double          OverLandFlow(double, double, double, double, double);
double          OvlFlowElemToElem(int, double, const elem_struct *, const elem_struct *);
double          OvlFlowElemToRiver(const river_struct *, elem_struct *);
void            ParseCmdLineParam(int, char *[], char []);
void            PIHM(double, pihm_struct *, cvode_struct *);
pihm_t_struct   PIHMTime(int);
void            PrintCVodeFinalStats(cvode_struct *);
void            PrintData(int, int, int, int, varctrl_struct *);
void            PrintInit(const char [], int, int, int, int, const elem_struct [], const river_struct []);
int             PrintNow(int, int, pihm_t_struct);
void            PrintPerf(int, int, double, double, double, FILE *, cvode_struct *);
void            PrintWaterBalance(int, int, int, const elem_struct [], const river_struct [], FILE *);
void            ProgressBar(double);
double          Psi(double, double, double);
double          PtfAlpha(double, double, double, double, int);
double          PtfBeta(double, double, double, double, int);
double          PtfKv(double, double, double, double, int);
double          PtfThetar(double, double);
double          PtfThetas(double, double, double, double, int);
double          Qtz(int);
void            ReadAlloc(pihm_struct *);
void            ReadAtt(const char [], atttbl_struct *);
void            ReadBc(const char [], const atttbl_struct *, forc_struct *);
void            ReadCalib(const char [], calib_struct *);
void            ReadMeteo(const char [], forc_struct *);
void            ReadIc(const char [], elem_struct [], river_struct []);
int             ReadKeyword(const char [], const char [], char, const char [], int, void *);
void            ReadLai(const char [], const atttbl_struct *, forc_struct *);
void            ReadLc(const char [], lctbl_struct *);
void            ReadMesh(const char [], meshtbl_struct *);
void            ReadPara(const char [], ctrl_struct *);
int             ReadPrintCtrl(const char [], const char [], const char [], int);
void            ReadRiver(const char [], rivtbl_struct *, shptbl_struct *, matltbl_struct *, forc_struct *);
void            ReadSoil(const char [], soiltbl_struct *);
int             ReadTs(const char [], int, int *, double *);
double          Recharge(const soil_struct *, const wstate_struct *, const wflux_struct *);
double          RiverCrossSectArea(int, double, double);
double          RiverEqWid(int, double, double);
void            RiverFlow(elem_struct [], river_struct []);
double          RiverPerim(int, double, double);
void            RiverToElem(river_struct *, elem_struct *, elem_struct *);
int             roundi(double);
#if defined(_OPENMP)
void            RunTime(double, double *, double *);
#else
void            RunTime (clock_t, double *, double *);
#endif
void            RelaxIc(elem_struct [], river_struct []);
void            SetCVodeParam(pihm_struct *, cvode_struct *);
int             SoilTex(double, double);
void            SolveCVode(double, const ctrl_struct *, int *, cvode_struct *);
void            Spinup(pihm_struct *, cvode_struct *);
void            StartupScreen(void);
int             StrTime(const char []);
double          SubsurfFlow(int, const elem_struct *, const elem_struct *);
void            UpdateVar(double, elem_struct [], river_struct [], cvode_struct *);
double          SurfH(double);
void            UpdatePrintVar(int, int, varctrl_struct *);
void            UpdPrintVarT(varctrl_struct *, int);
void            VerticalFlow(double, elem_struct []);
double          WiltingPoint(double, double, double, double);
void            WriteMetadata(const char []);

// DGW functions
#if defined(_DGW_)
double          DeepBoundFluxElem(int, int, const topo_struct *, const soil_struct *, const bc_struct *,
    const wstate_struct *);
double          DeepFlowElemToElem(double, double, const elem_struct *, const elem_struct *);
double          GeolInfil(const topo_struct *, const soil_struct *, const soil_struct *, const wstate_struct *);
double          GeolRecharge(const soil_struct *, const wstate_struct *, const wflux_struct *);
void            FreeGeoltbl(geoltbl_struct *);
void            InitGeol(const geoltbl_struct *, const calib_struct *, elem_struct []);
void            ReadBedrock(const char [], meshtbl_struct *, atttbl_struct *, ctrl_struct *);
void            ReadGeol(const char *, geoltbl_struct *);
#endif

// Noah functions
#if defined(_NOAH_)
void            AdjustSmcProfile(double, const double [], const soil_struct *, const phystate_struct *, wstate_struct *,
    wflux_struct *);
void            AlCalc(int, double, phystate_struct *);
void            CalcLateralFlux(const phystate_struct *, wflux_struct *);
void            CalcSlopeAspect(const meshtbl_struct *, elem_struct []);
void            CalHum(phystate_struct *, estate_struct *);
# if defined(_CYCLES_)
void            CanRes(const estate_struct *, phystate_struct *);
# else
void            CanRes(const soil_struct *, const epconst_struct *, const wstate_struct *, const estate_struct *,
    const eflux_struct *, phystate_struct *);
# endif
double          CSnow(double);
void            DefineSoilDepths(int, double, const double [], int *, double [],
    double []);
void            DEvap(const soil_struct *, const lc_struct *, const phystate_struct *, const wstate_struct *,
    wflux_struct *);
# if defined(_CYCLES_)
void            Evapo(const soil_struct *, const lc_struct *, const weather_struct *, const phystate_struct *,
    const estate_struct *es, const cstate_struct *, crop_struct [], wstate_struct *, wflux_struct *);
# else
void            Evapo(double, const soil_struct *, const lc_struct *, const phystate_struct *, const wstate_struct *,
    wflux_struct *);
# endif
int             FindLayer(double, int, const double []);
int             FindWaterTable(int, double, const double [], double []);
double          FrozRain(double, double);
double          GwTranspFrac(int, int, double, const double []);
void            HRT(double, double, double, double, const soil_struct *, const lc_struct *, const phystate_struct *,
    const estate_struct *, double [], double [], double [], double [], wstate_struct *);
void            HStep(int, double, double [], double [], double [], double [], estate_struct *);
void            IcePac(int, double, double, double, double, const soil_struct *, const lc_struct *, phystate_struct *,
    wstate_struct *, wflux_struct *, estate_struct *, eflux_struct *);
void            InitLsm(const char [], const ctrl_struct *, const noahtbl_struct *, const calib_struct *,
    elem_struct []);
double          Mod(double, double);
void            Noah(double, const lctbl_struct *, const calib_struct *, elem_struct []);
void            NoahHydrol(double, elem_struct []);
# if defined(_CYCLES_)
void            NoPac(double, double, const soil_struct *, const lc_struct *, const weather_struct *,
    const cstate_struct *, crop_struct [], phystate_struct *, wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *);
# else
void            NoPac(double, double, const soil_struct *, const lc_struct *, phystate_struct *, wstate_struct *,
    wflux_struct *, estate_struct *, eflux_struct *);
# endif
void            PcpDrp(double, double, const lc_struct *, wstate_struct *, wflux_struct *);
void            Penman(int, int, double, const estate_struct *, double *, phystate_struct *, wflux_struct *,
    eflux_struct *);
void            PenmanGlacial(int, int, double, const estate_struct *, double *, phystate_struct *, wflux_struct *,
    eflux_struct *);
double          Pslhs(double);
double          Pslhu(double);
double          Pslms(double);
double          Pslmu(double);
double          Psphs(double);
double          Psphu(double);
double          Pspms(double);
double          Pspmu(double);
void            ReadGlacierIce(const char [], double[]);
void            ReadLsm(const char [], ctrl_struct *, siteinfo_struct *, noahtbl_struct *);
void            ReadRad(const char [], forc_struct *);
void            RootDist(int, int, const double [], double []);
void            Rosr12(int, const double [], const double [], const double [], double [], double [], double []);
void            SfcDifOff(int, double, double, const lc_struct *, phystate_struct *);
# if defined(_CYCLES_)
void            SFlx(double, const weather_struct *, const cstate_struct *, soil_struct *, lc_struct *, crop_struct [],
    phystate_struct *, wstate_struct *, wflux_struct *, estate_struct *, eflux_struct *);
# else
void            SFlx(double, soil_struct *, lc_struct *, epconst_struct *, phystate_struct *, wstate_struct *,
    wflux_struct *, estate_struct *, eflux_struct *);
# endif
void            SFlxGlacial(double, soil_struct *, lc_struct *, phystate_struct *, wstate_struct *, wflux_struct *,
    estate_struct *, eflux_struct *);
void            ShFlx(double, double, double, double, const soil_struct *, const lc_struct *, const phystate_struct *,
    wstate_struct *, estate_struct *);
void            SmFlx(double, const soil_struct *, phystate_struct *, wstate_struct *, wflux_struct *);
double          SnFrac(double, double, double);
void            SnkSrc(int, double, double, double, double, const double [], const soil_struct *, double *, double *);
# if defined(_CYCLES_)
void            SnoPac(int, double, double, double, double, const soil_struct *, const lc_struct *,
    const weather_struct *, const cstate_struct *, crop_struct [], phystate_struct *, wstate_struct *,
    wflux_struct *, estate_struct *, eflux_struct *);
# else
void            SnoPac(int, double, double, double, double, const soil_struct *, const lc_struct *, phystate_struct *,
    wstate_struct *, wflux_struct *, estate_struct *, eflux_struct *);
# endif
void            SnowNew(double, const estate_struct *, phystate_struct *);
void            SnowPack(double, double, double, double, double *, double *);
double          Snowz0(double, double, double);
void            SRT(const soil_struct *, double [], double [], double [], double [], double [], phystate_struct *,
    wstate_struct *, wflux_struct *);
void            SStep(double, const soil_struct *, double [], double [], double [], double [], double [],
    phystate_struct *, wstate_struct *, wflux_struct *);
void            SunPos(int, const siteinfo_struct *, spa_data *);
double          TBnd(int, int, double, double, double, const double []);
double          TDfCnd(double, double, double, double, double);
double          TmpAvg(int, double, double, double, const double []);
double          TopoRadn(double, double, double, double, const topo_struct *);
void            Transp(const soil_struct *, const lc_struct *, const phystate_struct *, const wstate_struct *,
    wflux_struct *);
void            WDfCnd(double, double, const soil_struct *, double *, double *);
#endif

#if defined(_DAILY_)
void            DailyVar(int, int, elem_struct []);
void            InitDailyStruct(elem_struct []);
#endif

#if defined(_BGC_) || defined(_CYCLES_)
void            SetAbsTolArray(double, double, N_Vector);
#else
void            SetAbsTolArray(double, N_Vector);
#endif

#if defined(_BGC_) || defined(_CYCLES_)
double          Advection(double, double, double);
void            InitSolute(elem_struct []);
void            RiverElemSoluteFlow(int, int, elem_struct *, river_struct *);
void            SoluteTranspt(elem_struct [], river_struct []);
#endif

#if defined(_BGC_) || defined(_CYCLES_)
void            ReadAnnualFile(const char [], tsdata_struct *);
double          GetCO2(int, tsdata_struct *);
#endif

#if defined(_BGC_)
void            BackgroundLitterfall(const epconst_struct *, const cstate_struct *, epvar_struct *, cflux_struct *,
    nflux_struct *);
void            CanopyCond(const soil_struct *, const epconst_struct *, const daily_struct *, const phystate_struct *,
    const eflux_struct *, epvar_struct *);
void            CheckCarbonBalance(const cstate_struct *, double *);
void            CheckNitrogenBalance(const nstate_struct *, double *);
void            CSummary(const cstate_struct *, const cflux_struct *, summary_struct *);
void            DailyAllocation(const epconst_struct *, const cstate_struct *, const nstate_struct *, epvar_struct *,
    cflux_struct *, nflux_struct *, ntemp_struct *);
void            DailyBgc(int, pihm_struct *);
void            DailyCarbonStateUpdate(int, int, int, cstate_struct *, cflux_struct *);
void            DailyNitrogenStateUpdate(int, int, int, nstate_struct *, nflux_struct *, solute_struct *);
void            Decomp(double, const epconst_struct *, const cstate_struct *, const nstate_struct *, epvar_struct *,
    cflux_struct *, nflux_struct *, ntemp_struct *);
void            EvergreenPhenology(const epconst_struct *, const cstate_struct *, epvar_struct *);
void            FirstDay(const cninit_struct *, elem_struct [], river_struct []);
void            FRootLitFall(double, const epconst_struct *, cflux_struct *, nflux_struct *);
void            FreeEpctbl(epctbl_struct *);
double          GetNdep(int, tsdata_struct *);
void            GrowthResp(const epconst_struct *, cflux_struct *);
void            InitBgc(const epctbl_struct *, const calib_struct *, elem_struct []);
void            InitBgcVar(elem_struct [], river_struct [], cvode_struct *);
void            LeafLitFall(double, const epconst_struct *, cflux_struct *, nflux_struct *);
void            LivewoodTurnover(const epconst_struct *, const cstate_struct *, const nstate_struct *, epvar_struct *,
    cflux_struct *, nflux_struct *);
void            MaintResp(const epconst_struct *, const daily_struct *, const cstate_struct *, const nstate_struct *,
    epvar_struct *, cflux_struct *);
void            MakeZeroFluxStruct(cflux_struct *, nflux_struct *);
void            Mortality(const epconst_struct *, cstate_struct *, cflux_struct *, nstate_struct *, nflux_struct *);
void            OffsetLitterfall(const epconst_struct *, epvar_struct *, const cstate_struct *, cflux_struct *,
    nflux_struct *);
void            OnsetGrowth(const epconst_struct *, const epvar_struct *, const cstate_struct *, const nstate_struct *,
    cflux_struct *, nflux_struct *);
void            Phenology(const epconst_struct *, const daily_struct *, const cstate_struct *, const nstate_struct *,
    epvar_struct *, cflux_struct *, nflux_struct *);
void            Photosynthesis(psn_struct *);
void            PrecisionControl(cstate_struct *cs, nstate_struct *ns);
void            RadTrans(const cstate_struct *, const epconst_struct *, const daily_struct *, epvar_struct *,
    phystate_struct *, eflux_struct *);
void            ReadBgc(const char [], char [], char [], ctrl_struct *, co2control_struct *, ndepcontrol_struct *,
    cninit_struct *);
void            ReadBgcIc(const char [], elem_struct [], river_struct []);
void            ReadEpc(epctbl_struct *);
void            ResetSpinupStat(elem_struct []);
void            RestartInput(const bgcic_struct *, epvar_struct *, cstate_struct *, nstate_struct *);
void            RestartOutput(const epvar_struct *, const cstate_struct *, const nstate_struct *, bgcic_struct *);
void            SeasonDecidPhenology(const epconst_struct *, const daily_struct *, epvar_struct *);
void            SoilPsi(const soil_struct *, double, double *);
void            SoluteConc(elem_struct [], river_struct []);
void            TotalPhotosynthesis(const epconst_struct *, const daily_struct *, const phystate_struct *,
    epvar_struct *, cflux_struct *, psn_struct *, psn_struct *);
void            WriteBgcIc(const char [], elem_struct [], river_struct []);
void            ZeroSrcSnk(cstate_struct *, nstate_struct *, summary_struct *, solute_struct *);
#endif

#if defined(_CYCLES_)
double          AdjustClipThld(double, double);
void            AdjustThermalTime(crop_struct *);
double          Aeration(double);
double          AirMolarDensity(double, double);
void            ApplyDailyMeteoForcing(int, int, const siteinfo_struct *, forc_struct *, elem_struct []);
void            ApplyFert(const fert_struct *fixed_fert, const phystate_struct *phys, cstate_struct *cs,
    nstate_struct *ns, nflux_struct *nf);
void            AutoIrrig(int, const crop_struct [], const airrig_struct [], const soil_struct *, const wstate_struct *,
    const phystate_struct *, wflux_struct *);
double          BoundLayerCond(double, double, double, double);
double          BulkDensity(double, double, double);
void            BurnResidue(int, int, const tillage_struct *, wstate_struct *, cstate_struct *, nstate_struct *);
void            CalSnkSrc(int, const nflux_struct *, solute_struct []);
double          CNDestiny(double, double);
double          CO2FuncGrowth(int, double);
void            ColdDamage(int, int, const weather_struct *, crop_struct [], wstate_struct *, cstate_struct *,
    nstate_struct *, nflux_struct *, phystate_struct *);
double          ColdDamageFunc(double, double, double);
double          CommRadIntcp(const crop_struct []);
double          CommTotRadIntcp(const crop_struct []);
double          CommTransp(const crop_struct []);
void            ComposFactor(const soil_struct *, const wstate_struct *, const estate_struct *, phystate_struct *);
void            ComputeColdDamage(const weather_struct *, const phystate_struct *, crop_struct *, wstate_struct *,
    cstate_struct *, nstate_struct *);
double          ComputeHarvestIndex(const crop_struct *);
void            ConcWeight(const soil_struct *, const wstate_struct *, const phystate_struct *, double []);
int             CondPlant(int, int, const soil_struct *, const wstate_struct *, const estate_struct *,
    const phystate_struct *, plant_struct *);
double          CropGrowth(double, const soil_struct *, const weather_struct *, const wstate_struct *,
    const phystate_struct *, crop_struct *, cflux_struct *);
void            CropHarvest(int, int, crop_struct [], wstate_struct *, cstate_struct *, nstate_struct *, nflux_struct *,
    phystate_struct *);
void            CropNConc(double, const crop_struct *, double *, double *, double *, double *, double *, double *,
    double *);
void            CropNDemand(double, double, const crop_struct *, double *, double *, double *);
void            CropNStress(double, double, double, crop_struct *);
void            CropNUptake(int, const double [], const double [], const double [], const double [], const double [],
    const phystate_struct *, double [], double [], crop_struct [], nstate_struct *, nflux_struct *);
void            CropProcesses(int, crop_struct [], const soil_struct *, const weather_struct *, wstate_struct *,
    wflux_struct *, cstate_struct *, cflux_struct *, nstate_struct *, nflux_struct *, phystate_struct *);
void            CropStage(int, crop_struct []);
void            Cycles(int, const co2control_struct *, forc_struct *, elem_struct []);
void            DailyOper(int, int, int, weather_struct *, mgmt_struct *, crop_struct [], soil_struct *,
    wstate_struct *, wflux_struct *, estate_struct *, cstate_struct *, cflux_struct *, nstate_struct *, nflux_struct *,
    phystate_struct *);
void            Denitrification(const soil_struct *, const wstate_struct *, const cflux_struct *,
    const phystate_struct *, nstate_struct *, nflux_struct *);
double          DepthLimitToEvap(double);
void            DistRootDetritus(double, double, double, double, const crop_struct *, const phystate_struct *,
    cstate_struct *, nstate_struct *);
int             Doy(int, int, int);
void            Doy2Date(int, int, int *, int *);
void            EndRotation(elem_struct []);
void            ExecuteTillage(const tillage_struct *, const phystate_struct *, double [], soil_struct *,
    wstate_struct *, cstate_struct *, nstate_struct *);
void            FieldOper(int, int, mgmt_struct *, crop_struct [], soil_struct *, wstate_struct *, wflux_struct *,
    estate_struct *, cstate_struct *, nstate_struct *, nflux_struct *, phystate_struct *);
int             FindCrop(const char [], const crop_struct []);
double          FindIrrigVolume(int, double, const soil_struct *, const wstate_struct *, const wflux_struct *,
    const phystate_struct *);
int             FinalHarvestDate(int, double, double, double);
void            FirstDay(const soiltbl_struct *, const ctrl_struct *, elem_struct []);
void            FirstDOY(int, const cstate_struct *, mgmt_struct *);
void            FixedHarvest(int, int, const tillage_struct *, const phystate_struct *, crop_struct [], wstate_struct *,
    cstate_struct *, nstate_struct *, nflux_struct *);
void            ForageSeedHarvest(int, int, const phystate_struct *, crop_struct *, wstate_struct *, cstate_struct *,
    nstate_struct *, nflux_struct *);
int             ForcedClip(int, crop_struct []);
int             ForcedMaturity(int, int, int, int, int, int);
double          Fraction(double, double, double, double, double);
void            FreeAgtbl(agtbl_struct *);
void            FreeMgmttbl(int, mgmt_struct []);
double          GrainGrowth(const crop_struct *);
void            GrainHarvest(int, int, double, crop_struct *, wstate_struct *, cstate_struct *, nstate_struct *,
    nflux_struct *);
void            GrowingCrop(int, int, int, const soil_struct *, const weather_struct *, crop_struct [], wstate_struct *,
    wflux_struct *, cstate_struct *, cflux_struct *, nstate_struct *, nflux_struct *, phystate_struct *);
void            InitAgVar(elem_struct [], river_struct [], cvode_struct *);
void            InitCropStateVar(crop_struct *);
void            InitCycles(const calib_struct *, const agtbl_struct *, const mgmt_struct [], const crop_struct [],
    const soiltbl_struct *, elem_struct []);
void            InitMgmt(mgmt_struct *);
double          IntegRoot(double, double);
int             InTimeWindow(int, int, int);
int             IsLeapYear(int);
int             IsOperToday(int, int, const soil_struct *, const wstate_struct *, const estate_struct *, int *,
    mgmt_struct *, phystate_struct *);
void            KillCrop(int, int, const phystate_struct *, crop_struct *, wstate_struct *, cstate_struct *,
    nstate_struct *);
void            LateralNFlow(double, const soil_struct *, const wstate_struct *, const phystate_struct *,
    const double [], double, double, double []);
double          LinearEqmConc(double, double, double, double, double);
double          MaxAbgdHumifFactor(double);
double          MaxManureHumifFactor(double);
double          MaxRhizoHumifFactor(double);
double          MaxRootHumifFactor(double);
double          MobileNConc(double, const double [], const soil_struct *, const wstate_struct *,
    const phystate_struct *);
double          Moisture(double);
double          N2OFracNitrif(double);
void            Nitrification(const soil_struct *, const wstate_struct *, const estate_struct *,
    const phystate_struct *, nstate_struct *, nflux_struct *);
double          NMineral(double, double, double, double);
int             NumActiveCrop(const crop_struct []);
void            NXform(const crop_struct [], const weather_struct *, const soil_struct *, const wstate_struct *,
    const estate_struct *, const cstate_struct *, const cflux_struct *, const phystate_struct *, nstate_struct *,
    nflux_struct *);
double          PHFunction(double);
void            Phenology(const weather_struct *, crop_struct []);
void            PlantCrop(int, const plant_struct *, crop_struct *);
void            PotSoluteUptake(double, const double [], const soil_struct *, const wstate_struct *,
    const wflux_struct *, const phystate_struct *, double []);
double          Profile(int, const double []);
void            RadIntcp(double, crop_struct []);
void            ReadCrop(const char [], crop_struct []);
void            ReadCyclesCtrl(const char [], char [], agtbl_struct *, ctrl_struct *, co2control_struct *);
void            ReadCyclesIc(const char [], elem_struct []);
void            ReadMultOper(const agtbl_struct *, mgmt_struct [], crop_struct []);
void            ReadOper(const char [], int, int, mgmt_struct *, crop_struct []);
void            ReadSoilInit(const char [], soiltbl_struct *);
void            ResidueCover(const cstate_struct *, phystate_struct *);
void            ResetCrop(crop_struct *);
void            ResidueEvap(double, double, const crop_struct [], const cstate_struct *, const phystate_struct *,
    wstate_struct *, wflux_struct *);
void            ResidueWetting(double, const cstate_struct *, const phystate_struct *, wstate_struct *, wflux_struct *);
void            RootFrac(const crop_struct *, const phystate_struct *, double []);
double          SatVP(double);
double          ShootBiomassPartn(int, double, double, double);
double          SoilBufferPower(double, double, double);
void            SoilCarbonBalance(const double [], const soil_struct *, const crop_struct [], const phystate_struct *,
    wstate_struct *, cstate_struct *, cflux_struct *, nstate_struct *, nflux_struct *);
int             SoilCond(const plant_struct *, const soil_struct *, const wstate_struct *, const estate_struct *,
    const phystate_struct *);
void            SoilEvap(double, double, const crop_struct [], const soil_struct *, const phystate_struct *,
    wstate_struct *, wflux_struct *);
double          SoilTmpMovingAvg(double, const double []);
double          SoilWaterContent(double, double, double, double);
double          SoilWaterPot(double, double, double, double);
void            SoluteConc(double, elem_struct [], river_struct []);
void            SoluteTransp(double, double, const double [], const double [], const soil_struct *,
    const phystate_struct *, double []);
double          TextureFactor(double);
double          ThermalTime(double, double, double, double);
void            TillageFactor(double, const tillage_struct *, const soil_struct *, const phystate_struct *, double []);
void            TillageFactorSet(int , const double [], double , double []);
double          TmpFunc(double);
double          TmpFuncGrowth(double, double, double, double);
double          TmpLimit(double, double, double);
void            UpdateNProfile(double, const soil_struct *, const wstate_struct *, const nstate_struct *,
    const solute_struct [], double [], double [], phystate_struct *, nflux_struct *);
void            UpdateOperPtr(mgmt_struct *);
double          VolatilDepthFunc(double);
void            Volatilization(const weather_struct *, const crop_struct [], const soil_struct *, const wstate_struct *,
    const cstate_struct *, const phystate_struct *, nstate_struct *, nflux_struct *);
double          VolWCAt33Jkg(double, double, double);
double          VolWCAt1500Jkg(double, double, double);
double          WaterContentLimitToEvap(double, double, double);
void            WaterUptake(double, const soil_struct *, const weather_struct *, const phystate_struct *,
    crop_struct [], wstate_struct *, wflux_struct *);
void            WriteCyclesIc(const char [], const elem_struct []);
void            ZeroFluxes(wflux_struct *, cflux_struct *, nflux_struct *);
void            ZeroHarvest(crop_struct *);
void            NRT(double, double, double [], const soil_struct *, const wstate_struct *, const wstate_struct *,
    const phystate_struct *, double []);
#endif

#endif
