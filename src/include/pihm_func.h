#ifndef PIHM_FUNC_HEADER
#define PIHM_FUNC_HEADER

#define _ARITH_

/* State variables */
#define SURF(i)                 (i)
#define UNSAT(i)                (i + nelem)
#define GW(i)                   (i + 2 * nelem)
#define RIVER(i)                (i + 3 * nelem)

#if defined(_FBR_)
# define FBRUNSAT(i)            (i + 3 * nelem + nriver)
# define FBRGW(i)               (i + 4 * nelem + nriver)
#endif

#if defined(_BGC_) || defined(_CYCLES_) || defined(_RT_)
# if defined(_FBR_)
#  define SOLUTE_SOIL(i, j)     ((i) * nsolute + j + 5 * nelem + nriver)
#  define SOLUTE_RIVER(i, j)    ((i) * nsolute + j + (5 + nsolute) * nelem + nriver)
#  define SOLUTE_GEOL(i, j)     ((i) * nsolute + j + (5 + nsolute) * nelem + (1 + nsolute) * nriver)
# else
#  define SOLUTE_SOIL(i, j)     ((i) * nsolute + j + 3 * nelem + nriver)
#  define SOLUTE_RIVER(i, j)    ((i) * nsolute + j + (3 + nsolute) * nelem + nriver)
# endif
#endif

#if defined(_BGC_) && !defined(_LUMPED_)
# define SURFN(i)               (i + 3 * nelem + 2 * nriver)
# define SMINN(i)               (i + 4 * nelem + 2 * nriver)
# define STREAMN(i)             (i + 5 * nelem + 2 * nriver)
# define RIVBEDN(i)             (i + 5 * nelem + 3 * nriver)
#else
# define LUMPED_SMINN           (3 * nelem + 2 * nriver)
#endif

#define AvgElev(...)            _WsAreaElev(WS_ZMAX, __VA_ARGS__)
#define AvgZmin(...)            _WsAreaElev(WS_ZMIN, __VA_ARGS__)
#define TotalArea(...)          _WsAreaElev(WS_AREA, __VA_ARGS__)

/* CVode functions */
#if defined(_CVODE_OMP)
# define N_VNew(N)              N_VNew_OpenMP(N, nthreads)
# define NV_DATA                NV_DATA_OMP
# define NV_Ith                 NV_Ith_OMP
#else
# define N_VNew(N)              N_VNew_Serial(N)
# define NV_DATA                NV_DATA_S
# define NV_Ith                 NV_Ith_S
#endif

/* PIHM system function */
#define pihm_exit(...)          _custom_exit(__FILE__, __LINE__, __FUNCTION__, debug_mode,  __VA_ARGS__)
#define pihm_printf(...)        _custom_printf(__FILE__, __LINE__, __FUNCTION__, debug_mode, verbose_mode, __VA_ARGS__)
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
#define cycles_exit             pihm_exit
#define cycles_fopen            pihm_fopen
#define cycles_printf           pihm_printf
#endif

#define MIN(x, y)               (((x) < (y)) ? (x) : (y))
#define MAX(x, y)               (((x) > (y)) ? (x) : (y))

int             feenableexcept(int);

/*
 * Function Declarations
 */
void            _InitLc(elem_struct *, const lctbl_struct *,
    const calib_struct *);
double          _WsAreaElev(int, const elem_struct *);
void            AdjCVodeMaxStep(void *, ctrl_struct *);
#if defined(_RT_)
void            ApplyBc(const rttbl_struct *, forc_struct *, elem_struct *,
    river_struct *, int);
#else
void            ApplyBc(forc_struct *, elem_struct *, river_struct *, int);
#endif
#if defined(_RT_)
void            ApplyElemBc(const rttbl_struct *, forc_struct *, elem_struct *,
    int);
#else
void            ApplyElemBc(forc_struct *, elem_struct *, int);
#endif
#if defined(_RT_)
void            ApplyForc(forc_struct *, rttbl_struct *, elem_struct *, int,
    int, const siteinfo_struct *);
#elif defined(_NOAH_)
void            ApplyForc(forc_struct *, elem_struct *, int, int,
    const siteinfo_struct *);
#else
void            ApplyForc(forc_struct *, elem_struct *, int);
#endif
#if defined(_BGC_) || defined(_CYCLES_)
void            ApplyLai(elem_struct []);
#else
void            ApplyLai(forc_struct *, elem_struct *, int);
#endif
#if defined(_NOAH_)
void            ApplyMeteoForc(forc_struct *, elem_struct *, int, int,
    const siteinfo_struct *);
#else
void            ApplyMeteoForc(forc_struct *, elem_struct *, int);
#endif
void            ApplyRiverBc(forc_struct *, river_struct *, int);
double          AvgKv(const soil_struct *, double, double);
double          AvgH(double, double, double);
double          AvgHsurf(double, double, double);
void            BackupInput(const char *, const filename_struct *);
void            BoundFluxElem(int, int, const bc_struct *,
    const wstate_struct *, const topo_struct *, const soil_struct *,
    wflux_struct *);
double          BoundFluxRiver(int, const river_wstate_struct *,
    const river_topo_struct *, const shp_struct *, const matl_struct *,
    const river_bc_struct *bc);
void            CalcModelStep(ctrl_struct *);
double          ChanFlowElemToRiver(double, double, const river_struct *,
    elem_struct *);
double          ChanFlowRiverToRiver(const river_struct *, const river_struct *,
    int);
void            CheckCVodeFlag(int);
#if defined(_BGC_)
int             CheckSteadyState(const elem_struct *, double, int, int, int);
#else
int             CheckSteadyState(const elem_struct *, double, int, int);
#endif
void            CorrElev(elem_struct *, river_struct *);
void            CreateOutputDir(char *);
double          DhByDl(const double *, const double *, const double *);
double          EffKh(const soil_struct *, double);
double          EffKinf(const soil_struct *, double, double, double, double,
    double);
double          EffKv(const soil_struct *, double, int);
void            EtExtract(elem_struct *);
double          FieldCapacity(double, double, double, double);
void            FreeAtttbl(atttbl_struct *);
void            FreeCtrl(ctrl_struct *);
void            FreeForc(forc_struct *);
void            FreeLctbl(lctbl_struct *);
void            FreeMatltbl(matltbl_struct *);
void            FreeMeshtbl(meshtbl_struct *);
void            FreeMem(pihm_struct);
void            FreeRivtbl(rivtbl_struct *);
void            FreeShptbl(shptbl_struct *);
void            FreeSoiltbl(soiltbl_struct *);
void            FrictSlope(const elem_struct *, const river_struct *, double *,
    double *);
void            Hydrol(elem_struct *, river_struct *, const ctrl_struct *);
double          Infil(const wstate_struct *, const wstate_struct *,
    const wflux_struct *, const topo_struct *, const soil_struct *, double);
void            InitEFlux(eflux_struct *);
void            InitEState(estate_struct *);
#if defined(_RT_)
void            InitForc(elem_struct *, forc_struct *, const calib_struct *,
    const rttbl_struct *);
#else
void            InitForc(elem_struct *, forc_struct *, const calib_struct *);
#endif
void            Initialize(pihm_struct, N_Vector, void **);
void            InitLc(elem_struct *, const lctbl_struct *,
    const calib_struct *);
void            InitMesh(elem_struct *, const meshtbl_struct *);
#if defined(_TGM_) && defined(_RT_)
void            InitOutputFile(const char *, int, int, const chemtbl_struct [],
    const rttbl_struct *, print_struct *);
#else
void            InitOutputFile(const char *, int, int, print_struct *);
#endif
void            InitPrtVarCtrl(const char *, const char *, int, int, int,
    varctrl_struct *);
void            InitRiver(river_struct *, elem_struct *, const rivtbl_struct *,
    const shptbl_struct *, const matltbl_struct *, const meshtbl_struct *,
    const calib_struct *);
void            InitRiverWFlux(river_wflux_struct *);
void            InitRiverWState(river_wstate_struct *);
#if defined(_NOAH_)
void            InitSoil(elem_struct *, const soiltbl_struct *,
    const noahtbl_struct *, const calib_struct *);
#else
void            InitSoil(elem_struct *, const soiltbl_struct *,
    const calib_struct *);
#endif
void            InitSurfL(elem_struct *, const meshtbl_struct *);
void            InitTopo(elem_struct *, const meshtbl_struct *);
void            InitVar(elem_struct *, river_struct *, N_Vector);
void            InitWbFile(char *, char *, FILE *);
void            InitWFlux(wflux_struct *);
void            InitWState(wstate_struct *);
void            IntcpSnowEt(int, double, elem_struct *, const calib_struct *);
void            IntrplForc(tsdata_struct *, int, int, int);
double          KrFunc(double, double);
void            LateralFlow(elem_struct *, const river_struct *, int);
#if defined(_CYCLES_)
void            MapOutput(const int *, const crop_struct [],
    const elem_struct *, const river_struct *, const char *, print_struct *);
#elif defined(_RT_)
void            MapOutput(const int *, const chemtbl_struct [],
    const rttbl_struct *, const elem_struct *, const river_struct *,
    const char *, print_struct *);
#else
void            MapOutput(const int *, const elem_struct *,
    const river_struct *, const char *, print_struct *);
#endif
#if defined(_FBR_)
void            MassBalance(const wstate_struct *, const wstate_struct *,
    wflux_struct *, double *, const soil_struct *, const soil_struct *, double,
    double);
#else
void            MassBalance(const wstate_struct *, const wstate_struct *,
    wflux_struct *, double *, const soil_struct *, double, double);
#endif
double          MonthlyLai(int, int);
double          MonthlyMf(int);
double          MonthlyRl(int, int);
int             NumStateVar(void);
int             Ode(realtype, N_Vector, N_Vector, void *);
double          OutletFlux(int, const river_wstate_struct *,
    const river_topo_struct *, const shp_struct *, const matl_struct *,
    const river_bc_struct *);
double          OverLandFlow(double, double, double, double, double);
double          OvlFlowElemToElem(const elem_struct *, const elem_struct *, int,
    double, int);
double          OvlFlowElemToRiver(int, const river_struct *, elem_struct *);
void            ParseCmdLineParam(int, char *[], char *);
void            PIHM(pihm_struct, void *, N_Vector, double);
pihm_t_struct   PIHMTime(int);
void            PrintCVodeFinalStats(void *);
void            PrintData(varctrl_struct *, int, int, int, int);
void            PrintInit(const elem_struct *, const river_struct *,
    const char *, int, int, int, int);
int             PrintNow(int, int, const pihm_t_struct *);
void            PrintPerf(void *, int, int, double, double, double, FILE *);
void            PrintWaterBal(FILE *, int, int, int, const elem_struct *,
    const river_struct *);
void            ProgressBar(double);
double          Psi(double, double, double);
double          PtfAlpha(double, double, double, double, int);
double          PtfBeta(double, double, double, double, int);
double          PtfKv(double, double, double, double, int);
double          PtfThetar(double, double);
double          PtfThetas(double, double, double, double, int);
double          Qtz(int);
void            ReadAlloc(pihm_struct);
void            ReadAtt(const char *, atttbl_struct *);
#if defined(_RT_)
void            ReadBc(const char *, forc_struct *, const atttbl_struct *,
    const rttbl_struct *, const chemtbl_struct []);
#else
void            ReadBc(const char *, forc_struct *, const atttbl_struct *);
#endif
void            ReadCalib(const char *, calib_struct *);
void            ReadForc(const char *, forc_struct *);
void            ReadIc(const char *, elem_struct *, river_struct *);
int             ReadKeyword(const char *, const char *, void *, char,
    const char *, int);
void            ReadLai(const char *, forc_struct *, const atttbl_struct *);
void            ReadLc(const char *, lctbl_struct *);
void            ReadMesh(const char *, meshtbl_struct *);
void            ReadPara(const char *, ctrl_struct *);
int             ReadPrtCtrl(const char *, const char *, const char *, int);
void            ReadRiver(const char *, rivtbl_struct *, shptbl_struct *,
    matltbl_struct *, forc_struct *);
void            ReadSoil(const char *, soiltbl_struct *);
int             ReadTS(const char *, int *, double *, int);
double          Recharge(const wstate_struct *, const wflux_struct *,
    const soil_struct *);
double          RiverCroSectArea(int, double, double);
double          RiverEqWid(int, double, double);
void            RiverFlow(int, int, elem_struct *, river_struct *);
double          RiverPerim(int, double, double);
void            RiverToElem(int, river_struct *, elem_struct *, elem_struct *);
int             roundi(double);
#if defined(_OPENMP)
void            RunTime(double, double *, double *);
#else
void            RunTime (clock_t, double *, double *);
#endif
void            RelaxIc(elem_struct *, river_struct *);
void            SetCVodeParam(pihm_struct, void *, SUNLinearSolver *, N_Vector);
int             SoilTex(double, double);
void            SolveCVode(const ctrl_struct *, double, int *, void *,
    N_Vector);
void            Spinup(pihm_struct, N_Vector, void *, SUNLinearSolver *);
void            StartupScreen(void);
int             StrTime(const char *);
double          SubFlowElemToElem(const elem_struct *, const elem_struct *,
    int);
void            Summary(elem_struct *, river_struct *, N_Vector, double);
double          SurfH(double);
void            UpdPrintVar(varctrl_struct *, int, int);
void            UpdPrintVarT(varctrl_struct *, int);
void            VerticalFlow(elem_struct *, double);
double          WiltingPoint(double, double, double, double);

/*
 * Fractured bedrock functions
 */
#if defined(_FBR_)
double          FbrBoundFluxElem(int, int, const bc_struct *,
    const wstate_struct *, const topo_struct *, const soil_struct *);
double          FbrFlowElemToElem(const elem_struct *, const elem_struct *,
    double, double);
double          FbrInfil(const wstate_struct *, const soil_struct *,
    const soil_struct *, const topo_struct *);
double          FbrRecharge(const wstate_struct *, const wflux_struct *,
    const soil_struct *);
void            FreeGeoltbl(geoltbl_struct *);
void            InitGeol (elem_struct *, const geoltbl_struct *,
        const calib_struct *);
void            ReadBedrock(const char *, atttbl_struct *, meshtbl_struct *,
    ctrl_struct *);
void            ReadGeol(const char *, geoltbl_struct *);
#endif

/*
 * Noah functions
 */
#if defined(_NOAH_)
void            AdjSmProf(const soil_struct *, const phystate_struct *,
    const double *, double, wflux_struct *, wstate_struct *);
void            AlCalc(phystate_struct *, double, int);
void            CalcLatFlx(const phystate_struct *, wflux_struct *);
void            CalcSlopeAspect(elem_struct *, const meshtbl_struct *);
void            CalHum(phystate_struct *, estate_struct *);
# if defined(_CYCLES_)
void            CanRes(const estate_struct *, phystate_struct *);
# else
void            CanRes(const wstate_struct *, const estate_struct *,
    const eflux_struct *, phystate_struct *, const soil_struct *,
    const epconst_struct *);
# endif
double          CSnow(double);
void            DefSldpth(double *, int *, double *, double, const double *,
    int);
void            DEvap(const wstate_struct *, wflux_struct *,
    const phystate_struct *, const lc_struct *, const soil_struct *);
# if defined(_CYCLES_)
void            Evapo(const soil_struct *, const lc_struct *,
    const weather_struct *, const phystate_struct *, const estate_struct *es,
    const cstate_struct *, crop_struct [], wstate_struct *, wflux_struct *);
# else
void            Evapo(const wstate_struct *, wflux_struct *,
    const phystate_struct *, const lc_struct *, const soil_struct *, double);
# endif
int             FindLayer(const double *, int, double);
int             FindWaterTable(const double *, int, double, double *);
double          FrozRain(double, double);
double          GwTransp(double, const double *, int, int);
void            HRT(wstate_struct *, const estate_struct *,
    const phystate_struct *, const lc_struct *, const soil_struct *, double *,
    double, double, double, double, double *, double *, double *);
void            IcePac(wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, phystate_struct *, const lc_struct *, const soil_struct *,
    int, double, double, double, double);
void            InitLsm(elem_struct *, const char [], const ctrl_struct *,
    const noahtbl_struct *, const calib_struct *);
double          Mod(double, double);
void            Noah(elem_struct *, const lctbl_struct *, const calib_struct *,
    double);
void            NoahHydrol(elem_struct *, double);
# if defined(_CYCLES_)
void            NoPac(const soil_struct *, const lc_struct *,
    const weather_struct *, const cstate_struct *, double, double,
    crop_struct [], phystate_struct *, wstate_struct *, wflux_struct *,
    estate_struct *, eflux_struct *);
# else
void            NoPac(wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, phystate_struct *, const lc_struct *, const soil_struct *,
    double, double);
# endif
void            PcpDrp(wstate_struct *, wflux_struct *, const lc_struct *,
    double, double);
void            Penman(wflux_struct *, const estate_struct *, eflux_struct *,
    phystate_struct *, double *, double, int, int);
void            PenmanGlacial(wflux_struct *, const estate_struct *,
    eflux_struct *, phystate_struct *, double *, double, int, int);
double          Pslhs(double);
double          Pslhu(double);
double          Pslms(double);
double          Pslmu(double);
double          Psphs(double);
double          Psphu(double);
double          Pspms(double);
double          Pspmu(double);
void            ReadGlacierIce(const char [], double[]);
void            ReadLsm(const char *, siteinfo_struct *, ctrl_struct *,
    noahtbl_struct *);
void            ReadRad(const char *, forc_struct *);
void            RootDist(const double *, int, int, double *);
void            Rosr12(double *, const double *, const double *, double *,
    const double *, double *, int);
void            SfcDifOff(phystate_struct *, const lc_struct *, double, double,
    int);
# if defined(_CYCLES_)
void            SFlx(double, const weather_struct *, const cstate_struct *,
    soil_struct *, lc_struct *, crop_struct [], phystate_struct *,
    wstate_struct *, wflux_struct *, estate_struct *, eflux_struct *);
# else
void            SFlx(wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, phystate_struct *, lc_struct *, epconst_struct *,
    soil_struct *, double);
# endif
void            SFlxGlacial(wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, phystate_struct *, lc_struct *, soil_struct *, double);
void            ShFlx(wstate_struct *, estate_struct *, const phystate_struct *,
    const lc_struct *, const soil_struct *, double, double, double, double);
void            SmFlx(wstate_struct *, wflux_struct *, phystate_struct *,
    const soil_struct *, double);
double          SnFrac(double, double, double);
void            SnkSrc(double *, double, double, double *,
    const soil_struct *, const double *, double, int, double);
# if defined(_CYCLES_)
void            SnoPac(const soil_struct *, const lc_struct *,
    const weather_struct *, const cstate_struct *, int, double, double, double,
    double, crop_struct [], phystate_struct *, wstate_struct *, wflux_struct *,
    estate_struct *, eflux_struct *);
# else
void            SnoPac(wstate_struct *, wflux_struct *, estate_struct *,
    eflux_struct *, phystate_struct *, const lc_struct *, const soil_struct *,
    int, double, double, double, double);
# endif
void            SnowNew(const estate_struct *, double, phystate_struct *);
void            SnowPack(double, double, double *, double *, double, double);
double          Snowz0(double, double, double);
void            SRT(wstate_struct *, wflux_struct *, phystate_struct *,
    const soil_struct *, double *, double *, double *, double *, double *);
void            SStep(wstate_struct *, wflux_struct *, phystate_struct *,
    const soil_struct *, double *, double *, double *, double *, double *,
    double);
void            SunPos(const siteinfo_struct *, int, spa_data *);
double          TBnd(double, double, const double *, double, int, int);
double          TDfCnd(double, double, double, double, double);
double          TmpAvg(double, double, double, const double *, int);
double          TopoRadn(const topo_struct *, double, double, double, double);
void            Transp(const wstate_struct *, wflux_struct *,
    const phystate_struct *, const lc_struct *, const soil_struct *);
void            WDfCnd(double *, double *, double, double, const soil_struct *);
#endif

#if defined(_DAILY_)
void            DailyVar(int, int, elem_struct *);
void            InitDailyStruct(elem_struct *);
#endif

#if defined(_BGC_) || defined(_CYCLES_) || defined(_RT_)
void            SetAbsTolArray(double, double, N_Vector);
#else
void            SetAbsTolArray(double, N_Vector);
#endif

#if defined(_BGC_) || defined(_CYCLES_) || defined(_RT_)
double          AdvDiffDisp(double, double, double, double, double, double,
                    double, double, double);
void            InitSolute(elem_struct []);
# if defined(_FBR_) && defined(_TGM_)
void            RiverElemSoluteFlow(int, int, int, elem_struct *,
                    river_struct *);
# else
void            RiverElemSoluteFlow(int, int, elem_struct *, river_struct *);
# endif
void            SoluteTranspt(double, double, double, elem_struct [],
                    river_struct []);
#endif

#if defined(_BGC_)
void            BackgroundLitterfall(const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, nflux_struct *);
void            CanopyCond(const epconst_struct *, epvar_struct *,
    const eflux_struct *, const phystate_struct *, const soil_struct *,
    const daily_struct *);
void            CheckCarbonBalance(const cstate_struct *, double *);
void            CheckNitrogenBalance(const nstate_struct *, double *);
void            CSummary(const cflux_struct *, const cstate_struct *,
    summary_struct *);
void            DailyAllocation(cflux_struct *, const cstate_struct *,
    nflux_struct *, const nstate_struct *, const epconst_struct *,
    epvar_struct *, ntemp_struct *);
void            DailyBgc(pihm_struct, int);
void            DailyCarbonStateUpdate(cflux_struct *, cstate_struct *, int,
    int, int);
void            DailyNitrogenStateUpdate(nflux_struct *, nstate_struct *,
    solute_struct *, int, int, int);
void            Decomp(double, const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, const nstate_struct *,
    nflux_struct *, ntemp_struct *);
void            EvergreenPhenology(const epconst_struct *, epvar_struct *,
    const cstate_struct *);
void            FirstDay(elem_struct *, river_struct *, const cninit_struct *);
void            FRootLitFall(const epconst_struct *, double, cflux_struct *,
    nflux_struct *);
void            FreeEpctbl(epctbl_struct *);
double          GetCO2(tsdata_struct *, int);
double          GetNdep(tsdata_struct *, int);
void            GrowthResp(const epconst_struct *, cflux_struct *);
void            InitBgc(elem_struct *, const epctbl_struct *,
    const calib_struct *);
void            InitBgcVar(elem_struct *, river_struct *, N_Vector);
void            LeafLitFall(const epconst_struct *, double, cflux_struct *,
    nflux_struct *);
void            LivewoodTurnover(const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, const nstate_struct *,
    nflux_struct *);
void            MaintResp(const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, const nstate_struct *,
    const daily_struct *);
void            MakeZeroFluxStruct(cflux_struct *, nflux_struct *);
void            Mortality(const epconst_struct *, cstate_struct *,
    cflux_struct *, nstate_struct *, nflux_struct *);
# if defined(_LEACHING_)
void            NLeaching(elem_struct *);
# elif defined(_LUMPED_)
void            NLeachingLumped(elem_struct *, river_struct *);
# else
void            NTransport(elem_struct *, river_struct *);
# endif
void            OffsetLitterfall(const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, nflux_struct *);
void            OnsetGrowth(const epconst_struct *, const epvar_struct *,
    const cstate_struct *, cflux_struct *, const nstate_struct *,
    nflux_struct *);
void            Phenology(const epconst_struct *, epvar_struct *,
    const cstate_struct *, cflux_struct *, const nstate_struct *,
    nflux_struct *, const daily_struct *);
void            Photosynthesis(psn_struct *);
void            PrecisionControl(cstate_struct *cs, nstate_struct *ns);
void            RadTrans(const cstate_struct *, eflux_struct *,
    phystate_struct *, const epconst_struct *, epvar_struct *,
    const daily_struct *);
void            ReadAnnFile(tsdata_struct *, const char *);
void            ReadBgc(const char *, ctrl_struct *, co2control_struct *,
    ndepcontrol_struct *, cninit_struct *, char *, char *);
void            ReadBgcIc(const char *, elem_struct *, river_struct *);
void            ReadEpc(epctbl_struct *);
void            ResetSpinupStat(elem_struct *);
void            RestartInput(cstate_struct *, nstate_struct *,
    epvar_struct *, const bgcic_struct *);
void            RestartOutput(const cstate_struct *, const nstate_struct *,
    const epvar_struct *, bgcic_struct *);
void            SeasonDecidPhenology(const epconst_struct *, epvar_struct *,
    const daily_struct *);
void            SoilPsi(const soil_struct *, double, double *);
void            TotalPhotosynthesis(const epconst_struct *, epvar_struct *,
    const phystate_struct *, cflux_struct *, psn_struct *, psn_struct *,
    const daily_struct *);
void            WriteBgcIc(const char *, elem_struct *, river_struct *);
void            ZeroSrcSnk(cstate_struct *, nstate_struct *, summary_struct *,
    solute_struct *);
#endif

#if defined(_CYCLES_)
double          AdjustClipThld(double, double);
double          Aeration(double);
double          AirMolarDensity(double, double);
void            ApplyDailyMeteoForc(int, int, const siteinfo_struct *,
    forc_struct *, elem_struct []);
void            ApplyFert(const fert_struct *, cstate_struct *, nstate_struct *,
    nflux_struct *);
void            AutoIrrig(int, const crop_struct [], const airrig_struct [],
    const soil_struct *, const wstate_struct *, const phystate_struct *,
    wflux_struct *);
double          BoundLayerCond(double, double, double, double);
double          BulkDensity(double, double, double);
void            CalSnkSrc(int, const nflux_struct *, solute_struct []);
double          CNDestiny(double, double);
double          ColdDamage(double, double, double);
double          CommRadIntcp(const crop_struct []);
double          CommTotRadIntcp(const crop_struct []);
double          CommTransp(const crop_struct []);
void            ComposFactor(const soil_struct *, const wstate_struct *,
    const estate_struct *, phystate_struct *);
void            ComputeColdDamage(const weather_struct *,
    const phystate_struct *, crop_struct *, wstate_struct *, cstate_struct *,
    nstate_struct *);
double          ComputeHarvestIndex(double, double, double, double, double);
int             CondPlant(int, int, const plant_struct *, const soil_struct *,
    const wstate_struct *, const estate_struct *);
double          CropGrowth(double, const soil_struct *, const weather_struct *,
    const wstate_struct *, const phystate_struct *, crop_struct *);
void            CropNConc(double, const crop_struct *, double *, double *,
    double *, double *, double *, double *, double *);
void            CropNDemand(double, double, const crop_struct *, double *,
    double *, double *);
void            CropNStress(double, double, double, crop_struct *);
void            CropNUptake(int, double, double, const double [],
    const double [], const double [], const double [], const double [],
    const phystate_struct *, double [], double [], crop_struct [],
    nstate_struct *, nflux_struct *);
void            CropStage(int, crop_struct []);
void            Cycles(int, elem_struct []);
void            DailyOper(int, int, int, weather_struct *,
    mgmt_struct *, crop_struct [], soil_struct *, wstate_struct *,
    wflux_struct *, estate_struct *, cstate_struct *, cflux_struct *,
    nstate_struct *, nflux_struct *, phystate_struct *);
void            Denitrification(const soil_struct *, const wstate_struct *,
    const cflux_struct *, const phystate_struct *, nstate_struct *,
    nflux_struct *);
double          DepthLimitToEvap(double);
void            DistRootDetritus(double, double, double, double,
    const crop_struct *, const phystate_struct *, cstate_struct *,
    nstate_struct *);
int             Doy(int, int, int);
void            Doy2Date(int, int, int *, int *);
void            ExecuteTillage(const tillage_struct *, const phystate_struct *,
     double [], soil_struct *, wstate_struct *, cstate_struct *,
     nstate_struct *);
void            FieldOper(int, int, mgmt_struct *, crop_struct [],
    soil_struct *, wstate_struct *, wflux_struct *, estate_struct *,
    cstate_struct *, nstate_struct *, nflux_struct *, phystate_struct *);
int             FindCrop(const char [], const crop_struct []);
double          FindIrrigVolume(int, double, const soil_struct *,
    const wstate_struct *, const wflux_struct *, const phystate_struct *);
int             FinalHarvestDate(int, double, double, double);
void            FirstDay(const soiltbl_struct *, const ctrl_struct *,
    elem_struct []);
void            FirstDOY(int, const cstate_struct *, mgmt_struct *);
void            FixedHarvest(int, int, const tillage_struct *,
    const phystate_struct *, crop_struct [], wstate_struct *, cstate_struct *,
    nstate_struct *, nflux_struct *);
void            ForageSeedHarvest(int, int, const phystate_struct *,
    crop_struct *, wstate_struct *, cstate_struct *, nstate_struct *,
    nflux_struct *);
int             ForcedClip(int, crop_struct []);
int             ForcedMaturity(int, int, int, int, int, int);
double          Fraction(double, double, double, double, double);
void            GrainHarvest(int, int, double, crop_struct *, wstate_struct *,
    cstate_struct *, nstate_struct *);
void            GrowingCrop(int, int, int, const soil_struct *,
    const weather_struct *, crop_struct [], wstate_struct *,
    wflux_struct *, cstate_struct *, nstate_struct *, nflux_struct *,
    phystate_struct *);
void            InitAgVar(elem_struct [], river_struct [], N_Vector);
void            InitCropStateVar(crop_struct *);
void            InitCycles(const calib_struct *, const agtbl_struct *,
    const mgmt_struct [], const crop_struct [], const soiltbl_struct *,
    elem_struct []);
double          IntegRoot(double, double);
int             InTimeWindow(int, int, int);
int             IsLeapYear(int);
int             IsOperToday(int, int, const soil_struct *,
    const wstate_struct *, const estate_struct *, int *, mgmt_struct *);
void            KillCrop(int, int, const phystate_struct *, crop_struct *,
    wstate_struct *, cstate_struct *, nstate_struct *);
void            LateralNFlow(double, const soil_struct *,
    const wstate_struct *, const phystate_struct *, const double [], double,
    double, double []);
double          LinearEqmConc(double, double, double, double, double);
double          MaxAbgdHumifFactor(double);
double          MaxManureHumifFactor(double);
double          MaxRhizoHumifFactor(double);
double          MaxRootHumifFactor(double);
double          MobileNConc(double, const double [], const soil_struct *,
    const wstate_struct *, const phystate_struct *);
double          Moisture(double);
double          N2OFracNitrif(double);
void            Nitrification(const soil_struct *, const wstate_struct *,
    const estate_struct *, const phystate_struct *, nstate_struct *,
    nflux_struct *);
double          NMineral(double, double, double, double);
int             NumActiveCrop(const crop_struct []);
void            NXform(const crop_struct [], const weather_struct *,
    const soil_struct *, const wstate_struct *, const estate_struct *,
    const cstate_struct *, const cflux_struct *, const phystate_struct *,
    nstate_struct *, nflux_struct *);
double          PHFunction(double);
void            Phenology(const weather_struct *, crop_struct []);
void            PlantCrop(int, const plant_struct *, crop_struct *);
void            PotSoluteUptake(double, const double [], const soil_struct *,
    const wstate_struct *, const wflux_struct *, const phystate_struct *,
    double *, double []);
void            Processes(int, crop_struct [], const soil_struct *,
    const weather_struct *, wstate_struct *, wflux_struct *, cstate_struct *,
    nstate_struct *, nflux_struct *, phystate_struct *);
double          Profile(int, const double []);
void            RadIntcp(crop_struct []);
void            ReadCrop(const char [], crop_struct []);
void            ReadCyclesCtrl(const char [], agtbl_struct *, ctrl_struct *);
void            ReadMultOper(const agtbl_struct *, mgmt_struct [],
    crop_struct []);
void            ReadOper(const char [], int, int, mgmt_struct *,
    crop_struct []);
void            ReadSoilInit(const char [], soiltbl_struct *);
void            ResidueCover(const cstate_struct *, phystate_struct *);
void            ResetCrop(crop_struct *);
void            ResidueEvap(double, double, const crop_struct [],
    const cstate_struct *, const phystate_struct *, wstate_struct *,
    wflux_struct *);
void            ResidueWetting(double, const cstate_struct *,
    const phystate_struct *, wstate_struct *, wflux_struct *);
void            RootFrac(const crop_struct *, const phystate_struct *,
    double []);
double          SatVP(double);
double          ShootBiomassPartn(int, double, double, double);
double          SoilBufferPower(double, double, double);
void            SoilCarbonBalance(const double [], const soil_struct *,
    const crop_struct [], const phystate_struct *, wstate_struct *,
    cstate_struct *, cflux_struct *, nstate_struct *, nflux_struct *);
void            SoilEvap(double, double, const crop_struct [],
    const soil_struct *, const phystate_struct *, wstate_struct *,
    wflux_struct *);
double          SoilWaterContent(double, double, double, double);
double          SoilWaterPot(double, double, double, double);
void            SoluteConc(double, elem_struct [], river_struct []);
void            SoluteTransp(double, double, const double [], const double [],
    const soil_struct *, const phystate_struct *, double []);
double          TextureFactor(double);
double          ThermalTime(double, double, double, double);
void            TillageFactor(double, const tillage_struct *,
    const soil_struct *, const phystate_struct *, double []);
void            TillageFactorSet(int , const double [], double , double []);
double          TmpFunc(double);
double          TmpFuncGrowth(double, double, double, double);
double          TmpLimit(double, double, double);
void            UpdateNProfile(double, const soil_struct *,
    const wstate_struct *, const nstate_struct *, const solute_struct [],
    double [], double [], phystate_struct *);
void            UpdateOperPtr(mgmt_struct *);
double          VolatilDepthFunc(double);
void            Volatilization(const weather_struct *, const crop_struct [],
    const soil_struct *, const wstate_struct *, const cstate_struct *,
    const phystate_struct *, nstate_struct *, nflux_struct *);
double          VolWCAt33Jkg(double, double, double);
double          VolWCAt1500Jkg(double, double, double);
double          WaterContentLimitToEvap(double, double, double);
void            WaterUptake(double, const soil_struct *,
    const weather_struct *, const phystate_struct *,
    crop_struct [], wstate_struct *, wflux_struct *);
void            ZeroFluxes(wflux_struct *, cflux_struct *, nflux_struct *);
void            ZeroHarvest(crop_struct *);
void            NRT(double, double, double [], const soil_struct *,
    const wstate_struct *, const wstate_struct *, const phystate_struct *,
    double []);
#endif
#if defined(_CYCLES_OBSOLETE_)
void            ApplyFertilizer(const fixfert_struct *, cstate_struct *,
    nstate_struct *, nflux_struct *);
double          AvgSolConc(int, double, const double [],
    const double [], const double [], double, const double []);
#endif

#if defined(_RT_)
void            InitChem(const char [], const calib_struct *, forc_struct *forc,
    chemtbl_struct [], kintbl_struct [], rttbl_struct *, chmictbl_struct *,
    elem_struct []);
void            Reaction(double, const chemtbl_struct [], const kintbl_struct [],
    const rttbl_struct *, elem_struct []);
int             _React(double, const chemtbl_struct [], const kintbl_struct [],
    const rttbl_struct *, double, double, chmstate_struct *);
void            ReactControl(const chemtbl_struct [], const kintbl_struct [],
    const rttbl_struct *, double, double, double, chmstate_struct *, double []);
void            Lookup(FILE *, const calib_struct *, chemtbl_struct [],
    kintbl_struct [], rttbl_struct *);
void            Speciation(const chemtbl_struct [], const rttbl_struct *,
    river_struct []);
int             _Speciation(const chemtbl_struct [], const rttbl_struct *, int,
    chmstate_struct *);
int             SpeciesType(FILE *, const char []);
void            Unwrap(char *, const char *);
double          EqvUnsatH(double, double, double, double, double);
double          UnsatSatRatio(double, double, double);
void            SortChem(char[][MAXSTRING], const int [], int, chemtbl_struct []);
int             FindChem(const char [], const chemtbl_struct [], int);
void            ReadChem(const char[], const char[], chemtbl_struct [],
    kintbl_struct [], rttbl_struct *, forc_struct *, ctrl_struct *);
void            ReadPrep(const char[], const chemtbl_struct [],
    const rttbl_struct *, forc_struct *forc);
void            ReadCini(const char[], const chemtbl_struct *, int,
    atttbl_struct *, chmictbl_struct *);
int             ParseLocation(const char [], const char [], int);
void            ApplyPrcpConc(const  rttbl_struct *, forc_struct *,
    elem_struct [], int);
void            wrap(char *);
void            SoluteConc(const chemtbl_struct [], const rttbl_struct *,
    elem_struct [], river_struct []);
void            RTUpdate(const rttbl_struct *, elem_struct [], river_struct []);
void            InitRTVar(const chemtbl_struct [], const rttbl_struct *,
    elem_struct [], river_struct [], N_Vector);
int             MatchWrappedKey(const char [], const char []);
void            ReadTempPoints(const char [], double, int *, int *);
void            ReadDHParam(const char [], int, double *);
void            ReadPrimary(const char [], int, chemtbl_struct []);
void            ReadSecondary(const char [], int, int, chemtbl_struct [],
    rttbl_struct *);
void            ReadMinerals(const char [], int, int, double [][MAXSPS],
    double [], chemtbl_struct [], rttbl_struct *);
void            ReadAdsorption(const char [], int, int, chemtbl_struct [],
    rttbl_struct *);
void            ReadCationEchg(const char [], double, chemtbl_struct [],
    rttbl_struct *);
void            ReadMinKin(FILE *, int, double, int *, char [], chemtbl_struct [],
    kintbl_struct *);
void            InitChemS(const chemtbl_struct [], const rttbl_struct *,
    const rtic_struct *, double, double, chmstate_struct *);
void            ReadChemAtt(const char *, atttbl_struct *);
void            ReadRtIc(const char *, elem_struct []);
void            UpdatePConc(const rttbl_struct *, elem_struct [],
    river_struct []);
void            WriteRtIc(const char *, const chemtbl_struct [],
    const rttbl_struct *, elem_struct []);
double          SoilTempFactor(double);
#endif

#endif
