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
void            ReadIbc (char *, forc_struct *);
void            ReadPara (char *, ctrl_struct *);
void            ReadCalib (char *, calib_struct *);
void            ReadInit (char *, elem_struct *, int, river_struct *, int);
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
double          EffKH (int, double, double, double, double, double, double);

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
    es_struct *, ef_struct *, ps_struct *, lc_struct *, soil_struct *, int);
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
    const lc_struct *, const soil_struct *, const double *, double);
void            Transp (const ws_struct *, wf_struct *, const ps_struct *,
    const lc_struct *, const soil_struct *, const double *);
void            NoPac (ws_struct *, wf_struct *, const wf_struct *,
    es_struct *, ef_struct *, ps_struct *, lc_struct *,
    soil_struct *, const double *, double, double);
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
    const double *, double, double);
void            HRT (ws_struct *, es_struct *, ef_struct *, ps_struct *,
    const lc_struct *, const soil_struct *, double *,
    const double *, double, double, double, double, double *,
    double *, double *);
void            SRT (ws_struct *, wf_struct *, const wf_struct *, ps_struct *,
    const soil_struct *, double *, double *, double *,
    double *, double *, const double *, double);
void            SStep (ws_struct *, wf_struct *, ps_struct *,
    const soil_struct *, double *, double, const double *,
    double *, double *, double *, double *, double);
void            WDfCnd (double *, double *, double, double, double, int,
    const soil_struct *, const ps_struct *);
void            SnoPac (ws_struct *, wf_struct *, const wf_struct *,
    es_struct *, ef_struct *, ps_struct *, lc_struct *,
    const soil_struct *, int, const double *, double, double, double, double);
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
#endif
#endif
