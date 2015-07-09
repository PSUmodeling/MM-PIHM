#ifndef NOAH_HEADER
#define NOAH_HEADER

/*
 * Define constants 
 */

#define RD	    287.04
#define SIGMA	    5.67e-8
//#define CP	    1004.6
#define CPH2O	    4.218e3
#define CPICE	    2.106e3
#define LSUBF	    3.335e5
#define EMISSI_S    0.95
#define NLUS	    50
#define NSLTYPE	    32767
#define NSLOPE	    32767
#define TFREEZ	    273.15
//#define LVH2O	    2.501e6
#define LSUBS	    2.83e6
#define R	    287.04

#define MAXLYR     11

enum lsm_forcing_type {SOLAR_DIR_TS, SOLAR_DIF_TS};
enum lsm_print_type {T1_CTRL, STC_CTRL, SMC_CTRL, SH2O_CTRL, SNOWH_CTRL,
    ALBEDO_CTRL, LE_CTRL, SH_CTRL, G_CTRL, ETP_CTRL, ESNOW_CTRL,
    ROOTW_CTRL, SOILM_CTRL};

typedef struct genprmt_type
{
    double          sbeta_data;
    double          fxexp_data;
    double          csoil_data;
    double          salp_data;
    double          refdk_data;
    double          refkdt_data;
    double          frzk_data;
    double          zbot_data;
#ifdef _NOAH_
    double          tbot_data;
#endif
    double          smlow_data;
    double          smhigh_data;
    double          czil_data;
    double          lvcoef_data;
} genprmt_type;

typedef struct grid_type
{
    int             rdlai2d;    /* if rdlai2d == 1, then the xlai value that
                                 * we pass to sflx will be used;
                                 * if rdlai2d == 0, then xlai will be computed
                                 * within sflx, from table minimum and maximum
                                 * values in vegparm.tbl, and the current green
                                 * vegetation fraction. */
    int             usemonalb;  /* if usemonalb == 1, then the alb value passed to
                                 * sflx will be used as the background  snow-free
                                 * albedo term.
                                 * if usemonalb == 0, then alb will be computed
                                 * within sflx from minimum and maximum values in
                                 * vegparm.tbl, and the current green vegetation
                                 * fraction. */
    int             iz0tlnd;    /* option to turn on (iz0tlnd=1) or off (iz0tlnd=0)
                                 * the vegetation-category-dependent calculation of
                                 * the zilitinkivich coefficient czil in the sfcdif
                                 * subroutines. */

    double          rtdis[MAXLYR];

    /*
     * configuration information 
     */
    double          dt;         /* timestep (sec) (dt should not exceed 3600 secs,
                                 * recommend 1800 secs or less) */
    double          zlvl;       /* height (m) above ground of atmospheric forcing
                                 * variables */
    double          zlvl_wind;  /* height (m) above ground of wind observations */
    int             nsoil;      /* number of soil layers (at least 2, and not
                                 * greater than parameter nsold set below) */

    double          sldpth[MAXLYR];     /* the thickness of each soil layer (m) */

    /*
     * logical information 
     */
    int             lch;        /* exchange coefficient (ch) calculation flag
                                 * 0: using ch-routine sfcdif;
                                 * 1: ch is brought in */
    int             local;      /* flag for local-site simulation (where there is
                                 * no maps for albedo, veg fraction,
                                 * and roughness.
                                 * 1:  all lsm parameters (inluding albedo, veg
                                 * fraction and roughness length) will be defined
                                 * by three tables */
    int             llanduse;   /* (=usgs, using usgs landuse classification) */
    int             lsoil;      /* (=stas, using fao/statsgo soil texture
                                 * classification) */
    int             isurban;

    /*
     * forcing data 
     */
    double          longwave;   /* global lw downward radiation */
    double          lwdn;       /* lw downward radiation (w m-2; positive, not net
                                 * longwave) */
    double          soldn;      /* solar downward radiation (w m-2; positive, not net solar) */
    double          solnet;     /* net downward solar radiation ((w m-2; positive) */
    double          sfcprs;     /* pressure at height zlvl above ground (pascals) */
    double          prcp;       /* precip rate (kg m-2 s-1) (note, this is a rate) */
#ifdef _NOAH_
    double          pcpdrp;     /* combined prcp1 and drip (from cmc) that goes into the soil (m s-1) */
#endif
    double          sfctmp;     /* air temperature (k) at height zlvl above ground */
    double          th2;        /* air potential temperature (k) at height zlvl above ground */
    double          q2;         /* mixing ratio at height zlvl above ground (kg kg-1) */
    double          cosz;       /* solar zenith angle (not used for now) */
    double          prcprain;   /* liquid-precipitation rate (kg m-2 s-1) (not used) */
    double          solardirect;    /* direct component of downward solar radiation (w m-2) (not used) */
    double          ffrozp;     /* fraction of frozen precipitation */

    /*
     * other forcing (input) data 
     */
    double          sfcspd;     /* wind speed (m s-1) at height zlvl above ground */
    double          q2sat;      /* sat specific humidity at height zlvl above ground (kg kg-1) */
    double          dqsdt2;     /* slope of sat specific humidity curve at t=sfctmp (kg kg-1 k-1) */
    /*
     * canopy/soil characteristics 
     */
    int             vegtyp;     /* vegetation type (integer index) */
    int             soiltyp;    /* soil type (integer index) */
    int             slopetyp;   /* class of sfc slope (integer index) */

    double          shdfac;     /* areal fractional coverage of green vegetation (fraction= 0.0-1.0) */
    double          shdmin;     /* minimum areal fractional coverage of green vegetation (fraction= 0.0-1.0) <= shdfac */
    double          shdmax;     /* maximum areal fractional coverage of green vegetation (fraction= 0.0-1.0) >= shdfac */
    double          ptu;        /* photo thermal unit (plant phenology for annuals/crops) (not yet used, but passed to redprm for future use in veg parms) */
    double          alb;        /* backround snow-free surface albedo (fraction), for julian day of year (usually from temporal interpolation of monthly */
    /*
     * mean values' calling prog may or may not include diurnal sun angle effect) 
     */
    double          albedomax;  /* maximum background albedo */
    double          albedomin;  /* manimum background albedo */
    double          snoalb;     /* upper bound on maximum albedo over deep snow (e.g. from robinson and kukla, 1985, j. clim. & appl. meteor.) */
    double          tbot;       /* bottom soil temperature (local yearly-mean sfc air temperature) */
    double          z0brd;      /* background fixed roughness length (m) */
    double          z0;         /* time varying roughness length (m) as function of snow depth */
    double          z0max;      /* maximum roughness length (m) */
    double          z0min;      /* minimum roughness length (m) */
    double          embrd;      /* background surface emissivity (between 0 and 1) */
    double          emissi;     /* surface emissivity (between 0 and 1) */
    double          emissmax;   /* maximum emissivity */
    double          emissmin;   /* minimum emmisivity */
    double          laimax;     /* maximum lai */
    double          laimin;     /* minimum lai */

    double          snup;       /* threshold snow depth (in water equivalent m) that implies 100 percent snow cover */
    double          rgl;        /* parameters used in radiation stress function */
    double          hs;         /* parameter used in vapor pressure deficit function */
    double          topt;       /* optimum transpiration air temperature */
    double          cmcmax;     /* maximum canopy water capacity */
    double          cfactr;     /* parameter used in the canopy inteception calculation */
    double          rsmax;      /* max. stomatal resistance */

    double          vgalpha;    /* ys: van genuchten alpha */
    double          vgbeta;     /* ys: van genuchten beta */
    double          bexp;       /* b parameter */
    double          csoil;      /* soil heat capacity (j m-3 k-1) */

    double          refdk;      /* parameter in the surface runoff parameterization */
    double          refkdt;     /* parameter in the surface runoff parameterization */

    double          frzk;       /* frozen ground parameter */
    double          frzx;
    double          zbot;       /* depth (m) of lower boundary soil temperature */
    double          czil;       /* calculate roughness length of heat */

    double          dksat;      /* saturation soil conductivity */
#ifdef _NOAH_
    double          macksat;    /* ys: flux-pihm saturation macropore conductivity */
    double          areaf;      /* ys: flux-pihm fractional area of macropore */
    double          infil;        /* ys: flux-pihm infiltration */
    int             mac_status;
    int             nmacd;      /* ys: flux-pihm the layer where the macropore depth reaches */
    int             nwtbl;      /* ys: flux-pihm the layer where the water table is */
    double          cmcfactr;
#endif
    double          dwsat;      /* saturation soil diffusivity */
    double          f1;         /* soil thermal diffusivity/conductivity coefficient */
    double          fxexp;      /* soil evaporation exponent used in devap */
    double          psisat;     /* saturation soil potential */
    double          quartz;     /* soil quartz content */
    double          kdt;
    double          lvcoef;
    double          sbeta;      /* parameter used to calculate vegetation effect on soil heat */

    double          slope;      /* linear reservoir coefficient */
    double          salp;       /* shape parameter of distribution function of snow cover */
#ifdef _NOAH_
    double          aspect;     /* ys: surface aspect of grid */
    double          svf;        /* ys: sky view factor */
    double          h_phi[36];  /* unobstrcted angle */
#endif

    /*
     * history (state) variables 
     */
    double          cmc;        /* canopy moisture content (m) */
    double          t1;         /* ground/canopy/snowpack) effective skin temperature (k) */
    double          stc[MAXLYR];        /* soil temp (k) */
    double          smc[MAXLYR];        /* total soil moisture content (volumetric fraction) */
    double          sh2o[MAXLYR];       /* unfrozen soil moisture content (volumetric fraction). note: frozen soil moisture = smc - sh2o */
    double          snowh;      /* actual snow depth (m) */
    double          sneqv;      /* liquid water-equivalent snow depth (m). note: snow density = sneqv/snowh */
    double          albedo;     /* surface albedo including snow effect (unitless fraction) =snow-free albedo (alb) when sneqv=0, */
    /*
     * or =fct(msnoalb,alb,vegtyp,shdfac,shdmin) when sneqv>0 
     */
    double          ch;         /* surface exchange coefficient for heat and moisture (m s-1); */
    /*
     * note: ch is technically a conductance since it has been multiplied by wind speed. 
     */
    double          cm;         /* surface exchange coefficient for momentum (m s-1); */
    /*
     * note: cm is technically a conductance since it has been multiplied by wind speed. 
     */


    /*
     * output
     * output variables necessary for a coupled numerical weather prediction model,
     * e.g. noaa/nws/ncep mesoscale eta model.  for this application,
     * the remaining output/diagnostic/parameter blocks below are not
     * necessary.  other applications may require different output variables. 
     */
    double          eta;        /* actual latent heat flux (w m-2: negative, if up from surface) */
    double          eta_kinematic;  /* atctual latent heat flux in kg m-2 s-1 */
    double          sheat;      /* sensible heat flux (w m-2: negative, if upward from surface) */
    double          fdown;      /* radiation forcing at the surface (w m-2) = soldn*(1-alb)+lwdn */

    double          ec;         /* canopy water evaporation (w m-2) */
    double          edir;       /* direct soil evaporation (w m-2) */
    double          et[MAXLYR];         /* plant transpiration from a particular root (soil) layer (w m-2) */
    double          ett;        /* total plant transpiration (w m-2) */
    double          esnow;      /* sublimation from (or deposition to if <0) snowpack (w m-2) */
    double          drip;       /* through-fall of precip and/or dew in excess of canopy water-holding capacity (m) */
    double          dew;        /* dewfall (or frostfall for t<273.15) (m) */

    double          beta;       /* ratio of actual/potential evap (dimensionless) */
    double          etp;        /* potential evaporation (w m-2) */
    double          ssoil;      /* soil heat flux (w m-2: negative if downward from surface) */
    double          flx1;       /* precip-snow sfc (w m-2) */
    double          flx2;       /* freezing rain latent heat flux (w m-2) */
    double          flx3;       /* phase-change heat flux from snowmelt (w m-2) */

    double          snomlt;     /* snow melt (m) (water equivalent) */
    double          sncovr;     /* fractional snow cover (unitless fraction, 0-1) */

    double          runoff1;    /* surface runoff (m s-1), not infiltrating the surface */
    double          runoff2;    /* subsurface runoff (m s-1), drainage out bottom of last soil layer (baseflow) */
    double          runoff3;    /* numerical trunctation in excess of porosity (smcmax) for a given soil layer at the end of a time step (m s-1). note: the above runoff2 is actually the sum of runoff2 and runoff3 */

    double          rc;         /* canopy resistance (s m-1) */
    double          pc;         /* plant coefficient (unitless fraction, 0-1) where pc*etp = actual transp */
    double          xlai;       /* leaf area index (dimensionless) */
    double          rsmin;      /* minimum canopy resistance (s m-1) */
    double          rcs;        /* incoming solar rc factor (dimensionless) */
    double          rct;        /* air temperature rc factor (dimensionless) */
    double          rcq;        /* atmos vapor pressure deficit rc factor (dimensionless) */
    double          rcsoil;     /* soil moisture rc factor (dimensionless) */

    /*
     * diagnostic output 
     */
    double          soilw;      /* available soil moisture in root zone (unitless fraction between smcwlt and smcmax) */
    double          soilm;      /* total soil column moisture content (frozen+unfrozen) (m) */
    double          q1;         /* effective mixing ratio at surface (kg kg-1), used for diagnosing the mixing ratio at 2 meter for coupled model */
    double          smav[MAXLYR];       /* soil moisture availability for each layer, as a fraction between smcwlt and smcmax. */

    double          snotime1;   /* age of the snow on the ground */
    double          ribb;       /* bulk richardson number used to limit the dew/frost */

    /*
     * parameters 
     */
    double          smcwlt;     /* wilting point (volumetric) */
    double          smcdry;     /* dry soil moisture threshold where direct evap frm top layer ends (volumetric) */
    double          smcref;     /* soil moisture threshold where transpiration begins to stress (volumetric) */
    double          smcmax;     /* porosity, i.e. saturated value of soil moisture (volumetric) */
    double          smcmin;     /* ys: flux-pihm residual porosity */
    int             nroot;      /* number of root layers, a function of veg type, determined in subroutine redprm. */

    double          avgsubflux[3];
    double          avgrunoff;
    double          avginfil;

    double         *radn[2];
} grid_struct;

typedef struct radn_ts_struct
{
    int             nts;

    ts_struct      *ts;

    double         *radn[2];
} radn_ts_struct;

typedef struct lsm_ic_struct
{
    double         *t1;
    double         *snowh;
    double        **stc;
    double        **smc;
    double        **sh2o;
} lsm_ic_struct;

typedef struct lsm_struct
{
    genprmt_type    genprmt;
    grid_struct      *grid;
    lsm_ic_struct   ic;
#ifdef _NOAH_
    int             std_nsoil;
    double          std_sldpth[MAXLYR];
    int             rad_mode;   /* radiation mode; 1: topographic, 0: uniform */
    prtctrl_struct  prtctrl[NUM_PRINT];

    radn_ts_struct  forcing;

    int             nprint;
    int             prtvrbl[NUM_PRINT];
#endif
    double          longitude;
    double          latitude;
}              *lsm_struct;

void            SFlx (grid_struct *);

void            AlCalc (double *alb, double *snoalb, double *embrd, double *shdfac, double *shdmin, double *sncovr, double *tsnow, double *albedo, double *emissi, double *dt, int *snowng, double *snotime1, double *lvcoef);

void            CanRes (double *solar, double *ch, double *sfctmp, double *q2, double *sfcprs, double *smc, double *zsoil, int *nsoil, double *smcwlt, double *smcref, double *rsmin, double *rc, double *pc, int *nroot, double *q2sat, double *dqsdt2, double *topt, double *rsmax, double *rgl, double *hs, double *xlai, double *rcs, double *rct, double *rcq, double *rcsoil, double *emissi);

void            CSnow (double *sncond, double *dsnow);

#ifdef _NOAH_
void            DEvap (double *edir, double *etp1, double *smc, double *zsoil, double *shdfac, double *smcmax, double *dksat, double *smcdry, double *smcref, double *smcwlt, double *fxexp);
#else
void            DEvap (double *edir, double *etp1, double *smc, double *zsoil, double *shdfac, double *smcmax, double *bexp, double *dksat, double *dwsat, double *smcdry, double *smcref, double *smcwlt, double *fxexp);
#endif

#ifdef _NOAH_
void            Evapo (double *eta1, double *smc, int *nsoil, double *cmc, double *etp1, double *dt, double *zsoil, double *sh2o, double *smcmax, double *pc, double *smcwlt, double *dksat, double *smcref, double *shdfac, double *cmcmax, double *smcdry, double *cfactr, double *edir, double *ec, double *et, double *ett, double *sfctmp, double *q2, int *nroot, double *rtdis, double *fxexp);
#else
void            Evapo (double *eta1, double *smc, int *nsoil, double *cmc, double *etp1, double *dt, double *zsoil, double *sh2o, double *smcmax, double *bexp, double *pc, double *smcwlt, double *dksat, double *dwsat, double *smcref, double *shdfac, double *cmcmax, double *smcdry, double *cfactr, double *edir, double *ec, double *et, double *ett, double *sfctmp, double *q2, int *nroot, double *rtdis, double *fxexp);
#endif

void            Fac2Mit (double *smcmax, double *flimit);

#ifdef _NOAH_
void            FrH2O (double *free, double *tkelv, double *smc, double *sh2o, double *smcmax, double *smcmin, double *vgalpha, double *vgbeta);
#else
void            FrH2O (double *free, double *tkelv, double *smc, double *sh2o, double *smcmax, double *bexp, double *psis);
#endif

#ifdef _NOAH_
void            HRT (double *rhsts, double *stc, double *smc, double *smcmax, double *smcmin, int *nsoil, double *zsoil, double *yy, double *zz1, double *tbot, double *zbot, double *sh2o, double *dt, double *vgalpha, double *vgbeta, double *f1, double *df1, double *quartz, double *csoil, double *ai, double *bi, double *ci, int *vegtyp, int *isurban);
#else
void            HRT (double *rhsts, double *stc, double *smc, double *smcmax, int *nsoil, double *zsoil, double *yy, double *zz1, double *tbot, double *zbot, double *psisat, double *sh2o, double *dt, double *bexp, double *f1, double *df1, double *quartz, double *csoil, double *ai, double *bi, double *ci, int *vegtyp, int *isurban);
#endif

void            HStep (double *stcout, double *stcin, double *rhsts, double *dt, int *nsoil, double *ai, double *bi, double *ci);

#ifdef _NOAH_
void            NoPac (double *etp, double *eta, double *prcp, double *pcpdrp, double *smc, double *smcmax, double *smcmin, double *smcwlt, double *smcref, double *smcdry, double *cmc, double *cmcmax, int *nsoil, double *dt, double *shdfac, double *sbeta, double *q2, double *t1, double *sfctmp, double *t24, double *th2, double *fdown, double *f1, double *emissi, double *ssoil, double *stc, double *epsca, double *vgalpha, double *vgbeta, double *macksat, double *areaf, int *nmacd, int *mac_status, int *nwtbl, double *pc, double *rch, double *rr, double *cfactr, double *sh2o, double *frzfact, double *zsoil, double *dksat, double *tbot, double *zbot, double *inf, double *runoff2, double *runoff3, double *edir, double *ec, double *et, double *ett, int *nroot, double *rtdis, double *quartz, double *fxexp, double *csoil, double *beta, double *drip, double *dew, double *flx1, double *flx3, int *vegtyp, int *isurban);
#else
void            NoPac (double *etp, double *eta, double *prcp, double *smc, double *smcmax, double *smcwlt, double *smcref, double *smcdry, double *cmc, double *cmcmax, int *nsoil, double *dt, double *shdfac, double *sbeta, double *q2, double *t1, double *sfctmp, double *t24, double *th2, double *fdown, double *f1, double *emissi, double *ssoil, double *stc, double *epsca, double *bexp, double *pc, double *rch, double *rr, double *cfactr, double *sh2o, double *slope, double *kdt, double *frzfact, double *psisat, double *zsoil, double *dksat, double *dwsat, double *tbot, double *zbot, double *runoff1, double *runoff2, double *runoff3, double *edir, double *ec, double *et, double *ett, int *nroot, double *rtdis, double *quartz, double *fxexp, double *csoil, double *beta, double *drip, double *dew, double *flx1, double *flx3, int *vegtyp, int *isurban);
#endif

void            Penman (double *sfctmp, double *sfcprs, double *ch, double *t2v, double *th2, double *prcp, double *fdown, double *t24, double *ssoil, double *q2, double *q2sat, double *etp, double *rch, double *epsca, double *rr, int *snowng, int *frzgra, double *dqsdt2, double *flx2, double *emissi_in, double *sneqv, double *t1, double *sncovr);

void            RedPrm (grid_struct * noah, lsm_struct lsm, double *zsoil);

void            Rosr12 (double *p, double *a, double *b, double *c, double *d, double *delta, int *nsoil);

#ifdef _NOAH_
void            ShFlx (double *ssoil, double *stc, double *smc, double *smcmax, double *smcmin, int *nsoil, double *t1, double *dt, double *yy, double *zz1, double *zsoil, double *tbot, double *zbot, double *smcwlt, double *sh2o, double *vgalpha, double *vgbeta, double *f1, double *df1, double *quartz, double *csoil, int *vegtyp, int *isurban);
#else
void            ShFlx (double *ssoil, double *stc, double *smc, double *smcmax, int *nsoil, double *t1, double *dt, double *yy, double *zz1, double *zsoil, double *tbot, double *zbot, double *smcwlt, double *psisat, double *sh2o, double *bexp, double *f1, double *df1, double *quartz, double *csoil, int *vegtyp, int *isurban);
#endif

#ifdef _NOAH_
void            SmFlx (double *smc, int *nsoil, double *cmc, double *dt, double *prcp1, double *pcpdrp, double *zsoil, double *sh2o, double *frzfact, double *smcmax, double *smcmin, double *vgalpha, double *vgbeta, double *macksat, double *areaf, int *nmacd, int *mac_status, int *nwtbl, double *smcwlt, double *dksat, double *shdfac, double *cmcmax, double *inf, double *runoff2, double *runoff3, double *edir, double *ec, double *et, double *drip);
#else
void            SmFlx (double *smc, int *nsoil, double *cmc, double *dt, double *prcp1, double *zsoil, double *sh2o, double *slope, double *kdt, double *frzfact, double *smcmax, double *bexp, double *smcwlt, double *dksat, double *dwsat, double *shdfac, double *cmcmax, double *runoff1, double *runoff2, double *runoff3, double *edir, double *ec, double *et, double *drip);
#endif

void            SnFrac (double *sneqv, double *snup, double *salp, double *snowh, double *sncovr);

#ifdef _NOAH_
void            SnkSrc (double *tsnsr, double *tavg, double *smc, double *sh2o, double *zsoil, int *nsoil, double *smcmax, double *smcmin, double *vgalpha, double *vgbeta, double *dt, int k, double *qtot);
#else
void            SnkSrc (double *tsnsr, double *tavg, double *smc, double *sh2o, double *zsoil, int *nsoil, double *smcmax, double *psisat, double *bexp, double *dt, int k, double *qtot);
#endif

#ifdef _NOAH_
void            SnoPac (double *etp, double *eta, double *prcp, double *prcpf, double *pcpdrp, int *snowng, double *smc, double *smcmax, double *smcmin, double *smcwlt, double *smcref, double *smcdry, double *cmc, double *cmcmax, int *nsoil, double *dt, double *sbeta, double *df1, double *q2, double *t1, double *sfctmp, double *t24, double *th2, double *fdown, double *f1, double *ssoil, double *stc, double *epsca, double *sfcprs, double *vgalpha, double *vgbeta, double *macksat, double *areaf, int *nmacd, int *mac_status, int *nwtbl, double *pc, double *rch, double *rr, double *cfactr, double *sncovr, double *esd, double *sndens, double *snowh, double *sh2o, double *frzfact, double *zsoil, double *dksat, double *tbot, double *zbot, double *shdfac, double *inf, double *runoff2, double *runoff3, double *edir, double *ec, double *et, double *ett, int *nroot, double *snomlt, double *rtdis, double *quartz, double *fxexp, double *csoil, double *beta, double *drip, double *dew, double *flx1, double *flx2, double *flx3, double *esnow, double *etns, double *emissi, double *ribb, double *soldn, int *isurban, int *vegtyp);
#else
void            SnoPac (double *etp, double *eta, double *prcp, double *prcpf, int *snowng, double *smc, double *smcmax, double *smcwlt, double *smcref, double *smcdry, double *cmc, double *cmcmax, int *nsoil, double *dt, double *sbeta, double *df1, double *q2, double *t1, double *sfctmp, double *t24, double *th2, double *fdown, double *f1, double *ssoil, double *stc, double *epsca, double *sfcprs, double *bexp, double *pc, double *rch, double *rr, double *cfactr, double *sncovr, double *esd, double *sndens, double *snowh, double *sh2o, double *slope, double *kdt, double *frzfact, double *psisat, double *zsoil, double *dwsat, double *dksat, double *tbot, double *zbot, double *shdfac, double *runoff1, double *runoff2, double *runoff3, double *edir, double *ec, double *et, double *ett, int *nroot, double *snomlt, double *rtdis, double *quartz, double *fxexp, double *csoil, double *beta, double *drip, double *dew, double *flx1, double *flx2, double *flx3, double *esnow, double *etns, double *emissi, double *ribb, double *soldn, int *isurban, int *vegtyp);
#endif

void            SnowPack (double *esd, double *dtsec, double *snowh, double *sndens, double *tsnow, double *tsoil);

void            Snowz0 (double *sncovr, double *z0, double *z0brd, double *snowh);

void            SnowNew (double *temp, double *newsn, double *snowh, double *sndens);

#ifdef _NOAH_
void            SRT (double *rhstt, double *edir, double *et, double *sh2o, double *sh2oa, int *nsoil, int *nwtbl, double *pcpdrp, double *zsoil, double *dksat, double *smcmax, double *smcmin, double *vgalpha, double *vgbeta, double *macksat, double *areaf, int *nmacd, int *mac_status, double *inf, double *runoff2, double *dt, double *smcwlt, double *frzx, double *sice, double *ai, double *bi, double *ci);
#else
void            SRT (double *rhstt, double *edir, double *et, double *sh2o, double *sh2oa, int *nsoil, double *pcpdrp, double *zsoil, double *dwsat, double *dksat, double *smcmax, double *bexp, double *runoff1, double *runoff2, double *dt, double *smcwlt, double *slope, double *kdt, double *frzx, double *sice, double *ai, double *bi, double *ci);
#endif

#ifdef _NOAH_
void            SStep (double *sh2oout, double *sh2oin, double *cmc, double *rhstt, double *rhsct, double *dt, int *nsoil, double *smcmax, double *smcmin, double *cmcmax, double *runoff3, double *zsoil, double *smc, double *sice, double *ai, double *bi, double *ci);
#else
void            SStep (double *sh2oout, double *sh2oin, double *cmc, double *rhstt, double *rhsct, double *dt, int *nsoil, double *smcmax, double *cmcmax, double *runoff3, double *zsoil, double *smc, double *sice, double *ai, double *bi, double *ci);
#endif

void            TBnd (double *tu, double *tb, double *zsoil, double *zbot, int k, int *nsoil, double *tbnd1);

#ifdef _NOAH_
void            TDfCnd (double *df, double *smc, double *qz, double *smcmax, double *smcmin, double *sh2o);
#else
void            TDfCnd (double *df, double *smc, double *qz, double *smcmax, double *sh2o);
#endif

void            TmpAvg (double *tavg, double *tup, double *tm, double *tdn, double *zsoil, int *nsoil, int k);

void            Transp (double *et, int *nsoil, double *etp1, double *smc, double *cmc, double *zsoil, double *shdfac, double *smcwlt, double *cmcmax, double *pc, double *cfactr, double *smcref, double *sfctmp, double *q2, int *nroot, double *rtdis);

#ifdef _NOAH_
void            WDfCnd (double *wdf, double *wcnd, double *smc, double *smcmax, double *smcmin, double *vgalpha, double *vgbeta, double *dksat, double *macksat, double *areaf, int *mac_status, double *sicemax, double *dsmdz, int *macpore);
#else
void            WDfCnd (double *wdf, double *wcnd, double *smc, double *smcmax, double *bexp, double *dksat, double *dwsat, double *sicemax);
#endif

void            SfcDifOff (double *zlm, double *zlm_wind, double *z0, double *thz0, double *thlm, double *sfcspd, double *czil, double *akms, double *akhs, int *vegtyp, int *isurban, int *iz0tlnd);

#ifdef _NOAH_
double          EFFKV (double ksatfunc, double elemsatn, int status, double mackv, double kv, double areaf);
#endif

double          Pslmu (double zz);
double          Pslms (double zz);
double          Pslhu (double zz);
double          Pslhs (double zz);
double          Pspmu (double xx);
double          Pspms (double yy);
double          Psphu (double xx);
double          Psphs (double yy);

void LsmRead (char *project, lsm_struct noah, pihm_struct pihm);
void LsmInitialize (char *simulation, pihm_struct pihm, lsm_struct noah);
void MapLsmOutput (char *simulation, lsm_struct noah, int numele, char *outputdir);

void PIHMxNoah (int t, double stepsize, pihm_struct pihm, lsm_struct noah);

void AvgFlux (lsm_struct noah, pihm_struct pihm);

void DefSldpth (double *sldpth, int *nsoil, double total_depth, double *std_sldpth, int std_nsoil);
void CalcSlopeAspect (grid_struct *grid, elem_struct elem, pihm_struct pihm);
void LsmSaturationIC (lsm_ic_struct *ic, const grid_struct *grid, const elem_struct *elem, int numele);
void ReadLsmInit (char *project, char *simulation, lsm_ic_struct *ic, int numele);
void InitLsmVrbl (grid_struct *grid, elem_struct *elem, int numele, lsm_ic_struct ic);
int FindLayer (const double *sldpth, int nsoil, double depth);
double          mod (double, double);

double TopoRadiation (double sdir, double sdif, double zenith, double azimuth180, double slope, double aspect, double *h_phi, double svf);
void ApplyRadnForcing (radn_ts_struct *forcing, int t);
#endif
