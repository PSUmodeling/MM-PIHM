#ifndef NOAH_HEADER
#define NOAH_HEADER

#include "../pihm.h"
#include "../print.h"
#include "../forcing.h"
/*
#ifndef _FLUX_PIHM_
#define _FLUX_PIHM_
#endif
*/
/*
 * Define constants 
 */

#define RD		    287.04
#define SIGMA	    5.67e-8
#define CP		    1004.6
#define CPH2O	    4.218e3
#define CPICE	    2.106e3
#define LSUBF		3.335e5
#define EMISSI_S	0.95
#define NLUS		50
#define NSLTYPE		30
#define NSLOPE		30
#define TFREEZ		273.15
#define LVH2O		2.501e6
#define LSUBS		2.83e6
#define R		    287.04

enum lsm_forcing_type {SOLAR_DIR_TS, SOLAR_DIF_TS};
/*
 * VEGETATION PARAMETERS
 */

typedef struct VEGTBL_TYPE
{
    int             LUCATS;
    int             BARE;
    int             NATURAL;
    int             ISURBAN;
    char            LUTYPE[256];

    int             NROTBL[NLUS];
    double          SNUPTBL[NLUS];
    double          RSTBL[NLUS];
    double          RGLTBL[NLUS];
    double          HSTBL[NLUS];
    double          SHDTBL[NLUS];
    double          MAXALB[NLUS];
    double          EMISSMINTBL[NLUS];
    double          EMISSMAXTBL[NLUS];
    double          LAIMINTBL[NLUS];
    double          LAIMAXTBL[NLUS];
    double          Z0MINTBL[NLUS];
    double          Z0MAXTBL[NLUS];
    double          ALBEDOMINTBL[NLUS];
    double          ALBEDOMAXTBL[NLUS];
#ifdef _FLUX_PIHM_
    double          CMCFACTRTBL[NLUS];
#endif
    double          CMCMAX_DATA;
    double          TOPT_DATA;
    double          CFACTR_DATA;
    double          RSMAX_DATA;
} VEGTBL_TYPE;

/*
 * SOIL PARAMETERS 
 */
typedef struct SOILTBL_TYPE
{
    int             SLCATS;
    char            SLTYPE[256];
    double          BB[NSLTYPE];
#ifdef _FLUX_PIHM_
    double          VGA[NSLTYPE];   /* YS */
    double          VGB[NSLTYPE];   /* YS */
#endif
    double          DRYSMC[NSLTYPE];
    double          F11[NSLTYPE];
    double          MAXSMC[NSLTYPE];
#ifdef _FLUX_PIHM_
    double          MINSMC[NSLTYPE];    /* YS */
#endif
    double          REFSMC[NSLTYPE];
    double          SATPSI[NSLTYPE];
    double          SATDK[NSLTYPE];
    double          SATDW[NSLTYPE];
    double          WLTSMC[NSLTYPE];
    double          QTZ[NSLTYPE];
    double          MACKSAT[NSLTYPE];
    double          AREAF[NSLTYPE];
#ifdef _FLUX_PIHM_
    int             NMACD[NSLTYPE]; /* YS */
#endif
} SOILTBL_TYPE;


typedef struct GENPRMT_TYPE
{
    int             SLPCATS;
    double          SLOPE_DATA[NSLOPE];
    double          SBETA_DATA;
    double          FXEXP_DATA;
    double          CSOIL_DATA;
    double          SALP_DATA;
    double          REFDK_DATA;
    double          REFKDT_DATA;
    double          FRZK_DATA;
    double          ZBOT_DATA;
#ifdef _FLUX_PIHM_
    double          TBOT_DATA;
#endif
    double          SMLOW_DATA;
    double          SMHIGH_DATA;
    double          CZIL_DATA;
    double          LVCOEF_DATA;
} GENPRMT_TYPE;

typedef struct GRID_TYPE
{
    int             RDLAI2D;    /* If RDLAI2D == 1, then the XLAI value that
                                 * we pass to SFLX will be used;
                                 * If RDLAI2d == 0, then XLAI will be computed
                                 * within SFLX, from table minimum and maximum
                                 * values in VEGPARM.TBL, and the current Green
                                 * Vegetation Fraction. */
    int             USEMONALB;  /* If USEMONALB == 1, then the ALB value passed to
                                 * SFLX will be used as the background  snow-free
                                 * albedo term.
                                 * If USEMONALB == 0, then ALB will be computed
                                 * within SFLX from minimum and maximum values in
                                 * VEGPARM.TBL, and the current Green Vegetation
                                 * Fraction. */
    int             IZ0TLND;    /* Option to turn on (IZ0TLND=1) or off (IZ0TLND=0)
                                 * the vegetation-category-dependent calculation of
                                 * the Zilitinkivich coefficient CZIL in the SFCDIF
                                 * subroutines. */

    double         *RTDIS;

    /*
     * CONFIGURATION INFORMATION 
     */
    double          DT;         /* TIMESTEP (SEC) (DT SHOULD NOT EXCEED 3600 SECS,
                                 * RECOMMEND 1800 SECS OR LESS) */
    double          ZLVL;       /* HEIGHT (M) ABOVE GROUND OF ATMOSPHERIC FORCING
                                 * VARIABLES */
    double          ZLVL_WIND;  /* HEIGHT (M) ABOVE GROUND OF WIND OBSERVATIONS */
    int             NSOIL;      /* NUMBER OF SOIL LAYERS (AT LEAST 2, AND NOT
                                 * GREATER THAN PARAMETER NSOLD SET BELOW) */

    double         *SLDPTH;     /* THE THICKNESS OF EACH SOIL LAYER (M) */

    /*
     * LOGICAL INFORMATION 
     */
    int             LCH;        /* Exchange coefficient (Ch) calculation flag
                                 * 0: using ch-routine SFCDIF;
                                 * 1: Ch is brought in */
    int             LOCAL;      /* Flag for local-site simulation (where there is
                                 * no maps for albedo, veg fraction,
                                 * and roughness.
                                 * 1:  all LSM parameters (inluding albedo, veg
                                 * fraction and roughness length) will be defined
                                 * by three tables */
    int             LLANDUSE;   /* (=USGS, using USGS landuse classification) */
    int             LSOIL;      /* (=STAS, using FAO/STATSGO soil texture
                                 * classification) */
    int             ISURBAN;

    /*
     * FORCING DATA 
     */
    double          LONGWAVE;   /* GLOBAL LW DOWNWARD RADIATION */
    double          LWDN;       /* LW DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET
                                 * LONGWAVE) */
    double          SOLDN;      /* SOLAR DOWNWARD RADIATION (W M-2; POSITIVE, NOT NET SOLAR) */
    double          SOLNET;     /* NET DOWNWARD SOLAR RADIATION ((W M-2; POSITIVE) */
    double          SFCPRS;     /* PRESSURE AT HEIGHT ZLVL ABOVE GROUND (PASCALS) */
    double          PRCP;       /* PRECIP RATE (KG M-2 S-1) (NOTE, THIS IS A RATE) */
#ifdef _FLUX_PIHM_
    double          PCPDRP;     /* COMBINED PRCP1 AND DRIP (FROM CMC) THAT GOES INTO THE SOIL (M S-1) */
#endif
    double          SFCTMP;     /* AIR TEMPERATURE (K) AT HEIGHT ZLVL ABOVE GROUND */
    double          TH2;        /* AIR POTENTIAL TEMPERATURE (K) AT HEIGHT ZLVL ABOVE GROUND */
    double          Q2;         /* MIXING RATIO AT HEIGHT ZLVL ABOVE GROUND (KG KG-1) */
    double          COSZ;       /* Solar zenith angle (not used for now) */
    double          PRCPRAIN;   /* Liquid-precipitation rate (KG M-2 S-1) (not used) */
    double          SOLARDIRECT;    /* Direct component of downward solar radiation (W M-2) (not used) */
    double          FFROZP;     /* FRACTION OF FROZEN PRECIPITATION */

    /*
     * OTHER FORCING (INPUT) DATA 
     */
    double          SFCSPD;     /* WIND SPEED (M S-1) AT HEIGHT ZLVL ABOVE GROUND */
    double          Q2SAT;      /* SAT SPECIFIC HUMIDITY AT HEIGHT ZLVL ABOVE GROUND (KG KG-1) */
    double          DQSDT2;     /* SLOPE OF SAT SPECIFIC HUMIDITY CURVE AT T=SFCTMP (KG KG-1 K-1) */
    /*
     * CANOPY/SOIL CHARACTERISTICS 
     */
    int             VEGTYP;     /* VEGETATION TYPE (INTEGER INDEX) */
    int             SOILTYP;    /* SOIL TYPE (INTEGER INDEX) */
    int             SLOPETYP;   /* CLASS OF SFC SLOPE (INTEGER INDEX) */

    double          SHDFAC;     /* AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION (FRACTION= 0.0-1.0) */
    double          SHDMIN;     /* MINIMUM AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION (FRACTION= 0.0-1.0) <= SHDFAC */
    double          SHDMAX;     /* MAXIMUM AREAL FRACTIONAL COVERAGE OF GREEN VEGETATION (FRACTION= 0.0-1.0) >= SHDFAC */
    double          PTU;        /* PHOTO THERMAL UNIT (PLANT PHENOLOGY FOR ANNUALS/CROPS) (NOT YET USED, BUT PASSED TO REDPRM FOR FUTURE USE IN VEG PARMS) */
    double          ALB;        /* BACKROUND SNOW-FREE SURFACE ALBEDO (FRACTION), FOR JULIAN DAY OF YEAR (USUALLY FROM TEMPORAL INTERPOLATION OF MONTHLY */
    /*
     * MEAN VALUES' CALLING PROG MAY OR MAY NOT INCLUDE DIURNAL SUN ANGLE EFFECT) 
     */
    double          ALBEDOMAX;  /* Maximum background albedo */
    double          ALBEDOMIN;  /* Manimum background albedo */
    double          SNOALB;     /* UPPER BOUND ON MAXIMUM ALBEDO OVER DEEP SNOW (E.G. FROM ROBINSON AND KUKLA, 1985, J. CLIM. & APPL. METEOR.) */
    double          TBOT;       /* BOTTOM SOIL TEMPERATURE (LOCAL YEARLY-MEAN SFC AIR TEMPERATURE) */
    double          Z0BRD;      /* Background fixed roughness length (M) */
    double          Z0;         /* Time varying roughness length (M) as function of snow depth */
    double          Z0MAX;      /* Maximum roughness length (m) */
    double          Z0MIN;      /* Minimum roughness length (m) */
    double          EMBRD;      /* Background surface emissivity (between 0 and 1) */
    double          EMISSI;     /* Surface emissivity (between 0 and 1) */
    double          EMISSMAX;   /* Maximum emissivity */
    double          EMISSMIN;   /* Minimum emmisivity */
    double          LAIMAX;     /* Maximum LAI */
    double          LAIMIN;     /* Minimum LAI */

    double          SNUP;       /* Threshold snow depth (in water equivalent m) that implies 100 percent snow cover */
    double          RGL;        /* Parameters used in radiation stress function */
    double          HS;         /* Parameter used in vapor pressure deficit function */
    double          TOPT;       /* Optimum transpiration air temperature */
    double          CMCMAX;     /* Maximum canopy water capacity */
    double          CFACTR;     /* Parameter used in the canopy inteception calculation */
    double          RSMAX;      /* Max. stomatal resistance */

    double          VGALPHA;    /* YS: Van Genuchten alpha */
    double          VGBETA;     /* YS: Van Genuchten beta */
    double          BEXP;       /* B parameter */
    double          CSOIL;      /* Soil heat capacity (J M-3 K-1) */

    double          REFDK;      /* Parameter in the surface runoff parameterization */
    double          REFKDT;     /* Parameter in the surface runoff parameterization */

    double          FRZK;       /* Frozen ground parameter */
    double          FRZX;
    double          ZBOT;       /* Depth (M) of lower boundary soil temperature */
    double          CZIL;       /* Calculate roughness length of heat */

    double          DKSAT;      /* Saturation soil conductivity */
#ifdef _FLUX_PIHM_
    double          MACKSAT;    /* YS: Flux-PIHM Saturation macropore conductivity */
    double          AREAF;      /* YS: Flux-PIHM fractional area of macropore */
    double          INF;        /* YS: Flux-PIHM infiltration */
    int             NMACD;      /* YS: Flux-PIHM the layer where the macropore depth reaches */
    int             NWTBL;      /* YS: Flux-PIHM the layer where the water table is */
#endif
    double          DWSAT;      /* Saturation soil diffusivity */
    double          F1;         /* Soil thermal diffusivity/conductivity coefficient */
    double          FXEXP;      /* Soil evaporation exponent used in DEVAP */
    double          PSISAT;     /* Saturation soil potential */
    double          QUARTZ;     /* Soil quartz content */
    double          KDT;
    double          LVCOEF;
    double          SBETA;      /* Parameter used to calculate vegetation effect on soil heat */

    double          SLOPE;      /* Linear reservoir coefficient */
    double          SALP;       /* Shape parameter of distribution function of snow cover */
#ifdef _FLUX_PIHM_
    double          ASPECT;     /* YS: Surface aspect of grid */
    double          SVF;        /* YS: Sky view factor */
    double          H_PHI[36];  /* Unobstrcted angle */
#endif

    /*
     * HISTORY (STATE) VARIABLES 
     */
    double          CMC;        /* CANOPY MOISTURE CONTENT (M) */
    double          T1;         /* GROUND/CANOPY/SNOWPACK) EFFECTIVE SKIN TEMPERATURE (K) */
    double         *STC;        /* SOIL TEMP (K) */
    double         *SMC;        /* TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION) */
    double         *SH2O;       /* UNFROZEN SOIL MOISTURE CONTENT (VOLUMETRIC FRACTION). NOTE: FROZEN SOIL MOISTURE = SMC - SH2O */
    double          SNOWH;      /* ACTUAL SNOW DEPTH (M) */
    double          SNEQV;      /* LIQUID WATER-EQUIVALENT SNOW DEPTH (M). NOTE: SNOW DENSITY = SNEQV/SNOWH */
    double          ALBEDO;     /* SURFACE ALBEDO INCLUDING SNOW EFFECT (UNITLESS FRACTION) =SNOW-FREE ALBEDO (ALB) WHEN SNEQV=0, */
    /*
     * OR =FCT(MSNOALB,ALB,VEGTYP,SHDFAC,SHDMIN) WHEN SNEQV>0 
     */
    double          CH;         /* SURFACE EXCHANGE COEFFICIENT FOR HEAT AND MOISTURE (M S-1); */
    /*
     * NOTE: CH IS TECHNICALLY A CONDUCTANCE SINCE IT HAS BEEN MULTIPLIED BY WIND SPEED. 
     */
    double          CM;         /* SURFACE EXCHANGE COEFFICIENT FOR MOMENTUM (M S-1); */
    /*
     * NOTE: CM IS TECHNICALLY A CONDUCTANCE SINCE IT HAS BEEN MULTIPLIED BY WIND SPEED. 
     */


    /*
     * OUTPUT
     * OUTPUT VARIABLES NECESSARY FOR A COUPLED NUMERICAL WEATHER PREDICTION MODEL,
     * E.G. NOAA/NWS/NCEP MESOSCALE ETA MODEL.  FOR THIS APPLICATION,
     * THE REMAINING OUTPUT/DIAGNOSTIC/PARAMETER BLOCKS BELOW ARE NOT
     * NECESSARY.  OTHER APPLICATIONS MAY REQUIRE DIFFERENT OUTPUT VARIABLES. 
     */
    double          ETA;        /* ACTUAL LATENT HEAT FLUX (W m-2: NEGATIVE, IF UP FROM SURFACE) */
    double          ETA_KINEMATIC;  /* atctual latent heat flux in Kg m-2 s-1 */
    double          SHEAT;      /* SENSIBLE HEAT FLUX (W M-2: NEGATIVE, IF UPWARD FROM SURFACE) */
    double          FDOWN;      /* Radiation forcing at the surface (W m-2) = SOLDN*(1-alb)+LWDN */

    double          EC;         /* CANOPY WATER EVAPORATION (W m-2) */
    double          EDIR;       /* DIRECT SOIL EVAPORATION (W m-2) */
    double         *ET;         /* PLANT TRANSPIRATION FROM A PARTICULAR ROOT (SOIL) LAYER (W m-2) */
    double          ETT;        /* TOTAL PLANT TRANSPIRATION (W m-2) */
    double          ESNOW;      /* SUBLIMATION FROM (OR DEPOSITION TO IF <0) SNOWPACK (W m-2) */
    double          DRIP;       /* THROUGH-FALL OF PRECIP AND/OR DEW IN EXCESS OF CANOPY WATER-HOLDING CAPACITY (M) */
    double          DEW;        /* DEWFALL (OR FROSTFALL FOR T<273.15) (M) */

    double          BETA;       /* RATIO OF ACTUAL/POTENTIAL EVAP (DIMENSIONLESS) */
    double          ETP;        /* POTENTIAL EVAPORATION (W m-2) */
    double          SSOIL;      /* SOIL HEAT FLUX (W M-2: NEGATIVE IF DOWNWARD FROM SURFACE) */
    double          FLX1;       /* PRECIP-SNOW SFC (W M-2) */
    double          FLX2;       /* FREEZING RAIN LATENT HEAT FLUX (W M-2) */
    double          FLX3;       /* PHASE-CHANGE HEAT FLUX FROM SNOWMELT (W M-2) */

    double          SNOMLT;     /* SNOW MELT (M) (WATER EQUIVALENT) */
    double          SNCOVR;     /* FRACTIONAL SNOW COVER (UNITLESS FRACTION, 0-1) */

    double          RUNOFF1;    /* SURFACE RUNOFF (M S-1), NOT INFILTRATING THE SURFACE */
    double          RUNOFF2;    /* SUBSURFACE RUNOFF (M S-1), DRAINAGE OUT BOTTOM OF LAST SOIL LAYER (BASEFLOW) */
    double          RUNOFF3;    /* NUMERICAL TRUNCTATION IN EXCESS OF POROSITY (SMCMAX) FOR A GIVEN SOIL LAYER AT THE END OF A TIME STEP (M S-1). Note: the above RUNOFF2 is actually the sum of RUNOFF2 and RUNOFF3 */

    double          RC;         /* CANOPY RESISTANCE (S M-1) */
    double          PC;         /* PLANT COEFFICIENT (UNITLESS FRACTION, 0-1) WHERE PC*ETP = ACTUAL TRANSP */
    double          XLAI;       /* LEAF AREA INDEX (DIMENSIONLESS) */
    double          RSMIN;      /* MINIMUM CANOPY RESISTANCE (S M-1) */
    double          RCS;        /* INCOMING SOLAR RC FACTOR (DIMENSIONLESS) */
    double          RCT;        /* AIR TEMPERATURE RC FACTOR (DIMENSIONLESS) */
    double          RCQ;        /* ATMOS VAPOR PRESSURE DEFICIT RC FACTOR (DIMENSIONLESS) */
    double          RCSOIL;     /* SOIL MOISTURE RC FACTOR (DIMENSIONLESS) */

    /*
     * DIAGNOSTIC OUTPUT 
     */
    double          SOILW;      /* AVAILABLE SOIL MOISTURE IN ROOT ZONE (UNITLESS FRACTION BETWEEN SMCWLT AND SMCMAX) */
    double          SOILM;      /* TOTAL SOIL COLUMN MOISTURE CONTENT (FROZEN+UNFROZEN) (M) */
    double          Q1;         /* Effective mixing ratio at surface (kg kg-1), used for diagnosing the mixing ratio at 2 meter for coupled model */
    double         *SMAV;       /* Soil Moisture Availability for each layer, as a fraction between SMCWLT and SMCMAX. */

    double          SNOTIME1;   /* Age of the snow on the ground */
    double          RIBB;       /* Bulk Richardson number used to limit the dew/frost */

    /*
     * PARAMETERS 
     */
    double          SMCWLT;     /* WILTING POINT (VOLUMETRIC) */
    double          SMCDRY;     /* DRY SOIL MOISTURE THRESHOLD WHERE DIRECT EVAP FRM TOP LAYER ENDS (VOLUMETRIC) */
    double          SMCREF;     /* SOIL MOISTURE THRESHOLD WHERE TRANSPIRATION BEGINS TO STRESS (VOLUMETRIC) */
    double          SMCMAX;     /* POROSITY, I.E. SATURATED VALUE OF SOIL MOISTURE (VOLUMETRIC) */
    double          SMCMIN;     /* YS: Flux-PIHM RESIDUAL POROSITY */
    int             NROOT;      /* NUMBER OF ROOT LAYERS, A FUNCTION OF VEG TYPE, DETERMINED IN SUBROUTINE REDPRM. */
} GRID_TYPE;

typedef struct LSM_Print_Ctrl_struct
{
    char            name[100];
    int             Interval;
    int             NumVar;
    double        **PrintVar;
    double         *buffer;
} LSM_Print_ctrl;

typedef struct LSM_STRUCT
{
    SOILTBL_TYPE    SOILTBL;
    VEGTBL_TYPE     VEGTBL;
    GENPRMT_TYPE    GENPRMT;
    GRID_TYPE      *GRID;
#ifdef _FLUX_PIHM_
    int             STD_NSOIL;
    double         *STD_SLDPTH;
    int             RAD_MODE;   /* Radiation mode; 1: topographic, 0: uniform */
    Print_Ctrl      PCtrl[100];

    TSD            *TSD_rad;

    int             NPRINT;
    int             PRINT_T1;
    int             PRINT_STC;
    int             PRINT_SMC;
    int             PRINT_SH2O;
    int             PRINT_SNOWH;
    int             PRINT_ALBEDO;
    int             PRINT_LE;
    int             PRINT_SH;
    int             PRINT_G;
    int             PRINT_ETP;
#endif
    double          LONGITUDE;
    double          LATITUDE;
}              *LSM_STRUCT;

void            SFLX (GRID_TYPE *);

void            ALCALC (double *ALB, double *SNOALB, double *EMBRD,
   double *SHDFAC, double *SHDMIN, double *SNCOVR, double *TSNOW,
   double *ALBEDO, double *EMISSI, double *DT, int *SNOWNG, double *SNOTIME1,
   double *LVCOEF);

void            CANRES (double *SOLAR, double *CH, double *SFCTMP, double *Q2,
   double *SFCPRS, double *SMC, double *ZSOIL, int *NSOIL, double *SMCWLT,
   double *SMCREF, double *RSMIN, double *RC, double *PC, int *NROOT,
   double *Q2SAT, double *DQSDT2, double *TOPT, double *RSMAX, double *RGL,
   double *HS, double *XLAI, double *RCS, double *RCT, double *RCQ,
   double *RCSOIL, double *EMISSI);

void            CSNOW (double *SNCOND, double *DSNOW);

#ifdef _FLUX_PIHM_
void            DEVAP (double *EDIR, double *ETP1, double *SMC, double *ZSOIL,
   double *SHDFAC, double *SMCMAX, double *DKSAT, double *SMCDRY,
   double *SMCREF, double *SMCWLT, double *FXEXP);
#else
void            DEVAP (double *EDIR, double *ETP1, double *SMC, double *ZSOIL,
   double *SHDFAC, double *SMCMAX, double *BEXP, double *DKSAT, double *DWSAT,
   double *SMCDRY, double *SMCREF, double *SMCWLT, double *FXEXP);
#endif

#ifdef _FLUX_PIHM_
void            EVAPO (double *ETA1, double *SMC, int *NSOIL, double *CMC,
   double *ETP1, double *DT, double *ZSOIL, double *SH2O, double *SMCMAX,
   double *PC, double *SMCWLT, double *DKSAT, double *SMCREF, double *SHDFAC,
   double *CMCMAX, double *SMCDRY, double *CFACTR, double *EDIR, double *EC,
   double *ET, double *ETT, double *SFCTMP, double *Q2, int *NROOT,
   double *RTDIS, double *FXEXP);
#else
void            EVAPO (double *ETA1, double *SMC, int *NSOIL, double *CMC,
   double *ETP1, double *DT, double *ZSOIL, double *SH2O, double *SMCMAX,
   double *BEXP, double *PC, double *SMCWLT, double *DKSAT, double *DWSAT,
   double *SMCREF, double *SHDFAC, double *CMCMAX, double *SMCDRY,
   double *CFACTR, double *EDIR, double *EC, double *ET, double *ETT,
   double *SFCTMP, double *Q2, int *NROOT, double *RTDIS, double *FXEXP);
#endif

void            FAC2MIT (double *SMCMAX, double *FLIMIT);

#ifdef _FLUX_PIHM_
void            FRH2O (double *FREE, double *TKELV, double *SMC, double *SH2O,
   double *SMCMAX, double *SMCMIN, double *VGALPHA, double *VGBETA);
#else
void            FRH2O (double *FREE, double *TKELV, double *SMC, double *SH2O,
   double *SMCMAX, double *BEXP, double *PSIS);
#endif

#ifdef _FLUX_PIHM_
void            HRT (double *RHSTS, double *STC, double *SMC, double *SMCMAX,
   double *SMCMIN, int *NSOIL, double *ZSOIL, double *YY, double *ZZ1,
   double *TBOT, double *ZBOT, double *SH2O, double *DT, double *VGALPHA,
   double *VGBETA, double *F1, double *DF1, double *QUARTZ, double *CSOIL,
   double *AI, double *BI, double *CI, int *VEGTYP, int *ISURBAN);
#else
void            HRT (double *RHSTS, double *STC, double *SMC, double *SMCMAX,
   int *NSOIL, double *ZSOIL, double *YY, double *ZZ1, double *TBOT,
   double *ZBOT, double *PSISAT, double *SH2O, double *DT, double *BEXP,
   double *F1, double *DF1, double *QUARTZ, double *CSOIL, double *AI,
   double *BI, double *CI, int *VEGTYP, int *ISURBAN);
#endif

void            HSTEP (double *STCOUT, double *STCIN, double *RHSTS,
   double *DT, int *NSOIL, double *AI, double *BI, double *CI);

#ifdef _FLUX_PIHM_
void            NOPAC (double *ETP, double *ETA, double *PRCP, double *PCPDRP,
   double *SMC, double *SMCMAX, double *SMCMIN, double *SMCWLT,
   double *SMCREF, double *SMCDRY, double *CMC, double *CMCMAX, int *NSOIL,
   double *DT, double *SHDFAC, double *SBETA, double *Q2, double *T1,
   double *SFCTMP, double *T24, double *TH2, double *FDOWN, double *F1,
   double *EMISSI, double *SSOIL, double *STC, double *EPSCA, double *VGALPHA,
   double *VGBETA, double *MACKSAT, double *AREAF, int *NMACD, int *NWTBL,
   double *PC, double *RCH, double *RR, double *CFACTR, double *SH2O,
   double *FRZFACT, double *ZSOIL, double *DKSAT, double *TBOT, double *ZBOT,
   double *INF, double *RUNOFF2, double *RUNOFF3, double *EDIR, double *EC,
   double *ET, double *ETT, int *NROOT, double *RTDIS, double *QUARTZ,
   double *FXEXP, double *CSOIL, double *BETA, double *DRIP, double *DEW,
   double *FLX1, double *FLX3, int *VEGTYP, int *ISURBAN);
#else
void            NOPAC (double *ETP, double *ETA, double *PRCP, double *SMC,
   double *SMCMAX, double *SMCWLT, double *SMCREF, double *SMCDRY,
   double *CMC, double *CMCMAX, int *NSOIL, double *DT, double *SHDFAC,
   double *SBETA, double *Q2, double *T1, double *SFCTMP, double *T24,
   double *TH2, double *FDOWN, double *F1, double *EMISSI, double *SSOIL,
   double *STC, double *EPSCA, double *BEXP, double *PC, double *RCH,
   double *RR, double *CFACTR, double *SH2O, double *SLOPE, double *KDT,
   double *FRZFACT, double *PSISAT, double *ZSOIL, double *DKSAT,
   double *DWSAT, double *TBOT, double *ZBOT, double *RUNOFF1,
   double *RUNOFF2, double *RUNOFF3, double *EDIR, double *EC, double *ET,
   double *ETT, int *NROOT, double *RTDIS, double *QUARTZ, double *FXEXP,
   double *CSOIL, double *BETA, double *DRIP, double *DEW, double *FLX1,
   double *FLX3, int *VEGTYP, int *ISURBAN);
#endif

void            PENMAN (double *SFCTMP, double *SFCPRS, double *CH,
   double *T2V, double *TH2, double *PRCP, double *FDOWN, double *T24,
   double *SSOIL, double *Q2, double *Q2SAT, double *ETP, double *RCH,
   double *EPSCA, double *RR, int *SNOWNG, int *FRZGRA, double *DQSDT2,
   double *FLX2, double *EMISSI_IN, double *SNEQV, double *T1,
   double *SNCOVR);

void            REDPRM (GRID_TYPE * NOAH, LSM_STRUCT LSM, double *ZSOIL);

void            ROSR12 (double *P, double *A, double *B, double *C, double *D,
   double *DELTA, int *NSOIL);

#ifdef _FLUX_PIHM_
void            SHFLX (double *SSOIL, double *STC, double *SMC,
   double *SMCMAX, double *SMCMIN, int *NSOIL, double *T1, double *DT,
   double *YY, double *ZZ1, double *ZSOIL, double *TBOT, double *ZBOT,
   double *SMCWLT, double *SH2O, double *VGALPHA, double *VGBETA, double *F1,
   double *DF1, double *QUARTZ, double *CSOIL, int *VEGTYP, int *ISURBAN);
#else
void            SHFLX (double *SSOIL, double *STC, double *SMC,
   double *SMCMAX, int *NSOIL, double *T1, double *DT, double *YY,
   double *ZZ1, double *ZSOIL, double *TBOT, double *ZBOT, double *SMCWLT,
   double *PSISAT, double *SH2O, double *BEXP, double *F1, double *DF1,
   double *QUARTZ, double *CSOIL, int *VEGTYP, int *ISURBAN);
#endif

#ifdef _FLUX_PIHM_
void            SMFLX (double *SMC, int *NSOIL, double *CMC, double *DT,
   double *PRCP1, double *PCPDRP, double *ZSOIL, double *SH2O,
   double *FRZFACT, double *SMCMAX, double *SMCMIN, double *VGALPHA,
   double *VGBETA, double *MACKSAT, double *AREAF, int *NMACD, int *NWTBL,
   double *SMCWLT, double *DKSAT, double *SHDFAC, double *CMCMAX, double *INF,
   double *RUNOFF2, double *RUNOFF3, double *EDIR, double *EC, double *ET,
   double *DRIP);
#else
void            SMFLX (double *SMC, int *NSOIL, double *CMC, double *DT,
   double *PRCP1, double *ZSOIL, double *SH2O, double *SLOPE, double *KDT,
   double *FRZFACT, double *SMCMAX, double *BEXP, double *SMCWLT,
   double *DKSAT, double *DWSAT, double *SHDFAC, double *CMCMAX,
   double *RUNOFF1, double *RUNOFF2, double *RUNOFF3, double *EDIR,
   double *EC, double *ET, double *DRIP);
#endif

void            SNFRAC (double *SNEQV, double *SNUP, double *SALP,
   double *SNOWH, double *SNCOVR);

#ifdef _FLUX_PIHM_
void            SNKSRC (double *TSNSR, double *TAVG, double *SMC,
   double *SH2O, double *ZSOIL, int *NSOIL, double *SMCMAX, double *SMCMIN,
   double *VGALPHA, double *VGBETA, double *DT, int K, double *QTOT);
#else
void            SNKSRC (double *TSNSR, double *TAVG, double *SMC,
   double *SH2O, double *ZSOIL, int *NSOIL, double *SMCMAX, double *PSISAT,
   double *BEXP, double *DT, int K, double *QTOT);
#endif

#ifdef _FLUX_PIHM_
void            SNOPAC (double *ETP, double *ETA, double *PRCP, double *PRCPF,
   double *PCPDRP, int *SNOWNG, double *SMC, double *SMCMAX, double *SMCMIN,
   double *SMCWLT, double *SMCREF, double *SMCDRY, double *CMC,
   double *CMCMAX, int *NSOIL, double *DT, double *SBETA, double *DF1,
   double *Q2, double *T1, double *SFCTMP, double *T24, double *TH2,
   double *FDOWN, double *F1, double *SSOIL, double *STC, double *EPSCA,
   double *SFCPRS, double *VGALPHA, double *VGBETA, double *MACKSAT,
   double *AREAF, int *NMACD, int *NWTBL, double *PC, double *RCH, double *RR,
   double *CFACTR, double *SNCOVR, double *ESD, double *SNDENS, double *SNOWH,
   double *SH2O, double *FRZFACT, double *ZSOIL, double *DKSAT, double *TBOT,
   double *ZBOT, double *SHDFAC, double *INF, double *RUNOFF2,
   double *RUNOFF3, double *EDIR, double *EC, double *ET, double *ETT,
   int *NROOT, double *SNOMLT, double *RTDIS, double *QUARTZ, double *FXEXP,
   double *CSOIL, double *BETA, double *DRIP, double *DEW, double *FLX1,
   double *FLX2, double *FLX3, double *ESNOW, double *ETNS, double *EMISSI,
   double *RIBB, double *SOLDN, int *ISURBAN, int *VEGTYP);
#else
void            SNOPAC (double *ETP, double *ETA, double *PRCP, double *PRCPF,
   int *SNOWNG, double *SMC, double *SMCMAX, double *SMCWLT, double *SMCREF,
   double *SMCDRY, double *CMC, double *CMCMAX, int *NSOIL, double *DT,
   double *SBETA, double *DF1, double *Q2, double *T1, double *SFCTMP,
   double *T24, double *TH2, double *FDOWN, double *F1, double *SSOIL,
   double *STC, double *EPSCA, double *SFCPRS, double *BEXP, double *PC,
   double *RCH, double *RR, double *CFACTR, double *SNCOVR, double *ESD,
   double *SNDENS, double *SNOWH, double *SH2O, double *SLOPE, double *KDT,
   double *FRZFACT, double *PSISAT, double *ZSOIL, double *DWSAT,
   double *DKSAT, double *TBOT, double *ZBOT, double *SHDFAC, double *RUNOFF1,
   double *RUNOFF2, double *RUNOFF3, double *EDIR, double *EC, double *ET,
   double *ETT, int *NROOT, double *SNOMLT, double *RTDIS, double *QUARTZ,
   double *FXEXP, double *CSOIL, double *BETA, double *DRIP, double *DEW,
   double *FLX1, double *FLX2, double *FLX3, double *ESNOW, double *ETNS,
   double *EMISSI, double *RIBB, double *SOLDN, int *ISURBAN, int *VEGTYP);
#endif

void            SNOWPACK (double *ESD, double *DTSEC, double *SNOWH,
   double *SNDENS, double *TSNOW, double *TSOIL);

void            SNOWZ0 (double *SNCOVR, double *Z0, double *Z0BRD,
   double *SNOWH);

void            SNOW_NEW (double *TEMP, double *NEWSN, double *SNOWH,
   double *SNDENS);

#ifdef _FLUX_PIHM_
void            SRT (double *RHSTT, double *EDIR, double *ET, double *SH2O,
   double *SH2OA, int *NSOIL, double *PCPDRP, double *ZSOIL, double *DKSAT,
   double *SMCMAX, double *SMCMIN, double *VGALPHA, double *VGBETA,
   double *MACKSAT, double *AREAF, int *NMACD, double *INF, double *RUNOFF2,
   double *DT, double *SMCWLT, double *FRZX, double *SICE, double *AI,
   double *BI, double *CI);
#else
void            SRT (double *RHSTT, double *EDIR, double *ET, double *SH2O,
   double *SH2OA, int *NSOIL, double *PCPDRP, double *ZSOIL, double *DWSAT,
   double *DKSAT, double *SMCMAX, double *BEXP, double *RUNOFF1,
   double *RUNOFF2, double *DT, double *SMCWLT, double *SLOPE, double *KDT,
   double *FRZX, double *SICE, double *AI, double *BI, double *CI);
#endif

#ifdef _FLUX_PIHM_
void            SSTEP (double *SH2OOUT, double *SH2OIN, double *CMC,
   double *RHSTT, double *RHSCT, double *DT, int *NSOIL, double *SMCMAX,
   double *SMCMIN, double *CMCMAX, double *RUNOFF3, double *ZSOIL,
   double *SMC, double *SICE, double *AI, double *BI, double *CI);
#else
void            SSTEP (double *SH2OOUT, double *SH2OIN, double *CMC,
   double *RHSTT, double *RHSCT, double *DT, int *NSOIL, double *SMCMAX,
   double *CMCMAX, double *RUNOFF3, double *ZSOIL, double *SMC, double *SICE,
   double *AI, double *BI, double *CI);
#endif

void            TBND (double *TU, double *TB, double *ZSOIL, double *ZBOT,
   int K, int *NSOIL, double *TBND1);

#ifdef _FLUX_PIHM_
void            TDFCND (double *DF, double *SMC, double *QZ, double *SMCMAX,
   double *SMCMIN, double *SH2O);
#else
void            TDFCND (double *DF, double *SMC, double *QZ, double *SMCMAX,
   double *SH2O);
#endif

void            TMPAVG (double *TAVG, double *TUP, double *TM, double *TDN,
   double *ZSOIL, int *NSOIL, int K);

void            TRANSP (double *ET, int *NSOIL, double *ETP1, double *SMC,
   double *CMC, double *ZSOIL, double *SHDFAC, double *SMCWLT, double *CMCMAX,
   double *PC, double *CFACTR, double *SMCREF, double *SFCTMP, double *Q2,
   int *NROOT, double *RTDIS);

#ifdef _FLUX_PIHM_
void            WDFCND (double *WDF, double *WCND, double *SMC,
   double *SMCMAX, double *SMCMIN, double *VGALPHA, double *VGBETA,
   double *DKSAT, double *MACKSAT, double *AREAF, double *SICEMAX,
   double *DSMDZ, int *MACPORE);
#else
void            WDFCND (double *WDF, double *WCND, double *SMC,
   double *SMCMAX, double *BEXP, double *DKSAT, double *DWSAT,
   double *SICEMAX);
#endif

void            SFCDIF_off (double *ZLM, double *ZLM_WIND, double *Z0,
   double *THZ0, double *THLM, double *SFCSPD, double *CZIL, double *AKMS,
   double *AKHS, int *VEGTYP, int *ISURBAN, int *IZ0TLND);

#ifdef _FLUX_PIHM_
double          EFFKV (double KSATFUNC, double GRADY, double MACKV, double KV,
   double AREAF);
#endif

double          PSLMU (double ZZ);
double          PSLMS (double ZZ);
double          PSLHU (double ZZ);
double          PSLHS (double ZZ);
double          PSPMU (double XX);
double          PSPMS (double YY);
double          PSPHU (double XX);
double          PSPHS (double YY);

void            LSM_initialize (char *, Model_Data, Control_Data *,
   LSM_STRUCT);
void            LSM_read (char *, LSM_STRUCT);
void            LSM_initialize_output (char *, Model_Data, LSM_STRUCT,
   char *);
void            LSM_PrintInit (Model_Data, LSM_STRUCT, char *);
void            LSM_FreeData (Model_Data, LSM_STRUCT);

void            PIHM2Noah (realtype, realtype, Model_Data, LSM_STRUCT);
void            Noah2PIHM (Model_Data, LSM_STRUCT);

int             FindLayer (LSM_STRUCT, double);
double          mod (double, double);

double topo_radiation (double Sdir, double Sdif, double zenith, double azimuth180, double slope, double aspect, double *h_phi, double svf);
#endif
