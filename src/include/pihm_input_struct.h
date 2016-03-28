#ifndef PIHM_INPUT_STRUCT_HEADER
#define PIHM_INPUT_STRUCT_HEADER

typedef struct filename_struct
{
    char            riv[MAXSTRING];
    char            mesh[MAXSTRING];
    char            att[MAXSTRING];
    char            soil[MAXSTRING];
    char            geol[MAXSTRING];
    char            lc[MAXSTRING];
    char            forc[MAXSTRING];
    char            lai[MAXSTRING];
    char            ibc[MAXSTRING];
    char            para[MAXSTRING];
    char            calib[MAXSTRING];
    char            init[MAXSTRING];
#ifdef _NOAH_
    char            lsm[MAXSTRING];
    char            rad[MAXSTRING];
#endif
#ifdef _CYCLES_
    char            cycles[MAXSTRING];
    char            soilinit[MAXSTRING];
    char            crop[MAXSTRING];
    char            op[MAXOP][MAXSTRING];
#endif
} filename_struct;

typedef struct rivtbl_struct
{
    int             number;
    int            *fromnode;   /* Upstream Node no. */
    int            *tonode;     /* Dnstream Node no. */
    int            *down;       /* down stream segment */
    int            *leftele;    /* Left neighboring element */
    int            *rightele;   /* Right neighboring element */
    int            *shp;        /* shape type    */
    int            *matl;       /* material type */
    int            *bc;         /* BC type */
    int            *rsvr;
} rivtbl_struct;

typedef struct shptbl_struct
{
    int             number;
    double         *depth;      /* depth */
    int            *intrpl_ord; /* Interpolation order for river shape:
                                 * 1: rectangle,
                                 * 2: triangle,
                                 * 3: quadratic,
                                 * 4: cubic */
    double         *coeff;      /* Coefficient c in
                                 * D = c * pow(B / 2, interpOrd) */
} shptbl_struct;

typedef struct matltbl_struct
{
    int             number;     /* Number of River Bank/Bed Material */
    double         *rough;
    double         *cwr;        /* Weir Discharge Coefficient */
    double         *ksath;      /* Conductivity of river banks */
    double         *ksatv;      /* Conductivity of river bed */
    double         *bedthick;   /* thickeness of conductive river bed */
} matltbl_struct;

typedef struct meshtbl_struct
{
    int             numele;
    int             numnode;
    int           **node;
    int           **nabr;
    double         *x;
    double         *y;
    double         *zmin;
    double         *zmax;
} meshtbl_struct;

typedef struct atttbl_struct
{
    int            *soil;       /* soil type */
    int            *geol;
    int            *lc;         /* Land Cover type  */
    int           **bc;         /* Boundary condition type.
                                 * 0: Natural BC (no flow);
                                 * 1: Dirichlet BC;
                                 * 2:Neumann BC */
    int            *meteo;      /* precipitation (forcing) type */
    int            *lai;        /* LAI forcing type (0: use climatological
                                 * values; else: use user provided time
                                 *                                  * series */
    int            *source;     /* source (well) type */
    int            *macropore;
} atttbl_struct;

typedef struct soiltbl_struct
{
    int             number;     /* index */
    double         *silt;
    double         *clay;
    double         *om;
    double         *bd;
    double         *kinfv;
    double         *ksatv;      /* vertical saturated soil
                                 * conductivity */
    double         *ksath;
    double         *smcmax;     /* soil porosity */
    double         *smcmin;     /* soil moisture residual */
    double         *qtz;        /* ys: quartz content */
    double         *alpha;      /* soil curve parameter 1 */
    double         *beta;       /* soil curve parameter 2 */

    double         *areafh;     /* macroporous area fraction on
                                 * horizontal section */
    double         *areafv;
    double         *dmac;
    double         *smcref;
    double         *smcwlt;

    double          dinf;
    double          kmacv_ro;
    double          kmach_ro;
#ifdef _CYCLES_
    int            *totalLayers;
    double        **clay_lyr;
    double        **sand_lyr;
    double        **iom_lyr;
    double        **bd_lyr;
    double        **NO3_lyr;
    double        **NH4_lyr;
#endif
} soiltbl_struct;

typedef struct geoltbl_struct
{
    int             number;     /* index */
    double         *silt;
    double         *clay;
    double         *om;
    double         *bd;
    double         *ksath;      /* horizontal saturated geology
                                 * conductivity */
    double         *ksatv;      /* vertical saturated geology
                                 * conductivity */
    double         *smcmax;     /* geology porosity */
    double         *smcmin;     /* residual porosity */
    double         *alpha;      /* van genuchten parameter */
    double         *beta;       /* van genuchten parameter */
} geoltbl_struct;

typedef struct lctbl_struct
{
    int             number;     /* index */

    double         *laimax;     /* max lai */
    double         *laimin;     /* ys: min lai */
    double         *vegfrac;    /* canopy fracn */
    double         *albedomin;  /* ys: minimum albedo */
    double         *albedomax;  /* ys: maximum albedo */
    double         *emissmin;   /* ys: minimum emissivity */
    double         *emissmax;   /* ys: maximum emissivity */
    double         *z0min;      /* ys: minimum roughness length */
    double         *z0max;      /* ys: maximum roughness length */
    double         *hs;         /* ys: vapor pressure deficit stress
                                 * parameter */
    double         *snup;       /* ys */
    double         *rgl;        /* visible solar flux used in radiation
                                 * stress */
    double         *rsmin;      /* minimum stomatal resistance */
    double         *rough;      /* surface roughness factor  */
    double         *rzd;        /* rootzone depth */

    double          rsmax;      /* YS */
    int             bare;       /* YS */
    int             natural;
    double          cfactr;     /* YS */
    double          topt;       /* YS */
} lctbl_struct;

typedef struct tsdata_struct
{
    int             length;     /* length of time series */
    int            *ftime;
    double        **data;       /* 2D time series data */
    double         *value;
    double          zlvl_wind;
} tsdata_struct;

typedef struct forc_struct
{
    /* Forcing series */
    int             nbc;
    tsdata_struct  *bc;

    int             nmeteo;
    tsdata_struct  *meteo;

    int             nlai;
    tsdata_struct  *lai;

    int             nz0;
    tsdata_struct  *z0;

    int             nsource;
    tsdata_struct  *source;

    int             nmeltf;
    tsdata_struct  *meltf;

    int             nriverbc;
    tsdata_struct  *riverbc;
#ifdef _NOAH_
    int             nrad;
    tsdata_struct  *rad;
#endif
} forc_struct;

#ifdef _NOAH_
typedef struct noahtbl_struct
{
    double          sbeta;
    double          fxexp;
    double          csoil;
    double          salp;
    double          refdk;
    double          refkdt;
    double          frzk;
    double          zbot;
    double          tbot;
    double          smlow;
    double          smhigh;
    double          czil;
    double          lvcoef;
} noahtbl_struct;
#endif

#ifdef _CYCLES_
typedef struct agtbl_struct
{
    int            *op;
    int            *rotsz;
    int            *auto_N;
    int            *auto_P;
    int            *auto_S;
    int             nopfile;
    char            opfilen[MAXOP][MAXSTRING];
} agtbl_struct;

typedef struct croptbl_struct
{
    int             number;

    char          **cropName;
    int            *userSeedingDate;
    int            *userFloweringDate;
    int            *userMaturityDate;
    double         *userMaximumSoilCoverage;
    double         *userMaximumRootingDepth;
    double         *userExpectedYieldAvg;
    double         *userExpectedYieldMax;
    double         *userExpectedYieldMin;
    double         *userPercentMoistureInYield;
    double         *userFractionResidueStanding;
    double         *userFractionResidueRemoved;
    double         *userClippingBiomassThresholdUpper;
    double         *userClippingBiomassThresholdLower;
    double         *userClippingTiming;
    int            *userClippingDestiny;
    double         *userTranspirationMinTemperature;
    double         *userTranspirationThresholdTemperature;
    double         *userColdDamageMinTemperature;
    double         *userColdDamageThresholdTemperature;
    double         *userTemperatureBase;
    double         *userTemperatureOptimum;
    double         *userTemperatureMaximum;
    double         *userShootPartitionInitial;
    double         *userShootPartitionFinal;
    double         *userRadiationUseEfficiency;
    double         *userTranspirationUseEfficiency;
    double         *userHIx;
    double         *userHIo;    /* intercept harvest index */
    double         *userHIk;
    double         *userEmergenceTT;
    double         *userNMaxConcentration;
    double         *userNDilutionSlope;
    double         *userKc;
    int            *userAnnual;
    int            *userLegume;
    int            *userC3orC4;
    double         *userExtinctionCoefficient;
    double         *userPlantingDensity;
    int            *userClippingStart;
    int            *userClippingEnd;
    double         *LWP_StressOnset;
    double         *LWP_WiltingPoint;
    double         *transpirationMax;
} croptbl_struct;

typedef struct op_struct
{
    int             opYear;
    int             opDay;
    int             status;

    /* Planting Order */
    char            cropName[128];
    int             usesAutoIrrigation;
    int             usesAutoFertilization;
    int             plantID;
    double          plantingDensity;
    int             clippingStart;
    int             clippingEnd;

    /* Tillage */
    char            opToolName[MAXSTRING];
    double          opDepth;
    double          opSDR;
    double          opMixingEfficiency;
    char            cropNameT[128];
    double          fractionThermalTime;
    double          killEfficiency;
    int             grainHarvest;
    double          forageHarvest;

    /* Fixed Irrigation */
    double          opVolume;

    /* Fixed Fertilization */
    char            opSource[MAXSTRING];
    double          opMass;
    char            opForm[MAXSTRING];
    char            opMethod[MAXSTRING];
    int             opLayer;    /* Starting from 1 */
    double          opC_Organic;
    double          opC_Charcoal;
    double          opN_Organic;
    double          opN_Charcoal;
    double          opN_NH4;
    double          opN_NO3;
    double          opP_Organic;
    double          opP_Charcoal;
    double          opP_Inorganic;
    double          opK;
    double          opS;
} op_struct;

typedef struct autoirr_struct
{
    char            cropName[128];
    int             startDay;
    int             stopDay;
    double          waterDepletion;
    int             lastSoilLayer;
} autoirr_struct;
//
//typedef struct autoFertilizationStruct
//{
//    char            cropName[128];
//    int             startDay;
//    int             stopDay;
//    double          mass;
//    char            source[MAXSTRING];
//    char            form[MAXSTRING];
//    char            method[MAXSTRING];
//} autoFertilizationStruct;

typedef struct cropmgmt_struct
{
    int             yearsInRotation;
    int             adjustedYields;
    int             automaticNitrogen;
    int             automaticPhosphorus;
    int             automaticSulfur;
    int             rotationYear;

    op_struct      *FixedFertilization;
    int             numFertilization;

    op_struct      *FixedIrrigation;
    int             numIrrigation;

    op_struct      *Tillage;
    int             numTillage;
    double         *tillageFactor;

    op_struct      *plantingOrder;
    int             totalCropsPerRotation;

    autoirr_struct *autoIrrigation;
    int             numAutoIrrigation;

    int             usingAutoIrr;
    int             usingAutoFert;
} cropmgmt_struct;

typedef struct mgmttbl_struct
{
    int             number;
    cropmgmt_struct *cropmgmt;
} mgmttbl_struct;
#endif
#endif
