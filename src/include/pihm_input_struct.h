#ifndef PIHM_INPUT_STRUCT_HEADER
#define PIHM_INPUT_STRUCT_HEADER

/*****************************************************************************
 * Input file names
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * riv                      char[]      river input file name
 * mesh                     char[]      mesh structure file name
 * att                      char[]      element attribute file name
 * soil                     char[]      soil property file name
 * geol                     char[]      geology property file name
 * lc                       char[]      land cover property file name
 * meteo                    char[]      meteorological forcing file name
 * lai                      char[]      lai forcing file name
 * bc                       char[]      boundary condition file name
 * para                     char[]      control parameter file name
 * calib                    char[]      calibration file name
 * ic                       char[]      initial condition file name
 * lsm                      char[]      land surface module control file name
 * rad                      char[]      radiation forcing file name
 * bgc                      char[]      bgc module control file name
 * co2                      char[]      co2 forcing file name
 * ndep                     char[]      nitrogen deposition forcing file name
 * bgcic                    char[]      bgc module initial condition file name
 * cyclesic                 char[]      cycles module initial condition file
 *                                        name
 ****************************************************************************/
typedef struct filename_struct
{
    char            riv[MAXSTRING];
    char            mesh[MAXSTRING];
    char            att[MAXSTRING];
    char            soil[MAXSTRING];
    char            geol[MAXSTRING];
    char            lc[MAXSTRING];
    char            meteo[MAXSTRING];
    char            lai[MAXSTRING];
    char            bc[MAXSTRING];
    char            para[MAXSTRING];
    char            calib[MAXSTRING];
    char            ic[MAXSTRING];
#ifdef _NOAH_
    char            lsm[MAXSTRING];
    char            rad[MAXSTRING];
#endif
#ifdef _CYCLES_
    char            cycles[MAXSTRING];
    char            soilinit[MAXSTRING];
    char            crop[MAXSTRING];
    char            op[MAXOP][MAXSTRING];
    char            cyclesic[MAXSTRING];
#endif
#ifdef _BGC_
    char            bgc[MAXSTRING];
    char            co2[MAXSTRING];
    char            ndep[MAXSTRING];
    char            bgcic[MAXSTRING];
#endif
} filename_struct;

/*****************************************************************************
 * River input structure
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * number                   int*        number of river segments
 * fromnode                 int*        upstream node id
 * tonode                   int*        downstream node id
 * down                     int*        downstream channel id
 * leftele                  int*        left bank id
 * rightele                 int*        right bank id
 * shp                      int*        river shape type
 * matl                     int*        material type
 * bc                       int*        boundary condition type
 * rsvr                     int*        reservoir type
 ****************************************************************************/
typedef struct rivtbl_struct
{
    int            *fromnode;
    int            *tonode;
    int            *down;
    int            *leftele;
    int            *rightele;
    int            *shp;
    int            *matl;
    int            *bc;
    int            *rsvr;
} rivtbl_struct;

/*****************************************************************************
 * River shape parameters
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * number                   int         number of shape types
 * depth                    double*     river channel depth
 * intrpl_ord               int*        interpolation order (shape of channel)
 *                                        1: rectangle
 *                                        2: triangle
 *                                        3: quadratic
 *                                        4: cubic
 * coeff                    double*     width coefficient
 ****************************************************************************/
typedef struct shptbl_struct
{
    int             number;
    double         *depth;
    int            *intrpl_ord;
    double         *coeff;
} shptbl_struct;

/*****************************************************************************
 * River channel matierial parameters
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * number                   int         number of bank/bed material types
 * rough                    double*     river channel roughness [s m-1/3]
 * cwr                      double*     discharge coefficient [-]
 * ksath                    double*     bank hydraulic conductivity [m s-1]
 * ksatv                    double*     bed hydraulic conductivity [m s-1]
 * bedthick                 double*     bed thickness [m]
 ****************************************************************************/
typedef struct matltbl_struct
{
    int             number;
    double         *rough;
    double         *cwr;
    double         *ksath;
    double         *ksatv;
    double         *bedthick;
} matltbl_struct;

/*****************************************************************************
 * Mesh structure
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * numnode                  int         number of nodes
 * node                     int**       nodes of element
 * nabr                     int**       neighbors of element
 * x                        double*     x of node [m]
 * y                        double*     y of node [m]
 * zmin                     double*     bedrock elevation of node [m]
 * zmax                     double*     surface elevation of node [m]
 ****************************************************************************/
typedef struct meshtbl_struct
{
    int             numnode;
    int           **node;
    int           **nabr;
    double         *x;
    double         *y;
    double         *zmin;
    double         *zmax;
} meshtbl_struct;

/*****************************************************************************
 * Element attribute
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * soil                     int*        element soil type
 * geol                     int*        element geology type
 * lc                       int*        element land cover type
 * bc                       int**       element boudnary condition type
 * meteo                    int*        element meteorological forcing type
 * lai_type                 int*        element leaf area index forcing type
 *                                        0: use climatological values;
 *                                        else: use forcing file
 * source                   int*        element source forcing type
 ****************************************************************************/
typedef struct atttbl_struct
{
    int            *soil;
    int            *geol;
    int            *lc;
    int           **bc;
    int            *meteo;
    int            *lai;
    int            *source;
} atttbl_struct;

/*****************************************************************************
 * Soil parameter
 * ---------------------------------------------------------------------------
 *
 * Variables                Type        Description
 * ==========               ==========  ====================
 * number                   int         number of soil types
 * silt                     double*     silt percentage [%]
 * clay                     double*     clay percentage [%]
 * om                       double*     organic matter percentage [%]
 * bd                       double*     bulk density [g cm-3]
 * kinfv                    double*     saturated infiltration conductivity
 *                                        [m s-1]
 * ksatv                    double*     vertical saturated hydraulic
 *                                        conductivity [m s-1]
 * ksath                    double*     horizontal saturated hydraulic
 *                                        conductivity [m s-1]
 * smcmax                   double*     maximum soil moisture content [m3 m-3]
 * smcmin                   double*     residual soil moisture content
 *                                        [m3 m-3]
 * smcwlt                   double*     wilting point [m3 m-3]
 * smcref                   double*     soil moisture threshold where
 *                                        transpiration begins to stress
 *                                        [m3 m-3]
 * qtz                      double*     soil quartz content [-]
 * alpha                    double*     alpha from van Genuchten eqn [m-1]
 * beta                     double*     beta (n) from van Genuchten eqn [-]
 * areafh                   double*     macropore area fraction on a
 *                                        horizontal cross-section [m2 m-2]
 * areafv                   double*     macropore area fraction on a vertical
 *                                        cross-section [m2 m-2]
 * dmac                     double*     macropore depth [m]
 * dinf                     double      depth from ground surface accross which
 *                                        head gradient is calculated for
 *                                        infiltration [m]
 * kmacv_ro                 double      ratio between vertical macropore
 *                                        hydarulic conductivity and vertical
 *                                        saturated infiltration hydarulic
 *                                        conductivity [-]
 * kmach_ro                 double      ratio between horizontal macropore
 *                                        hydarulic conductivity and
 *                                        horizontal saturated hydarulic
 *                                        conductivity [-]
 ****************************************************************************/
typedef struct soiltbl_struct
{
    int             number;
    double         *silt;
    double         *clay;
    double         *om;
    double         *bd;
    double         *kinfv;
    double         *ksatv;
    double         *ksath;
    double         *smcmax;
    double         *smcmin;
    double         *smcwlt;
    double         *smcref;
    double         *qtz;
    double         *alpha;
    double         *beta;
    double         *areafh;
    double         *areafv;
    double         *dmac;
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

/*****************************************************************************
 * Geology parameter
 * ---------------------------------------------------------------------------
 *
 * Variables                Type        Description
 * ==========               ==========  ====================
 * number                   int         number of soil types
 * silt                     double*     silt percentage [%]
 * clay                     double*     clay percentage [%]
 * om                       double*     organic matter percentage [%]
 * bd                       double*     bulk density [g cm-3]
 * ksath                    double*     horizontal saturated hydraulic
 *                                        conductivity [m s-1]
 * ksatv                    double*     vertical saturated hydraulic
 *                                        conductivity [m s-1]
 * smcmax                   double*     maximum soil moisture content [m3 m-3]
 * smcmin                   double*     residual soil moisture content
 *                                        [m3 m-3]
 * alpha                    double*     alpha from van Genuchten eqn [m-1]
 * beta                     double*     beta (n) from van Genuchten eqn [-]
 ****************************************************************************/
typedef struct geoltbl_struct
{
    int             number;
    double         *silt;
    double         *clay;
    double         *om;
    double         *bd;
    double         *ksath;
    double         *ksatv;
    double         *smcmax;
    double         *smcmin;
    double         *alpha;
    double         *beta;
} geoltbl_struct;

/*****************************************************************************
 * Land cover parameters
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * number                   int         number of land cover types
 * laimax                   double      maximum LAI across all seasons for a
 *                                        vegetation type [m2 m-2]
 * laimin                   double      minimum LAI across all seasons for a
 *                                        vegetation type [m2 m-2]
 * vegfrac                  double      areal fractional coverage of green
 *                                        vegetation (0.0-1.0) [-]
 * albedomin                double      minimum background albedo [-]
 * albedomax                double      maximum background albedo [-]
 * emissmin                 double      minimum emissivity [-]
 * emissmax                 double      maximum emissivity [-]
 * z0min                    double      minimum roughness length [m]
 * z0max                    double      maximum roughness length [m]
 * hs                       double      parameter used in vapor pressure
 *                                        deficit function [-]
 * snup                     double      threshold snow depth (in water
 *                                        equivalent) that implies 100% snow
 *                                        cover [m]
 * rgl                      double      reference incoming solar flux for
 *                                        photosynthetically active canopy
 *                                        [W m-2]
 * rsmin                    double      minimum canopy resistance [s m-1]
 * rough                    double      surface roughness (Manning's n)
 *                                        [s m-1/3]
 * rzd                      double      rooting depth [m]
 * rsmax                    double      cuticular resistance [s m-1]
 * bare                     int         the land-use category representing
 *                                        bare ground
 * natural                  int         the land-use category representing
 *                                        non-urban portion of urban land-use
 *                                        points
 * cfactr                   double      parameter used in the canopy
 *                                        inteception calculation [-]
 * topt                     double      optimum transpiration air temperature
 *                                        [K]
 ****************************************************************************/
typedef struct lctbl_struct
{
    int             number;
    double         *laimax;
    double         *laimin;
    double         *vegfrac;
    double         *albedomin;
    double         *albedomax;
    double         *emissmin;
    double         *emissmax;
    double         *z0min;
    double         *z0max;
    double         *hs;
    double         *snup;
    double         *rgl;
    double         *rsmin;
    double         *rough;
    double         *rzd;
    double          rsmax;
    int             bare;
    int             natural;
    double          cfactr;
    double          topt;
} lctbl_struct;

/*****************************************************************************
 * Time series data structure
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * length                   int         length of time series
 * ftime                    double*     forcing time
 * data                     double**    forcing values at forcing time
 * value                    double*     forcing values at model time t
 * zlvl_wind                double      height above groundof wind
 *                                        observations [m]
 ****************************************************************************/
typedef struct tsdata_struct
{
    int             length;
    int            *ftime;
    double        **data;
    double         *value;
    double          zlvl_wind;
} tsdata_struct;

/*****************************************************************************
 * Forcing structure
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * nbc                      int         number of boundary condition series
 * bc                       tsdata_struct*
 *                                      boundary condition time series
 * nmeteo                   int         number of meteorological forcing
 *                                        series
 * meteo                    tsdata_struct*
 *                                      meteorological forcing series
 * nlai                     int         number of lai series
 * lai                      tsdata_struct*
 *                                      lai forcing series
 * nsource                  int         number of source forcing series
 * source                   tsdata_struct*
 *                                      source forcing series
 * nriverbc                 int         number of river boundary conditions
 * riverbc                  tsdata_struct*
 *                                      river boundary condition series
 * nrad                     int         number of radiation forcing series
 * rad                      tsdata_struct*
 *                                      radiation forcing series
 * co2                      tsdata_struct*
 *                                      co2 forcing series
 * ndep                     tsdata_struct*
 *                                      nitrogen deposition forcing series
 ****************************************************************************/
typedef struct forc_struct
{
    int             nbc;
    tsdata_struct  *bc;
    int             nmeteo;
    tsdata_struct  *meteo;
    int             nlai;
    tsdata_struct  *lai;
    int             nsource;
    tsdata_struct  *source;
    int             nriverbc;
    tsdata_struct  *riverbc;
#ifdef _NOAH_
    int             nrad;
    tsdata_struct  *rad;
#endif
#ifdef _BGC_
    tsdata_struct  *co2;
    tsdata_struct  *ndep;
#endif
} forc_struct;

#ifdef _NOAH_

/*****************************************************************************
 * Land surface parameters
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * sbeta                    double      parameter used to calculate vegetation
 *                                        effect on soil heat [-]
 * fxexp                    double      soil evaporation exponent used in
 *                                        direct evaporation [-]
 * csoil                    double      soil heat capacity [J m-3 K-1]
 * salp                     double      shape parameter of distribution
 *                                        function of snow cover [-]
 * frzk                     double      frozen ground parameter [-]
 * zbot                     double      depth of lower boundary soil
 *                                        temperature [m]
 * tbot                     double      bottom soil temperature (local yearly-
 *                                        mean sfc air temperature) [K]
 * czil                     double      zilitinkevich constant [-]
 * lvcoef                   double      parameter controls surface snow albedo
 *                                        in the presence of snowcover [-]
 ****************************************************************************/
typedef struct noahtbl_struct
{
    double          sbeta;
    double          fxexp;
    double          csoil;
    double          salp;
    double          frzk;
    double          zbot;
    double          tbot;
    double          czil;
    double          lvcoef;
} noahtbl_struct;
#endif

#ifdef _BGC_
/*****************************************************************************
 * Ecophysiological parameters
 * ---------------------------------------------------------------------------
 * Variables                Type        Description
 * ==========               ==========  ====================
 * woody                    int*        flag: 1 = woody, 0 = non-woody
 * evergreen                int*        flag: 1 = evergreen, 0 = deciduous
 * c3_flag                  int*        flag: 1 = C3,  0 = C4
 * phenology_flag           int*        flag: 1 = phenology model, 0 = user
 *                                        defined
 * onday                    int*        day of year when leaves on
 * offday                   int*        day of year yearday leaves off
 * transfer_days            int*        growth period for transfer [day]
 * litfall_days             int*        growth period for litfall [day]
 * leaf_turnover            double*     annual leaf turnover fraction [yr-1]
 * froot_turnover           double*     annual fine root turnover fraction
 *                                        [yr-1]
 * livewood_turnover        double*     annual live wood turnover fraction
 *                                        [yr-1]
 * daily_mortality_turnover double*     daily mortality turnover [day-1]
 * daily_fire_turnover      double*     daily fire turnover [day-1]
 * alloc_frootc_leafc       double*     new fine root C to new leaf C [-]
 * alloc_newstemc_newleafc  double*     new stem C to new leaf C [-]
 * alloc_newlivewoodc_newwoodc
 *                          double*     new livewood C:new wood C [-]
 * alloc_crootc_stemc       double*     new live croot C to new live stem C
 *                                        [-]
 * alloc_prop_curgrowth     double*     daily allocation to current growth [-]
 * avg_proj_sla             double*     canopy average projected SLA
 *                                        [m2 kgC-1]
 * sla_ratio                double*     ratio of shaded to sunlit projected
 *                                        SLA [-]
 * lai_ratio                double*     ratio of (all-sided LA / one-sided LA)
 *                                        [-]
 * ext_coef                 double*     canopy light extinction coefficient
 *                                        [-]
 * flnr                     double*     leaf N in Rubisco [kgNRub kgNleaf-1]
 * psi_open                 double*     psi at start of conductance reduction
 *                                        [MPa]
 * psi_close                double*     psi at complete conductance reduction
 *                                        [MPa]
 * vpd_open                 double*     vpd at start of conductance reduction
 *                                        [Pa]
 * vpd_close                double*     vpd at complete conductance reduction
 *                                        [Pa]
 * froot_cn                 double*     C:N for fine roots [kgC kgN-1]
 * leaf_cn                  double*     C:N for leaves [kgC kgN-1]
 * livewood_cn              double*     C:N for live wood [kgC kgN-1]
 * deadwood_cn              double*     C:N for dead wood [kgC kgN-1]
 * leaflitr_cn              double*     constant C:N for leaf litter
 *                                        [kgC kgN-1]
 * leaflitr_flab            double*     leaf litter labile fraction [-]
 * leaflitr_fucel           double*     leaf litter unshielded cellulose
 *                                        fraction [-]
 * leaflitr_fscel           double*     leaf litter shielded cellulose
 *                                        fraction [-]
 * leaflitr_flig            double*     leaf litter lignin fraction [-]
 * frootlitr_flab           double*     fine root litter labile fraction [-]
 * frootlitr_fucel          double*     fine root litter unshielded cellulose
 *                                        fraction [-]
 * frootlitr_fscel          double*     fine root litter shielded cellulose
 *                                        fraction [-]
 * frootlitr_flig           double*     fine root litter lignin fraction [-]
 * deadwood_fucel           double*     dead wood unshileded cellulose
 *                                        fraction [-]
 * deadwood_fscel           double*     dead wood shielded cellulose fraction
 *                                        [-]
 * deadwood_flig            double*     dead wood lignin fraction [-]
 ****************************************************************************/
typedef struct epctbl_struct
{
    int            *woody;
    int            *evergreen;
    int            *c3_flag;
    int            *phenology_flag;
    int            *onday;
    int            *offday;
    int            *transfer_days;
    int            *litfall_days;
    double         *leaf_turnover;
    double         *froot_turnover;
    double         *livewood_turnover;
    double         *daily_mortality_turnover;
    double         *daily_fire_turnover;
    double         *alloc_frootc_leafc;
    double         *alloc_newstemc_newleafc;
    double         *alloc_newlivewoodc_newwoodc;
    double         *alloc_crootc_stemc;
    double         *alloc_prop_curgrowth;
    double         *avg_proj_sla;
    double         *sla_ratio;
    double         *lai_ratio;
    double         *ext_coef;
    double         *flnr;
    double         *psi_open;
    double         *psi_close;
    double         *vpd_open;
    double         *vpd_close;
    double         *froot_cn;
    double         *leaf_cn;
    double         *livewood_cn;
    double         *deadwood_cn;
    double         *leaflitr_cn;
    double         *leaflitr_flab;
    double         *leaflitr_fucel;
    double         *leaflitr_fscel;
    double         *leaflitr_flig;
    double         *frootlitr_flab;
    double         *frootlitr_fucel;
    double         *frootlitr_fscel;
    double         *frootlitr_flig;
    double         *deadwood_fucel;
    double         *deadwood_fscel;
    double         *deadwood_flig;
} epctbl_struct;
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
    double         *userFloweringTT;
    double         *userMaturityTT;
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

typedef struct plant_struct
{
    /* Planting */
    int             opYear;
    int             opDay;
    char            cropName[128];
    int             usesAutoIrrigation;
    int             usesAutoFertilization;
    int             plantID;
    double          plantingDensity;
    int             clippingStart;
    int             clippingEnd;
} plant_struct;

typedef struct tillage_struct
{
    /* Tillage */
    int             opYear;
    int             opDay;
    char            opToolName[MAXSTRING];
    double          opDepth;
    double          opSDR;
    double          opMixingEfficiency;
    char            cropNameT[128];
    double          fractionThermalTime;
    double          killEfficiency;
    int             grainHarvest;
    double          forageHarvest;
} tillage_struct;

typedef struct fixirr_struct
{
    /* Fixed Irrigation */
    int             opYear;
    int             opDay;
    double          opVolume;
} fixirr_struct;

typedef struct fixfert_struct
{
    /* Fixed Fertilization */
    int             opYear;
    int             opDay;
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
} fixfert_struct;

typedef struct autoirr_struct
{
    char            cropName[MAXSTRING];
    int             startDay;
    int             stopDay;
    double          waterDepletion;
    int             lastSoilLayer;
} autoirr_struct;

typedef struct cropmgmt_struct
{
    int             yearsInRotation;
    //int             adjustedYields;
    int             automaticNitrogen;
    int             automaticPhosphorus;
    int             automaticSulfur;
    int             rotationYear;

    fixfert_struct *FixedFertilization;
    int             numFertilization;

    fixirr_struct  *FixedIrrigation;
    int             numIrrigation;

    tillage_struct *Tillage;
    int             numTillage;
    double          tillageFactor[MAXLYR];

    plant_struct   *plantingOrder;
    int             totalCropsPerRotation;

    autoirr_struct *autoIrrigation;
    int             numAutoIrrigation;

    int            *op_status[4];

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
