#ifndef PIHM_CONST_HEADER
#define PIHM_CONST_HEADER

// Physical parameters
#define PI                      3.14159265  // pi
#define DAYINSEC                86400       // number of seconds in a day
#define GRAV                    9.80665     // gravity constant (m s-2)
#define CP                      1004.0      // specific heat capacity of air (J kg-1 K-1)
#define RHOH2O                  1000.0      // water density (kg m-3)
#define LVH2O                   2.501E6     // latent heat of vaporization (J kg-1)
#define SIGMA                   5.67E-8     // Stefan-Boltzmann constant (W m-2 K-4)
#define RD                      287.04      // gas constant for dry air (J kg-1 K-1)
#define RV                      461.5       // gas constant for water vapor (J kg-1 K-1)
#define CPH2O                   4.218E3     // specific heat capacity of water (J kg-1 K-1)
#define CPICE                   2.106E3     // specific heat capacity of ice (J kg-1 K-1)
#define LSUBF                   3.335E5     // latent heat of fusion (J kg-1)
#define EMISSI_S                0.95        // emissivity of snow (-)
#define TFREEZ                  273.15      // freezing point (K)
#define LSUBS                   2.83E6      // latent heat of sublimation (J kg-1)
#define ICEDENS                 0.9         // glacier ice density (g cm-3, or dimensionless fraction of H2O density

#define F_OK                    0

// Simulation mode
#define NORMAL_MODE             0
#define SPINUP_MODE             1
#define ACC_SPINUP_MODE         2

// Default bad value
#define BADVAL                  -999

// Model steps
#define HYDROL_STEP             0
#define LS_STEP                 1
#define CN_STEP                 2

// Maximum number of output files
#define MAXPRINT                1024

// Meteorological forcing related
#define NUM_METEO_VAR           7           // number of meteorological forcing variables
#define PRCP_TS                 0           // index of precipitation forcing
#define SFCTMP_TS               1           // index of air temperature forcing
#define RH_TS                   2           // index of RH forcing
#define SFCSPD_TS               3           // index of wind speed forcing
#define SOLAR_TS                4           // index of solar radiation forcing
#define LONGWAVE_TS             5           // index of longwave radiation forcing
#define PRES_TS                 6           // index of surface pressure forcing

// Radiation forcing variables
#define UNIF_SOL                0           // use global solar radiation
#define TOPO_SOL                1           // use topographic solar radiation
#define SOLDIR_TS               0           // index of direct solar radiation forcing
#define SOLDIF_TS               1           // index of diffused solar radiation forcing

// Number of edges of an element
#define NUM_EDGE                3

#define WS_ZMAX                 0
#define WS_ZMIN                 1
#define WS_AREA                 2

// Hydrology parameters
#define PSIMIN                  -70.0       // minimum psi allowed (m)
#define DEPRSTG                 1E-4        // depression storage (m)
#define GRADMIN                 5E-8        // minimum hydraulic gradient (m m-1)
#define SATMIN                  0.1         // minimum saturation ratio (-)
#define RIVDPTHMIN              0.05        // minimum river depth (m)
#define RIVGRADMIN              0.05        // minimum river hydraulic gradient (m m-1)
#define CMCFACTR                2E-4        // canopy water capacity per LAI (m)
#define SH2OMIN                 0.02        // minimum swc (m3 m-3)

// Maximum of soil layers in Flux-PIHM
#define MAXLYR                  11

// Land cover parameters
#define NLCTYPE                 40          // number of land cover types

// Land cover types
#define IGBP_ENF                1
#define IGBP_EBF                2
#define IGBP_DNF                3
#define IGBP_DBF                4
#define IGBP_MIXF               5
#define IGBP_CLOSE_SHRUB        6
#define IGBP_OPEN_SHRUB         7
#define IGBP_WOODY_SAVANNA      8
#define IGBP_SAVANNA            9
#define IGBP_GRASS              10
#define IGBP_WETLAND            11
#define IGBP_CROP               12
#define IGBP_URBAN_BUILDUP      13
#define IGBP_CROP_NATURAL       14
#define IGBP_SNOW_ICE           15
#define IGBP_BARREN             16
#define IGBP_WATER              17
#define IGBP_UNCLASSIFIED1      18
#define IGBP_UNCLASSIFIED2      19
#define IGBP_UNCLASSIFIED3      20
#define NLCD40_WATER            21
#define NLCD40_SNOW_ICE         22
#define NLCD40_DEVELOPED_OPEN   23
#define NLCD40_DEVELOPED_LOW    24
#define NLCD40_DEVELOPED_MID    25
#define NLCD40_DEVELOPED_HIGH   26
#define NLCD40_BARREN           27
#define NLCD40_DECIDUOUS        28
#define NLCD40_EVERGREEN        29
#define NLCD40_MIXF             30
#define NLCD40_DWARF_SCRUB      31
#define NLCD40_SHRUB            32
#define NLCD40_GRASS            33
#define NLCD40_SEDGE            34
#define NLCD40_LICHENS          35
#define NLCD40_MOSS             36
#define NLCD40_PASTURE          37
#define NLCD40_CROP             38
#define NLCD40_WOODY_WETLAND    39
#define NLCD40_HERB_WETLAND     40

// Soil textures
#define SAND                    0
#define LOAMY_SAND              1
#define SANDY_LOAM              2
#define LOAM                    3
#define SILT_LOAM               4
#define SILT                    5
#define SANDY_CLAY_LOAM         6
#define CLAY_LOAM               7
#define SILTY_CLAY_LOAM         8
#define SANDY_CLAY              9
#define SILTY_CLAY              10
#define CLAY                    11

// Number of river fluxes of a river segment
# define NUM_RIVFLX             6

// River fluxes
#define UPSTREAM                0
#define DOWNSTREAM              1
#define SURF_LEFT               2
#define SURF_RIGHT              3
#define AQUIFER_LEFT            4
#define AQUIFER_RIGHT           5

// River boundary condition types
#define OUTLET_DIRICHLET        -1
#define OUTLET_NEUMANN          -2
#define ZERO_DPTH_GRAD          -3
#define CRIT_DPTH               -4

// Element boundary condition types
#define NO_FLOW                 0
#define DIRICHLET               1
#define NEUMANN                 2

// River segment interpolation order
#define RECTANGLE               1
#define TRIANGLE                2
#define QUADRATIC               3
#define CUBIC                   4

// Approximation
#define KINEMATIC               1
#define DIFF_WAVE               2

// Initialization type
#define RELAX                   0
#define RST_FILE                1

// Average flux
#define SUM                     0
#define AVG                     1

// Interpolate forcing time series
#define NO_INTRPL               0
#define INTRPL                  1

// Maximum allowable difference between simulation cycles in subsurface water storage at steady-state (m)
#define SPINUP_W_TOLERANCE      0.01

// Ecosystem constants
#define RAD2PAR                 0.45        // ratio PAR / SWtotal (-)
#define EPAR                    4.55        // (umol/J) PAR photon energy ratio
#define SOIL1_CN                12.0        // C:N for fast microbial recycling pool
#define SOIL2_CN                12.0        // C:N for slow microbial recycling pool
#define SOIL3_CN                10.0        // C:N for recalcitrant SOM pool (humus)
#define SOIL4_CN                10.0        // C:N for recalcitrant SOM pool (humus)
#define GRPERC                  0.3         // growth resp per unit of C grown (-)
#define GRPNOW                  1.0         // proportion of storage growth resp at fixation (-)
#define PPFD50                  75.0        // PPFD for 1/2 stomatal closure (umol m-2 s-1)
#define DENITRIF_PROPORTION     0.01        // fraction of mineralization to volatile
# define MOBILEN_PROPORTION     0.1         // fraction mineral N avail for leaching

// Respiration fractions for fluxes between compartments (-)
#define RFL1S1                  0.39        // transfer from litter 1 to soil 1
#define RFL2S2                  0.55        // transfer from litter 2 to soil 2
#define RFL4S3                  0.29        // transfer from litter 4 to soil 3
#define RFS1S2                  0.28        // transfer from soil 1 to soil 2
#define RFS2S3                  0.46        // transfer from soil 2 to soil 3
#define RFS3S4                  0.55        // transfer from soil 3 to soil 4

// Base decomposition rate constants (day -1)
#define KL1_BASE                0.7         // labile litter pool
#define KL2_BASE                0.07        // cellulose litter pool
#define KL4_BASE                0.014       // lignin litter pool
#define KS1_BASE                0.07        // fast microbial recycling pool
#define KS2_BASE                0.014       // medium microbial recycling pool
#define KS3_BASE                0.0014      // slow microbial recycling pool
#define KS4_BASE                0.0001      // recalcitrant SOM (humus) pool
#define KFRAG_BASE              0.001       // physical fragmentation of coarse woody debris

// Decomposition acceleration terms
#define KS1_ACC                 1.0
#define KS2_ACC                 1.0
#define KS3_ACC                 5.0
#define KS4_ACC                 70.0

// This constant determines the lower limit of state variables before they are set to 0.0 to control rounding and
// overflow errors
#define CRIT_PREC               1E-20

// This constant is used in if conditions where floating point values are compared
#define FLT_COND_TOL            1E-10

// Maximum allowable trend in slow soil carbon at steady-state (kgC m-2 yr-1)
#define SPINUP_C_TOLERANCE      0.0005

// Allocation parameters
#define DAYSNDEPLOY             365.0
#define DAYSCRECOVER            365.0
#define BULK_DENITRIF_PROPORTION 0.5

// Output variables
#define YEARLY_OUTPUT           -1
#define MONTHLY_OUTPUT          -2
#define DAILY_OUTPUT            -3
#define HOURLY_OUTPUT           -4

#define FORCING_BC              0
#define FORCING_METEO           1
#define FORCING_RAD             2
#define FORCING_LAI             3
#define FORCING_RIVERBC         4
#define FORCING_CO2             5
#define FORCING_NDEP            6

// Output variable types
enum output_var{
    SURF_CTRL,
    UNSAT_CTRL,
    GW_CTRL,
    STAGE_CTRL,
    SNOW_CTRL,
    CMC_CTRL,
    INFIL_CTRL,
    RECHARGE_CTRL,
    EC_CTRL,
    ETT_CTRL,
    EDIR_CTRL,
    RIVFLX0_CTRL,
    RIVFLX1_CTRL,
    RIVFLX2_CTRL,
    RIVFLX3_CTRL,
    RIVFLX4_CTRL,
    RIVFLX5_CTRL,
    SUBFLX_CTRL,
    SURFFLX_CTRL,
    // System
    IC_CTRL,
    WB_CTRL,
    // Noah
    T1_CTRL,
    STC_CTRL,
    SMC_CTRL,
    SH2O_CTRL,
    SNOWH_CTRL,
    ALBEDO_CTRL,
    LE_CTRL,
    SH_CTRL,
    G_CTRL,
    ETP_CTRL,
    ESNOW_CTRL,
    ROOTW_CTRL,
    SOILM_CTRL,
    SOLAR_CTRL,
    CH_CTRL,
    // BGC or Cycles
    LAI_CTRL,
    N_PROFILE_CTRL,
    N_RIVER_CTRL,
    LEACHING_CTRL,
    DENITRIF_CTRL,
    // Cycles
    YIELD_CTRL,
    BIOMASS_CTRL,
    RADNINTCP_CTRL,
    WATER_STS_CTRL,
    N_STS_CTRL,
    CROP_TR_CTRL,
    CROP_POTTR_CTRL,
    NITRIF_CTRL,
    IMMOBIL_CTRL,
    MINERAL_CTRL,
    VOLATIL_CTRL,
    SOC_CTRL,
    N2O_CTRL,
    N_HARVEST_CTRL,
    N_INPUT_CTRL,
    // BGC
    NPP_CTRL,
    NEP_CTRL,
    NEE_CTRL,
    GPP_CTRL,
    MR_CTRL,
    GR_CTRL,
    HR_CTRL,
    FIRE_CTRL,
    LITFALLC_CTRL,
    VEGC_CTRL,
    AGC_CTRL,
    LITRC_CTRL,
    SOILC_CTRL,
    TOTALC_CTRL,
    SMINN_CTRL,
    SURFTEC_CTRL,
    UNSATTEC_CTRL,
    GWTEC_CTRL,
    RIVSTGTEC_CTRL,
    RIVGWTEC_CTRL,
};

#if defined(_CYCLES_)
#define MAXOP                   100         // maximum number of operations
#define MAXCROP                 100         // maximum number of crops

// Dimensions in solute arrays
#define NO3                     0
#define NH4                     1

#define NUM_MA_DAYS             7           // length of moving window for soil temperature average used by conditional
                                            // planting
#define COND_DEPTH              0.075       // depth at which soil temperature and moisture are checked planting

// Clipping biomass destiny
#define REMOVE_CLIPPING         0           // clipped biomass harvested
#define RETURN_CLIPPING         1           // clipped biomass returned to soil surface
#define GRAZING_CLIPPING        2           // clipped biomass consumed by livestock

// Operation types
#define OPER_TYPES              4           // number of operation types
#define PLANT_OP                0           // planting
#define TILLAGE_OP              1           // tillage
#define FIXIRR_OP               2           // scheduled irrigation
#define FIXFERT_OP              3           // scheduled fertilization

// Tillage types
#define TILLAGE                 0           // tillage
#define GRAIN_HARVEST           1           // grain harvest
#define FORAGE_HARVEST          2           // forage harvest
#define KILL_CROP               3           // kill crops
#define BURN_RESIDUE            4           // burn residue
#define ALL_CROPS               -1          // flag to kill all crops

// Crop life cycle types
#define ANNUAL_CROP             1           // annual crop
#define PERENNIAL_CROP          0           // perennial crop

// Crop photosynthesis pathways
#define C3_CROP                 1           // C3 crop
#define C4_CROP                 0           // C4 crop

// Crop growth stages
#define NOT_USED                -999        // not being used in the simulation
#define KILLED                  -1          // killed
#define NO_CROP                 0           // not being planted
#define PRE_EMERGENCE           1           // pre-emergence
#define VEGETATIVE_GROWTH       2           // vegetative growth
#define PERENNIAL               3           // perennial
#define REPRODUCTIVE_GROWTH     4           // reproductive growth
#define MATURITY                5           // maturity
#define CLIPPING                6           // clipping
#define PLANTING                7           // planting

#define SPINUP_TOLERANCE        0.01        // spin-up tolerance for change in soil organic carbon (Mg ha-1 year -1)

#define STAN_RESIDUE_SA         4.0         // standing residue area to mass ratio (m2/kg)
#define FLAT_RESIDUE_SA         4.0         // flat residue area to mass ratio (m2/kg)
#define STAN_RESIDUE_K          0.25        // standing residue extinction coefficient
#define FLAT_RESIDUE_K          1.0         // flat residue extinction

#define MAX_SOC_DECOMP_RATE     1.5E-4      // maximum soil organic carbon decomposition rate (day-1)
#define MAX_RESIDUE_DECOMP_RATE 0.05        // maximum residue decomposition rate (day-1)
#define MAX_ROOT_DECOMP_RATE    0.05        // maximum root decomposition rate (day-1)
#define MAX_RHIZO_DECOMP_RATE   0.1         // maximum rhizome decomposition rate (day-1)
#define MAX_MANURE_DECOMP_RATE  0.05        // maximum manure decomposition rate (day-1)
#define MAX_MICROB_DECOMP_RATE  1.0         // maximum microbe decomposition rate (calculated internally) (day-1)
#define C_FRAC_PLANT            0.43        // C fraction in plant
#define C_FRAC_RHIZO            0.43        // C fraction in rhizome
#define C_FRAC_MANURE           0.4         // C fraction in manure
#define SOC_HUMIF_POWER         6.0         // soil organic carbon humification exponent

#define WATER_DENSITY           1000.0      // water density (kg m-3)

#define CO2_DEFAULT             400.0       // default atmospheric CO2 concentration (ppm)

#define KD_NO3                  0.0         // adsorption coefficient for NO3 (cm3 g-1)
#define KD_NH4                  5.6         // adsorption coefficient for NH4 (cm3 g-1)

#define NITRIF_CONST            0.2         // nitrification rate (day-1)
#define POT_DENITRIF            3.2E-5      // potential denitrification rate (kg N kg-1 soil day-1)
#define DENITRIF_HALF_RATE      6.0E-5      // half saturation constant for denitrification (kg N kg-1 soil)
#define DECOMP_HALF_RESP        0.22        // decomposition half response to saturation
#define DECOMP_RESP_POWER       3.0         // decomposition exponential response to saturation

#define ALPHA_C3                0.0033      // RUE and TUE adjustment factor for C3 crops
#define BETA_C3                 1.28        // RUE and TUE adjustment factor for C3 crops
#define THETA_C3                0.94        // RUE and TUE adjustment factor for C3 crops
#define ALPHA_C4                0.006       // RUE and TUE adjustment factor for C4 crops
#define BETA_C4                 1.05        // RUE and TUE adjustment factor for C4 crops
#define THETA_C4                0.95        // RUE and TUE adjustment factor for C4 crops

#define TB                      0           // normalized time at beginning of growth period
#define TM                      0.52        // normalized time at maximum growth rate
#define TE                      1           // normalized time at end of growth period
#endif

// Both macro NSOLUTE and global variable nsolute are needed. NSOLUTE is used for declare a large enough array size and
// nsolute is for loops in the code
#if defined(_BGC_)
# define NSOLUTE                1
#elif defined(_CYCLES_)
# define NSOLUTE                2
#endif

// External variable
extern int     verbose_mode;
extern int     debug_mode;
extern int     append_mode;
extern int     corr_mode;
extern int     spinup_mode;
extern int     fixed_length;
extern char    project[MAXSTRING];
extern int     nelem;
extern int     nriver;
#if defined(_BGC_)
extern int     first_balance;
#endif
#if defined(_OPENMP)
extern int     nthreads;
#endif
#if defined(_TRANSPORT_)
int            nsolute;
#endif

#endif
