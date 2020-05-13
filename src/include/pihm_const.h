#ifndef PIHM_CONST_HEADER
#define PIHM_CONST_HEADER

/* Physical parameters */
#define PI                      3.14159265  /* pi */
#define DAYINSEC                86400       /* number of seconds in a day */
#define GRAV                    9.80665     /* gravity constant (m s-2) */
#define CP                      1004.0      /* specific heat capacity of air
                                             * (J kg-1 K-1) */
#define RHOH2O                  1000.0      /* water density (kg m-3) */
#define LVH2O                   2.501E6     /* latent heat of vaporization
                                             * (J kg-1) */
#define SIGMA                   5.67E-8     /* Stefan-Boltzmann constant
                                             * (W m-2 K-4) */
#define RD                      287.04      /* gas constant for dry air
                                             * (J kg-1 K-1) */
#define RV                      461.5       /* gas constant for water vapor
                                             * (J kg-1 K-1) */
#define CPH2O                   4.218E3     /* specific heat capacity of water
                                             * (J kg-1 K-1) */
#define CPICE                   2.106E3     /* specific heat capacity of ice
                                             * (J kg-1 K-1) */
#define LSUBF                   3.335E5     /* latent heat of fusion (J kg-1) */
#define EMISSI_S                0.95        /* emissivity of snow (-) */
#define TFREEZ                  273.15      /* freezing point (K) */
#define LSUBS                   2.83E6      /* latent heat of sublimation
                                             * (J kg-1) */
#define ICEDENS                 0.9         /* glacier ice density (g cm-3, or
                                             * dimensionless fraction of H2O
                                             * density */

#define F_OK                    0

#if defined(_LUMPED_)
# define LUMPED    nelem
#endif

/* Simulation mode */
#define NORMAL_MODE             0
#define SPINUP_MODE             1
#define ACC_SPINUP_MODE         2

/* Default bad value */
#define BADVAL                  -999

/* Model steps */
#define HYDROL_STEP             0
#define LS_STEP                 1
#define CN_STEP                 2
#define RT_STEP                 3

/* Maximum number of output files */
#define MAXPRINT                1024

/* Meteorological forcing related */
#define NUM_METEO_VAR           7           /* number of meteorological forcing
                                             * variables */
#define PRCP_TS                 0           /* index of precipitation forcing */
#define SFCTMP_TS               1           /* index of air temperature forcing
                                             */
#define RH_TS                   2           /* index of RH forcing */
#define SFCSPD_TS               3           /* index of wind speed forcing */
#define SOLAR_TS                4           /* index of solar radiation forcing
                                             */
#define LONGWAVE_TS             5           /* index of longwave radiation
                                             * forcing */
#define PRES_TS                 6           /* index of surface pressure forcing
                                             */

/* Radiation forcing variables */
#define UNIF_SOL                0           /* use global solar radiation */
#define TOPO_SOL                1           /* use topographic solar radiation*/
#define SOLDIR_TS               0           /* index of direct solar radiation
                                             * forcing */
#define SOLDIF_TS               1           /* index of diffused solar radiation
                                             * forcing */

/* Number of edges of an element */
#define NUM_EDGE                3

#define WS_ZMAX                 0
#define WS_ZMIN                 1
#define WS_AREA                 2

/* Hydrology parameters */
#define PSIMIN                  -70.0       /* minimum psi allowed (m) */
#define DEPRSTG                 1E-4        /* depression storage (m) */
#define GRADMIN                 5E-8        /* minimum hydraulic gradient
                                             * (m m-1) */
#define SATMIN                  0.1         /* minimum saturation ratio (-) */
#define RIVDPTHMIN              0.05        /* minimum river depth (m) */
#define RIVGRADMIN              0.05        /* minimum river hydraulic gradient
                                             * (m m-1) */
#define CMCFACTR                2E-4        /* canopy water capacity per LAI (m)
                                             */
#define SH2OMIN                 0.02        /* minimum sh2o (m3 m-3) */

/* Maximum of soil layers in Flux-PIHM */
#define MAXLYR                  11

/* Land cover parameters */
#define NLCTYPE                 40          /* number of land cover types */

/* Land cover types */
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


/* Soil textures */
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

/* Number of river fluxes of a river segment */
#if defined(_FBR_) && defined(_TGM_)
# define NUM_RIVFLX             8
#else
# define NUM_RIVFLX             6
#endif

/* River fluxes */
#define UP_CHANL2CHANL          0
#define DOWN_CHANL2CHANL        1
#define LEFT_SURF2CHANL         2
#define RIGHT_SURF2CHANL        3
#define LEFT_AQUIF2CHANL        4
#define RIGHT_AQUIF2CHANL       5
#if defined(_FBR_) && defined(_TGM_)
# define LEFT_FBR2CHANL         6
# define RIGHT_FBR2CHANL        7
#endif

/* River boundary condition types */
#define OUTLET_DIRICHLET        -1
#define OUTLET_NEUMANN          -2
#define ZERO_DPTH_GRAD          -3
#define CRIT_DPTH               -4

/* Element boundary condition types */
#define NO_FLOW                 0
#define DIRICHLET               1
#define NEUMANN                 2

/* River segment interpolation order */
#define RECTANGLE               1
#define TRIANGLE                2
#define QUADRATIC               3
#define CUBIC                   4

/* Approximation */
#define KINEMATIC               1
#define DIFF_WAVE               2

/* Initialization type */
#define RELAX                   0
#define RST_FILE                1

/* Average flux */
#define SUM                     0
#define AVG                     1

/* Interpolate forcing time series */
#define NO_INTRPL               0
#define INTRPL                  1

/* Maximum allowable difference between simulation cycles in subsurface water
 * storage at steady-state (m) */
#define SPINUP_W_TOLERANCE      0.01

/* Ecosystem constants */
#define RAD2PAR                 0.45        /* ratio PAR / SWtotal (-) */
#define EPAR                    4.55        /* (umol/J) PAR photon energy ratio
                                             */
#define SOIL1_CN                12.0        /* C:N for fast microbial recycling
                                             * pool */
#define SOIL2_CN                12.0        /* C:N for slow microbial recycling
                                             * pool */
#define SOIL3_CN                10.0        /* C:N for recalcitrant SOM pool
                                             * (humus) */
#define SOIL4_CN                10.0        /* C:N for recalcitrant SOM pool
                                             * (humus) */
#define GRPERC                  0.3         /* growth resp per unit of C grown
                                             * (-) */
#define GRPNOW                  1.0         /* proportion of storage growth resp
                                             * at fixation (-) */
#define PPFD50                  75.0        /* PPFD for 1/2 stomatal closure
                                             * (umol m-2 s-1) */
#define DENITRIF_PROPORTION     0.01        /* fraction of mineralization to
                                             * volatile */
#if defined(_CYCLES_OBSOLETE_)
# define MOBILEN_PROPORTION     1.0         /* fraction mineral N avail for
                                             * leaching */
#else
# define MOBILEN_PROPORTION     0.1         /* fraction mineral N avail for
                                             * leaching */
#endif

/* Respiration fractions for fluxes between compartments (-) */
#define RFL1S1                  0.39        /* transfer from litter 1 to soil 1
                                             */
#define RFL2S2                  0.55        /* transfer from litter 2 to soil 2
                                             */
#define RFL4S3                  0.29        /* transfer from litter 4 to soil 3
                                             */
#define RFS1S2                  0.28        /* transfer from soil 1 to soil 2 */
#define RFS2S3                  0.46        /* transfer from soil 2 to soil 3 */
#define RFS3S4                  0.55        /* transfer from soil 3 to soil 4 */

/* Base decomposition rate constants (day -1) */
#define KL1_BASE                0.7         /* labile litter pool */
#define KL2_BASE                0.07        /* cellulose litter pool */
#define KL4_BASE                0.014       /* lignin litter pool */
#define KS1_BASE                0.07        /* fast microbial recycling pool */
#define KS2_BASE                0.014       /* medium microbial recycling pool*/
#define KS3_BASE                0.0014      /* slow microbial recycling pool */
#define KS4_BASE                0.0001      /* recalcitrant SOM (humus) pool */
#define KFRAG_BASE              0.001       /* physical fragmentation of coarse
                                             * woody debris */

/* Decomposition acceleration terms */
#define KS1_ACC                 1.0
#define KS2_ACC                 1.0
#define KS3_ACC                 5.0
#define KS4_ACC                 70.0

/* This constant determines the lower limit of state variables before they are
 * set to 0.0 to control rounding and overflow errors */
#define CRIT_PREC               1E-20

/* This constant is used in if conditions where floating point values are
 * compared */
#define FLT_COND_TOL            1E-10

/* Maximum allowable trend in slow soil carbon at steady-state (kgC m-2 yr-1) */
#define SPINUP_C_TOLERANCE      0.0005

/* Allocation parameters */
#define DAYSNDEPLOY             365.0
#define DAYSCRECOVER            365.0
#define BULK_DENITRIF_PROPORTION 0.5

/*
 * RT constants
 */
/* Maximum number of species */
#define MAXSPS                  20

/* Maximum number of dependece, monod, and inhibition terms */
#define MAXDEP                  4

#define ZERO_CONC               1.0E-20

/* Threshold of water storage when deposition occurs */
#define DEPTHR                  1.0E-5

/* RT simulation mode */
#define KIN_REACTION            0
#define TRANSPORT_ONLY          1

/* RT primary species types */
#define AQUEOUS                 1
#define ADSORPTION              2
#define CATION_ECHG             3
#define MINERAL                 4

/* RT mass action types */
#define IMMOBILE_MA             0
#define MOBILE_MA               1
#define MIXED_MA                2

/* RT kinetic reaction types */
#define TST                     1
#define PRCP_ONLY               2
#define DISS_ONLY               3
#define MONOD                   4

/* RT volumes in each model grid */
#if defined(_FBR_)
# define NCHMVOL                2
#else
# define NCHMVOL                1
#endif
#define SOIL_CHMVOL             0
#define GEOL_CHMVOL             1

/* Output variables */
#define YEARLY_OUTPUT           -1
#define MONTHLY_OUTPUT          -2
#define DAILY_OUTPUT            -3
#define HOURLY_OUTPUT           -4

/* Output variable types */
#define ELEMVAR                 0
#define RIVERVAR                1

#define SURF_CTRL               0
#define UNSAT_CTRL              1
#define GW_CTRL                 2
#define RIVSTG_CTRL             3
#define RIVGW_CTRL              4
#define SNOW_CTRL               5
#define CMC_CTRL                6
#define INFIL_CTRL              7
#define RECHARGE_CTRL           8
#define EC_CTRL                 9
#define ETT_CTRL                10
#define EDIR_CTRL               11
#define RIVFLX0_CTRL            12
#define RIVFLX1_CTRL            13
#define RIVFLX2_CTRL            14
#define RIVFLX3_CTRL            15
#define RIVFLX4_CTRL            16
#define RIVFLX5_CTRL            17
#define RIVFLX6_CTRL            18
#define RIVFLX7_CTRL            19
#define RIVFLX8_CTRL            20
#define RIVFLX9_CTRL            21
#define RIVFLX10_CTRL           22
#define SUBFLX_CTRL             23
#define SURFFLX_CTRL            24
/* System */
#define IC_CTRL                 25
#define WB_CTRL                 26
/* Noah */
#define T1_CTRL                 27
#define STC_CTRL                28
#define SMC_CTRL                29
#define SH2O_CTRL               30
#define SNOWH_CTRL              31
#define ALBEDO_CTRL             32
#define LE_CTRL                 33
#define SH_CTRL                 34
#define G_CTRL                  35
#define ETP_CTRL                36
#define ESNOW_CTRL              37
#define ROOTW_CTRL              38
#define SOILM_CTRL              39
#define SOLAR_CTRL              40
#define CH_CTRL                 41
/* BGC or Cycles */
#define LAI_CTRL                42
#define N_PROFILE_CTRL          43
#define N_RIVER_CTRL            44
#define LEACHING_CTRL           45
/* Cycles */
#define BIOMASS_CTRL            46
#define RADNINTCP_CTRL          47
#define WATER_STS_CTRL          48
#define N_STS_CTRL              49
#define CROP_TR_CTRL            50
#define CROP_POTTR_CTRL         51
#define RES_EVAP_CTRL           52
#define DENITRIF_CTRL           53
/* BGC */
#define NPP_CTRL                54
#define NEP_CTRL                55
#define NEE_CTRL                56
#define GPP_CTRL                57
#define MR_CTRL                 58
#define GR_CTRL                 59
#define HR_CTRL                 60
#define FIRE_CTRL               61
#define LITFALLC_CTRL           62
#define VEGC_CTRL               63
#define AGC_CTRL                64
#define LITRC_CTRL              65
#define SOILC_CTRL              66
#define TOTALC_CTRL             67
#define SMINN_CTRL              68
#define SURFTEC_CTRL            69
#define UNSATTEC_CTRL           70
#define GWTEC_CTRL              71
#define RIVSTGTEC_CTRL          72
#define RIVGWTEC_CTRL           73
/* DGW */
#define FBRUNSAT_CTRL           74
#define FBRGW_CTRL              75
#define FBRINFIL_CTRL           76
#define FBRRECHG_CTRL           77
#define FBRFLOW_CTRL            78
/* RT */
#define CHEM_CTRL               79

#if defined(_CYCLES_)
#define MAXOP                   100
#define MAXCROP                 100

#define REMOVE_CLIPPING         0
#define RETURN_CLIPPING         1
#define GRAZING_CLIPPING        2

#define NOT_USED                -999
#define KILLED                  -1
#define NO_CROP                 0
#define PRE_EMERGENCE           1
#define VEGETATIVE_GROWTH       2
#define PERENNIAL               3
#define REPRODUCTIVE_GROWTH     4
#define MATURITY                5
#define CLIPPING                6
#define PLANTING                7

#define OPER_TYPES              4
#define PLANT_OP                0
#define TILLAGE_OP              1
#define FIXIRR_OP               2
#define FIXFERT_OP              3

#define ALL_CROPS                -1

#define FC_BOUND                -1
#define PWP_BOUND               -2

#endif
#if defined(_CYCLES_OBSOLETE_)
#define STAN_RESIDUE_SA         4.0         /* standing residue area to mass
                                             * ratio (m2 kg-1) */
#define FLAT_RESIDUE_SA         4.0         /* flat residue area to mass ratio
                                             * (m2 kg-1) */
#define STAN_RESIDUE_K          0.25        /* standing residue extinction
                                             * coefficient (-) */
#define FLAT_RESIDUE_K          1.0         /* flat residue extinction (-) */

#define MAXIMUM_UNDISTURBED_SOC_DECOMPOSITION_RATE  0.00015 /* (day-1)*/
#define MAXIMUM_RESIDUE_DECOMPOSITION_RATE          0.05    /* (day-1) */
#define MAXIMUM_ROOT_DECOMPOSITION_RATE             0.05    /* (day-1) */
#define MAXIMUM_RHIZO_DECOMPOSITION_RATE            0.1     /* (day-1) */
#define MAXIMUM_MANURE_DECOMPOSITION_RATE           0.05    /* (day-1) */
#define MAXIMUM_MICROBIAL_DECOMPOSITION_RATE        1.0     /* (day-1) */
#define FRACTION_CARBON_PLANT   0.43
#define FRACTION_CARBON_RIZHO   0.43
#define FRACTION_CARBON_MANURE  0.4

#define SOC_DECOMPOSITION_POWER 0.5
#define SOC_HUMIFICATION_POWER  6.0

#define WATER_DENSITY           1000.0      /* (kg m-3) */

#define THRESHOLD_TEMPERATURE_SNOWFALL  1   /* (degree C) */
#define THRESHOLD_TEMPERATURE_SNOWMELT  -1  /* (degree C) */
#define SNOWMELT_RATE                   2.5 /* (mm C-1 day-1 or degree day) */
#define NITRIFICATION_CONSTANT          0.2 /* (1 day-1) */
#define POTENTIAL_DENITRIFICATION       0.000032    /* (kgN kgSoil-1 day-1) */
#define DENITRIFICATION_HALF_RATE       0.00006     /* (kgN kgSoil-1) */
#define NITRIFICATION_NO3_NH4_RATIO     8           /* (NO3 / NH4) */

enum stage
{NO_CROP, PRE_EMERGENCE, VEGETATIVE_GROWTH, PERENNIAL, REPRODUCTIVE_GROWTH,
    MATURITY, CLIPPING, PLANTING
};
#endif

/* Both macro NSOLUTE and global variable nsolute are needed. NSOLUTE is used
 * for declare a large enough array size and nsolute is for loops in the code */
#if defined(_BGC_)
# define NSOLUTE                        1
#elif defined(_CYCLES)
# define NSOLUTE                        2
#elif defined(_RT_)
# define NSOLUTE                        MAXSPS
#endif

/* External variable */
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
#if defined(_BGC_) || defined(_CYCLES_OBSOLETE_) || defined(_RT_)
int             nsolute;
#endif

#endif
