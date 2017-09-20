#ifndef PIHM_CONST_HEADER
#define PIHM_CONST_HEADER

/* Physical parameters */
#define PI          3.14159265
#define DAYINSEC    86400       /* number of seconds in a day */
#define GRAV        9.80665     /* gravity constant [m s-2] */
#define CP          1004.0      /* specific heat capacity of air
                                 * [J kg-1 K-1] */
#define LVH2O       2.501e6     /* latent heat of vaporization [J kg-1] */
#define SIGMA       5.67e-8     /* stefan-Boltzmann constant [W m-2 K-4] */
#define RD          287.04      /* gas constant for dry air [J kg-1 K-1] */
#define RV          461.5       /* gas constant for water vapor [J kg-1 K-1] */
#define CPH2O	    4.218e3     /* specific heat capacity of water
                                 * [J kg-1 K-1] */
#define CPICE	    2.106e3     /* specific heat capacity of ice
                                 * [J kg-1 K-1] */
#define LSUBF	    3.335e5     /* latent heat of fusion [J kg-1] */
#define EMISSI_S    0.95        /* emissivity of snow [-] */
#define TFREEZ      273.15      /* freezing point [K] */
#define LSUBS       2.83e6      /* latent heat of sublimation [J kg-1] */

#define F_OK            0

/* Simulation mode */
#define NORMAL_MODE     0
#define SPINUP_MODE     1
#define ACC_SPINUP_MODE 2

/* Default bad value */
#define BADVAL	    -999

/* Verbosity level */
#define VL_ERROR    -1
#define VL_NORMAL   0
#define VL_VERBOSE  1

/* Model steps */
#define HYDROL_STEP 0
#define LS_STEP     1
#define CN_STEP     2

/* Maximum string length */
#define MAXSTRING   1024

/* Maximum number of output files */
#define MAXPRINT    1024

/* Meteorological forcing related */
#define NUM_METEO_VAR   7       /* number of meteo forcing variables */
#define PRCP_TS         0       /* index of precipitation forcing */
#define SFCTMP_TS       1       /* index of air temperature forcing */
#define RH_TS           2       /* index of RH forcing */
#define SFCSPD_TS       3       /* index of wind speed forcing */
#define SOLAR_TS        4       /* index of solar radiation forcing */
#define LONGWAVE_TS     5       /* index of longwave radiation forcing */
#define PRES_TS         6       /* index of surface pressure forcing */

/* Radiation forcing variables */
#define UNIF_SOL    0
#define TOPO_SOL    1

#define SOLDIR_TS   0           /* index of direct solar radiation forcing */
#define SOLDIF_TS   1           /* index of diffused solar radiation forcing */

/* Number of edges of an element */
#define NUM_EDGE    3

#define RIVER_WDTH  0
#define RIVER_AREA  1
#define RIVER_PERIM 2

/* Number of river fluxes of a river segment */
#define NUM_RIVFLX  11

/* Hydrology parameters */
#define PSIMIN	    -70.0       /* minimum psi allowed [m] */
#define DEPRSTG     1e-4        /* depressio storage [m] */
#define GRADMIN     5e-8        /* minimum hydraulic gradient [m m-1] */
#define SATMIN      0.1         /* minimum saturation ratio [-] */
#define RIVDPTHMIN  0.05        /* minimum river depth [m] */
#define RIVGRADMIN  0.05        /* minimum river hydarulic gradient [m m-1] */
#define CMCFACTR    2E-4        /* canopy water capacity per LAI [m] */

/* Maximum of soil layers in Flux-PIHM */
#define MAXLYR      11

/* Land cover parameters */
#define NLCTYPE     40          /* number of land cover types */
#define ISURBAN     13          /* land cover type representing urban */

/* Land cover types */
#define ENF             1
#define EBF             2
#define DNF             3
#define DBF             4
#define MIXF            5
#define CLOSE_SHRUB     6
#define OPEN_SHRUB      7
#define WOODY_SAVANNA   8
#define SAVANNA         9
#define GRASS           10
#define PWL             11
#define CROP            12
#define URBAN_BUILDUP   13
#define CROP_NATURAL    14
#define SNOW_ICE        15
#define BARREN          16
#define WATER           17
#define WOOD_TUNDRA     18
#define MIX_TUNDRA      19
#define BARREN_TUNDRA   20

/* Soil textures */
#define SAND            0
#define LOAMY_SAND      1
#define SANDY_LOAM      2
#define LOAM            3
#define SILT_LOAM       4
#define SILT            5
#define SANDY_CLAY_LOAM 6
#define CLAY_LOAM       7
#define SILTY_CLAY_LOAM 8
#define SANDY_CLAY      9
#define SILTY_CLAY      10
#define CLAY            11

/* Macropore status */
#define MTX_CTRL    0           /* matrix control */
#define APP_CTRL    1           /* application control */
#define MAC_CTRL    2           /* macropore control */

/* River fluxes */
#define UP_CHANL2CHANL      0
#define DOWN_CHANL2CHANL    1
#define LEFT_SURF2CHANL     2
#define RIGHT_SURF2CHANL    3
#define LEFT_AQUIF2CHANL    4
#define RIGHT_AQUIF2CHANL   5
#define CHANL_LKG           6
#define LEFT_AQUIF2AQUIF    7
#define RIGHT_AQUIF2AQUIF   8
#define DOWN_AQUIF2AQUIF    9
#define UP_AQUIF2AQUIF      10



/* River segment interpolation order */
#define RECTANGLE           1
#define TRIANGLE            2
#define QUADRATIC           3
#define CUBIC               4

/* Approximation */
#define KINEMATIC           1
#define DIFF_WAVE           2

/* Initialization type */
#define RELAX               0
#define RST_FILE            1

/* Average flux */
#define SUM                 0
#define AVG                 1

/* Ecosystem constants */
#define RAD2PAR             0.45    /* ratio PAR / SWtotal [-] */
#define EPAR                4.55    /* (umol/J) PAR photon energy ratio */
#define SOIL1_CN            12.0    /* C:N for fast microbial recycling pool */
#define SOIL2_CN            12.0    /* C:N for slow microbial recycling pool */
#define SOIL3_CN            10.0    /* C:N for recalcitrant SOM pool (humus) */
#define SOIL4_CN            10.0    /* C:N for recalcitrant SOM pool (humus) */
#define GRPERC              0.3 /* growth resp per unit of C grown [-] */
#define GRPNOW              1.0 /* proportion of storage growth resp at
                                 * fixation [-] */
#define PPFD50              75.0    /* PPFD for 1/2 stomatal closure
                                     * [umol m-2 s-1] */
#define DENITRIF_PROPORTION 0.01    /* fraction of mineralization to
                                     * volatile */
#define MOBILEN_PROPORTION  0.1 /* fraction mineral N avail for leaching */

/* Respiration fractions for fluxes between compartments [-] */
#define	RFL1S1	    0.39        /* transfer from litter 1 to soil 1 */
#define	RFL2S2	    0.55        /* transfer from litter 2 to soil 2 */
#define	RFL4S3	    0.29        /* transfer from litter 4 to soil 3 */
#define	RFS1S2	    0.28        /* transfer from soil 1 to soil 2 */
#define	RFS2S3	    0.46        /* transfer from soil 2 to soil 3 */
#define	RFS3S4	    0.55        /* transfer from soil 3 to soil 4 */

/* Base decomposition rate constants [day -1] */
#define KL1_BASE    0.7         /* labile litter pool */
#define KL2_BASE    0.07        /* cellulose litter pool */
#define KL4_BASE    0.014       /* lignin litter pool */
#define KS1_BASE    0.07        /* fast microbial recycling pool */
#define KS2_BASE    0.014       /* medium microbial recycling pool */
#define KS3_BASE    0.0014      /* slow microbial recycling pool */
#define KS4_BASE    0.0001      /* recalcitrant SOM (humus) pool */
#define KFRAG_BASE  0.001       /* physical fragmentation of coarse woody
                                 * debris */

/* Decomposition acceleration terms */
#define KS1_ACC     1.0
#define KS2_ACC     1.0
#define KS3_ACC     5.0
#define KS4_ACC     70.0

/* This constant determines the lower limit of state variables before they are
 * set to 0.0 to control rounding and overflow errors */
#define CRIT_PREC   1e-20

/* This constant is used in if conditions where floating point values are
 * compared */
#define FLT_COND_TOL    1e-10

/* Maximum allowable trend in slow soil carbon at steady-state
 * [kgC m-2 yr-1] */
#define SPINUP_TOLERANCE    0.0005

/* Allocation parameters */
#define DAYSNDEPLOY                 365.0
#define DAYSCRECOVER                365.0
#define BULK_DENITRIF_PROPORTION    0.5

/* Output variables */
#define YEARLY_OUTPUT           -1
#define MONTHLY_OUTPUT          -2
#define DAILY_OUTPUT            -3
#define HOURLY_OUTPUT           -4

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
#define T1_CTRL                 25
#define STC_CTRL                26
#define SMC_CTRL                27
#define SH2O_CTRL               28
#define SNOWH_CTRL              29
#define ALBEDO_CTRL             30
#define LE_CTRL                 31
#define SH_CTRL                 32
#define G_CTRL                  33
#define ETP_CTRL                34
#define ESNOW_CTRL              35
#define ROOTW_CTRL              36
#define SOILM_CTRL              37
#define SOLAR_CTRL              38
#define CH_CTRL                 39
#define BIOMASS_CTRL            40
#define RADNINTCP_CTRL          41
#define WATER_STS_CTRL          42
#define N_STS_CTRL              43
#define CROP_TR_CTRL            44
#define CROP_POTTR_CTRL         45
#define RES_EVAP_CTRL           46
#define NO3_PROF_CTRL           47
#define NO3_RIVER_CTRL          48
#define NH4_PROF_CTRL           49
#define NH4_RIVER_CTRL          50
#define NO3_DENIT_CTRL          51
#define NO3_LEACH_CTRL          52
#define NH4_LEACH_CTRL          53
#define NO3_LEACH_RIVER_CTRL    54
#define NH4_LEACH_RIVER_CTRL    55
#define N_LEACH_CTRL            56
#define LAI_CTRL                57
#define VEGC_CTRL               58
#define LITRC_CTRL              59
#define SOILC_CTRL              60
#define TOTALC_CTRL             61
#define NPP_CTRL                62
#define NEP_CTRL                63
#define NEE_CTRL                64
#define GPP_CTRL                65
#define SMINN_CTRL              66
#define LEAFC_CTRL              67
#define LIVESTEMC_CTRL          68
#define DEADSTEMC_CTRL          69
#define SURFTEC_CTRL            70
#define UNSATTEC_CTRL           71
#define GWTEC_CTRL              72
#define RIVSTGTEC_CTRL          73
#define RIVGWTEC_CTRL           74
#define IC_CTRL					75
#define WB_CTRL					76
#ifdef _CYCLES_
#define MAXOP               100

#define PLANT_OP    0
#define TILLAGE_OP  1
#define FIXIRR_OP   2
#define FIXFERT_OP  3

#define REMOVE_CLIPPING     0
#define RETURN_CLIPPING     1
#define GRAZING_CLIPPING    2

#define STAN_RESIDUE_SA     4.0 /* Standing residue area to mass ratio
                                 * (m2/kg) */
#define FLAT_RESIDUE_SA     4.0 /* Flat residue area to mass ratio (m2/kg) */
#define STAN_RESIDUE_K      0.25    /* Standing residue extinction coefficient */
#define FLAT_RESIDUE_K      1.0 /* flat residue extinction */

#define MAXIMUM_UNDISTURBED_SOC_DECOMPOSITION_RATE  0.00015 /* (1 + 0.056) ^ (1 / 365) - 1  ' 1/day (1/5 for Urbana) */
#define MAXIMUM_RESIDUE_DECOMPOSITION_RATE          0.05    /* 1/day */
#define MAXIMUM_ROOT_DECOMPOSITION_RATE             0.05    /* 1/day */
#define MAXIMUM_RHIZO_DECOMPOSITION_RATE            0.1 /*  1/day */
#define MAXIMUM_MANURE_DECOMPOSITION_RATE           0.05    /* 1/day */
#define MAXIMUM_MICROBIAL_DECOMPOSITION_RATE        1.0 /* calculated internaly
                                                         * 1/day */
#define FRACTION_CARBON_PLANT               0.43
#define FRACTION_CARBON_RIZHO               0.43
#define FRACTION_CARBON_MANURE              0.4

#define SOC_DECOMPOSITION_POWER             0.5
#define SOC_HUMIFICATION_POWER              6.0

#define WATER_DENSITY                       1000.0  /* kg/m3 */

#define THRESHOLD_TEMPERATURE_SNOWFALL      1   /* degree C */
#define THRESHOLD_TEMPERATURE_SNOWMELT      -1  /* degree C */
#define SNOWMELT_RATE                       2.5 /* mm/(C day) or degree day
                                                 * melting factor */

#define NITRIFICATION_CONSTANT              0.2 /* 1/day */
#define POTENTIAL_DENITRIFICATION           0.000032  /* kg N / kg soil / day */
#define DENITRIFICATION_HALF_RATE           0.00006 /* kg N / kg Soil */
#define NITRIFICATION_NO3_NH4_RATIO         8   /* NO3-N / NH4-N */

enum stage
{ NO_CROP, PRE_EMERGENCE, VEGETATIVE_GROWTH, PERENNIAL, REPRODUCTIVE_GROWTH,
    MATURITY, CLIPPING, PLANTING
};

#endif

/* External variable */
extern int      verbose_mode;
extern int      debug_mode;
extern int      corr_mode;
extern int      spinup_mode;
extern int      tecplot;
extern char     project[MAXSTRING];
extern int      nelem;
extern int      nriver;
#ifdef _OPENMP
extern int      nthreads;
#endif
#if defined(_BGC_) || defined (_CYCLES_)
extern int      first_balance;
#endif

#endif
