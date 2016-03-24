#ifndef PIHM_CONST_HEADER
#define PIHM_CONST_HEADER

/* Physical parameters */
#define PI		3.14159265
#define GRAV		9.80665
#define CP              1004.0  /* Specific heat capacity of air (J/kg/K) */
#define LVH2O           2.501e6 /* Latent heat of vaporization (J/kg) */
#define SIGMA           5.67e-8 /* Stefan-Boltzmann constant (W/m2/K4) */
#define RD              287.04  /* Gas constant for dry air (J/kg/K) */
#define RV              461.5   /* Gas constant for water vapor (J/kg/K) */
#define CPH2O	        4.218e3 /* Specific heat capacity of water (J/kg/K) */
#define CPICE	        2.106e3
#define LSUBF	        3.335e5
#define EMISSI_S        0.95
#define TFREEZ          273.15
#define LSUBS           2.83e6

#define BADVAL		-999
#define MAXSTRING	1024
#define ISURBAN         13

#define MAXLYR          11
#define NUM_METEO_VAR   7
#define NUM_PRINT       200

#define PSIMIN		-70.0
#define DEPRSTG         1E-4
#define GRADMIN         5E-8
#define SATMIN          0.1
#define RIVDPTHMIN      0.05
#define RIVGRADMIN      0.05

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
#define MTX_CTRL        0
#define APP_CTRL        1
#define MAC_CTRL        2

/* River segment interpolation order */
#define RECTANGLE       1
#define TRIANGLE        2
#define QUADRATIC       3
#define CUBIC           4

/* Approximation */
#define KINEMATIC       1
#define DIFF_WAVE       2

/* Initialization type */
#define RELAX           0
#define RST_FILE        3

/* Average flux */
#define SUM             0
#define AVG             1

/* Output variables */
#define SURF_CTRL       0
#define UNSAT_CTRL      1
#define GW_CTRL         2
#define RIVSTG_CTRL     3
#define RIVGW_CTRL      4
#define SNOW_CTRL       5
#define CMC_CTRL        6
#define INFIL_CTRL      7
#define RECHARGE_CTRL   8
#define EC_CTRL         9
#define ETT_CTRL        10
#define EDIR_CTRL       11
#define RIVFLX0_CTRL    12
#define RIVFLX1_CTRL    13
#define RIVFLX2_CTRL    14
#define RIVFLX3_CTRL    15
#define RIVFLX4_CTRL    16
#define RIVFLX5_CTRL    17
#define RIVFLX6_CTRL    18
#define RIVFLX7_CTRL    19
#define RIVFLX8_CTRL    20
#define RIVFLX9_CTRL    21
#define RIVFLX10_CTRL   22
#define SUBFLX_CTRL     23
#define SURFFLX_CTRL    24

#ifdef _NOAH_
#define T1_CTRL         25
#define STC_CTRL        26
#define SMC_CTRL        27
#define SH2O_CTRL       28
#define SNOWH_CTRL      29
#define ALBEDO_CTRL     30
#define LE_CTRL         31
#define SH_CTRL         32
#define G_CTRL          33
#define ETP_CTRL        34
#define ESNOW_CTRL      35
#define ROOTW_CTRL      36
#define SOILM_CTRL      37
#define SOLAR_CTRL      38
#endif

extern int      verbose_mode;
extern int      debug_mode;
extern char     project[MAXSTRING];

/* Enumrate type for forcing time series */
enum meteo_forcing_type
{ PRCP_TS, SFCTMP_TS, RH_TS, SFCSPD_TS, SOLAR_TS, LONGWAVE_TS, PRES_TS };

enum rad_forcing_type
{ SOLAR_DIR_TS, SOLAR_DIF_TS };

#ifdef _CYCLES_
#define REMOVE_CLIPPING     0
#define RETURN_CLIPPING     1
#define GRAZING_CLIPPING    2

#define STAN_RESIDUE_SA 4.0     /* Standing residue area to mass ratio
                                 * (m2/kg) */
#define FLAT_RESIDUE_SA 4.0     /* Flat residue area to mass ratio (m2/kg) */
#define STAN_RESIDUE_K  0.25    /* Standing residue extinction coefficient */
#define FLAT_RESIDUE_K  1.0     /* flat residue extinction */

#define MAXIMUM_UNDISTURBED_SOC_DECOMPOSITION_RATE  0.00015     /* (1 + 0.056) ^ (1 / 365) - 1  ' 1/day (1/5 for Urbana) */
#define MAXIMUM_RESIDUE_DECOMPOSITION_RATE          0.05        /* 1/day */
#define MAXIMUM_ROOT_DECOMPOSITION_RATE             0.05        /* 1/day */
#define MAXIMUM_RHIZO_DECOMPOSITION_RATE            0.1 /*  1/day */
#define MAXIMUM_MANURE_DECOMPOSITION_RATE           0.05        /* 1/day */
#define MAXIMUM_MICROBIAL_DECOMPOSITION_RATE        1.0 /* calculated internaly
                                                         * 1/day */
#define FRACTION_CARBON_PLANT               0.43
#define FRACTION_CARBON_RIZHO               0.43
#define FRACTION_CARBON_MANURE              0.4

#define SOC_DECOMPOSITION_POWER             0.5
#define SOC_HUMIFICATION_POWER              6.0

#define WATER_DENSITY                       1000.0      /* kg/m3 */

#define THRESHOLD_TEMPERATURE_SNOWFALL      1   /* degree C */
#define THRESHOLD_TEMPERATURE_SNOWMELT      -1  /* degree C */
#define SNOWMELT_RATE                       2.5 /* mm/(C day) or degree day
                                                 * melting factor */

#define NITRIFICATION_CONSTANT              0.2 /* 1/day */
#define POTENTIAL_DENITRIFICATION           0.000032    /* kg N / kg soil / day */
#define DENITRIFICATION_HALF_RATE           0.00006     /* kg N / kg Soil */
#define NITRIFICATION_NO3_NH4_RATIO         8   /* NO3-N / NH4-N */

enum stage
{ NO_CROP, PRE_EMERGENCE, VEGETATIVE_GROWTH, PERENNIAL, REPRODUCTIVE_GROWTH,
        MATURITY, CLIPPING, PLANTING };

#endif
#ifdef _ENKF_
#define MAXPARAM        100
#define MAXVAR          100
#define MAXINT          2147483647
#define CORRMAX         0.25
#define SUCCESS_TAG     2
#define CYCLE_TAG       1
#define PARAM_TAG       3
#define LOG_TYPE        1

enum prmt_type
{ KSATH, KSATV, KINF, KMACH, KMACV, DINF, RZD, DMAC, POROSITY,
    ALPHA, BETA, AREAFV, AREAFH, VEGFRAC, ALBEDO, ROUGH, PRCP, SFCTMP,
    EC, ETT, EDIR, RIVROUGH, RIVKSATH, RIVKSATV, RIVBEDTHICK, RIVDEPTH,
    RIVSHPCOEFF, DRIP, INTCP, RSMIN, CZIL, FXEXP, CFACTR, RGL, HS, THETAREF,
    THETAW
};

enum obs_type
{ RUNOFF_OBS, TSKIN_OBS };

#endif
#endif
