#ifndef CYCLES_CONST_HEADER
#define CYCLES_CONST_HEADER

#define PI              3.14159265358979
#define MAXSTRING       1024
#define BADVAL          -999

#define REMOVE_CLIPPING     0
#define RETURN_CLIPPING     1
#define GRAZING_CLIPPING    2

#define STAN_RESIDUE_SA 4.0     /* Standing residue area to mass ratio
                                 * (m2/kg) */
#define FLAT_RESIDUE_SA 4.0     /* Flat residue area to mass ratio (m2/kg) */
#define STAN_RESIDUE_K  0.25    /* Standing residue extinction coefficient */
#define FLAT_RESIDUE_K  1.0     /* flat residue extinction */

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
#define POTENTIAL_DENITRIFICATION           0.000032    /* kg N / kg soil / day */
#define DENITRIFICATION_HALF_RATE           0.00006 /* kg N / kg Soil */
#define NITRIFICATION_NO3_NH4_RATIO         8   /* NO3-N / NH4-N */

enum stage
{ NO_CROP, PRE_EMERGENCE, VEGETATIVE_GROWTH, PERENNIAL, REPRODUCTIVE_GROWTH, MATURITY, CLIPPING, PLANTING };

extern int      verbose_mode;
extern int      debug_mode;
#endif
