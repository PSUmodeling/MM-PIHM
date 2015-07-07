/* 
 * get_co2.c
 * retrieve the appropriate CO2 concentration for the current simulation year
 * 
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 * Biome-BGC version 4.2 (final release)
 * See copyright.txt for Copyright information
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 */

#include "bgc.h"

double GetCO2 (ts_struct co2_ts, int t)
{
    double          co2;

    IntrplForcing (co2_ts, t, 1, &co2);

    return (co2);
}
