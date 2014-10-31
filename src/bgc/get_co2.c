/* 
get_co2.c
retrieve the appropriate CO2 concentration for the current simulation year

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

double get_co2(TSD co2_ts, double t)
{
    return Interpolation (&co2_ts, t);
}
