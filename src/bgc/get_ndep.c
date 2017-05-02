
/*
 * get_ndep.c
 * retrieve the appropriate nitrogen deposition value for the current simulation year
 *
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 * Biome-BGC version 4.2 (final release)
 * See copyright.txt for Copyright information
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 */

#include "pihm.h"

double GetNdep (tsdata_struct ndep_ts, int t)
{
    double          ndep;

    IntrplForcing (&ndep_ts, t, 1);

    ndep = ndep_ts.value[0];

    return (ndep);
}
