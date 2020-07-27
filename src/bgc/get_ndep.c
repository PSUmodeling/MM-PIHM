#include "pihm.h"

double GetNdep(tsdata_struct *ndep_ts, int t)
{
    double          ndep;

    IntrplForcing(t, 1, INTRPL, ndep_ts);

    ndep = ndep_ts->value[0];

    return ndep;
}
