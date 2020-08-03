#include "pihm.h"

double GetNdep(int t, tsdata_struct *ndep_ts)
{
    double          ndep;

    IntrplForcing(t, 1, INTRPL, ndep_ts);

    ndep = ndep_ts->value[0];

    return ndep;
}
