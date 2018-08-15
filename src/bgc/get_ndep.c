#include "pihm.h"

double GetNdep(tsdata_struct *ndep_ts, int t)
{
    double          ndep;

    IntrplForc(t, 1, ndep_ts);

    ndep = ndep_ts->value[0];

    return ndep;
}
