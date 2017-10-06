#include "pihm.h"

double GetNdep(tsdata_struct *ndep_ts, int t)
{
    double          ndep;

    IntrplForc(ndep_ts, t, 1);

    ndep = ndep_ts->value[0];

    return ndep;
}
