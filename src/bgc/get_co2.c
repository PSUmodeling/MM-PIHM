#include "pihm.h"

double GetCO2(tsdata_struct co2_ts, int t)
{
    double          co2;

    IntrplForcing(&co2_ts, t, 1);

    co2 = co2_ts.value[0];

    return co2;
}
