#include "pihm.h"

double GetCO2(tsdata_struct *co2_ts, int t)
{
    double          co2;

    IntrplForc(co2_ts, t, 0, 1, INTRPL);

    co2 = co2_ts->value[0];

    return co2;
}
