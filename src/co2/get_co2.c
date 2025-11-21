#include "pihm.h"

double GetCO2(int t, tsdata_struct *co2_ts)
{
    double          co2;

    IntrplForcing(t, 1, INTRPL, FORCING_CO2, co2_ts);

    co2 = co2_ts->value[0];

    return co2;
}
