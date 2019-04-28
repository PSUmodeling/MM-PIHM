#include "pihm.h"

void RTUpdate(const rttbl_struct *rttbl, elem_struct elem[],
    river_struct river[])
{
    int             i;
    double          riv;

    for (i = 0; i < rttbl->NumBTC; i++)
    {
        int             k;

        for (k = 0; k < rttbl->NumStc; k++)
        {
            if (rttbl->BTC_loc[i] < 0)
            {
                if (-rttbl->BTC_loc[i] >= -rttbl->pumps[0].Pump_Location &&
                    k == rttbl->pumps[0].Position_Species)
                {
                    river[-rttbl->BTC_loc[i] - 1].chms_stream.btcv_pconc[k] =
                        (river[-rttbl->BTC_loc[i] - 1].chms_stream.p_conc[k] * riv +
                        rttbl->pumps[0].Injection_conc * 1.0E-3 * rttbl->pumps[0].flow_rate) /
                        (riv + rttbl->pumps[0].flow_rate);
                }
                else
                {
                    river[-rttbl->BTC_loc[i] - 1].chms_stream.btcv_pconc[k] =
                        river[-rttbl->BTC_loc[i] - 1].chms_stream.p_conc[k];
                }
            }
            else if (rttbl->BTC_loc[i] < nelem)
            {
                elem[rttbl->BTC_loc[i] - 1].chms_unsat.btcv_pconc[k] =
                    elem[rttbl->BTC_loc[i] - 1].chms_unsat.p_conc[k];
            }
            else
            {
                elem[rttbl->BTC_loc[i] - nelem - 1].chms_gw.btcv_pconc[k] =
                    elem[rttbl->BTC_loc[i] - 1].chms_gw.p_conc[k];
            }
        }
    }
}
