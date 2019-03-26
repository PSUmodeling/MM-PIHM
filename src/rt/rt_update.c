#include "pihm.h"

void RTUpdate(const rttbl_struct *rttbl, elem_struct elem[],
    river_struct river[])
{
    int             i;
    double          riv;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             k;

        for (k = 0; k < rttbl->NumStc; k++)
        {
            elem[i].chms_unsat.log10_pconc[k] =
                log10(elem[i].chms_unsat.p_conc[k]);

            elem[i].chms_gw.log10_pconc[k] =
                log10(elem[i].chms_gw.p_conc[k]);
        }
        for (k = 0; k < rttbl->NumSsc; k++)
        {
            elem[i].chms_unsat.log10_sconc[k] =
                log10(elem[i].chms_unsat.s_conc[k]);

            elem[i].chms_gw.log10_sconc[k] =
                log10(elem[i].chms_gw.s_conc[k]);
        }
    }

    /* Flux for RIVER flow */
    for (i = 0; i < nriver; i++)
    {
        int             k;

        for (k = 0; k < rttbl->NumStc; k++)
        {
            river[i].chms_stream.log10_pconc[k] =
                log10(river[i].chms_stream.p_conc[k]);
        }
        for (k = 0; k < rttbl->NumSsc; k++)
        {
            river[i].chms_stream.log10_sconc[k] =
                log10(river[i].chms_stream.s_conc[k]);
        }

        if (river[i].down < 0)
        {
            riv = river[i].wf.rivflow[1];
        }
    }

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
                        log10((river[-rttbl->BTC_loc[i] - 1].chms_stream.p_conc[k] * riv +
                        rttbl->pumps[0].Injection_conc * 1.0E-3 * rttbl->pumps[0].flow_rate) /
                        (riv + rttbl->pumps[0].flow_rate));
                }
                else
                {
                    river[-rttbl->BTC_loc[i] - 1].chms_stream.btcv_pconc[k] =
                        log10(river[-rttbl->BTC_loc[i] - 1].chms_stream.p_conc[k]);
                }
            }
            else if (rttbl->BTC_loc[i] < nelem)
            {
                elem[rttbl->BTC_loc[i] - 1].chms_unsat.btcv_pconc[k] =
                    elem[rttbl->BTC_loc[i] - 1].chms_unsat.log10_pconc[k];
            }
            else
            {
                elem[rttbl->BTC_loc[i] - nelem - 1].chms_gw.btcv_pconc[k] =
                    elem[rttbl->BTC_loc[i] - 1].chms_gw.log10_pconc[k];
            }
        }
    }
}
