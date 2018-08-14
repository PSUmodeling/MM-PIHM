#include "pihm.h"

void UpdateNProf(int nsoil, double dt, const nflux_struct *nf,
    nstate_struct *ns)
{
    int             k;

    ns->no3[0] += nf->surplusn;
    ns->nh4[0] += nf->urine;

    for (k = 0; k < nsoil; k++)
    {
        ns->no3[k] += (nf->uptake_no3[k] + nf->fert_no3[k] + nf->immob_no3[k] +
            nf->nitrif_nh4_to_no3[k] + nf->denitn[k] + nf->till_no3[k]) *
            (dt / DAYINSEC);

        ns->nh4[k] += (nf->uptake_nh4[k] + nf->fert_nh4[k] + nf->till_nh4[k] +
            nf->immob_nh4[k] - nf->nitrif_nh4_to_no3[k] -
            nf->nitrif_nh4_to_n2o[k] - nf->nh4volat[k]) * (dt / DAYINSEC);
    }
}
