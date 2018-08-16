#include "pihm.h"

void UpdNProf(double dt, const soil_struct *soil, const wstate_struct *ws,
    const nstate_struct *ns0, const nflux_struct *nf, const nprof_struct *np,
    pstate_struct *ps, nstate_struct *ns)
{
    int             k;

    /*
     * Add source sink terms to NO3 and NH4 mass at different layers
     */
    ns->no3[0] = ns0->no3[0] + nf->surplusn * (dt / DAYINSEC);
    ns->nh4[0] = ns0->nh4[0] + nf->urine * (dt / DAYINSEC);

    for (k = 0; k < ps->nsoil; k++)
    {
        ns->no3[k] = ns0->no3[k] +
            (nf->uptake_no3[k] + nf->fert_no3[k] + nf->immob_no3[k] +
            nf->nitrif_nh4_to_no3[k] + nf->denitn[k] + nf->till_no3[k]) *
            (dt / DAYINSEC);

        ns->nh4[k] = ns0->nh4[k] +
            (nf->uptake_nh4[k] + nf->fert_nh4[k] + nf->till_nh4[k] +
            nf->immob_nh4[k] - nf->nitrif_nh4_to_no3[k] -
            nf->nitrif_nh4_to_n2o[k] - nf->nh4volat[k]) * (dt / DAYINSEC);
    }

    /*
     * Add lateral transport fluxes to NO3 and NH4 mass
     */
    ps->nwtbl = FindWaterTable(ps->sldpth, ps->nsoil, ws->gw, ps->satdpth);

    CalcLatNFlux(0.0, ps->nsoil, ps->sldpth, ps->satdpth, soil->bd, ws->smc,
        ns0->no3, np->no3, ns->no3);
    CalcLatNFlux(5.6, ps->nsoil, ps->sldpth, ps->satdpth, soil->bd, ws->smc,
        ns0->nh4, np->nh4, ns->nh4);
}

void CalcLatNFlux(double kd, int nsoil, const double sldpth[],
    const double satdpth[], const double bd[], const double swc[],
    const double solute0[], double np, double solute[])
{
    int             k;
    double          np0 = 0.0;
    double          totwght = 0.0;
    double          weight[MAXLYR];

    for (k = 0; k < nsoil; k++)
    {
        np0 += solute[k];

        if (satdpth[k] > 0.0)
        {
            weight[k] = (solute[k] > 0.0) ?
                LinearEquilibriumConcentration(kd, bd[k], sldpth[k], swc[k],
                    solute0[k]) * satdpth[k] * swc[k] * RHOH2O : 0.0;
            totwght += weight[k];
        }
        else
        {
            weight[k] = 0.0;
        }
    }

    if (totwght <= 0.0)
    {
        totwght = 1.0;
        weight[nsoil - 1] = 1.0;
    }

    for (k = 0; k < nsoil; k++)
    {
        weight[k] /= totwght;

        solute[k] += weight[k] * (np - np0);
    }

    for (k = nsoil - 1; k > 0; k--)
    {
        if (solute[k] < 0.0)
        {
            solute[k - 1] += solute[k];
            solute[k] = 0.0;
        }
    }

    for (k = 0; k < nsoil - 1; k++)
    {
        if (solute[k] < 0.0)
        {
            solute[k + 1] += solute[k];
            solute[k] = 0.0;
        }
    }
}
