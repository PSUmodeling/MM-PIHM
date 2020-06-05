#include "pihm.h"

void SoluteConc(double dt, elem_struct elem[], river_struct river[])
{
    int             i;
    /*
     * Calculate solute N concentrations
     */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             k;
        double          no3[MAXLYR];
        double          nh4[MAXLYR];

        for (k = 0; k < MAXLYR; k++)
        {
            no3[k] = 0.0;
            nh4[k] = 0.0;
        }

        UpdateNProfile(dt, &elem[i].soil, &elem[i].ws, &elem[i].ns,
            elem[i].solute, no3, nh4, &elem[i].ps);

        /* Calculate NO3 and NH4 concentrations */
        elem[i].solute[NO3].conc_surf = 0.0;
        elem[i].solute[NO3].conc = MobileNConc(KD_NO3, no3, &elem[i].soil,
            &elem[i].ws, &elem[i].ps);

        elem[i].solute[NH4].conc_surf = 0.0;
        elem[i].solute[NH4].conc = MobileNConc(KD_NH4, nh4, &elem[i].soil,
            &elem[i].ws, &elem[i].ps);
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river[i].solute[NO3].conc = (river[i].ws.stage > 0.0) ?
            river[i].ns.no3 / river[i].ws.stage : 0.0;
        river[i].solute[NO3].conc = MAX(river[i].solute[NO3].conc, 0.0);

        river[i].solute[NH4].conc = (river[i].ws.stage > 0.0) ?
            river[i].ns.nh4 / river[i].ws.stage : 0.0;
        river[i].solute[NH4].conc = MAX(river[i].solute[NH4].conc, 0.0);
    }
}

void UpdateNProfile(double dt, const soil_struct *soil,
    const wstate_struct *ws, const nstate_struct *ns,
    const solute_struct solute[], double no3[], double nh4[],
    phystate_struct *ps)
{
    int             k;

    /*
    * Add source sink terms to NO3 and NH4 mass at different layers
    */
    for (k = 0; k < ps->nlayers; k++)
    {
        no3[k] = ns->no3[k] + solute[NO3].snksrc[k] * dt;

        nh4[k] = ns->nh4[k] + solute[NH4].snksrc[k] * dt;
    }

    /*
    * Add lateral transport fluxes to NO3 and NH4 mass
    */
    FindWaterTable(ps->soil_depth, ps->nlayers, ws->gw, ps->satdpth);

    LateralNFlow(KD_NO3, soil, ws, ps, ns->no3, Profile(ps->nlayers, no3),
        ps->no3, no3);

    LateralNFlow(KD_NH4, soil, ws, ps, ns->nh4, Profile(ps->nlayers, nh4),
        ps->nh4, nh4);
}

void LateralNFlow(double kd, const soil_struct *soil,
    const wstate_struct *ws, const phystate_struct *ps,
    const double solute0[], double profile0, double profile, double solute[])
{
    int             k;
    double          total_weight = 0.0;
    double          weight[MAXLYR];

    for (k = 0; k < ps->nlayers; k++)
    {
        weight[k] = (ps->satdpth[k] > 0.0 && solute[k] > 0.0) ?
            LinearEqmConc(kd, soil->bd[k], ps->soil_depth[k], ws->swc[k],
            solute[k]) * ps->satdpth[k] * ws->swc[k] * RHOH2O : 0.0;
    }

    total_weight = Profile(ps->nlayers, weight);

    if (total_weight <= 0.0)
    {
        total_weight = 1.0;
        weight[ps->nlayers - 1] = 1.0;
    }

    for (k = 0; k < ps->nlayers; k++)
    {
        weight[k] /= total_weight;
        solute[k] += weight[k] * (profile - profile0);
    }

    for (k = ps->nlayers - 1; k > 0; k--)
    {
        if (solute[k] < 0.0)
        {
            solute[k - 1] += solute[k];
            solute[k] = 0.0;
        }
    }

    for (k = 0; k < ps->nlayers - 1; k++)
    {
        if (solute[k] < 0.0)
        {
            solute[k + 1] += solute[k];
            solute[k] = 0.0;
        }
    }
}

double MobileNConc(double kd, const double solute[], const soil_struct *soil,
    const wstate_struct *ws, const phystate_struct *ps)
{
    int             k;
    double          conc;
    double          avg_conc = 0.0;
    double          storage = 0.0;

    for (k = 0; k < ps->nlayers; k++)
    {
        if (ps->satdpth[k] > 0.0)
        {
            conc = (solute[k] > 0.0) ?
                LinearEqmConc(kd, soil->bd[k], ps->soil_depth[k], ws->swc[k],
                    solute[k]) * RHOH2O : 0.0;

            avg_conc += ps->satdpth[k] * conc;
            storage += ps->satdpth[k];
        }
    }

    avg_conc /= storage;

    return avg_conc;
}

#if NOT_YET_IMPLEMETED
void NRT(double kd, double dt, double wflux[], const soil_struct *soil,
    const wstate_struct *ws0, const wstate_struct *ws,
    const phystate_struct *ps, double solute[])
{
    int             k;
    double          adsorb[MAXLYR];
    double          conc[MAXLYR];
    double          conc0[MAXLYR];
    double          ai[MAXLYR], bi[MAXLYR], ci[MAXLYR], rhstt[MAXLYR];
    double          ciin[MAXLYR], rhsttin[MAXLYR];

    for (k = 0; k < ps->nlayers; k++)
    {
        conc[k] = (solute[k] > 0.0) ? LinearEqmConc(kd, soil->bd[k],
            ps->soil_depth[k], ws0->smc[k], solute[k]) : 0.0;
        adsorb[k] = solute[k] -
            ps->soil_depth[k] * RHOH2O * ws0->smc[k] * conc[k];

        ai[k] = -0.5 * wflux[k] / ps->soil_depth[k] / RHOH2O;
        bi[k] = ws->smc[k] / dt +
            0.5 * wflux[k + 1] / ps->soil_depth[k] / RHOH2O;
        ci[k] = 0.0;
        rhstt[k] = -(0.5 * wflux[k + 1] / ps->soil_depth[k] / RHOH2O -
            ws0->smc[k] / dt) * conc[k] + 0.0 / ps->soil_depth[k] / RHOH2O;
        rhstt[k] += (k == 0) ? 0.0 :
            0.5 * wflux[k] / ps->soil_depth[k] / RHOH2O * conc[k - 1];

        rhsttin[k] = rhstt[k];
        ciin[k] = ci[k];
    }

    Rosr12(ci, ai, bi, ciin, rhsttin, rhstt, ps->nlayers);

    for (k = 0; k < ps->nlayers; k++)
    {
        conc[k] = ci[k];

        solute[k] = adsorb[k] +
            ps->soil_depth[k] * RHOH2O * ws->smc[k] * conc[k];
    }

    for (k = ps->nlayers - 1; k > 0; k--)
    {
        if (solute[k] < 0.0)
        {
            solute[k - 1] += solute[k];
            solute[k] = 0.0;
        }
    }

    for (k = 0; k < ps->nlayers - 1; k++)
    {
        if (solute[k] < 0.0)
        {
            solute[k + 1] += solute[k];
            solute[k] = 0.0;
        }
    }
}
#endif
