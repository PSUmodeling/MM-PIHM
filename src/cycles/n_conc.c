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

        UpdateNProfile(dt, &elem[i].soil, &elem[i].ws, &elem[i].ns0,
            &elem[i].nf, &elem[i].ns, &elem[i].ps);

        /* Calculate NO3 and NH4 concentrations */
        elem[i].solute[NO3].conc_surf = 0.0;
        elem[i].solute[NO3].conc = MobileNConc(KD_NO3, elem[i].ns.no3,
            &elem[i].soil, &elem[i].ws, &elem[i].ps);

        elem[i].solute[NH4].conc_surf = 0.0;
        elem[i].solute[NH4].conc = MobileNConc(KD_NH4, elem[i].ns.nh4,
            &elem[i].soil, &elem[i].ws, &elem[i].ps);
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river[i].solute[NO3].conc = (river[i].ws.stage > DEPTHR) ?
            river[i].ns.no3 / river[i].ws.stage : 0.0;
        river[i].solute[NO3].conc = MAX(river[i].solute[NO3].conc, 0.0);

        river[i].solute[NH4].conc = (river[i].ws.stage > DEPTHR) ?
            river[i].ns.nh4 / river[i].ws.stage : 0.0;
        river[i].solute[NH4].conc = MAX(river[i].solute[NH4].conc, 0.0);
    }
}

void UpdateNProfile(double dt, const soil_struct *soil,
    const wstate_struct *ws, const nstate_struct *ns0, const nflux_struct *nf,
    nstate_struct *ns, phystate_struct *ps)
{
    int             k;

    /*
    * Add source sink terms to NO3 and NH4 mass at different layers
    */
    ns->no3[0] = ns0->no3[0] + nf->surplus / DAYINSEC * dt;
    ns->nh4[0] = ns0->nh4[0] + nf->urine / DAYINSEC * dt;

    for (k = 0; k < ps->nlayers; k++)
    {
        ns->no3[k] = ns0->no3[k] +
            (nf->nitrif[k] - nf->n2o_from_nitrif[k] -
            nf->denitrif[k] - nf->no3_uptake[k] +
            nf->no3_fert[k] + nf->no3_immobil[k]) /
            DAYINSEC * dt;

        ns->nh4[k] = ns0->nh4[k] +
            (-nf->nitrif[k] - nf->volatil[k] -
            nf->nh4_uptake[k] + nf->nh4_fert[k] +
            nf->nh4_immobil[k] + nf->mineral[k]) /
            DAYINSEC * dt;
    }

    /*
    * Add lateral transport fluxes to NO3 and NH4 mass
    */
    FindWaterTable(ps->soil_depth, ps->nlayers, ws->gw, ps->satdpth);

    LateralNFlow(KD_NO3, soil, ws, ps, ns0->no3, ns0->no3_profile,
        ns->no3_profile, ns->no3);

    LateralNFlow(KD_NH4, soil, ws, ps, ns0->nh4, ns0->nh4_profile,
        ns->nh4_profile, ns->nh4);
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
            solute0[k]) * ps->satdpth[k] * ws->swc[k] * RHOH2O : 0.0;
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
                    solute[k]) : 0.0;

            avg_conc += ps->satdpth[k] * conc;
            storage += ps->satdpth[k];
        }
    }

    avg_conc /= storage;

    return avg_conc;
}
