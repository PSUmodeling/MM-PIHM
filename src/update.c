#include "pihm.h"

void UpdateVar(double stepsize, elem_struct elem[], river_struct river[], cvode_struct *cvode)
{
    double         *y;
    int             i;

    y = NV_DATA(cvode->CV_Y);

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        double          subrunoff;

        elem[i].ws.surf = y[SURF(i)];
        elem[i].ws.unsat = y[UNSAT(i)];
        elem[i].ws.gw = y[GW(i)];

        AdjustFluxes(elem[i].topo.area, stepsize, &elem[i].soil, &elem[i].ws, &elem[i].ws0, &subrunoff, &elem[i].wf);

#if defined(_NOAH_)
        elem[i].wf.runoff2 = subrunoff;

        elem[i].ps.nwtbl = FindWaterTable(elem[i].ps.nlayers, elem[i].ws.gw, elem[i].ps.soil_depth, elem[i].ps.satdpth);

        CalcLateralFlux(&elem[i].ps, &elem[i].wf);
#endif

        elem[i].ws0 = elem[i].ws;

#if defined(_TRANSPORT_) && !defined(_BGC_) && !defined(_CYCLES_)
        int             k;

        for (k = 0; k < nsolute; k++)
        {
            elem[i].solute[k].amount = MAX(y[SOLUTE_SOIL(i, k)], 0.0);
        }
#endif

#if defined(_BGC_)
        elem[i].ns.sminn = MAX(y[SOLUTE_SOIL(i, 0)], 0.0);

        elem[i].ns.nleached_snk += elem[i].nt.sminn0 - elem[i].ns.sminn + elem[i].solute[0].snksrc * stepsize;

        elem[i].nt.sminn0 = elem[i].ns.sminn;
#endif

#if defined(_CYCLES_)
        elem[i].ps.no3 = MAX(y[SOLUTE_SOIL(i, NO3)], 0.0);
        elem[i].ps.nh4 = MAX(y[SOLUTE_SOIL(i, NH4)], 0.0);

        UpdateNProfile(stepsize, &elem[i].soil, &elem[i].ws, &elem[i].ns, elem[i].solute, elem[i].ns.no3, elem[i].ns.nh4, &elem[i].ps, &elem[i].nf);

        elem[i].ns0 = elem[i].ns;
        elem[i].ps.no3_prev = elem[i].ps.no3;
        elem[i].ps.nh4_prev = elem[i].ps.nh4;
#endif
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river[i].ws.stage = y[RIVER(i)];

#if defined(_TRANSPORT_) && !defined(_BGC_) && !defined(_CYCLES_)
        int             j;

        for (j = 0; j < nsolute; j++)
        {
            river[i].solute[j].amount = MAX(y[SOLUTE_RIVER(i, j)], 0.0);
        }
#endif
#if defined(_BGC_)
        river[i].ns.streamn = MAX(y[SOLUTE_RIVER(i, 0)], 0.0);
#endif

#if defined(_CYCLES_)
        river[i].ns.no3 = MAX(y[SOLUTE_RIVER(i, NO3)], 0.0);
        river[i].ns.nh4 = MAX(y[SOLUTE_RIVER(i, NH4)], 0.0);
#endif
    }
}

// The AdjustFluxes subroutine adjusts model infiltration rate and recharge rate based on mass balance in the soil
// column. This subroutine also calculates an equivalent infiltration rate, which is used as the boundary condition
// in the Noah soil moisture calculation. Unlike the infiltration rate, the equivalent infiltration rate is calculated
// without considering oversturation, and is based on the mass balance of the whole soil column, instead of just the
// unsaturated zone.
void AdjustFluxes(double area, double stepsize, const soil_struct *soil, const wstate_struct *ws, const wstate_struct *ws0, double *subrunoff, wflux_struct *wf)
{
    int             j;
    double          soilw0, soilw1;

    // Adjust recharge
    soilw0 = ws0->gw;
    soilw1 = ws->gw;

    *subrunoff = 0.0;
    for (j = 0; j < NUM_EDGE; j++)
    {
        *subrunoff += wf->subsurf[j] / area;
    }

    wf->recharge = (soilw1 - soilw0) * soil->porosity / stepsize + *subrunoff + wf->edir_gw + wf->ett_gw;

    // Adjust infiltration
    soilw0 = ws0->unsat;
    soilw1 = ws->unsat;

    wf->infil = (soilw1 - soilw0) * soil->porosity / stepsize + wf->recharge + wf->edir_unsat + wf->ett_unsat;

    // Calculate equivalent infiltration based on mass conservation
    soilw0 = ws0->gw + ws0->unsat;
    soilw0 = MIN(soilw0, soil->depth);
    soilw0 = MAX(soilw0, 0.0);

    soilw1 = ws->gw + ws->unsat;
    soilw1 = MIN(soilw1, soil->depth);
    soilw1 = MAX(soilw1, 0.0);

    wf->eqv_infil = (soilw1 - soilw0) * soil->porosity / stepsize + *subrunoff + wf->edir_unsat + wf->edir_gw + wf->ett_unsat + wf->ett_gw;

    if (wf->eqv_infil < 0.0)
    {
        *subrunoff -= wf->eqv_infil;
        wf->eqv_infil = 0.0;
    }
}
