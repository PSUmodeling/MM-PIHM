#include "pihm.h"

void Summary (pihm_struct pihm, N_Vector CV_Y, double stepsize)
{
    double         *y;
    int             i;
    double          subrunoff;

    y = NV_DATA (CV_Y);

#ifdef _OPENMP
#pragma omp parallel for private(subrunoff)
#endif
    for (i = 0; i < nelem; i++)
    {
        pihm->elem[i].ws.surf = y[SURF (i)];
        pihm->elem[i].ws.unsat = y[UNSAT (i)];
        pihm->elem[i].ws.gw = y[GW (i)];

        MassBalance (&pihm->elem[i].ws, &pihm->elem[i].ws0, &pihm->elem[i].wf,
            &subrunoff, &pihm->elem[i].soil, pihm->elem[i].topo.area,
            stepsize);

#ifdef _NOAH_
        pihm->elem[i].wf.runoff2 = subrunoff;

        pihm->elem[i].ps.nwtbl = FindWT (pihm->elem[i].ps.sldpth,
            pihm->elem[i].ps.nsoil, pihm->elem[i].ws.gw,
            pihm->elem[i].ps.satdpth);

        CalcLatFlx (&pihm->elem[i].ps, &pihm->elem[i].wf,
            pihm->elem[i].topo.area);
#endif

        pihm->elem[i].ws0 = pihm->elem[i].ws;

#ifdef _BGC_
        pihm->elem[i].ns.surfn = (y[SURFN (i)] > 0.0) ? y[SURFN (i)] : 0.0;
        pihm->elem[i].ns.sminn = (y[SMINN (i)] > 0.0) ? y[SMINN (i)] : 0.0;

        pihm->elem[i].ns.nleached_snk +=
            (pihm->elem[i].nt.surfn0 + pihm->elem[i].nt.sminn0) -
            (pihm->elem[i].ns.surfn + pihm->elem[i].ns.sminn) +
            pihm->elem[i].nf.ndep_to_sminn * stepsize +
            pihm->elem[i].nf.nfix_to_sminn * stepsize +
            pihm->elem[i].nsol.snksrc * stepsize;

        pihm->elem[i].nt.surfn0 = pihm->elem[i].ns.surfn;
        pihm->elem[i].nt.sminn0 = pihm->elem[i].ns.sminn;
#endif
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        pihm->riv[i].ws.stage = y[RIVSTG (i)];
        pihm->riv[i].ws.gw = y[RIVGW (i)];

        pihm->riv[i].ws0 = pihm->riv[i].ws;
    }
}

void MassBalance (wstate_struct *ws, wstate_struct *ws0, wflux_struct *wf,
    double *subrunoff, const soil_struct *soil, double area, double stepsize)
{
    int             j;
    double          wtd0, wtd1;
    double          elemsatn0, elemsatn1;
    double          realunsat0, realunsat1;
    double          realgw0, realgw1;
    double          recharge;
    double          soilw0, soilw1;

    /*
     * Calculate infiltration based on mass conservation
     */
    wtd0 = soil->depth - ((ws0->gw > 0.0) ? ws0->gw : 0.0);
    wtd0 = (wtd0 < 0.0) ? 0.0 : wtd0;

    if (wtd0 <= 0.0)
    {
        elemsatn0 = 1.0;
    }
    else if (ws0->unsat < 0.0)
    {
        elemsatn0 = 0.0;
    }
    else
    {
        elemsatn0 = ws0->unsat / wtd0;
    }

    elemsatn0 = (elemsatn0 > 1.0) ? 1.0 : elemsatn0;
    elemsatn0 = (elemsatn0 < 0.0) ? 0.0 : elemsatn0;

    realunsat0 = elemsatn0 * wtd0;

    realgw0 = ws0->gw;
    realgw0 = (realgw0 > soil->depth) ? soil->depth : realgw0;
    realgw0 = (realgw0 < 0.0) ? 0.0 : realgw0;

    soilw0 = ws0->gw + ws0->unsat;
    soilw0 = (soilw0 > soil->depth) ? soil->depth : soilw0;
    soilw0 = (soilw0 < 0.0) ? 0.0 : soilw0;

    wtd1 = soil->depth - ((ws->gw > 0.0) ? ws->gw : 0.0);
    wtd1 = (wtd1 < 0.0) ? 0.0 : wtd1;

    if (wtd1 <= 0.0)
    {
        elemsatn1 = 1.0;
    }
    else if (ws->unsat < 0.0)
    {
        elemsatn1 = 0.0;
    }
    else
    {
        elemsatn1 = ws->unsat / wtd1;
    }

    elemsatn1 = (elemsatn1 > 1.0) ? 1.0 : elemsatn1;
    elemsatn1 = (elemsatn1 > 0.0) ? 0.0 : elemsatn1;

    realunsat1 = elemsatn1 * wtd1;

    realgw1 =
        (ws->gw > soil->depth) ? soil->depth : ws->gw;

    soilw1 = ws->gw + ws->unsat;
    soilw1 = (soilw1 > soil->depth) ? soil->depth : soilw1;
    soilw1 = (soilw1 < 0.0) ? 0.0 : soilw1;

    /*
     * Subsurface runoff rate
     */
    *subrunoff = 0.0;
    for (j = 0; j < NUM_EDGE; j++)
    {
        *subrunoff += wf->subsurf[j] / area;
    }

    recharge = (realgw1 - realgw0) * soil->porosity / stepsize +
        *subrunoff + wf->edir_gw + wf->ett_gw;

    wf->infil = (soilw1 - soilw0) * soil->porosity / stepsize +
        *subrunoff + wf->edir_unsat + wf->edir_gw +
        wf->ett_unsat + wf->ett_gw;

    if (wf->infil < 0.0)
    {
        *subrunoff -= wf->infil;
        wf->infil = 0.0;
    }
}
