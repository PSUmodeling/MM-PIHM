#include "pihm.h"

void Summary (pihm_struct pihm, N_Vector CV_Y, double stepsize)
{
    double         *y;
    double          wtd0, wtd1;
    double          elemsatn0, elemsatn1;
    double          realunsat0, realunsat1;
    double          realgw0, realgw1;
    double          recharge;
    double          soilw0, soilw1;
    double          subrunoff;
    int             i, j;
    elem_struct    *elem;
    river_struct   *riv;

    y = NV_DATA_S (CV_Y);

    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        elem->ws.surf = y[SURF (i)];
        elem->ws.unsat = y[UNSAT (i)];
        elem->ws.gw = y[GW (i)];

        /*
         * Calculate infiltration based on mass conservation
         */
        wtd0 = elem->soil.depth - ((elem->ws0.gw > 0.0) ? elem->ws0.gw : 0.0);
        wtd0 = (wtd0 < 0.0) ? 0.0 : wtd0;
        if (wtd0 <= 0.0)
        {
            elemsatn0 = 1.0;
        }
        else if (elem->ws0.unsat < 0.0)
        {
            elemsatn0 = 0.0;
        }
        else
        {
            elemsatn0 = elem->ws0.unsat / wtd0;
        }
        elemsatn0 = (elemsatn0 > 1.0) ? 1.0 : elemsatn0;
        elemsatn0 = (elemsatn0 < 0.0) ? 0.0 : elemsatn0;
        realunsat0 = elemsatn0 * wtd0;

        realgw0 = elem->ws0.gw;
        realgw0 = (realgw0 > elem->soil.depth) ? elem->soil.depth : realgw0;
        realgw0 = (realgw0 < 0.0) ? 0.0 : realgw0;

        soilw0 = elem->ws0.gw + elem->ws0.unsat;
        soilw0 = (soilw0 > elem->soil.depth) ? elem->soil.depth : soilw0;
        soilw0 = (soilw0 < 0.0) ? 0.0 : soilw0;

        wtd1 = elem->soil.depth - ((elem->ws.gw > 0.0) ? elem->ws.gw : 0.0);
        wtd1 = (wtd1 < 0.0) ? 0.0 : wtd1;
        if (wtd1 <= 0.0)
        {
            elemsatn1 = 1.0;
        }
        else if (elem->ws.unsat < 0.0)
        {
            elemsatn1 = 0.0;
        }
        else
        {
            elemsatn1 = elem->ws.unsat / wtd1;
        }
        elemsatn1 = (elemsatn1 > 1.0) ? 1.0 : elemsatn1;
        elemsatn1 = (elemsatn1 > 0.0) ? 0.0 : elemsatn1;
        realunsat1 = elemsatn1 * wtd1;

        realgw1 =
            (elem->ws.gw > elem->soil.depth) ? elem->soil.depth : elem->ws.gw;

        soilw1 = elem->ws.gw + elem->ws.unsat;
        soilw1 = (soilw1 > elem->soil.depth) ? elem->soil.depth : soilw1;
        soilw1 = (soilw1 < 0.0) ? 0.0 : soilw1;

        /* 
         * Subsurface runoff rate
         */
        subrunoff = 0.0;
        for (j = 0; j < 3; j++)
        {
            subrunoff += elem->wf.subsurf[j] / elem->topo.area;
        }

        recharge = (realgw1 - realgw0) * elem->soil.porosity / stepsize +
            subrunoff + elem->wf.edir_gw + elem->wf.ett_gw;

        elem->wf.infil = (soilw1 - soilw0) * elem->soil.porosity / stepsize +
            subrunoff + elem->wf.edir_unsat + elem->wf.edir_gw +
            elem->wf.ett_unsat + elem->wf.ett_gw;

        if (elem->wf.infil < 0.0)
        {
            subrunoff -= elem->wf.infil;
            elem->wf.infil = 0.0;
        }

#ifdef _NOAH_
        elem->wf.runoff2 = subrunoff;

        elem->ps.nwtbl = FindWT (elem->ps.sldpth, elem->ps.nsoil,
            elem->ws.gw, elem->ps.satdpth);

        CalcLatFlx (&elem->ws, &elem->ps, &elem->wf);
#endif

        elem->ws0 = elem->ws;
    }

    for (i = 0; i < pihm->numriv; i++)
    {
        riv = &pihm->riv[i];

        riv->ws.stage = y[RIVSTG (i)];
        riv->ws.gw = y[RIVGW (i)];

        riv->ws0 = riv->ws;
    }

#ifdef _NOAH_
    AvgFlux (pihm->elem, pihm->numele, SUM);
#endif
}
