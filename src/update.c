
/*****************************************************************************
 * File		: update.c
 * Function	: Update variables during simulation
 ****************************************************************************/

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
    double          totalw0;
    int             i, j;
    elem_struct    *elem;
    river_struct   *riv;

    y = NV_DATA_S (CV_Y);


    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        elem->surf = y[SURF(i)];
        elem->unsat = y[UNSAT(i)];
        elem->gw = y[GW(i)];

        /* calculate infiltration based on mass conservation */
        wtd0 = elem->soil.depth - ((elem->gw0 > 0.0) ? elem->gw0 : 0.0);
        wtd0 = (wtd0 < 0.0) ? 0.0 : wtd0;
        if (wtd0 <= 0.0)
        {
            elemsatn0 = 1.0;
        }
        else if (elem->unsat0 < 0.0)
        {
            elemsatn0 = 0.0;
        }
        else
        {
            elemsatn0 = elem->unsat0 / wtd0;
        }
        elemsatn0 = (elemsatn0 > 1.0) ? 1.0 : elemsatn0;
        elemsatn0 = (elemsatn0 < 0.0) ? 0.0 : elemsatn0;
        realunsat0 = elemsatn0 * wtd0;

        realgw0 = elem->gw0;
        realgw0 = (realgw0 > elem->soil.depth) ? elem->soil.depth : realgw0;
        realgw0 = (realgw0 < 0.0) ? 0.0 : realgw0;

        soilw0 = elem->gw0 + elem->unsat0;
        soilw0 = (soilw0 > elem->soil.depth) ? elem->soil.depth : soilw0;
        soilw0 = (soilw0 < 0.0) ? 0.0 : soilw0;

        wtd1 = elem->soil.depth - ((elem->gw > 0.0) ? elem->gw : 0.0);
        wtd1 = (wtd1 < 0.0) ? 0.0 : wtd1;
        if (wtd1 <= 0.0)
        {
            elemsatn1 = 1.0;
        }
        else if (elem->unsat < 0.0)
        {
            elemsatn1 = 0.0;
        }
        else
        {
            elemsatn1 = elem->unsat / wtd1;
        }
        elemsatn1 = (elemsatn1 > 1.0) ? 1.0 : elemsatn1;
        elemsatn1 = (elemsatn1 > 0.0) ? 0.0 : elemsatn1;
        realunsat1 = elemsatn1 * wtd1;

        realgw1 = (elem->gw > elem->soil.depth) ? elem->soil.depth : elem->gw;

        soilw1 = elem->gw + elem->unsat;
        soilw1 = (soilw1 > elem->soil.depth) ? elem->soil.depth : soilw1;
        soilw1 = (soilw1 < 0.0) ? 0.0 : soilw1;
        /* subsurface runoff rate */
        elem->runoff = 0.0;
        for (j = 0; j < 3; j++)
        {
            elem->runoff += elem->fluxsub[j] / elem->topo.area;
        }

        recharge = (realgw1 - realgw0) * elem->soil.porosity / stepsize +
            elem->runoff + elem->edir[2] + elem->ett[2];
        //elem->infil = (realunsat1 - realunsat0) * elem->soil.porosity / stepsize + recharge + (1.0 - elem->et_from_sat) * elem->et[1] + elem->et[2];
        elem->infil = (soilw1 - soilw0) * elem->soil.porosity / stepsize +
            elem->runoff + elem->edir[1] +elem->edir[2] + elem->ett[1] + elem->ett[2];

        if (elem->infil < 0.0)
        {
            elem->runoff -= elem->infil;
            elem->infil = 0.0;
        }

        totalw0 = elem->surf0 + elem->unsat0 * elem->soil.porosity +
            elem->gw0 * elem->soil.porosity + stepsize *
            (elem->netprcp - elem->et[0] - elem->et[1] - elem->et[2]);
        for (j = 0; j < 3; j++)
        {
            totalw0 -= (elem->fluxsurf[j] + elem->fluxsub[j]) / elem->topo.area;
        }

        elem->totalw = elem->surf + elem->unsat * elem->soil.porosity + elem->gw * elem->soil.porosity;

        elem->mbc = elem->totalw / totalw0;
    }

    for (i = 0; i < pihm->numriv; i++)
    {
        riv = &pihm->riv[i];

        riv->stage = y[RIVSTG(i)];
        riv->gw = y[RIVGW(i)];

        totalw0 = riv->stage0 + riv->gw0 * riv->matl.porosity;

        for (j = 0; j < 7; j++)
        {
            totalw0 -= riv->fluxriv[j] / (riv->shp.length * EqWid (riv->shp.intrpl_ord, riv->shp.depth, riv->shp.coeff));
        }
        totalw0 += (-riv->fluxriv[7] - riv->fluxriv[8] - riv->fluxriv[9] - riv->fluxriv[10] + riv->fluxriv[6]) / (riv->shp.length * EqWid (riv->shp.intrpl_ord, riv->shp.depth, riv->shp.coeff));

        riv->totalw = riv->stage + riv->gw * riv->matl.porosity;

        riv->mbc = riv->totalw / totalw0;
    }

    for (i = 0; i < pihm->numele; i++)
    {
        pihm->elem[i].surf0 = y[SURF(i)];
        pihm->elem[i].unsat0 = y[UNSAT(i)];
        pihm->elem[i].gw0 = y[GW(i)];
    }
    for (i = 0; i < pihm->numriv; i++)
    {
        pihm->riv[i].stage0 = y[RIVSTG(i)];
        pihm->riv[i].gw0 = y[RIVGW(i)];
    }
}
