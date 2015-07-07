/*****************************************************************************
 * File		: update.c
 * Function	: Update variables during simulation
 ****************************************************************************/

#include "pihm.h"

void Summary (pihm_struct pihm, N_Vector CV_Y, double stepsize)
{
    double       *y;
    double        wtd0, wtd1, elemsatn0, elemsatn1, realunsat0, realunsat1, realgw0, realgw1, recharge, runoff;
    double        h, aquiferdepth;
    int             i, j;
    elem_struct  *elem;

    y = NV_DATA_S (CV_Y);


    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        elem->unsat = (y[i + pihm->numele] >= 0.0) ? y[i + pihm->numele] : 0.0;
        elem->gw = (y[i + 2 * pihm->numele] >= 0.0) ? y[i + 2 * pihm->numele] : 0.0;

        h = elem->gw;
        /* calculate infiltration based on mass conservation */
        aquiferdepth = elem->topo.zmax - elem->topo.zmin;
        wtd0 = aquiferdepth - (elem->gw0 > 0.0 ? elem->gw0 : 0.0);
        wtd0 = wtd0 < 0.0 ? 0.0 : wtd0;
        elemsatn0 = (wtd0 <= 0.0) ? 1.0 : (elem->unsat0 < 0.0 ? 0.0 : elem->unsat0 / wtd0);
        elemsatn0 = elemsatn0 > 1.0 ? 1.0 : (elemsatn0 < 0.0 ? 0.0 : elemsatn0);
        realunsat0 = elemsatn0 * wtd0;
        realgw0 = elem->gw0 > aquiferdepth ? aquiferdepth : (elem->gw0 < 0.0 ? 0.0 : elem->gw0);
        wtd1 = aquiferdepth - elem->gw;
        wtd1 = wtd1 < 0.0 ? 0.0 : wtd1;
        elemsatn1 = (wtd1 <= 0.0) ? 1.0 : (elem->unsat < 0.0 ? 0.0 : elem->unsat / wtd1);
        elemsatn1 = elemsatn1 > 1.0 ? 1.0 : (elemsatn1 < 0.0 ? 0.0 : elemsatn1);
        realunsat1 = elemsatn1 * wtd1;
        realgw1 = elem->gw > aquiferdepth ? aquiferdepth : elem->gw;

        /* subsurface runoff rate */
        runoff = 0.0;
        for (j = 0; j < 3; j++)
        {
            runoff = runoff + elem->fluxsub[j] / elem->topo.area;
        }
#ifdef _NOAH_
        recharge = (realgw1 - realgw0) * elem->soil.porosity / stepsize + runoff + elem->et_from_sat * elem->et[1];
        elem->infil = (realunsat1 - realunsat0) * elem->soil.porosity / stepsize + recharge + (1.0 - elem->et_from_sat) * elem->et[1] + elem->et[2];
#else
        recharge = (realgw1 - realgw0) * elem->soil.porosity / stepsize + runoff + ((elem->gw0 > aquiferdepth - elem->lc.rzd) ? elem->et[1] : 0.0);
        elem->infil = (realunsat1 - realunsat0) * elem->soil.porosity / stepsize + recharge + (elem->surf0 < EPS / 100.0 ? elem->et[2] : 0.0) + ((elem->gw0 <= aquiferdepth - elem->lc.rzd) ? elem->et[1] : 0.0);
#endif
        elem->infil = elem->infil > 0.0 ? elem->infil : 0.0;
    }
    for (i = 0; i < pihm->numele; i++)
    {
        pihm->elem[i].surf0 = y[i];
        pihm->elem[i].unsat0 = y[i + pihm->numele];
        pihm->elem[i].gw0 = y[i + 2 * pihm->numele];
    }
    for (i = 0; i < pihm->numriv; i++)
    {
        pihm->riv[i].stage0 = y[i + 3 * pihm->numele];
        pihm->riv[i].gw0 = y[i + 3 * pihm->numele + pihm->numriv];
    }
}
