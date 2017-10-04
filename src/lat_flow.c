#include "pihm.h"

void LateralFlow(elem_struct *elem, river_struct *rivseg, int surf_mode)
{
    int             i;
    double         *dhbydx;
    double         *dhbydy;

    dhbydx = (double *)malloc(nelem * sizeof(double));
    dhbydy = (double *)malloc(nelem * sizeof(double));

    FrictSlope(elem, rivseg, surf_mode, dhbydx, dhbydy);

#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;
        double          avg_sf;
        elem_struct    *nabr;

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].nabr[j] > 0)
            {
                nabr = &elem[elem[i].nabr[j] - 1];

                /* Subsurface flow between triangular elements */
                elem[i].wf.subsurf[j] =
                    SubFlowElemToElem(&elem[i].ws, &elem[i].topo,
                    &elem[i].soil, j, &nabr->ws, &nabr->topo, &nabr->soil);

                /* Surface flux between triangular elements */
                avg_sf = 0.5 *
                    (sqrt(dhbydx[i] * dhbydx[i] + dhbydy[i] * dhbydy[i]) +
                     sqrt(dhbydx[nabr->ind - 1] * dhbydx[nabr->ind - 1] +
                     dhbydy[nabr->ind - 1] * dhbydy[nabr->ind - 1]));
                elem[i].wf.ovlflow[j] =
                    OvlFlowElemToElem(&elem[i].ws, &elem[i].topo, &elem[i].lc,
                    j, &nabr->ws, &nabr->topo, &nabr->lc, avg_sf, surf_mode);
            }
            else if (elem[i].nabr[j] < 0)
            {
                /* Do nothing. River-element interactions are calculated
                 * in river_flow.c */
            }
            else    /* Boundary condition flux */
            {
                BoundFluxElem(elem[i].attrib.bc_type[j], j, &elem[i].bc,
                    &elem[i].ws, &elem[i].topo, &elem[i].soil, &elem[i].wf);
            }
        }    /* End of neighbor loop */
    }    /* End of element loop */

    free(dhbydx);
    free(dhbydy);
}

void FrictSlope(elem_struct *elem, river_struct *rivseg, int surf_mode,
    double *dhbydx, double *dhbydy)
{
    int             i;
#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;
        double          surfh[NUM_EDGE];
        elem_struct    *nabr;
        river_struct   *rivnabr;

        if (surf_mode == DIFF_WAVE)
        {
            for (j = 0; j < NUM_EDGE; j++)
            {
                if (elem[i].nabr[j] > 0)
                {
                    nabr = &elem[elem[i].nabr[j] - 1];
                    surfh[j] = nabr->topo.zmax + nabr->ws.surfh;
                }
                else if (elem[i].nabr[j] < 0)
                {
                    rivnabr = &rivseg[-elem[i].nabr[j] - 1];

                    if (rivnabr->ws.stage > rivnabr->shp.depth)
                    {
                        surfh[j] = rivnabr->topo.zbed + rivnabr->ws.stage;
                    }
                    else
                    {
                        surfh[j] = rivnabr->topo.zmax;
                    }
                }
                else
                {
                    if (elem[i].attrib.bc_type[j] == 0)
                    {
                        surfh[j] = elem[i].topo.zmax + elem[i].ws.surfh;
                    }
                    else
                    {
                        surfh[j] = elem[i].bc.head[j];
                    }
                }
            }

            dhbydx[i] =
                DhByDl(elem[i].topo.nabrdist_y, elem[i].topo.nabrdist_x, surfh);
            dhbydy[i] =
                DhByDl(elem[i].topo.nabrdist_x, elem[i].topo.nabrdist_y, surfh);
        }
    }
}

double AvgHsurf(double diff, double yi, double yinabr)
{
    double          avg_y;

    if (diff > 0.0)
    {
        if (yi > DEPRSTG)
        {
            avg_y = 1.0 * (yi - DEPRSTG);
        }
        else
        {
            avg_y = 0.0;
        }
    }
    else
    {
        if (yinabr > DEPRSTG)
        {
            avg_y = 1.0 * (yinabr - DEPRSTG);
        }
        else
        {
            avg_y = 0.0;
        }
    }

    return avg_y;
}

double AvgH(double diff, double yi, double yinabr)
{
    double          avg_y = 0.0;

    if (diff > 0.0)
    {
        if (yi > 0.0)
        {
            avg_y = yi;
        }
    }
    else
    {
        if (yinabr > 0.0)
        {
            avg_y = yinabr;
        }
    }

    return avg_y;
}

double DhByDl(double *l1, double *l2, double *surfh)
{
    return -1.0 *
        (l1[2] * (surfh[1] - surfh[0]) + l1[1] * (surfh[0] - surfh[2]) +
        l1[0] * (surfh[2] - surfh[1])) /
        (l2[2] * (l1[1] - l1[0]) + l2[1] * (l1[0] - l1[2]) +
        l2[0] * (l1[2] - l1[1]));
}

double EffKh(double tmpy, double aqdepth, double macd, double macksath,
    double areaf, double ksath)
{
    tmpy = (tmpy > 0.0) ? tmpy : 0.0;

    if (tmpy > aqdepth - macd)
    {
        if (tmpy > aqdepth)
        {
            return (macksath * macd * areaf + ksath * (aqdepth -
                macd * areaf)) / aqdepth;
        }
        else
        {
            return (macksath * (tmpy - (aqdepth - macd)) * areaf +
                ksath * (aqdepth - macd + (tmpy - (aqdepth - macd)) *
                (1.0 - areaf))) / tmpy;
        }
    }
    else
    {
        return ksath;
    }
}


double OverLandFlow(double avg_y, double grad_y, double avg_sf, double crossa,
    double avg_rough)
{
    return crossa * pow(avg_y, 0.6666667) * grad_y / (sqrt(avg_sf) * avg_rough);
}

double SubFlowElemToElem(const wstate_struct *ws, const topo_struct *topo,
    const soil_struct *soil, int j, const wstate_struct *nabr_ws,
    const topo_struct *nabr_topo, const soil_struct *nabr_soil)
{
    double          diff_h;
    double          avg_h;
    double          grad_h;
    double          effk, effk_nabr;
    double          avg_ksat;

    /*
     * Subsurface lateral flux calculation between triangular
     * elements
     */
    diff_h = (ws->gw + topo->zmin) - (nabr_ws->gw + nabr_topo->zmin);
    avg_h = AvgH(diff_h, ws->gw, nabr_ws->gw);
    grad_h = diff_h / topo->nabrdist[j];

    /* Take into account macropore effect */
    effk = EffKh(ws->gw, soil->depth, soil->dmac, soil->kmach, soil->areafv,
        soil->ksath);
    effk_nabr = EffKh(nabr_ws->gw, nabr_soil->depth, nabr_soil->dmac,
        nabr_soil->kmach, nabr_soil->areafv, nabr_soil->ksath);
    avg_ksat = 0.5 * (effk + effk_nabr);

    /* Groundwater flow modeled by Darcy's Law */
    return avg_ksat * grad_h * avg_h * topo->edge[j];
}

double OvlFlowElemToElem(const wstate_struct *ws, const topo_struct *topo,
    const lc_struct *lc, int j, const wstate_struct *nabr_ws,
    const topo_struct *nabr_topo, const lc_struct *nabr_lc, double avg_sf,
    int surf_mode)
{
    double          diff_h;
    double          avg_h;
    double          grad_h;
    double          avg_rough;
    double          crossa;

    diff_h = (surf_mode == KINEMATIC) ?
        topo->zmax - nabr_topo->zmax :
        (ws->surfh + topo->zmax) - (nabr_ws->surfh + nabr_topo->zmax);
    avg_h = AvgHsurf(diff_h, ws->surfh, nabr_ws->surfh);
    grad_h = diff_h / topo->nabrdist[j];
    if (surf_mode == KINEMATIC)
    {
        avg_sf = (grad_h > 0.0) ? grad_h : GRADMIN;
    }
    else
    {
        avg_sf = (avg_sf > GRADMIN) ? avg_sf : GRADMIN;
    }
    /* Weighting needed */
    avg_rough = 0.5 * (lc->rough + nabr_lc->rough);
    crossa = avg_h * topo->edge[j];

    return OverLandFlow(avg_h, grad_h, avg_sf, crossa, avg_rough);
}

void BoundFluxElem(int bc_type, int j, const bc_struct *bc,
    const wstate_struct *ws, const topo_struct *topo,
    const soil_struct *soil, wflux_struct *wf)
{
    double          diff_h;
    double          avg_h;
    double          effk;
    double          avg_ksat;
    double          grad_h;

    /* No flow (natural) boundary condition is default */
    if (bc_type == 0)
    {
        wf->ovlflow[j] = 0.0;
        wf->subsurf[j] = 0.0;
    }
    /* Note: ideally different boundary conditions need to be
     * incorporated for surf and subsurf respectively */
    else if (bc_type > 0)
    {
        /* Note: the formulation assumes only Dirichlet TS right now */
        /* note the assumption here is no flow for surface */
        wf->ovlflow[j] = 0.0;

        diff_h = ws->gw + topo->zmin - bc->head[j];
        avg_h = AvgH(diff_h, ws->gw, bc->head[j] - topo->zmin);
        /* Minimum distance from circumcenter to the edge of the triangle
         * on which boundary condition is defined */
        effk = EffKh(ws->gw, soil->depth, soil->dmac, soil->kmach,
            soil->areafv, soil->ksath);
        avg_ksat = effk;
        grad_h = diff_h / topo->nabrdist[j];
        wf->subsurf[j] = avg_ksat * grad_h * avg_h * topo->edge[j];
    }
    else
    {
        /* Neumann bc (note: md->ele[i].bc[j] value has to be
         * = 2+(index of Neumann boundary ts) */
        wf->ovlflow[j] = 0.0;
        wf->subsurf[j] = bc->flux[j];
    }
}
