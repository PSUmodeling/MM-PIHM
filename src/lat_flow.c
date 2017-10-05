#include "pihm.h"

void LateralFlow(elem_struct *elem, river_struct *river, int surf_mode)
{
    int             i;
    double         *dhbydx;
    double         *dhbydy;

    dhbydx = (double *)malloc(nelem * sizeof(double));
    dhbydy = (double *)malloc(nelem * sizeof(double));

    FrictSlope(elem, river, surf_mode, dhbydx, dhbydy);

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
                    SubFlowElemToElem(&elem[i], nabr, j);

                /* Surface flux between triangular elements */
                avg_sf = 0.5 *
                    (sqrt(dhbydx[i] * dhbydx[i] + dhbydy[i] * dhbydy[i]) +
                     sqrt(dhbydx[nabr->ind - 1] * dhbydx[nabr->ind - 1] +
                     dhbydy[nabr->ind - 1] * dhbydy[nabr->ind - 1]));
                elem[i].wf.ovlflow[j] =
                    OvlFlowElemToElem(&elem[i], nabr, j, avg_sf, surf_mode);
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

void FrictSlope(elem_struct *elem, river_struct *river, int surf_mode,
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
                    rivnabr = &river[-elem[i].nabr[j] - 1];

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

double AvgHsurf(double diff, double hsurf, double hnabr)
{
    double          avg_h;

    if (diff > 0.0)
    {
        if (hsurf > DEPRSTG)
        {
            avg_h = 1.0 * (hsurf - DEPRSTG);
        }
        else
        {
            avg_h = 0.0;
        }
    }
    else
    {
        if (hnabr > DEPRSTG)
        {
            avg_h = 1.0 * (hnabr - DEPRSTG);
        }
        else
        {
            avg_h = 0.0;
        }
    }

    return avg_h;
}

double AvgH(double diff, double hsub, double hnabr)
{
    double          avg_h = 0.0;

    if (diff > 0.0)
    {
        if (hsub > 0.0)
        {
            avg_h = hsub;
        }
    }
    else
    {
        if (hnabr > 0.0)
        {
            avg_h = hnabr;
        }
    }

    return avg_h;
}

double DhByDl(double *l1, double *l2, double *surfh)
{
    return -1.0 *
        (l1[2] * (surfh[1] - surfh[0]) + l1[1] * (surfh[0] - surfh[2]) +
        l1[0] * (surfh[2] - surfh[1])) /
        (l2[2] * (l1[1] - l1[0]) + l2[1] * (l1[0] - l1[2]) +
        l2[0] * (l1[2] - l1[1]));
}

double EffKh(const soil_struct *soil, double gw)
{
    double          k1, k2;
    double          d1, d2;

    gw = (gw > 0.0) ? gw : 0.0;

    if (gw > soil->depth - soil->dmac)
    {
        k1 = soil->kmach * soil->areafv + soil->ksath * (1.0 - soil->areafv);
        k2 = soil->ksath;

        if (gw > soil->depth)
        {
            d1 = soil->dmac;
            d2 = soil->depth - soil->dmac;
        }
        else
        {
            d1 = gw - (soil->depth - soil->dmac);
            d2 = soil->depth - soil->dmac;
        }

        return (k1 * d1 + k2 * d2) / (d1 + d2);
    }
    else
    {
        return soil->ksath;
    }
}

double OverLandFlow(double avg_h, double grad_h, double avg_sf, double crossa,
    double avg_rough)
{
    return crossa * pow(avg_h, 0.6666667) * grad_h / (sqrt(avg_sf) * avg_rough);
}

double SubFlowElemToElem(const elem_struct *elem, const elem_struct *nabr,
    int j)
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
    diff_h = (elem->ws.gw + elem->topo.zmin) - (nabr->ws.gw + nabr->topo.zmin);
    avg_h = AvgH(diff_h, elem->ws.gw, nabr->ws.gw);
    grad_h = diff_h / elem->topo.nabrdist[j];

    /* Take into account macropore effect */
    effk = EffKh(&elem->soil, elem->ws.gw);
    effk_nabr = EffKh(&nabr->soil, nabr->ws.gw);
    avg_ksat = 0.5 * (effk + effk_nabr);

    /* Groundwater flow modeled by Darcy's Law */
    return avg_ksat * grad_h * avg_h * elem->topo.edge[j];
}

double OvlFlowElemToElem(const elem_struct *elem, const elem_struct *nabr,
    int j, double avg_sf, int surf_mode)
{
    double          diff_h;
    double          avg_h;
    double          grad_h;
    double          avg_rough;
    double          crossa;

    diff_h = (surf_mode == KINEMATIC) ?
        elem->topo.zmax - nabr->topo.zmax :
        (elem->ws.surfh + elem->topo.zmax) - (nabr->ws.surfh + nabr->topo.zmax);
    avg_h = AvgHsurf(diff_h, elem->ws.surfh, nabr->ws.surfh);
    grad_h = diff_h / elem->topo.nabrdist[j];
    if (surf_mode == KINEMATIC)
    {
        avg_sf = (grad_h > 0.0) ? grad_h : GRADMIN;
    }
    else
    {
        avg_sf = (avg_sf > GRADMIN) ? avg_sf : GRADMIN;
    }
    /* Weighting needed */
    avg_rough = 0.5 * (elem->lc.rough + nabr->lc.rough);
    crossa = avg_h * elem->topo.edge[j];

    return OverLandFlow(avg_h, grad_h, avg_sf, crossa, avg_rough);
}

void BoundFluxElem(int bc_type, int j, const bc_struct *bc,
    const wstate_struct *ws, const topo_struct *topo, const soil_struct *soil,
    wflux_struct *wf)
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
        effk = EffKh(soil, ws->gw);
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
