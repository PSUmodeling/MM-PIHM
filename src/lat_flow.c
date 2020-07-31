#include "pihm.h"

void LateralFlow(int surf_mode, const river_struct river[], elem_struct elem[])
{
    int             i;
    double         *dh_dx;
    double         *dh_dy;

    if (surf_mode == DIFF_WAVE)
    {
        dh_dx = (double *)malloc(nelem * sizeof(double));
        dh_dy = (double *)malloc(nelem * sizeof(double));

        FrictionSlope(elem, river, dh_dx, dh_dy);
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;
        double          avg_sf;
        elem_struct    *nabr;

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].nabr[j] == 0)       /* Boundary condition flux */
            {
                BoundFluxElem(elem[i].attrib.bc_type[j], j, &elem[i].topo,
                    &elem[i].soil, &elem[i].bc, &elem[i].ws, &elem[i].wf);
            }
            else
            {
                nabr = &elem[elem[i].nabr[j] - 1];

                /* Subsurface flow between triangular elements */
                elem[i].wf.subsurf[j] = SubsurfFlow(j, &elem[i], nabr);

                /* Surface flow between triangular elements */
                if (elem[i].nabr_river[j] == 0)
                {
                    /* Surface flux between triangular elements */
                    /* avg_sf not needed in kinematic mode */
                    avg_sf = (surf_mode == DIFF_WAVE) ?
                        0.5 * (sqrt(dh_dx[i] * dh_dx[i] + dh_dy[i] * dh_dy[i]) +
                        sqrt(dh_dx[nabr->ind - 1] * dh_dx[nabr->ind - 1] +
                        dh_dy[nabr->ind - 1] * dh_dy[nabr->ind - 1])) : 0.0;

                    elem[i].wf.ovlflow[j] = OvlFlowElemToElem(j, avg_sf,
                        surf_mode, &elem[i], nabr);
                }
            }
        }    /* End of neighbor loop */
    }    /* End of element loop */

    if (surf_mode == DIFF_WAVE)
    {
        free(dh_dx);
        free(dh_dy);
    }

#if defined(_FBR_)
    /*
     * Lateral fractured bedrock flow
     */
# if defined(_OPENMP)
#  pragma omp parallel for
# endif
    for (i = 0; i < nelem; i++)
    {
        int             j;
        elem_struct    *nabr;

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].nabr[j] == 0)
            {
                elem[i].wf.fbrflow[j] =
                    DeepBoundFluxElem(elem[i].attrib.fbrbc_type[j], j,
                    &elem[i].topo, &elem[i].geol, &elem[i].fbr_bc, &elem[i].ws);
            }
            else
            {
                nabr = &elem[elem[i].nabr[j] - 1];

                /* Groundwater flow modeled by Darcy's Law */
                elem[i].wf.fbrflow[j] =
                    DeepFlowElemToElem(elem[i].topo.nabrdist[j],
                        elem[i].topo.edge[j], &elem[i], nabr);
            }
        }
    }
#endif
}

void FrictionSlope(const elem_struct elem[], const river_struct river[],
    double dh_dx[], double dh_dy[])
{
    int             i;
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;
        double          surfh[NUM_EDGE];
        double          x[NUM_EDGE];
        double          y[NUM_EDGE];
        const elem_struct *nabr;
        const river_struct *river_nabr;

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].nabr[j] == 0)
            {
                surfh[j] = elem[i].topo.zmax + elem[i].ws.surfh;
                x[j] = elem[i].topo.nabr_x[j];
                y[j] = elem[i].topo.nabr_y[j];
            }
            else if (elem[i].nabr_river[j] == 0)
            {
                nabr = &elem[elem[i].nabr[j] - 1];

                surfh[j] = nabr->topo.zmax + nabr->ws.surfh;
                x[j] = elem[i].topo.nabr_x[j];
                y[j] = elem[i].topo.nabr_y[j];
            }
            else
            {
                river_nabr = &river[elem[i].nabr_river[j] - 1];

                surfh[j] = (river_nabr->ws.stage > river_nabr->shp.depth) ?
                    river_nabr->topo.zbed + river_nabr->ws.stage :
                    river_nabr->topo.zmax;
                x[j] = river[elem[i].nabr_river[j] - 1].topo.x;
                y[j] = river[elem[i].nabr_river[j] - 1].topo.y;
            }
        }

        dh_dx[i] = DhByDl(y, x, surfh);
        dh_dy[i] = DhByDl(x, y, surfh);
    }
}

double AvgHsurf(double diff, double h_surf, double h_nabr)
{
    return (diff > 0.0) ?
        MAX(h_surf - DEPRSTG, 0.0) : MAX(h_nabr- DEPRSTG, 0.0);
}

double AvgH(double diff, double h_sub, double h_nabr)
{
    return (diff > 0.0) ? MAX(h_sub, 0.0) : MAX(h_nabr, 0.0);
}

double DhByDl(const double l1[], const double l2[], const double surfh[])
{
    return -(l1[2] * (surfh[1] - surfh[0]) + l1[1] * (surfh[0] - surfh[2]) +
        l1[0] * (surfh[2] - surfh[1])) /
        (l2[2] * (l1[1] - l1[0]) + l2[1] * (l1[0] - l1[2]) +
        l2[0] * (l1[2] - l1[1]));
}

double EffKh(double gw, const soil_struct *soil)
{
    double          k1, k2;
    double          d1, d2;

    gw = MAX(gw, 0.0);

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

double OverLandFlow(double avg_h, double grad_h, double avg_sf,
    double cross_area, double avg_rough)
{
    return cross_area * pow(avg_h, 0.6666667) * grad_h /
        (sqrt(avg_sf) * avg_rough);
}

double SubsurfFlow(int j, const elem_struct *elem_ptr, const elem_struct *nabr)
{
    double          diff_h;
    double          avg_h;
    double          grad_h;
    double          effk, effk_nabr;
    double          avg_ksat;

    /*
     * Subsurface lateral flux calculation between triangular elements
     */
    diff_h = (elem_ptr->ws.gw + elem_ptr->topo.zmin) -
        (nabr->ws.gw + nabr->topo.zmin);
    avg_h = AvgH(diff_h, elem_ptr->ws.gw, nabr->ws.gw);
    grad_h = diff_h / elem_ptr->topo.nabrdist[j];

    /* Take into account macropore effect */
    effk = EffKh(elem_ptr->ws.gw, &elem_ptr->soil);
    effk_nabr = EffKh(nabr->ws.gw, &nabr->soil);
    avg_ksat = 0.5 * (effk + effk_nabr);

    /* Groundwater flow modeled by Darcy's Law */
    return avg_ksat * grad_h * avg_h * elem_ptr->topo.edge[j];
}

double OvlFlowElemToElem(int j, double avg_sf, int surf_mode,
    const elem_struct *elem_ptr, const elem_struct *nabr)
{
    double          diff_h;
    double          avg_h;
    double          grad_h;
    double          avg_rough;
    double          cross_area;

    diff_h = (surf_mode == KINEMATIC) ?
        elem_ptr->topo.zmax - nabr->topo.zmax :
        (elem_ptr->ws.surfh + elem_ptr->topo.zmax) -
        (nabr->ws.surfh + nabr->topo.zmax);
    avg_h = AvgHsurf(diff_h, elem_ptr->ws.surfh, nabr->ws.surfh);
    grad_h = MAX(diff_h / elem_ptr->topo.nabrdist[j], GRADMIN);
    avg_sf = (surf_mode == KINEMATIC) ? grad_h : MAX(avg_sf, GRADMIN);
    avg_rough = 0.5 * (elem_ptr->lc.rough + nabr->lc.rough);
    cross_area = avg_h * elem_ptr->topo.edge[j];

    return OverLandFlow(avg_h, grad_h, avg_sf, cross_area, avg_rough);
}

void BoundFluxElem(int bc_type, int j, const topo_struct *topo,
    const soil_struct *soil, const bc_struct *bc, const wstate_struct *ws,
    wflux_struct *wf)
{
    double          diff_h;
    double          avg_h;
    double          effk;
    double          avg_ksat;
    double          grad_h;

    /* Assume no flow for surface */
    wf->ovlflow[j] = 0.0;

    /* No flow (natural) boundary condition is default */
    if (bc_type == NO_FLOW)
    {
        wf->subsurf[j] = 0.0;
    }
    /* Note: ideally different boundary conditions need to be
     * incorporated for surf and subsurf respectively */
    else if (bc_type > 0)
    {
        diff_h = ws->gw + topo->zmin - bc->head[j];
        avg_h = AvgH(diff_h, ws->gw, bc->head[j] - topo->zmin);
        /* Minimum distance from circumcenter to the edge of the triangle
         * on which boundary condition is defined */
        effk = EffKh(ws->gw, soil);
        avg_ksat = effk;
        grad_h = diff_h / topo->nabrdist[j];
        wf->subsurf[j] = avg_ksat * grad_h * avg_h * topo->edge[j];
    }
    else
    {
        /* Negative sign is added so the positive numbers in forcing time series
         * represents source */
        wf->subsurf[j] = -bc->flux[j];
    }
}

#if defined(_FBR_)
double DeepFlowElemToElem(double distance, double edge,
    const elem_struct *elem_ptr, const elem_struct *nabr)
{
    double          diff_h;
    double          avg_h;
    double          grad_h;
    double          effk, effk_nabr;
    double          avg_ksat;

    diff_h = (elem_ptr->ws.fbr_gw + elem_ptr->topo.zbed) -
        (nabr->ws.fbr_gw + nabr->topo.zbed);
    avg_h = AvgH(diff_h, elem_ptr->ws.fbr_gw, nabr->ws.fbr_gw);
    grad_h = diff_h / distance;

    effk = EffKh(elem_ptr->ws.fbr_gw, &elem_ptr->geol);
    effk_nabr = EffKh(nabr->ws.fbr_gw, &nabr->geol);
    avg_ksat = 0.5 * (effk + effk_nabr);

    return avg_ksat * grad_h * avg_h * edge;
}

double DeepBoundFluxElem(int bc_type, int j, const topo_struct *topo,
    const soil_struct *geol, const bc_struct *bc, const wstate_struct *ws)
{
    double          diff_h;
    double          avg_h;
    double          effk;
    double          grad_h;
    double          flux;

    /* No flow (natural) boundary condition is default */
    if (bc_type == NO_FLOW)
    {
        flux = 0.0;
    }
    else if (bc_type > 0)
    {
        /* Dirichlet boundary conditions */
        diff_h = ws->fbr_gw + topo->zbed - bc->head[j];
        avg_h = AvgH(diff_h, ws->fbr_gw, bc->head[j] - topo->zbed);
        /* Minimum distance from circumcenter to the edge of the triangle
         * on which boundary condition is defined */
        effk = geol->ksath;
        grad_h = diff_h / topo->nabrdist[j];
        flux = effk * grad_h * avg_h * topo->edge[j];
    }
    else
    {
        /* Neumann boundary conditions */
        flux = -bc->flux[j];
    }

    return flux;
}
#endif
