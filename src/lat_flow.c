#include "pihm.h"

void LateralFlow (pihm_struct pihm)
{
    int             i, j;
    double          dif_y_sub;
    double          avg_y_sub;
    double          distance;
    double          grad_y_sub;
    double          surfh[3];
    double         *dhbydx;
    double         *dhbydy;
    double          effk;
    double          effk_nabr;
    double          avg_ksat;
    double          dif_y_surf;
    double          avg_y_surf;
    double          grad_y_surf;
    double          avg_sf;
    double          avg_rough;
    double          crossa;

    elem_struct    *elem;
    elem_struct    *nabr;
    river_struct   *riv;

    dhbydx = (double *)malloc (pihm->numele * sizeof (double));
    dhbydy = (double *)malloc (pihm->numele * sizeof (double));

    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        if (pihm->ctrl.surf_mode == DIFF_WAVE)
        {
            for (j = 0; j < 3; j++)
            {
                if (elem->nabr[j] > 0)
                {
                    nabr = &pihm->elem[elem->nabr[j] - 1];
                    surfh[j] = nabr->topo.zmax + nabr->ws.surf;
                }
                else if (elem->nabr[j] < 0)
                {
                    riv = &pihm->riv[-elem->nabr[j] - 1];

                    if (riv->ws.stage > riv->shp.depth)
                    {
                        surfh[j] = riv->topo.zbed + riv->ws.stage;
                    }
                    else
                    {
                        surfh[j] = riv->topo.zmax;
                    }
                }
                else
                {
                    if (elem->attrib.bc_type[j] == 0)
                    {
                        surfh[j] = elem->topo.zmax + elem->ws.surf;
                    }
                    else
                    {
                        surfh[j] = elem->bc.head[j];
                    }
                }
            }

            dhbydx[i] = DhByDl (elem->topo.surfy, elem->topo.surfx, surfh);
            dhbydy[i] = DhByDl (elem->topo.surfx, elem->topo.surfy, surfh);
        }
    }

    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        for (j = 0; j < 3; j++)
        {
            if (elem->nabr[j] != 0)
            {
                if (elem->nabr[j] > 0)
                {
                    nabr = &pihm->elem[elem->nabr[j] - 1];
                }
                else
                {
                    riv = &pihm->riv[-elem->nabr[j] - 1];
                    nabr = (i == riv->leftele - 1) ?
                        &pihm->elem[riv->rightele - 1] :
                        &pihm->elem[riv->leftele - 1];
                }

                /*
                 * Subsurface lateral flux calculation between triangular
                 * elements
                 */
                dif_y_sub =
                    (elem->ws.gw + elem->topo.zmin) - (nabr->ws.gw +
                    nabr->topo.zmin);
                avg_y_sub = AvgY (dif_y_sub, elem->ws.gw, nabr->ws.gw);
                distance =
                    sqrt (pow (elem->topo.x - nabr->topo.x,
                        2) + pow (elem->topo.y - nabr->topo.y, 2));
                grad_y_sub = dif_y_sub / distance;
                /* Take into account macropore effect */
                effk =
                    EffKH (elem->ws.gw, elem->soil.depth, elem->soil.dmac,
                    elem->soil.kmach, elem->soil.areafv, elem->soil.ksath);
                effk_nabr =
                    EffKH (nabr->ws.gw, nabr->soil.depth, nabr->soil.dmac,
                    nabr->soil.kmach, nabr->soil.areafv, nabr->soil.ksath);
                avg_ksat = 0.5 * (effk + effk_nabr);
                /* Groundwater flow modeled by Darcy's Law */
                elem->wf.fluxsub[j] =
                    avg_ksat * grad_y_sub * avg_y_sub * elem->topo.edge[j];

                /* 
                 * Surface lateral flux calculation between triangular
                 * elements
                 */
                if (pihm->ctrl.surf_mode == KINEMATIC)
                {
                    dif_y_surf = elem->topo.zmax - nabr->topo.zmax;
                }
                else
                {
                    dif_y_surf =
                        (elem->ws.surf + elem->topo.zmax) - (nabr->ws.surf +
                        nabr->topo.zmax);
                }
                avg_y_surf = AvgY (dif_y_surf, elem->ws.surf, nabr->ws.surf);
                grad_y_surf = dif_y_surf / distance;
                avg_sf =
                    0.5 * (sqrt (pow (dhbydx[i], 2) + pow (dhbydy[i],
                            2)) + sqrt (pow (dhbydx[nabr->ind - 1],
                            2) + pow (dhbydy[nabr->ind - 1], 2)));
                if (pihm->ctrl.surf_mode == KINEMATIC)
                {
                    avg_sf = (grad_y_surf > 0.0) ? grad_y_surf : GRADMIN;
                }
                else
                {
                    avg_sf = (avg_sf > GRADMIN) ? avg_sf : GRADMIN;
                }
                /* Weighting needed */
                avg_rough = 0.5 * (elem->lc.rough + nabr->lc.rough);
                crossa = avg_y_surf * elem->topo.edge[j];
                elem->wf.fluxsurf[j] =
                    OverlandFlow (avg_y_surf, grad_y_surf, avg_sf, crossa,
                    avg_rough);
            }
            else                /* Boundary condition flux */
            {
                /* No flow (natural) boundary condition is default */
                if (elem->attrib.bc_type[j] == 0)
                {
                    elem->wf.fluxsurf[j] = 0.0;
                    elem->wf.fluxsub[j] = 0.0;
                }
                /* Note: ideally different boundary conditions need to be
                 * incorporated for surf and subsurf respectively */
                else if (elem->attrib.bc_type[j] > 0)
                {
                    /* Note: the formulation assumes only Dirichlet TS right
                     * now */
                    /* note the assumption here is no flow for surface */
                    elem->wf.fluxsurf[j] = 0.0;
                    dif_y_sub =
                        elem->ws.gw + elem->topo.zmin - elem->bc.head[j];
                    avg_y_sub =
                        AvgY (dif_y_sub, elem->ws.gw,
                        elem->bc.head[j] - elem->topo.zmin);
                    /* Minimum distance from circumcenter to the edge of the
                     * triangle on which boundary condition is defined */
                    distance =
                        sqrt (pow (elem->topo.edge[0] * elem->topo.edge[1] *
                            elem->topo.edge[2] / (4.0 * elem->topo.area),
                            2) - pow (elem->topo.edge[j] / 2.0, 2));
                    effk =
                        EffKH (elem->ws.gw, elem->soil.depth, elem->soil.dmac,
                        elem->soil.kmach, elem->soil.areafv,
                        elem->soil.ksath);
                    avg_ksat = effk;
                    grad_y_sub = dif_y_sub / distance;
                    elem->wf.fluxsub[j] =
                        avg_ksat * grad_y_sub * avg_y_sub *
                        elem->topo.edge[j];
                }
                else
                {
                    /* Neumann bc (note: md->ele[i].bc[j] value has to be
                     * = 2+(index of neumann boundary ts) */
                    elem->wf.fluxsurf[j] = elem->bc.flux[j];
                    elem->wf.fluxsub[j] = elem->bc.flux[j];
                }
            }                   /* End of specified boundary condition */
        }                       /* End of neighbor loop */
    }                           /* End of element loop */

    free (dhbydx);
    free (dhbydy);
}

double AvgYsfc (double diff, double yi, double yinabr)
{
    double          avg_y;

    if (diff > 0.0)
    {
        if (yi > DEPRSTG)
        {
            avg_y = 1.0 * yi;
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
            avg_y = 1.0 * yinabr;
        }
        else
        {
            avg_y = 0.0;
        }
    }

    return (avg_y);
}

double AvgY (double diff, double yi, double yinabr)
{
    double          avg_y;

    if (diff > 0.0)
    {
        if (yi > 0.0)
        {
            avg_y = 1.0 * yi;
        }
        else
        {
            avg_y = 0.0;
        }
    }
    else
    {
        if (yinabr > 0.0)
        {
            avg_y = 1.0 * yinabr;
        }
        else
        {
            avg_y = 0.0;
        }
    }

    return (avg_y);
}


double DhByDl (double *l1, double *l2, double *surfh)
{
    return (-1.0 * (l1[2] * (surfh[1] - surfh[0]) + l1[1] * (surfh[0] -
                surfh[2]) + l1[0] * (surfh[2] - surfh[1])) / (l2[2] * (l1[1] -
                l1[0]) + l2[1] * (l1[0] - l1[2]) + l2[0] * (l1[2] - l1[1])));
}

double EffKH (double tmpy, double aqdepth, double macd, double macksath,
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
                ksath * (aqdepth - macd + (tmpy - (aqdepth -
                            macd)) * (1.0 - areaf))) / tmpy;
        }
    }
    else
    {
        return (ksath);
    }
}


double OverlandFlow (double avg_y, double grad_y, double avg_sf,
    double crossa, double avg_rough)
{
    return (crossa * pow (avg_y,
            2.0 / 3.0) * grad_y / (sqrt (fabs (avg_sf)) * avg_rough));
}
