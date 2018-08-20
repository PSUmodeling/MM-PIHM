#include "pihm.h"

void RiverFlow(elem_struct *elem, river_struct *river, int riv_mode)
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river_struct   *down;
        elem_struct    *left;
        elem_struct    *right;
        double          effk_nabr;
        double          effk;

        if (river[i].down > 0)
        {
            /*
             * Boundary conditions
             *
             * When a downstream segment is present, boundary conditions are
             * always applied to the upstream node
             */
            if (river[i].attrib.riverbc_type != 0)
            {
                river[i].wf.rivflow[UP_CHANL2CHANL] +=
                    BoundFluxRiver(river[i].attrib.riverbc_type, &river[i].ws,
                    &river[i].topo, &river[i].shp, &river[i].matl,
                    &river[i].bc);
            }

            down = &river[river[i].down - 1];

            /*
             * Channel flow between river-river segments
             */
            river[i].wf.rivflow[DOWN_CHANL2CHANL] =
                ChanFlowRiverToRiver(&river[i], down, riv_mode);

            /*
             * Subsurface flow between river-river segments
             */
            left = &elem[river[i].leftele - 1];
            right = &elem[river[i].rightele - 1];
            effk = 0.5 *
                (EffKh(left->ws.gw, &left->soil) +
                EffKh(right->ws.gw, &right->soil));
            left = &elem[down->leftele - 1];
            right = &elem[down->rightele - 1];
            effk_nabr = 0.5 *
                (EffKh(left->ws.gw, &left->soil) +
                EffKh(right->ws.gw, &right->soil));

            river[i].wf.rivflow[DOWN_AQUIF2AQUIF] =
                SubFlowRiverToRiver(&river[i], effk, down, effk_nabr);
        }
        else
        {
            /*
             * Outlet flux
             */
            river[i].wf.rivflow[DOWN_CHANL2CHANL] =
                OutletFlux(river[i].down, &river[i].ws, &river[i].topo,
                &river[i].shp, &river[i].matl, &river[i].bc);
            /* Note: boundary condition for subsurface element can be changed.
             * Assumption: no flow condition */
            river[i].wf.rivflow[DOWN_AQUIF2AQUIF] = 0.0;
        }

        /*
         * Flux between river segments and triangular elements
         */
        left = &elem[river[i].leftele - 1];
        right = &elem[river[i].rightele - 1];

        RiverToElem(&river[i], left, right);

        /*
         * Flux between river channel and subsurface
         */
        river[i].wf.rivflow[CHANL_LKG] =
            ChanLeak(&river[i].ws, &river[i].topo, &river[i].shp,
            &river[i].matl);
    }

    /*
     * Accumulate to get in-flow for down segments
     *
     * NOTE: Upstream flux summation must be calculated outside OMP to avoid
     * different threads accessing the same variable at the same time
     */
    for (i = 0; i < nriver; i++)
    {
        river_struct   *down;

        if (river[i].down > 0)
        {
            down = &river[river[i].down - 1];

            down->wf.rivflow[UP_CHANL2CHANL] -=
                river[i].wf.rivflow[DOWN_CHANL2CHANL];

            down->wf.rivflow[UP_AQUIF2AQUIF] -=
                river[i].wf.rivflow[DOWN_AQUIF2AQUIF];
        }
    }
}

void RiverToElem(river_struct *river, elem_struct *left, elem_struct *right)
{
    double          effk_left, effk_right;
    int             j;

    /* Lateral surface flux calculation between river-triangular element */
    if (river->leftele > 0)
    {
        river->wf.rivflow[LEFT_SURF2CHANL] = OvlFlowElemToRiver(left, river);
    }
    if (river->rightele > 0)
    {
        river->wf.rivflow[RIGHT_SURF2CHANL] = OvlFlowElemToRiver(right, river);
    }

    effk_left = EffKh(left->ws.gw, &left->soil);
    effk_right = EffKh(right->ws.gw, &right->soil);

    /* Lateral subsurface flux calculation between river-triangular element */
    if (river->leftele > 0)
    {
        river->wf.rivflow[LEFT_AQUIF2CHANL] =
            ChanFlowElemToRiver(left, effk_left, river, river->topo.dist_left);
    }
    if (river->rightele > 0)
    {
        river->wf.rivflow[RIGHT_AQUIF2CHANL] =
            ChanFlowElemToRiver(right, effk_right, river,
            river->topo.dist_right);
    }

    /* Lateral flux between rectangular element (beneath river) and triangular
     * element */
    if (river->leftele > 0)
    {
        river->wf.rivflow[LEFT_AQUIF2AQUIF] =
            SubFlowElemToRiver(left, effk_left, river,
            0.5 * (effk_left + effk_right), river->topo.dist_left);
    }
    if (river->rightele > 0)
    {
        river->wf.rivflow[RIGHT_AQUIF2AQUIF] =
            SubFlowElemToRiver(right, effk_right, river,
            0.5 * (effk_left + effk_right), river->topo.dist_right);
    }

    /* Replace flux term */
    /* Left */
    for (j = 0; j < NUM_EDGE; j++)
    {
        if (left->nabr[j] == -river->ind)
        {
            left->wf.ovlflow[j] = -river->wf.rivflow[LEFT_SURF2CHANL];
            left->wf.subsurf[j] = -(river->wf.rivflow[LEFT_AQUIF2CHANL] +
                river->wf.rivflow[LEFT_AQUIF2AQUIF]);
            break;
        }
    }

    /* Right */
    for (j = 0; j < NUM_EDGE; j++)
    {
        if (right->nabr[j] == -river->ind)
        {
            right->wf.ovlflow[j] = -river->wf.rivflow[RIGHT_SURF2CHANL];
            right->wf.subsurf[j] = -(river->wf.rivflow[RIGHT_AQUIF2CHANL] +
                river->wf.rivflow[RIGHT_AQUIF2AQUIF]);
            break;
        }
    }
}

double OvlFlowElemToRiver(const elem_struct *elem, const river_struct *river)
{
    double          zbank;
    double          flux;
    double          elem_h;
    double          rivseg_h;

    zbank = (river->topo.zmax > elem->topo.zmax) ?
        river->topo.zmax : elem->topo.zmax;

    elem_h = elem->topo.zmax + elem->ws.surfh;
    rivseg_h = river->topo.zbed + river->ws.stage;

    /*
     * Panday and Hyakorn 2004 AWR Eqs. (23) and (24)
     */
    if (rivseg_h > elem_h)
    {
        if (elem_h > zbank)
        {
            /* Submerged weir */
            flux = river->matl.cwr * 2.0 * sqrt(2.0 * GRAV) *
                river->shp.length * sqrt(rivseg_h - elem_h) *
                (rivseg_h - zbank) / 3.0;
        }
        else
        {
            if (zbank < rivseg_h)
            {
                /* Free-flowing weir */
                flux = river->matl.cwr * 2.0 * sqrt(2.0 * GRAV) *
                    river->shp.length * sqrt(rivseg_h - zbank) *
                    (rivseg_h - zbank) / 3.0;
            }
            else
            {
                flux = 0.0;
            }
        }
    }
    else if (elem->ws.surfh > DEPRSTG)
    {
        if (rivseg_h > zbank)
        {
            /* Submerged weir */
            flux = -river->matl.cwr * 2.0 * sqrt(2.0 * GRAV) *
                river->shp.length * sqrt(elem_h - rivseg_h) *
                (elem_h - zbank) / 3.0;
        }
        else
        {
            if (zbank < elem_h)
            {
                /* Free-flowing weir */
                flux = -river->matl.cwr * 2.0 * sqrt(2.0 * GRAV) *
                    river->shp.length * sqrt(elem_h - zbank) *
                    (elem_h - zbank) / 3.0;
            }
            else
            {
                flux = 0.0;
            }
        }
    }
    else
    {
        flux = 0.0;
    }

    return flux;
}

double ChanFlowRiverToRiver(const river_struct *river, const river_struct *down,
    int riv_mode)
{
    double          total_h;
    double          perim;
    double          total_h_down;
    double          perim_down;
    double          avg_perim;
    double          avg_rough;
    double          distance;
    double          diff_h;
    double          grad_h;
    double          avg_sf;
    double          crossa;
    double          crossa_down;
    double          avg_crossa;
    double          avg_h;

    total_h = river->ws.stage + river->topo.zbed;
    perim =
        RiverPerim(river->shp.intrpl_ord, river->ws.stage, river->shp.coeff);

    total_h_down = down->ws.stage + down->topo.zbed;
    perim_down =
        RiverPerim(down->shp.intrpl_ord, down->ws.stage, down->shp.coeff);

    avg_perim = (perim + perim_down) / 2.0;
    avg_rough = (river->matl.rough + down->matl.rough) / 2.0;
    distance = 0.5 * (river->shp.length + down->shp.length);

    diff_h = (riv_mode == KINEMATIC) ?
        (river->topo.zbed - down->topo.zbed) : (total_h - total_h_down);
    grad_h = diff_h / distance;
    avg_sf = (grad_h > 0.0) ? grad_h : RIVGRADMIN;
    crossa =
        RiverCroSectArea(river->shp.intrpl_ord, river->ws.stage,
            river->shp.coeff);
    crossa_down =
        RiverCroSectArea(down->shp.intrpl_ord, down->ws.stage, down->shp.coeff);
    avg_crossa = 0.5 * (crossa + crossa_down);
    avg_h = (avg_perim == 0.0) ? 0.0 : (avg_crossa / avg_perim);

    return OverLandFlow(avg_h, grad_h, avg_sf, crossa, avg_rough);
}

double SubFlowRiverToRiver(const river_struct *river, double effk,
    const river_struct *down, double effk_nabr)
{
    double          total_h;
    double          total_h_down;
    double          avg_wid;
    double          diff_h;
    double          avg_h;
    double          distance;
    double          grad_h;
    double          aquifer_depth;
    double          avg_ksat;

    /* Lateral flux calculation between element beneath river (ebr) * and ebr */
    total_h = river->ws.gw + river->topo.zmin;
    total_h_down = down->ws.gw + down->topo.zmin;
    avg_wid = (river->shp.width + down->shp.width) / 2.0;
    diff_h = total_h - total_h_down;
    avg_h = AvgH(diff_h, river->ws.gw, down->ws.gw);
    distance = 0.5 * (river->shp.length + down->shp.length);
    grad_h = diff_h / distance;
    aquifer_depth = river->topo.zbed - river->topo.zmin;
#if defined(_ARITH_)
    avg_ksat = 0.5 * (effk + effk_nabr);
#else
    avg_ksat = 2.0 / (1.0 / effk + 1.0 / effk_nabr);
#endif
    /* Groundwater flow modeled by Darcy's law */
    return avg_ksat * grad_h * avg_h * avg_wid;
}

double OutletFlux(int down, const river_wstate_struct *ws,
    const river_topo_struct *topo, const shp_struct *shp,
    const matl_struct *matl, const river_bc_struct *bc)
{
    double          total_h;
    double          total_h_down;
    double          distance;
    double          grad_h;
    double          avg_h;
    double          avg_perim;
    double          crossa;
    double          discharge = 0.0;

    switch (down)
    {
        case DIRICHLET:
            /* Dirichlet boundary condition */
            total_h = ws->stage + topo->zbed;
            total_h_down = bc->head;
            distance = 0.5 * shp->length;
            grad_h = (total_h - total_h_down) / distance;
            avg_h = AvgH(grad_h, ws->stage,
                ((bc->head - (topo->node_zmax - shp->depth) > 0.0) ?
                bc->head - (topo->node_zmax - shp->depth) : 0.0));
            avg_perim = RiverPerim(shp->intrpl_ord, ws->stage, shp->coeff);
            crossa = RiverCroSectArea(shp->intrpl_ord, ws->stage, shp->coeff);
            avg_h = (avg_perim == 0.0) ? 0.0 : (crossa / avg_perim);
            discharge =
                OverLandFlow(avg_h, grad_h, grad_h, crossa, matl->rough);
            break;
        case NEUMANN:
            /* Neumann boundary condition */
            discharge = -bc->flux;
            break;
        case ZERO_DPTH_GRAD:
            /* Zero-depth-gradient boundary conditions */
            distance = 0.5 * shp->length;
            grad_h = (topo->zbed - (topo->node_zmax - shp->depth)) / distance;
            avg_h = ws->stage;
            avg_perim = RiverPerim(shp->intrpl_ord, ws->stage, shp->coeff);
            crossa = RiverCroSectArea(shp->intrpl_ord, ws->stage, shp->coeff);
            discharge = sqrt(grad_h) * crossa * ((avg_perim > 0.0) ?
                pow(crossa / avg_perim, 2.0 / 3.0) : 0.0) / matl->rough;
            break;
        case CRIT_DPTH:
            /* Critical depth boundary conditions */
            crossa = RiverCroSectArea(shp->intrpl_ord, ws->stage, shp->coeff);
            discharge = crossa * sqrt(GRAV * ws->stage);
            break;
        default:
            PIHMprintf(VL_ERROR,
                "Error: River routing boundary condition type (%d) "
                "is not recognized.\n", down);
            PIHMexit(EXIT_FAILURE);
    }

    return discharge;
}

double BoundFluxRiver(int riverbc_type, const river_wstate_struct *ws,
    const river_topo_struct *topo, const shp_struct *shp,
    const matl_struct *matl, const river_bc_struct *bc)
{
    double          total_h;
    double          total_h_down;
    double          distance;
    double          grad_h;
    double          avg_h;
    double          avg_perim;
    double          crossa;
    double          flux = 0.0;

    if (riverbc_type > 0)
    {
        /* Dirichlet boundary condition */
        total_h = ws->stage + topo->zbed;
        total_h_down = bc->head;
        distance = 0.5 * shp->length;
        grad_h = (total_h - total_h_down) / distance;
        avg_h = AvgH(grad_h, ws->stage,
            ((bc->head - topo->zbed > 0.0) ? (bc->head - topo->zbed) : 0.0));
        avg_perim = RiverPerim(shp->intrpl_ord, ws->stage, shp->coeff);
        crossa = RiverCroSectArea(shp->intrpl_ord, ws->stage, shp->coeff);
        avg_h = (avg_perim == 0.0) ? 0.0 : (crossa / avg_perim);
        flux =
            OverLandFlow(avg_h, grad_h, grad_h, crossa, matl->rough);
    }
    else if (riverbc_type < 0)
    {
        /* Neumann boundary condition */
        flux = -bc->flux;
    }

    return flux;
}

double ChanFlowElemToRiver(const elem_struct *elem, double effk,
    const river_struct *river, double distance)
{
    double          diff_h;
    double          avg_h;
    double          grad_h;
    double          avg_ksat;

    diff_h = (river->ws.stage + river->topo.zbed) -
        (elem->ws.gw + elem->topo.zmin);

    /* This is head in neighboring cell representation */
    if (elem->topo.zmin > river->topo.zbed)
    {
        avg_h = elem->ws.gw;
    }
    else if (elem->topo.zmin + elem->ws.gw > river->topo.zbed)
    {
        avg_h = elem->topo.zmin + elem->ws.gw - river->topo.zbed;
    }
    else
    {
        avg_h = 0.0;
    }
    avg_h = AvgH(diff_h, river->ws.stage, avg_h);

    grad_h = diff_h / distance;

    avg_ksat = 0.5 * (effk + river->matl.ksath);

    return  river->shp.length * avg_ksat * grad_h * avg_h;
}

double SubFlowElemToRiver(const elem_struct *elem, double effk,
    const river_struct *river, double effk_riv, double distance)
{
    double          diff_h;
    double          avg_h;
    double          aquifer_depth;
    double          avg_ksat;
    double          grad_h;

    diff_h = (river->ws.gw + river->topo.zmin) -
        (elem->ws.gw + elem->topo.zmin);

    /* This is head in neighboring cell represention */
    if (elem->topo.zmin > river->topo.zbed)
    {
        avg_h = 0.0;
    }
    else if (elem->topo.zmin + elem->ws.gw > river->topo.zbed)
    {
        avg_h = river->topo.zbed - elem->topo.zmin;
    }
    else
    {
        avg_h = elem->ws.gw;
    }
    avg_h = AvgH(diff_h, river->ws.gw, avg_h);
    aquifer_depth = river->topo.zbed - river->topo.zmin;

#if defined(_ARITH_)
    avg_ksat = 0.5 * (effk + effk_riv);
#else
    avg_ksat = 2.0 / (1.0 / effk + 1.0 / effk_riv);
#endif
    grad_h = diff_h / distance;

    return river->shp.length * avg_ksat * grad_h * avg_h;
}

double ChanLeak(const river_wstate_struct *ws, const river_topo_struct *topo,
    const shp_struct *shp, const matl_struct *matl)
{
    double          diff_h;
    double          grad_h;

    if (topo->zbed - (ws->gw + topo->zmin) > 0.0)
    {
        diff_h = ws->stage;
    }
    else
    {
        diff_h = ws->stage + topo->zbed - (ws->gw + topo->zmin);
    }

    grad_h = diff_h / matl->bedthick;

    return matl->ksatv * shp->width * shp->length * grad_h;
}

double RiverCroSectArea(int order, double depth, double coeff)
{
    double          cs_area = 0.0;

    depth = (depth > 0.0) ? depth : 0.0;

    switch (order)
    {
        case RECTANGLE:
            cs_area = depth * coeff;
            break;
        case TRIANGLE:
            cs_area = depth * depth / coeff;
            break;
        case QUADRATIC:
            cs_area = 4.0 * depth * sqrt(depth) / (3.0 * sqrt(coeff));
            break;
        case CUBIC:
            cs_area =
                3.0 * pow(depth, 4.0 / 3.0) / (2.0 * pow(coeff, 1.0 / 3.0));
            break;
        default:
            PIHMprintf(VL_ERROR, "Error: River order %d is not defined.\n",
                order);
            PIHMexit(EXIT_FAILURE);
    }

    return cs_area;
}

double RiverPerim(int order, double depth, double coeff)
{
    double          perim = 0.0;

    depth = (depth > 0.0) ? depth : 0.0;

    switch (order)
    {
        case RECTANGLE:
            perim = 2.0 * depth + coeff;
            break;
        case TRIANGLE:
            perim = 2.0 * depth * sqrt(1.0 + coeff * coeff) / coeff;
            break;
        case QUADRATIC:
            perim = sqrt(depth * (1.0 + 4.0 * coeff * depth) / coeff) +
                log(2.0 * sqrt(coeff * depth) +
                sqrt(1.0 + 4.0 * coeff * depth)) / (2.0 * coeff);
            break;
        case CUBIC:
            perim = 2.0 * ((pow(depth * (1.0 + 9.0 * pow(coeff, 2.0 / 3.0) *
                depth), 0.5) / 3.0) +
                (log(3.0 * pow(coeff, 1.0 / 3.0) * sqrt(depth) +
                pow(1.0 + 9.0 * pow(coeff, 2.0 / 3.0) * depth, 0.5)) /
                (9.0 * pow(coeff, 1.0 / 3.0))));
            break;
        default:
            PIHMprintf(VL_ERROR, "Error: River order %d is not defined.\n",
                order);
            PIHMexit(EXIT_FAILURE);
    }

    return perim;
}
