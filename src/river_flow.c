#include "pihm.h"

void RiverFlow(elem_struct *elem, river_struct *rivseg, int riv_mode)
{
    int             i;

#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river_struct   *down;
        elem_struct    *left;
        elem_struct    *right;
        double          dif_y;
        double          grad_y;
        double          effk_nabr;
        double          effk;

        if (rivseg[i].down > 0)
        {
            down = &rivseg[rivseg[i].down - 1];

            /*
             * Channel flow between river-river segments
             */
            rivseg[i].wf.rivflow[DOWN_CHANL2CHANL] =
                ChanFlowRiverToRiver(&rivseg[i].ws, &rivseg[i].topo,
                &rivseg[i].shp, &rivseg[i].matl, &down->ws, &down->topo,
                &down->shp, &down->matl, riv_mode);
            /* Accumulate to get in-flow for down segments */
            down->wf.rivflow[UP_CHANL2CHANL] -=
                rivseg[i].wf.rivflow[DOWN_CHANL2CHANL];

            /*
             * Subsurface flow between river-river segments
             */
            left = &elem[rivseg[i].leftele - 1];
            right = &elem[rivseg[i].rightele - 1];
            effk = 0.5 *
                (EffKh(left->ws.gw, left->soil.depth, left->soil.dmac,
                left->soil.kmach, left->soil.areafv, left->soil.ksath) +
                EffKh(right->ws.gw, right->soil.depth, right->soil.dmac,
                right->soil.kmach, left->soil.areafv, right->soil.ksath));
            left = &elem[down->leftele - 1];
            right = &elem[down->rightele - 1];
            effk_nabr = 0.5 *
                (EffKh(left->ws.gw, left->soil.depth, left->soil.dmac,
                left->soil.kmach, left->soil.areafv, left->soil.ksath) +
                EffKh(right->ws.gw, right->soil.depth, right->soil.dmac,
                right->soil.kmach, left->soil.areafv, right->soil.ksath));

            rivseg[i].wf.rivflow[DOWN_AQUIF2AQUIF] =
                SubFlowRiverToRiver(&rivseg[i].ws, &rivseg[i].topo,
                &rivseg[i].shp, effk, &down->ws, &down->topo, &down->shp,
                effk_nabr);
            /* Accumulate to get in-flow for down segments */
            down->wf.rivflow[UP_AQUIF2AQUIF] -=
                rivseg[i].wf.rivflow[DOWN_AQUIF2AQUIF];
        }
        else
        {
            /*
             * Outlet flux
             */
            /* Note: boundary condition for subsurface element can be changed.
             * Assumption: no flow condition */
            rivseg[i].wf.rivflow[DOWN_CHANL2CHANL] =
                OutletFlux(rivseg[i].down, &rivseg[i].ws, &rivseg[i].topo,
                &rivseg[i].shp, &rivseg[i].matl, &rivseg[i].bc);

            rivseg[i].wf.rivflow[DOWN_AQUIF2AQUIF] = 0.0;
        }

        /*
         * Flux between river segments and triangular elements
         */
        left = &elem[rivseg[i].leftele - 1];
        right = &elem[rivseg[i].rightele - 1];

        if (rivseg[i].leftele > 0)
        {
            RiverToElem(&rivseg[i], left, right, i + 1,
                rivseg[i].topo.dist_left,
                &rivseg[i].wf.rivflow[LEFT_SURF2CHANL],
                &rivseg[i].wf.rivflow[LEFT_AQUIF2CHANL],
                &rivseg[i].wf.rivflow[LEFT_AQUIF2AQUIF]);
        }

        if (rivseg[i].rightele > 0)
        {
            RiverToElem(&rivseg[i], right, left, i + 1,
                rivseg[i].topo.dist_right,
                &rivseg[i].wf.rivflow[RIGHT_SURF2CHANL],
                &rivseg[i].wf.rivflow[RIGHT_AQUIF2CHANL],
                &rivseg[i].wf.rivflow[RIGHT_AQUIF2AQUIF]);
        }

        if (rivseg[i].topo.zbed - (rivseg[i].ws.gw + rivseg[i].topo.zmin) > 0.0)
        {
            dif_y = rivseg[i].ws.stage;
        }
        else
        {
            dif_y = rivseg[i].ws.stage + rivseg[i].topo.zbed -
                (rivseg[i].ws.gw + rivseg[i].topo.zmin);
        }
        grad_y = dif_y / rivseg[i].matl.bedthick;
        rivseg[i].wf.rivflow[CHANL_LKG] = rivseg[i].matl.ksatv *
            rivseg[i].shp.width * rivseg[i].shp.length * grad_y;
    }
}

void RiverToElem(river_struct *rivseg, elem_struct *elem, elem_struct *oppbank,
    int ind, double distance, double *fluxsurf, double *fluxriv,
    double *fluxsub)
{
    double          dif_y_sub;
    double          avg_y_sub;
    double          effk;
    double          grad_y_sub;
    double          effk_nabr;
    double          avg_ksat;
    double          aquifer_depth;
    double          zbank;
    int             j;

    /* Lateral surface flux calculation between river-triangular element */
    zbank = (rivseg->topo.zmax > elem->topo.zmax) ?
        rivseg->topo.zmax : elem->topo.zmax;
    *fluxsurf = OvlFlowElemToRiver(elem->topo.zmax, elem->ws.surfh, zbank,
        rivseg->ws.stage, rivseg->topo.zbed, rivseg->matl.cwr,
        rivseg->shp.length);

    /* Lateral subsurface flux calculation between river-triangular element */
    dif_y_sub =
        (rivseg->ws.stage + rivseg->topo.zbed) -
        (elem->ws.gw + elem->topo.zmin);
    /* This is head in neighboring cell representation */
    if (elem->topo.zmin > rivseg->topo.zbed)
    {
        avg_y_sub = elem->ws.gw;
    }
    else if (elem->topo.zmin + elem->ws.gw > rivseg->topo.zbed)
    {
        avg_y_sub = elem->topo.zmin + elem->ws.gw - rivseg->topo.zbed;
    }
    else
    {
        avg_y_sub = 0.0;
    }
    avg_y_sub = AvgH(dif_y_sub, rivseg->ws.stage, avg_y_sub);
    effk = rivseg->matl.ksath;
    grad_y_sub = dif_y_sub / distance;
    /* Take into account macropore effect */
    effk_nabr = EffKh(elem->ws.gw, elem->soil.depth, elem->soil.dmac,
        elem->soil.kmach, elem->soil.areafv, elem->soil.ksath);
    avg_ksat = 0.5 * (effk + effk_nabr);
    *fluxriv = rivseg->shp.length * avg_ksat * grad_y_sub * avg_y_sub;

    /* Lateral flux between rectangular element (beneath river) and triangular
     * element */
    dif_y_sub =
        (rivseg->ws.gw + rivseg->topo.zmin) - (elem->ws.gw + elem->topo.zmin);
    /* this is head in neighboring cell represention */
    if (elem->topo.zmin > rivseg->topo.zbed)
    {
        avg_y_sub = 0.0;
    }
    else if (elem->topo.zmin + elem->ws.gw > rivseg->topo.zbed)
    {
        avg_y_sub = rivseg->topo.zbed - elem->topo.zmin;
    }
    else
    {
        avg_y_sub = elem->ws.gw;
    }
    avg_y_sub = AvgH(dif_y_sub, rivseg->ws.gw, avg_y_sub);
    aquifer_depth = rivseg->topo.zbed - rivseg->topo.zmin;
    effk = 0.5 *
        (EffKh(elem->ws.gw, elem->soil.depth, elem->soil.dmac,
        elem->soil.kmach, elem->soil.areafv, elem->soil.ksath) +
        EffKh(oppbank->ws.gw, oppbank->soil.depth, oppbank->soil.dmac,
        oppbank->soil.kmach, oppbank->soil.areafv, oppbank->soil.ksath));
    effk_nabr =
        EffKh(elem->ws.gw, elem->soil.depth, elem->soil.dmac, elem->soil.kmach,
        elem->soil.areafv, elem->soil.ksath);
#ifdef _ARITH_
    avg_ksat = 0.5 * (effk + effk_nabr);
#else
    avg_ksat = 2.0 / (1.0 / effk + 1.0 / effk_nabr);
#endif
    grad_y_sub = dif_y_sub / distance;
    *fluxsub = rivseg->shp.length * avg_ksat * grad_y_sub * avg_y_sub;

    /* Replace flux term */
    for (j = 0; j < NUM_EDGE; j++)
    {
        if (elem->nabr[j] == -ind)
        {
            elem->wf.ovlflow[j] = -(*fluxsurf);
            elem->wf.subsurf[j] = -(*fluxriv + *fluxsub);
            break;
        }
    }
}

double OvlFlowElemToRiver(double zmax, double surfh, double zbank,
    double stage, double zbed, double cwr, double length)
{
    double          flux;
    double          elem_h = zmax + surfh;
    double          rivseg_h = stage + zbed;

    /*
     * Panday and Hyakorn 2004 AWR Eqs. (23) and (24)
     */
    if (rivseg_h > elem_h)
    {
        if (elem_h > zbank)
        {
            /* Submerged weir */
            flux = cwr * 2.0 * sqrt(2.0 * GRAV) * length *
                sqrt(rivseg_h - elem_h) * (rivseg_h - zbank) / 3.0;
        }
        else
        {
            if (zbank < rivseg_h)
            {
                /* Free-flowing weir */
                flux = cwr * 2.0 * sqrt(2.0 * GRAV) * length *
                    sqrt(rivseg_h - zbank) * (rivseg_h - zbank) / 3.0;
            }
            else
            {
                flux = 0.0;
            }
        }
    }
    else if (surfh > DEPRSTG)
    {
        if (rivseg_h > zbank)
        {
            /* Submerged weir */
            flux = -cwr * 2.0 * sqrt(2.0 * GRAV) * length *
                sqrt(elem_h - rivseg_h) * (elem_h - zbank) / 3.0;
        }
        else
        {
            if (zbank < elem_h)
            {
                /* Free-flowing weir */
                flux = -cwr * 2.0 * sqrt(2.0 * GRAV) * length *
                    sqrt(elem_h - zbank) * (elem_h - zbank) / 3.0;
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

double _RiverWdthAreaPerim(int type, int riv_order, double riv_depth,
    double riv_coeff)
{
    double          eq_wid = 0.0;
    double          riv_area = 0.0;
    double          riv_perim = 0.0;
    double          ans;

    riv_depth = (riv_depth > 0.0) ? riv_depth : 0.0;

    switch (riv_order)
    {
        case RECTANGLE:
            eq_wid = riv_coeff;
            riv_area = riv_depth * riv_coeff;
            riv_perim = 2.0 * riv_depth + riv_coeff;
            break;
        case TRIANGLE:
            eq_wid = 2.0 *
                pow(riv_depth + RIVDPTHMIN, 1.0 / (riv_order - 1)) /
                pow(riv_coeff, 1.0 / (riv_order - 1));
            riv_area = riv_depth * riv_depth / riv_coeff;
            riv_perim =
                2.0 * riv_depth * sqrt(1.0 + riv_coeff * riv_coeff) / riv_coeff;
            break;
        case QUADRATIC:
            eq_wid = 2.0 *
                pow(riv_depth + RIVDPTHMIN, 1.0 / (riv_order - 1)) /
                pow(riv_coeff, 1.0 / (riv_order - 1));
            riv_area =
                4.0 * riv_depth * sqrt(riv_depth) / (3.0 * sqrt(riv_coeff));
            riv_perim =
                sqrt(riv_depth * (1.0 + 4.0 * riv_coeff * riv_depth) /
                riv_coeff) +
                log(2.0 * sqrt(riv_coeff * riv_depth) +
                sqrt(1.0 + 4.0 * riv_coeff * riv_depth)) / (2.0 * riv_coeff);
            break;
        case CUBIC:
            eq_wid = 2.0 *
                pow(riv_depth + RIVDPTHMIN, 1.0 / (riv_order - 1)) /
                pow(riv_coeff, 1.0 / (riv_order - 1));
            riv_area = 3.0 * pow(riv_depth, 4.0 / 3.0) /
                (2.0 * pow(riv_coeff, 1.0 / 3.0));
            riv_perim = 2.0 * (
                (pow(riv_depth * (1.0 + 9.0 * pow(riv_coeff, 2.0 / 3.0) *
                riv_depth), 0.5) / 3.0) +
                (log(3.0 * pow(riv_coeff, 1.0 / 3.0) *
                sqrt(riv_depth) + pow(1.0 + 9.0 *
                pow(riv_coeff, 2.0 / 3.0) * riv_depth, 0.5)) /
                (9.0 * pow(riv_coeff, 1.0 / 3.0))));
            break;
        default:
            PIHMprintf(VL_ERROR, "Error: River order %d is not defined.\n",
                riv_order);
            PIHMexit(EXIT_FAILURE);
    }

    switch (type)
    {
        case RIVER_WDTH:
            ans = eq_wid;
            break;
        case RIVER_AREA:
            ans = riv_area;
            break;
        case RIVER_PERIM:
            ans = riv_perim;
            break;
        default:
            ans = BADVAL;
            PIHMprintf(VL_ERROR,
                "Error: Return value type %d id not defined.\n", type);
            PIHMexit(EXIT_FAILURE);
    }

    return ans;
}

double ChanFlowRiverToRiver(const river_wstate_struct *ws,
    const river_topo_struct *topo, const shp_struct *shp,
    const matl_struct *matl, const river_wstate_struct *down_ws,
    const river_topo_struct *down_topo, const shp_struct *down_shp,
    const matl_struct *down_matl, int riv_mode)
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

    total_h = ws->stage + topo->zbed;
    perim = RiverPerim(shp->intrpl_ord, ws->stage, shp->coeff);

    total_h_down = down_ws->stage + down_topo->zbed;
    perim_down = RiverPerim(down_shp->intrpl_ord, down_ws->stage,
        down_shp->coeff);
    avg_perim = (perim + perim_down) / 2.0;
    avg_rough = (matl->rough + down_matl->rough) / 2.0;
    distance = 0.5 * (shp->length + down_shp->length);
    diff_h = (riv_mode == 1) ?
        (topo->zbed - down_topo->zbed) : (total_h - total_h_down);
    grad_h = diff_h / distance;
    avg_sf = (grad_h > 0.0) ? grad_h : RIVGRADMIN;
    crossa = RiverArea(shp->intrpl_ord, ws->stage, shp->coeff);
    crossa_down =
        RiverArea(down_shp->intrpl_ord, down_ws->stage, down_shp->coeff);
    avg_crossa = 0.5 * (crossa + crossa_down);
    avg_h = (avg_perim == 0.0) ? 0.0 : (avg_crossa / avg_perim);
    return OverLandFlow(avg_h, grad_h, avg_sf, crossa, avg_rough);
}

double SubFlowRiverToRiver(const river_wstate_struct *ws,
    const river_topo_struct *topo, const shp_struct *shp, double effk,
    const river_wstate_struct *down_ws, const river_topo_struct *down_topo,
    const shp_struct *down_shp, double effk_nabr)
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

    /* Lateral flux calculation between element beneath river (ebr)
     * and ebr */
    total_h = ws->gw + topo->zmin;
    total_h_down = down_ws->gw + down_topo->zmin;
    avg_wid = (shp->width + down_shp->width) / 2.0;
    diff_h = total_h - total_h_down;
    avg_h = AvgH(diff_h, ws->gw, down_ws->gw);
    distance = 0.5 * (shp->length + down_shp->length);
    grad_h = diff_h / distance;
    aquifer_depth = topo->zbed - topo->zmin;
#ifdef _ARITH_
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
            total_h = ws->gw + topo->zmin;
            total_h_down = bc->head + (topo->node_zmax - shp->depth);
            distance = 0.5 * shp->length;
            grad_h = (total_h - total_h_down) / distance;
            avg_h = AvgH(grad_h, ws->stage, bc->head);
            avg_perim = RiverPerim(shp->intrpl_ord, ws->stage, shp->coeff);
            crossa = RiverArea(shp->intrpl_ord, ws->stage, shp->coeff);
            avg_h = (avg_perim == 0.0) ? 0.0 : (crossa / avg_perim);
            discharge =
                OverLandFlow(avg_h, grad_h, grad_h, crossa, matl->rough);
            break;
        case NEUMANN:
            /* Neumann boundary condition */
            discharge = bc->flux;
            break;
        case ZERO_DPTH_GRAD:
            /* Zero-depth-gradient boundary conditions */
            distance = 0.5 * shp->length;
            grad_h = (topo->zbed - (topo->node_zmax - shp->depth)) / distance;
            avg_h = ws->stage;
            avg_perim = RiverPerim(shp->intrpl_ord, ws->stage, shp->coeff);
            crossa = RiverArea(shp->intrpl_ord, ws->stage, shp->coeff);
            discharge = sqrt(grad_h) * crossa * ((avg_perim > 0.0) ?
                pow(crossa / avg_perim, 2.0 / 3.0) : 0.0) / matl->rough;
            break;
        case CRIT_DPTH:
            /* Critical depth boundary conditions */
            crossa = RiverArea(shp->intrpl_ord, ws->stage, shp->coeff);
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
