#include "pihm.h"

void RiverFlow(elem_struct *elem, river_struct *riv, int riv_mode)
{
    int             i;

    /*
     * Lateral flux calculation between river-river and river-triangular
     * elements
     */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river_struct   *down;
        elem_struct    *left;
        elem_struct    *right;
        double          total_y;
        double          perim;
        double          total_y_down;
        double          perim_down;
        double          avg_perim;
        double          avg_rough;
        double          distance;
        double          dif_y;
        double          grad_y;
        double          avg_sf;
        double          crossa;
        double          crossa_down;
        double          avg_crossa;
        double          avg_y;
        double          avg_ksat;
        double          effk_nabr;
        double          effk;
        double          aquifer_depth;
        double          grad_y_sub;
        double          avg_y_sub;
        double          dif_y_sub;
        double          avg_wid;

        total_y = riv[i].ws.stage + riv[i].topo.zbed;
        perim = RivPerim(riv[i].shp.intrpl_ord, riv[i].ws.stage,
            riv[i].shp.coeff);

        if (riv[i].down > 0)
        {
            down = &riv[riv[i].down - 1];

            /* Lateral flux calculation between river-river element */
            total_y_down = down->ws.stage + down->topo.zbed;
            perim_down =
                RivPerim(down->shp.intrpl_ord, down->ws.stage, down->shp.coeff);
            avg_perim = (perim + perim_down) / 2.0;
            avg_rough = (riv[i].matl.rough + down->matl.rough) / 2.0;
            distance = 0.5 * (riv[i].shp.length + down->shp.length);
            dif_y = (riv_mode == 1) ?
                (riv[i].topo.zbed - down->topo.zbed) : (total_y - total_y_down);
            grad_y = dif_y / distance;
            avg_sf = (grad_y > 0.0) ? grad_y : RIVGRADMIN;
            crossa =
                RivArea(riv[i].shp.intrpl_ord, riv[i].ws.stage,
                riv[i].shp.coeff);
            crossa_down =
                RivArea(down->shp.intrpl_ord, down->ws.stage, down->shp.coeff);
            avg_crossa = 0.5 * (crossa + crossa_down);
            avg_y = (avg_perim == 0.0) ? 0.0 : (avg_crossa / avg_perim);
            riv[i].wf.rivflow[DOWN_CHANL2CHANL] =
                OverlandFlow(avg_y, grad_y, avg_sf, crossa, avg_rough);
            /* Accumulate to get in-flow for down segments */
            down->wf.rivflow[UP_CHANL2CHANL] -=
                riv[i].wf.rivflow[DOWN_CHANL2CHANL];

            /* Lateral flux calculation between element beneath river (ebr)
             * and ebr */
            total_y = riv[i].ws.gw + riv[i].topo.zmin;
            total_y_down = down->ws.gw + down->topo.zmin;
            avg_wid = (riv[i].shp.width + down->shp.width) / 2.0;
            dif_y_sub = total_y - total_y_down;
            avg_y_sub = AvgY(dif_y_sub, riv[i].ws.gw, down->ws.gw);
            grad_y_sub = dif_y_sub / distance;
            aquifer_depth = riv[i].topo.zbed - riv[i].topo.zmin;
            left = &elem[riv[i].leftele - 1];
            right = &elem[riv[i].rightele - 1];
            effk = 0.5 *
                (EffKH(left->ws.gw, left->soil.depth, left->soil.dmac,
                left->soil.kmach, left->soil.areafv, left->soil.ksath) +
                EffKH(right->ws.gw, right->soil.depth, right->soil.dmac,
                right->soil.kmach, left->soil.areafv, right->soil.ksath));
            left = &elem[down->leftele - 1];
            right = &elem[down->rightele - 1];
            effk_nabr = 0.5 *
                (EffKH(left->ws.gw, left->soil.depth, left->soil.dmac,
                left->soil.kmach, left->soil.areafv, left->soil.ksath) +
                EffKH(right->ws.gw, right->soil.depth, right->soil.dmac,
                right->soil.kmach, left->soil.areafv, right->soil.ksath));
#ifdef _ARITH_
            avg_ksat = 0.5 * (effk + effk_nabr);
#else
            avg_ksat = 2.0 / (1.0 / effk + 1.0 / effk_nabr);
#endif
            /* Groundwater flow modeled by Darcy's law */
            riv[i].wf.rivflow[DOWN_AQUIF2AQUIF] =
                avg_ksat * grad_y_sub * avg_y_sub * avg_wid;
            /* Accumulate to get in-flow for down segments */
            down->wf.rivflow[UP_AQUIF2AQUIF] -=
                riv[i].wf.rivflow[DOWN_AQUIF2AQUIF];
        }
        else
        {
            switch (riv[i].down)
            {
                case DIRICHLET:
                    /* Dirichlet boundary condition */
                    total_y_down = riv[i].bc.head +
                        (riv[i].topo.node_zmax - riv[i].shp.depth);
                    distance = 0.5 * riv[i].shp.length;
                    grad_y = (total_y - total_y_down) / distance;
                    avg_sf = grad_y;
                    avg_rough = riv[i].matl.rough;
                    avg_y = AvgY(grad_y, riv[i].ws.stage, riv[i].bc.head);
                    avg_perim = perim;
                    crossa = RivArea(riv[i].shp.intrpl_ord, riv[i].ws.stage,
                        riv[i].shp.coeff);
                    avg_y = (avg_perim == 0.0) ? 0.0 : (crossa / avg_perim);
                    riv[i].wf.rivflow[DOWN_CHANL2CHANL] =
                        OverlandFlow(avg_y, grad_y, avg_sf, crossa, avg_rough);
                    break;
                case NEUMANN:
                    /* Neumann boundary condition */
                    riv[i].wf.rivflow[DOWN_CHANL2CHANL] = riv[i].bc.flux;
                    break;
                case ZERO_DPTH_GRAD:
                    /* Zero-depth-gradient boundary conditions */
                    distance = 0.5 * riv[i].shp.length;
                    grad_y = (riv[i].topo.zbed -
                        (riv[i].topo.node_zmax - riv[i].shp.depth)) / distance;
                    avg_rough = riv[i].matl.rough;
                    avg_y = riv[i].ws.stage;
                    avg_perim = perim;
                    crossa = RivArea(riv[i].shp.intrpl_ord, riv[i].ws.stage,
                        riv[i].shp.coeff);
                    riv[i].wf.rivflow[DOWN_CHANL2CHANL] =
                        sqrt(grad_y) * crossa *
                        ((avg_perim > 0.0) ?
                        pow(crossa / avg_perim, 2.0 / 3.0) : 0.0) / avg_rough;
                    break;
                case CRIT_DPTH:
                    /* Critical depth boundary conditions */
                    crossa = RivArea(riv[i].shp.intrpl_ord, riv[i].ws.stage,
                        riv[i].shp.coeff);
                    riv[i].wf.rivflow[DOWN_CHANL2CHANL] =
                        crossa * sqrt(GRAV * riv[i].ws.stage);
                    break;
                default:
                    PIHMprintf(VL_ERROR,
                        "Error: River routing boundary condition type (%d) "
                        "is not recognized.\n", riv[i].down);
                    PIHMexit(EXIT_FAILURE);
            }
            /* Note: boundary condition for subsurface element can be changed.
             * Assumption: no flow condition */
            riv[i].wf.rivflow[DOWN_AQUIF2AQUIF] = 0.0;
        }

        left = &elem[riv[i].leftele - 1];
        right = &elem[riv[i].rightele - 1];

        if (riv[i].leftele > 0)
        {
            RiverToEle(&riv[i], left, right, i + 1, riv[i].topo.dist_left,
                &riv[i].wf.rivflow[LEFT_SURF2CHANL],
                &riv[i].wf.rivflow[LEFT_AQUIF2CHANL],
                &riv[i].wf.rivflow[LEFT_AQUIF2AQUIF]);
        }

        if (riv[i].rightele > 0)
        {
            RiverToEle(&riv[i], right, left, i + 1, riv[i].topo.dist_right,
                &riv[i].wf.rivflow[RIGHT_SURF2CHANL],
                &riv[i].wf.rivflow[RIGHT_AQUIF2CHANL],
                &riv[i].wf.rivflow[RIGHT_AQUIF2AQUIF]);
        }

        if (riv[i].topo.zbed - (riv[i].ws.gw + riv[i].topo.zmin) > 0.0)
        {
            dif_y = riv[i].ws.stage;
        }
        else
        {
            dif_y = riv[i].ws.stage + riv[i].topo.zbed -
                (riv[i].ws.gw + riv[i].topo.zmin);
        }
        grad_y = dif_y / riv[i].matl.bedthick;
        riv[i].wf.rivflow[CHANL_LKG] =
            riv[i].matl.ksatv * riv[i].shp.width * riv[i].shp.length * grad_y;
    }
}

void RiverToEle(river_struct *riv, elem_struct *elem, elem_struct *oppbank,
    int ind, double distance, double *fluxsurf, double *fluxriv,
    double *fluxsub)
{
    double          total_y;
    double          dif_y_sub;
    double          avg_y_sub;
    double          effk;
    double          grad_y_sub;
    double          effk_nabr;
    double          avg_ksat;
    double          aquifer_depth;
    int             j;

    total_y = riv->ws.stage + riv->topo.zbed;

    /* Lateral surface flux calculation between river-triangular element */
    *fluxsurf = OLFEleToRiv(elem->ws.surfh + elem->topo.zmax, elem->topo.zmax,
        riv->matl.cwr, riv->topo.zmax, total_y, riv->shp.length);

    /* Lateral subsurface flux calculation between river-triangular element */
    dif_y_sub =
        (riv->ws.stage + riv->topo.zbed) - (elem->ws.gw + elem->topo.zmin);
    /* This is head in neighboring cell represention */
    if (elem->topo.zmin > riv->topo.zbed)
    {
        avg_y_sub = elem->ws.gw;
    }
    else if (elem->topo.zmin + elem->ws.gw > riv->topo.zbed)
    {
        avg_y_sub = elem->topo.zmin + elem->ws.gw - riv->topo.zbed;
    }
    else
    {
        avg_y_sub = 0.0;
    }
    avg_y_sub = AvgY(dif_y_sub, riv->ws.stage, avg_y_sub);
    effk = riv->matl.ksath;
    grad_y_sub = dif_y_sub / distance;
    /* Take into account macropore effect */
    effk_nabr = EffKH(elem->ws.gw, elem->soil.depth, elem->soil.dmac,
        elem->soil.kmach, elem->soil.areafv, elem->soil.ksath);
    avg_ksat = 0.5 * (effk + effk_nabr);
    *fluxriv = riv->shp.length * avg_ksat * grad_y_sub * avg_y_sub;

    /* Lateral flux between rectangular element (beneath river) and triangular
     * element */
    dif_y_sub = (riv->ws.gw + riv->topo.zmin) - (elem->ws.gw + elem->topo.zmin);
    /* this is head in neighboring cell represention */
    if (elem->topo.zmin > riv->topo.zbed)
    {
        avg_y_sub = 0.0;
    }
    else if (elem->topo.zmin + elem->ws.gw > riv->topo.zbed)
    {
        avg_y_sub = riv->topo.zbed - elem->topo.zmin;
    }
    else
    {
        avg_y_sub = elem->ws.gw;
    }
    avg_y_sub = AvgY(dif_y_sub, riv->ws.gw, avg_y_sub);
    aquifer_depth = riv->topo.zbed - riv->topo.zmin;
    effk = 0.5 *
        (EffKH(elem->ws.gw, elem->soil.depth, elem->soil.dmac,
        elem->soil.kmach, elem->soil.areafv, elem->soil.ksath) +
        EffKH(oppbank->ws.gw, oppbank->soil.depth, oppbank->soil.dmac,
        oppbank->soil.kmach, oppbank->soil.areafv, oppbank->soil.ksath));
    effk_nabr =
        EffKH(elem->ws.gw, elem->soil.depth, elem->soil.dmac, elem->soil.kmach,
        elem->soil.areafv, elem->soil.ksath);
#ifdef _ARITH_
    avg_ksat = 0.5 * (effk + effk_nabr);
#else
    avg_ksat = 2.0 / (1.0 / effk + 1.0 / effk_nabr);
#endif
    grad_y_sub = dif_y_sub / distance;
    *fluxsub = riv->shp.length * avg_ksat * grad_y_sub * avg_y_sub;

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

double OLFEleToRiv(double eleytot, double elez, double cwr, double rivzmax,
    double rivytot, double length)
{
    double          flux;
    double          zbank;

    /*
     * Panday and Hyakorn 2004 AWR Eqs. (23) and (24)
     */
    zbank = (rivzmax < elez) ? elez : rivzmax;

    if (rivytot > eleytot)
    {
        if (eleytot > zbank)
        {
            /* Submerged weir */
            flux = cwr * 2.0 * sqrt(2.0 * GRAV) * length *
                sqrt(rivytot - eleytot) * (rivytot - zbank) / 3.0;
        }
        else
        {
            if (zbank < rivytot)
            {
                /* Free-flowing weir */
                flux = cwr * 2.0 * sqrt(2.0 * GRAV) * length *
                    sqrt(rivytot - zbank) * (rivytot - zbank) / 3.0;
            }
            else
            {
                flux = 0.0;
            }
        }
    }
    else if (eleytot - elez > DEPRSTG)
    {
        if (rivytot > zbank)
        {
            /* Submerged weir */
            flux = -cwr * 2.0 * sqrt(2.0 * GRAV) * length *
                sqrt(eleytot - rivytot) * (eleytot - zbank) / 3.0;
        }
        else
        {
            if (zbank < eleytot)
            {
                /* Free-flowing weir */
                flux = -cwr * 2.0 * sqrt(2.0 * GRAV) * length *
                    sqrt(eleytot - zbank) * (eleytot - zbank) / 3.0;
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

double _RivWdthAreaPerim(int type, int riv_order, double riv_depth,
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
