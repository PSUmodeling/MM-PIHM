#include "pihm.h"

void RiverFlow (pihm_struct pihm)
{
    river_struct   *riv;
    river_struct   *down;
    elem_struct    *left;
    elem_struct    *right;

    int             i;
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
    double          wid;
    double          avg_ksat;
    double          effk_nabr;
    double          effk;
    double          aquifer_depth;
    double          grad_y_sub;
    double          avg_y_sub;
    double          dif_y_sub;
    double          wid_down;
    double          avg_wid;
    double          dt;

    dt = (double) pihm->ctrl.stepsize;

    /*
     * Lateral flux calculation between river-river and river-triangular
     * elements
     */
    for (i = 0; i < pihm->numriv; i++)
    {
        riv = &pihm->riv[i];

        total_y = riv->ws.stage + riv->topo.zbed;
        perim = RivPerim (riv->shp.intrpl_ord, riv->ws.stage, riv->shp.coeff);

        if (riv->down > 0)
        {
            down = &pihm->riv[riv->down - 1];

            /* Lateral flux calculation between river-river element */
            total_y_down = down->ws.stage + down->topo.zbed;
            perim_down =
                RivPerim (down->shp.intrpl_ord, down->ws.stage, down->shp.coeff);
            avg_perim = (perim + perim_down) / 2.0;
            avg_rough = (riv->matl.rough + down->matl.rough) / 2.0;
            distance = 0.5 * (riv->shp.length + down->shp.length);
            dif_y =
                (pihm->ctrl.riv_mode ==
                1) ? (riv->topo.zbed - down->topo.zbed) : (total_y -
                total_y_down);
            grad_y = dif_y / distance;
            avg_sf = (grad_y > 0.0) ? grad_y : RIVGRADMIN;
            crossa =
                RivArea (riv->shp.intrpl_ord, riv->ws.stage, riv->shp.coeff);
            crossa_down =
                RivArea (down->shp.intrpl_ord, down->ws.stage, down->shp.coeff);
            avg_crossa = 0.5 * (crossa + crossa_down);
            avg_y = (avg_perim == 0.0) ? 0.0 : (avg_crossa / avg_perim);
            riv->wf.fluxriv[1] =
                OverlandFlow (avg_y, grad_y, avg_sf, crossa, avg_rough);
            /* Accumulate to get in-flow for down segments:
             * [0] for inflow, [1] for outflow */
            down->wf.fluxriv[0] -= riv->wf.fluxriv[1];

            /* Lateral flux calculation between element beneath river (ebr)
             * and ebr */
            total_y = riv->ws.gw + riv->topo.zmin;
            total_y_down = down->ws.gw + down->topo.zmin;
            wid = EqWid (riv->shp.intrpl_ord, riv->shp.depth, riv->shp.coeff);
            wid_down =
                EqWid (down->shp.intrpl_ord, down->shp.depth,
                down->shp.coeff);
            avg_wid = (wid + wid_down) / 2.0;
            distance = 0.5 * (riv->shp.length + down->shp.length);
            dif_y_sub = total_y - total_y_down;
            avg_y_sub = AvgY (dif_y_sub, riv->ws.gw, down->ws.gw);
            grad_y_sub = dif_y_sub / distance;
            /* take care of macropore effect */
            aquifer_depth = riv->topo.zbed - riv->topo.zmin;
            left = &pihm->elem[riv->leftele - 1];
            right = &pihm->elem[riv->rightele - 1];
            effk =
                0.5 * (EffKH (left->soil.macropore, left->ws.gw,
                    left->soil.depth, left->soil.dmac, left->soil.kmach,
                    left->soil.areafv,
                    left->soil.ksath) + EffKH (right->soil.macropore,
                    right->ws.gw, right->soil.depth, right->soil.dmac,
                    right->soil.kmach, left->soil.areafv, right->soil.ksath));
            left = &pihm->elem[down->leftele - 1];
            right = &pihm->elem[down->rightele - 1];
            effk_nabr =
                0.5 * (EffKH (left->soil.macropore, left->ws.gw,
                    left->soil.depth, left->soil.dmac, left->soil.kmach,
                    left->soil.areafv,
                    left->soil.ksath) + EffKH (right->soil.macropore,
                    right->ws.gw, right->soil.depth, right->soil.dmac,
                    right->soil.kmach, left->soil.areafv, right->soil.ksath));
#ifdef _ARITH_
            avg_ksat = 0.5 * (effk + effk_nabr);
#else
            avg_ksat = 2.0 / (1.0 / effk + 1.0 / effk_nabr);
#endif
            /* Groundwater flow modeled by Darcy's law */
            riv->wf.fluxriv[9] = avg_ksat * grad_y_sub * avg_y_sub * avg_wid;
            /* Accumulate to get in-flow for down segments:
             * [10] for inflow, [9] for outflow */
            down->wf.fluxriv[10] -= riv->wf.fluxriv[9];
        }
        else
        {
            switch (riv->down)
            {
                case -1:
                    /* Dirichlet boundary condition */
                    total_y_down =
                        *riv->forc.riverbc + (riv->topo.node_zmax -
                        riv->shp.depth);
                    distance = 0.5 * riv->shp.length;
                    grad_y = (total_y - total_y_down) / distance;
                    avg_sf = grad_y;
                    avg_rough = riv->matl.rough;
                    avg_y = AvgY (grad_y, riv->ws.stage, *riv->forc.riverbc);
                    avg_perim = perim;
                    crossa =
                        RivArea (riv->shp.intrpl_ord, riv->ws.stage,
                        riv->shp.coeff);
                    avg_y = (avg_perim == 0.0) ? 0.0 : (crossa / avg_perim);
                    riv->wf.fluxriv[1] =
                        OverlandFlow (avg_y, grad_y, avg_sf, crossa,
                        avg_rough);
                    break;
                case -2:
                    /* Neumann boundary condition */
                    riv->wf.fluxriv[1] = *riv->forc.riverbc;
                    break;
                case -3:
                    /* Zero-depth-gradient boundary conditions */
                    distance = 0.5 * riv->shp.length;
                    grad_y =
                        (riv->topo.zbed - (riv->topo.node_zmax -
                            riv->shp.depth)) / distance;
                    avg_rough = riv->matl.rough;
                    avg_y = riv->ws.stage;
                    avg_perim = perim;
                    crossa =
                        RivArea (riv->shp.intrpl_ord, riv->ws.stage,
                        riv->shp.coeff);
                    riv->wf.fluxriv[1] =
                        sqrt (grad_y) * crossa * ((avg_perim >
                            0.0) ? pow (crossa / avg_perim,
                            2.0 / 3.0) : 0.0) / avg_rough;
                    break;
                case -4:
                    /* Critical depth boundary conditions */
                    crossa =
                        RivArea (riv->shp.intrpl_ord, riv->ws.stage,
                        riv->shp.coeff);
                    riv->wf.fluxriv[1] = crossa * sqrt (GRAV * riv->ws.stage);
                    break;
                default:
                    printf
                        ("Error: river routing boundary condition type is wrong!\n");
                    PihmExit (1);
            }
            /* Note: boundary condition for subsurface element can be changed.
             * Assumption: no flow condition */
            riv->wf.fluxriv[9] = 0.0;
        }

        left = &pihm->elem[riv->leftele - 1];
        right = &pihm->elem[riv->rightele - 1];

        if (riv->leftele > 0)
        {
            RiverToEle (riv, left, right, i + 1, &riv->wf.fluxriv[2],
                &riv->wf.fluxriv[4], &riv->wf.fluxriv[7], dt);
        }

        if (riv->rightele > 0)
        {
            RiverToEle (riv, right, left, i + 1, &riv->wf.fluxriv[3],
                &riv->wf.fluxriv[5], &riv->wf.fluxriv[8], dt);
        }

        avg_wid = EqWid (riv->shp.intrpl_ord, riv->ws.stage, riv->shp.coeff);
        if (riv->topo.zbed - (riv->ws.gw + riv->topo.zmin) > 0.0)
        {
            dif_y = riv->ws.stage;
        }
        else
        {
            dif_y = riv->ws.stage + riv->topo.zbed - (riv->ws.gw + riv->topo.zmin);
        }
        grad_y = dif_y / riv->matl.bedthick;
        riv->wf.fluxriv[6] =
            riv->matl.ksatv * avg_wid * riv->shp.length * grad_y;
    }
}

void RiverToEle (river_struct *riv, elem_struct *elem, elem_struct *oppbank,
    int ind, double *fluxsurf, double *fluxriv, double *fluxsub, double dt)
{
    double          total_y;
    double          dif_y_sub;
    double          avg_y_sub;
    double          effk;
    double          distance;
    double          grad_y_sub;
    double          effk_nabr;
    double          avg_ksat;
    double          aquifer_depth;
    int             j;

    total_y = riv->ws.stage + riv->topo.zbed;

    /* Lateral surface flux calculation between river-triangular element */
    *fluxsurf =
        OLFEleToRiv (elem->ps.surfavail + elem->topo.zmax, elem->topo.zmax,
        riv->matl.cwr, riv->topo.zmax, total_y, riv->shp.length);

    /* Lateral subsurface flux calculation between river-triangular element */
    dif_y_sub = (riv->ws.stage + riv->topo.zbed) - (elem->ws.gw + elem->topo.zmin);
    /* This is head at river edge representation */
    //avg_y_sub = ((md->riv[i].zmax-(md->ele[md->riv[i].leftele-1].zmax-md->ele[md->riv[i].leftele-1].zmin)+md->dummyy[md->riv[i].leftele-1 + 2*md->numele])>md->riv[i].zmin)?((md->riv[i].zmax-(md->ele[md->riv[i].leftele-1].zmax-md->ele[md->riv[i].leftele-1].zmin)+md->dummyy[md->riv[i].leftele-1 + 2*md->numele])-md->riv[i].zmin):0;
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
    avg_y_sub = AvgY (dif_y_sub, riv->ws.stage, avg_y_sub);
    effk = riv->matl.ksath;
    distance =
        sqrt (pow (riv->topo.x - elem->topo.x,
            2) + pow (riv->topo.y - elem->topo.y, 2));
    grad_y_sub = dif_y_sub / distance;
    /* Take into account macropore effect */
    effk_nabr =
        EffKH (elem->soil.macropore, elem->ws.gw, elem->soil.depth,
        elem->soil.dmac, elem->soil.kmach, elem->soil.areafv,
        elem->soil.ksath);
    avg_ksat = 0.5 * (effk + effk_nabr);
    *fluxriv = riv->shp.length * avg_ksat * grad_y_sub * avg_y_sub;

    /* Lateral flux between rectangular element (beneath river) and triangular
     * element */
    dif_y_sub = (riv->ws.gw + riv->topo.zmin) - (elem->ws.gw + elem->topo.zmin);
    /* This is head at river edge representation */
    //avg_y_sub = ((md->riv[i].zmax-(md->ele[md->riv[i].leftele-1].zmax-md->ele[md->riv[i].leftele-1].zmin)+md->dummyy[md->riv[i].leftele-1 + 2*md->numele])>md->riv[i].zmin)?md->riv[i].zmin-(md->riv[i].zmax-(md->ele[md->riv[i].leftele-1].zmax-md->ele[md->riv[i].leftele-1].zmin)):md->dummyy[md->riv[i].leftele-1 + 2*md->numele];
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
    avg_y_sub = AvgY (dif_y_sub, riv->ws.gw, avg_y_sub);
    aquifer_depth = riv->topo.zbed - riv->topo.zmin;
    effk =
        0.5 * (EffKH (elem->soil.macropore, elem->ws.gw, elem->soil.depth,
            elem->soil.dmac, elem->soil.kmach, elem->soil.areafv,
            elem->soil.ksath) + EffKH (oppbank->soil.macropore, oppbank->ws.gw,
            oppbank->soil.depth, oppbank->soil.dmac, oppbank->soil.kmach,
            oppbank->soil.areafv, oppbank->soil.ksath));
    effk_nabr =
        EffKH (elem->soil.macropore, elem->ws.gw, elem->soil.depth,
        elem->soil.dmac, elem->soil.kmach, elem->soil.areafv,
        elem->soil.ksath);
#ifdef _ARITH_
    avg_ksat = 0.5 * (effk + effk_nabr);
#else
    avg_ksat = 2.0 / (1.0 / effk + 1.0 / effk_nabr);
#endif
    grad_y_sub = dif_y_sub / distance;
    /* Take into account macropore effect */
    *fluxsub = riv->shp.length * avg_ksat * grad_y_sub * avg_y_sub;

    /* Replace flux term */
    for (j = 0; j < 3; j++)
    {
        if (elem->nabr[j] == - ind)
        {
            elem->wf.fluxsurf[j] = -(*fluxsurf);
            elem->wf.fluxsub[j] = -(*fluxriv + *fluxsub);
            break;
        }
    }
}

double EqWid (int riv_order, double riv_depth, double riv_coeff)
{
    double          eq_wid;

    riv_depth = (riv_depth > 0.0) ? riv_depth : 0.0;

    switch (riv_order)
    {
        case RECTANGLE:
            eq_wid = riv_coeff;
            break;
        case TRIANGLE:
            eq_wid = 2.0 * pow (riv_depth + RIVDPTHMIN, 1.0 / (riv_order - 1)) /
                pow (riv_coeff, 1.0 / (riv_order - 1));
            break;
        case QUADRATIC:
            eq_wid = 2.0 * pow (riv_depth + RIVDPTHMIN, 1.0 / (riv_order - 1)) /
                pow (riv_coeff, 1.0 / (riv_order - 1));
            break;
        case CUBIC:
            eq_wid = 2.0 * pow (riv_depth + RIVDPTHMIN, 1.0 / (riv_order - 1)) /
                pow (riv_coeff, 1.0 / (riv_order - 1));
            break;
        default:
            printf ("Error: River order %d is not defined!\n", riv_order);
            PihmExit (1);
    }
    return (eq_wid);
}

double OLFEleToRiv (double eleytot, double elez, double cwr, double rivzmax,
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
            flux =
                cwr * 2.0 * sqrt (2.0 * GRAV) * length * sqrt (rivytot -
                eleytot) * (rivytot - zbank) / 3.0;
        }
        else
        {
            if (zbank < rivytot)
            {
                /* Free-flowing weir */
                flux =
                    cwr * 2.0 * sqrt (2.0 * GRAV) * length * sqrt (rivytot -
                    zbank) * (rivytot - zbank) / 3.0;
            }
            else
            {
                flux = 0.0;
            }
        }
    }
    else
    {
        if (rivytot > zbank)
        {
            /* Submerged weir */
            flux =
                -cwr * 2.0 * sqrt (2.0 * GRAV) * length * sqrt (eleytot -
                rivytot) * (eleytot - zbank) / 3.0;
        }
        else
        {
            if (zbank < eleytot)
            {
                /* Free-flowing weir */
                flux =
                    -cwr * 2.0 * sqrt (2.0 * GRAV) * length * sqrt (eleytot -
                    zbank) * (eleytot - zbank) / 3.0;
            }
            else
            {
                flux = 0.0;
            }
        }
    }

    return (flux);
}

double RivArea (int riv_order, double riv_depth, double riv_coeff)
{
    double          riv_area;

    riv_depth = (riv_depth > 0.0) ? riv_depth : 0.0;

    switch (riv_order)
    {
        case RECTANGLE:
            riv_area = riv_depth * riv_coeff;
            break;
        case TRIANGLE:
            riv_area = pow (riv_depth, 2) / riv_coeff;
            break;
        case QUADRATIC:
            riv_area =
                4.0 * pow (riv_depth, 1.5) / (3.0 * pow (riv_coeff, 0.5));
            break;
        case CUBIC:
            riv_area =
                3.0 * pow (riv_depth, 4.0 / 3.0) / (2.0 * pow (riv_coeff,
                    1.0 / 3.0));
            break;
        default:
            printf ("Error: River order %d is not defined!\n", riv_order);
            PihmExit (1);
    }

    return (riv_area);
}

double RivPerim (int riv_order, double riv_depth, double riv_coeff)
{
    double          riv_perim;

    riv_depth = (riv_depth > 0.0) ? riv_depth : 0.0;

    switch (riv_order)
    {
        case RECTANGLE:
            riv_perim = 2.0 * riv_depth + riv_coeff;
            break;
        case TRIANGLE:
            riv_perim =
                2.0 * riv_depth * pow (1.0 + pow (riv_coeff, 2),
                0.5) / riv_coeff;
            break;
        case QUADRATIC:
            riv_perim =
                (pow (riv_depth * (1.0 +
                        4.0 * riv_coeff * riv_depth) / riv_coeff,
                    0.5)) + (log (2.0 * pow (riv_coeff * riv_depth,
                        0.5) + pow (1.0 + 4.0 * riv_coeff * riv_depth,
                        0.5)) / (2.0 * riv_coeff));
            break;
        case CUBIC:
            riv_perim =
                2.0 * ((pow (riv_depth * (1.0 + 9.0 * pow (riv_coeff,
                                2.0 / 3.0) * riv_depth),
                        0.5) / 3.0) + (log (3.0 * pow (riv_coeff,
                            1.0 / 3.0) * pow (riv_depth,
                            0.5) + pow (1.0 + 9.0 * pow (riv_coeff,
                                2.0 / 3.0) * riv_depth,
                            0.5)) / (9.0 * pow (riv_coeff, 1.0 / 3.0))));
            break;
        default:
            printf ("Error: River order %d is not defined!\n", riv_order);
            PihmExit (1);
    }
    return (riv_perim);
}
