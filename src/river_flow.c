#include "pihm.h"

void RiverFlow(elem_struct elem[], river_struct river[])
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

        InitRiverWFlux(&river[i].wf);

        if (river[i].down > 0)
        {
            // Boundary conditions
            // When a downstream segment is present, boundary conditions are always applied to the upstream node
            if (river[i].attrib.riverbc_type != 0)
            {
                river[i].wf.rivflow[UPSTREAM] += BoundFluxRiver(river[i].attrib.riverbc_type, &river[i].topo, &river[i].shp, &river[i].matl, &river[i].bc, &river[i].ws);
            }

            // Channel flow between river-river segments
            down = &river[river[i].down - 1];

            river[i].wf.rivflow[DOWNSTREAM] = ChannelFlowRiverToRiver(&river[i], down);
        }
        else
        {
            // Outlet flux
            river[i].wf.rivflow[DOWNSTREAM] = OutletFlux(river[i].down, &river[i].topo, &river[i].shp, &river[i].matl, &river[i].bc, &river[i].ws);
        }

        // Flux between river segments and triangular elements
        left = &elem[river[i].left - 1];
        right = &elem[river[i].right - 1];

        RiverToElem(&river[i], left, right);
    }

    // Accumulate to get in-flow for down segments
    // NOTE: Upstream flux summation must be calculated outside OMP to avoid different threads accessing the same
    // variable at the same time
    for (i = 0; i < nriver; i++)
    {
        river_struct   *down;

        if (river[i].down > 0)
        {
            down = &river[river[i].down - 1];

            down->wf.rivflow[UPSTREAM] -= river[i].wf.rivflow[DOWNSTREAM];
        }
    }
}

void RiverToElem(river_struct *river_ptr, elem_struct *left, elem_struct *right)
{
    if (river_ptr->left > 0)
    {
        double          effk_left;

        river_ptr->wf.rivflow[SURF_LEFT] = OvlFlowElemToRiver(river_ptr, left);

        effk_left = EffKh(left->ws.gw, &left->soil);
        river_ptr->wf.rivflow[AQUIFER_LEFT] = ChannelFlowElemToRiver(effk_left, river_ptr->topo.dist_left, river_ptr, left);
    }

    if (river_ptr->right > 0)
    {
        double          effk_right;

        river_ptr->wf.rivflow[SURF_RIGHT] = OvlFlowElemToRiver(river_ptr, right);

        effk_right = EffKh(right->ws.gw, &right->soil);
        river_ptr->wf.rivflow[AQUIFER_RIGHT] = ChannelFlowElemToRiver(effk_right, river_ptr->topo.dist_right, river_ptr, right);
    }
}

double OvlFlowElemToRiver(const river_struct *river_ptr, elem_struct *bank)
{
    int             j;
    double          z_bank;
    double          flux;
    double          bank_h;
    double          river_h;

    z_bank = (river_ptr->topo.zmax > bank->topo.zmax) ? river_ptr->topo.zmax : bank->topo.zmax;

    bank_h = bank->topo.zmax + bank->ws.surfh;
    river_h = river_ptr->topo.zbed + river_ptr->ws.stage;

    // Panday and Hyakorn 2004 AWR Eqs. (23) and (24)
    if (river_h > bank_h)
    {
        if (bank_h > z_bank)
        {
            // Submerged weir
            flux = river_ptr->matl.cwr * 2.0 * sqrt(2.0 * GRAV) * river_ptr->shp.length * sqrt(river_h - bank_h) * (river_h - z_bank) / 3.0;
        }
        else
        {
            // Free-flowing weir
            flux = (z_bank < river_h) ? river_ptr->matl.cwr * 2.0 * sqrt(2.0 * GRAV) * river_ptr->shp.length * sqrt(river_h - z_bank) * (river_h - z_bank) / 3.0 : 0.0;
        }
    }
    else if (bank->ws.surfh > DEPRSTG)
    {
        if (river_h > z_bank)
        {
            // Submerged weir
            flux = -river_ptr->matl.cwr * 2.0 * sqrt(2.0 * GRAV) * river_ptr->shp.length * sqrt(bank_h - river_h) * (bank_h - z_bank) / 3.0;
        }
        else
        {
            // Free-flowing weir
            flux = (z_bank < bank_h) ? -river_ptr->matl.cwr * 2.0 * sqrt(2.0 * GRAV) * river_ptr->shp.length * sqrt(bank_h - z_bank) * (bank_h - z_bank) / 3.0 : 0.0;
        }
    }
    else
    {
        flux = 0.0;
    }

    for (j = 0; j < NUM_EDGE; j++)
    {
        if (bank->nabr_river[j] == river_ptr->ind)
        {
            bank->wf.overland[j] = -flux;
            break;
        }
    }

    return flux;
}

double ChannelFlowRiverToRiver(const river_struct *river_ptr, const river_struct *down)
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

    total_h = river_ptr->ws.stage + river_ptr->topo.zbed;
    perim = RiverPerim(river_ptr->shp.intrpl_ord, river_ptr->ws.stage, river_ptr->shp.coeff);

    total_h_down = down->ws.stage + down->topo.zbed;
    perim_down = RiverPerim(down->shp.intrpl_ord, down->ws.stage, down->shp.coeff);

    avg_perim = 0.5 * (perim + perim_down);
    avg_rough = 0.5 * (river_ptr->matl.rough + down->matl.rough);
    distance = 0.5 * (river_ptr->shp.length + down->shp.length);

    diff_h = total_h - total_h_down;
    grad_h = diff_h / distance;
    avg_sf = (grad_h > 0.0) ? grad_h : RIVGRADMIN;
    crossa = RiverCrossSectArea(river_ptr->shp.intrpl_ord, river_ptr->ws.stage, river_ptr->shp.coeff);
    crossa_down = RiverCrossSectArea(down->shp.intrpl_ord, down->ws.stage, down->shp.coeff);
    avg_crossa = 0.5 * (crossa + crossa_down);
    avg_h = (avg_perim == 0.0) ? 0.0 : avg_crossa / avg_perim;

    return OverLandFlow(avg_h, grad_h, avg_sf, crossa, avg_rough);
}

double OutletFlux(int down, const river_topo_struct *topo, const shp_struct *shp, const matl_struct *matl, const river_bc_struct *bc, const river_wstate_struct *ws)
{
    double          total_h;
    double          total_h_down;
    double          distance;
    double          grad_h;
    double          avg_h;
    double          avg_perim;
    double          cross_area;
    double          discharge = 0.0;

    switch (down)
    {
        case OUTLET_DIRICHLET:
            // Dirichlet boundary condition
            total_h = ws->stage + topo->zbed;
            total_h_down = bc->head;
            distance = 0.5 * shp->length;
            grad_h = (total_h - total_h_down) / distance;
            avg_perim = RiverPerim(shp->intrpl_ord, ws->stage, shp->coeff);
            cross_area = RiverCrossSectArea(shp->intrpl_ord, ws->stage, shp->coeff);
            avg_h = (avg_perim == 0.0) ? 0.0 : cross_area / avg_perim;
            discharge = OverLandFlow(avg_h, grad_h, grad_h, cross_area, matl->rough);
            break;
        case OUTLET_NEUMANN:
            // Neumann boundary condition
            discharge = -bc->flux;
            break;
        case ZERO_DPTH_GRAD:
            // Zero-depth-gradient boundary conditions
            distance = 0.5 * shp->length;
            grad_h = (topo->zbed - (topo->node_zmax - shp->depth)) / distance;
            avg_h = ws->stage;
            avg_perim = RiverPerim(shp->intrpl_ord, ws->stage, shp->coeff);
            cross_area = RiverCrossSectArea(shp->intrpl_ord, ws->stage, shp->coeff);
            discharge = sqrt(grad_h) * cross_area * ((avg_perim > 0.0) ? pow(cross_area / avg_perim, 2.0 / 3.0) : 0.0) / matl->rough;
            break;
        case CRIT_DPTH:
            // Critical depth boundary conditions
            cross_area = RiverCrossSectArea(shp->intrpl_ord, ws->stage, shp->coeff);
            discharge = cross_area * sqrt(GRAV * ws->stage);
            break;
        default:
            pihm_printf(VL_ERROR, "Error: River routing boundary condition type (%d) is not recognized.\n", down);
            pihm_exit(EXIT_FAILURE);
    }

    return discharge;
}

double BoundFluxRiver(int riverbc_type, const river_topo_struct *topo, const shp_struct *shp, const matl_struct *matl, const river_bc_struct *bc, const river_wstate_struct *ws)
{
    double          total_h;
    double          total_h_down;
    double          distance;
    double          grad_h;
    double          avg_h;
    double          avg_perim;
    double          cross_area;
    double          flux = 0.0;

    if (riverbc_type > 0)
    {
        // Dirichlet boundary condition
        total_h = ws->stage + topo->zbed;
        total_h_down = bc->head;
        distance = 0.5 * shp->length;
        grad_h = (total_h - total_h_down) / distance;
        avg_h = AvgH(grad_h, ws->stage, (bc->head - topo->zbed > 0.0) ? (bc->head - topo->zbed) : 0.0);
        avg_perim = RiverPerim(shp->intrpl_ord, ws->stage, shp->coeff);
        cross_area = RiverCrossSectArea(shp->intrpl_ord, ws->stage, shp->coeff);
        avg_h = (avg_perim == 0.0) ? 0.0 : cross_area / avg_perim;
        flux = OverLandFlow(avg_h, grad_h, grad_h, cross_area, matl->rough);
    }
    else if (riverbc_type < 0)
    {
        // Neumann boundary condition
        flux = -bc->flux;
    }

    return flux;
}

double ChannelFlowElemToRiver(double effk, double distance, const river_struct *river_ptr, elem_struct *bank)
{
    int             j;
    double          diff_h;
    double          avg_h;
    double          grad_h;
    double          avg_ksat;
    double          flux;

    diff_h = (river_ptr->ws.stage + river_ptr->topo.zbed) - (bank->ws.gw + bank->topo.zmin);

    // This is head in neighboring cell representation
    if (bank->topo.zmin > river_ptr->topo.zbed)
    {
        avg_h = bank->ws.gw;
    }
    else if (bank->topo.zmin + bank->ws.gw > river_ptr->topo.zbed)
    {
        avg_h = bank->topo.zmin + bank->ws.gw - river_ptr->topo.zbed;
    }
    else
    {
        avg_h = 0.0;
    }
    avg_h = AvgH(diff_h, river_ptr->ws.stage, avg_h);

    grad_h = diff_h / distance;

    avg_ksat = 0.5 * (effk + river_ptr->matl.ksath);

    flux = river_ptr->shp.length * avg_ksat * grad_h * avg_h;

    for (j = 0; j < NUM_EDGE; j++)
    {
        if (bank->nabr_river[j] == river_ptr->ind)
        {
            bank->wf.subsurf[j] -= flux;
            break;
        }
    }

    return flux;
}

double RiverCrossSectArea(int order, double depth, double coeff)
{
    double          cs_area = 0.0;

    depth = MAX(depth, 0.0);

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
            cs_area = 3.0 * pow(depth, 4.0 / 3.0) / (2.0 * pow(coeff, 1.0 / 3.0));
            break;
        default:
            pihm_printf(VL_ERROR, "Error: River order %d is not defined.\n", order);
            pihm_exit(EXIT_FAILURE);
    }

    return cs_area;
}

double RiverPerim(int order, double depth, double coeff)
{
    double          perim = 0.0;

    depth = MAX(depth, 0.0);

    switch (order)
    {
        case RECTANGLE:
            perim = 2.0 * depth + coeff;
            break;
        case TRIANGLE:
            perim = 2.0 * depth * sqrt(1.0 + coeff * coeff) / coeff;
            break;
        case QUADRATIC:
            perim = sqrt(depth * (1.0 + 4.0 * coeff * depth) / coeff) + log(2.0 * sqrt(coeff * depth) +
                sqrt(1.0 + 4.0 * coeff * depth)) / (2.0 * coeff);
            break;
        case CUBIC:
            perim = 2.0 * (pow(depth * (1.0 + 9.0 * pow(coeff, 2.0 / 3.0) * depth), 0.5) / 3.0 + log(3.0 * pow(coeff, 1.0 / 3.0) * sqrt(depth) +
                pow(1.0 + 9.0 * pow(coeff, 2.0 / 3.0) * depth, 0.5)) / (9.0 * pow(coeff, 1.0 / 3.0)));
            break;
        default:
            pihm_printf(VL_ERROR, "Error: River order %d is not defined.\n", order);
            pihm_exit(EXIT_FAILURE);
    }

    return perim;
}
