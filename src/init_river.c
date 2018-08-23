#include "pihm.h"

void InitRiver(river_struct *river, elem_struct *elem,
    const rivtbl_struct *rivtbl, const shptbl_struct *shptbl,
    const matltbl_struct *matltbl, const meshtbl_struct *meshtbl,
    const calib_struct *cal)
{
    int             i;

    for (i = 0; i < nriver; i++)
    {
        int             j;

        river[i].ind = i + 1;
        river[i].leftele = rivtbl->leftele[i];
        river[i].rightele = rivtbl->rightele[i];
        river[i].fromnode = rivtbl->fromnode[i];
        river[i].tonode = rivtbl->tonode[i];
        river[i].down = rivtbl->down[i];

        for (j = 0; j < NUM_EDGE; j++)
        {
            /* Note: use element nabr < 0 for river identification */
            if (elem[river[i].leftele - 1].nabr[j] == river[i].rightele)
            {
                elem[river[i].leftele - 1].nabr[j] = -(i + 1);
            }
            if (elem[river[i].rightele - 1].nabr[j] == river[i].leftele)
            {
                elem[river[i].rightele - 1].nabr[j] = -(i + 1);
            }
        }

        river[i].topo.x = 0.5 *
            (meshtbl->x[river[i].fromnode - 1] +
            meshtbl->x[river[i].tonode - 1]);
        river[i].topo.y = 0.5 *
            (meshtbl->y[river[i].fromnode - 1] +
            meshtbl->y[river[i].tonode - 1]);
        river[i].topo.zmax = 0.5 *
            (meshtbl->zmax[river[i].fromnode - 1] +
            meshtbl->zmax[river[i].tonode - 1]);
        river[i].topo.zmin = river[i].topo.zmax -
            (0.5 * (elem[river[i].leftele - 1].topo.zmax +
            elem[river[i].rightele - 1].topo.zmax) -
            0.5 * (elem[river[i].leftele - 1].topo.zmin +
            elem[river[i].rightele - 1].topo.zmin));
        river[i].topo.node_zmax = meshtbl->zmax[river[i].tonode - 1];
        river[i].topo.dist_left = sqrt(
            (river[i].topo.x - elem[river[i].leftele - 1].topo.x) *
            (river[i].topo.x - elem[river[i].leftele - 1].topo.x) +
            (river[i].topo.y - elem[river[i].leftele - 1].topo.y) *
            (river[i].topo.y - elem[river[i].leftele - 1].topo.y));
        river[i].topo.dist_right = sqrt(
            (river[i].topo.x - elem[river[i].rightele - 1].topo.x) *
            (river[i].topo.x - elem[river[i].rightele - 1].topo.x) +
            (river[i].topo.y - elem[river[i].rightele - 1].topo.y) *
            (river[i].topo.y - elem[river[i].rightele - 1].topo.y));

        river[i].shp.depth = cal->rivdepth * shptbl->depth[rivtbl->shp[i] - 1];
        river[i].shp.intrpl_ord = shptbl->intrpl_ord[rivtbl->shp[i] - 1];
        river[i].shp.coeff =
            cal->rivshpcoeff * shptbl->coeff[rivtbl->shp[i] - 1];
        river[i].shp.length = sqrt(
            pow(meshtbl->x[river[i].fromnode - 1] -
            meshtbl->x[river[i].tonode - 1], 2) +
            pow(meshtbl->y[river[i].fromnode - 1] -
            meshtbl->y[river[i].tonode - 1], 2));
        river[i].shp.width = RiverEqWid(river[i].shp.intrpl_ord,
            river[i].shp.depth, river[i].shp.coeff);

        river[i].topo.zbed = river[i].topo.zmax - river[i].shp.depth;

        river[i].matl.rough =
            cal->rivrough * matltbl->rough[rivtbl->matl[i] - 1];
        river[i].matl.cwr = matltbl->cwr[rivtbl->matl[i] - 1];
        river[i].matl.ksath =
            cal->rivksath * matltbl->ksath[rivtbl->matl[i] - 1];
        river[i].matl.ksatv =
            cal->rivksatv * matltbl->ksatv[rivtbl->matl[i] - 1];
        river[i].matl.bedthick =
            cal->rivbedthick * matltbl->bedthick[rivtbl->matl[i] - 1];
        river[i].matl.porosity = 0.5 *
            (elem[river[i].leftele - 1].soil.porosity +
            elem[river[i].rightele - 1].soil.porosity);
        river[i].matl.smcmin = 0.5 *
            (elem[river[i].leftele - 1].soil.smcmin +
            elem[river[i].rightele - 1].soil.smcmin);

        river[i].topo.area = river[i].shp.length *
            RiverEqWid(river[i].shp.intrpl_ord, river[i].shp.depth,
            river[i].shp.coeff);
    }
}

double RiverEqWid(int order, double depth, double coeff)
{
    double          eq_wid = 0.0;

    depth = (depth > 0.0) ? depth : 0.0;

    switch (order)
    {
        case RECTANGLE:
            eq_wid = coeff;
            break;
        case TRIANGLE:
        case QUADRATIC:
        case CUBIC:
            eq_wid = 2.0 *
                pow(depth + RIVDPTHMIN, 1.0 / ((double)order - 1.0)) /
                pow(coeff, 1.0 / ((double)order - 1.0));
            break;
        default:
            PIHMprintf(VL_ERROR, "Error: River order %d is not defined.\n",
                order);
            PIHMexit(EXIT_FAILURE);
    }

    return eq_wid;
}
