#include "pihm.h"

void InitRiver(const meshtbl_struct *meshtbl, const rivtbl_struct *rivtbl,
    const shptbl_struct *shptbl, const matltbl_struct *matltbl,
    const calib_struct *calib, elem_struct elem[], river_struct river[])
{
    int             i;

    for (i = 0; i < nriver; i++)
    {
        int             j;

        river[i].ind = i + 1;
        river[i].left = rivtbl->left[i];
        river[i].right = rivtbl->right[i];
        river[i].from = rivtbl->from[i];
        river[i].to = rivtbl->to[i];
        river[i].down = rivtbl->down[i];

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[river[i].left - 1].nabr[j] == river[i].right)
            {
                elem[river[i].left - 1].nabr_river[j] = i + 1;
            }
            if (elem[river[i].right - 1].nabr[j] == river[i].left)
            {
                elem[river[i].right - 1].nabr_river[j] = i + 1;
            }
        }

        river[i].topo.x = 0.5 *
            (meshtbl->x[river[i].from - 1] + meshtbl->x[river[i].to - 1]);
        river[i].topo.y = 0.5 *
            (meshtbl->y[river[i].from - 1] + meshtbl->y[river[i].to - 1]);
        river[i].topo.zmax = 0.5 *
            (meshtbl->zmax[river[i].from - 1] + meshtbl->zmax[river[i].to - 1]);
        river[i].topo.node_zmax = meshtbl->zmax[river[i].to - 1];
        river[i].topo.dist_left = sqrt(
            (river[i].topo.x - elem[river[i].left - 1].topo.x) *
            (river[i].topo.x - elem[river[i].left - 1].topo.x) +
            (river[i].topo.y - elem[river[i].left - 1].topo.y) *
            (river[i].topo.y - elem[river[i].left - 1].topo.y));
        river[i].topo.dist_right = sqrt(
            (river[i].topo.x - elem[river[i].right - 1].topo.x) *
            (river[i].topo.x - elem[river[i].right - 1].topo.x) +
            (river[i].topo.y - elem[river[i].right - 1].topo.y) *
            (river[i].topo.y - elem[river[i].right - 1].topo.y));

        river[i].shp.depth = calib->rivdepth *
            shptbl->depth[rivtbl->shp[i] - 1];
        river[i].shp.intrpl_ord = shptbl->intrpl_ord[rivtbl->shp[i] - 1];
        river[i].shp.coeff = calib->rivshpcoeff *
            shptbl->coeff[rivtbl->shp[i] - 1];
        river[i].shp.length = sqrt(
            (meshtbl->x[river[i].from - 1] - meshtbl->x[river[i].to - 1]) *
            (meshtbl->x[river[i].from - 1] - meshtbl->x[river[i].to - 1]) +
            (meshtbl->y[river[i].from - 1] - meshtbl->y[river[i].to - 1]) *
            (meshtbl->y[river[i].from - 1] - meshtbl->y[river[i].to - 1]));
        river[i].shp.width = RiverEqWid(river[i].shp.intrpl_ord,
            river[i].shp.depth, river[i].shp.coeff);

        river[i].topo.zbed = river[i].topo.zmax - river[i].shp.depth;

        river[i].matl.rough = calib->rivrough *
            matltbl->rough[rivtbl->matl[i] - 1];
        river[i].matl.cwr = matltbl->cwr[rivtbl->matl[i] - 1];
        river[i].matl.ksath = calib->rivksath *
            matltbl->ksath[rivtbl->matl[i] - 1];

        river[i].topo.area = river[i].shp.length *
            RiverEqWid(river[i].shp.intrpl_ord, river[i].shp.depth,
            river[i].shp.coeff);
    }
}

double RiverEqWid(int order, double depth, double coeff)
{
    double          eq_wid = 0.0;

    depth = MAX(depth, 0.0);

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
            pihm_printf(VL_ERROR, "Error: River order %d is not defined.\n",
                order);
            pihm_exit(EXIT_FAILURE);
    }

    return eq_wid;
}
