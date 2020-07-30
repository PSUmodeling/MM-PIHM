#include "pihm.h"

void Hydrol(const ctrl_struct *ctrl, elem_struct elem[], river_struct river[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        /* Calculate actual surface water depth */
        elem[i].ws.surfh = SurfH(elem[i].ws.surf);
    }

    /* Determine which layers does ET extract water from */
    EtUptake(elem);

    /* Water flow */
    LateralFlow(ctrl->surf_mode, river, elem);

    VerticalFlow((double)ctrl->stepsize, elem);

    RiverFlow(ctrl->surf_mode, ctrl->riv_mode, elem, river);
}

void EtUptake(elem_struct elem[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        /* Source of direct evaporation */
#if defined(_NOAH_)
        if (elem[i].ws.gw > elem[i].soil.depth - elem[i].soil.dinf)
        {
            elem[i].wf.edir_surf = 0.0;
            elem[i].wf.edir_unsat = 0.0;
            elem[i].wf.edir_gw = elem[i].wf.edir;
        }
        else
        {
            elem[i].wf.edir_surf = 0.0;
            elem[i].wf.edir_unsat = elem[i].wf.edir;
            elem[i].wf.edir_gw = 0.0;
        }
#else
        if (elem[i].ws.surfh >= DEPRSTG)
        {
            elem[i].wf.edir_surf = elem[i].wf.edir;
            elem[i].wf.edir_unsat = 0.0;
            elem[i].wf.edir_gw = 0.0;
        }
        else if (elem[i].ws.gw > elem[i].soil.depth - elem[i].soil.dinf)
        {
            elem[i].wf.edir_surf = 0.0;
            elem[i].wf.edir_unsat = 0.0;
            elem[i].wf.edir_gw = elem[i].wf.edir;
        }
        else
        {
            elem[i].wf.edir_surf = 0.0;
            elem[i].wf.edir_unsat = elem[i].wf.edir;
            elem[i].wf.edir_gw = 0.0;
        }
#endif

        /* Source of transpiration */
#if defined(_NOAH_)
        elem[i].ps.gwet = GwTransp(elem[i].wf.ett, elem[i].wf.et,
            elem[i].ps.nwtbl, elem[i].ps.nroot);
        elem[i].wf.ett_unsat = (1.0 - elem[i].ps.gwet) * elem[i].wf.ett;
        elem[i].wf.ett_gw = elem[i].ps.gwet * elem[i].wf.ett;
#else
        if (elem[i].ws.gw > elem[i].soil.depth - elem[i].ps.rzd)
        {
            elem[i].wf.ett_unsat = 0.0;
            elem[i].wf.ett_gw = elem[i].wf.ett;
        }
        else
        {
            elem[i].wf.ett_unsat = elem[i].wf.ett;
            elem[i].wf.ett_gw = 0.0;
        }
#endif
    }
}

double SurfH(double surf_eqv)
{
    /*
     * Following Panday and Huyakorn (2004) AWR:
     * Use a parabolic curve to express the equivalent surface water depth
     * (surf_eqv) in terms of actual flow depth (surfh) when the actual flow
     * depth is below depression storage; assume that
     * d(surf_eqv) / d(surfh) = 1.0 when surfh = DEPRSTG. Thus
     *   surf_eqv = (1 / 2 * DEPRSTG) * surfh ^ 2, i.e.
     *   surfh = sqrt(2 * DEPRSTG * surf_eqv)
     */
    double          surfh;

    if (DEPRSTG == 0.0)
    {
        surfh = surf_eqv;
    }
    else
    {
        if (surf_eqv < 0.0)
        {
            surfh = 0.0;
        }
        else if (surf_eqv <= 0.5 * DEPRSTG)
        {
            surfh = sqrt(2.0 * DEPRSTG * surf_eqv);
        }
        else
        {
            surfh = DEPRSTG + (surf_eqv - 0.5 * DEPRSTG);
        }
    }

    return surfh;
}
