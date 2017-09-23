#include "pihm.h"

void Hydrol(elem_struct *elem, river_struct *riv, ctrl_struct *ctrl)
{
    int             i;

#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        /* Calculate actual surface water depth */
        elem[i].ws.surfh = SurfH(elem[i].ws.surf);
    }

    /* Determine which layers does ET extract water from */
    ETExtract(elem);

    /* Water flow */
    VerticalFlow(elem, (double)ctrl->stepsize);

    LateralFlow(elem, riv, ctrl->surf_mode);

    RiverFlow(elem, riv, ctrl->riv_mode);
}

void ETExtract(elem_struct *elem)
{
    int             i;

#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        /* Source of direct evaporation */
#ifdef _NOAH_
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
#ifdef _NOAH_
        elem[i].ps.gwet = GWTransp(elem[i].wf.ett, elem[i].wf.et, elem[i].ps.nwtbl,
            elem[i].ps.nroot);
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

double SurfH(double surfeqv)
{
    /*
     * Following Panday and Huyakorn (2004) AWR:
     * Use a parabolic curve to express the equivalent surface water depth
     * (surfeqv) in terms of actual flow depth (surfh) when the actual flow
     * depth is below depression storage; assume that
     * d(surfeqv) / d(surfh) = 1.0 when surfh = DEPRSTG
     */
    double          surfh;

    if (surfeqv < 0.0)
    {
        surfh = 0.0;
    }
    else if (surfeqv <= 0.5 * DEPRSTG)
    {
        surfh = sqrt(2.0 * DEPRSTG * surfeqv);
    }
    else
    {
        surfh = DEPRSTG + (surfeqv - 0.5 * DEPRSTG);
    }

    return surfh;
}
