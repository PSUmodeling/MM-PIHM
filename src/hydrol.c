#include "pihm.h"

void Hydrol(elem_struct *elem, river_struct *river, const ctrl_struct *ctrl)
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
    EtExtract(elem);

    /* Water flow */
    LateralFlow(elem, river, ctrl->surf_mode);

    VerticalFlow(elem, (double)ctrl->stepsize);

    RiverFlow(elem, river, ctrl->riv_mode);
}

void EtExtract(elem_struct *elem)
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

double SurfH(double surfeqv)
{
    /*
     * Following Panday and Huyakorn (2004) AWR:
     * Use a parabolic curve to express the equivalent surface water depth
     * (surfeqv) in terms of actual flow depth (surfh) when the actual flow
     * depth is below depression storage; assume that
     * d(surfeqv) / d(surfh) = 1.0 when surfh = DEPRSTG. Thus
     *   surfeqv = (1 / 2 * DEPRSTG) * surfh ^ 2, i.e.
     *   surfh = sqrt(2 * DEPRSTG * surfeqv)
     */
    double          surfh;

    if (DEPRSTG == 0.0)
    {
        return surfeqv;
    }
    else
    {
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
}
