#include "pihm.h"

void Hydrol(const ctrl_struct *ctrl, elem_struct elem[], river_struct river[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        /*
         * Calculate actual surface water depth
         */
        elem[i].ws.surfh = SurfH(elem[i].ws.surf);

        /*
         * Determine which layers does ET extract water from
         */
        EtExtract(&elem[i].soil, &elem[i].ws, &elem[i].ps, &elem[i].wf);
    }


    /*
     * Water flow
     */
    LateralFlow(elem, river, ctrl->surf_mode);

    VerticalFlow(elem, (double)ctrl->stepsize);

    RiverFlow(elem, river, ctrl->riv_mode);
}

void EtExtract(const soil_struct *soil, const wstate_struct *ws,
    pstate_struct *ps, wflux_struct *wf)
{
    /*
     * Source of direct evaporation
     */
#if defined(_NOAH_)
    if (ws->gw > soil->depth - soil->dinf)
    {
        wf->edir_surf = 0.0;
        wf->edir_unsat = 0.0;
        wf->edir_gw = wf->edir;
    }
    else
    {
        wf->edir_surf = 0.0;
        wf->edir_unsat = wf->edir;
        wf->edir_gw = 0.0;
    }
#else
    if (ws->surfh >= DEPRSTG)
    {
        wf->edir_surf = wf->edir;
        wf->edir_unsat = 0.0;
        wf->edir_gw = 0.0;
    }
    else if (ws->gw > soil->depth - soil->dinf)
    {
        wf->edir_surf = 0.0;
        wf->edir_unsat = 0.0;
        wf->edir_gw = wf->edir;
    }
    else
    {
        wf->edir_surf = 0.0;
        wf->edir_unsat = wf->edir;
        wf->edir_gw = 0.0;
    }
#endif

    /*
     * Source of transpiration
     */
#if defined(_NOAH_)
    ps->gwet = GwTransp(wf->ett, wf->et, ps->nwtbl, ps->nroot);
    wf->ett_unsat = (1.0 - ps->gwet) * wf->ett;
    wf->ett_gw = ps->gwet * wf->ett;
#else
    if (ws->gw > soil->depth - ps->rzd)
    {
        wf->ett_unsat = 0.0;
        wf->ett_gw = wf->ett;
    }
    else
    {
        wf->ett_unsat = wf->ett;
        wf->ett_gw = 0.0;
    }
#endif
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
        surfh = surfeqv;
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
    }

    return surfh;
}
