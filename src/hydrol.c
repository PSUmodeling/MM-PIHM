#include "pihm.h"

int Hydrol (realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *pihm_data)
{
    int             i, j;
    double         *y;
    double         *dy;
    double          dt;
    pihm_struct     pihm;
    elem_struct    *elem;
    river_struct   *riv;

    y = NV_DATA_S (CV_Y);
    dy = NV_DATA_S (CV_Ydot);
    pihm = (pihm_struct)pihm_data;

    dt = (double)pihm->ctrl.stepsize;

    /*
     * Initialization of temporary state variables 
     */
    for (i = 0; i < 3 * pihm->numele + 2 * pihm->numriv; i++)
    {
        dy[i] = 0.0;
    }

    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        elem->ws.surf = (y[SURF (i)] >= 0.0) ? y[SURF (i)] : 0.0;
        elem->ws.unsat = (y[UNSAT (i)] >= 0.0) ? y[UNSAT (i)] : 0.0;
        elem->ws.gw = (y[GW (i)] >= 0.0) ? y[GW (i)] : 0.0;
    }

    for (i = 0; i < pihm->numriv; i++)
    {
        riv = &pihm->riv[i];

        riv->ws.stage = (y[RIVSTG (i)] >= 0.0) ? y[RIVSTG (i)] : 0.0;
        riv->ws.gw = (y[RIVGW (i)] >= 0.0) ? y[RIVGW (i)] : 0.0;

        riv->wf.river[UP_CHANL2CHANL] = 0.0;
        riv->wf.river[UP_AQUIF2AQUIF] = 0.0;
    }

    /*
     * Determine source of ET
     */
    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        /* Source of direct evaporation */
#ifdef _NOAH_
        if (elem->ws.gw > elem->soil.depth - elem->soil.dinf)
        {
            elem->wf.edir_surf = 0.0;
            elem->wf.edir_unsat = 0.0;
            elem->wf.edir_gw = elem->wf.edir;
        }
        else
        {
            elem->wf.edir_surf = 0.0;
            elem->wf.edir_unsat = elem->wf.edir;
            elem->wf.edir_gw = 0.0;
        }
#else
        if (elem->ws.surf >= DEPRSTG)
        {
            elem->wf.edir_surf = elem->wf.edir;
            elem->wf.edir_unsat = 0.0;
            elem->wf.edir_gw = 0.0;
        }
        else if (elem->ws.gw > elem->soil.depth - elem->soil.dinf)
        {
            elem->wf.edir_surf = 0.0;
            elem->wf.edir_unsat = 0.0;
            elem->wf.edir_gw = elem->wf.edir;
        }
        else
        {
            elem->wf.edir_surf = 0.0;
            elem->wf.edir_unsat = elem->wf.edir;
            elem->wf.edir_gw = 0.0;
        }
#endif

        /* Source of transpiration */
#ifdef _NOAH_
        elem->ps.gwet = GWTransp (elem->wf.ett, elem->wf.et, elem->ps.nwtbl,
            elem->ps.nroot);
        elem->wf.ett_unsat = (1.0 - elem->ps.gwet) * elem->wf.ett;
        elem->wf.ett_gw = elem->ps.gwet * elem->wf.ett;
#else
        if (elem->ws.gw > elem->soil.depth - elem->ps.rzd)
        {
            elem->wf.ett_unsat = 0.0;
            elem->wf.ett_gw = elem->wf.ett;
        }
        else
        {
            elem->wf.ett_unsat = elem->wf.ett;
            elem->wf.ett_gw = 0.0;
        }
#endif
    }

    /*
     * Water flow
     */
    VerticalFlow (pihm);

    LateralFlow (pihm);

    RiverFlow (pihm);

    /*
     * RHS of ODEs
     */
    for (i = 0; i < pihm->numele; i++)
    {
        elem = &pihm->elem[i];

        dy[SURF (i)] += elem->wf.pcpdrp - elem->wf.infil - elem->wf.edir_surf;
        dy[UNSAT (i)] += elem->wf.infil - elem->wf.rechg -
            elem->wf.edir_unsat - elem->wf.ett_unsat;
        dy[GW (i)] += elem->wf.rechg - elem->wf.edir_gw - elem->wf.ett_gw;

        for (j = 0; j < 3; j++)
        {
            dy[SURF (i)] -= elem->wf.surf[j] / elem->topo.area;
            dy[GW (i)] -= elem->wf.subsurf[j] / elem->topo.area;
        }

        dy[UNSAT (i)] /= elem->soil.porosity;
        dy[GW (i)] /= elem->soil.porosity;

        if (isnan (dy[SURF (i)]))
        {
            printf
                ("ERROR: NAN error for Element %d (surface water) at %lf\n",
                i + 1, t);
            PihmExit (1);
        }
        if (isnan (dy[UNSAT (i)]))
        {
            printf ("ERROR: NAN error for Element %d (unsat water) at %lf\n",
                i + 1, t);
            PihmExit (1);
        }
        if (isnan (dy[GW (i)]))
        {
            printf ("ERROR: NAN error for Element %d (groundwater) at %lf\n",
                i + 1, t);
            PihmExit (1);
        }
    }

    for (i = 0; i < pihm->numriv; i++)
    {
        riv = &(pihm->riv[i]);

        for (j = 0; j <= 6; j++)
        {
            /* Note the limitation due to
             * d(v)/dt=a*dy/dt+y*da/dt
             * for cs other than rectangle */
            dy[RIVSTG (i)] -= riv->wf.river[j] / riv->topo.area;
        }
        dy[RIVGW (i)] +=
            -riv->wf.river[LEFT_AQUIF2AQUIF] - riv->wf.river[RIGHT_AQUIF2AQUIF] - riv->wf.river[DOWN_AQUIF2AQUIF] -
            riv->wf.river[UP_AQUIF2AQUIF] + riv->wf.river[CHANL_LKG];
        dy[RIVGW (i)] /= riv->matl.porosity * riv->topo.area;

        if (isnan (dy[RIVSTG (i)]))
        {
            printf ("ERROR: NAN error for River Segment %d (river) at %lf\n",
                i + 1, t);
            PihmExit (1);
        }
        if (isnan (dy[RIVGW (i)]))
        {
            printf
                ("ERROR: NAN error for River Segment %d (groundwater) at"
                "%lf\n", i + 1, t);
            PihmExit (1);
        }
    }

    return (0);
}
