#include "pihm.h"

void NTransport (pihm_struct pihm)
{
    int             i;

    /*
     * Calculate solute N concentrantions
     */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem_struct *elem;
        double      strg;

        elem = &pihm->elem[i];

        strg = elem->ws.surf +
            (elem->ws.unsat + elem->ws.gw) * elem->soil.porosity +
            elem->soil.depth * elem->soil.smcmin;
        elem->nsol.conc = (strg > 0.0) ?
            MOBILEN_PROPORTION * elem->ns.sminn / strg / 1000.0 : 0.0;
        elem->nsol.conc = (elem->nsol.conc > 0.0) ?
            elem->nsol.conc : 0.0;
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river_struct *riv;
        double      strg;

        riv = &pihm->riv[i];

        /* River stream */
        strg = riv->ws.stage + riv->ws.gw * riv->matl.porosity +
            riv->matl.bedthick * riv->matl.smcmin;
        riv->nsol.conc = (strg > 0.0) ?
            MOBILEN_PROPORTION * riv->ns.rivern / strg / 1000.0 : 0.0;
        riv->nsol.conc = (riv->nsol.conc > 0.0) ?
            riv->nsol.conc : 0.0;
    }

    /*
     * Calculate solute fluxes
     */
    for (i = 0; i < nelem; i++)
    {
        elem_struct *elem;
        elem_struct *nabr;
        double      wflux;
        int         j;

        elem = &pihm->elem[i];

        /* Element to element */
        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem->nabr[j] > 0)
            {
                nabr = &pihm->elem[elem->nabr[j] - 1];

                wflux = elem->wf.ovlflow[j] + elem->wf.subsurf[j];

                elem->nsol.flux[j] = wflux * 1000.0 *
                    ((wflux > 0.0) ?  elem->nsol.conc : nabr->nsol.conc);
            }
            else if (elem->nabr[j] < 0)
            {
                /* Do nothing. River-element interactions are calculated
                 * later */
            }
            else                        /* Boundary condition flux */
            {
                elem->nsol.flux[j] = 0.0;
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river_struct *riv;
        river_struct *down;
        elem_struct  *left;
        elem_struct  *right;
        double      wflux;
        int         j;

        riv = &pihm->riv[i];

        /*
         * Initialize solute fluxes
         */
        for (j = 0; j < 4; j++)
        {
            riv->nsol.flux[j] = 0.0;
        }

        /* Downstream and upstream */
        wflux = riv->wf.rivflow[DOWN_CHANL2CHANL] +
            riv->wf.rivflow[DOWN_AQUIF2AQUIF];

        if (riv->down > 0)
        {
            down = &pihm->riv[riv->down - 1];

            riv->nsol.flux[DOWN] = wflux * 1000.0 *
                ((wflux > 0.0) ? riv->nsol.conc : down->nsol.conc);

            /* Uptream */
            down->nsol.flux[UP] = -riv->nsol.flux[DOWN];
        }
        else
        {
            /* Outlet */
            riv->nsol.flux[DOWN] = wflux * 1000.0 * riv->nsol.conc;
        }

        /* Left and right banks */
        left = &pihm->elem[riv->leftele - 1];
        right = &pihm->elem[riv->rightele - 1];

        if (riv->leftele > 0)
        {
            wflux = riv->wf.rivflow[LEFT_SURF2CHANL] +
                riv->wf.rivflow[LEFT_AQUIF2CHANL] +
                riv->wf.rivflow[LEFT_AQUIF2AQUIF];

            riv->nsol.flux[LEFT] = wflux * 1000.0 *
                ((wflux > 0.0) ? riv->nsol.conc: left->nsol.conc);

            for (j = 0; j < NUM_EDGE; j++)
            {
                if (left->nabr[j] == - (i + 1))
                {
                    left->nsol.flux[j] = -riv->nsol.flux[LEFT];
                    break;
                }
            }

        }

        if (riv->rightele > 0)
        {
            wflux = riv->wf.rivflow[RIGHT_SURF2CHANL] +
                riv->wf.rivflow[RIGHT_AQUIF2CHANL] +
                riv->wf.rivflow[RIGHT_AQUIF2AQUIF];

            riv->nsol.flux[RIGHT] = wflux * 1000.0 *
                ((wflux > 0.0) ?  riv->nsol.conc: right->nsol.conc);

            for (j = 0; j < NUM_EDGE; j++)
            {
                if (right->nabr[j] == - (i + 1))
                {
                    right->nsol.flux[j] = -riv->nsol.flux[RIGHT];
                    break;
                }
            }
        }
    }
}
