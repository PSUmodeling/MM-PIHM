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
        strg = riv->ws.stage;
        riv->nsol.conc_stream = (strg > 0.0) ?
            riv->ns.streamn / strg / 1000.0 : 0.0;
        riv->nsol.conc_stream = (riv->nsol.conc_stream > 0.0) ?
            riv->nsol.conc_stream : 0.0;

        /* River bed */
        strg = riv->ws.gw * riv->matl.porosity;
        riv->nsol.conc_bed = (strg > 0.0) ?
            MOBILEN_PROPORTION * riv->ns.sminn / strg / 1000.0 : 0.0;
        riv->nsol.conc_bed = (riv->nsol.conc_bed > 0.0) ?
            riv->nsol.conc_bed : 0.0;
    }

    /*
     * Calculate solute fluxes
     */
    for (i = 0; i < nelem; i++)
    {
        elem_struct *elem;
        elem_struct *nabr;
        int         j;

        elem = &pihm->elem[i];

        /* Element to element */
        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem->nabr[j] > 0)
            {
                nabr = &pihm->elem[elem->nabr[j] - 1];

                elem->nsol.flux[j] =
                    (elem->wf.ovlflow[j] + elem->wf.subsurf[j]) * 1000.0 *
                    ((elem->wf.ovlflow[j] + elem->wf.subsurf[j] > 0.0) ?
                    elem->nsol.conc : nabr->nsol.conc);
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
        int         j;

        riv = &pihm->riv[i];

        /* Downstream and upstream */
        if (riv->down > 0)
        {
            down = &pihm->riv[riv->down - 1];

            /* Stream */
            riv->nsol.flux[DOWN_CHANL2CHANL] =
                riv->wf.rivflow[DOWN_CHANL2CHANL] * 1000.0 *
                ((riv->wf.rivflow[DOWN_CHANL2CHANL] > 0.0) ?
                riv->nsol.conc_stream : down->nsol.conc_stream);

            down->nsol.flux[UP_CHANL2CHANL] =
                -riv->nsol.flux[DOWN_CHANL2CHANL];

            /* Bed */
            riv->nsol.flux[DOWN_AQUIF2AQUIF] =
                riv->wf.rivflow[DOWN_AQUIF2AQUIF] * 1000.0 *
                ((riv->wf.rivflow[DOWN_AQUIF2AQUIF] > 0.0) ?
                riv->nsol.conc_bed : down->nsol.conc_bed);

            down->nsol.flux[UP_AQUIF2AQUIF] =
                -riv->nsol.flux[DOWN_AQUIF2AQUIF];
        }
        else
        {
            riv->nsol.flux[DOWN_CHANL2CHANL] =
                riv->wf.rivflow[DOWN_CHANL2CHANL] * 1000.0 *
                riv->nsol.conc_stream;

            riv->nsol.flux[DOWN_AQUIF2AQUIF] = 0.0;
        }

        /* Left and right banks */
        left = &pihm->elem[riv->leftele - 1];
        right = &pihm->elem[riv->rightele - 1];

        if (riv->leftele > 0)
        {
            riv->nsol.flux[LEFT_SURF2CHANL] =
                riv->wf.rivflow[LEFT_SURF2CHANL] * 1000.0 *
                ((riv->wf.rivflow[LEFT_SURF2CHANL] > 0.0) ?
                riv->nsol.conc_stream : left->nsol.conc);

            riv->nsol.flux[LEFT_AQUIF2CHANL] =
                riv->wf.rivflow[LEFT_AQUIF2CHANL] * 1000.0 *
                ((riv->wf.rivflow[LEFT_AQUIF2CHANL] > 0.0) ?
                riv->nsol.conc_stream : left->nsol.conc);

            riv->nsol.flux[LEFT_AQUIF2AQUIF] =
                riv->wf.rivflow[LEFT_AQUIF2AQUIF] * 1000.0 *
                ((riv->wf.rivflow[LEFT_AQUIF2AQUIF] > 0.0) ?
                riv->nsol.conc_bed : left->nsol.conc);

            for (j = 0; j < NUM_EDGE; j++)
            {
                if (left->nabr[j] == - (i + 1))
                {
                    left->nsol.flux[j] = -(riv->nsol.flux[LEFT_SURF2CHANL] +
                        riv->nsol.flux[LEFT_AQUIF2CHANL] +
                        riv->nsol.flux[LEFT_AQUIF2AQUIF]);
                    break;
                }
            }

        }

        if (riv->rightele > 0)
        {
            riv->nsol.flux[RIGHT_SURF2CHANL] =
                riv->wf.rivflow[RIGHT_SURF2CHANL] * 1000.0 *
                ((riv->wf.rivflow[RIGHT_SURF2CHANL] > 0.0) ?
                riv->nsol.conc_stream : right->nsol.conc);

            riv->nsol.flux[RIGHT_AQUIF2CHANL] =
                riv->wf.rivflow[RIGHT_AQUIF2CHANL] * 1000.0 *
                ((riv->wf.rivflow[RIGHT_AQUIF2CHANL] > 0.0) ?
                riv->nsol.conc_stream : right->nsol.conc);

            riv->nsol.flux[RIGHT_AQUIF2AQUIF] =
                riv->wf.rivflow[RIGHT_AQUIF2AQUIF] * 1000.0 *
                ((riv->wf.rivflow[RIGHT_AQUIF2AQUIF] > 0.0) ?
                riv->nsol.conc_bed : right->nsol.conc);

            for (j = 0; j < NUM_EDGE; j++)
            {
                if (right->nabr[j] == - (i + 1))
                {
                    right->nsol.flux[j] =
                        -(riv->nsol.flux[RIGHT_SURF2CHANL] +
                        riv->nsol.flux[RIGHT_AQUIF2CHANL] +
                        riv->nsol.flux[RIGHT_AQUIF2AQUIF]);
                    break;
                }
            }
        }

        riv->nsol.flux[CHANL_LKG] = riv->wf.rivflow[CHANL_LKG] * 1000.0 *
            ((riv->wf.rivflow[CHANL_LKG] > 0.0) ?
            riv->nsol.conc_stream : riv->nsol.conc_bed);
    }
}
