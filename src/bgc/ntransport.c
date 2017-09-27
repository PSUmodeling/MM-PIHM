#include "pihm.h"

void NTransport(elem_struct *elem, river_struct *riv)
{
    int             i;

    /*
     * Calculate solute N concentrations
     */
#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;
        double          strg;

        /* Initialize N fluxes */
        elem[i].nsol.infilflux = 0.0;
        for (j = 0; j < NUM_EDGE; j++)
        {
            elem[i].nsol.ovlflux[j] = 0.0;
            elem[i].nsol.subflux[j] = 0.0;
        }

        /* Element surface */
        strg = elem[i].ws.surf;
        elem[i].nsol.conc_surf = (strg > 0.0) ?
            elem[i].ns.surfn / strg / 1000.0 : 0.0;
        elem[i].nsol.conc_surf = (elem[i].nsol.conc_surf > 0.0) ?
            elem[i].nsol.conc_surf : 0.0;

        /* Element subsurface */
        strg = (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity +
            elem[i].soil.depth * elem[i].soil.smcmin;
        elem[i].nsol.conc_subsurf = (strg > 0.0) ?
            elem[i].ns.sminn / strg / 1000.0 : 0.0;
        elem[i].nsol.conc_subsurf = (elem[i].nsol.conc_subsurf > 0.0) ?
            elem[i].nsol.conc_subsurf : 0.0;
    }

#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             j;
        double          strg;

        /* Initialize N fluxes */
        for (j = 0; j < NUM_RIVFLX; j++)
        {
            riv[i].nsol.flux[j] = 0.0;
        }

        /* River stream */
        strg = riv[i].ws.stage;
        riv[i].nsol.conc_stream = (strg > 0.0) ?
            riv[i].ns.streamn / strg / 1000.0 : 0.0;
        riv[i].nsol.conc_stream = (riv[i].nsol.conc_stream > 0.0) ?
            riv[i].nsol.conc_stream : 0.0;

        /* River bed */
        strg = riv[i].ws.gw * riv[i].matl.porosity +
            riv[i].matl.bedthick * riv[i].matl.smcmin;
        riv[i].nsol.conc_bed = (strg > 0.0) ? riv[i].ns.sminn / strg / 1000.0 : 0.0;
        riv[i].nsol.conc_bed = (riv[i].nsol.conc_bed > 0.0) ?
            riv[i].nsol.conc_bed : 0.0;
    }

    /*
     * Calculate solute fluxes
     */
#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem_struct    *nabr;
        int             j;

        /* Infiltration */
        elem[i].nsol.infilflux = elem[i].wf.infil * 1000.0 *
            ((elem[i].wf.infil > 0.0) ? elem[i].nsol.conc_surf :
            MOBILEN_PROPORTION * elem[i].nsol.conc_subsurf);

        /* Element to element */
        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].nabr[j] > 0)
            {
                nabr = &elem[elem[i].nabr[j] - 1];

                elem[i].nsol.subflux[j] = elem[i].wf.subsurf[j] * 1000.0 *
                    ((elem[i].wf.subsurf[j] > 0.0) ?
                    MOBILEN_PROPORTION * elem[i].nsol.conc_subsurf :
                    MOBILEN_PROPORTION * nabr->nsol.conc_subsurf);

                elem[i].nsol.ovlflux[j] = elem[i].wf.ovlflow[j] * 1000.0 *
                    ((elem[i].wf.ovlflow[j] > 0.0) ?
                    MOBILEN_PROPORTION * elem[i].nsol.conc_surf :
                    MOBILEN_PROPORTION * nabr->nsol.conc_surf);
            }
            else if (elem[i].nabr[j] < 0)
            {
                /* Do nothing. River-element interactions are calculated
                 * later */
            }
            else    /* Boundary condition flux */
            {
                elem[i].nsol.ovlflux[j] = 0.0;
                elem[i].nsol.subflux[j] = 0.0;
            }
        }
    }

#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river_struct   *down;
        elem_struct    *left;
        elem_struct    *right;
        int             j;

        /* Downstream and upstream */
        if (riv[i].down > 0)
        {
            down = &riv[riv[i].down - 1];

            /* Stream */
            riv[i].nsol.flux[DOWN_CHANL2CHANL] =
                riv[i].wf.rivflow[DOWN_CHANL2CHANL] * 1000.0 *
                ((riv[i].wf.rivflow[DOWN_CHANL2CHANL] > 0.0) ?
                riv[i].nsol.conc_stream : down->nsol.conc_stream);

            down->nsol.flux[UP_CHANL2CHANL] = -riv[i].nsol.flux[DOWN_CHANL2CHANL];

            /* Bed */
            riv[i].nsol.flux[DOWN_AQUIF2AQUIF] =
                riv[i].wf.rivflow[DOWN_AQUIF2AQUIF] * 1000.0 *
                ((riv[i].wf.rivflow[DOWN_AQUIF2AQUIF] > 0.0) ?
                MOBILEN_PROPORTION * riv[i].nsol.conc_bed :
                MOBILEN_PROPORTION * down->nsol.conc_bed);

            down->nsol.flux[UP_AQUIF2AQUIF] = -riv[i].nsol.flux[DOWN_AQUIF2AQUIF];
        }
        else
        {
            riv[i].nsol.flux[DOWN_CHANL2CHANL] =
                riv[i].wf.rivflow[DOWN_CHANL2CHANL] * 1000.0 *
                riv[i].nsol.conc_stream;

            riv[i].nsol.flux[DOWN_AQUIF2AQUIF] = 0.0;
        }

        /* Left and right banks */
        left = &elem[riv[i].leftele - 1];
        right = &elem[riv[i].rightele - 1];

        if (riv[i].leftele > 0)
        {
            riv[i].nsol.flux[LEFT_SURF2CHANL] =
                riv[i].wf.rivflow[LEFT_SURF2CHANL] * 1000.0 *
                ((riv[i].wf.rivflow[LEFT_SURF2CHANL] > 0.0) ?
                riv[i].nsol.conc_stream :
                MOBILEN_PROPORTION * left->nsol.conc_surf);

            riv[i].nsol.flux[LEFT_AQUIF2CHANL] =
                riv[i].wf.rivflow[LEFT_AQUIF2CHANL] * 1000.0 *
                ((riv[i].wf.rivflow[LEFT_AQUIF2CHANL] > 0.0) ?
                riv[i].nsol.conc_stream :
                MOBILEN_PROPORTION * left->nsol.conc_subsurf);

            riv[i].nsol.flux[LEFT_AQUIF2AQUIF] =
                riv[i].wf.rivflow[LEFT_AQUIF2AQUIF] * 1000.0 *
                ((riv[i].wf.rivflow[LEFT_AQUIF2AQUIF] > 0.0) ?
                MOBILEN_PROPORTION * riv[i].nsol.conc_bed :
                MOBILEN_PROPORTION * left->nsol.conc_subsurf);

            for (j = 0; j < NUM_EDGE; j++)
            {
                if (left->nabr[j] == -(i + 1))
                {
                    left->nsol.ovlflux[j] = -riv[i].nsol.flux[LEFT_SURF2CHANL];
                    left->nsol.subflux[j] =
                        -(riv[i].nsol.flux[LEFT_AQUIF2CHANL] +
                        riv[i].nsol.flux[LEFT_AQUIF2AQUIF]);
                    break;
                }
            }

        }

        if (riv[i].rightele > 0)
        {
            riv[i].nsol.flux[RIGHT_SURF2CHANL] =
                riv[i].wf.rivflow[RIGHT_SURF2CHANL] * 1000.0 *
                ((riv[i].wf.rivflow[RIGHT_SURF2CHANL] > 0.0) ?
                riv[i].nsol.conc_stream :
                MOBILEN_PROPORTION * right->nsol.conc_surf);

            riv[i].nsol.flux[RIGHT_AQUIF2CHANL] =
                riv[i].wf.rivflow[RIGHT_AQUIF2CHANL] * 1000.0 *
                ((riv[i].wf.rivflow[RIGHT_AQUIF2CHANL] > 0.0) ?
                riv[i].nsol.conc_stream :
                MOBILEN_PROPORTION * right->nsol.conc_subsurf);

            riv[i].nsol.flux[RIGHT_AQUIF2AQUIF] =
                riv[i].wf.rivflow[RIGHT_AQUIF2AQUIF] * 1000.0 *
                ((riv[i].wf.rivflow[RIGHT_AQUIF2AQUIF] > 0.0) ?
                MOBILEN_PROPORTION * riv[i].nsol.conc_bed :
                MOBILEN_PROPORTION * right->nsol.conc_subsurf);

            for (j = 0; j < NUM_EDGE; j++)
            {
                if (right->nabr[j] == -(i + 1))
                {
                    right->nsol.ovlflux[j] =
                        -riv[i].nsol.flux[RIGHT_SURF2CHANL];
                    right->nsol.subflux[j] =
                        -(riv[i].nsol.flux[RIGHT_AQUIF2CHANL] +
                        riv[i].nsol.flux[RIGHT_AQUIF2AQUIF]);
                    break;
                }
            }
        }

        riv[i].nsol.flux[CHANL_LKG] = riv[i].wf.rivflow[CHANL_LKG] * 1000.0 *
            ((riv[i].wf.rivflow[CHANL_LKG] > 0.0) ?
            riv[i].nsol.conc_stream :
            MOBILEN_PROPORTION * riv[i].nsol.conc_bed);
    }
}
