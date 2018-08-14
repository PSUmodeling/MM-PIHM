#include "pihm.h"

void NTransport(elem_struct *elem, river_struct *river)
{
    int             i;
    /*
     * Calculate solute N concentrations
     */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        /* Initialize N fluxes */
        for (j = 0; j < NUM_EDGE; j++)
        {
            elem[i].no3sol.flux[j] = 0.0;
            elem[i].nh4sol.flux[j] = 0.0;
        }

        /* Calculate NO3 and NH4 average concentrations in saturated zone */
        elem[i].no3sol.conc = AvgSolConc(elem[i].ps.nsoil, 0.0, elem[i].soil.bd,
            elem[i].ps.sldpth, elem[i].ws.smc, elem[i].ws.gw, elem[i].ns.no3);

        elem[i].nh4sol.conc = AvgSolConc(elem[i].ps.nsoil, 5.6, elem[i].soil.bd,
            elem[i].ps.sldpth, elem[i].ws.smc, elem[i].ws.gw, elem[i].ns.no3);
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             j;
        double          strg;

        /* Initialize N fluxes */
        for (j = 0; j < NUM_RIVFLX; j++)
        {
            river[i].no3sol.flux[j] = 0.0;
            river[i].nh4sol.flux[j] = 0.0;
        }

        /* River stream */
        if (river[i].ws.stage > 0.0)
        {
            river[i].no3sol.conc_stream = (river[i].ns.streamno3 > 0.0) ?
                LinearEquilibriumConcentration(0.0, 0.0, river[i].ws.stage, 1.0,
                river[i].ns.streamno3) : 0.0;

            river[i].nh4sol.conc_stream = (river[i].ns.streamnh4 > 0.0) ?
                LinearEquilibriumConcentration(5.6, 0.0, river[i].ws.stage, 1.0,
                river[i].ns.streamnh4) : 0.0;
        }
        else
        {
            river[i].no3sol.conc_stream = 0.0;
            river[i].nh4sol.conc_stream = 0.0;
        }

        /* River bed */
        if (river[i].ws.gw > 0.0)
        {
            river[i].no3sol.conc_bed = (river[i].ns.bedno3 > 0.0) ?
                LinearEquilibriumConcentration(0.0, river[i].matl.bd,
                river[i].ws.gw, river[i].matl.porosity, river[i].ns.bedno3) :
                0.0;

            river[i].nh4sol.conc_bed = (river[i].ns.bednh4 > 0.0) ?
                LinearEquilibriumConcentration(0.0, river[i].matl.bd,
                river[i].ws.gw, river[i].matl.porosity, river[i].ns.bednh4) :
                0.0;
        }
        else
        {
            river[i].no3sol.conc_bed = 0.0;
            river[i].nh4sol.conc_bed = 0.0;
        }
    }

    /*
     * Calculate solute fluxes
     */
//#if defined(_OPENMP)
//# pragma omp parallel for
//#endif
//    for (i = 0; i < nelem; i++)
//    {
//        elem_struct    *nabr;
//        int             j;
//
//        /* Element to element */
//        for (j = 0; j < NUM_EDGE; j++)
//        {
//            if (elem[i].nabr[j] > 0)
//            {
//                nabr = &elem[elem[i].nabr[j] - 1];
//
//                elem[i].nsol.subflux[j] = elem[i].wf.subsurf[j] * 1000.0 *
//                    ((elem[i].wf.subsurf[j] > 0.0) ?
//                    MOBILEN_PROPORTION * elem[i].nsol.conc_subsurf :
//                    MOBILEN_PROPORTION * nabr->nsol.conc_subsurf);
//            }
//            else if (elem[i].nabr[j] < 0)
//            {
//                /* Do nothing. River-element interactions are calculated
//                 * later */
//            }
//            else    /* Boundary condition flux */
//            {
//                elem[i].nsol.subflux[j] = 0.0;
//            }
//        }
//    }
//
//#if defined(_OPENMP)
//# pragma omp parallel for
//#endif
//    for (i = 0; i < nriver; i++)
//    {
//        river_struct   *down;
//        elem_struct    *left;
//        elem_struct    *right;
//        int             j;
//
//        /* Downstream and upstream */
//        if (river[i].down > 0)
//        {
//            down = &river[river[i].down - 1];
//
//            /* Stream */
//            river[i].nsol.flux[DOWN_CHANL2CHANL] =
//                river[i].wf.rivflow[DOWN_CHANL2CHANL] * 1000.0 *
//                ((river[i].wf.rivflow[DOWN_CHANL2CHANL] > 0.0) ?
//                river[i].nsol.conc_stream : down->nsol.conc_stream);
//
//            /* Bed */
//            river[i].nsol.flux[DOWN_AQUIF2AQUIF] =
//                river[i].wf.rivflow[DOWN_AQUIF2AQUIF] * 1000.0 *
//                ((river[i].wf.rivflow[DOWN_AQUIF2AQUIF] > 0.0) ?
//                MOBILEN_PROPORTION * river[i].nsol.conc_bed :
//                MOBILEN_PROPORTION * down->nsol.conc_bed);
//        }
//        else
//        {
//            river[i].nsol.flux[DOWN_CHANL2CHANL] =
//                river[i].wf.rivflow[DOWN_CHANL2CHANL] * 1000.0 *
//                river[i].nsol.conc_stream;
//
//            river[i].nsol.flux[DOWN_AQUIF2AQUIF] = 0.0;
//        }
//
//        /* Left and right banks */
//        left = &elem[river[i].leftele - 1];
//        right = &elem[river[i].rightele - 1];
//
//        if (river[i].leftele > 0)
//        {
//            river[i].nsol.flux[LEFT_SURF2CHANL] =
//                river[i].wf.rivflow[LEFT_SURF2CHANL] * 1000.0 *
//                ((river[i].wf.rivflow[LEFT_SURF2CHANL] > 0.0) ?
//                river[i].nsol.conc_stream :
//                MOBILEN_PROPORTION * left->nsol.conc_surf);
//
//            river[i].nsol.flux[LEFT_AQUIF2CHANL] =
//                river[i].wf.rivflow[LEFT_AQUIF2CHANL] * 1000.0 *
//                ((river[i].wf.rivflow[LEFT_AQUIF2CHANL] > 0.0) ?
//                river[i].nsol.conc_stream :
//                MOBILEN_PROPORTION * left->nsol.conc_subsurf);
//
//            river[i].nsol.flux[LEFT_AQUIF2AQUIF] =
//                river[i].wf.rivflow[LEFT_AQUIF2AQUIF] * 1000.0 *
//                ((river[i].wf.rivflow[LEFT_AQUIF2AQUIF] > 0.0) ?
//                MOBILEN_PROPORTION * river[i].nsol.conc_bed :
//                MOBILEN_PROPORTION * left->nsol.conc_subsurf);
//
//            for (j = 0; j < NUM_EDGE; j++)
//            {
//                if (left->nabr[j] == -(i + 1))
//                {
//                    left->nsol.ovlflux[j] =
//                        -river[i].nsol.flux[LEFT_SURF2CHANL];
//                    left->nsol.subflux[j] =
//                        -(river[i].nsol.flux[LEFT_AQUIF2CHANL] +
//                        river[i].nsol.flux[LEFT_AQUIF2AQUIF]);
//                    break;
//                }
//            }
//
//        }
//
//        if (river[i].rightele > 0)
//        {
//            river[i].nsol.flux[RIGHT_SURF2CHANL] =
//                river[i].wf.rivflow[RIGHT_SURF2CHANL] * 1000.0 *
//                ((river[i].wf.rivflow[RIGHT_SURF2CHANL] > 0.0) ?
//                river[i].nsol.conc_stream :
//                MOBILEN_PROPORTION * right->nsol.conc_surf);
//
//            river[i].nsol.flux[RIGHT_AQUIF2CHANL] =
//                river[i].wf.rivflow[RIGHT_AQUIF2CHANL] * 1000.0 *
//                ((river[i].wf.rivflow[RIGHT_AQUIF2CHANL] > 0.0) ?
//                river[i].nsol.conc_stream :
//                MOBILEN_PROPORTION * right->nsol.conc_subsurf);
//
//            river[i].nsol.flux[RIGHT_AQUIF2AQUIF] =
//                river[i].wf.rivflow[RIGHT_AQUIF2AQUIF] * 1000.0 *
//                ((river[i].wf.rivflow[RIGHT_AQUIF2AQUIF] > 0.0) ?
//                MOBILEN_PROPORTION * river[i].nsol.conc_bed :
//                MOBILEN_PROPORTION * right->nsol.conc_subsurf);
//
//            for (j = 0; j < NUM_EDGE; j++)
//            {
//                if (right->nabr[j] == -(i + 1))
//                {
//                    right->nsol.ovlflux[j] =
//                        -river[i].nsol.flux[RIGHT_SURF2CHANL];
//                    right->nsol.subflux[j] =
//                        -(river[i].nsol.flux[RIGHT_AQUIF2CHANL] +
//                        river[i].nsol.flux[RIGHT_AQUIF2AQUIF]);
//                    break;
//                }
//            }
//        }
//
//        river[i].nsol.flux[CHANL_LKG] =
//            river[i].wf.rivflow[CHANL_LKG] * 1000.0 *
//            ((river[i].wf.rivflow[CHANL_LKG] > 0.0) ?
//            river[i].nsol.conc_stream :
//            MOBILEN_PROPORTION * river[i].nsol.conc_bed);
//    }
//
//    /*
//     * Accumulate to get in-flow for down segments
//     */
//    for (i = 0; i < nriver; i++)
//    {
//        river_struct   *down;
//
//        if (river[i].down > 0)
//        {
//            down = &river[river[i].down - 1];
//
//            down->nsol.flux[UP_CHANL2CHANL] -=
//                river[i].nsol.flux[DOWN_CHANL2CHANL];
//
//            down->nsol.flux[UP_AQUIF2AQUIF] -=
//                river[i].nsol.flux[DOWN_AQUIF2AQUIF];
//        }
//    }
}

double AvgSolConc(int nsoil, double kd, const double bd[],
    const double sldpth[], const double swc[], double gw,
    const double solute[])
{
    int             k;
    double          conc_lyr;
    double          satdpth[MAXLYR];
    double          avg_conc = 0.0;
    double          sattot = 0.0;

    FindWaterTable(sldpth, nsoil, gw, satdpth);

    for (k = 0; k < nsoil; k++)
    {
        if (satdpth[k] > 0.0)
        {
            conc_lyr = (solute[k] > 0.0) ?
                LinearEquilibriumConcentration(kd, bd[k], sldpth[k], swc[k],
                    solute[k]) : 0.0;

            avg_conc += satdpth[k] * conc_lyr;
            sattot += satdpth[k];
        }
    }

    avg_conc /= sattot;

    return avg_conc;
}
