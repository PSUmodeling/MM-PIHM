#include "pihm.h"

void NTransport(double dt, elem_struct elem[], river_struct river[])
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

        UpdNProf(dt, &elem[i].soil, &elem[i].ws, &elem[i].ns0,
            &elem[i].nf, &elem[i].np, &elem[i].ps, &elem[i].ns);

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
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem_struct    *nabr;
        int             j;

        /* Element to element */
        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].nabr[j] > 0)
            {
                nabr = &elem[elem[i].nabr[j] - 1];

                elem[i].no3sol.flux[j] = elem[i].wf.subsurf[j] * RHOH2O *
                    ((elem[i].wf.subsurf[j] > 0.0) ?
                    MOBILEN_PROPORTION * elem[i].no3sol.conc :
                    MOBILEN_PROPORTION * nabr->no3sol.conc);

                elem[i].nh4sol.flux[j] = elem[i].wf.subsurf[j] * RHOH2O *
                    ((elem[i].wf.subsurf[j] > 0.0) ?
                    MOBILEN_PROPORTION * elem[i].nh4sol.conc :
                    MOBILEN_PROPORTION * nabr->nh4sol.conc);
            }
            else if (elem[i].nabr[j] < 0)
            {
                /* Do nothing. River-element interactions are calculated
                 * later */
            }
            else    /* Boundary condition flux */
            {
                elem[i].no3sol.flux[j] = 0.0;
                elem[i].nh4sol.flux[j] = 0.0;
            }
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river_struct   *down;
        elem_struct    *left;
        elem_struct    *right;
        int             j;

        /* Downstream and upstream */
        if (river[i].down > 0)
        {
            down = &river[river[i].down - 1];

            /* Stream */
            river[i].no3sol.flux[DOWN_CHANL2CHANL] =
                river[i].wf.rivflow[DOWN_CHANL2CHANL] * RHOH2O *
                ((river[i].wf.rivflow[DOWN_CHANL2CHANL] > 0.0) ?
                river[i].no3sol.conc_stream : down->no3sol.conc_stream);

            river[i].nh4sol.flux[DOWN_CHANL2CHANL] =
                river[i].wf.rivflow[DOWN_CHANL2CHANL] * RHOH2O *
                ((river[i].wf.rivflow[DOWN_CHANL2CHANL] > 0.0) ?
                river[i].nh4sol.conc_stream : down->nh4sol.conc_stream);

            /* Bed */
            river[i].no3sol.flux[DOWN_AQUIF2AQUIF] =
                river[i].wf.rivflow[DOWN_AQUIF2AQUIF] * RHOH2O *
                ((river[i].wf.rivflow[DOWN_AQUIF2AQUIF] > 0.0) ?
                MOBILEN_PROPORTION * river[i].no3sol.conc_bed :
                MOBILEN_PROPORTION * down->no3sol.conc_bed);

            river[i].nh4sol.flux[DOWN_AQUIF2AQUIF] =
                river[i].wf.rivflow[DOWN_AQUIF2AQUIF] * RHOH2O *
                ((river[i].wf.rivflow[DOWN_AQUIF2AQUIF] > 0.0) ?
                MOBILEN_PROPORTION * river[i].nh4sol.conc_bed :
                MOBILEN_PROPORTION * down->nh4sol.conc_bed);
        }
        else
        {
            river[i].no3sol.flux[DOWN_CHANL2CHANL] =
                river[i].wf.rivflow[DOWN_CHANL2CHANL] * RHOH2O *
                river[i].no3sol.conc_stream;

            river[i].nh4sol.flux[DOWN_CHANL2CHANL] =
                river[i].wf.rivflow[DOWN_CHANL2CHANL] * RHOH2O *
                river[i].nh4sol.conc_stream;

            river[i].no3sol.flux[DOWN_AQUIF2AQUIF] = 0.0;
            river[i].nh4sol.flux[DOWN_AQUIF2AQUIF] = 0.0;
        }

        /* Left and right banks */
        left = &elem[river[i].leftele - 1];
        right = &elem[river[i].rightele - 1];

        if (river[i].leftele > 0)
        {
            river[i].no3sol.flux[LEFT_SURF2CHANL] = 0.0;
            river[i].nh4sol.flux[LEFT_SURF2CHANL] = 0.0;

            river[i].no3sol.flux[LEFT_AQUIF2CHANL] =
                river[i].wf.rivflow[LEFT_AQUIF2CHANL] * RHOH2O *
                ((river[i].wf.rivflow[LEFT_AQUIF2CHANL] > 0.0) ?
                river[i].no3sol.conc_stream :
                MOBILEN_PROPORTION * left->no3sol.conc);

            river[i].nh4sol.flux[LEFT_AQUIF2CHANL] =
                river[i].wf.rivflow[LEFT_AQUIF2CHANL] * RHOH2O *
                ((river[i].wf.rivflow[LEFT_AQUIF2CHANL] > 0.0) ?
                river[i].nh4sol.conc_stream :
                MOBILEN_PROPORTION * left->nh4sol.conc);

            river[i].no3sol.flux[LEFT_AQUIF2AQUIF] =
                river[i].wf.rivflow[LEFT_AQUIF2AQUIF] * RHOH2O *
                ((river[i].wf.rivflow[LEFT_AQUIF2AQUIF] > 0.0) ?
                MOBILEN_PROPORTION * river[i].no3sol.conc_bed :
                MOBILEN_PROPORTION * left->no3sol.conc);

            river[i].nh4sol.flux[LEFT_AQUIF2AQUIF] =
                river[i].wf.rivflow[LEFT_AQUIF2AQUIF] * RHOH2O *
                ((river[i].wf.rivflow[LEFT_AQUIF2AQUIF] > 0.0) ?
                MOBILEN_PROPORTION * river[i].nh4sol.conc_bed :
                MOBILEN_PROPORTION * left->nh4sol.conc);

            for (j = 0; j < NUM_EDGE; j++)
            {
                if (left->nabr[j] == -(i + 1))
                {
                    left->no3sol.flux[j] =
                        -(river[i].no3sol.flux[LEFT_AQUIF2CHANL] +
                        river[i].no3sol.flux[LEFT_AQUIF2AQUIF]);

                    left->nh4sol.flux[j] =
                        -(river[i].nh4sol.flux[LEFT_AQUIF2CHANL] +
                        river[i].nh4sol.flux[LEFT_AQUIF2AQUIF]);
                    break;
                }
            }

        }

        if (river[i].rightele > 0)
        {
            river[i].no3sol.flux[RIGHT_SURF2CHANL] = 0.0;
            river[i].nh4sol.flux[RIGHT_SURF2CHANL] = 0.0;

            river[i].no3sol.flux[RIGHT_AQUIF2CHANL] =
                river[i].wf.rivflow[RIGHT_AQUIF2CHANL] * RHOH2O *
                ((river[i].wf.rivflow[RIGHT_AQUIF2CHANL] > 0.0) ?
                river[i].no3sol.conc_stream :
                MOBILEN_PROPORTION * right->no3sol.conc);

            river[i].nh4sol.flux[RIGHT_AQUIF2CHANL] =
                river[i].wf.rivflow[RIGHT_AQUIF2CHANL] * RHOH2O *
                ((river[i].wf.rivflow[RIGHT_AQUIF2CHANL] > 0.0) ?
                river[i].nh4sol.conc_stream :
                MOBILEN_PROPORTION * right->nh4sol.conc);

            river[i].no3sol.flux[RIGHT_AQUIF2AQUIF] =
                river[i].wf.rivflow[RIGHT_AQUIF2AQUIF] * RHOH2O *
                ((river[i].wf.rivflow[RIGHT_AQUIF2AQUIF] > 0.0) ?
                MOBILEN_PROPORTION * river[i].no3sol.conc_bed :
                MOBILEN_PROPORTION * right->no3sol.conc);

            river[i].nh4sol.flux[RIGHT_AQUIF2AQUIF] =
                river[i].wf.rivflow[RIGHT_AQUIF2AQUIF] * RHOH2O *
                ((river[i].wf.rivflow[RIGHT_AQUIF2AQUIF] > 0.0) ?
                MOBILEN_PROPORTION * river[i].nh4sol.conc_bed :
                MOBILEN_PROPORTION * right->nh4sol.conc);

            for (j = 0; j < NUM_EDGE; j++)
            {
                if (right->nabr[j] == -(i + 1))
                {
                    right->no3sol.flux[j] =
                        -(river[i].no3sol.flux[RIGHT_AQUIF2CHANL] +
                        river[i].no3sol.flux[RIGHT_AQUIF2AQUIF]);

                    right->nh4sol.flux[j] =
                        -(river[i].nh4sol.flux[RIGHT_AQUIF2CHANL] +
                        river[i].nh4sol.flux[RIGHT_AQUIF2AQUIF]);
                    break;
                }
            }
        }

        river[i].no3sol.flux[CHANL_LKG] =
            river[i].wf.rivflow[CHANL_LKG] * RHOH2O *
            ((river[i].wf.rivflow[CHANL_LKG] > 0.0) ?
            river[i].no3sol.conc_stream :
            MOBILEN_PROPORTION * river[i].no3sol.conc_bed);

        river[i].nh4sol.flux[CHANL_LKG] =
            river[i].wf.rivflow[CHANL_LKG] * RHOH2O *
            ((river[i].wf.rivflow[CHANL_LKG] > 0.0) ?
            river[i].nh4sol.conc_stream :
            MOBILEN_PROPORTION * river[i].nh4sol.conc_bed);
    }

    /*
     * Accumulate to get in-flow for down segments
     */
    for (i = 0; i < nriver; i++)
    {
        river_struct   *down;

        if (river[i].down > 0)
        {
            down = &river[river[i].down - 1];

            down->no3sol.flux[UP_CHANL2CHANL] -=
                river[i].no3sol.flux[DOWN_CHANL2CHANL];

            down->no3sol.flux[UP_AQUIF2AQUIF] -=
                river[i].no3sol.flux[DOWN_AQUIF2AQUIF];

            down->nh4sol.flux[UP_CHANL2CHANL] -=
                river[i].nh4sol.flux[DOWN_CHANL2CHANL];

            down->nh4sol.flux[UP_AQUIF2AQUIF] -=
                river[i].nh4sol.flux[DOWN_AQUIF2AQUIF];
        }
    }
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
