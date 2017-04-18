/*
 * nleaching.c
 * daily nitrogen leaching to groundwater
 *
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 * Biome-BGC version 4.2 (final release)
 * See copyright.txt for Copyright information
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 */

#include "pihm.h"

void NTransport (elem_struct *elem, int numele, river_struct *riv, int numriv)
{
    int             i;

    /* N leaching flux is calculated after all the other nfluxes are
     * reconciled to avoid the possibility of removing more N than is there.
     * This follows the implicit logic of precedence for soil mineral N
     * resources:
     * 1) microbial processes and plant uptake (competing)
     * 2) leaching
     *
     * leaching happens when there is outflow, as a function of the presumed
     * proportion of the soil mineral N pool which is soluble (nitrates), the
     * soil water content, and the outflow */

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < numele; i++)
    {
        int         k;
        double      totwater;

        totwater = elem[i].daily.avg_surf;

        for (k = 0; k < elem[i].ps.nsoil; k++)
        {
            totwater += elem[i].daily.avg_smc[k] * elem[i].ps.sldpth[k];
        }

        totwater *= 1000.0;

        elem[i].nsol.conc = MOBILEN_PROPORTION * elem[i].ns.sminn / totwater;
        elem[i].nsol.conc = (elem[i].nsol.conc > 0.0) ?
            elem[i].nsol.conc : 0.0;
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < numriv; i++)
    {
        double      totwater;

        totwater = (riv[i].daily.avg_stage +
            riv[i].daily.avg_gw * riv[i].matl.porosity) * 1000.0;
        riv[i].nsol.conc = MOBILEN_PROPORTION * riv[i].ns.sminn / totwater;
        riv[i].nsol.conc = (riv[i].nsol.conc > 0.0) ?
            riv[i].nsol.conc : 0.0;
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < numele; i++)
    {
        int         j;
        double      wflux;

        elem[i].nf.sminn_leached = 0.0;

        for (j = 0; j < NUM_EDGE; j++)
        {
            wflux = (elem[i].daily.avg_subsurf[j] +
                elem[i].daily.avg_ovlflow[j]) *
                24.0 * 3600.0 * 1000.0 / elem[i].topo.area;

            if (wflux > 0.0)
            {
                elem[i].nf.sminn_leached +=
                    wflux * elem[i].nsol.conc;
            }
            else
            {
                if (elem[i].nabr[j] > 0)
                {
                    elem[i].nf.sminn_leached +=
                        wflux * elem[elem[i].nabr[j] - 1].nsol.conc;
                }
                else if (elem[i].nabr[j] < 0)
                {
                    elem[i].nf.sminn_leached +=
                        wflux * riv[-elem[i].nabr[j] - 1].nsol.conc;
                }
            }
        }

        elem[i].nf.sminn_leached =
            (elem[i].nf.sminn_leached < elem[i].ns.sminn) ?
            elem[i].nf.sminn_leached : elem[i].ns.sminn;

        elem[i].ns.nleached_snk += elem[i].nf.sminn_leached;
        elem[i].ns.sminn -= elem[i].nf.sminn_leached;
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < numriv; i++)
    {
        double      wflux;

        riv[i].nf.sminn_leached = 0.0;

        /* Upstream */
        wflux = (riv[i].daily.avg_rivflow[UP_CHANL2CHANL] +
            riv[i].daily.avg_rivflow[UP_AQUIF2AQUIF]) *
            1000.0 * 24.0 * 3600.0 / riv[i].topo.area;

        if (wflux > 0.0)
        {
            riv[i].nf.sminn_leached +=
                wflux * riv[i].nsol.conc;
        }
        else
        {
            if (riv[i].up > 0)
            {
                riv[i].nf.sminn_leached +=
                    wflux * riv[riv[i].up - 1].nsol.conc;
            }
        }

        /* Downstream */
        wflux = (riv[i].daily.avg_rivflow[DOWN_CHANL2CHANL] +
            riv[i].daily.avg_rivflow[DOWN_AQUIF2AQUIF]) *
            1000.0 * 24.0 * 3600.0 / riv[i].topo.area;

        if (wflux > 0.0)
        {
            riv[i].nf.sminn_leached +=
                wflux * riv[i].nsol.conc;
        }
        else
        {
            if (riv[i].down > 0)
            {
                riv[i].nf.sminn_leached +=
                    wflux * riv[riv[i].down - 1].nsol.conc;
            }
        }

        /* Left bank */
        wflux = (riv[i].daily.avg_rivflow[LEFT_SURF2CHANL] +
            riv[i].daily.avg_rivflow[LEFT_AQUIF2CHANL] +
            riv[i].daily.avg_rivflow[LEFT_AQUIF2AQUIF]) *
            1000.0 * 24.0 * 3600.0 / riv[i].topo.area;

        if (wflux > 0.0)
        {
            riv[i].nf.sminn_leached +=
                wflux * riv[i].nsol.conc;
        }
        else
        {
            riv[i].nf.sminn_leached +=
                wflux * elem[riv[i].leftele - 1].nsol.conc;
        }

        /* Right bank */
        wflux = (riv[i].daily.avg_rivflow[RIGHT_SURF2CHANL] +
            riv[i].daily.avg_rivflow[RIGHT_AQUIF2CHANL] +
            riv[i].daily.avg_rivflow[RIGHT_AQUIF2AQUIF]) *
            1000.0 * 24.0 * 3600.0 / riv[i].topo.area;

        if (wflux > 0.0)
        {
            riv[i].nf.sminn_leached +=
                wflux * riv[i].nsol.conc;
        }
        else
        {
            riv[i].nf.sminn_leached +=
                wflux * elem[riv[i].rightele - 1].nsol.conc;
        }

        if (riv[i].nf.sminn_leached > riv[i].ns.sminn)
        {
            riv[i].nf.sminn_leached = riv[i].ns.sminn;
        }

        riv[i].ns.sminn -= riv[i].nf.sminn_leached;
    }
}

//void nleaching_nt (bgc_grid *grid, int numele, bgc_river *riv, int numriv)
//{
//    double          nconc;
//    int             i;
//    double          total_soilwater = 0.0;
//    double          total_sminn = 0.0;
//    double          nleached;
//    double          outflow;
//
//    for (i = 0; i < numele; i++)
//    {
//        total_soilwater += grid[i].ws.soilw;
//        total_sminn += grid[i].ns.sminn;
//    }
//
//    nconc = total_sminn / total_soilwater;
//
//    outflow = riv[numriv - 1].metv.latflux[1];
//
//    nleached = MOBILEN_PROPORTION * nconc * outflow;
//
//    for (i = 0; i < numele; i++)
//    {
//        grid[i].nf.sminn_leached = nleached / numele;
//
//        grid[i].nf.sminn_leached = (grid[i].nf.sminn_leached < grid[i].ns.sminn) ? grid[i].nf.sminn_leached : grid[i].ns.sminn;
//
//        grid[i].ns.nleached_snk += grid[i].nf.sminn_leached;
//        grid[i].ns.sminn -= grid[i].nf.sminn_leached;
//    }
//}
