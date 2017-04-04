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
    double         *nconc;
    int             i;
    const int       UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3;

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

    nconc = (double *)malloc ((numele + numriv) * sizeof (double));

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

        nconc[i] = elem[i].ns.sminn / totwater;
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < numriv; i++)
    {
        double      totwater;

        totwater = (riv[i].daily.avg_stage +
            riv[i].daily.avg_gw * riv[i].matl.porosity) * 1000.0;
        nconc[i + numele] = riv[i].ns.sminn / totwater;
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < numele; i++)
    {
        int         j;
        double      nabr_nconc[4];
        double      latflux;

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].nabr[j] > 0)
            {
                nabr_nconc[j] = nconc[elem[i].nabr[j] - 1];
            }
            else if (elem[i].nabr[j] < 0)
            {
                nabr_nconc[j] = nconc[numele - elem[i].nabr[j] - 1];
            }
            else
            {
                nabr_nconc[j] = 0.0;
            }

            nabr_nconc[j] = (nabr_nconc[j] > 0.0) ? nabr_nconc[j] : 0.0;
        }

        elem[i].nf.sminn_leached = 0.0;

        for (j = 0; j < NUM_EDGE; j++)
        {
            latflux = (elem[i].daily.avg_subsurf[j] +
                elem[i].daily.avg_ovlflow[j]) *
                24.0 * 3600.0 * 1000.0 / elem[i].topo.area;

            if (latflux > 0.0)
            {
                elem[i].nf.sminn_leached +=
                    MOBILEN_PROPORTION * nconc[i] * latflux;
            }
            else
            {
                elem[i].nf.sminn_leached +=
                    MOBILEN_PROPORTION * nabr_nconc[j] * latflux;
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
        double      nleached = 0.0;
        double      nabr_nconc[4];
        double      latflux;

        /* Upstream */
        if (riv[i].up > 0)
        {
            nabr_nconc[UP] = nconc[numele + riv[i].up - 1];
        }
        else
        {
            nabr_nconc[UP] = 0.0;
        }
        nabr_nconc[UP] = (nabr_nconc[UP] > 0.0) ? nabr_nconc[UP] : 0.0;

        latflux = (riv[i].daily.avg_rivflow[UP_CHANL2CHANL] +
            riv[i].daily.avg_rivflow[UP_AQUIF2AQUIF]) *
            1000.0 * 24.0 * 3600.0 / riv[i].topo.area;

        if (latflux > 0.0)
        {
            nleached += MOBILEN_PROPORTION * nconc[numele + i] * latflux;
        }
        else
        {
            nleached += MOBILEN_PROPORTION * nabr_nconc[UP] * latflux;
        }

        /* Downstream */
        if (riv[i].down > 0)
        {
            nabr_nconc[DOWN] = nconc[numele + riv[i].down - 1];
        }
        else
        {
            nabr_nconc[DOWN] = 0.0;
        }
        nabr_nconc[DOWN] = (nabr_nconc[DOWN] > 0.0) ? nabr_nconc[DOWN] : 0.0;

        latflux = (riv[i].daily.avg_rivflow[DOWN_CHANL2CHANL] +
            riv[i].daily.avg_rivflow[DOWN_AQUIF2AQUIF]) *
            1000.0 * 24.0 * 3600.0 / riv[i].topo.area;

        if (latflux > 0.0)
        {
            nleached += MOBILEN_PROPORTION * nconc[numele + i] * latflux;
        }
        else
        {
            nleached += MOBILEN_PROPORTION * nabr_nconc[DOWN] * latflux;
        }

        /* Left bank */
        if (riv[i].leftele > 0)
        {
            nabr_nconc[LEFT] = nconc[riv[i].leftele - 1];
        }
        else
        {
            nabr_nconc[LEFT] = 0.0;
        }
        nabr_nconc[LEFT] = (nabr_nconc[LEFT] > 0.0) ? nabr_nconc[LEFT] : 0.0;

        latflux = (riv[i].daily.avg_rivflow[LEFT_SURF2CHANL] +
            riv[i].daily.avg_rivflow[LEFT_AQUIF2CHANL] +
            riv[i].daily.avg_rivflow[LEFT_AQUIF2AQUIF]) *
            1000.0 * 24.0 * 3600.0 / riv[i].topo.area;

        if (latflux > 0.0)
        {
            nleached += MOBILEN_PROPORTION * nconc[numele + i] * latflux;
        }
        else
        {
            nleached += MOBILEN_PROPORTION * nabr_nconc[LEFT] * latflux;
        }

        if (riv[i].rightele > 0)
        {
            nabr_nconc[RIGHT] = nconc[riv[i].rightele - 1];
        }
        else
        {
            nabr_nconc[RIGHT] = 0.0;
        }
        nabr_nconc[RIGHT] =
            (nabr_nconc[RIGHT] > 0.0) ? nabr_nconc[RIGHT] : 0.0;

        latflux = (riv[i].daily.avg_rivflow[RIGHT_SURF2CHANL] +
            riv[i].daily.avg_rivflow[RIGHT_AQUIF2CHANL] +
            riv[i].daily.avg_rivflow[RIGHT_AQUIF2AQUIF]) *
            1000.0 * 24.0 * 3600.0 / riv[i].topo.area;

        if (latflux > 0.0)
        {
            nleached += MOBILEN_PROPORTION * nconc[numele + i] * latflux;
        }
        else
        {
            nleached += MOBILEN_PROPORTION * nabr_nconc[RIGHT] * latflux;
        }

        if (nleached > riv[i].ns.sminn)
        {
            nleached = riv[i].ns.sminn;
        }

        riv[i].nf.sminn_leached = nleached;
        riv[i].ns.sminn -= nleached;
    }

    free (nconc);
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
