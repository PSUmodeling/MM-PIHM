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

const double        MINSTRG = 1.0E-4;
const int           UP = 0;
const int           DOWN = 1;
const int           LEFT = 2;
const int           RIGHT = 3;

void NTransport (elem_struct *elem, river_struct *riv, double dt)
{
    int             i;
    int             nsteps = 1;

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

    /*
     * Copy information to solute structures
     */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int         j, k;

        elem[i].nsol.prev_strg = elem[i].daily.prev_surf +
            (elem[i].daily.prev_unsat + elem[i].daily.prev_gw) *
            elem[i].soil.porosity;
        elem[i].nsol.strg = elem[i].daily.surf +
            (elem[i].daily.unsat + elem[i].daily.gw) * elem[i].soil.porosity;
        //elem[i].nsol.prev_strg = elem[i].daily.prev_surf;
        //elem[i].nsol.strg = elem[i].daily.surf;

        //for (k = 0; k < elem[i].ps.nsoil; k++)
        //{
        //    elem[i].nsol.prev_strg +=
        //        elem[i].daily.prev_smc[k] * elem[i].ps.sldpth[k];
        //    elem[i].nsol.strg +=
        //        elem[i].daily.smc[k] * elem[i].ps.sldpth[k];
        //}

        elem[i].nsol.prev_strg = (elem[i].nsol.prev_strg > MINSTRG) ?
            elem[i].nsol.prev_strg : MINSTRG;
        elem[i].nsol.prev_strg *= 1000.0;
        elem[i].nsol.strg = (elem[i].nsol.strg > MINSTRG) ?
            elem[i].nsol.strg : MINSTRG;
        elem[i].nsol.strg *= 1000.0;

        for (j = 0; j < NUM_EDGE; j++)
        {
            elem[i].nsol.wflux[j] = (elem[i].daily.avg_ovlflow[j] +
                elem[i].daily.avg_subsurf[j]) *
                1000.0 / elem[i].topo.area;
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        riv[i].nsol.prev_strg = riv[i].daily.prev_stage +
            riv[i].daily.prev_gw * riv[i].matl.porosity;
        riv[i].nsol.strg = riv[i].daily.stage +
            riv[i].daily.gw * riv[i].matl.porosity;

        riv[i].nsol.prev_strg = (riv[i].nsol.prev_strg > MINSTRG) ?
            riv[i].nsol.prev_strg : MINSTRG;
        riv[i].nsol.prev_strg *= 1000.0;
        riv[i].nsol.strg = (riv[i].nsol.strg > MINSTRG) ?
            riv[i].nsol.strg : MINSTRG;
        riv[i].nsol.strg *= 1000.0;

        riv[i].nsol.wflux[UP] = (riv[i].daily.avg_rivflow[UP_CHANL2CHANL] +
            riv[i].daily.avg_rivflow[UP_AQUIF2AQUIF]) *
            1000.0 / riv[i].topo.area;

        riv[i].nsol.wflux[DOWN] =
            (riv[i].daily.avg_rivflow[DOWN_CHANL2CHANL] +
            riv[i].daily.avg_rivflow[DOWN_AQUIF2AQUIF]) *
            1000.0 / riv[i].topo.area;

        riv[i].nsol.wflux[LEFT] =
            (riv[i].daily.avg_rivflow[LEFT_SURF2CHANL] +
            riv[i].daily.avg_rivflow[LEFT_AQUIF2CHANL] +
            riv[i].daily.avg_rivflow[LEFT_AQUIF2AQUIF]) *
            1000.0 / riv[i].topo.area;

        riv[i].nsol.wflux[RIGHT] =
            (riv[i].daily.avg_rivflow[RIGHT_SURF2CHANL] +
            riv[i].daily.avg_rivflow[RIGHT_AQUIF2CHANL] +
            riv[i].daily.avg_rivflow[RIGHT_AQUIF2AQUIF]) *
            1000.0 / riv[i].topo.area;
    }

    for (i = 0; i < nelem; i++)
    {
        int         j;
        double      totsnk;
        int         ratio;

        totsnk = 0.0;

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].nsol.wflux[j] > 0.0)
            {
                totsnk += elem[i].nsol.wflux[j] * dt;
            }
        }

        if (totsnk > elem[i].nsol.prev_strg)
        {
            ratio = (int)ceil (totsnk / elem[i].nsol.prev_strg);
            nsteps = (ratio > nsteps) ? ratio : nsteps;
        }
    }

    for (i = 0; i < nriver; i++)
    {
        int         j;
        double      totsnk;
        int         ratio;

        totsnk = 0.0;

        for (j = 0; j < 4; j++)
        {
            if (riv[i].nsol.wflux[j] > 0.0)
            {
                totsnk += riv[i].nsol.wflux[j] * dt;
            }
        }

        if (totsnk > riv[i].nsol.prev_strg)
        {
            ratio = (int)ceil (totsnk / riv[i].nsol.prev_strg);
            nsteps = (ratio > nsteps) ? ratio : nsteps;
        }
    }

    AdptStepTrnsp (elem, riv, dt, nsteps);
}

void AdptStepTrnsp (elem_struct *elem, river_struct *riv, double dt,
    int nsteps)
{
    int             i;
    int             step;
    double          frac;
    double          adpt_dt;

    adpt_dt = dt / (double)nsteps;

    /*
     * Initialize nf.sminn_leached fluxes
     */
    for (i = 0; i < nelem; i++)
    {
        elem[i].nf.sminn_leached = 0.0;
    }

    for (i = 0; i < nriver; i++)
    {
        riv[i].nf.sminn_leached = 0.0;
    }

    /*
     * Local time step loops
     */
    for (step = 0; step < nsteps; step++)
    {
        frac = (double)(step + 1) / (double)nsteps;

        /*
         * Calculate soluble N concentrations
         */
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < nelem; i++)
        {
            double  local_strg;

            local_strg = elem[i].nsol.prev_strg +
                (elem[i].nsol.strg - elem[i].nsol.prev_strg) * frac;

            elem[i].nsol.conc =
                MOBILEN_PROPORTION * elem[i].ns.sminn / local_strg;

            elem[i].nsol.conc = (elem[i].nsol.conc > 0.0) ?
                elem[i].nsol.conc : 0.0;
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < nriver; i++)
        {
            double      local_strg;

            local_strg = riv[i].nsol.prev_strg +
                (riv[i].nsol.strg - riv[i].nsol.prev_strg) * frac;

            riv[i].nsol.conc =
                MOBILEN_PROPORTION * riv[i].ns.sminn / local_strg;

            riv[i].nsol.conc = (riv[i].nsol.conc > 0.0) ?
                riv[i].nsol.conc : 0.0;
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < nelem; i++)
        {
            int         j;

            elem[i].nsol.trnsp_flux = 0.0;

            for (j = 0; j < NUM_EDGE; j++)
            {
                if (elem[i].nabr[j] > 0)
                {
                    elem[i].nsol.trnsp_flux += (elem[i].nsol.wflux[j] > 0.0) ?
                        elem[i].nsol.conc * elem[i].nsol.wflux[j] * adpt_dt :
                        elem[elem[i].nabr[j] - 1].nsol.conc *
                        elem[i].nsol.wflux[j] * adpt_dt;
                }
                else if (elem[i].nabr[j] < 0)
                {
                    elem[i].nsol.trnsp_flux += (elem[i].nsol.wflux[j] > 0.0) ?
                        elem[i].nsol.conc * elem[i].nsol.wflux[j] * adpt_dt :
                        riv[-elem[i].nabr[j] - 1].nsol.conc *
                        elem[i].nsol.wflux[j] * adpt_dt;
                }
            }

            elem[i].nsol.trnsp_flux =
                (elem[i].nsol.trnsp_flux < elem[i].ns.sminn) ?
                elem[i].nsol.trnsp_flux : elem[i].ns.sminn;

            elem[i].nf.sminn_leached += elem[i].nsol.trnsp_flux;
            elem[i].ns.nleached_snk += elem[i].nsol.trnsp_flux;
            elem[i].ns.sminn -= elem[i].nsol.trnsp_flux;
#ifdef _DEBUG_
            if (fabs(elem[i].nsol.trnsp_flux) > 1.0e3)
            {
                printf ("Error\n");
            }
#endif
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < nriver; i++)
        {
            riv[i].nsol.trnsp_flux = 0.0;

            /* Upstream */
            if (riv[i].up > 0)
            {
                riv[i].nsol.trnsp_flux += (riv[i].nsol.wflux[UP] > 0.0) ?
                    riv[i].nsol.conc * riv[i].nsol.wflux[UP] * adpt_dt :
                    riv[riv[i].up - 1].nsol.conc *
                    riv[i].nsol.wflux[UP] * adpt_dt;
            }

            /* Downstream */
            if (riv[i].down > 0)
            {
                riv[i].nsol.trnsp_flux += (riv[i].nsol.wflux[DOWN] > 0.0) ?
                    riv[i].nsol.conc * riv[i].nsol.wflux[DOWN] * adpt_dt :
                    riv[riv[i].down - 1].nsol.conc *
                    riv[i].nsol.wflux[DOWN] * adpt_dt;
            }
            else
            {
                riv[i].nsol.trnsp_flux +=
                    riv[i].nsol.conc * riv[i].nsol.wflux[DOWN] * adpt_dt;
            }

            /* Left bank */
            riv[i].nsol.trnsp_flux += (riv[i].nsol.wflux[LEFT] > 0.0) ?
                riv[i].nsol.conc * riv[i].nsol.wflux[LEFT] * adpt_dt :
                elem[riv[i].leftele - 1].nsol.conc *
                riv[i].nsol.wflux[LEFT] * adpt_dt;

            /* Right bank */
            riv[i].nsol.trnsp_flux += (riv[i].nsol.wflux[RIGHT] > 0.0) ?
                riv[i].nsol.conc * riv[i].nsol.wflux[RIGHT] * adpt_dt :
                elem[riv[i].rightele - 1].nsol.conc *
                riv[i].nsol.wflux[RIGHT] * adpt_dt;

            riv[i].nsol.trnsp_flux =
                (riv[i].nsol.trnsp_flux < riv[i].ns.sminn) ?
                riv[i].nsol.trnsp_flux : riv[i].ns.sminn;

            riv[i].nf.sminn_leached += riv[i].nsol.trnsp_flux;
            riv[i].ns.sminn -= riv[i].nsol.trnsp_flux;
        }
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
