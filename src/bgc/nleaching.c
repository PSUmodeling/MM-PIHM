/* 
 * nleaching.c
 * daily nitrogen leaching to groundwater
 * 
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 * Biome-BGC version 4.2 (final release)
 * See copyright.txt for Copyright information
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 */

#include "bgc.h"

void nleaching (bgc_grid *grid, int numele, bgc_river *riv, int numriv)
{
    double         *soilwater_nconc;
    int             i, j;
    double          nleached;
    double          nabr_nconc[4];
    siteconst_struct *sitec;
    metvar_struct  *metv;

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

    soilwater_nconc = (double *)malloc ((numele + numriv) * sizeof (double));

    for (i = 0; i < numele; i++)
    {
        soilwater_nconc[i] = grid[i].ns.sminn / grid[i].ws.soilw;
        //printf ("%lf\t", soilwater_nconc[i]);
    }

    for (i = 0; i < numriv; i++)
    {
        soilwater_nconc[i + numele] = riv[i].sminn / riv[i].soilw;
        //printf ("%lf\t", soilwater_nconc[i + numele]);
    }
    //printf ("\n");

    for (i = 0; i < numele; i++)
    {
        sitec = &grid[i].sitec;
        metv = &grid[i].metv;

        for (j = 0; j < 3; j++)
        {
            if (sitec->nabr[j] > 0)
            {
                nabr_nconc[j] = soilwater_nconc[sitec->nabr[j] - 1];
            }
            else if (sitec->nabr[j] < 0)
            {
                nabr_nconc[j] = soilwater_nconc[numele - sitec->nabr[j] - 1];
            }
            else
            {
                nabr_nconc[j] = 0.0;
            }

            nabr_nconc[j] = (nabr_nconc[j] > 0.0) ? nabr_nconc[j] : 0.0;
        }

        grid[i].nf.sminn_leached = 0.0;

        for (j = 0; j < 3; j++)
        {
            if (metv->latflux[j] > 0.0)
            {
                grid[i].nf.sminn_leached += MOBILEN_PROPORTION * soilwater_nconc[i] * metv->latflux[j];
            }
            else
            {
                grid[i].nf.sminn_leached += MOBILEN_PROPORTION * nabr_nconc[j] * metv->latflux[j];
            }
        }

        grid[i].ns.nleached_snk += grid[i].nf.sminn_leached;
        grid[i].ns.sminn -= grid[i].nf.sminn_leached;
    }

    for (i = 0; i < numriv; i++)
    {
        metv = &riv[i].metv;

        for (j = 0; j < 4; j++)
        {
            if (riv[i].nabr[j] > 0)
            {
                nabr_nconc[j] = soilwater_nconc[riv[i].nabr[j] - 1];
            }
            else if (riv[i].nabr[j] < 0)
            {
                nabr_nconc[j] = soilwater_nconc[numele - riv[i].nabr[j] - 1];
            }
            else
            {
                nabr_nconc[j] = 0.0;
            }

            nabr_nconc[j] = (nabr_nconc[j] > 0.0) ? nabr_nconc[j] : 0.0;
        }

        nleached = 0.0;

        for (j = 0; j < 4; j++)
        {
            if (metv->latflux[j] > 0.0)
            {
                nleached += MOBILEN_PROPORTION * soilwater_nconc[numele + i] * metv->latflux[j];
            }
            else
            {
                nleached += MOBILEN_PROPORTION * nabr_nconc[j] * metv->latflux[j];
            }
        }

        nleached = (nleached < riv[i].sminn) ? nleached : riv[i].sminn;

        riv[i].sminn -= nleached;
    }
    fflush (stdout);

    free (soilwater_nconc);
}
