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

void nleaching (nstate_struct * ns, nflux_struct * nf, wstate_struct * ws, const double *nabr_nconc, const metvar_struct *metv)
{
    double          soilwater_nconc;
    int             k;

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
    //if (wf->soilw_outflow)
    //{
    //    soilwater_nconc = MOBILEN_PROPORTION * ns->sminn / ws->soilw;
    //    nf->sminn_leached = soilwater_nconc * wf->soilw_outflow;
    //    /* update state variables */
    //    ns->nleached_snk += nf->sminn_leached;
    //    ns->sminn -= nf->sminn_leached;
    //}

    soilwater_nconc = MOBILEN_PROPORTION * ns->sminn / ws->soilw;

    nf->sminn_leached = 0.0;
    for (k = 0; k < 3; k++)
    {
        if (metv->subflux[k] > 0.0)
        {
            nf->sminn_leached += soilwater_nconc * metv->subflux[k];
        }
        else
            nf->sminn_leached += nabr_nconc[k] * metv->subflux[k];

        ns->nleached_snk += nf->sminn_leached;
        ns->sminn -= nf->sminn_leached;
    }
}
