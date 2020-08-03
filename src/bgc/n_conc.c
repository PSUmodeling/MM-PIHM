#include "pihm.h"

#if !defined(_LUMPEDBGC_) && !defined(_LEACHING_)
void SoluteConc(elem_struct elem[], river_struct river[])
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
        double          storage;

        /* Element surface */
        elem[i].solute[0].conc_surf = 0.0;

        /* Element subsurface */
        storage = (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity +
            elem[i].soil.depth * elem[i].soil.smcmin;
        elem[i].solute[0].conc = (storage > 0.0) ?
            MOBILEN_PROPORTION * elem[i].ns.sminn / storage : 0.0;
        elem[i].solute[0].conc = MAX(elem[i].solute[0].conc, 0.0);
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             j;
        double          storage;

        storage = river[i].ws.stage;
        river[i].solute[0].conc = (storage > 0.0) ?
            river[i].ns.streamn / storage : 0.0;
        river[i].solute[0].conc = MAX(river[i].solute[0].conc, 0.0);
    }
}
#endif

#if defined(_LUMPEDBGC_)
void NLeachingLumped(elem_struct *elem, river_struct *river)
{
    int             i;
    double          strg = 0.0;      /* Total water storage (m3 m-2) */
    double          runoff = 0.0;    /* Total runoff (kg m-2 s-1) */

    for (i = 0; i < nelem; i++)
    {
        int             j;

        strg += ((elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity +
            elem[i].soil.depth * elem[i].soil.smcmin) * elem[i].topo.area;
    }

    for (i = 0; i < nriver; i++)
    {
        if (river[i].down < 0)
        {
           runoff += river[i].wf.rivflow[DOWN_CHANL2CHANL] * 1000.0;
        }
    }

    strg /= elem[LUMPEDBGC].topo.area;
    runoff /= elem[LUMPEDBGC].topo.area;

    elem[LUMPEDBGC].nf.sminn_leached = (runoff > 0.0) ?
        runoff * MOBILEN_PROPORTION * elem[LUMPEDBGC].ns.sminn / strg / 1000.0 :
        0.0;
}
#endif

#if defined(_LEACHING_)
void NLeaching(elem_struct *elem)
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
        double          strg;

        /* Initialize N fluxes */
        for (j = 0; j < NUM_EDGE; j++)
        {
            elem[i].nsol.subflux[j] = 0.0;
        }

        /* Element subsurface */
        strg = (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity +
            elem[i].soil.depth * elem[i].soil.smcmin;
        elem[i].nsol.conc_subsurf = (strg > 0.0) ?
            elem[i].ns.sminn / strg / 1000.0 : 0.0;
        elem[i].nsol.conc_subsurf = (elem[i].nsol.conc_subsurf > 0.0) ?
            elem[i].nsol.conc_subsurf : 0.0;
    }

    /*
     * Calculate solute fluxes
     */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        /* Element to element */
        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].nabr[j] > 0)
            {
                elem[i].nsol.subflux[j] = elem[i].wf.subsurf[j] * 1000.0 *
                    ((elem[i].wf.subsurf[j] > 0.0) ?
                    MOBILEN_PROPORTION * elem[i].nsol.conc_subsurf : 0.0);
            }
            else
            {
                elem[i].nsol.subflux[j] = 0.0;
            }
        }
    }
}
#endif
