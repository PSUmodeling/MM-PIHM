#include "pihm.h"

void SoluteConc(elem_struct elem[], river_struct river[])
{
    int             i;

    // Calculate solute N concentrations
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        double          storage;

        // Element surface
        elem[i].solute[0].conc_surf = 0.0;

        // Element subsurface
        storage = (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity + elem[i].soil.depth * elem[i].soil.smcmin;
        elem[i].solute[0].conc = (storage > 0.0) ? MOBILEN_PROPORTION * elem[i].ns.sminn / storage : 0.0;
        elem[i].solute[0].conc = MAX(elem[i].solute[0].conc, 0.0);
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        double          storage;

        storage = river[i].ws.stage;
        river[i].solute[0].conc = (storage > 0.0) ? river[i].ns.streamn / storage : 0.0;
        river[i].solute[0].conc = MAX(river[i].solute[0].conc, 0.0);
    }
}
