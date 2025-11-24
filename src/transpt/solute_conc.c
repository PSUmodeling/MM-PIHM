#include "pihm.h"

void SoluteConc(elem_struct elem[], river_struct river[])
{
    int             kelem, kriver, ksolute;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (kelem = 0; kelem < nelem; kelem++)
    {
        double          storage;

        storage = (elem[kelem].ws.unsat + elem[kelem].ws.gw) * elem[kelem].soil.porosity + elem[kelem].soil.depth * elem[kelem].soil.smcmin;

        for (ksolute = 0; ksolute < nsolute; ksolute++)
        {
            elem[kelem].solute[ksolute].conc_surf = 0.0;
            elem[kelem].solute[ksolute].conc = (storage > 0.0) ? elem[kelem].solute[ksolute].amount / storage : 0.0;
            elem[kelem].solute[ksolute].conc = MAX(elem[kelem].solute[ksolute].conc, 0.0);
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (kriver = 0; kriver < nriver; kriver++)
    {
        double          storage;

        storage = river[kriver].ws.stage;

        for (ksolute = 0; ksolute < nsolute; ksolute++)
        {
            river[kriver].solute[ksolute].conc = (storage > 0.0) ? river[kriver].solute[ksolute].amount / storage : 0.0;
            river[kriver].solute[ksolute].conc = MAX(river[kriver].solute[ksolute].conc, 0.0);
        }
    }
}
