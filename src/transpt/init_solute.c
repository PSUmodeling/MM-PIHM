#include "pihm.h"

void InitSolute(elem_struct elem[])
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j, k;

        for (k = 0; k < nsolute; k++)
        {
            elem[i].solute[k].conc_surf      = 0.0;
            elem[i].solute[k].conc           = 0.0;

            elem[i].solute[k].infil          = 0.0;
            for (j = 0; j < NUM_EDGE; j++)
            {
                elem[i].solute[k].subflux[j] = 0.0;
            }
#if defined(_CYCLES_)
            for (j = 0; j < MAXLYR; j++)
            {
                elem[i].solute[k].snksrc[j]  = 0.0;
            }
#else
            elem[i].solute[k].snksrc         = 0.0;
#endif

#if defined(_FBR_)
            elem[i].solute[k].conc_geol      = 0.0;
            elem[i].solute[k].fbr_infil      = 0.0;
            for (j = 0; j < NUM_EDGE; j++)
            {
                elem[i].solute[k].fbrflow[j] = 0.0;
            }
            elem[i].solute[k].snksrc_geol    = 0.0;
# if defined(_LUMPED_)
            elem[i].solute[k].fbr_discharge  = 0.0;
# endif
#endif
        }
    }
}
