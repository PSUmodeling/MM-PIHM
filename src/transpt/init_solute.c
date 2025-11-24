#include "pihm.h"

#if !defined(_BGC_) && !defined(_CYCLES_)
void InitSolute(elem_struct elem[], river_struct river[], solutetbl_struct *solutetbl)
#else
void InitSolute(elem_struct elem[], river_struct river[])
#endif
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j, k;

        elem[i].solute = (solute_struct *)malloc(nsolute * sizeof(solute_struct));

        for (k = 0; k < nsolute; k++)
        {
#if !defined(_BGC_) && !defined(_CYCLES_)
            elem[i].solute[k].amount = solutetbl->amount[k][i];
#endif
            elem[i].solute[k].conc_surf = 0.0;
            elem[i].solute[k].conc = 0.0;

            elem[i].solute[k].infil = 0.0;
            for (j = 0; j < NUM_EDGE; j++)
            {
                elem[i].solute[k].subflux[j] = 0.0;
            }
#if defined(_CYCLES_)
            for (j = 0; j < MAXLYR; j++)
            {
                elem[i].solute[k].snksrc[j] = 0.0;
            }
#else
            elem[i].solute[k].snksrc = 0.0;
#endif
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             j, k;

        river[i].solute = (river_solute_struct *)malloc(nsolute * sizeof(river_solute_struct));

        for (k = 0; k < nsolute; k++)
        {
#if !defined(_BGC_) && !defined(_CYCLES_)
            river[i].solute[k].amount = 0.0;
#endif
            river[i].solute[k].conc = 0.0;

            for (j = 0; j < NUM_RIVFLX; j++)
            {
                river[i].solute[k].flux[j] = 0.0;
            }
        }
    }
}
