#include "pihm.h"

void InitMesh(elem_struct *elem, const meshtbl_struct *meshtbl)
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        elem[i].ind = i + 1;

        for (j = 0; j < NUM_EDGE; j++)
        {
            elem[i].node[j] = meshtbl->node[i][j];
            elem[i].nabr[j] = meshtbl->nabr[i][j];
            elem[i].nabr_river[j] = 0;      /* initialize to 0 */
        }
    }
}
