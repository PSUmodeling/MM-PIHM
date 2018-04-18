#include "pihm.h"

void InitMesh(elem_struct *elem, const meshtbl_struct *meshtbl)
{
    int             i, j;

    for (i = 0; i < nelem; i++)
    {
        elem[i].ind = i + 1;

        for (j = 0; j < NUM_EDGE; j++)
        {
            elem[i].node[j] = meshtbl->node[i][j];
            elem[i].nabr[j] = meshtbl->nabr[i][j];
        }
    }
}
