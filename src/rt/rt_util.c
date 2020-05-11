#include "pihm.h"

int FindChem(const char chemn[MAXSTRING], const chemtbl_struct  chemtbl[],
    int nsps)
{
    int             i;
    int             ind = -1;

    for (i = 0; i < nsps; i++)
    {
        if (strcmp(chemn, chemtbl[i].name) == 0)
        {
            ind = i;
            break;
        }
    }

    return ind;
}
