#include "pihm.h"

double Profile(int nlayers, double x[])
{
    int             k;
    double          sum = 0.0;

    for (k = 0; k < nlayers; k++)
    {
        sum += x[k];
    }

    return sum;
}