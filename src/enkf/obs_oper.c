#include "pihm.h"
#include "enkf.h"

void DisOper (obs_struct *obs, pihm_struct pihm)
{
    double          dist;
    double          dist_min = 999999999.9;
    int             ind_min;
    int             i;

    for (i = 0; i < pihm->numriv; i++)
    {
        dist = sqrt ((pihm->riv[i].topo.x - obs->x) * (pihm->riv[i].topo.x - obs->x) +
            (pihm->riv[i].topo.y - obs->y) * (pihm->riv[i].topo.y - obs->y));
        if (dist < dist_min)
        {
            ind_min = i;
            dist_min = dist;
        }
    }

    obs->var_ind = (int *) malloc (sizeof (int));
    obs->weight = (double *) malloc (pihm->numriv * sizeof (double));
    obs->nctrl = 1;
    obs->k = (double *) malloc (obs->nctrl * sizeof (double));
    obs->b = (double *) malloc (obs->nctrl * sizeof (double));

    obs->var_ind[0] = 10;

    for (i = 0; i < pihm->numriv; i++)
    {
        obs->weight[i] = 0.0;
    }
    obs->weight[ind_min] = 24.0 * 3600.0;

    obs->k[0] = 1.0;
    obs->b[0] = 0.0;
}

