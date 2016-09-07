#include "pihm.h"

void InitOper (pihm_struct pihm, enkf_struct ens)
{
    int             i;

    for (i = 0; i < ens->nobs; i++)
    {
        if (strcasecmp (ens->obs[i].name, "discharge") == 0)
        {
            ens->obs[i].type = RUNOFF_OBS;
            DisOper (&ens->obs[i], ens->var, pihm);
        }
        else if (strcasecmp (ens->obs[i].name, "skin_temperature") == 0)
        {
            ens->obs[i].type = TSKIN_OBS;
            LandSfcTmpOper (&ens->obs[i], ens->var, pihm);
        }
        else
        {
            PIHMprintf (VL_ERROR,
                "Error finding observation operator for %s.\n",
                ens->obs[i].name);
            PIHMexit (EXIT_FAILURE);
        }
    }
}

void DisOper (obs_struct *obs, var_struct *var, pihm_struct pihm)
{
    double          dist;
    double          dist_min = 999999999.9;
    int             ind_min;
    int             i;

    obs->var_ind = (int *)malloc (sizeof (int));
    obs->var_ind[0] = FindVar (var, "rivflx1");

    obs->dim = pihm->numriv;

    obs->weight = (double *)malloc (obs->dim * sizeof (double));

    obs->nlyr = 1;
    obs->k = (double **)malloc (obs->dim * sizeof (double *));
    obs->b = (double **)malloc (obs->dim * sizeof (double *));
    for (i = 0; i < obs->dim; i++)
    {
        obs->k[i] = (double *)malloc (sizeof (double));
        obs->b[i] = (double *)malloc (sizeof (double));
    }

    for (i = 0; i < pihm->numriv; i++)
    {
        dist =
            sqrt ((pihm->riv[i].topo.x - obs->x) * (pihm->riv[i].topo.x -
                obs->x) + (pihm->riv[i].topo.y -
                obs->y) * (pihm->riv[i].topo.y - obs->y));
        if (dist < dist_min)
        {
            ind_min = i;
            dist_min = dist;
        }
    }

    for (i = 0; i < pihm->numriv; i++)
    {
        obs->weight[i] = 0.0;
        obs->k[i][0] = 1.0;
        obs->b[i][0] = 0.0;
    }

    obs->weight[ind_min] = 24.0 * 3600.0;
}

void LandSfcTmpOper (obs_struct *obs, var_struct *var, pihm_struct pihm)
{
    double          dist;
    int             i;
    double          total_area = 0.0;

    obs->var_ind = (int *)malloc (sizeof (int));
    obs->var_ind[0] = FindVar (var, "t1");

    obs->dim = pihm->numele;

    obs->weight = (double *)malloc (obs->dim * sizeof (double));

    obs->nlyr = 1;
    obs->k = (double **)malloc (obs->dim * sizeof (double));
    obs->b = (double **)malloc (obs->dim * sizeof (double));
    for (i = 0; i < obs->dim; i++)
    {
        obs->k[i] = (double *)malloc (sizeof (double));
        obs->b[i] = (double *)malloc (sizeof (double));
    }

    for (i = 0; i < pihm->numele; i++)
    {
        dist =
            sqrt ((pihm->elem[i].topo.x - obs->x) * (pihm->elem[i].topo.x -
                obs->x) + (pihm->elem[i].topo.y -
                obs->y) * (pihm->elem[i].topo.y - obs->y));

        obs->k[i][0] = 1.0;
        obs->b[i][0] = 0.0;

        if (dist < obs->rad)
        {
            obs->weight[i] = pihm->elem[i].topo.area;
            total_area += pihm->elem[i].topo.area;
        }
        else
        {
            obs->weight[i] = 0.0;
        }
    }

    for (i = 0; i < pihm->numele; i++)
    {
        obs->weight[i] /= total_area;
    }
}

int FindVar (var_struct *var, char *varname)
{
    int             i;
    int             id = -999;

    for (i = 0; i < MAXVAR; i++)
    {
        if (strcasecmp (varname, var[i].name) == 0)
        {
            id = i;
            break;
        }
    }

    if (id == -999)
    {
        PIHMprintf (VL_ERROR,
            "Cannot find variable \"%s\".\n", varname);
        PIHMexit (EXIT_FAILURE);
    }

    return (id);
}

void COSMOSOper (obs_struct *obs, var_struct *var, pihm_struct pihm)
{
    double          dist;
    int             i;
    double          total_area = 0.0;

    obs->var_ind = (int *)malloc (sizeof (int));
    obs->var_ind[0] = FindVar (var, "smc0");

    obs->dim = pihm->numele;

    obs->weight = (double *)malloc (obs->dim * sizeof (double));

    obs->nlyr = 1;
    obs->k = (double **)malloc (obs->dim * sizeof (double));
    obs->b = (double **)malloc (obs->dim * sizeof (double));
    for (i = 0; i < obs->dim; i++)
    {
        obs->k[i] = (double *)malloc (sizeof (double));
        obs->b[i] = (double *)malloc (sizeof (double));
    }

    for (i = 0; i < pihm->numele; i++)
    {
        dist =
            sqrt ((pihm->elem[i].topo.x - obs->x) * (pihm->elem[i].topo.x -
                obs->x) + (pihm->elem[i].topo.y -
                obs->y) * (pihm->elem[i].topo.y - obs->y));

        obs->k[i][0] = 1.0;
        obs->b[i][0] = 0.0;

        if (dist < obs->rad)
        {
            obs->weight[i] = pihm->elem[i].topo.area;
            total_area += pihm->elem[i].topo.area;
        }
        else
        {
            obs->weight[i] = 0.0;
        }
    }

    for (i = 0; i < pihm->numele; i++)
    {
        obs->weight[i] /= total_area;
    }
}
