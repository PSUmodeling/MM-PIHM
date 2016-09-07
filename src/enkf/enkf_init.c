#include "pihm.h"

void InitEns (enkf_struct ens)
{
    int             ne;
    pihm_struct     pihm;
    N_Vector        CV_Y;       /* State Variables Vector */
    int             nsv;
    int             i, j;
    char            outputdir[MAXSTRING];

    outputdir[0] = '\0';

    pihm = (pihm_struct)malloc (sizeof *pihm);

    ReadAlloc (project, pihm);

    /* problem size */
    nsv = 3 * pihm->numele + 2 * pihm->numriv;

    CV_Y = N_VNew_Serial (nsv);

    Initialize (pihm, CV_Y);

    MapOutput (project, pihm, outputdir);

    ens->numele = pihm->numele;
    ens->numriv = pihm->numriv;
    ens->ascii = pihm->ctrl.ascii;

    /* Initialize ensemble members */
    ne = ens->ne;

    ens->member = (ensmbr_struct *)malloc (ne * sizeof (ensmbr_struct));

    /*
     * Define variable controls: vairable names, variable dimension, etc.
     */
    MapVar (ens->var, ens->numele, ens->numriv);

    InitOper (pihm, ens);

    PIHMprintf (VL_NORMAL, "Ensemble members: %d\n", ne);
    PIHMprintf (VL_NORMAL, "Default observation cycle: %-d hour(s)\n",
        ens->interval / 3600);
    PIHMprintf (VL_NORMAL, "Observations:");
    if (ens->nobs == 0)
    {
        PIHMprintf (VL_NORMAL, " none");
    }
    else
    {
        for (i = 0; i < ens->nobs - 1; i++)
        {
            PIHMprintf (VL_NORMAL, " %s,", ens->obs[i].name);
        }
        PIHMprintf (VL_NORMAL, " %s\n", ens->obs[ens->nobs - 1].name);
    }

    for (i = 0; i < ne; i++)
    {
        for (j = 0; j < MAXVAR; j++)
        {
            if (ens->var[j].dim > 0)
            {
                ens->member[i].var[j] =
                    (double *)malloc (ens->var[j].dim * sizeof (double));
            }
        }

    }

    N_VDestroy_Serial (CV_Y);

    FreeData (pihm);
    free (pihm);
}
