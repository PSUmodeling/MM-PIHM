#include "pihm.h"

void BGCSpinup (char *simulation, pihm_struct pihm, char *outputdir)
{
    FILE           *soilc_file;
    FILE           *vegc_file;
    FILE           *spinyr_file;
    FILE           *restart_file;
    FILE           *sminn_file;
    char            fn[MAXSTRING];
    char            restart_fn[MAXSTRING];
    int             i, j;
    struct tm      *timestamp;
    int             metyr;
    int             t;
    int             first_balance = 1;
    int             metyears;
    int            *steady;
    int             metcycle = 0;
    int             spinyears = 0;
    double         *tally1, *tally1b;
    double         *tally2, *tally2b;
    double          t1;
    int             total_complete;

    timestamp = (struct tm *)malloc (sizeof (struct tm));

    steady = (int *)malloc (pihm->numele * sizeof (int));
    tally1 = (double *)malloc (pihm->numele * sizeof (double));
    tally1b = (double *)malloc (pihm->numele * sizeof (double));
    tally2 = (double *)malloc (pihm->numele * sizeof (double));
    tally2b = (double *)malloc (pihm->numele * sizeof (double));

    metyears = pihm->ctrl.spinupendyear - pihm->ctrl.spinupstartyear + 1;
    metyr = 0;

    spinyears = 0;
    metcycle = 0;

    do
    {
        if (metyr == metyears)
        {
            metyr = 0;
        }

        PIHMprintf (VL_NORMAL, "Year: %6d\n", spinyears);

        for (j = 0;
            j < (pihm->ctrl.spinupend - pihm->ctrl.spinupstart) / 24 / 3600;
            j++)
        {
            t = pihm->ctrl.spinupstart + (j + 1) * 24 * 3600;

            DailyBgc (pihm, t, pihm->ctrl.spinupstart, first_balance);

            first_balance = 0;

            for (i = 0; i < pihm->numele; i++)
            {
                tally2[i] += pihm->elem[i].summary.soilc;
                tally2b[i] += pihm->elem[i].summary.totalc;
            }
        }

        metyr++;

        spinyears = spinyears + j / 365;

        for (i = 0; i < pihm->numele; i++)
        {
            if (metcycle > 0)
            {
                /* convert tally1 and tally2 to average daily soilc */
                tally2[i] /= (double)metyears * 365.0;
                t1 = (tally2[i] - tally1[i]) / (double)metyears;
                steady[i] = (fabs (t1) < SPINUP_TOLERANCE);
            }
            else
            {
                steady[i] = 0;
            }

            tally1[i] = tally2[i];
        }

        metcycle++;

        total_complete = 0;
        for (i = 0; i < pihm->numele; i++)
        {
            total_complete += steady[i];
        }

        PIHMprintf (VL_NORMAL,
            "%d elements steady, %d elements to go\n",
            total_complete, pihm->numele - total_complete);

        if (total_complete == pihm->numele)
        {
            PIHMprintf (VL_NORMAL, "All elements steady after %d year.\n",
                spinyears);
        }
    } while (spinyears < pihm->ctrl.maxspinyears);

    sprintf (fn, "%ssoilc.dat", outputdir);
    soilc_file = fopen (fn, "w");

    sprintf (fn, "%svegc.dat", outputdir);
    vegc_file = fopen (fn, "w");

    sprintf (fn, "%sspinyr.dat", outputdir);
    spinyr_file = fopen (fn, "w");

    sprintf (fn, "%ssminn.dat", outputdir);
    sminn_file = fopen (fn, "w");

    sprintf (restart_fn, "input/%s/%s.bgcic", simulation, simulation);
    restart_file = fopen (restart_fn, "wb");
    for (i = 0; i < pihm->numele; i++)
    {
        RestartOutput (&pihm->elem[i].cs, &pihm->elem[i].ns,
            &pihm->elem[i].epv, &pihm->elem[i].restart_output);
        fwrite (&(pihm->elem[i].restart_output), sizeof (bgcic_struct), 1,
            restart_file);
        fprintf (soilc_file, "%lf\t", pihm->elem[i].summary.soilc);
        fprintf (vegc_file, "%lf\t", pihm->elem[i].summary.vegc);
        fprintf (sminn_file, "%lf\t", pihm->elem[i].ns.sminn);
    }
    for (i = 0; i < pihm->numriv; i++)
    {
        fwrite (&pihm->riv[i].ns.sminn, sizeof (double), 1, restart_file);
        fprintf (sminn_file, "%lf\t", pihm->riv[i].ns.sminn);
    }

    fclose (sminn_file);
    fclose (soilc_file);
    fclose (vegc_file);
    fclose (restart_file);
    fclose (spinyr_file);
    free (steady);
    free (tally1);
    free (tally1b);
    free (tally2);
    free (tally2b);
}
