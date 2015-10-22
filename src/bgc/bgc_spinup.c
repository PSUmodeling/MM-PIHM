#include "bgc.h"

void BgcSpinup (char *simulation, bgc_struct bgc, pihm_struct pihm, lsm_struct noah)
{
    FILE           *soilc_file;
    FILE           *vegc_file;
    FILE           *spinyr_file;
    FILE           *restart_file;
    FILE           *sminn_file;
    char            restart_fn[MAXSTRING];
    int             i, j;
    struct tm      *timestamp;
    int             metyr;
    int             t;
    int             first_balance = 1;
    int             metyears;
    int            *steady1, *steady2;
    int            *rising;
    int             metcycle = 0;
    int             spinyears = 0;
    double         *tally1, *tally1b;
    double         *tally2, *tally2b;
    double          t1;
    int            *spinup_complete;
    int            *spinup_year;
    int             total_complete;
    double         *naddfrac;

    timestamp = (struct tm *)malloc (sizeof (struct tm));

    steady1 = (int *)malloc (pihm->numele * sizeof (int));
    steady2 = (int *)malloc (pihm->numele * sizeof (int));
    rising = (int *)malloc (pihm->numele * sizeof (int));

    spinup_complete = (int *)malloc (pihm->numele * sizeof (int));
    spinup_year = (int *)malloc (pihm->numele * sizeof (int));

    tally1 = (double *)malloc (pihm->numele * sizeof (double));
    tally1b = (double *)malloc (pihm->numele * sizeof (double));
    tally2 = (double *)malloc (pihm->numele * sizeof (double));
    tally2b = (double *)malloc (pihm->numele * sizeof (double));

    naddfrac = (double *)malloc (pihm->numele * sizeof (double));

    metyears = bgc->ctrl.spinupendyear - bgc->ctrl.spinupstartyear + 1;
    metyr = 0;

    spinyears = 0;
    metcycle = 0;

    for (i = 0; i < pihm->numele; i++)
    {
        spinup_complete[i] = 0;
        steady1[i] = 0;
        steady2[i] = 0;
        rising[i] = 1;
    }

    do
    {
        if (metcycle == 0)
        {
            for (i = 0; i < pihm->numele; i++)
            {
                tally1[i] = 0.0;
                tally1b[i] = 0.0;
                tally2[i] = 0.0;
                tally2b[i] = 0.0;
            }
        }

        if (metyr == metyears)
        {
            metyr = 0;
        }

        printf ("Year: %6d\n", spinyears);

        for (j = 0; j < (bgc->ctrl.spinupend - bgc->ctrl.spinupstart) / 24 / 3600; j++)
        {
            t = bgc->ctrl.spinupstart + (j + 1) * 24 * 3600;

            for (i = 0; i < pihm->numele; i++)
            {
                naddfrac[i] = 1.0;

                if (!steady1[i] && rising[i] && metcycle == 0)
                {
                    naddfrac[i] = 1.0 - ((double) j / (double) metyears / 365.0);
                }
                else
                {
                    naddfrac[i] = 0.0;
                }

                naddfrac[i] = 0.0;

                //daymet (&bgc->grid[i].metarr, &bgc->grid[i].metv, j);
                //bgc->grid[i].ws.soilw = bgc->grid[i].metv.soilw;
                //bgc->grid[i].epv.vwc = bgc->grid[i].metv.swc;
            }

            //for (i = 0; i < pihm->numriv; i++)
            //{
            //    bgc->riv[i].soilw = bgc->riv[i].metarr.soilw[j];

            //    for (k = 0; k < 4; k++)
            //    {
            //        bgc->riv[i].metv.latflux[k] = bgc->riv[i].metarr.latflux[k][j];
            //    }
            //}

            DailyBgc (bgc, pihm->numele, pihm->numriv, t, bgc->ctrl.spinupstart, naddfrac, first_balance);

            for (i = 0; i < pihm->numele; i++)
            {
                if (metcycle == 1)
                {
                    tally1[i] += bgc->grid[i].summary.soilc;
                    tally1b[i] += bgc->grid[i].summary.totalc;
                }
                if (metcycle == 2)
                {
                    tally2[i] += bgc->grid[i].summary.soilc;
                    tally2b[i] += bgc->grid[i].summary.totalc;
                }
            }
            first_balance = 0;
        }

        spinyears = spinyears + j / 365;

        for (i = 0; i < pihm->numele; i++)
        {

            if (!steady1[i] && metcycle == 2)
            {
                /* convert tally1 and tally2 to average daily soilc */
                tally1[i] /= (double) metyears * 365.0;
                tally2[i] /= (double) metyears * 365.0;
                rising[i] = (tally2[i] > tally1[i]);
                t1 = (tally2[i] - tally1[i]) / (double) metyears;
                steady1[i] = (fabs (t1) < SPINUP_TOLERANCE);
            }
            /* second block is after supplemental N turned off */
            else if (steady1[i] && metcycle == 2)
            {
                /* convert tally1 and tally2 to average daily soilc */
                tally1[i] /= (double) metyears * 365.0;
                tally2[i] /= (double) metyears * 365.0;
                t1 = (tally2[i] - tally1[i]) / (double)metyears;
                steady2[i] = (fabs (t1) < SPINUP_TOLERANCE);

                /* if rising above critical rate, back to steady1=0 */
                if (t1 > SPINUP_TOLERANCE)
                {
                    steady1[i] = 0;
                    rising[i] = 1;
                }
            }
            if (steady1[i] && steady2[i])
            {
                if (spinup_complete[i] == 0)
                {
                    printf ("Ele %d spinup %d Avg daily soilc = %lf (%lf)\n", i, steady1[i] && steady2[i], tally1[i], bgc->grid[i].summary.soilc);
                    spinup_year[i] = spinyears;
                }
                spinup_complete[i] = 1;
            }
        }

        if (metcycle == 2)
        {
            metcycle = 0;
        }
        else
        {
            metcycle++;
        }

        total_complete = 0;
        for (i = 0; i < pihm->numele; i++)
        {
            total_complete += spinup_complete[i];
        }

        printf ("%d elements completed spin-up, %d elements to go\n", total_complete, pihm->numele - total_complete);

        /* spinup control */
        /* if this is the third pass through metcycle, do comparison */
        /* first block is during the rising phase */

        /* end of do block, test for steady state */
    } while (spinyears < bgc->ctrl.maxspinyears || metcycle != 0);// || total_complete < PIHM->NumEle);

    soilc_file = fopen ("soilc.dat", "w");
    vegc_file = fopen("vegc.dat", "w");
    spinyr_file = fopen ("spinyr.dat", "w");
    sminn_file = fopen ("sminn.dat", "w");
    sprintf (restart_fn, "input/%s/%s.bgcinit", simulation, simulation);
    restart_file = fopen (restart_fn, "wb");
    for (i = 0; i < pihm->numele; i++)
    {
        restart_output (&bgc->grid[i].ws, &bgc->grid[i].cs, &bgc->grid[i].ns, &bgc->grid[i].epv, &bgc->grid[i].restart_output);
        fwrite(&(bgc->grid[i].restart_output), sizeof(restart_data_struct), 1, restart_file);
        fprintf (soilc_file, "%lf\t", bgc->grid[i].summary.soilc);
        fprintf (vegc_file, "%lf\t", bgc->grid[i].summary.vegc);
        fprintf (spinyr_file, "%d\t", spinup_year[i]);
        fprintf (sminn_file, "%lf\t", bgc->grid[i].ns.sminn);
    }
    for (i = 0; i < pihm->numriv; i++)
    {
        fwrite (&bgc->riv[i].sminn, sizeof (double), 1, restart_file);
        fprintf (sminn_file, "%lf\t", bgc->riv[i].sminn);
    }

    fclose (sminn_file);
    fclose (soilc_file);
    fclose (vegc_file);
    fclose (restart_file);
    fclose (spinyr_file);

    free (steady1);
    free (steady2);
    free (rising);
    free (spinup_complete);
    free (spinup_year);
    free (tally1);
    free (tally1b);
    free (tally2);
    free (tally2b);
    free (naddfrac);
}
