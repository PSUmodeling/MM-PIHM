#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#include "../pihm.h"            /* Data Model and Variable Declarations */
#include "../noah/noah.h"
#include "bgc.h"

void bgc_spinup (char *filename, bgc_struct BGCM, Model_Data PIHM, LSM_STRUCT LSM)
{
    FILE           *soilc_file;
    FILE           *vegc_file;
    FILE           *restart_file;
    char            restart_fn[100];
    int             i, j, k;
    struct tm      *timestamp;
    time_t         *rawtime;
    int             simyr, yday, metyr, metday;
    int             ind_simyr;
    int             nmetdays;
    double          spinup_starttime;
    double          spinup_endtime;
    double          t;
    int             tmpyears;
    int             first_balance = 1;
    int             ntimesmet, nblock;
    int             metyears;
    int             steady1[PIHM->NumEle], steady2[PIHM->NumEle], rising[PIHM->NumEle], metcycle = 0, spinyears = 0;
    double          tally1[PIHM->NumEle], tally1b[PIHM->NumEle], tally2[PIHM->NumEle], tally2b[PIHM->NumEle], t1;
    int             spinup_complete[PIHM->NumEle];
    int             total_complete;
    double          naddfrac = 1.0;

    rawtime = (time_t *) malloc (sizeof (time_t));
    timestamp = (struct tm *)malloc (sizeof (struct tm));

    printf ("\n\nRunning BGC model in spin-up mode using prescribed soil moisture and soil temperature\n");
    printf ("PIHM is not running.\n");
    timestamp->tm_year = BGCM->ctrl.spinupstart;
    timestamp->tm_mon = 1;
    timestamp->tm_mday = 1;
    timestamp->tm_hour = 0;
    timestamp->tm_min = 0;
    timestamp->tm_sec = 0;
    timestamp->tm_year = timestamp->tm_year - 1900;
    timestamp->tm_mon = timestamp->tm_mon - 1;
    *rawtime = timegm (timestamp);
    spinup_starttime = (realtype) * rawtime;

    timestamp->tm_year = BGCM->ctrl.spinupend + 1;
    //printf ("year = %d\n", timestamp->tm_year);
    timestamp->tm_year = timestamp->tm_year - 1900;
    *rawtime = timegm (timestamp);
    spinup_endtime = (realtype) * rawtime;

    metarr_init (BGCM, PIHM, LSM, spinup_starttime, spinup_endtime);

    metyears = BGCM->ctrl.spinupend - BGCM->ctrl.spinupstart + 1;
    metyr = 0;

    spinyears = 0;
    metcycle = 0;
    for (i = 0; i < PIHM->NumEle; i++)
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
            for (i = 0; i < PIHM->NumEle; i++)
            {
                tally1[i] = 0.0;
                tally1b[i] = 0.0;
                tally2[i] = 0.0;
                tally2b[i] = 0.0;
            }
        }

        if (metyr == metyears)
            metyr = 0;

        printf ("Year: %6d\n", spinyears);

        for (j = 0; j < (int)((spinup_endtime - spinup_starttime) / 24. / 3600.); j++)
        {
            t = spinup_starttime + j * 24. * 3600.;
            *rawtime = t;


            //printf ("%lf %lf\n", spinup_starttime, spinup_endtime);
            for (i = 0; i < PIHM->NumEle; i++)
                //for (i = 0; i < 10; i++)
            {
                if (!steady1[i] && rising[i] && metcycle == 0)
                    naddfrac = 1. - ((double)j / (double)metyears / 365.);
                else
                    naddfrac = 0.;


                //printf ("metday = %d ELE %d\n", j + 1, i + 1);
                daymet (&BGCM->grid[i].metarr, &BGCM->grid[i].metv, j);
                //printf ("prcp %lf tnight %lf tmax %lf tmin %lf tavg %lf tday %lf tsoil %lf swc %lf vpd %lf swavgfd %lf par %lf dayl %lf prev_dayl %lf pa %lf\n", BGCM->grid[i].metv.prcp, BGCM->grid[i].metv.tnight, BGCM->grid[i].metv.tmax, BGCM->grid[i].metv.tmin, BGCM->grid[i].metv.tavg, BGCM->grid[i].metv.tday, BGCM->grid[i].metv.tsoil, BGCM->grid[i].metv.swc, BGCM->grid[i].metv.vpd, BGCM->grid[i].metv.swavgfd, BGCM->grid[i].metv.par, BGCM->grid[i].metv.dayl, BGCM->grid[i].metv.prev_dayl, BGCM->grid[i].metv.pa);
                //printf ("tnight = %lf in main \n", BGCM->grid[0].metarr.tnight[1]);
                daily_bgc (BGCM, &BGCM->grid[i], t, naddfrac, first_balance);
                //printf ("leafc = %lf\n", BGCM->grid[i].cs.leafc);
                //printf ("tnight = %lf in main \n", BGCM->grid[0].metarr.tnight[1]);

                if (metcycle == 1)
                {
                    tally1[i] += BGCM->grid[i].summary.soilc;
                    tally1b[i] += BGCM->grid[i].summary.totalc;
                }
                if (metcycle == 2)
                {
                    tally2[i] += BGCM->grid[i].summary.soilc;
                    tally2b[i] += BGCM->grid[i].summary.totalc;
                }
            }
            first_balance = 0;
        }

        spinyears = spinyears + j / 365;

        for (i = 0; i < PIHM->NumEle; i++)
        {

            if (!steady1[i] && metcycle == 2)
            {
                /* convert tally1 and tally2 to average daily soilc */
                tally1[i] /= (double)metyears *365.0;
                tally2[i] /= (double)metyears *365.0;
                rising[i] = (tally2[i] > tally1[i]);
                t1 = (tally2[i] - tally1[i]) / (double)metyears;
                steady1[i] = (fabs (t1) < SPINUP_TOLERANCE);

                //("spinyears = %d rising = %d steady1 = %d\n", spinyears, rising, steady1);
                //bgc_printf (BV_DIAG, "metcycle = %d tally1 = %lf tally2 = %lf pdif = %lf\n\n", metcycle, tally1, tally2, t1);
            }
            /* second block is after supplemental N turned off */
            else if (steady1[i] && metcycle == 2)
            {
                /* convert tally1 and tally2 to average daily soilc */
                tally1[i] /= (double)metyears *365.0;
                tally2[i] /= (double)metyears *365.0;
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
                    printf ("Ele %d spinup %d Avg daily soilc = %lf (%lf)\n", i, steady1[i] && steady2[i], tally1[i], BGCM->grid[i].summary.soilc);
                spinup_complete[i] = 1;
            }
        }

        if (metcycle == 2)
            metcycle = 0;
        else
            metcycle++;

        total_complete = 0;
        for (i = 0; i < PIHM->NumEle; i++)
            total_complete = total_complete + spinup_complete[i];

        printf ("%d elements completed spin-up, %d elements to go\n", total_complete, PIHM->NumEle - total_complete);

        /* spinup control */
        /* if this is the third pass through metcycle, do comparison */
        /* first block is during the rising phase */

        /* end of do block, test for steady state */
    } while (spinyears < BGCM->ctrl.maxspinyears || metcycle != 0);// || total_complete < PIHM->NumEle);

    soilc_file = fopen ("soilc.dat", "w");
    vegc_file = fopen("vegc.dat", "w");
    sprintf (restart_fn, "input/%s/%s.bgcinit", filename, filename);
    restart_file = fopen (restart_fn, "wb");
    for (i = 0; i < PIHM->NumEle; i++)
    {
        restart_output (&BGCM->grid[i].ws, &BGCM->grid[i].cs, &BGCM->grid[i].ns, &BGCM->grid[i].epv, &BGCM->grid[i].restart_output);
	fwrite(&(BGCM->grid[i].restart_output), sizeof(restart_data_struct), 1, restart_file);
        fprintf (soilc_file, "%lf\t", BGCM->grid[i].summary.soilc);
        fprintf (vegc_file, "%lf\t", BGCM->grid[i].summary.vegc);
    }
    fclose (soilc_file);
    fclose (vegc_file);
    fclose (restart_file);
}
