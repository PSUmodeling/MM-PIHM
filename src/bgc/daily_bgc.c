/*
bgc.c
Core BGC model logic

Includes in-line output handling routines that write to daily and annual
output files. This is the only library module that has external
I/O connections, and so it is the only module that includes bgc_io.h.

*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
Biome-BGC version 4.2 (final release)
See copyright.txt for Copyright information

Revisions since 4.1.2
	Merged spinup_bgc.c with bgc.c to eliminate
	code duplication
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*/

#include "bgc.h"

void daily_bgc(bgc_struct BGCM, bgc_grid *grid, double t)
{
    metvar_struct   *metv;
    co2control_struct *co2;
    ndepcontrol_struct *ndepctrl;
    control_struct  *ctrl;    
    struct tm      *timestamp;
    time_t          *rawtime;

    /* miscelaneous variables for program control in main */
    int simyr, yday, metyr, metday;
    int first_balance;
    int annual_alloc;
    int outv;
    int i, nmetdays;
    double tair_avg, tdiff;
    int dayout;

    double daily_ndep, daily_nfix, ndep_scalar, ndep_diff, ndep;
    int ind_simyr;

    metv = &grid->metv;
    co2 = &BGCM->co2;
    ndepctrl = &BGCM->ndepctrl;
    ctrl = &BGCM->ctrl;

    printf ("BGC daily cycle\n");
    rawtime = (time_t *) malloc (sizeof (time_t));
    *rawtime = (int)t;
    timestamp = gmtime (rawtime);

    printf ("BGC Time = %4.4d-%2.2d-%2.2d %2.2d:%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);

    /* Get co2 and ndep */
    if (ctrl->spinup == 1)      /* Spinup mode */
        metv->co2 = co2->co2ppm;
    else                /* Model mode */
    {
        /* atmospheric CO2 and Ndep handling */
        if (!(co2->varco2))
        {
            /* constant CO2, constant Ndep */
            metv->co2 = co2->co2ppm;
            daily_ndep = ndepctrl->ndep / 365.0;
            daily_nfix = ndepctrl->nfix / 365.0;
        }
        else 
        {
            /* when varco2 = 1, use file for co2 */
            if (co2->varco2 == 1)
                metv->co2 = get_co2(BGCM->Forcing[CO2_TS][0], t);
            if (metv->co2 < -999)
            {
                printf ("Error finding CO2 value on %4.4d-%2.2d-%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday);
                exit (1);
            }

            /* when varco2 = 2, use the constant CO2 value, but can vary Ndep */
            if (co2->varco2 == 2)
                metv->co2 = co2->co2ppm;

            if (ndepctrl->varndep == 0)
            {
                /* increasing CO2, constant Ndep */
                daily_ndep = ndepctrl->ndep / 365.0;
                daily_nfix = ndepctrl->nfix / 365.0;
            }
            else
            {
                daily_ndep = get_ndep(BGCM->Forcing[NDEP_TS][0], t);
                daily_nfix = ndepctrl->nfix / 365.0;
                if(daily_ndep < -999)
                {
                    printf("Error finding NDEP %4.4d-%2.2d-%2.2d\n", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday);
                    exit (1);
                }
                else
                {
                    daily_ndep = daily_ndep / 365.0;	
                }
            }
        }
    }
    printf ("co2 = %lf, ndep = %lf, nfix = %lf\n", metv->co2, daily_ndep, daily_nfix);
}
