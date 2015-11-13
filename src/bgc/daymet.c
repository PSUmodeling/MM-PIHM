#include "bgc.h"

void daymet (const metarr_struct * metarr, metvar_struct * metv, int metday)
{
    int             k;

    if (metarr->flag[metday] == 0)
    {
        printf ("ERROR: BGC forcing of the %dth day is not available!\n", metday + 1);
        fflush (stdout);
        PihmExit (1);
    }

    /* convert prcp from cm --> kg/m2 */
    metv->prcp = metarr->prcp[metday];

    /* air temperature calculations (all temperatures deg C) */
    //printf ("tnight %lf\n", metv->tnight);
    //printf ("metday %d\n", metday);
    //printf ("metday in daymet = %d, %lf\n", metarr->tnight[metday]);
    metv->tnight = metarr->tnight[metday];
    metv->tmax = metarr->tmax[metday];
    metv->tmin = metarr->tmin[metday];
    metv->tavg = metarr->tavg[metday];
    metv->tday = metarr->tday[metday];
    //printf ("metday in daymet = %d, %lf\n", metarr->tnight[metday]);
    //printf ("metday in daymet = %d, %lf\n", metarr->tnight[metday]);

    metv->tsoil = metarr->tsoil[metday];
    metv->swc = metarr->swc[metday];
    metv->soilw = metarr->soilw[metday];
    metv->sw_alb = metarr->sw_alb[metday];
    metv->gl_bl = metarr->gl_bl[metday];

    /* daylight average vapor pressure deficit (Pa) */
    metv->vpd = metarr->vpd[metday];
    metv->q2d = metarr->q2d[metday];

    /* daylight average shortwave flux density (W/m2) */
    metv->swavgfd = metarr->swavgfd[metday];

    /* PAR (W/m2) */
    metv->par = metarr->par[metday];

    /* daylength (s) */
    metv->dayl = metarr->dayl[metday];
    metv->prev_dayl = metarr->prev_dayl[metday];

    metv->pa = metarr->pa[metday];

    for (k = 0; k < 3; k++)
    {
        metv->latflux[k] = metarr->latflux[k][metday];
    }
    metv->latflux[3] = 0.0;

    //printf ("prcp %lf tmax %lf tmin %lf tavg %lf tday %lf tnight %lf tsoil %lf swc %lf vpd %lf swavgfd %lf par %lf dayl %lf prev_dayl %lf pa %lf\n", metv->prcp, metv->tmax, metv->tmin, metv->tavg, metv->tday, metv->tnight, metv->tsoil, metv->swc, metv->vpd, metv->swavgfd, metv->par, metv->dayl, metv->prev_dayl, metv->pa);
}
