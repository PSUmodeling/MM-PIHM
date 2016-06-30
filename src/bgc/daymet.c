#include "pihm.h"

void daymet (const metarr_struct *metarr, wstate_struct *ws, wflux_struct *wf, estate_struct *es, eflux_struct *ef, pstate_struct *ps, int metday)
{
    int             k;

    if (metarr->flag[metday] == 0)
    {
        printf ("ERROR: BGC forcing of the %dth day is not available!\n", metday + 1);
        fflush (stdout);
        PihmExit (1);
    }

    /* convert prcp from cm --> kg/m2 */
    wf->prcp = metarr->prcp[metday];

    es->tnight = metarr->tnight[metday];
    es->tmax = metarr->tmax[metday];
    es->tmin = metarr->tmin[metday];
    es->sfctmp = metarr->tavg[metday];
    es->tday = metarr->tday[metday];

    es->stc[0] = metarr->tsoil[metday];
    ws->sh2o[0] = metarr->swc[metday];
    ws->soilw = metarr->soilw[metday];
    ps->albedo = metarr->sw_alb[metday];
    ps->ch = metarr->gl_bl[metday];

    /* daylight average vapor pressure deficit (Pa) */
    ps->q2d = metarr->q2d[metday];

    /* daylight average shortwave flux density (W/m2) */
    ef->soldn = metarr->swavgfd[metday];

    /* PAR (W/m2) */
    ef->par = metarr->par[metday];

    /* daylength (s) */
    ps->dayl = metarr->dayl[metday];
    ps->prev_dayl = metarr->prev_dayl[metday];

    ps->sfcprs = metarr->pa[metday];

    for (k = 0; k < 3; k++)
    {
        wf->fluxlat[k] = metarr->latflux[k][metday];
    }
    wf->fluxlat[3] = 0.0;
}
