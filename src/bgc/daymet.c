#include "pihm.h"

void daymet (const stor_struct *stor, wstate_struct *ws, wflux_struct *wf, estate_struct *es, eflux_struct *ef, pstate_struct *ps, int metday)
{
    int             k;

    if (stor->flag[metday] == 0)
    {
        printf ("ERROR: BGC forcing of the %dth day is not available!\n", metday + 1);
        fflush (stdout);
        PihmExit (1);
    }

    /* daylength (s) */
    ps->dayl = stor->dayl[metday];
    ps->prev_dayl = stor->prev_dayl[metday];

    es->tmax = stor->tmax[metday];
    es->tmin = stor->tmin[metday];
    es->sfctmp = stor->sfctmp[metday];
    es->tday = stor->tday[metday];
    es->tnight = stor->tnight[metday];

    /* daylight average vapor pressure deficit (Pa) */
    ps->q2d = stor->q2d[metday];

    ps->sfcprs = stor->sfcprs[metday];

    /* daylight average shortwave flux density (W/m2) */
    ef->soldn = stor->soldn[metday];

    /* PAR (W/m2) */
    ef->par = stor->par[metday];

    for (k = 0; k < MAXLYR; k++)
    {
        es->stc[k] = stor->stc[k][metday];
        ws->sh2o[k] = stor->sh2o[k][metday];
    }
    ws->surf = stor->surf[metday];
    ws->stage = stor->stage[metday];
    ws->unsat = stor->unsat[metday];
    ws->gw = stor->gw[metday];
    ps->albedo = stor->albedo[metday];
    ps->ch = stor->ch[metday];

    for (k = 0; k < 3; k++)
    {
        wf->surf[k] = stor->surfflx[k][metday];
        wf->subsurf[k] = stor->subsurfflx[k][metday];
    }

    for (k = 0; k < 11; k++)
    {
        wf->river[k] = stor->riverflx[k][metday];
    }
}
