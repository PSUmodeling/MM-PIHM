#include "pihm.h"

void DayMet (const stor_struct *stor, daily_wstate_struct *daily_ws,
    daily_wflux_struct *daily_wf, daily_estate_struct *daily_es,
    daily_eflux_struct *daily_ef, daily_pstate_struct *daily_ps, int metday)
{
    int             k;

    if (stor->flag[metday] == 0)
    {
        printf ("ERROR: BGC forcing of the %dth day is not available!\n", metday + 1);
        fflush (stdout);
        PihmExit (1);
    }

    /* daylength (s) */
    daily_ps->dayl = stor->dayl[metday];
    daily_ps->prev_dayl = stor->prev_dayl[metday];
    daily_ps->avg_q2d = stor->q2d[metday];
    daily_ps->avg_sfcprs = stor->sfcprs[metday];
    daily_ps->avg_albedo = stor->albedo[metday];
    daily_ps->avg_ch = stor->ch[metday];

    daily_es->tmax = stor->tmax[metday];
    daily_es->tmin = stor->tmin[metday];
    daily_es->avg_sfctmp = stor->sfctmp[metday];
    daily_es->tday = stor->tday[metday];
    daily_es->tnight = stor->tnight[metday];
    for (k = 0; k < MAXLYR; k++)
    {
        daily_es->avg_stc[k] = stor->stc[k][metday];
    }

    daily_ef->avg_soldn = stor->soldn[metday];
    //daily_ef->par = stor->par[metday];

    daily_ws->avg_surf = stor->surf[metday];
    daily_ws->avg_unsat = stor->unsat[metday];
    daily_ws->avg_gw = stor->gw[metday];
    for (k = 0; k < MAXLYR; k++)
    {
        daily_ws->avg_sh2o[k] = stor->sh2o[k][metday];
    }

    for (k = 0; k < NUM_EDGE; k++)
    {
        daily_wf->avg_ovlflow[k] = stor->surfflx[k][metday];
        daily_wf->avg_subsurf[k] = stor->subsurfflx[k][metday];
    }
}

void RiverDayMet (const river_stor_struct *stor, river_daily_wstate_struct *daily_ws, river_daily_wflux_struct *daily_wf, int metday)
{
    int             k;

    if (stor->flag[metday] == 0)
    {
        printf ("ERROR: BGC forcing of the %dth day is not available!\n", metday + 1);
        fflush (stdout);
        PihmExit (1);
    }

    daily_ws->avg_stage = stor->stage[metday];
    daily_ws->avg_gw = stor->gw[metday];

    for (k = 0; k < NUM_RIVFLX; k++)
    {
        daily_wf->avg_rivflow[k] = stor->rivflow[k][metday];
    }
}
