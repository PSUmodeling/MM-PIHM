#include "pihm.h"

void DayMet (const stor_struct *stor, daily_struct *daily, int metday)
{
    int             k;

    if (stor->flag[metday] == 0)
    {
        printf ("ERROR: BGC forcing of the %dth day is not available!\n", metday + 1);
        fflush (stdout);
        PihmExit (1);
    }

    daily->dayl = stor->dayl[metday];
    daily->prev_dayl = stor->prev_dayl[metday];
    daily->avg_q2d = stor->q2d[metday];
    daily->avg_sfcprs = stor->sfcprs[metday];
    daily->avg_albedo = stor->albedo[metday];
    daily->avg_ch = stor->ch[metday];

    daily->tmax = stor->tmax[metday];
    daily->tmin = stor->tmin[metday];
    daily->avg_sfctmp = stor->sfctmp[metday];
    daily->tday = stor->tday[metday];
    daily->tnight = stor->tnight[metday];
    for (k = 0; k < MAXLYR; k++)
    {
        daily->avg_stc[k] = stor->stc[k][metday];
    }

    daily->avg_soldn = stor->soldn[metday];

    daily->avg_surf = stor->surf[metday];
    daily->avg_unsat = stor->unsat[metday];
    daily->avg_gw = stor->gw[metday];
    for (k = 0; k < MAXLYR; k++)
    {
        daily->avg_sh2o[k] = stor->sh2o[k][metday];
    }

    for (k = 0; k < NUM_EDGE; k++)
    {
        daily->avg_ovlflow[k] = stor->surfflx[k][metday];
        daily->avg_subsurf[k] = stor->subsurfflx[k][metday];
    }
}

void RiverDayMet (const river_stor_struct *stor, river_daily_struct *daily, int metday)
{
    int             k;

    if (stor->flag[metday] == 0)
    {
        printf ("ERROR: BGC forcing of the %dth day is not available!\n", metday + 1);
        fflush (stdout);
        PihmExit (1);
    }

    daily->avg_stage = stor->stage[metday];
    daily->avg_gw = stor->gw[metday];

    for (k = 0; k < NUM_RIVFLX; k++)
    {
        daily->avg_rivflow[k] = stor->rivflow[k][metday];
    }
}
