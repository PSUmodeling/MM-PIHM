#include "pihm.h"

void DayMet (const stor_struct *stor, daily_struct *daily, int metday)
{
    int             k;

    if (stor->flag[metday] == 0)
    {
        PIHMprintf (VL_ERROR,
            "Error finding BGC forcing of the %dth.\n", metday + 1);
        PIHMexit (EXIT_FAILURE);
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

    daily->prev_surf = stor->prev_surf[metday];
    daily->avg_surf = stor->avg_surf[metday];
    daily->surf = stor->surf[metday];

    daily->prev_unsat = stor->prev_unsat[metday];
    daily->avg_unsat = stor->avg_unsat[metday];
    daily->unsat = stor->unsat[metday];

    daily->prev_gw = stor->prev_gw[metday];
    daily->avg_gw = stor->avg_gw[metday];
    daily->gw = stor->gw[metday];

    for (k = 0; k < MAXLYR; k++)
    {
        daily->avg_sh2o[k] = stor->sh2o[k][metday];
        daily->prev_smc[k] = stor->prev_smc[k][metday];
        daily->avg_smc[k] = stor->avg_smc[k][metday];
        daily->smc[k] = stor->smc[k][metday];
    }

    for (k = 0; k < NUM_EDGE; k++)
    {
        daily->avg_ovlflow[k] = stor->surfflx[k][metday];
        daily->avg_subsurf[k] = stor->subsurfflx[k][metday];
    }
}

void RiverDayMet (const river_stor_struct *stor, river_daily_struct *daily,
    int metday)
{
    int             k;

    if (stor->flag[metday] == 0)
    {
        PIHMprintf (VL_ERROR,
            "Error finding BGC forcing of the %dth day.\n", metday + 1);
        PIHMexit (EXIT_FAILURE);
    }

    daily->prev_stage = stor->prev_stage[metday];
    daily->avg_stage = stor->avg_stage[metday];
    daily->stage = stor->stage[metday];

    daily->prev_gw = stor->prev_gw[metday];
    daily->avg_gw = stor->avg_gw[metday];
    daily->gw = stor->gw[metday];

    for (k = 0; k < NUM_RIVFLX; k++)
    {
        daily->avg_rivflow[k] = stor->rivflow[k][metday];
    }
}
