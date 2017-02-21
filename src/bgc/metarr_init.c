#include "pihm.h"

void InitElemStor (stor_struct *stor, int start_time, int end_time)
{
    /*
     * Generate meteorological forcing array for spin-up
     */
    int             j, k;
    int             length;

    length = (end_time - start_time) / 24 / 3600;

    stor->dayl = (double *)malloc (length * sizeof (double));
    stor->prev_dayl = (double *)malloc (length * sizeof (double));
    stor->tmax = (double *)malloc (length * sizeof (double));
    stor->tmin = (double *)malloc (length * sizeof (double));
    stor->sfctmp = (double *)malloc (length * sizeof (double));
    stor->tday = (double *)malloc (length * sizeof (double));
    stor->tnight = (double *)malloc (length * sizeof (double));
    stor->q2d = (double *)malloc (length * sizeof (double));
    stor->sfcprs = (double *)malloc (length * sizeof (double));
    stor->soldn = (double *)malloc (length * sizeof (double));
    //stor->par = (double *)malloc (length * sizeof (double));
    for (k = 0; k < MAXLYR; k++)
    {
        stor->stc[k] = (double *)malloc (length * sizeof (double));
        stor->sh2o[k] = (double *)malloc (length * sizeof (double));
        stor->smc[k] = (double *)malloc (length * sizeof (double));
    }
    stor->surf = (double *)malloc (length * sizeof (double));
    stor->unsat = (double *)malloc (length * sizeof (double));
    stor->gw = (double *)malloc (length * sizeof (double));
    stor->albedo = (double *)malloc (length * sizeof (double));
    stor->ch = (double *)malloc (length * sizeof (double));
    for (j = 0; j < 3; j++)
    {
        stor->subsurfflx[j] = (double *)malloc (length * sizeof (double));
        stor->surfflx[j] = (double *)malloc (length * sizeof (double));
    }
    stor->flag = (int *)malloc (length * sizeof (int));

    for (j = 0; j < length; j++)
    {
        stor->flag[j] = 0;
    }
}

void InitRiverStor (river_stor_struct *stor, int start_time, int end_time)
{
    /*
     * Generate meteorological forcing array for spin-up
     */
    int             j;
    int             length;

    length = (end_time - start_time) / 24 / 3600;

    stor->stage = (double *)malloc (length * sizeof (double));
    stor->gw = (double *)malloc (length * sizeof (double));
    for (j = 0; j < 11; j++)
    {
        stor->rivflow[j] = (double *)malloc (length * sizeof (double));
    }
    stor->flag = (int *)malloc (length * sizeof (int));

    for (j = 0; j < length; j++)
    {
        stor->flag[j] = 0;
    }
}

void Save2Stor (pihm_struct pihm, int t, int start_time, int end_time)
{
    int             i, k;
    int             ind;
    stor_struct    *stor;
    daily_struct   *daily;

    if (t > start_time)
    {
        ind = (t - start_time) / 24 / 3600 - 1; /* t is already the end of day */

        for (i = 0; i < pihm->numele; i++)
        {
            stor = &(pihm->elem[i].stor);
            daily = &(pihm->elem[i].daily);

            stor->dayl[ind] = daily->dayl;
            stor->prev_dayl[ind] = daily->prev_dayl;

            stor->tmax[ind] = daily->tmax;
            stor->tmin[ind] = daily->tmin;
            stor->sfctmp[ind] = daily->avg_sfctmp;
            stor->tday[ind] = daily->tday;
            stor->tnight[ind] = daily->tnight;
            stor->q2d[ind] = daily->avg_q2d;
            stor->sfcprs[ind] = daily->avg_sfcprs;
            stor->soldn[ind] = daily->avg_soldn;
            //stor->par[ind] = stor->soldn[ind] * RAD2PAR;
            for (k = 0; k < MAXLYR; k++)
            {
                stor->stc[k][ind] = daily->avg_stc[k];
                stor->sh2o[k][ind] = daily->avg_sh2o[k];
                stor->smc[k][ind] = daily->avg_sh2o[k];
            }
            stor->surf[ind] = daily->avg_surf;
            stor->unsat[ind] = daily->avg_unsat;
            stor->gw[ind] = daily->avg_gw;
            stor->albedo[ind] = daily->avg_albedo;
            stor->ch[ind] = daily->avg_ch;

            for (k = 0; k < 3; k++)
            {
                stor->subsurfflx[k][ind] = daily->avg_subsurf[k];
                stor->surfflx[k][ind] = daily->avg_ovlflow[k];
            }

            stor->flag[ind] = 1;
        }

        for (i = 0; i < pihm->numriv; i++)
        {
            pihm->riv[i].stor.stage[ind] = pihm->riv[i].daily.avg_stage;
            pihm->riv[i].stor.gw[ind] = pihm->riv[i].daily.avg_gw;

            for (k = 0; k < 11; k++)
            {
                pihm->riv[i].stor.rivflow[k][ind] =
                    pihm->riv[i].daily.avg_rivflow[k];
            }

            pihm->riv[i].stor.flag[ind] = 1;
        }

    }
}
