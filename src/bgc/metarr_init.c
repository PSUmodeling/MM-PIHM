#include "pihm.h"

void InitElemStor (elem_stor_struct *stor, int start_time, int end_time)
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
    stor->par = (double *)malloc (length * sizeof (double));
    for (k = 0; k < MAXLYR; k++)
    {
        stor->stc[k] = (double *)malloc (length * sizeof (double));
        stor->sh2o[k] = (double *)malloc (length * sizeof (double));
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
        stor->riverflx[j] = (double *)malloc (length * sizeof (double));
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
    elem_stor_struct    *stor;
    elem_daily_struct   *daily;

    if (t > start_time)
    {
        ind = (t - start_time) / 24 / 3600 - 1; /* t is already the end of day */

        for (i = 0; i < pihm->numele; i++)
        {
            stor = &(pihm->elem[i].stor);
            daily = &(pihm->elem[i].daily);

            stor->dayl[ind] = daily->ps.dayl;
            stor->prev_dayl[ind] = daily->ps.prev_dayl;

            stor->tmax[ind] = daily->es.tmax;
            stor->tmin[ind] = daily->es.tmin;
            stor->sfctmp[ind] = daily->es.sfctmp;
            stor->tday[ind] = daily->es.tday;
            stor->tnight[ind] = daily->es.tnight;
            stor->q2d[ind] = daily->ps.q2d;
            stor->sfcprs[ind] = daily->ps.sfcprs;
            stor->soldn[ind] = daily->ef.soldn;
            stor->par[ind] = stor->soldn[ind] * RAD2PAR;
            for (k = 0; k < MAXLYR; k++)
            {
                stor->stc[k][ind] = daily->es.stc[k];
                stor->sh2o[k][ind] = daily->ws.sh2o[k];
            }
            stor->surf[ind] = daily->ws.surf;
            stor->unsat[ind] = daily->ws.unsat;
            stor->gw[ind] = daily->ws.gw;
            stor->albedo[ind] = daily->ps.albedo;
            stor->ch[ind] = daily->ps.ch;

            for (k = 0; k < 3; k++)
            {
                stor->subsurfflx[k][ind] = daily->wf.surf[k];
                stor->surfflx[k][ind] = daily->wf.surf[k];
            }

            //stor->tmax[ind] = daily->es.tmax - 273.15;
            //stor->tmin[ind] = daily->es.tmin - 273.15;
            //stor->tavg[ind] = daily->es.sfctmp - 273.15;
            //stor->tday[ind] = daily->es.tday - 273.15;
            //stor->tnight[ind] = daily->es.tnight - 273.15;

            //stor->q2d[ind] = daily->ps.q2d;
            //stor->pa[ind] = daily->ps.sfcprs;
            //stor->swavgfd[ind] = daily->ef.soldn;
            //stor->par[ind] = stor->swavgfd[ind] * RAD2PAR;

            //stor->tsoil[ind] = daily->es.stc[0] - 273.15;
            //stor->rootw[ind] = daily->ws.sh2o[0];
            //stor->totalw[ind] = (daily->ws.surf
            //    + daily->ws.unsat * pihm->elem[i].soil.porosity
            //    + daily->ws.gw * pihm->elem[i].soil.porosity)
            //    * 1000.0;

            //stor->sw_alb[ind] = daily->ps.albedo;
            //stor->gl_bl[ind] = daily->ps.ch;

            //for (k = 0; k < 3; k++)
            //{
            //    /* Convert from m3/s to kg/m2/d */
            //    stor->latflux[k][ind] = 1000.0
            //        * (daily->wf.subsurf[k] + daily->wf.surf[k])
            //        * 24.0 * 3600.0 / pihm->elem[i].topo.area;
            //}
            //stor->latflux[3][ind] = 0.0;

            stor->flag[ind] = 1;
        }

        for (i = 0; i < pihm->numriv; i++)
        {
            pihm->riv[i].stor.stage[ind] = pihm->riv[i].daily.ws.stage;
            pihm->riv[i].stor.gw[ind] = pihm->riv[i].daily.ws.gw;

            /* Convert from m3/s to kg/m2/d */
            for (k = 0; k < 11; k++)
            {
                pihm->riv[i].stor.riverflx[k][ind] = pihm->riv[i].daily.wf.river[k];
            }

            pihm->riv[i].stor.flag[ind] = 1;
        }

    }
}
