#include "pihm.h"

void metarr_init (metarr_struct *metarr, int start_time, int end_time)
{
    /*
     * Generate meteorological forcing array for spin-up
     */
    int             j, k;
    int             length;

    length = (end_time - start_time) / 24 / 3600;

    metarr->tmax = (double *)malloc (length * sizeof (double));
    metarr->tmin = (double *)malloc (length * sizeof (double));
    metarr->prcp = (double *)malloc (length * sizeof (double));
    metarr->vpd = (double *)malloc (length * sizeof (double));
    metarr->q2d = (double *)malloc (length * sizeof (double));
    metarr->swavgfd = (double *)malloc (length * sizeof (double));
    metarr->par = (double *)malloc (length * sizeof (double));
    metarr->dayl = (double *)malloc (length * sizeof (double));
    metarr->prev_dayl = (double *)malloc (length * sizeof (double));
    metarr->tavg = (double *)malloc (length * sizeof (double));
    metarr->tday = (double *)malloc (length * sizeof (double));
    metarr->tnight = (double *)malloc (length * sizeof (double));
    metarr->tsoil = (double *)malloc (length * sizeof (double));
    metarr->swc = (double *)malloc (length * sizeof (double));
    metarr->pa = (double *)malloc (length * sizeof (double));
    metarr->soilw = (double *)malloc (length * sizeof (double));
    metarr->sw_alb = (double *)malloc (length * sizeof (double));
    metarr->gl_bl = (double *)malloc (length * sizeof (double));
    metarr->flag = (int *)malloc (length * sizeof (int));

    for (k = 0; k < 4; k++)
    {
        metarr->latflux[k] = (double *)malloc (length * sizeof (double));
    }

    for (j = 0; j < length; j++)
    {
        metarr->flag[j] = 0;
    }
}

void Save2MetArr (pihm_struct pihm, int t, int start_time, int end_time)
{
    int             i, k;
    int             ind;
    metarr_struct  *metarr;
    daily_struct   *daily;

    if (t > start_time)
    {
        ind = (t - start_time) / 24 / 3600 - 1; /* t is already the end of day */

        for (i = 0; i < pihm->numele; i++)
        {
            metarr = &(pihm->elem[i].metarr);
            daily = &(pihm->elem[i].daily);

            metarr->dayl[ind] = daily->ps.dayl;
            metarr->prev_dayl[ind] = daily->ps.prev_dayl;

            metarr->tmax[ind] = daily->es.tmax - 273.15;
            metarr->tmin[ind] = daily->es.tmin - 273.15;
            metarr->tavg[ind] = daily->es.sfctmp - 273.15;
            metarr->tday[ind] = daily->es.tday - 273.15;
            metarr->tnight[ind] = daily->es.tnight - 273.15;

            metarr->q2d[ind] = daily->ps.q2d;
            metarr->pa[ind] = daily->ps.sfcprs;
            metarr->swavgfd[ind] = daily->ef.soldn;
            metarr->par[ind] = metarr->swavgfd[ind] * RAD2PAR;

            metarr->tsoil[ind] = daily->es.stc[0] - 273.15;
            metarr->swc[ind] = daily->ws.sh2o[0];
            metarr->soilw[ind] = (daily->ws.surf
                + daily->ws.unsat * pihm->elem[i].soil.porosity
                + daily->ws.gw * pihm->elem[i].soil.porosity)
                * 1000.0;

            metarr->sw_alb[ind] = daily->ps.albedo;
            metarr->gl_bl[ind] = daily->ps.ch;

            for (k = 0; k < 3; k++)
            {
                /* Convert from m3/s to kg/m2/d */
                metarr->latflux[k][ind] = 1000.0
                    * (daily->wf.fluxsub[k] + daily->wf.fluxsurf[k])
                    * 24.0 * 3600.0 / pihm->elem[i].topo.area;
            }
            metarr->latflux[3][ind] = 0.0;

            metarr->flag[ind] = 1;
        }

        for (i = 0; i < pihm->numriv; i++)
        {
            metarr = &pihm->riv[i].metarr;
            daily = &pihm->riv[i].daily;

            metarr->soilw[ind] = (daily->ws.surf
                + daily->ws.gw * pihm->riv[i].matl.porosity) * 1000.0;

            /* Convert from m3/s to kg/m2/d */
            for (k = 0; k < 4; k++)
            {
                metarr->latflux[k][ind] = 1000.0 * (daily->wf.fluxsurf[k] + daily->wf.fluxsub[k]) * 24.0 * 3600.0 / pihm->riv[i].topo.area;
            }

            metarr->flag[ind] = 1;
        }

    }
}
