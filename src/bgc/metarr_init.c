#include "bgc.h"

void metarr_init (bgc_struct bgc, pihm_struct pihm, int start_time, int end_time)
{
    /*
     * Generate meteorological forcing array for spin-up
     */
    int             i, j, k;
    int             length;

    printf ("Initialize meteorological forcing array...\n");

    length = (end_time - start_time) / 24 / 3600;

    for (i = 0; i < pihm->numele; i++)
    {
        bgc->grid[i].metarr.tmax = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.tmin = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.prcp = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.vpd = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.q2d = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.swavgfd = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.par = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.dayl = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.prev_dayl = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.tavg = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.tday = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.tnight = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.tsoil = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.swc = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.pa = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.soilw = (double *)malloc (length * sizeof (double));
        bgc->grid[i].metarr.flag = (int *)malloc (length * sizeof (int));
        for (k = 0; k < 4; k++)
        {
            bgc->grid[i].metarr.latflux[k] = (double *)malloc (length * sizeof (double));
        }

        for (j = 0; j < length; j++)
        {
            bgc->grid[i].metarr.flag[j] = 0;
        }
    }

    for (i = 0; i < pihm->numriv; i++)
    {
        bgc->riv[i].metarr.flag = (int *)malloc (length * sizeof (int));
        bgc->riv[i].metarr.soilw = (double *)malloc (length * sizeof (double));

        for (k = 0; k < 4; k++)
        {
            bgc->riv[i].metarr.latflux[k] = (double *)malloc (length * sizeof (double));
        }

        for (j = 0; j < length; j++)
        {
            bgc->riv[i].metarr.flag[j] = 0;
        }
    }
}

void pihm2metarr (bgc_struct bgc, pihm_struct pihm, int t, int start_time, int end_time)
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
            metarr = &(bgc->grid[i].metarr);
            daily = &(pihm->elem[i].daily);

            metarr->dayl[ind] = daily->dayl;
            metarr->prev_dayl[ind] = daily->prev_dayl;

            metarr->tmax[ind] = daily->tmax - 273.15;
            metarr->tmin[ind] = daily->tmin - 273.15;
            metarr->tavg[ind] = daily->sfctmp - 273.15;
            metarr->tday[ind] = daily->tday - 273.15;
            metarr->tnight[ind] = daily->tnight - 273.15;

            metarr->q2d[ind] = daily->q2d;
            metarr->pa[ind] = daily->sfcprs;
            metarr->swavgfd[ind] = daily->solar;
            metarr->par[ind] = metarr->swavgfd[ind] * RAD2PAR;

            metarr->tsoil[ind] = daily->stc - 273.15;
            metarr->swc[ind] = daily->sh2o;
            metarr->soilw[ind] = (daily->surf
                + daily->unsat * pihm->elem[i].soil.porosity
                + daily->gw * pihm->elem[i].soil.porosity)
                * 1000.0;

            for (k = 0; k < 3; k++)
            {
                /* Convert from m3/s to kg/m2/d */
                metarr->latflux[k][ind] = 1000.0
                    * (daily->fluxsub[k] + daily->fluxsurf[k])
                    * 24.0 * 3600.0 / pihm->elem[i].topo.area;
            }
            metarr->latflux[3][ind] = 0.0;

            metarr->flag[ind] = 1;
        }

        for (i = 0; i < pihm->numriv; i++)
        {
            metarr = &(bgc->riv[i].metarr);
            daily = &(pihm->riv[i].daily);

            metarr->soilw[ind] = (daily->surf
                + daily->gw * pihm->riv[i].matl.porosity) * 1000.0;

            /* Convert from m3/s to kg/m2/d */
            for (k = 0; k < 4; k++)
            {
                metarr->latflux[k][ind] = 1000.0 * (daily->fluxsurf[k] + daily->fluxsub[k]) * 24.0 * 3600.0 / pihm->riv[i].topo.area;
            }

            metarr->flag[ind] = 1;
        }

    }
}
