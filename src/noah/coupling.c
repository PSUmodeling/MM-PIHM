/*****************************************************************************
 * File		: coupling.c
 * Function	: Coupling between PIHM and Noah LSM
 * Version	: August, 2014
 ****************************************************************************/

#include "pihm.h"
#include "noah.h"
#include "spa.h"

void PIHMxNoah (int t, double stepsize, pihm_struct pihm, lsm_struct noah)
{
    grid_struct      *grid;

    const double    A2 = 17.67;
    const double    A3 = 273.15;
    const double    A4 = 29.65;
    const double    ELWV = 2.501e6;
    double          a23m4;
    const double    E0 = 611.0;
    const double    RVV = 461.0;
    const double    EPSILON = 0.622;
    double          e;
    double          svp;
    const double    SVP1 = 611.2;
    const double    SVP2 = 17.67;
    const double    SVP3 = 29.65;
    const double    SVPT0 = 273.15;

    double          t1v;
    double          th2v;
    double          t2v;
    double          rho;
    double          es;
    double          rh;
    int             i, j;
    int             spa_result;

    double          etsat;
    double          fcr, acrt, dice, sum;
    int             ialp1;
    int             jj, k;
    int             cvfrz = 3;

    double          soldown, sdir, sdif;
    time_t          rawtime;
    struct tm      *timestamp;

    spa_data        spa;
    rawtime = (int)t;
    timestamp = gmtime (&rawtime);

    spa.year = timestamp->tm_year + 1900;
    spa.month = timestamp->tm_mon + 1;
    spa.day = timestamp->tm_mday;
    spa.hour = timestamp->tm_hour;
    spa.minute = timestamp->tm_min;
    spa.second = timestamp->tm_sec;
    spa.timezone = 0;

    spa.delta_t = 67;
    spa.delta_ut1 = 0;
    spa.atmos_refract = 0.5667;

    spa.longitude = noah->longitude;
    spa.latitude = noah->latitude;

    spa.elevation = 0.0;
    for (i = 0; i < pihm->numele; i++)
    {
        spa.elevation = spa.elevation + pihm->elem[i].topo.zmax;
    }
    spa.elevation = spa.elevation / (double)pihm->numele;

    /* Calculate surface pressure based on fao 1998 method (narasimhan 2002) */
    spa.pressure = 1013.25 * pow ((293.0 - 0.0065 * spa.elevation) / 293.0, 5.26);
    spa.temperature = noah->genprmt.tbot_data;

    spa.function = SPA_ZA;
    spa_result = spa_calculate (&spa);

    if (spa_result != 0)
    {
        printf ("spa error code: %d\n", spa_result);
        exit (1);
    }

    spa.azimuth180 = mod ((360.0 + spa.azimuth180), 360.0);

    ApplyRadnForcing (&noah->forcing, t);

    for (i = 0; i < pihm->numele; i++)
    {
        grid = &(noah->grid[i]);

        /* Read time step */
        grid->dt = (double) stepsize;

        /* Read forcing */
        grid->sfcspd = *pihm->elem[i].forc.meteo[SFCSPD_TS];
        grid->sfctmp = *pihm->elem[i].forc.meteo[SFCTMP_TS];
        rh = *pihm->elem[i].forc.meteo[RH_TS];
        grid->sfcprs = *pihm->elem[i].forc.meteo[PRES_TS] / 100.0;
        grid->longwave = *pihm->elem[i].forc.meteo[LONGWAVE_TS];
        grid->prcp = *pihm->elem[i].forc.meteo[PRCP_TS];

        /* Calculate solar radiation */
        if (noah->rad_mode > 0)
        {
            sdir = *grid->radn[SOLAR_DIR_TS];
            sdif = *grid->radn[SOLAR_DIF_TS];

            soldown = TopoRadiation (sdir, sdif, spa.zenith, spa.azimuth180, grid->slope, grid->aspect, grid->h_phi, grid->svf);
        }
        else
        {
            soldown = *pihm->elem[i].forc.meteo[SOLAR_TS];
        }

        grid->soldn = soldown;

        grid->sfcprs = grid->sfcprs * 1.0e2;

        /* Initiate lsm variables */

        rh = rh / 100.0;

        svp = SVP1 * exp (SVP2 * (grid->sfctmp - SVPT0) / (grid->sfctmp - SVP3));
        e = rh * svp;

        grid->q2 = (0.622 * e) / (grid->sfcprs - (1.0 - 0.622) * e);

        if (grid->prcp > 0.0 && grid->sfctmp < 273.15)
        {
            grid->ffrozp = 1.0;
        }
        else
        {
            grid->ffrozp = 0.0;
        }

        grid->th2 = grid->sfctmp + (0.0098 * grid->zlvl);
        t1v = grid->t1 * (1.0 + 0.61 * grid->q2);
        th2v = grid->th2 * (1.0 + 0.61 * grid->q2);
        t2v = grid->sfctmp * (1.0 + 0.61 * grid->q2);
        rho = grid->sfcprs / (RD * t2v);

        a23m4 = A2 * (A3 - A4);

        es = E0 * exp (ELWV / RVV * (1.0 / A3 - 1.0 / grid->sfctmp));
        grid->q2sat = EPSILON * es / (grid->sfcprs - (1.0 - EPSILON) * es);

        grid->dqsdt2 = grid->q2sat * a23m4 / pow (grid->sfctmp - A4, 2);

        if (grid->usemonalb)
        {
            grid->alb = 0.18;
        }
        else
        {
            grid->alb = BADVAL;
        }

        if (grid->rdlai2d)
        {
            grid->xlai = 2.0;
        }
        //else
        //    grid->xlai = badval;

        grid->shdfac = pihm->elem[i].lc.vegfrac;

#ifndef _BGC_
        if (pihm->elem[i].forc.lai_type > 0)
        {
            grid->xlai = *pihm->elem[i].forc.lai;
        }
        else
        {
            grid->xlai = MonthlyLAI (t, pihm->elem[i].lc.type);
        }
#endif

        grid->cmcmax = grid->cmcfactr * grid->xlai;

        if (grid->q1 == BADVAL)
        {
            grid->q1 = grid->q2;
        }

        SfcDifOff (&(grid->zlvl), &(grid->zlvl_wind), &(grid->z0), &t1v, &th2v, &(grid->sfcspd), &(grid->czil), &(grid->cm), &(grid->ch), &(grid->vegtyp), &(grid->isurban), &(grid->iz0tlnd));

        grid->solnet = grid->soldn * (1.0 - grid->albedo);
        grid->lwdn = grid->longwave * grid->emissi;

        grid->avginfil /= (double) (grid->dt / pihm->dt);
        for (j = 0; j < 3; j++)
        {
            grid->avgsubflux[j] /= (double) (grid->dt / pihm->dt);
        }

        grid->runoff2 = 0.0;
        for (j = 0; j < 3; j++)
        {
            grid->runoff2 += grid->avgsubflux[j] / pihm->elem[i].topo.area;
        }
        grid->infil = grid->avginfil;
        grid->nwtbl = FindLayer (grid->sldpth, grid->nsoil, pihm->elem[i].soil.depth - pihm->elem[i].gw0);
        if (grid->nwtbl > grid->nsoil)
        {
            grid->nwtbl = grid->nsoil;
        }

        grid->mac_status = pihm->elem[i].macpore_status;

        SFlx (grid);

        grid->drip = 1.0e3 * grid->drip / grid->dt;  /* convert drip from m/s to kg m{-2} s{-1} (mm/s) */

        /*
         * Transfer Noah variables to PIHM
         */
#ifdef _rt_
	pihm->ele[i].temp   = (realtype) grid->stc[2];
#endif
        pihm->elem[i].netprcp =  grid->pcpdrp;
        /*
         * ET: convert from w m-2 to m s-1 
         */
        pihm->elem[i].et[0] = grid->ec / LVH2O / 1000.0 ;
        pihm->elem[i].et[1] = grid->ett / LVH2O / 1000.0 ;
        pihm->elem[i].et[2] = grid->edir / LVH2O / 1000.0 ;

        /* Calculate transpiration from saturated zone */
        pihm->elem[i].et_from_sat = 0.0;
        etsat = 0.0;
        if (grid->ett > 0.0)
        {
            if (grid->nwtbl <= grid->nroot)
            {
                for (j = (grid->nwtbl <= 0 ? 0 : grid->nwtbl - 1); j < grid->nroot; j++)
                {
                    etsat += grid->et[j];
                }
                pihm->elem[i].et_from_sat =  etsat / grid->ett;
                pihm->elem[i].et_from_sat = (pihm->elem[i].et_from_sat > 1.0) ? 1.0 : pihm->elem[i].et_from_sat;
                pihm->elem[i].et_from_sat = (pihm->elem[i].et_from_sat < 0.0) ? 0.0 : pihm->elem[i].et_from_sat;
            }
        }
        pihm->elem[i].snow = grid->sneqv;
        pihm->elem[i].intcp = grid->cmc;

        /* Calculate surface saturation ratio for pihm infiltration */
        pihm->elem[i].sfcsat = (grid->sh2o[0] - grid->smcmin) / (grid->smcmax - grid->smcmin);
        pihm->elem[i].sfcsat = (pihm->elem[i].sfcsat > 0.0) ? pihm->elem[i].sfcsat : 0.0;
        pihm->elem[i].sfcsat = (pihm->elem[i].sfcsat < 1.0) ? pihm->elem[i].sfcsat : 1.0;

        /* Calculate infiltration reduction factor due to frozen soil */
        pihm->elem[i].fcr = 1.0;
        dice = 0.0;
        for (j = 0; j < grid->nsoil; j++)
        {
            dice = dice + grid->sldpth[j] * (grid->smc[j] - grid->sh2o[j]);
        }
        fcr = 1.0;
        if (dice > 1.0e-2)
        {
            acrt = (double) cvfrz * grid->frzx / dice;
            sum = 1.0;
            ialp1 = cvfrz - 1;
            for (j = 1; j < ialp1 + 1; j++)
            {
                k = 1;
                for (jj = j + 1; jj < ialp1; jj++)
                    k = k * jj;
                sum = sum + pow (acrt, (double)(cvfrz - j)) / (double)k;
            }
            fcr = 1.0 - exp (-acrt) * sum;
        }
        pihm->elem[i].fcr = fcr;

        grid->avginfil = 0.0;
        for (j = 0; j < 3; j++)
        {
            grid->avgsubflux[j] = 0.0;
        }
    }

}

void ApplyRadnForcing (radn_ts_struct *forcing, int t)
{
    int         i, j;
    double      radn[2];

    for (i = 0; i < forcing->nts; i++)
    {
        IntrplForcing (forcing->ts[i], t, 2, radn);

        for (j = 0; j < 2; j++)
        {
            forcing->radn[j][i] = radn[j];
        }
    }
}

void AvgFlux (lsm_struct noah, pihm_struct pihm)
{
    int         i, j;
    for (i = 0; i < pihm->numele; i++)
    {
        noah->grid[i].avginfil += pihm->elem[i].infil;
        for (j = 0; j < 3; j++)
        {
            noah->grid[i].avgsubflux[j] += pihm->elem[i].fluxsub[j];
        }
    }
}
