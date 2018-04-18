#include "pihm.h"

int FindWaterTable(const double *sldpth, int nsoil, double gw, double *satdpth)
{
    int             layer = -999;
    int             j;
    double          dsum = 0.0;
    double          depth;

    for (j = 0; j < MAXLYR; j++)
    {
        satdpth[j] = 0.0;
    }

    depth = 0.0;
    for (j = 0; j < nsoil; j++)
    {
        depth += sldpth[j];
    }

    if (gw <= 0.0)
    {
        layer = nsoil;
        satdpth[nsoil - 1] = 1.0E-3;
    }
    else if (gw > depth)
    {
        layer = 0;
        for (j = 0; j < nsoil; j++)
        {
            satdpth[j] = sldpth[j];
        }
    }
    else
    {
        for (j = nsoil - 1; j >= 0; j--)
        {
            if (dsum + sldpth[j] > gw)
            {
                satdpth[j] = gw - dsum;
                layer = j + 1;
                break;
            }
            else
            {
                satdpth[j] = sldpth[j];
                dsum += sldpth[j];
            }
        }
    }

    return layer;
}

int FindLayer(const double *sldpth, int nsoil, double depth)
{
    int             layer;
    int             j = 0;
    int             ind = 0;
    double          dsum = 0.0;

    if (depth <= 0.0)
    {
        layer = 0;
    }
    else
    {
        while (dsum < depth)
        {
            if (sldpth[j] < 0.0)
            {
                break;
            }
            dsum += sldpth[j];
            ind = j;
            j++;
        }
        layer = ind + 1;
        layer = (layer > nsoil) ? nsoil : layer;
    }
    return layer;
}

double TopoRadn(const topo_struct *topo, double sdir, double sdif,
    double zenith, double azimuth180)
{
    double          incidence;    /* Sun incidence angle (degree) */
    double          tcf;          /* terrain configuration factor (-) */
    double          soldown;

    azimuth180 = Mod((360.0 + azimuth180), 360.0);

    /* If the Sun is blocked, set direct solar radiation to 0.0 */
    if (zenith > topo->h_phi[(int)floor(azimuth180 / 10.0)])
    {
        sdir = 0.0;
    }

    /* Calculate Sun incidence angle */
    incidence = acos(cos(zenith * PI / 180.0) * cos(topo->slope * PI / 180.0) +
        sin(zenith * PI / 180.0) * sin(topo->slope * PI / 180.0) *
        cos((azimuth180 - topo->aspect) * PI / 180.0));
    incidence *= 180.0 / PI;
    incidence = (incidence > 90.0) ? 90.0 : incidence;

    /* Calculate terrain configuration factor
     * Dozier and Frew 1990, IEEE Transactions on Geoscience and Remote Sensing,
     * 28(5), 963--969 */
    tcf = (1.0 + cos(topo->slope * PI / 180.0)) / 2.0 - topo->svf;
    tcf = (tcf < 0.0) ? 0.0 : tcf;

    soldown = sdir * cos(incidence * PI / 180.0) +
        topo->svf * sdif + 0.2 * tcf * (sdir * cos(zenith * PI / 180.0) + sdif);
    soldown = (soldown < 0.0) ? 0.0 : soldown;

    return soldown;
}

void DefSldpth(double *sldpth, int *nsoil, double *zsoil, double total_depth,
    const double *std_sldpth, int std_nsoil)
{
    int             j, k;
    double          std_zsoil[MAXLYR];

    std_zsoil[0] = std_sldpth[0];

    for (j = 1; j < MAXLYR; j++)
    {
        std_zsoil[j] = std_zsoil[j - 1] + std_sldpth[j];
    }

    if (total_depth <= std_zsoil[0])
    {
        sldpth[0] = total_depth;
        *nsoil = 1;
        for (j = 1; j < MAXLYR; j++)
        {
            sldpth[j] = BADVAL;
        }
    }
    else if (total_depth <= std_zsoil[std_nsoil - 1])
    {
        for (j = 1; j < std_nsoil + 1; j++)
        {
            if (total_depth <= std_zsoil[j])
            {
                for (k = 0; k < j; k++)
                {
                    sldpth[k] = std_sldpth[k];
                }
                sldpth[j] = total_depth - std_zsoil[j - 1];
                *nsoil = j + 1;

                /* The following calculations guarantee that each layer is
                 * thicker than the layer on top */
                if (sldpth[j] < sldpth[j - 1])
                {
                    sldpth[j - 1] += sldpth[j];
                    sldpth[j] = BADVAL;
                    *nsoil -= 1;
                }
                for (k = j + 1; k < MAXLYR; k++)
                {
                    sldpth[k] = BADVAL;
                }
                break;
            }
        }
    }
    else
    {
        for (j = 0; j < std_nsoil; j++)
        {
            sldpth[j] = std_sldpth[j];
        }
        sldpth[std_nsoil] = total_depth - std_zsoil[std_nsoil - 1];
        *nsoil = std_nsoil + 1;
        if (sldpth[std_nsoil] < sldpth[std_nsoil - 1])
        {
            sldpth[std_nsoil - 1] += sldpth[std_nsoil];
            sldpth[std_nsoil] = BADVAL;
            *nsoil -= 1;
        }
    }

    /* Calculate depth (negative) below ground from top skin sfc to bottom of
     * each soil layer. Note: sign of zsoil is negative (denoting below ground)
     */
    zsoil[0] = -sldpth[0];
    for (k = 1; k < *nsoil; k++)
    {
        zsoil[k] = -sldpth[k] + zsoil[k - 1];
    }
}

double GwTransp(double ett, const double *et, int nwtbl, int nroot)
{
    /* Calculate transpiration from saturated zone */
    int             j;
    double          gw_transp = 0.0;

    if (ett > 0.0)
    {
        if (nwtbl <= nroot)
        {
            for (j = (nwtbl <= 0 ? 0 : nwtbl - 1); j < nroot; j++)
            {
                gw_transp += et[j];
            }

            gw_transp = gw_transp / ett;
            gw_transp = (gw_transp > 1.0) ? 1.0 : gw_transp;
            gw_transp = (gw_transp < 0.0) ? 0.0 : gw_transp;
        }
    }

    return gw_transp;
}

void RootDist(const double *sldpth, int nsoil, int nroot, double *rtdis)
{
    /* Calculate root distribution.
     * Present version assumes uniform distribution based on soil layer depths.
     */
    double          zsoil[MAXLYR];
    int             j, kz;

    zsoil[0] = -sldpth[0];
    for (kz = 1; kz < nsoil; kz++)
    {
        zsoil[kz] = -sldpth[kz] + zsoil[kz - 1];
    }

    for (j = 0; j < nroot; j++)
    {
        rtdis[j] = -sldpth[j] / zsoil[nroot - 1];
    }
}

void SunPos(const siteinfo_struct *siteinfo, int t, spa_data *spa)
{
    int             spa_result;
    pihm_t_struct   pihm_time;

    pihm_time = PIHMTime(t);

    spa->year = pihm_time.year;
    spa->month = pihm_time.month;
    spa->day = pihm_time.day;
    spa->hour = pihm_time.hour;
    spa->minute = pihm_time.minute;
    spa->second = 0;
    spa->timezone = 0;

    spa->delta_t = 67;
    spa->delta_ut1 = 0;
    spa->atmos_refract = 0.5667;

    spa->longitude = siteinfo->longitude;
    spa->latitude = siteinfo->latitude;
    spa->elevation = siteinfo->zmax;

    /* Calculate surface pressure based on FAO 1998 method (Narasimhan 2002) */
    spa->pressure =
        1013.25 * pow((293.0 - 0.0065 * spa->elevation) / 293.0, 5.26);
    spa->temperature = siteinfo->tavg;

    spa->function = SPA_ZA_RTS;
    spa_result = spa_calculate(spa);

    if (spa_result != 0)
    {
        PIHMprintf(VL_ERROR, "Error with spa error code: %d.\n", spa_result);
        PIHMexit(EXIT_FAILURE);
    }

}

void CalHum(pstate_struct *ps, estate_struct *es)
{
    const double    A2 = 17.67;
    const double    A3 = 273.15;
    const double    A4 = 29.65;
    const double    ELWV = 2.501e6;
    double          a23m4;
    const double    E0 = 611.0;
    const double    RVV = 461.0;
    const double    EPSILON = 0.622;
    double          e;
    double          esat;
    double          svp;
    const double    SVP1 = 611.2;
    const double    SVP2 = 17.67;
    const double    SVP3 = 29.65;
    const double    SVPT0 = 273.15;

    double          t2v;
    double          rho;
    double          rh;

    rh = ps->rh / 100.0;

    svp = SVP1 * exp(SVP2 * (es->sfctmp - SVPT0) / (es->sfctmp - SVP3));
    e = rh * svp;

    ps->q2 = (0.622 * e) / (ps->sfcprs - (1.0 - 0.622) * e);

    es->th2 = es->sfctmp + (0.0098 * ps->zlvl);
    //es->t1v = es->t1 * (1.0 + 0.61 * ps->q2);
    //es->th2v = es->th2 * (1.0 + 0.61 * ps->q2);
    t2v = es->sfctmp * (1.0 + 0.61 * ps->q2);
    rho = ps->sfcprs / (RD * t2v);

    a23m4 = A2 * (A3 - A4);

    esat = E0 * exp(ELWV / RVV * (1.0 / A3 - 1.0 / es->sfctmp));

    ps->q2sat = EPSILON * esat / (ps->sfcprs - (1.0 - EPSILON) * esat);

    ps->dqsdt2 = ps->q2sat * a23m4 / ((es->sfctmp - A4) * (es->sfctmp - A4));
}

double FrozRain(double prcp, double sfctmp)
{
    double          ffrozp;

    if (prcp > 0.0 && sfctmp < TFREEZ)
    {
        ffrozp = 1.0;
    }
    else
    {
        ffrozp = 0.0;
    }

    return ffrozp;
}

void CalcLatFlx(const pstate_struct *ps, wflux_struct *wf, double area)
{
    double          sattot;
    int             ks;

#if defined(_CYCLES_)
    int             k;

    for (k = 0; k < NUM_EDGE; k++)
    {
        for (ks = 0; ks < MAXLYR; ks++)
        {
            wf->smflxh[k][ks] = 0.0;
        }
    }
#endif

    /* Determine runoff from each layer */
    sattot = 0.0;
    for (ks = 0; ks < ps->nsoil; ks++)
    {
        sattot += ps->satdpth[ks];
    }

    if (sattot <= 0.0)
    {
        wf->runoff2_lyr[ps->nsoil - 1] = wf->runoff2;

#if defined(_CYCLES_)
        for (k = 0; k < NUM_EDGE; k++)
        {
            wf->smflxh[k][ps->nsoil - 1] = wf->subsurf[k] / area;
        }
#endif
    }
    else
    {
        for (ks = 0; ks < ps->nsoil; ks++)
        {
            wf->runoff2_lyr[ks] = ps->satdpth[ks] / sattot * wf->runoff2;

#if defined(_CYCLES_)
            for (k = 0; k < NUM_EDGE; k++)
            {
                wf->smflxh[k][ks] =
                    ps->satdpth[ks] / sattot * wf->subsurf[k] / area;
            }
#endif
        }
    }
}
