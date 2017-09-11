#include "pihm.h"

int FindWT(const double *sldpth, int nsoil, double gw, double *satdpth)
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

double Mod(double a, double n)
{
    return a - n * floor(a / n);
}

double TopoRadn(double sdir, double sdif, double zenith, double azimuth180,
    double slope, double aspect, const double *h_phi, double svf)
{
    double          incidence;
    double          gvf;
    double          soldown;

    azimuth180 = Mod((360.0 + azimuth180), 360.0);

    if (zenith > h_phi[(int)floor(azimuth180 / 10.0)])
    {
        sdir = 0.0;
    }

    incidence = acos(cos(zenith * PI / 180.0) * cos(slope * PI / 180.0) +
        sin(zenith * PI / 180.0) * sin(slope * PI / 180.0) *
        cos((azimuth180 - aspect) * PI / 180.0));
    incidence *= 180.0 / PI;
    incidence = incidence > 90.0 ? 90.0 : incidence;

    gvf = (1.0 + cos(slope * PI / 180.0)) / 2.0 - svf;
    gvf = (gvf < 0.0) ? 0.0 : gvf;

    soldown = sdir * cos(incidence * PI / 180.0) +
        svf * sdif + 0.2 * gvf * (sdir * cos(zenith * PI / 180.0) + sdif);
    soldown = soldown < 0.0 ? 0.0 : soldown;

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

                /* The following calculations gurantee that each layer is
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
     * each soil layer.  note:  sign of zsoil is negative (denoting below
     * ground) */
    zsoil[0] = -sldpth[0];
    for (k = 1; k < *nsoil; k++)
    {
        zsoil[k] = -sldpth[k] + zsoil[k - 1];
    }
}

void CalcSlopeAspect(elem_struct *elem, const meshtbl_struct *meshtbl)
{
    const int       XCOMP = 0;
    const int       YCOMP = 1;
    const int       ZCOMP = 2;
    double          x[NUM_EDGE];
    double          y[NUM_EDGE];
    double          zmax[NUM_EDGE];
    double          edge_vector[2][NUM_EDGE];
    double          normal_vector[NUM_EDGE];
    double          vector[NUM_EDGE];
    double          h, c;
    double          se, ce;
    int             nodes[2];
    double          x1, y1, z1, x2, y2, z2, xc, yc, zc;
    double          c1, c2, ce1, ce2, se1, se2, phi1, phi2;
    double          integrable;
    int             ind, ind1, ind2;
    int             i, j, k;

    for (i = 0; i < nelem; i++)
    {
        for (j = 0; j < NUM_EDGE; j++)
        {
            x[j] = meshtbl->x[elem[i].node[j] - 1];
            y[j] = meshtbl->y[elem[i].node[j] - 1];
            zmax[j] = meshtbl->zmax[elem[i].node[j] - 1];
        }

        edge_vector[0][XCOMP] = x[0] - x[2];
        edge_vector[0][YCOMP] = y[0] - y[2];
        edge_vector[0][ZCOMP] = zmax[0] - zmax[2];

        edge_vector[1][XCOMP] = x[1] - x[2];
        edge_vector[1][YCOMP] = y[1] - y[2];
        edge_vector[1][ZCOMP] = zmax[1] - zmax[2];

        /* Calculate normal vector */
        normal_vector[XCOMP] = edge_vector[0][YCOMP] * edge_vector[1][ZCOMP] -
            edge_vector[0][ZCOMP] * edge_vector[1][YCOMP];
        normal_vector[YCOMP] = edge_vector[0][ZCOMP] * edge_vector[1][XCOMP] -
            edge_vector[0][XCOMP] * edge_vector[1][ZCOMP];
        normal_vector[ZCOMP] = edge_vector[0][XCOMP] * edge_vector[1][YCOMP] -
            edge_vector[0][YCOMP] * edge_vector[1][XCOMP];

        if (normal_vector[ZCOMP] < 0.0)
        {
            normal_vector[XCOMP] = -normal_vector[XCOMP];
            normal_vector[YCOMP] = -normal_vector[YCOMP];
            normal_vector[ZCOMP] = -normal_vector[ZCOMP];
        }

        /* Calculate slope */
        c = sqrt(normal_vector[XCOMP] * normal_vector[XCOMP] +
            normal_vector[YCOMP] * normal_vector[YCOMP]);
        elem[i].topo.slope = atan(c / normal_vector[ZCOMP]) * 180.0 / PI;

        /* Calculte aspect */
        ce = normal_vector[XCOMP] / c;
        se = normal_vector[YCOMP] / c;
        elem[i].topo.aspect = acos(ce) * 180.0 / PI;

        if (se < 0.0)
        {
            elem[i].topo.aspect = 360.0 - elem[i].topo.aspect;
        }

        elem[i].topo.aspect = Mod(360.0 - elem[i].topo.aspect + 270.0, 360.0);

        /*
         * Calculate sky view factor (Dozier and Frew 1990)
         */
        elem[i].topo.svf = 0.0;

        /* Calculate unobstructed angle for every 10 degrees */
        for (j = 0; j < 36; j++)
        {
            elem[i].topo.h_phi[j] = 90.0;
        }

        /* Consider every edge of every triangular grid */
        for (j = 0; j < nelem; j++)
        {
            for (k = 0; k < NUM_EDGE; k++)
            {
                switch (k)
                {
                    case 0:
                        nodes[0] = 1;
                        nodes[1] = 2;
                        break;
                    case 1:
                        nodes[0] = 0;
                        nodes[1] = 2;
                        break;
                    case 2:
                        nodes[0] = 0;
                        nodes[1] = 1;
                        break;
                }
                x1 = meshtbl->x[elem[j].node[nodes[0]] - 1];
                y1 = meshtbl->y[elem[j].node[nodes[0]] - 1];
                z1 = meshtbl->zmax[elem[j].node[nodes[0]] - 1];
                x2 = meshtbl->x[elem[j].node[nodes[1]] - 1];
                y2 = meshtbl->y[elem[j].node[nodes[1]] - 1];
                z2 = meshtbl->zmax[elem[j].node[nodes[1]] - 1];

                xc = 0.5 * (x1 + x2);
                yc = 0.5 * (y1 + y2);
                zc = 0.5 * (z1 + z2);

                vector[XCOMP] = xc - elem[i].topo.x;
                vector[YCOMP] = yc - elem[i].topo.y;
                vector[ZCOMP] = zc - elem[i].topo.zmax;
                c = sqrt(vector[XCOMP] * vector[XCOMP] +
                    vector[YCOMP] * vector[YCOMP]);
                /* Unobstructed angle of the kth edge of the jth grid */
                h = atan(c / vector[ZCOMP]) * 180.0 / PI;
                h = (h < 0.0) ? 90.0 : h;

                /* Find out which directions are blocked */
                edge_vector[0][XCOMP] = x1 - elem[i].topo.x;
                edge_vector[0][YCOMP] = y1 - elem[i].topo.y;
                edge_vector[0][ZCOMP] = z1 - elem[i].topo.zmax;
                edge_vector[1][XCOMP] = x2 - elem[i].topo.x;
                edge_vector[1][YCOMP] = y2 - elem[i].topo.y;
                edge_vector[1][ZCOMP] = z2 - elem[i].topo.zmax;

                c1 = sqrt(edge_vector[0][XCOMP] * edge_vector[0][XCOMP] +
                    edge_vector[0][YCOMP] * edge_vector[0][YCOMP]);
                c2 = sqrt(edge_vector[1][XCOMP] * edge_vector[1][XCOMP] +
                    edge_vector[1][YCOMP] * edge_vector[1][YCOMP]);

                ce1 = edge_vector[0][XCOMP] / c1;
                se1 = edge_vector[0][YCOMP] / c1;
                phi1 = acos(ce1) * 180.0 / PI;
                if (se1 < 0.0)
                {
                    phi1 = 360.0 - phi1;
                }
                phi1 = Mod(360.0 - phi1 + 270.0, 360.0);

                ce2 = edge_vector[1][XCOMP] / c2;
                se2 = edge_vector[1][YCOMP] / c2;
                phi2 = acos(ce2) * 180.0 / PI;
                if (se2 < 0.0)
                {
                    phi2 = 360.0 - phi2;
                }
                phi2 = Mod(360.0 - phi2 + 270.0, 360.0);

                if (fabs(phi1 - phi2) > 180.0)
                {
                    ind1 = 0;
                    ind2 = (int)floor((phi1 < phi2 ? phi1 : phi2) / 10.0);
                    for (ind = ind1; ind <= ind2; ind++)
                    {
                        if (h < elem[i].topo.h_phi[ind])
                        {
                            elem[i].topo.h_phi[ind] = h;
                        }
                    }

                    ind1 = (int)floor((phi1 > phi2 ? phi1 : phi2) / 10.0);
                    ind2 = 35;
                    for (ind = ind1; ind <= ind2; ind++)
                    {
                        if (h < elem[i].topo.h_phi[ind])
                        {
                            elem[i].topo.h_phi[ind] = h;
                        }
                    }
                }
                else
                {
                    ind1 = (int)floor((phi1 < phi2 ? phi1 : phi2) / 10.0);
                    ind2 = (int)floor((phi1 > phi2 ? phi1 : phi2) / 10.0);
                    for (ind = ind1; ind <= ind2; ind++)
                    {
                        if (h < elem[i].topo.h_phi[ind])
                        {
                            elem[i].topo.h_phi[ind] = h;
                        }
                    }
                }
            }
        }

        /* Calculate sky view factor (Eq. 7b) */
        for (ind = 0; ind < 36; ind++)
        {
            integrable = sin(elem[i].topo.slope * PI / 180.0) *
                cos((ind * 10.0 + 5.0 - elem[i].topo.aspect) * PI / 180.0);
            integrable *= elem[i].topo.h_phi[ind] * PI / 180.0 -
                sin(elem[i].topo.h_phi[ind] * PI / 180.0) *
                cos(elem[i].topo.h_phi[ind] * PI / 180.0);
            integrable += cos(elem[i].topo.slope * PI / 180.0) *
                pow(sin(elem[i].topo.h_phi[ind] * PI / 180.0), 2);

            elem[i].topo.svf += 0.5 / PI * integrable * 10.0 / 180.0 * PI;
        }
    }
}

double GWTransp(double ett, double *et, int nwtbl, int nroot)
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

void SunPos(int t, double latitude, double longitude, double elevation,
    double tmp, spa_data * spa)
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

    spa->longitude = longitude;
    spa->latitude = latitude;
    spa->elevation = elevation;

    /* Calculate surface pressure based on FAO 1998 method (Narasimhan 2002) */
    spa->pressure =
        1013.25 * pow((293.0 - 0.0065 * spa->elevation) / 293.0, 5.26);
    spa->temperature = tmp;

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

    ps->dqsdt2 = ps->q2sat * a23m4 / pow(es->sfctmp - A4, 2);
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

double AvgElev(elem_struct *elem)
{
    double          elev = 0.0;
    int             i;

    for (i = 0; i < nelem; i++)
    {
        elev += elem[i].topo.zmax;
    }

    elev /= (double)nelem;

    return elev;
}

double TotalArea(elem_struct *elem)
{
    double          area = 0.0;
    int             i;

    for (i = 0; i < nelem; i++)
    {
        area += elem[i].topo.area;
    }

    return area;
}

void CalcLatFlx(const pstate_struct *ps, wflux_struct *wf, double area)
{
    double          sattot;
    int             ks;

#ifdef _CYCLES_
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
    sattot = 0;
    for (ks = 0; ks < ps->nsoil; ks++)
    {
        sattot += ps->satdpth[ks];
    }

    if (sattot <= 0.0)
    {
        wf->runoff2_lyr[ps->nsoil - 1] = wf->runoff2;

#ifdef _CYCLES_
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

#ifdef _CYCLES_
            for (k = 0; k < NUM_EDGE; k++)
            {
                wf->smflxh[k][ks] =
                    ps->satdpth[ks] / sattot * wf->subsurf[k] / area;
            }
#endif
        }
    }
}
