#include "pihm.h"

int FindWaterTable(int nlayers, double gw, const double soil_depth[],
    double sat_depth[])
{
    int             layer = BADVAL;
    int             kz;
    double          dsum = 0.0;
    double          depth;

    for (kz = 0; kz < MAXLYR; kz++)
    {
        sat_depth[kz] = 0.0;
    }

    depth = 0.0;
    for (kz = 0; kz < nlayers; kz++)
    {
        depth += soil_depth[kz];
    }

    if (gw <= 0.0)
    {
        layer = nlayers;
        sat_depth[nlayers - 1] = 1.0E-3;
    }
    else if (gw > depth)
    {
        layer = 0;
        for (kz = 0; kz < nlayers; kz++)
        {
            sat_depth[kz] = soil_depth[kz];
        }
    }
    else
    {
        for (kz = nlayers - 1; kz >= 0; kz--)
        {
            if (dsum + soil_depth[kz] > gw)
            {
                sat_depth[kz] = gw - dsum;
                layer = kz + 1;
                break;
            }
            else
            {
                sat_depth[kz] = soil_depth[kz];
                dsum += soil_depth[kz];
            }
        }
    }

    return layer;
}

int FindLayer(int nlayers, double depth, const double soil_depth[])
{
    int             layer;
    int             kz = 0;
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
            if (soil_depth[kz] < 0.0)
            {
                break;
            }
            dsum += soil_depth[kz];
            ind = kz;
            kz++;
        }
        layer = ind + 1;
        layer = MIN(layer, nlayers);
    }

    return layer;
}

void DefineSoilDepths(int nsoil_std, double total_depth,
    const double soil_depth_std[], int *nlayers, double soil_depth[],
    double zsoil[])
{
    int             j, k;
    double          zsoil_std[MAXLYR];

    zsoil_std[0] = soil_depth_std[0];

    for (j = 1; j < MAXLYR; j++)
    {
        zsoil_std[j] = zsoil_std[j - 1] + soil_depth_std[j];
    }

    if (total_depth <= zsoil_std[0])
    {
        soil_depth[0] = total_depth;
        *nlayers = 1;
        for (j = 1; j < MAXLYR; j++)
        {
            soil_depth[j] = BADVAL;
        }
    }
    else if (total_depth <= zsoil_std[nsoil_std - 1])
    {
        for (j = 1; j < nsoil_std + 1; j++)
        {
            if (total_depth <= zsoil_std[j])
            {
                for (k = 0; k < j; k++)
                {
                    soil_depth[k] = soil_depth_std[k];
                }
                soil_depth[j] = total_depth - zsoil_std[j - 1];
                *nlayers = j + 1;

                /* The following calculations guarantee that each layer is
                 * thicker than the layer on top */
                if (soil_depth[j] < soil_depth[j - 1])
                {
                    soil_depth[j - 1] += soil_depth[j];
                    soil_depth[j] = BADVAL;
                    (*nlayers)--;
                }
                for (k = j + 1; k < MAXLYR; k++)
                {
                    soil_depth[k] = BADVAL;
                }
                break;
            }
        }
    }
    else
    {
        for (j = 0; j < nsoil_std; j++)
        {
            soil_depth[j] = soil_depth_std[j];
        }
        soil_depth[nsoil_std] = total_depth - zsoil_std[nsoil_std - 1];
        *nlayers = nsoil_std + 1;
        if (soil_depth[nsoil_std] < soil_depth[nsoil_std - 1])
        {
            soil_depth[nsoil_std - 1] += soil_depth[nsoil_std];
            soil_depth[nsoil_std] = BADVAL;
            (*nlayers)--;
        }
    }

    /* Calculate depth (negative) below ground from top skin sfc to bottom of
     * each soil layer. Note: sign of zsoil is negative (denoting below ground)
     */
    zsoil[0] = -soil_depth[0];
    for (k = 1; k < *nlayers; k++)
    {
        zsoil[k] = -soil_depth[k] + zsoil[k - 1];
    }
}

double GwTranspFrac(int nwtbl, int nroot, double ett, const double et[])
{
    /* Calculate transpiration from saturated zone */
    int             j;
    double          gw_transp = 0.0;

    if (ett > 0.0)
    {
        if (nwtbl <= nroot)
        {
            for (j = MAX(nwtbl - 1, 0); j < nroot; j++)
            {
                gw_transp += et[j];
            }

            gw_transp = gw_transp / ett;
            gw_transp = MIN(gw_transp, 1.0);
            gw_transp = MAX(gw_transp, 0.0);
        }
    }

    return gw_transp;
}

void RootDist(int nlayers, int nroot, const double soil_depth[],
    double root_dist[])
{
    /*
     * Calculate root distribution.
     * Present version assumes uniform distribution based on soil layer depths.
     */
    double          zsoil[MAXLYR];
    int             kz;

    zsoil[0] = -soil_depth[0];
    for (kz = 1; kz < nlayers; kz++)
    {
        zsoil[kz] = -soil_depth[kz] + zsoil[kz - 1];
    }

    for (kz = 0; kz < nroot; kz++)
    {
        root_dist[kz] = -soil_depth[kz] / zsoil[nroot - 1];
    }
}

void CalcLateralFlux(const phystate_struct *ps, wflux_struct *wf)
{
    double          sattot;
    int             ks;

    /* Determine runoff from each layer */
    sattot = 0.0;
    for (ks = 0; ks < ps->nlayers; ks++)
    {
        sattot += ps->satdpth[ks];
    }

    if (sattot <= 0.0)
    {
        wf->runoff2_lyr[ps->nlayers - 1] = wf->runoff2;
    }
    else
    {
        for (ks = 0; ks < ps->nlayers; ks++)
        {
            wf->runoff2_lyr[ks] = ps->satdpth[ks] / sattot * wf->runoff2;
        }
    }
}
