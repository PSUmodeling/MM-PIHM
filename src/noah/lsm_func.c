#include "pihm.h"

int FindWaterTable(int nlayers, double gw, const double soil_depth[], double sat_depth[])
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

void DefineSoilDepths(int nsoil_std, double total_depth, const double soil_depth_std[], int *nlayers,
    double soil_depth[], double zsoil[])
{
    int             kz, k;
    double          zsoil_std[MAXLYR];

    zsoil_std[0] = soil_depth_std[0];

    for (kz = 1; kz < MAXLYR; kz++)
    {
        zsoil_std[kz] = zsoil_std[kz - 1] + soil_depth_std[kz];
    }

    if (total_depth <= zsoil_std[0])
    {
        soil_depth[0] = total_depth;
        *nlayers = 1;
        for (kz = 1; kz < MAXLYR; kz++)
        {
            soil_depth[kz] = BADVAL;
        }
    }
    else if (total_depth <= zsoil_std[nsoil_std - 1])
    {
        for (kz = 1; kz < nsoil_std + 1; kz++)
        {
            if (total_depth <= zsoil_std[kz])
            {
                for (k = 0; k < kz; k++)
                {
                    soil_depth[k] = soil_depth_std[k];
                }
                soil_depth[kz] = total_depth - zsoil_std[kz - 1];
                *nlayers = kz + 1;

                // The following calculations guarantee that each layer is thicker than the layer on top
                if (soil_depth[kz] < soil_depth[kz - 1])
                {
                    soil_depth[kz - 1] += soil_depth[kz];
                    soil_depth[kz] = BADVAL;
                    (*nlayers)--;
                }
                for (k = kz + 1; k < MAXLYR; k++)
                {
                    soil_depth[k] = BADVAL;
                }
                break;
            }
        }
    }
    else
    {
        for (kz = 0; kz < nsoil_std; kz++)
        {
            soil_depth[kz] = soil_depth_std[kz];
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

    // Calculate depth (negative) below ground from top skin sfc to bottom of each soil layer. Note: sign of zsoil is
    // negative (denoting below ground)
    zsoil[0] = -soil_depth[0];
    for (kz = 1; kz < *nlayers; kz++)
    {
        zsoil[kz] = -soil_depth[kz] + zsoil[kz - 1];
    }
}

// Calculate transpiration from saturated zone
double GwTranspFrac(int nwtbl, int nroot, double ett, const double et[])
{
    int             kz;
    double          gw_transp = 0.0;

    if (ett > 0.0)
    {
        if (nwtbl <= nroot)
        {
            for (kz = MAX(nwtbl - 1, 0); kz < nroot; kz++)
            {
                gw_transp += et[kz];
            }

            gw_transp = gw_transp / ett;
            gw_transp = MIN(gw_transp, 1.0);
            gw_transp = MAX(gw_transp, 0.0);
        }
    }

    return gw_transp;
}

// Calculate root distribution.
// Present version assumes uniform distribution based on soil layer depths.
void RootDist(int nlayers, int nroot, const double soil_depth[], double root_dist[])
{
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

// Determine runoff from each layer
void CalcLateralFlux(const phystate_struct *ps, wflux_struct *wf)
{
    double          sattot;
    int             kz;

    sattot = 0.0;
    for (kz = 0; kz < ps->nlayers; kz++)
    {
        sattot += ps->satdpth[kz];
    }

    if (sattot <= 0.0)
    {
        wf->runoff2_lyr[ps->nlayers - 1] = wf->runoff2;
    }
    else
    {
        for (kz = 0; kz < ps->nlayers; kz++)
        {
            wf->runoff2_lyr[kz] = ps->satdpth[kz] / sattot * wf->runoff2;
        }
    }
}
