#include "pihm.h"

int FindWaterTable(const double *soil_depth, int nlayers, double gw, double *satdpth)
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
    for (j = 0; j < nlayers; j++)
    {
        depth += soil_depth[j];
    }

    if (gw <= 0.0)
    {
        layer = nlayers;
        satdpth[nlayers - 1] = 1.0E-3;
    }
    else if (gw > depth)
    {
        layer = 0;
        for (j = 0; j < nlayers; j++)
        {
            satdpth[j] = soil_depth[j];
        }
    }
    else
    {
        for (j = nlayers - 1; j >= 0; j--)
        {
            if (dsum + soil_depth[j] > gw)
            {
                satdpth[j] = gw - dsum;
                layer = j + 1;
                break;
            }
            else
            {
                satdpth[j] = soil_depth[j];
                dsum += soil_depth[j];
            }
        }
    }

    return layer;
}

int FindLayer(const double *soil_depth, int nlayers, double depth)
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
            if (soil_depth[j] < 0.0)
            {
                break;
            }
            dsum += soil_depth[j];
            ind = j;
            j++;
        }
        layer = ind + 1;
        layer = (layer > nlayers) ? nlayers : layer;
    }
    return layer;
}

void DefSldpth(double *soil_depth, int *nlayers, double *zsoil, double total_depth,
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
        soil_depth[0] = total_depth;
        *nlayers = 1;
        for (j = 1; j < MAXLYR; j++)
        {
            soil_depth[j] = BADVAL;
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
                    soil_depth[k] = std_sldpth[k];
                }
                soil_depth[j] = total_depth - std_zsoil[j - 1];
                *nlayers = j + 1;

                /* The following calculations guarantee that each layer is
                 * thicker than the layer on top */
                if (soil_depth[j] < soil_depth[j - 1])
                {
                    soil_depth[j - 1] += soil_depth[j];
                    soil_depth[j] = BADVAL;
                    *nlayers -= 1;
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
        for (j = 0; j < std_nsoil; j++)
        {
            soil_depth[j] = std_sldpth[j];
        }
        soil_depth[std_nsoil] = total_depth - std_zsoil[std_nsoil - 1];
        *nlayers = std_nsoil + 1;
        if (soil_depth[std_nsoil] < soil_depth[std_nsoil - 1])
        {
            soil_depth[std_nsoil - 1] += soil_depth[std_nsoil];
            soil_depth[std_nsoil] = BADVAL;
            *nlayers -= 1;
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

void RootDist(const double *soil_depth, int nlayers, int nroot, double *rtdis)
{
    /* Calculate root distribution.
     * Present version assumes uniform distribution based on soil layer depths.
     */
    double          zsoil[MAXLYR];
    int             j, kz;

    zsoil[0] = -soil_depth[0];
    for (kz = 1; kz < nlayers; kz++)
    {
        zsoil[kz] = -soil_depth[kz] + zsoil[kz - 1];
    }

    for (j = 0; j < nroot; j++)
    {
        rtdis[j] = -soil_depth[j] / zsoil[nroot - 1];
    }
}

void CalcLatFlx(const phystate_struct *ps, wflux_struct *wf)
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
