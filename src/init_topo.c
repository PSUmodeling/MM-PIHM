#include "pihm.h"

void InitTopo(const meshtbl_struct *meshtbl, elem_struct elem[])
{
    int             i, j;
    double          x[NUM_EDGE];
    double          y[NUM_EDGE];
    double          zmin[NUM_EDGE];
    double          zmax[NUM_EDGE];
#if defined(_DGW_)
    double          zbed[NUM_EDGE];
#endif

    for (i = 0; i < nelem; i++)
    {
        for (j = 0; j < NUM_EDGE; j++)
        {
            x[j] = meshtbl->x[elem[i].node[j] - 1];
            y[j] = meshtbl->y[elem[i].node[j] - 1];
            zmin[j] = meshtbl->zmin[elem[i].node[j] - 1];
            zmax[j] = meshtbl->zmax[elem[i].node[j] - 1];
#if defined(_DGW_)
            zbed[j] = meshtbl->zbed[elem[i].node[j] - 1];
#endif
        }

        elem[i].topo.area = 0.5 *
            ((x[1] - x[0]) * (y[2] - y[0]) - (y[1] - y[0]) * (x[2] - x[0]));
        /* Calculate centroid of triangle */
        elem[i].topo.x = (x[0] + x[1] + x[2]) / 3.0;
        elem[i].topo.y = (y[0] + y[1] + y[2]) / 3.0;

        elem[i].topo.zmin = (zmin[0] + zmin[1] + zmin[2]) / 3.0;
        elem[i].topo.zmax = (zmax[0] + zmax[1] + zmax[2]) / 3.0;
#if defined(_DGW_)
        elem[i].topo.zbed = (zbed[0] + zbed[1] + zbed[2]) / 3.0;
#endif
        elem[i].topo.edge[0] = sqrt(
            ((x[1] - x[2]) * (x[1] - x[2]) + (y[1] - y[2]) * (y[1] - y[2])));
        elem[i].topo.edge[1] = sqrt(
            ((x[2] - x[0]) * (x[2] - x[0]) + (y[2] - y[0]) * (y[2] - y[0])));
        elem[i].topo.edge[2] = sqrt(
            ((x[0] - x[1]) * (x[0] - x[1]) + (y[0] - y[1]) * (y[0] - y[1])));
    }

#if defined(_NOAH_)
    CalcSlopeAspect(meshtbl, elem);
#endif
}

#if defined(_NOAH_)
void CalcSlopeAspect(const meshtbl_struct *meshtbl, elem_struct elem[])
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

        /* Calculate aspect */
        if (c == 0.0)
        {
            elem[i].topo.aspect = 0.0;
        }
        else
        {
            ce = normal_vector[XCOMP] / c;
            se = normal_vector[YCOMP] / c;
            elem[i].topo.aspect = acos(ce) * 180.0 / PI;
            elem[i].topo.aspect = (se < 0.0) ?
                360.0 - elem[i].topo.aspect : elem[i].topo.aspect;
            elem[i].topo.aspect =
                Mod(360.0 - elem[i].topo.aspect + 270.0, 360.0);
        }

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
                h = (vector[ZCOMP] <= 0.0) ?
                    90.0 : atan(c / vector[ZCOMP]) * 180.0 / PI;

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
                phi1 = (se1 < 0.0) ? 360.0 - phi1 : phi1;
                phi1 = Mod(360.0 - phi1 + 270.0, 360.0);

                ce2 = edge_vector[1][XCOMP] / c2;
                se2 = edge_vector[1][YCOMP] / c2;
                phi2 = acos(ce2) * 180.0 / PI;
                phi2 = (se2 < 0.0) ? 360.0 - phi2 : phi2;
                phi2 = Mod(360.0 - phi2 + 270.0, 360.0);

                if (fabs(phi1 - phi2) > 180.0)
                {
                    ind1 = 0;
                    ind2 = (int)floor((phi1 < phi2 ? phi1 : phi2) / 10.0);
                    for (ind = ind1; ind <= ind2; ind++)
                    {
                        elem[i].topo.h_phi[ind] =
                            MIN(elem[i].topo.h_phi[ind], h);
                    }

                    ind1 = (int)floor((phi1 > phi2 ? phi1 : phi2) / 10.0);
                    ind2 = 35;
                    for (ind = ind1; ind <= ind2; ind++)
                    {
                        elem[i].topo.h_phi[ind] =
                            MIN(elem[i].topo.h_phi[ind], h);
                    }
                }
                else
                {
                    ind1 = (int)floor((phi1 < phi2 ? phi1 : phi2) / 10.0);
                    ind2 = (int)floor((phi1 > phi2 ? phi1 : phi2) / 10.0);
                    for (ind = ind1; ind <= ind2; ind++)
                    {
                        elem[i].topo.h_phi[ind] =
                            MIN(elem[i].topo.h_phi[ind], h);
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

double Mod(double a, double n)
{
    return a - n * floor(a / n);
}
#endif
