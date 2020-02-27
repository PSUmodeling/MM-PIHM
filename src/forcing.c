#include "pihm.h"

#if defined(_RT_)
void ApplyBc(const rttbl_struct *rttbl, forc_struct *forc, elem_struct *elem,
    river_struct *river, int t)
#else
void ApplyBc(forc_struct *forc, elem_struct *elem, river_struct *river, int t)
#endif
{
    /* Element boundary conditions */
    if (forc->nbc > 0)
    {
#if defined(_RT_)
        ApplyElemBc(rttbl, forc, elem, t);
#else
        ApplyElemBc(forc, elem, t);
#endif
    }

    /* River boundary condition */
    if (forc->nriverbc > 0)
    {
        ApplyRiverBc(forc, river, t);
    }
}

#if defined(_RT_)
void ApplyForc(forc_struct *forc, rttbl_struct *rttbl, elem_struct *elem,
    int t, int stepsize, int rad_mode, const siteinfo_struct *siteinfo)
#elif defined(_NOAH_)
void ApplyForc(forc_struct *forc, elem_struct *elem, int t, int stepsize,
    int rad_mode, const siteinfo_struct *siteinfo)
#else
void ApplyForc(forc_struct *forc, elem_struct *elem, int t, int stepsize)
#endif
{
    /* Meteorological forcing */
#if defined(_NOAH_)
    ApplyMeteoForc(forc, elem, t, stepsize, rad_mode, siteinfo);
#else
    ApplyMeteoForc(forc, elem, t, stepsize);
#endif

    /* LAI forcing */
#if defined(_BGC_) || defined(_CYCLES_)
    ApplyLai(elem);
#else
    ApplyLai(forc, elem, t);
#endif

#if defined(_RT_)
    /* Precipitation concentration */
    ApplyPrcpConc(rttbl, forc, elem, t);
#endif
}

#if defined(_RT_)
void ApplyElemBc(const rttbl_struct *rttbl, forc_struct *forc,
    elem_struct *elem, int t)
#else
void ApplyElemBc(forc_struct *forc, elem_struct *elem, int t)
#endif
{
    int             i, k;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (k = 0; k < forc->nbc; k++)
    {
#if defined(_RT_)
        IntrplForc(&forc->bc[k], t, 0, 1 + rttbl->NumStc, INTRPL);
#else
        IntrplForc(&forc->bc[k], t, 0, 1, INTRPL);
#endif
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             ind;
        int             j;
#if defined(_RT_)
        int             k;
#endif

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].attrib.bc_type[j] > 0)
            {
                ind = elem[i].attrib.bc_type[j] - 1;
                elem[i].bc.head[j] = forc->bc[ind].value[0];
#if defined(_RT_)
                for (k = 0; k < rttbl->NumStc; k++)
                {
                    elem[i].bc.conc[j][k] = forc->bc[ind].value[k + 1];
                }
#endif
            }
            else if (elem[i].attrib.bc_type[j] < 0)
            {
                ind = -elem[i].attrib.bc_type[j] - 1;
                elem[i].bc.flux[j] = forc->bc[ind].value[0];
#if defined(_RT_)
                for (k = 0; k < rttbl->NumStc; k++)
                {
                    elem[i].bc.conc[j][k] = forc->bc[ind].value[k + 1];
                }
#endif
            }


#if defined(_FBR_)
            if (elem[i].attrib.fbrbc_type[j] > 0)
            {
                ind = elem[i].attrib.fbrbc_type[j] - 1;
                elem[i].fbr_bc.head[j] = forc->bc[ind].value[0];
# if defined(_RT_)
                for (k = 0; k < rttbl->NumStc; k++)
                {
                    elem[i].fbr_bc.conc[j][k] = forc->bc[ind].value[k + 1];
                }
# endif
            }
            else if (elem[i].attrib.fbrbc_type[j] < 0)
            {
                ind = -elem[i].attrib.fbrbc_type[j] - 1;
                elem[i].fbr_bc.flux[j] = forc->bc[ind].value[0];
# if defined(_RT_)
                for (k = 0; k < rttbl->NumStc; k++)
                {
                    elem[i].fbr_bc.conc[j][k] = forc->bc[ind].value[k + 1];
                }
# endif
            }
#endif
        }
    }
}

#if defined(_NOAH_)
void ApplyMeteoForc(forc_struct *forc, elem_struct *elem, int t, int stepsize,
    int rad_mode, const siteinfo_struct *siteinfo)
#else
void ApplyMeteoForc(forc_struct *forc, elem_struct *elem, int t, int stepsize)
#endif
{
    int             i, k;
#if defined(_NOAH_)
    spa_data        spa;
#endif

    /*
     * Meteorological forcing for PIHM
     */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (k = 0; k < forc->nmeteo; k++)
    {
        IntrplForc(&forc->meteo[k], t, stepsize, NUM_METEO_VAR, INTRPL);
    }

#if defined(_NOAH_)
    /*
     * Topographic radiation for Noah
     */
    if (rad_mode > 0)
    {
        if (forc->nrad > 0)
        {
# if defined(_OPENMP)
#  pragma omp parallel for
# endif
            for (k = 0; k < forc->nrad; k++)
            {
                IntrplForc(&forc->rad[k], t, stepsize, 2, INTRPL);
            }
        }

        /* Calculate Sun position for topographic solar radiation */
        SunPos(siteinfo, t, &spa);
    }
#endif

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             ind;

        ind = elem[i].attrib.meteo_type - 1;

        elem[i].wf.prcp = forc->meteo[ind].value[PRCP_TS] / 1000.0;
        elem[i].es.sfctmp = forc->meteo[ind].value[SFCTMP_TS];
        elem[i].ps.rh = forc->meteo[ind].value[RH_TS];
        elem[i].ps.sfcspd = forc->meteo[ind].value[SFCSPD_TS];
        elem[i].ef.soldn = forc->meteo[ind].value[SOLAR_TS];
        elem[i].ef.soldn = MAX(elem[i].ef.soldn, 0.0);
#if defined(_NOAH_)
        elem[i].ef.longwave = forc->meteo[ind].value[LONGWAVE_TS];
#endif
        elem[i].ps.sfcprs = forc->meteo[ind].value[PRES_TS];

#if defined(_NOAH_)
        /* Calculate solar radiation */
        if (rad_mode > 0)
        {
            if (forc->nrad > 0)
            {
                elem[i].ef.soldir = forc->rad[ind].value[SOLDIR_TS];
                elem[i].ef.soldif = forc->rad[ind].value[SOLDIF_TS];
            }

            elem[i].ef.soldn = TopoRadn(&elem[i].topo, elem[i].ef.soldir,
                elem[i].ef.soldif, spa.zenith, spa.azimuth180);
            elem[i].ef.soldn = MAX(elem[i].ef.soldn, 0.0);
        }
#endif
    }
}

#if defined(_BGC_) || defined(_CYCLES_)
void ApplyLai(elem_struct *elem)
#else
void ApplyLai(forc_struct *forc, elem_struct *elem, int t)
#endif
{
    int             i;

#if defined(_CYCLES_)
    /*
     * Cycles coupling
     */
# if defined(_OPENMP)
#  pragma omp parallel for
# endif
    for (i = 0; i < nelem; i++)
    {
        double          ksolar;
        double          tau;

        if (CommRadIntcp(elem[i].crop) > 0.0)
        {
            ksolar = 0.5;

            tau = 1.0 - MIN(CommRadIntcp(elem[i].crop), 0.98);

            elem[i].ps.proj_lai = -log(tau) / ksolar;
        }
        else
        {
            elem[i].ps.proj_lai = 0.0;
        }
    }
#elif  _BGC_
    /*
     * BGC coupling
     */
# if defined(_OPENMP)
#  pragma omp parallel for
# endif
    for (i = 0; i < nelem; i++)
    {
# if defined(_LUMPED_)
        elem[i].ps.proj_lai =
            elem[LUMPED].cs.leafc * elem[LUMPED].epc.avg_proj_sla;
# else
        elem[i].ps.proj_lai = elem[i].cs.leafc * elem[i].epc.avg_proj_sla;
# endif
    }
#else
    /*
     * Use LAI forcing
     */
    int             k;

    if (forc->nlai > 0)
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (k = 0; k < forc->nlai; k++)
        {
            IntrplForc(&forc->lai[k], t, 0, 1, INTRPL);
        }
    }

# if defined(_OPENMP)
#  pragma omp parallel for
# endif
    for (i = 0; i < nelem; i++)
    {
        int             ind;

        if (elem[i].attrib.lai_type > 0)
        {
            ind = elem[i].attrib.lai_type - 1;

            elem[i].ps.proj_lai = forc->lai[ind].value[0];
        }
        else
        {
            elem[i].ps.proj_lai = MonthlyLai(t, elem[i].attrib.lc_type);
        }
    }
#endif
}

#if defined(_RT_)
void ApplyPrcpConc(const rttbl_struct *rttbl, forc_struct *forc,
    elem_struct elem[], int t)
{
    int             i, j;

    if (forc->PrpFlg == 2)
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (j = 0; j < forc->nprcpc; j++)
        {
            IntrplForc(&forc->prcpc[j], t, 0, NumSpc, NO_INTRPL);
        }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (i = 0; i < nelem; i++)
        {
            int             k;
            int             ind;

            ind = elem[i].attrib.prcpc_type - 1;

            for (k = 0; k < NumSpc; k++)
            {
                elem[i].prcps.t_conc[k] = forc->prcpc[ind].value[k];
            }
        }
    }
    else if (forc->PrpFlg == 1)
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (i = 0; i < nelem; i++)
        {
            int             k;

            for (k = 0; k < NumSpc; k++)
            {
                elem[i].prcps.t_conc[k] = rttbl->prcp_conc[k];
            }
        }
    }
    else if (forc->PrpFlg == 0)
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (i = 0; i < nelem; i++)
        {
            int             k;

            for (k = 0; k < NumSpc; k++)
            {
                elem[i].prcps.t_conc[k] = 0.0;
            }
        }
    }
}
#endif

void ApplyRiverBc(forc_struct *forc, river_struct *river, int t)
{
    int             i, k;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (k = 0; k < forc->nriverbc; k++)
    {
        IntrplForc(&forc->riverbc[k], t, 0, 1, INTRPL);
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             ind;

        if (river[i].attrib.riverbc_type > 0)
        {
            ind = river[i].attrib.riverbc_type - 1;
            river[i].bc.head = forc->riverbc[ind].value[0];
        }
        else if (river[i].attrib.riverbc_type < 0)
        {
            ind = -river[i].attrib.riverbc_type - 1;
            river[i].bc.flux = forc->riverbc[ind].value[0];
        }
    }
}

void IntrplForc(tsdata_struct *ts, int t, int interval, int nvrbl, int intrpl)
{
    int             j;
    int             first, middle, last;

    if (t < ts->ftime[0])
    {
        PIHMprintf(VL_ERROR, "Error finding forcing for current time step.\n");
        PIHMprintf(VL_ERROR, "Please check your forcing file.\n");
        PIHMexit(EXIT_FAILURE);
    }
    else if (t > ts->ftime[ts->length - 1])
    {
        PIHMprintf(VL_ERROR, "Error finding forcing for current time step.\n");
        PIHMprintf(VL_ERROR, "Please check your forcing file.\n");
        PIHMexit(EXIT_FAILURE);
    }

    first = 1;
    last = ts->length - 1;

    while (first <= last)
    {
        middle = (first + last) / 2;
        if (t >= ts->ftime[middle - 1] && t < ts->ftime[middle])
        {
            if (intrpl)
            {
                if (middle + 1 <= ts->length - 1 &&
                    ts->ftime[middle + 1] <= t + interval)
                {
                    int             k;
                    int             n = 0;

                    for (j = 0; j < nvrbl; j++)
                    {
                        ts->value[j] = 0.0;
                    }

                    k = (t == ts->ftime[middle - 1]) ? middle - 1 : middle;

                    while (ts->ftime[k] < t + interval && k <= ts->length - 1)
                    {
                        for (j = 0; j < nvrbl; j++)
                        {
                            ts->value[j] += ts->data[k][j];
                        }
                        n++;
                        k++;
                    }

                    for (j = 0; j < nvrbl; j++)
                    {
                        ts->value[j] /= (double)n;
                    }

                }
                else                        /* Interpolation w/o aggregation */
                {
                    for (j = 0; j < nvrbl; j++)
                    {
                        ts->value[j] =
                            ((double)(ts->ftime[middle] - t) *
                            ts->data[middle - 1][j] +
                            (double)(t - ts->ftime[middle - 1]) *
                            ts->data[middle][j]) /
                            (double)(ts->ftime[middle] - ts->ftime[middle - 1]);
                    }
                }
            }
            else                            /* No interpolation */
            {
                for (j = 0; j < nvrbl; j++)
                {
                    ts->value[j] = ts->data[middle - 1][j];
                }
            }

            break;
        }
        else if (ts->ftime[middle] > t)
        {
            last = middle - 1;
        }
        else
        {
            first = middle + 1;
        }
    }
}

double MonthlyLai(int t, int lc_type)
{
    /*
     * Monthly LAI data come from WRF MPTABLE.TBL for Noah MODIS land
     * cover categories
     */
    pihm_t_struct   pihm_time;

    double          lai_tbl[40][12] = {
        /* Evergreen Needleleaf Forest */
        {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0},
        /* Evergreen Broadleaf Forest */
        {4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5},
        /* Deciduous Needleleaf Forest */
        {0.0, 0.0, 0.0, 0.6, 1.2, 2.0, 2.6, 1.7, 1.0, 0.5, 0.2, 0.0},
        /* Deciduous Broadleaf Forest */
        {0.0, 0.0, 0.3, 1.2, 3.0, 4.7, 4.5, 3.4, 1.2, 0.3, 0.0, 0.0},
        /* Mixed Forest */
        {2.0, 2.0, 2.2, 2.6, 3.5, 4.3, 4.3, 3.7, 2.6, 2.2, 2.0, 2.0},
        /* Closed Shrubland */
        {0.0, 0.0, 0.3, 0.9, 2.2, 3.5, 3.5, 2.5, 0.9, 0.3, 0.0, 0.0},
        /* Open Shrubland */
        {0.0, 0.0, 0.2, 0.6, 1.5, 2.3, 2.3, 1.7, 0.6, 0.2, 0.0, 0.0},
        /* Woody Savanna */
        {0.2, 0.2, 0.4, 1.0, 2.4, 4.1, 4.1, 2.7, 1.0, 0.4, 0.2, 0.2},
        /* Savanna */
        {0.3, 0.3, 0.5, 0.8, 1.8, 3.6, 3.8, 2.1, 0.9, 0.5, 0.3, 0.3},
        /* Grassland */
        {0.4, 0.5, 0.6, 0.7, 1.2, 3.0, 3.5, 1.5, 0.7, 0.6, 0.5, 0.4},
        /* Permanent Wetland */
        {0.2, 0.3, 0.3, 0.5, 1.5, 2.9, 3.5, 2.7, 1.2, 0.3, 0.3, 0.2},
        /* Cropland */
        {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 1.5, 0.0, 0.0, 0.0},
        /* Urban and Built-Up */
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* Cropland/Natural Veg. Mosaic */
        {0.2, 0.3, 0.3, 0.4, 1.1, 2.5, 3.2, 2.2, 1.1, 0.3, 0.3, 0.2},
        /* Permanent Snow */
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* Barren/Sparsely Vegetated */
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* IGBP Water */
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* Unclassified */
        {999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999},
        /* Fill Value */
        {999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999},
        /* Unclassified */
        {999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999},
        /* Open Water */
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        /* Perennial Ice/Snow */
        {1.06, 1.11, 1.20, 1.29, 1.45, 1.67, 1.71, 1.63, 1.45, 1.27, 1.13,
            1.06},
        /* Developed Open Space */
        {0.83, 0.94, 1.06, 1.18, 1.86, 3.72, 4.81, 4.26, 2.09, 1.29, 1.06,
            0.94},
        /* Developed Low Intensity */
        {0.96, 1.07, 1.20, 1.35, 2.04, 3.83, 4.87, 4.35, 2.30, 1.46, 1.18,
            1.06},
        /* Developed Medium Intensity */
        {1.11, 1.22, 1.36, 1.54, 2.26, 3.97, 4.94, 4.46, 2.55, 1.66, 1.34,
            1.20},
        /* Developed High Intensity */
        {1.24, 1.34, 1.50, 1.71, 2.45, 4.09, 5.00, 4.54, 2.76, 1.82, 1.47,
            1.32},
        /* Barren Land */
        {0.03, 0.03, 0.03, 0.02, 0.02, 0.03, 0.04, 0.06, 0.09, 0.06, 0.04,
            0.03},
        /* Deciduous Forest */
        {0.62, 0.67, 0.92, 1.71, 3.42, 5.53, 6.22, 5.60, 3.83, 1.79, 0.92,
            0.67},
        /* Evergreen Forest */
        {3.38, 3.43, 3.47, 3.52, 3.78, 4.54, 4.98, 4.76, 3.87, 3.56, 3.47,
            3.43},
        /* Mixed Forest */
        {3.10, 3.26, 3.61, 4.11, 5.17, 6.73, 7.21, 6.71, 5.34, 4.09, 3.41,
            3.14},
        /* Dwarf Scrub */
        {0.24, 0.24, 0.19, 0.13, 0.15, 0.20, 0.26, 0.48, 0.70, 0.48, 0.30,
            0.24},
        /* Shrub/Scrub */
        {0.35, 0.38, 0.38, 0.38, 0.55, 1.06, 1.53, 1.53, 1.04, 0.58, 0.44,
            0.38},
        /* Grassland/Herbaceous */
        {0.70, 0.80, 0.90, 1.00, 1.60, 3.30, 4.30, 3.80, 1.80, 1.10, 0.90,
            0.80},
        /* Sedge/Herbaceous */
        {0.70, 0.80, 0.90, 1.00, 1.60, 3.30, 4.30, 3.80, 1.80, 1.10, 0.90,
            0.80},
        /* Lichens */
        {0.70, 0.80, 0.90, 1.00, 1.60, 3.30, 4.30, 3.80, 1.80, 1.10, 0.90,
            0.80},
        /* Moss */
        {0.70, 0.80, 0.90, 1.00, 1.60, 3.30, 4.30, 3.80, 1.80, 1.10, 0.90,
            0.80},
        /* Pasture/Hay */
        {0.47, 0.54, 0.60, 0.67, 1.07, 2.20, 2.87, 2.54, 1.20, 0.74, 0.60,
            0.54},
        /* Cultivated Crops */
        {0.47, 0.54, 0.60, 0.67, 1.07, 2.20, 2.87, 2.54, 1.20, 0.74, 0.60,
            0.54},
        /* Woody Wetland */
        {0.35, 0.38, 0.38, 0.38, 0.55, 1.06, 1.53, 1.53, 1.04, 0.58, 0.44,
            0.38},
        /* Emergent Herbaceous Wetland */
        {0.24, 0.24, 0.19, 0.13, 0.15, 0.20, 0.26, 0.48, 0.70, 0.48, 0.30,
            0.24}
    };

    pihm_time = PIHMTime(t);

    return lai_tbl[lc_type - 1][pihm_time.month - 1];
}

double MonthlyRl(int t, int lc_type)
{
    /*
     * Monthly roughness length data are calculated using monthly LAI
     * data above with max/min LAI and max/min roughness length data
     * in the vegprmt.tbl
     */
    pihm_t_struct   pihm_time;

    double          rl_tbl[40][12] = {
        /* Evergreen Needleleaf Forest */
        {0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500,
            0.500, 0.500},
        /* Evergreen Broadleaf Forest */
        {0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500,
            0.500, 0.500},
        /* Deciduous Needleleaf Forest */
        {0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500,
            0.500, 0.500},
        /* Deciduous Broadleaf Forest */
        {0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500,
            0.500, 0.500},
        /* Mixed Forest */
        {0.200, 0.200, 0.200, 0.200, 0.278, 0.367, 0.367, 0.300, 0.200, 0.200,
            0.200, 0.200},
        /* Closed Shrubland */
        {0.010, 0.010, 0.010, 0.015, 0.032, 0.048, 0.048, 0.035, 0.015, 0.010,
            0.010, 0.010},
        /* Open Shrubland */
        {0.010, 0.010, 0.010, 0.010, 0.033, 0.053, 0.053, 0.038, 0.010, 0.010,
            0.010, 0.010},
        /* Woody Savanna */
        {0.010, 0.010, 0.010, 0.016, 0.034, 0.050, 0.050, 0.038, 0.016, 0.010,
            0.010, 0.010},
        /* Savanna */
        {0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150,
            0.150, 0.150},
        /* Grassland */
        {0.100, 0.100, 0.101, 0.102, 0.106, 0.120, 0.120, 0.108, 0.102, 0.101,
            0.100, 0.100},
        /* Permanent Wetland */
        {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
            0.000, 0.000},
        /* Cropland */
        {0.050, 0.050, 0.050, 0.050, 0.050, 0.061, 0.085, 0.085, 0.050, 0.050,
            0.050, 0.050},
        /* Urban and Built-Up */
        {0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500,
            0.500, 0.500},
        /* Cropland/Natural Veg. Mosaic */
        {0.050, 0.050, 0.050, 0.050, 0.050, 0.059, 0.091, 0.050, 0.050, 0.050,
            0.050, 0.050},
        /* Permanent Snow */
        {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
            0.001, 0.001},
        /* Barren/Sparsely Vegetated */
        {0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010,
            0.010, 0.010},
        /* IGBP Water */
        {0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001,
            0.0001, 0.0001, 0.0001, 0.0001},
        /* Unclassified */
        {0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300,
            0.300, 0.300},
        /* Fill Value */
        {0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150,
            0.150, 0.150},
        /* Unclassified */
        {0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050,
            0.050, 0.050},
        /* Open Water */
        {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
            0.000, 0.000},
        /* Perennial Ice/Snow */
        {0.152, 0.151, 0.155, 0.165, 0.169, 0.169, 0.169, 0.170, 0.170, 0.166,
            0.157, 0.152},
        /* Developed Open Space */
        {0.089, 0.089, 0.091, 0.093, 0.095, 0.094, 0.093, 0.094, 0.095, 0.094,
            0.091, 0.089},
        /* Developed Low Intensity */
        {0.119, 0.119, 0.123, 0.132, 0.137, 0.136, 0.135, 0.136, 0.137, 0.133,
            0.124, 0.119},
        /* Developed Medium Intensity */
        {0.154, 0.153, 0.163, 0.179, 0.187, 0.187, 0.186, 0.187, 0.188, 0.180,
            0.163, 0.154},
        /* Developed High Intensity */
        {0.183, 0.183, 0.195, 0.218, 0.229, 0.229, 0.228, 0.229, 0.230, 0.220,
            0.196, 0.183},
        /* Barren Land */
        {0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.013, 0.014, 0.014,
            0.013, 0.013},
        /* Deciduous Forest */
        {0.343, 0.343, 0.431, 0.577, 0.650, 0.657, 0.656, 0.653, 0.653, 0.581,
            0.431, 0.343},
        /* Evergreen Forest */
        {1.623, 1.623, 1.623, 1.623, 1.623, 1.623, 1.622, 1.622, 1.623, 1.623,
            1.623, 1.623},
        /* Mixed Forest */
        {0.521, 0.518, 0.557, 0.629, 0.663, 0.664, 0.665, 0.665, 0.667, 0.633,
            0.562, 0.521},
        /* Dwarf Scrub */
        {0.022, 0.022, 0.021, 0.020, 0.020, 0.020, 0.020, 0.023, 0.025, 0.024,
            0.023, 0.022},
        /* Shrub/Scrub */
        {0.034, 0.034, 0.033, 0.033, 0.033, 0.032, 0.032, 0.034, 0.035, 0.035,
            0.034, 0.034},
        /* Grassland/Herbaceous */
        {0.070, 0.070, 0.070, 0.070, 0.070, 0.069, 0.068, 0.069, 0.070, 0.070,
            0.070, 0.070},
        /* Sedge/Herbaceous */
        {0.070, 0.070, 0.070, 0.070, 0.070, 0.069, 0.068, 0.069, 0.070, 0.070,
            0.070, 0.070},
        /* Lichens */
        {0.070, 0.070, 0.070, 0.070, 0.070, 0.069, 0.068, 0.069, 0.070, 0.070,
            0.070, 0.070},
        /* Moss */
        {0.070, 0.070, 0.070, 0.070, 0.070, 0.069, 0.068, 0.069, 0.070, 0.070,
            0.070, 0.070},
        /* Pasture/Hay */
        {0.047, 0.047, 0.047, 0.047, 0.047, 0.046, 0.046, 0.046, 0.047, 0.047,
            0.047, 0.047},
        /* Cultivated Crops */
        {0.047, 0.047, 0.047, 0.047, 0.047, 0.046, 0.046, 0.046, 0.047, 0.047,
            0.047, 0.047},
        /* Woody Wetland */
        {0.038, 0.038, 0.038, 0.037, 0.037, 0.037, 0.037, 0.039, 0.040, 0.039,
            0.039, 0.038},
        /* Emergent Herbaceous Wetland */
        {0.027, 0.027, 0.026, 0.024, 0.024, 0.024, 0.024, 0.028, 0.029, 0.029,
            0.027, 0.027}
    };

    pihm_time = PIHMTime(t);

    return rl_tbl[lc_type - 1][pihm_time.month - 1];
}

double MonthlyMf(int t)
{
    pihm_t_struct   pihm_time;

    double          mf_tbl[12] = {
        0.001308019, 0.001633298, 0.002131198, 0.002632776, 0.003031171,
        0.003197325, 0.003095839, 0.002745240, 0.002260213, 0.001759481,
        0.001373646, 0.001202083
    };

    pihm_time = PIHMTime(t);

    return mf_tbl[pihm_time.month - 1];
}
