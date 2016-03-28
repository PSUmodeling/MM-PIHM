#include "pihm.h"

void Initialize (pihm_struct pihm, N_Vector CV_Y)
{
    if (verbose_mode)
    {
        printf ("\n\nInitialize data structure\n");
    }

    pihm->elem = (elem_struct *)malloc (pihm->numele * sizeof (elem_struct));
    pihm->riv = (river_struct *)malloc (pihm->numriv * sizeof (river_struct));

    InitMeshStruct (pihm->elem, pihm->numele, pihm->meshtbl);

    InitTopo (pihm->elem, pihm->numele, pihm->meshtbl);

#ifdef _NOAH_
    /* Calculate average elevation of model domain */
    pihm->elevation = AvgElev (pihm->elem, pihm->numele);
#endif

    InitSoil (pihm->elem, pihm->numele, pihm->atttbl, pihm->soiltbl,
#ifdef _NOAH_
        pihm->noahtbl,
#endif
        pihm->cal);

    InitLC (pihm->elem, pihm->numele, pihm->atttbl, pihm->lctbl, pihm->cal);

    InitForcing (pihm->elem, pihm->numele, pihm->riv, pihm->numriv,
        pihm->atttbl, pihm->rivtbl, &pihm->forc, pihm->cal);

    InitRiver (pihm->riv, pihm->numriv, pihm->elem, pihm->rivtbl,
        pihm->shptbl, pihm->matltbl, pihm->meshtbl, pihm->cal);

    if (debug_mode)
    {
        CorrectElevation (pihm->elem, pihm->numele, pihm->riv, pihm->numriv);
    }

    InitSurfL (pihm->elem, pihm->numele, pihm->riv, pihm->meshtbl);

#ifdef _NOAH_
    InitLsm (pihm->elem, pihm->numele, pihm->ctrl, pihm->noahtbl, pihm->cal);
#endif

#ifdef _CYCLES_
    InitCycles (pihm->elem, pihm->numele, &pihm->mgmttbl, &pihm->agtbl);
#endif

    if (pihm->ctrl.init_type == RELAX)
    {
#ifdef _NOAH_
        ApplyForcing (&pihm->forc, pihm->ctrl.starttime);
#endif
        SaturationIC (pihm->elem, pihm->numele, pihm->riv, pihm->numriv);
    }
    else if (pihm->ctrl.init_type == RST_FILE)
    {
        ReadIC (pihm->filename.ic, pihm->elem, pihm->numele, pihm->riv,
            pihm->numriv);
    }

    InitVar (pihm->elem, pihm->numele, pihm->riv, pihm->numriv, CV_Y);

    CalcModelStep (&pihm->ctrl);

//#ifdef _DAILY_
//    InitDailyStruct (pihm);
//#endif
}

void InitMeshStruct (elem_struct *elem, int numele, meshtbl_struct meshtbl)
{
    int             i, j;

    for (i = 0; i < numele; i++)
    {
        elem[i].ind = i + 1;

        for (j = 0; j < 3; j++)
        {
            elem[i].node[j] = meshtbl.node[i][j];
            elem[i].nabr[j] = meshtbl.nabr[i][j];
        }
    }
}

void InitTopo (elem_struct *elem, int numele, meshtbl_struct meshtbl)
{
    int             i, j;
    double          x[3];
    double          y[3];
    double          zmin[3];
    double          zmax[3];

    for (i = 0; i < numele; i++)
    {
        for (j = 0; j < 3; j++)
        {
            x[j] = meshtbl.x[elem[i].node[j] - 1];
            y[j] = meshtbl.y[elem[i].node[j] - 1];
            zmin[j] = meshtbl.zmin[elem[i].node[j] - 1];
            zmax[j] = meshtbl.zmax[elem[i].node[j] - 1];
        }

        elem[i].topo.area =
            0.5 * ((x[1] - x[0]) * (y[2] - y[0]) - (y[1] - y[0]) * (x[2] -
                x[0]));
        /* Calculate centroid of triangle */
        elem[i].topo.x = (x[0] + x[1] + x[2]) / 3.0;
        elem[i].topo.y = (y[0] + y[1] + y[2]) / 3.0;

        elem[i].topo.zmin = (zmin[0] + zmin[1] + zmin[2]) / 3.0;
        elem[i].topo.zmax = (zmax[0] + zmax[1] + zmax[2]) / 3.0;
        elem[i].topo.edge[0] =
            sqrt (pow ((x[1] - x[2]), 2) + pow ((y[1] - y[2]), 2));
        elem[i].topo.edge[1] =
            sqrt (pow ((x[2] - x[0]), 2) + pow ((y[2] - y[0]), 2));
        elem[i].topo.edge[2] =
            sqrt (pow ((x[0] - x[1]), 2) + pow ((y[0] - y[1]), 2));
    }

#ifdef _NOAH_
    CalcSlopeAspect (elem, numele, meshtbl);
#endif
}

void InitSoil (elem_struct *elem, int numele, atttbl_struct atttbl,
    soiltbl_struct soiltbl,
#ifdef _NOAH_
    noahtbl_struct noahtbl,
#endif
    calib_struct cal)
{
    int             i;
    int             soil_ind;

    for (i = 0; i < numele; i++)
    {
        soil_ind = atttbl.soil[i] - 1;

        elem[i].soil.dinf = cal.dinf * soiltbl.dinf;

        elem[i].soil.depth = elem[i].topo.zmax - elem[i].topo.zmin;

        elem[i].soil.ksath = cal.ksath * soiltbl.ksath[soil_ind];
        elem[i].soil.ksatv = cal.ksatv * soiltbl.ksatv[soil_ind];
        elem[i].soil.kinfv = cal.kinfv * soiltbl.kinfv[soil_ind];

        elem[i].soil.smcmin = cal.porosity * soiltbl.smcmin[soil_ind];
        elem[i].soil.smcmax = cal.porosity * soiltbl.smcmax[soil_ind];
        elem[i].soil.porosity = elem[i].soil.smcmax - elem[i].soil.smcmin;
        if (elem[i].soil.porosity > 1.0 || elem[i].soil.porosity <= 0.0)
        {
            printf ("Warning: Porosity value out of bounds for element %d",
                i + 1);
            PihmExit (1);
        }
        elem[i].soil.alpha = cal.alpha * soiltbl.alpha[soil_ind];
        elem[i].soil.beta = cal.beta * soiltbl.beta[soil_ind];


        /* Calculate field capacity and wilting point following Chan and
         * Dudhia 2001 MWR, but replacing Campbell with van Genuchten */
        elem[i].soil.smcwlt = cal.porosity * soiltbl.smcwlt[soil_ind];
#ifdef _NOAH_
        elem[i].soil.smcwlt *= cal.thetaw;
#endif
        elem[i].soil.smcref = cal.porosity * soiltbl.smcref[soil_ind];
#ifdef _NOAH_
        elem[i].soil.smcref *= cal.thetaref;
#endif

        elem[i].soil.dmac = cal.dmac * soiltbl.dmac[soil_ind];
        elem[i].soil.dmac = (elem[i].soil.dmac > elem[i].soil.depth) ?
            elem[i].soil.depth : elem[i].soil.dmac;

        elem[i].soil.areafh = cal.areafh * soiltbl.areafh[soil_ind];
        elem[i].soil.areafv = cal.areafv * soiltbl.areafv[soil_ind];
        elem[i].soil.macropore = atttbl.macropore[i];

        elem[i].soil.kmacv =
            cal.kmacv * soiltbl.kmacv_ro * soiltbl.kinfv[soil_ind];
        elem[i].soil.kmach =
            cal.kmach * soiltbl.kmach_ro * soiltbl.ksath[soil_ind];
#ifdef _NOAH_
        elem[i].soil.csoil = noahtbl.csoil;
        elem[i].soil.quartz = soiltbl.qtz[soil_ind];
        elem[i].soil.smcdry = elem[i].soil.smcwlt;
#endif
    }
}

double FieldCapacity (double alpha, double beta, double kv, double smcmax,
    double smcmin)
{
    /*
     * Solve field capacity using Newton's method
     * Field capacity is defined as (Chen and Dudhia 2001 MWR)
     * Theta_ref = Theta_s * (1 / 3 + 2 / 3 * satn_ref)
     * where satn_ref is the soil saturation ratio when the soil hydraulic
     * conductivity reaches 5.70E-9 m/s
     */
    double          satn;
    double          satnk;
    double          denom;
    double          df;
    double          dsatn;
    double          mx;
    double          ftheta;
    int             n;
    double          smcref;
    const double    KFC = 5.79E-9;
    const double    ERROR = 1.0E-3;

    mx = 1.0 - 1.0 / beta;

    n = 0;

    satn = 0.75;

    while (n < 10 && dsatn > ERROR)
    {
        n++;

        df = KrFunc (alpha, beta, satn) - KFC / kv;

        ftheta = 1.0 - pow (satn, 1.0 / mx);

        denom = 0.5 * pow (satn, -0.5) * pow (1.0 - pow (ftheta, mx), 2.0) +
            2.0 * pow (satn, 1.0 / mx - 0.5) * (pow (ftheta,
                mx - 1.0) - pow (ftheta, 2.0 * mx - 1.0));

        satnk = satn - df / denom;
        satnk = (satnk > 1.0 - 1.0E-3) ? 1.0 - 1.0E-3 : satnk;
        satnk = (satnk < SATMIN) ? SATMIN : satnk;

        dsatn = fabs (satnk - satn);

        satn = satnk;
    }

    if (dsatn > ERROR)
    {
        satn = 0.75;
    }

    smcref = (1.0 / 3.0 + 2.0 / 3.0 * satn) * (smcmax - smcmin) + smcmin;

    return (smcref);
}

double WiltingPoint (double smcmax, double smcmin, double alpha, double beta)
{
    const double    PSIW = 200.0;

    return (0.5 * (smcmax - smcmin) * pow (1.0 / (1.0 + pow (PSIW * alpha,
                    beta)), 1.0 - 1.0 / beta) + smcmin);
}

void InitLC (elem_struct *elem, int numele, atttbl_struct atttbl,
    lctbl_struct lctbl, calib_struct cal)
{
    int             i;
    int             lc_ind;

    for (i = 0; i < numele; i++)
    {
        lc_ind = atttbl.lc[i] - 1;
        elem[i].lc.type = lc_ind + 1;
        elem[i].lc.shdfac = cal.vegfrac * lctbl.vegfrac[lc_ind];
        elem[i].lc.rzd = cal.rzd * lctbl.rzd[lc_ind];
        elem[i].lc.rsmin = lctbl.rsmin[lc_ind];
        elem[i].lc.rgl = lctbl.rgl[lc_ind];
        elem[i].lc.hs = lctbl.hs[lc_ind];
        elem[i].lc.laimin = lctbl.laimin[lc_ind];
        elem[i].lc.laimax = lctbl.laimax[lc_ind];
        elem[i].lc.emissmin = lctbl.emissmin[lc_ind];
        elem[i].lc.emissmax = lctbl.emissmax[lc_ind];
        elem[i].lc.albedomin = cal.albedo * lctbl.albedomin[lc_ind];
        elem[i].lc.albedomax = cal.albedo * lctbl.albedomax[lc_ind];
        elem[i].lc.z0min = lctbl.z0min[lc_ind];
        elem[i].lc.z0max = lctbl.z0max[lc_ind];
        elem[i].lc.rough = cal.rough * lctbl.rough[lc_ind];
        elem[i].lc.cmcfactr = 0.0002;
        elem[i].lc.rsmax = lctbl.rsmax;
        elem[i].lc.cfactr = lctbl.cfactr;
        elem[i].lc.topt = lctbl.topt;
        elem[i].lc.bare = (elem[i].lc.type == lctbl.bare) ? 1 : 0;
        elem[i].lc.shdfac = (elem[i].lc.bare == 1) ? 0.0 : elem[i].lc.shdfac;
#ifdef _NOAH_
        elem[i].lc.ptu = 0.0;   /* Not yet used */
        elem[i].lc.snup = lctbl.snup[lc_ind];
        elem[i].lc.isurban = (elem[i].lc.type == ISURBAN) ? 1 : 0;
        elem[i].lc.shdmin = 0.01;
        elem[i].lc.shdmax = 0.96;
#endif

#ifdef _NOAH_
        elem[i].lc.rsmin *= cal.rsmin;
        elem[i].lc.cmcfactr *= cal.intcp;
        elem[i].lc.cfactr *= cal.cfactr;
        elem[i].lc.rgl *= cal.rgl;
        elem[i].lc.hs *= cal.hs;
#endif
    }
}

void InitRiver (river_struct *riv, int numriv, elem_struct *elem,
    rivtbl_struct rivtbl, shptbl_struct shptbl,
    matltbl_struct matltbl, meshtbl_struct meshtbl, calib_struct cal)
{
    int             i, j;

    for (i = 0; i < numriv; i++)
    {
        riv[i].leftele = rivtbl.leftele[i];
        riv[i].rightele = rivtbl.rightele[i];
        riv[i].fromnode = rivtbl.fromnode[i];
        riv[i].tonode = rivtbl.tonode[i];
        riv[i].down = rivtbl.down[i];

        for (j = 0; j < 3; j++)
        {
            /* Note: Strategy to use BC < -4 for river identification */
            if (elem[riv[i].leftele - 1].nabr[j] == riv[i].rightele)
            {
                elem[riv[i].leftele - 1].nabr[j] = 0 - (i + 1);
            }
            if (elem[riv[i].rightele - 1].nabr[j] == riv[i].leftele)
            {
                elem[riv[i].rightele - 1].nabr[j] = 0 - (i + 1);
            }
        }

        riv[i].topo.x =
            (meshtbl.x[riv[i].fromnode - 1] + meshtbl.x[riv[i].tonode -
                1]) / 2.0;
        riv[i].topo.y =
            (meshtbl.y[riv[i].fromnode - 1] + meshtbl.y[riv[i].tonode -
                1]) / 2.0;
        riv[i].topo.zmax =
            (meshtbl.zmax[riv[i].fromnode - 1] +
            meshtbl.zmax[riv[i].tonode - 1]) / 2.0;
        riv[i].topo.zmin =
            riv[i].topo.zmax - (0.5 * (elem[riv[i].leftele - 1].topo.zmax +
                elem[riv[i].rightele - 1].topo.zmax) -
            0.5 * (elem[riv[i].leftele - 1].topo.zmin + elem[riv[i].rightele -
                    1].topo.zmin));
        riv[i].topo.node_zmax = meshtbl.zmax[riv[i].tonode - 1];

        riv[i].shp.depth = cal.rivdepth * shptbl.depth[rivtbl.shp[i] - 1];
        riv[i].shp.intrpl_ord = shptbl.intrpl_ord[rivtbl.shp[i] - 1];
        riv[i].shp.coeff = cal.rivshpcoeff * shptbl.coeff[rivtbl.shp[i] - 1];
        riv[i].shp.length =
            sqrt (pow (meshtbl.x[riv[i].fromnode - 1] -
                meshtbl.x[riv[i].tonode - 1],
                2) + pow (meshtbl.y[riv[i].fromnode - 1] -
                meshtbl.y[riv[i].tonode - 1], 2));

        riv[i].topo.zbed = riv[i].topo.zmax - riv[i].shp.depth;

        riv[i].matl.rough = cal.rivrough * matltbl.rough[rivtbl.matl[i] - 1];
        riv[i].matl.cwr = matltbl.cwr[rivtbl.matl[i] - 1];
        riv[i].matl.ksath = cal.rivksath * matltbl.ksath[rivtbl.matl[i] - 1];
        riv[i].matl.ksatv = cal.rivksatv * matltbl.ksatv[rivtbl.matl[i] - 1];
        riv[i].matl.bedthick =
            cal.rivbedthick * matltbl.bedthick[rivtbl.matl[i] - 1];
        riv[i].matl.porosity =
            0.5 * (elem[riv[i].leftele - 1].soil.porosity +
            elem[riv[i].rightele - 1].soil.porosity);

        riv[i].topo.area = riv[i].shp.length *
            EqWid (riv[i].shp.intrpl_ord, riv[i].shp.depth, riv[i].shp.coeff);
    }
}

void InitForcing (elem_struct *elem, int numele, river_struct *riv,
    int numriv, atttbl_struct atttbl, rivtbl_struct rivtbl,
    forc_struct *forc, calib_struct cal)
{
    int             i, j;

    for (i = 0; i < forc->nmeteo; i++)
    {
        for (j = 0; j < forc->meteo[i].length; j++)
        {
            forc->meteo[i].data[j][PRCP_TS] *= cal.prcp;
            forc->meteo[i].data[j][SFCTMP_TS] *= cal.sfctmp;
        }
    }

    if (forc->nbc > 0)
    {
        for (i = 0; i < forc->nbc; i++)
        {
            forc->bc[i].value = (double *)malloc (sizeof (double));
        }
    }
    if (forc->nmeteo > 0)
    {
        for (i = 0; i < forc->nmeteo; i++)
        {
            forc->meteo[i].value =
                (double *)malloc (NUM_METEO_VAR * sizeof (double));
        }
    }
    if (forc->nlai > 0)
    {
        for (i = 0; i < forc->nlai; i++)
        {
            forc->lai[i].value = (double *)malloc (sizeof (double));
        }
    }
    if (forc->nriverbc > 0)
    {
        for (i = 0; i < forc->nriverbc; i++)
        {
            forc->riverbc[i].value = (double *)malloc (sizeof (double));
        }
    }
    if (forc->nsource > 0)
    {
        for (i = 0; i < forc->nsource; i++)
        {
            forc->source[i].value = (double *)malloc (sizeof (double));
        }
    }
#ifdef _NOAH_
    if (forc->nrad > 0)
    {
        for (i = 0; i < forc->nrad; i++)
        {
            forc->rad[i].value = (double *)malloc (2 * sizeof (double));
        }
    }
#endif

    for (i = 0; i < numele; i++)
    {
        for (j = 0; j < 3; j++)
        {
            elem[i].forc.bc_type[j] = atttbl.bc[i][j];
            if (elem[i].forc.bc_type[j] > 0)
            {
                elem[i].forc.bc[j] =
                    forc->bc[elem[i].forc.bc_type[j] - 1].value;
            }
            else if (elem[i].forc.bc_type[j] < 0)
            {
                elem[i].forc.bc[j] =
                    forc->bc[0 - elem[i].forc.bc_type[j] - 1].value;
            }
        }

        for (j = 0; j < NUM_METEO_VAR; j++)
        {
            elem[i].forc.meteo[j] =
                &forc->meteo[atttbl.meteo[i] - 1].value[j];
        }
        elem[i].ps.zlvl_wind = forc->meteo[atttbl.meteo[i] - 1].zlvl_wind;

        elem[i].forc.lai_type = atttbl.lai[i];
        if (elem[i].forc.lai_type > 0)
        {
            elem[i].forc.lai = forc->lai[atttbl.lai[i] - 1].value;
        }

        if (atttbl.source[i] > 0)
        {
            elem[i].forc.source = forc->source[atttbl.source[i] - 1].value;
        }

#ifdef _NOAH_
        if (forc->nrad > 0)
        {
            for (j = 0; j < 2; j++)
            {
                elem[i].forc.rad[j] =
                    &forc->rad[atttbl.meteo[i] - 1].value[j];
            }
        }
#endif
    }

    for (i = 0; i < numriv; i++)
    {
        if (rivtbl.bc[i] > 0)
        {
            riv[i].forc.riverbc = forc->riverbc[rivtbl.bc[i] - 1].value;
        }
    }
}

void CorrectElevation (elem_struct *elem, int numele, river_struct *riv,
    int numriv)
{
    int             i, j;
    int             sink;
    int             bedrock_flag;
    int             river_flag;
    double          new_elevation;

    for (i = 0; i < numele; i++)
    {
        /* Correction of Surf Elev (artifacts due to coarse scale
         * discretization). Not needed if there is lake feature. */
        sink = 1;

        for (j = 0; j < 3; j++)
        {
            if (elem[i].nabr[j] != 0)
            {
                new_elevation = (elem[i].nabr[j] > 0) ?
                    elem[elem[i].nabr[j] - 1].topo.zmax :
                    riv[0 - elem[i].nabr[j] - 1].topo.zmax;
                if (elem[i].topo.zmax - new_elevation >= 0)
                {
                    sink = 0;
                    break;
                }
            }
        }

        if (sink == 1)
        {
            printf ("Ele %d (surface) is sink", i + 1);
            /* Note: Following correction is being applied for debug==1
             * case only */
            printf ("\tBefore: %lf Corrected using: ", elem[i].topo.zmax);
            new_elevation = 1.0e7;
            for (j = 0; j < 3; j++)
            {
                if (elem[i].nabr[j] != 0)
                {
                    elem[i].topo.zmax = (elem[i].nabr[j] > 0) ?
                        elem[elem[i].nabr[j] - 1].topo.zmax :
                        riv[0 - elem[i].nabr[j] - 1].topo.zmax;
                    new_elevation = (new_elevation > elem[i].topo.zmax) ?
                        elem[i].topo.zmax : new_elevation;
                    printf ("(%d)%lf  ", j + 1,
                        (elem[i].nabr[j] > 0) ?
                        elem[elem[i].nabr[j] - 1].topo.zmax :
                        riv[0 - elem[i].nabr[j] - 1].topo.zmax);
                }
            }
            elem[i].topo.zmax = new_elevation;
            printf ("=(New)%lf\n", elem[i].topo.zmax);
        }
    }
    /* Correction of BedRck Elev. Is this needed? */

    bedrock_flag = 1;
    if (bedrock_flag == 1)
    {
        for (i = 0; i < numele; i++)
        {
            sink = 1;
            for (j = 0; j < 3; j++)
            {
                if (elem[i].nabr[j] != 0)
                {
                    new_elevation = (elem[i].nabr[j] > 0) ?
                        elem[elem[i].nabr[j] - 1].topo.zmin :
                        riv[-elem[i].nabr[j] - 1].topo.zmin;
                    if (elem[i].topo.zmin - new_elevation >= 0.0)
                    {
                        sink = 0;
                        break;
                    }
                }
            }
            if (sink == 1)
            {
                printf ("Ele %d (bedrock) is sink", i + 1);
                /* Note: Following correction is being applied for debug==1
                 * case only */
                printf ("\tBefore: %lf Corrected using:", elem[i].topo.zmin);
                new_elevation = 1.0e7;
                for (j = 0; j < 3; j++)
                {
                    if (elem[i].nabr[j] != 0)
                    {
                        elem[i].topo.zmin = (elem[i].nabr[j] > 0) ?
                            elem[elem[i].nabr[j] - 1].topo.zmin :
                            riv[0 - elem[i].nabr[j] - 1].topo.zmin;
                        new_elevation = (new_elevation > elem[i].topo.zmin) ?
                            elem[i].topo.zmin : new_elevation;
                        printf ("(%d)%lf  ", j + 1,
                            (elem[i].nabr[j] > 0) ?
                            elem[elem[i].nabr[j] - 1].topo.zmin :
                            riv[0 - elem[i].nabr[j] - 1].topo.zmin);
                    }
                }
                elem[i].topo.zmin = new_elevation;
                printf ("=(New)%lf\n", elem[i].topo.zmin);
            }
        }
    }

    for (i = 0; i < numriv; i++)
    {
        if (riv[i].down > 0)
        {
            if (riv[i].topo.zbed < riv[riv[i].down - 1].topo.zbed)
            {
                river_flag = 1;
                printf ("Riv %d is lower than downstream Riv %d\n",
                    i + 1, riv[i].down);
            }
        }
    }
    if (river_flag == 1)
    {
        printf ("\n\tRiver elevation correction needed\n");
        getchar ();
    }
}

void InitSurfL (elem_struct *elem, int numele, river_struct *riv,
    meshtbl_struct meshtbl)
{
    int             i, j;
    double          x[3];
    double          y[3];
    double          zmin[3];
    double          zmax[3];
    double          distx;
    double          disty;

    for (i = 0; i < numele; i++)
    {
        for (j = 0; j < 3; j++)
        {
            x[j] = meshtbl.x[elem[i].node[j] - 1];
            y[j] = meshtbl.y[elem[i].node[j] - 1];
            zmin[j] = meshtbl.zmin[elem[i].node[j] - 1];
            zmax[j] = meshtbl.zmax[elem[i].node[j] - 1];
        }

        for (j = 0; j < 3; j++)
        {
            /* Note: Assumption here is that the forumulation is circumcenter
             * based */
            switch (j)
            {
                case 0:
                    distx = (elem[i].topo.x - 0.5 * (x[1] + x[2]));
                    disty = (elem[i].topo.y - 0.5 * (y[1] + y[2]));
                    break;
                case 1:
                    distx = (elem[i].topo.x - 0.5 * (x[2] + x[0]));
                    disty = (elem[i].topo.y - 0.5 * (y[2] + y[0]));
                    break;
                case 2:
                    distx = (elem[i].topo.x - 0.5 * (x[0] + x[1]));
                    disty = (elem[i].topo.y - 0.5 * (y[0] + y[1]));
                    break;
            }

            if (elem[i].nabr[j] == 0)
            {
                elem[i].topo.surfx[j] = elem[i].topo.x - 2.0 * distx;
                elem[i].topo.surfy[j] = elem[i].topo.y - 2.0 * disty;
            }
            else
            {
                elem[i].topo.surfx[j] = (elem[i].nabr[j] > 0) ?
                    elem[elem[i].nabr[j] - 1].topo.x :
                    riv[0 - elem[i].nabr[j] - 1].topo.x;
                elem[i].topo.surfy[j] = (elem[i].nabr[j] > 0) ?
                    elem[elem[i].nabr[j] - 1].topo.y :
                    riv[0 - elem[i].nabr[j] - 1].topo.y;
            }
        }
    }
}

void SaturationIC (elem_struct *elem, int numele, river_struct *riv,
    int numriv)
{
    int             i;
#ifdef _NOAH_
    int             j;
    double          sfctmp;
#endif

    for (i = 0; i < numele; i++)
    {
        elem[i].ic.intcp = 0.0;
        elem[i].ic.sneqv = 0.0;
        elem[i].ic.surf = 0.0;
        elem[i].ic.unsat = 0.1;
        elem[i].ic.gw = elem[i].soil.depth - 0.1;

#ifdef _NOAH_
        sfctmp = *elem[i].forc.meteo[SFCTMP_TS];

        elem[i].ic.t1 = sfctmp;

        elem[i].ic.stc[0] = sfctmp +
            (sfctmp - elem[i].ps.tbot) / elem[i].ps.zbot *
            elem[i].ps.sldpth[0] * 0.5;

        for (j = 1; j < MAXLYR; j++)
        {
            if (elem[i].ps.sldpth[j] > 0.0)
            {
                elem[i].ic.stc[j] =
                    elem[i].ic.stc[j - 1] +
                    (sfctmp - elem[i].ps.tbot) / elem[i].ps.zbot *
                    (elem[i].ps.sldpth[j - 1] + elem[i].ps.sldpth[j]) * 0.5;
            }
            else
            {
                elem[i].ic.stc[j] = BADVAL;
            }
        }

        for (j = 0; j < MAXLYR; j++)
        {
            if (elem[i].ps.sldpth[j] > 0.0)
            {
                elem[i].ic.smc[j] = elem[i].soil.smcmax;
                elem[i].ic.sh2o[j] = elem[i].soil.smcmax;
            }
            else
            {
                elem[i].ic.smc[j] = BADVAL;
                elem[i].ic.sh2o[j] = BADVAL;
            }
        }
        elem[i].ic.snowh = 0.0;
#endif
    }

    for (i = 0; i < numriv; i++)
    {
        riv[i].ic.stage = 0.0;
        riv[i].ic.gw = riv[i].topo.zbed - riv[i].topo.zmin - 0.1;
    }
}

void InitVar (elem_struct *elem, int numele, river_struct *riv,
    int numriv, N_Vector CV_Y)
{
    int             i;
#ifdef _NOAH_
    int             j;
#endif

    /* State variables (initial conditions) */
    for (i = 0; i < numele; i++)
    {
        elem[i].ws.cmc = elem[i].ic.intcp;
        elem[i].ws.sneqv = elem[i].ic.sneqv;

        elem[i].ws.surf = elem[i].ic.surf;
        elem[i].ws.unsat = elem[i].ic.unsat;
        elem[i].ws.gw = elem[i].ic.gw;

        NV_Ith_S (CV_Y, i) = elem[i].ic.surf;
        NV_Ith_S (CV_Y, i + numele) = elem[i].ic.unsat;
        NV_Ith_S (CV_Y, i + 2 * numele) = elem[i].ic.gw;
#ifdef _NOAH_
        elem[i].es.t1 = elem[i].ic.t1;
        elem[i].ps.snowh = elem[i].ic.snowh;

        for (j = 0; j < MAXLYR; j++)
        {
            elem[i].es.stc[j] = elem[i].ic.stc[j];
            elem[i].ws.smc[j] = elem[i].ic.smc[j];
            elem[i].ws.sh2o[j] = elem[i].ic.sh2o[j];
        }
#endif

        elem[i].ws0 = elem[i].ws;
    }

    for (i = 0; i < numriv; i++)
    {
        riv[i].ws.stage = riv[i].ic.stage;
        riv[i].ws.gw = riv[i].ic.gw;

        NV_Ith_S (CV_Y, i + 3 * numele) = riv[i].ic.stage;
        NV_Ith_S (CV_Y, i + 3 * numele + numriv) = riv[i].ic.gw;

        riv[i].ws0 = riv[i].ws;
    }

    /* Other variables */
    for (i = 0; i < numele; i++)
    {
        ZeroWaterFlux (&elem[i].wf);

#ifdef _NOAH_
        ZeroWaterFlux (&elem[i].avgwf);

        elem[i].ps.snotime1 = 0.0;
        elem[i].ps.ribb = 0.0;
        elem[i].ps.fcr = 1.0;
        elem[i].ps.snoalb = 0.75;
        elem[i].ps.zlvl = 3.0;
        elem[i].ps.emissi = 0.96;
        elem[i].ps.albedo = 0.18;
        elem[i].ps.z0 = 0.1;
        elem[i].ps.ch = 1.0e-4;
        elem[i].ps.cm = 1.0e-4;
        elem[i].ps.beta = BADVAL;
        elem[i].ps.sncovr = BADVAL;
        elem[i].ps.rc = BADVAL;
        elem[i].ps.pc = BADVAL;
        elem[i].ps.rcs = BADVAL;
        elem[i].ps.rct = BADVAL;
        elem[i].ps.rcsoil = BADVAL;
        elem[i].ps.q1 = BADVAL;
        elem[i].ps.cosz = BADVAL;
        elem[i].ps.z0brd = BADVAL;
        elem[i].ps.eta_kinematic = BADVAL;

        elem[i].ef.sheat = BADVAL;
        elem[i].ef.eta = BADVAL;
        elem[i].ef.fdown = BADVAL;
        elem[i].ef.ec = BADVAL;
        elem[i].ef.edir = BADVAL;
        elem[i].ef.ett = BADVAL;
        elem[i].ef.etp = BADVAL;
        elem[i].ef.ssoil = BADVAL;
        elem[i].ef.flx1 = BADVAL;
        elem[i].ef.flx2 = BADVAL;
        elem[i].ef.flx3 = BADVAL;
        elem[i].ef.solardirect = BADVAL;
        elem[i].ef.esnow = BADVAL;

        elem[i].wf.runoff1 = BADVAL;
        elem[i].wf.runoff2 = BADVAL;
        elem[i].wf.runoff3 = BADVAL;
        elem[i].wf.pcpdrp = 0.0;
        elem[i].wf.drip = 0.0;
        elem[i].wf.prcprain = BADVAL;
        elem[i].wf.dew = BADVAL;
        elem[i].wf.snomlt = BADVAL;
        elem[i].wf.esnow = BADVAL;

        elem[i].ws.soilw = BADVAL;
        elem[i].ws.soilm = BADVAL;
#endif
    }
}

void CalcModelStep (ctrl_struct *ctrl)
{
    int             i;

    ctrl->nstep = (ctrl->endtime - ctrl->starttime) / ctrl->stepsize;

    ctrl->tout = (int *)malloc ((ctrl->nstep + 1) * sizeof (int));

    for (i = 0; i < ctrl->nstep + 1; i++)
    {
        if (i == 0)
        {
            ctrl->tout[i] = ctrl->starttime;
        }
        else
        {
            ctrl->tout[i] = ctrl->tout[i - 1] + ctrl->stepsize;
        }
    }

    if (ctrl->tout[ctrl->nstep] < ctrl->endtime)
    {
        ctrl->tout[ctrl->nstep] = ctrl->endtime;
    }
}

void ZeroWaterFlux (wf_struct *wf)
{
    int             j;

#ifdef _NOAH_
    wf->runoff2 = 0.0;
#endif
    wf->infil = 0.0;
    wf->prcp = 0.0;
    wf->netprcp = 0.0;
    wf->rechg = 0.0;
    wf->drip = 0.0;
    wf->edir = 0.0;
    wf->ett = 0.0;
    wf->ec = 0.0;
    wf->etp = 0.0;
    wf->eta = 0.0;
    wf->edir_sfc = 0.0;
    wf->edir_unsat = 0.0;
    wf->edir_gw = 0.0;
    wf->ett_unsat = 0.0;
    wf->ett_gw = 0.0;

    for (j = 0; j < 3; j++)
    {
        wf->fluxsurf[j] = 0.0;
        wf->fluxsub[j] = 0.0;
    }

    for (j = 0; j < MAXLYR; j++)
    {
        wf->et[j] = 0.0;
    }

    for (j = 0; j < 11; j++)
    {
        wf->fluxriv[j] = 0.0;
    }
}
