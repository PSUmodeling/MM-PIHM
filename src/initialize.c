#include "pihm.h"

void Initialize(pihm_struct pihm, N_Vector CV_Y, void **cvode_mem)
{
    int             i, j;

    PIHMprintf(VL_VERBOSE, "\n\nInitialize data structure\n");

    /* Allocate memory for solver */
    *cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (*cvode_mem == NULL)
    {
        PIHMprintf(VL_ERROR, "Error in allocating memory for solver.\n");
        PIHMexit(EXIT_FAILURE);
    }

    /*
     * Initialize PIHM structure
     */
    pihm->elem = (elem_struct *)malloc(nelem * sizeof(elem_struct));
    pihm->rivseg = (river_struct *)malloc(nriver * sizeof(river_struct));

    for (i = 0; i < nelem; i++)
    {
        pihm->elem[i].attrib.soil_type = pihm->atttbl.soil[i];
        pihm->elem[i].attrib.lc_type = pihm->atttbl.lc[i];
        for (j = 0; j < NUM_EDGE; j++)
        {
            pihm->elem[i].attrib.bc_type[j] = pihm->atttbl.bc[i][j];
        }
        pihm->elem[i].attrib.meteo_type = pihm->atttbl.meteo[i];
        pihm->elem[i].attrib.lai_type = pihm->atttbl.lai[i];
    }

    for (i = 0; i < nriver; i++)
    {
        pihm->rivseg[i].attrib.riverbc_type = pihm->rivtbl.bc[i];
    }

    /* Initialize element mesh structures */
    InitMesh(pihm->elem, &pihm->meshtbl);

    /* Initialize element topography */
    InitTopo(pihm->elem, &pihm->meshtbl);

#ifdef _NOAH_
    /* Calculate average elevation and total area of model domain */
    pihm->siteinfo.elevation = AvgElev(pihm->elem);
    pihm->siteinfo.area = TotalArea(pihm->elem);
#endif

    /* Initialize element soil properties */
#ifdef _NOAH_
    InitSoil(pihm->elem, &pihm->soiltbl, &pihm->noahtbl, &pihm->cal);
#else
    InitSoil(pihm->elem, &pihm->soiltbl, &pihm->cal);
#endif

    /* Initialize element land cover properties */
    InitLc(pihm->elem, &pihm->lctbl, &pihm->cal);

    /* Initialize element forcing */
#ifdef _BGC_
    InitForc(pihm->elem, &pihm->forc, &pihm->cal, pihm->co2.varco2,
        pihm->ndepctrl.varndep);
#else
    InitForc(pihm->elem, &pihm->forc, &pihm->cal);
#endif

    /* Initialize river segment properties */
    InitRiver(pihm->rivseg, pihm->elem, &pihm->rivtbl, &pihm->shptbl,
        &pihm->matltbl, &pihm->meshtbl, &pihm->cal);

    /* Correct element elevations to avoid sinks */
    if (corr_mode)
    {
        CorrElev(pihm->elem, pihm->rivseg);
    }

    /* Calculate distances between elements */
    InitSurfL(pihm->elem, pihm->rivseg, &pihm->meshtbl);

#ifdef _NOAH_
    /* Initialize land surface module (Noah) */
    InitLsm(pihm->elem, &pihm->ctrl, &pihm->noahtbl, &pihm->cal);
#endif

#ifdef _CYCLES_
    /* Initialize Cycles modeule */
    if (pihm->ctrl.read_cycles_restart)
    {
        ReadCyclesIC(pihm->filename.cyclesic, pihm->elem, pihm->rivseg);
    }

    InitCycles(pihm->elem, pihm->rivseg, &pihm->ctrl, &pihm->mgmttbl,
        &pihm->agtbl, &pihm->croptbl, &pihm->soiltbl);
#endif

#ifdef _BGC_
    /* Initialize CN (Biome-BGC) module */
    InitBgc(pihm->elem, &pihm->epctbl);
#endif

    /*
     * Create hydrological and land surface initial conditions
     */
    if (pihm->ctrl.init_type == RELAX)
    {
        /* Relaxation mode */
#ifdef _NOAH_
        /* Noah initialization needs air temperature thus forcing is applied */
        ApplyForc(&pihm->forc, pihm->elem, pihm->ctrl.starttime, &pihm->ctrl,
            &pihm->siteinfo);
#endif
        RelaxIc(pihm->elem, pihm->rivseg);
    }
    else if (pihm->ctrl.init_type == RST_FILE)
    {
        /* Hot start (using .ic file) */
        ReadIc(pihm->filename.ic, pihm->elem, pihm->rivseg);
    }

    /* Initialize state variables */
    InitVar(pihm->elem, pihm->rivseg, CV_Y);

#ifdef _BGC_
    /* Initialize CN variables */
    if (pihm->ctrl.read_bgc_restart)
    {
        ReadBgcIc(pihm->filename.bgcic, pihm->elem, pihm->rivseg);
    }
    else
    {
        FirstDay(pihm->elem, pihm->rivseg, &pihm->cninit);
    }

    InitBgcVar(pihm->elem, pihm->rivseg, CV_Y);
#endif

    /* Calculate model time steps */
    CalcModelStep(&pihm->ctrl);

#ifdef _DAILY_
    InitDailyStruct(pihm->elem, pihm->rivseg);
#endif
}

void InitMesh(elem_struct *elem, const meshtbl_struct *meshtbl)
{
    int             i, j;

    for (i = 0; i < nelem; i++)
    {
        elem[i].ind = i + 1;

        for (j = 0; j < NUM_EDGE; j++)
        {
            elem[i].node[j] = meshtbl->node[i][j];
            elem[i].nabr[j] = meshtbl->nabr[i][j];
        }
    }
}

void InitTopo(elem_struct *elem, const meshtbl_struct *meshtbl)
{
    int             i, j;
    double          x[NUM_EDGE];
    double          y[NUM_EDGE];
    double          zmin[NUM_EDGE];
    double          zmax[NUM_EDGE];

    for (i = 0; i < nelem; i++)
    {
        for (j = 0; j < NUM_EDGE; j++)
        {
            x[j] = meshtbl->x[elem[i].node[j] - 1];
            y[j] = meshtbl->y[elem[i].node[j] - 1];
            zmin[j] = meshtbl->zmin[elem[i].node[j] - 1];
            zmax[j] = meshtbl->zmax[elem[i].node[j] - 1];
        }

        elem[i].topo.area = 0.5 *
            ((x[1] - x[0]) * (y[2] - y[0]) - (y[1] - y[0]) * (x[2] - x[0]));
        /* Calculate centroid of triangle */
        elem[i].topo.x = (x[0] + x[1] + x[2]) / 3.0;
        elem[i].topo.y = (y[0] + y[1] + y[2]) / 3.0;

        elem[i].topo.zmin = (zmin[0] + zmin[1] + zmin[2]) / 3.0;
        elem[i].topo.zmax = (zmax[0] + zmax[1] + zmax[2]) / 3.0;
        elem[i].topo.edge[0] =
            sqrt(pow((x[1] - x[2]), 2) + pow((y[1] - y[2]), 2));
        elem[i].topo.edge[1] =
            sqrt(pow((x[2] - x[0]), 2) + pow((y[2] - y[0]), 2));
        elem[i].topo.edge[2] =
            sqrt(pow((x[0] - x[1]), 2) + pow((y[0] - y[1]), 2));
    }

#ifdef _NOAH_
    CalcSlopeAspect(elem, meshtbl);
#endif
}

#ifdef _NOAH_
void InitSoil(elem_struct *elem, const soiltbl_struct *soiltbl,
    const noahtbl_struct *noahtbl, const calib_struct *cal)
#else
void InitSoil(elem_struct *elem, const soiltbl_struct *soiltbl,
    const calib_struct *cal)
#endif
{
    int             i;
    int             soil_ind;

    for (i = 0; i < nelem; i++)
    {
        soil_ind = elem[i].attrib.soil_type - 1;

        elem[i].soil.dinf = cal->dinf * soiltbl->dinf;

        elem[i].soil.depth = elem[i].topo.zmax - elem[i].topo.zmin;

        elem[i].soil.ksath = cal->ksath * soiltbl->ksath[soil_ind];
        elem[i].soil.ksatv = cal->ksatv * soiltbl->ksatv[soil_ind];
        elem[i].soil.kinfv = cal->kinfv * soiltbl->kinfv[soil_ind];

        elem[i].soil.smcmin = cal->porosity * soiltbl->smcmin[soil_ind];
        elem[i].soil.smcmax = cal->porosity * soiltbl->smcmax[soil_ind];
        elem[i].soil.porosity = elem[i].soil.smcmax - elem[i].soil.smcmin;
        if (elem[i].soil.porosity > 1.0 || elem[i].soil.porosity <= 0.0)
        {
            PIHMprintf(VL_ERROR,
                "Error: Porosity value out of bounds for Element %d", i + 1);
            PIHMexit(EXIT_FAILURE);
        }
        elem[i].soil.alpha = cal->alpha * soiltbl->alpha[soil_ind];
        elem[i].soil.beta = cal->beta * soiltbl->beta[soil_ind];

        /* Calculate field capacity and wilting point following Chan and
         * Dudhia 2001 MWR, but replacing Campbell with van Genuchten */
        elem[i].soil.smcwlt = cal->porosity * soiltbl->smcwlt[soil_ind];
#ifdef _NOAH_
        elem[i].soil.smcwlt *= cal->smcwlt;
#endif
        elem[i].soil.smcref = cal->porosity * soiltbl->smcref[soil_ind];
#ifdef _NOAH_
        elem[i].soil.smcref *= cal->smcref;
#endif

        elem[i].soil.dmac = cal->dmac * soiltbl->dmac[soil_ind];
        elem[i].soil.dmac = (elem[i].soil.dmac > elem[i].soil.depth) ?
            elem[i].soil.depth : elem[i].soil.dmac;

        elem[i].soil.areafh = cal->areafh * soiltbl->areafh[soil_ind];
        elem[i].soil.areafv = cal->areafv * soiltbl->areafv[soil_ind];

        elem[i].soil.kmacv =
            cal->kmacv * soiltbl->kmacv_ro * soiltbl->kinfv[soil_ind];
        elem[i].soil.kmach =
            cal->kmach * soiltbl->kmach_ro * soiltbl->ksath[soil_ind];
#ifdef _NOAH_
        elem[i].soil.csoil = noahtbl->csoil;
        elem[i].soil.quartz = soiltbl->qtz[soil_ind];
        elem[i].soil.smcdry = elem[i].soil.smcwlt;
#endif
    }
}

void InitLc(elem_struct *elem, const lctbl_struct *lctbl,
    const calib_struct *cal)
{
    int             i;
    int             lc_ind;

    for (i = 0; i < nelem; i++)
    {
        lc_ind = elem[i].attrib.lc_type - 1;

        elem[i].ps.rzd = cal->rzd * lctbl->rzd[lc_ind];
        elem[i].epc.rsmin = lctbl->rsmin[lc_ind];
        elem[i].epc.rgl = lctbl->rgl[lc_ind];
        elem[i].epc.hs = lctbl->hs[lc_ind];
        elem[i].epc.rsmax = lctbl->rsmax;
        elem[i].epc.topt = lctbl->topt;
        elem[i].lc.shdfac = cal->vegfrac * lctbl->vegfrac[lc_ind];
        elem[i].lc.laimin = lctbl->laimin[lc_ind];
        elem[i].lc.laimax = lctbl->laimax[lc_ind];
        elem[i].lc.emissmin = lctbl->emissmin[lc_ind];
        elem[i].lc.emissmax = lctbl->emissmax[lc_ind];
        elem[i].lc.albedomin = cal->albedo * lctbl->albedomin[lc_ind];
        elem[i].lc.albedomax = cal->albedo * lctbl->albedomax[lc_ind];
        elem[i].lc.z0min = lctbl->z0min[lc_ind];
        elem[i].lc.z0max = lctbl->z0max[lc_ind];
        elem[i].lc.rough = cal->rough * lctbl->rough[lc_ind];
        elem[i].lc.cmcfactr = CMCFACTR;
        elem[i].lc.cfactr = lctbl->cfactr;
        elem[i].lc.bare = (elem[i].attrib.lc_type == lctbl->bare) ? 1 : 0;
        elem[i].lc.shdfac = (elem[i].lc.bare == 1) ? 0.0 : elem[i].lc.shdfac;
#ifdef _NOAH_
        elem[i].lc.snup = lctbl->snup[lc_ind];
        elem[i].lc.isurban = (elem[i].attrib.lc_type == ISURBAN) ? 1 : 0;
        elem[i].lc.shdmin = 0.01;
        elem[i].lc.shdmax = 0.96;
#endif

#ifdef _NOAH_
        elem[i].epc.rgl *= cal->rgl;
        elem[i].epc.hs *= cal->hs;
        elem[i].epc.rsmin *= cal->rsmin;
        elem[i].lc.cmcfactr *= cal->cmcmax;
        elem[i].lc.cfactr *= cal->cfactr;
#endif
    }
}

void InitRiver(river_struct *rivseg, elem_struct *elem,
    const rivtbl_struct *rivtbl, const shptbl_struct *shptbl,
    const matltbl_struct *matltbl, const meshtbl_struct *meshtbl,
    const calib_struct *cal)
{
    int             i, ii, j;

    for (i = 0; i < nriver; i++)
    {
        rivseg[i].leftele = rivtbl->leftele[i];
        rivseg[i].rightele = rivtbl->rightele[i];
        rivseg[i].fromnode = rivtbl->fromnode[i];
        rivseg[i].tonode = rivtbl->tonode[i];
        rivseg[i].down = rivtbl->down[i];
        rivseg[i].up = 0;
        for (ii = 0; ii < nriver; ii++)
        {
            if (rivtbl->down[ii] == i + 1)
            {
                rivseg[i].up = ii + 1;
            }
        }

        for (j = 0; j < NUM_EDGE; j++)
        {
            /* Note: Strategy to use BC < -4 for river identification */
            if (elem[rivseg[i].leftele - 1].nabr[j] == rivseg[i].rightele)
            {
                elem[rivseg[i].leftele - 1].nabr[j] = -(i + 1);
            }
            if (elem[rivseg[i].rightele - 1].nabr[j] == rivseg[i].leftele)
            {
                elem[rivseg[i].rightele - 1].nabr[j] = -(i + 1);
            }
        }

        rivseg[i].topo.x = 0.5 *
            (meshtbl->x[rivseg[i].fromnode - 1] +
            meshtbl->x[rivseg[i].tonode - 1]);
        rivseg[i].topo.y = 0.5 *
            (meshtbl->y[rivseg[i].fromnode - 1] +
            meshtbl->y[rivseg[i].tonode - 1]);
        rivseg[i].topo.zmax = 0.5 *
            (meshtbl->zmax[rivseg[i].fromnode - 1] +
            meshtbl->zmax[rivseg[i].tonode - 1]);
        rivseg[i].topo.zmin = rivseg[i].topo.zmax -
            (0.5 * (elem[rivseg[i].leftele - 1].topo.zmax +
            elem[rivseg[i].rightele - 1].topo.zmax) -
            0.5 * (elem[rivseg[i].leftele - 1].topo.zmin +
            elem[rivseg[i].rightele - 1].topo.zmin));
        rivseg[i].topo.node_zmax = meshtbl->zmax[rivseg[i].tonode - 1];
        rivseg[i].topo.dist_left = sqrt(
            (rivseg[i].topo.x - elem[rivseg[i].leftele - 1].topo.x) *
            (rivseg[i].topo.x - elem[rivseg[i].leftele - 1].topo.x) +
            (rivseg[i].topo.y - elem[rivseg[i].leftele - 1].topo.y) *
            (rivseg[i].topo.y - elem[rivseg[i].leftele - 1].topo.y));
        rivseg[i].topo.dist_right = sqrt(
            (rivseg[i].topo.x - elem[rivseg[i].rightele - 1].topo.x) *
            (rivseg[i].topo.x - elem[rivseg[i].rightele - 1].topo.x) +
            (rivseg[i].topo.y - elem[rivseg[i].rightele - 1].topo.y) *
            (rivseg[i].topo.y - elem[rivseg[i].rightele - 1].topo.y));

        rivseg[i].shp.depth =
            cal->rivdepth * shptbl->depth[rivtbl->shp[i] - 1];
        rivseg[i].shp.intrpl_ord = shptbl->intrpl_ord[rivtbl->shp[i] - 1];
        rivseg[i].shp.coeff =
            cal->rivshpcoeff * shptbl->coeff[rivtbl->shp[i] - 1];
        rivseg[i].shp.length = sqrt(
            pow(meshtbl->x[rivseg[i].fromnode - 1] -
            meshtbl->x[rivseg[i].tonode - 1], 2) +
            pow(meshtbl->y[rivseg[i].fromnode - 1] -
            meshtbl->y[rivseg[i].tonode - 1], 2));
        rivseg[i].shp.width = RiverEqWid(rivseg->shp.intrpl_ord,
            rivseg->shp.depth, rivseg->shp.coeff);

        rivseg[i].topo.zbed = rivseg[i].topo.zmax - rivseg[i].shp.depth;

        rivseg[i].matl.rough =
            cal->rivrough * matltbl->rough[rivtbl->matl[i] - 1];
        rivseg[i].matl.cwr = matltbl->cwr[rivtbl->matl[i] - 1];
        rivseg[i].matl.ksath =
            cal->rivksath * matltbl->ksath[rivtbl->matl[i] - 1];
        rivseg[i].matl.ksatv =
            cal->rivksatv * matltbl->ksatv[rivtbl->matl[i] - 1];
        rivseg[i].matl.bedthick =
            cal->rivbedthick * matltbl->bedthick[rivtbl->matl[i] - 1];
        rivseg[i].matl.porosity = 0.5 *
            (elem[rivseg[i].leftele - 1].soil.porosity +
            elem[rivseg[i].rightele - 1].soil.porosity);
        rivseg[i].matl.smcmin = 0.5 *
            (elem[rivseg[i].leftele - 1].soil.smcmin +
            elem[rivseg[i].rightele - 1].soil.smcmin);

        rivseg[i].topo.area = rivseg[i].shp.length *
            RiverEqWid(rivseg[i].shp.intrpl_ord, rivseg[i].shp.depth,
            rivseg[i].shp.coeff);
    }
}

#ifdef _BGC_
void InitForc(elem_struct *elem, forc_struct *forc, const calib_struct *cal,
    int varco2, int varndep)
#else
void InitForc(elem_struct *elem, forc_struct *forc, const calib_struct *cal)
#endif
{
    int             i, j;

    /* Apply scenarios */
    for (i = 0; i < forc->nmeteo; i++)
    {
        for (j = 0; j < forc->meteo[i].length; j++)
        {
            forc->meteo[i].data[j][PRCP_TS] *= cal->prcp;
            forc->meteo[i].data[j][SFCTMP_TS] += cal->sfctmp;
        }
    }

    if (forc->nbc > 0)
    {
        for (i = 0; i < forc->nbc; i++)
        {
            forc->bc[i].value = (double *)malloc(sizeof(double));
        }
    }
    if (forc->nmeteo > 0)
    {
        for (i = 0; i < forc->nmeteo; i++)
        {
            forc->meteo[i].value =
                (double *)malloc(NUM_METEO_VAR * sizeof(double));
        }
    }
    if (forc->nlai > 0)
    {
        for (i = 0; i < forc->nlai; i++)
        {
            forc->lai[i].value = (double *)malloc(sizeof(double));
        }
    }
    if (forc->nriverbc > 0)
    {
        for (i = 0; i < forc->nriverbc; i++)
        {
            forc->riverbc[i].value = (double *)malloc(sizeof(double));
        }
    }
    if (forc->nsource > 0)
    {
        for (i = 0; i < forc->nsource; i++)
        {
            forc->source[i].value = (double *)malloc(sizeof(double));
        }
    }
#ifdef _NOAH_
    if (forc->nrad > 0)
    {
        for (i = 0; i < forc->nrad; i++)
        {
            forc->rad[i].value = (double *)malloc(2 * sizeof(double));
        }
    }
#endif

#ifdef _BGC_
    if (!spinup_mode && varco2)
    {
        forc->co2[0].value = (double *)malloc(sizeof(double));
    }

    if (!spinup_mode && varndep)
    {
        forc->ndep[0].value = (double *)malloc(sizeof(double));
    }
#endif

    for (i = 0; i < nelem; i++)
    {
        elem[i].ps.zlvl_wind =
            forc->meteo[elem[i].attrib.meteo_type - 1].zlvl_wind;
    }
}

void CorrElev(elem_struct *elem, river_struct *rivseg)
{
    int             i, j;
    int             sink;
    int             bedrock_flag = 1;
    int             river_flag = 0;
    double          nabr_zmax;
    double          new_elevation;

    PIHMprintf(VL_VERBOSE, "Correct surface elevation.\n");

    for (i = 0; i < nelem; i++)
    {
        /* Correction of surface elevation (artifacts due to coarse scale
         * discretization). Not needed if there is lake feature. */
        sink = 1;

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].nabr[j] != 0)
            {
                nabr_zmax = (elem[i].nabr[j] > 0) ?
                    elem[elem[i].nabr[j] - 1].topo.zmax :
                    rivseg[0 - elem[i].nabr[j] - 1].topo.zmax;
                if (elem[i].topo.zmax >= nabr_zmax)
                {
                    sink = 0;
                    break;
                }
            }
        }

        if (sink == 1)
        {
            PIHMprintf(VL_NORMAL, "Element %d is a sink\n", i + 1);

            /* Note: Following correction is being applied for correction
             * mode only */
            PIHMprintf(VL_NORMAL, "\tBefore: surface %lf, "
                "bedrock %lf. Neighbors surface:",
                elem[i].topo.zmax, elem[i].topo.zmin);

            new_elevation = 1.0e7;
            for (j = 0; j < NUM_EDGE; j++)
            {
                if (elem[i].nabr[j] != 0)
                {
                    nabr_zmax = (elem[i].nabr[j] > 0) ?
                        elem[elem[i].nabr[j] - 1].topo.zmax :
                        rivseg[0 - elem[i].nabr[j] - 1].topo.zmax;
                    new_elevation = (nabr_zmax < new_elevation) ?
                        nabr_zmax : new_elevation;
                    PIHMprintf(VL_NORMAL, " (%d)%lf", j + 1,
                        (elem[i].nabr[j] > 0) ?
                        elem[elem[i].nabr[j] - 1].topo.zmax :
                        rivseg[0 - elem[i].nabr[j] - 1].topo.zmax);
                }
            }

            /* Lift bedrock elevation by the same amount */
            elem[i].topo.zmin += (new_elevation - elem[i].topo.zmax);

            /* Apply new surface elevation */
            elem[i].topo.zmax = new_elevation;

            PIHMprintf(VL_NORMAL, ". Corrected = %lf, %lf\n",
                elem[i].topo.zmax, elem[i].topo.zmin);
        }
    }

#ifdef OBSOLETE
    /* Correction of BedRck Elev. Is this needed? */
    if (bedrock_flag == 1)
    {
        for (i = 0; i < numele; i++)
        {
            sink = 1;
            for (j = 0; j < NUM_EDGE; j++)
            {
                if (elem[i].nabr[j] != 0)
                {
                    new_elevation = (elem[i].nabr[j] > 0) ?
                        elem[elem[i].nabr[j] - 1].topo.zmin :
                        rivseg[-elem[i].nabr[j] - 1].topo.zmin;
                    if (elem[i].topo.zmin - new_elevation >= 0.0)
                    {
                        sink = 0;
                        break;
                    }
                }
            }
            if (sink == 1)
            {
                PIHMprintf (VL_NORMAL, "Ele %d (bedrock) is sink", i + 1);
                /* Note: Following correction is being applied for debug==1
                 * case only */
                PIHMprintf (VL_NORMAL, "\tBefore: %lf Corrected using:",
                    elem[i].topo.zmin);
                new_elevation = 1.0e7;
                for (j = 0; j < NUM_EDGE; j++)
                {
                    if (elem[i].nabr[j] != 0)
                    {
                        elem[i].topo.zmin = (elem[i].nabr[j] > 0) ?
                            elem[elem[i].nabr[j] - 1].topo.zmin :
                            rivseg[0 - elem[i].nabr[j] - 1].topo.zmin;
                        new_elevation = (new_elevation > elem[i].topo.zmin) ?
                            elem[i].topo.zmin : new_elevation;
                        PIHMprintf (VL_NORMAL, "(%d)%lf  ", j + 1,
                            (elem[i].nabr[j] > 0) ?
                            elem[elem[i].nabr[j] - 1].topo.zmin :
                            rivseg[0 - elem[i].nabr[j] - 1].topo.zmin);
                    }
                }
                elem[i].topo.zmin = new_elevation;
                PIHMprintf (VL_NORMAL, "=(New)%lf\n", elem[i].topo.zmin);
            }
        }
    }
#endif

    for (i = 0; i < nriver; i++)
    {
        if (rivseg[i].down > 0)
        {
            if (rivseg[i].topo.zbed < rivseg[rivseg[i].down - 1].topo.zbed)
            {
                river_flag = 1;
                PIHMprintf(VL_NORMAL,
                    "River %d is lower than downstream River %d.\n",
                    i + 1, rivseg[i].down);
            }
        }
        else
        {
            if (rivseg[i].topo.zbed <=
                rivseg[i].topo.node_zmax - rivseg[i].shp.depth)
            {
                river_flag = 1;
                PIHMprintf(VL_NORMAL,
                    "River outlet is higher than the channel (River %d).\n",
                    i + 1);
            }
        }
    }

    if (river_flag == 1)
    {
        PIHMprintf(VL_NORMAL, "\nRiver elevation correction needed "
            "but PIHM will continue without fixing river elevation.\n\n");
    }

    sleep(5);
}

void InitSurfL(elem_struct *elem, river_struct *rivseg,
    const meshtbl_struct *meshtbl)
{
    int             i, j;
    double          x[NUM_EDGE];
    double          y[NUM_EDGE];
    double          zmin[NUM_EDGE];
    double          zmax[NUM_EDGE];
    double          distx;
    double          disty;

    for (i = 0; i < nelem; i++)
    {
        for (j = 0; j < NUM_EDGE; j++)
        {
            x[j] = meshtbl->x[elem[i].node[j] - 1];
            y[j] = meshtbl->y[elem[i].node[j] - 1];
            zmin[j] = meshtbl->zmin[elem[i].node[j] - 1];
            zmax[j] = meshtbl->zmax[elem[i].node[j] - 1];
        }

        for (j = 0; j < NUM_EDGE; j++)
        {
            /* Note: Assumption here is that the formulation is circumcenter
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
                elem[i].topo.nabrdist_x[j] = elem[i].topo.x - 2.0 * distx;
                elem[i].topo.nabrdist_y[j] = elem[i].topo.y - 2.0 * disty;
                elem[i].topo.nabrdist[j] =
                    sqrt(pow(elem->topo.edge[0] * elem->topo.edge[1] *
                    elem->topo.edge[2] / (4.0 * elem->topo.area), 2) -
                    pow(elem->topo.edge[j] / 2.0, 2));
            }
            else
            {
                elem[i].topo.nabrdist_x[j] = (elem[i].nabr[j] > 0) ?
                    elem[elem[i].nabr[j] - 1].topo.x :
                    rivseg[0 - elem[i].nabr[j] - 1].topo.x;
                elem[i].topo.nabrdist_y[j] = (elem[i].nabr[j] > 0) ?
                    elem[elem[i].nabr[j] - 1].topo.y :
                    rivseg[0 - elem[i].nabr[j] - 1].topo.y;
                elem[i].topo.nabrdist[j] =
                    (elem[i].topo.x - elem[i].topo.nabrdist_x[j]) *
                    (elem[i].topo.x - elem[i].topo.nabrdist_x[j]);
                elem[i].topo.nabrdist[j] +=
                    (elem[i].topo.y - elem[i].topo.nabrdist_y[j]) *
                    (elem[i].topo.y - elem[i].topo.nabrdist_y[j]);
                elem[i].topo.nabrdist[j] = sqrt(elem[i].topo.nabrdist[j]);
            }
        }
    }
}

void RelaxIc(elem_struct *elem, river_struct *rivseg)
{
    int             i;
#ifdef _NOAH_
    int             j;
    double          sfctmp;
#endif

    for (i = 0; i < nelem; i++)
    {
        elem[i].ic.cmc = 0.0;
        elem[i].ic.sneqv = 0.0;
        elem[i].ic.surf = 0.0;
        elem[i].ic.unsat = 0.1;
        elem[i].ic.gw = elem[i].soil.depth - 0.1;

#ifdef _NOAH_
        sfctmp = elem[i].es.sfctmp;

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

    for (i = 0; i < nriver; i++)
    {
        rivseg[i].ic.stage = 0.0;
        rivseg[i].ic.gw = rivseg[i].topo.zbed - rivseg[i].topo.zmin - 0.1;
    }
}

void InitVar(elem_struct *elem, river_struct *rivseg, N_Vector CV_Y)
{
    int             i;
#ifdef _NOAH_
    int             j;
#endif

    /* State variables (initial conditions) */
    for (i = 0; i < nelem; i++)
    {
        elem[i].ws.cmc = elem[i].ic.cmc;
        elem[i].ws.sneqv = elem[i].ic.sneqv;

        elem[i].ws.surf = elem[i].ic.surf;
        elem[i].ws.unsat = elem[i].ic.unsat;
        elem[i].ws.gw = elem[i].ic.gw;

        NV_Ith(CV_Y, SURF(i)) = elem[i].ic.surf;
        NV_Ith(CV_Y, UNSAT(i)) = elem[i].ic.unsat;
        NV_Ith(CV_Y, GW(i)) = elem[i].ic.gw;

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

    for (i = 0; i < nriver; i++)
    {
        rivseg[i].ws.stage = rivseg[i].ic.stage;
        rivseg[i].ws.gw = rivseg[i].ic.gw;

        NV_Ith(CV_Y, RIVSTG(i)) = rivseg[i].ic.stage;
        NV_Ith(CV_Y, RIVGW(i)) = rivseg[i].ic.gw;

        rivseg[i].ws0 = rivseg[i].ws;
    }

    /* Other variables */
    for (i = 0; i < nelem; i++)
    {
        InitWFlux(&elem[i].wf);

#ifdef _NOAH_
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
        elem[i].ef.esnow = BADVAL;

        elem[i].wf.runoff2 = BADVAL;
        elem[i].wf.runoff3 = BADVAL;
        elem[i].wf.pcpdrp = 0.0;
        elem[i].wf.drip = 0.0;
        elem[i].wf.dew = BADVAL;
        elem[i].wf.snomlt = BADVAL;
        elem[i].wf.esnow = BADVAL;

        elem[i].ws.soilm = BADVAL;
#endif
    }
}

void CalcModelStep(ctrl_struct *ctrl)
{
    int             i;

    ctrl->nstep = (ctrl->endtime - ctrl->starttime) / ctrl->stepsize;

    ctrl->tout = (int *)malloc((ctrl->nstep + 1) * sizeof(int));

    for (i = 0; i < ctrl->nstep + 1; i++)
    {
        ctrl->tout[i] = (i == 0) ?
            ctrl->starttime : ctrl->tout[i - 1] + ctrl->stepsize;
    }

    if (ctrl->tout[ctrl->nstep] < ctrl->endtime)
    {
        ctrl->tout[ctrl->nstep] = ctrl->endtime;
    }
}

void InitWFlux(wflux_struct *wf)
{
    int             j;

    for (j = 0; j < NUM_EDGE; j++)
    {
        wf->ovlflow[j] = 0.0;
        wf->subsurf[j] = 0.0;
    }
    wf->prcp = 0.0;
    wf->pcpdrp = 0.0;
    wf->infil = 0.0;
    wf->rechg = 0.0;
    wf->drip = 0.0;
    wf->edir = 0.0;
    wf->ett = 0.0;
    wf->ec = 0.0;
    wf->etp = 0.0;
    wf->eta = 0.0;
    wf->edir_surf = 0.0;
    wf->edir_unsat = 0.0;
    wf->edir_gw = 0.0;
    wf->ett_unsat = 0.0;
    wf->ett_gw = 0.0;

#ifdef _NOAH_
    int             k;

    for (k = 0; k < MAXLYR; k++)
    {
        wf->et[k] = 0.0;
    }
    wf->runoff2 = 0.0;
    for (k = 0; k < MAXLYR; k++)
    {
        wf->runoff2_lyr[k] = 0.0;
    }
    wf->runoff3 = 0.0;
    for (k = 0; k < MAXLYR; k++)
    {
        wf->smflxv[k] = 0.0;
    }
    for (j = 0; j < NUM_EDGE; j++)
    {
        for (k = 0; k < MAXLYR; k++)
        {
            wf->smflxh[j][k] = 0.0;
        }
    }
    wf->dew = 0.0;
    wf->snomlt = 0.0;
    wf->esnow = 0.0;
    wf->etns = 0.0;
#endif
#ifdef _CYCLES_
    wf->eres = 0.0;
#endif
}

void InitRiverWFlux(river_wflux_struct *wf)
{
    int             j;

    for (j = 0; j < NUM_RIVFLX; j++)
    {
        wf->rivflow[j] = 0.0;
    }
}

void InitEFlux(eflux_struct *ef)
{
    ef->soldn = 0.0;

#ifdef _NOAH_
    int             k;

    ef->solnet = 0.0;
    ef->etp = 0.0;
    ef->ssoil = 0.0;
    ef->eta = 0.0;
    ef->sheat = 0.0;
    ef->fdown = 0.0;
    ef->lwdn = 0.0;
    ef->ec = 0.0;
    ef->edir = 0.0;
    for (k = 0; k < MAXLYR; k++)
    {
        ef->et[k] = 0.0;
    }
    ef->ett = 0.0;
    ef->esnow = 0.0;
    ef->soldir = 0.0;
    ef->soldif = 0.0;
    ef->longwave = 0.0;
    ef->flx1 = 0.0;
    ef->flx2 = 0.0;
    ef->flx3 = 0.0;
#endif
#ifdef _BGC_
    ef->swabs_per_plaisun = 0.0;
    ef->swabs_per_plaishade = 0.0;
#endif
}
