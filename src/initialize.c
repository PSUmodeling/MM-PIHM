#include "pihm.h"

#define MAX_TYPE    100

void Initialize(pihm_struct pihm, N_Vector CV_Y, void **cvode_mem)
{
    int             i, j;
#if defined(_LUMPED_)
    int             soil_counter[MAX_TYPE];
    int             lc_counter[MAX_TYPE];

    for (i = 0; i < MAX_TYPE; i++)
    {
        soil_counter[i] = 0;
        lc_counter[i] = 0;
    }
#endif

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
#if defined(_LUMPED_)
    pihm->elem = (elem_struct *)malloc((nelem + 1) * sizeof(elem_struct));
#else
    pihm->elem = (elem_struct *)malloc(nelem * sizeof(elem_struct));
#endif
    pihm->river = (river_struct *)malloc(nriver * sizeof(river_struct));

    for (i = 0; i < nelem; i++)
    {
        pihm->elem[i].attrib.soil_type = pihm->atttbl.soil[i];
#if defined(_FBR_)
        pihm->elem[i].attrib.geol_type = pihm->atttbl.geol[i];
#endif
        pihm->elem[i].attrib.lc_type = pihm->atttbl.lc[i];
#if defined(_CYCLES_)
        pihm->elem[i].attrib.op_type = pihm->agtbl.op[i];
#endif

#if defined(_LUMPED_)
        soil_counter[pihm->elem[i].attrib.soil_type]++;
        lc_counter[pihm->elem[i].attrib.lc_type]++;
#endif
        for (j = 0; j < NUM_EDGE; j++)
        {
            pihm->elem[i].attrib.bc_type[j] = pihm->atttbl.bc[i][j];
#if defined(_FBR_)
            pihm->elem[i].attrib.fbrbc_type[j] = pihm->atttbl.fbr_bc[i][j];
#endif
        }
        pihm->elem[i].attrib.meteo_type = pihm->atttbl.meteo[i];
        pihm->elem[i].attrib.lai_type = pihm->atttbl.lai[i];
    }

#if defined(_LUMPED_)
    /* Use the soil type (land cover type) that covers the most number of model
     * grids for the lumped grid */
    pihm->elem[LUMPED].attrib.soil_type = 0;
    pihm->elem[LUMPED].attrib.lc_type = 0;
    for (i = 0; i < MAX_TYPE; i++)
    {
        pihm->elem[LUMPED].attrib.soil_type =
            (soil_counter[i] >
            soil_counter[pihm->elem[LUMPED].attrib.soil_type]) ?
            i : pihm->elem[LUMPED].attrib.soil_type;
        pihm->elem[LUMPED].attrib.lc_type =
            (lc_counter[i] >
            lc_counter[pihm->elem[LUMPED].attrib.lc_type]) ?
            i : pihm->elem[LUMPED].attrib.lc_type;
    }
#endif

    for (i = 0; i < nriver; i++)
    {
        pihm->river[i].attrib.riverbc_type = pihm->rivtbl.bc[i];
    }

    /* Initialize element mesh structures */
    InitMesh(pihm->elem, &pihm->meshtbl);

    /* Initialize element topography */
    InitTopo(pihm->elem, &pihm->meshtbl);

    /* Calculate average elevation and total area of model domain */
    pihm->siteinfo.zmax = AvgElev(pihm->elem);
    pihm->siteinfo.zmin = AvgZmin(pihm->elem);
    pihm->siteinfo.area = TotalArea(pihm->elem);
#if defined(_LUMPED_)
    pihm->elem[LUMPED].topo.zmax = pihm->siteinfo.zmax;
    pihm->elem[LUMPED].topo.zmin = pihm->siteinfo.zmin;
    pihm->elem[LUMPED].topo.area = pihm->siteinfo.area;
#endif

    /* Initialize element soil properties */
#if defined(_NOAH_)
    InitSoil(pihm->elem, &pihm->soiltbl, &pihm->noahtbl, &pihm->cal);
#else
    InitSoil(pihm->elem, &pihm->soiltbl, &pihm->cal);
#endif

#if defined(_FBR_)
    /* Initialize element geol properties */
    InitGeol(pihm->elem, &pihm->geoltbl, &pihm->cal);
#endif

    /* Initialize element land cover properties */
    InitLc(&pihm->lctbl, &pihm->cal, pihm->elem);

    /* Initialize element forcing */
    InitForc(&pihm->cal, &pihm->forc, pihm->elem);

    /* Initialize river segment properties */
    InitRiver(pihm->river, pihm->elem, &pihm->rivtbl, &pihm->shptbl,
        &pihm->matltbl, &pihm->meshtbl, &pihm->cal);

    /* Correct element elevations to avoid sinks */
    if (corr_mode)
    {
        CorrElev(pihm->elem, pihm->river);
    }

    /* Calculate distances between elements */
    InitSurfL(pihm->elem, pihm->river, &pihm->meshtbl);

#if defined(_NOAH_)
    /* Initialize land surface module (Noah) */
    InitLsm(pihm->elem, &pihm->ctrl, &pihm->noahtbl, &pihm->cal);
#endif

#if defined(_CYCLES_)
    InitCycles(&pihm->agtbl, &pihm->soiltbl, pihm->epctbl, pihm->elem,
        pihm->river);
#endif

#if defined(_BGC_)
    /* Initialize CN (Biome-BGC) module */
    InitBgc(pihm->elem, &pihm->epctbl);
#endif

    /*
     * Create hydrological and land surface initial conditions
     */
    if (pihm->ctrl.init_type == RELAX)
    {
        /* Relaxation mode */
#if defined(_NOAH_)
        /* Noah initialization needs air temperature thus forcing is applied */
        ApplyForc(pihm->ctrl.starttime, pihm->ctrl.rad_mode, &pihm->siteinfo,
            &pihm->forc, pihm->elem);
#endif

        RelaxIc(pihm->elem, pihm->river);
    }
    else if (pihm->ctrl.init_type == RST_FILE)
    {
        /* Hot start (using .ic file) */
        ReadIc(pihm->filename.ic, pihm->elem, pihm->river);
    }

    /* Initialize state variables */
    InitVar(pihm->elem, pihm->river, CV_Y);

#if defined(_CYCLES_)
# if TEMP_DISABLED
    /* Initialize Cycles module */
    if (pihm->ctrl.read_cycles_restart)
    {
        ReadCyclesIC(pihm->filename.cyclesic, pihm->elem, pihm->river);
    }
# endif
    FirstDay(&pihm->soiltbl, pihm->elem, pihm->river);

    InitCyclesVar(pihm->elem, pihm->river, CV_Y);
#endif

#if defined(_BGC_)
    /* Initialize CN variables */
    if (pihm->ctrl.read_bgc_restart)
    {
        ReadBgcIc(pihm->filename.bgcic, pihm->elem, pihm->river);
    }
    else
    {
        FirstDay(pihm->elem, pihm->river, &pihm->cninit);
    }

    InitBgcVar(pihm->elem, pihm->river, CV_Y);
#endif

    /* Calculate model time steps */
    CalcModelStep(&pihm->ctrl);

#if defined(_DAILY_)
    InitDailyStruct(pihm->elem);
#endif
}

void CorrElev(elem_struct *elem, river_struct *river)
{
    int             i, j;
    int             sink;
#if OBSOLETE
    int             bedrock_flag = 1;
#endif
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
                    river[-elem[i].nabr[j] - 1].topo.zmax;
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
                        river[-elem[i].nabr[j] - 1].topo.zmax;
                    new_elevation = (nabr_zmax < new_elevation) ?
                        nabr_zmax : new_elevation;
                    PIHMprintf(VL_NORMAL, " (%d)%lf", j + 1,
                        (elem[i].nabr[j] > 0) ?
                        elem[elem[i].nabr[j] - 1].topo.zmax :
                        river[-elem[i].nabr[j] - 1].topo.zmax);
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

#if OBSOLETE
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
                        river[-elem[i].nabr[j] - 1].topo.zmin;
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
                            river[0 - elem[i].nabr[j] - 1].topo.zmin;
                        new_elevation = (new_elevation > elem[i].topo.zmin) ?
                            elem[i].topo.zmin : new_elevation;
                        PIHMprintf (VL_NORMAL, "(%d)%lf  ", j + 1,
                            (elem[i].nabr[j] > 0) ?
                            elem[elem[i].nabr[j] - 1].topo.zmin :
                            river[0 - elem[i].nabr[j] - 1].topo.zmin);
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
        if (river[i].down > 0)
        {
            if (river[i].topo.zbed < river[river[i].down - 1].topo.zbed)
            {
                river_flag = 1;
                PIHMprintf(VL_NORMAL,
                    "River %d is lower than downstream River %d.\n",
                    i + 1, river[i].down);
            }
        }
        else
        {
            if (river[i].topo.zbed <=
                river[i].topo.node_zmax - river[i].shp.depth)
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

void InitSurfL(elem_struct *elem, const river_struct *river,
    const meshtbl_struct *meshtbl)
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;
        double          x[NUM_EDGE];
        double          y[NUM_EDGE];
        double          zmin[NUM_EDGE];
        double          zmax[NUM_EDGE];
        double          distx;
        double          disty;

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

            if (0 == elem[i].nabr[j])
            {
                elem[i].topo.nabr_x[j] = elem[i].topo.x - 2.0 * distx;
                elem[i].topo.nabr_y[j] = elem[i].topo.y - 2.0 * disty;
                elem[i].topo.nabrdist[j] =
                    sqrt(pow(elem->topo.edge[0] * elem->topo.edge[1] *
                    elem->topo.edge[2] / (4.0 * elem->topo.area), 2) -
                    pow(elem->topo.edge[j] / 2.0, 2));
            }
            else
            {
                elem[i].topo.nabr_x[j] = (elem[i].nabr[j] > 0) ?
                    elem[elem[i].nabr[j] - 1].topo.x :
                    river[0 - elem[i].nabr[j] - 1].topo.x;
                elem[i].topo.nabr_y[j] = (elem[i].nabr[j] > 0) ?
                    elem[elem[i].nabr[j] - 1].topo.y :
                    river[0 - elem[i].nabr[j] - 1].topo.y;
                elem[i].topo.nabrdist[j] =
                    (elem[i].topo.x - elem[i].topo.nabr_x[j]) *
                    (elem[i].topo.x - elem[i].topo.nabr_x[j]);
                elem[i].topo.nabrdist[j] +=
                    (elem[i].topo.y - elem[i].topo.nabr_y[j]) *
                    (elem[i].topo.y - elem[i].topo.nabr_y[j]);
                elem[i].topo.nabrdist[j] = sqrt(elem[i].topo.nabrdist[j]);
            }
        }
    }
}

double _WsAreaElev(int type, const elem_struct *elem)
{
    double          ans = 0.0;
    int             i;

    for (i = 0; i < nelem; i++)
    {
        switch (type)
        {
            case WS_ZMAX:
                ans += elem[i].topo.zmax;
                break;
            case WS_ZMIN:
                ans += elem[i].topo.zmin;
                break;
            case WS_AREA:
                ans += elem[i].topo.area;
                break;
            default:
                ans = BADVAL;
                PIHMprintf(VL_ERROR,
                    "Error: Return value type %d id not defined.\n", type);
                PIHMexit(EXIT_FAILURE);
        }
    }

    if (type == WS_AREA)
    {
        return ans;
    }
    else
    {
        return ans / (double)nelem;
    }
}

void RelaxIc(elem_struct *elem, river_struct *river)
{
    int             i;
    const double    INIT_UNSAT = 0.1;
#if defined(_FBR_)
    const double    INIT_FBR_GW = 5.0;
#endif

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem[i].ic.cmc = 0.0;
        elem[i].ic.sneqv = 0.0;
        elem[i].ic.surf = 0.0;
        elem[i].ic.unsat = INIT_UNSAT;
        elem[i].ic.gw = elem[i].soil.depth - INIT_UNSAT;

#if defined(_FBR_)
        elem[i].ic.fbr_gw = (elem[i].geol.depth > INIT_FBR_GW) ?
            INIT_FBR_GW : elem[i].geol.depth;
        elem[i].ic.fbr_unsat = 0.5 * (elem[i].geol.depth - elem[i].ic.fbr_gw);
#endif

#if defined(_NOAH_)
        int             j;
        double          sfctmp;

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

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river[i].ic.stage = 0.0;
        river[i].ic.gw = river[i].topo.zbed - river[i].topo.zmin - 0.1;
    }
}

void InitVar(elem_struct *elem, river_struct *river, N_Vector CV_Y)
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
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

#if defined(_FBR_)
        elem[i].ws.fbr_unsat = elem[i].ic.fbr_unsat;
        elem[i].ws.fbr_gw = elem[i].ic.fbr_gw;

        NV_Ith(CV_Y, FBRUNSAT(i)) = elem[i].ic.fbr_unsat;
        NV_Ith(CV_Y, FBRGW(i)) = elem[i].ic.fbr_gw;
#endif

#if defined(_NOAH_)
        int             j;

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

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        river[i].ws.stage = river[i].ic.stage;
        river[i].ws.gw = river[i].ic.gw;

        NV_Ith(CV_Y, RIVSTG(i)) = river[i].ic.stage;
        NV_Ith(CV_Y, RIVGW(i)) = river[i].ic.gw;

        river[i].ws0 = river[i].ws;
    }

    /* Other variables */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        InitWFlux(&elem[i].wf);

#if defined(_NOAH_)
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

    ctrl->tout[ctrl->nstep] = (ctrl->tout[ctrl->nstep] < ctrl->endtime) ?
        ctrl->endtime : ctrl->tout[ctrl->nstep];
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
    wf->esnow = 0.0;

#if defined(_FBR_)
    wf->fbr_infil = 0.0;
    wf->fbr_rechg = 0.0;
    for (j = 0; j < NUM_EDGE; j++)
    {
        wf->fbrflow[j] = 0.0;
    }
#endif

#if defined(_NOAH_)
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
#if defined(_CYCLES_)
    wf->eres = 0.0;
    wf->irrigationVol = 0.0;
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

#if defined(_NOAH_)
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
#if defined(_BGC_)
    ef->swabs_per_plaisun = 0.0;
    ef->swabs_per_plaishade = 0.0;
#endif
}
