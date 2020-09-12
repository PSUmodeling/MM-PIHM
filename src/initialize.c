#include "pihm.h"

#define MAX_TYPE    100

void Initialize(pihm_struct pihm, N_Vector CV_Y, void **cvode_mem)
{
    int             i, j;
    int             bc;
#if defined(_LUMPEDBGC_)
    int             soil_counter[MAX_TYPE];
    int             lc_counter[MAX_TYPE];

    for (i = 0; i < MAX_TYPE; i++)
    {
        soil_counter[i] = 0;
        lc_counter[i] = 0;
    }
#endif

    pihm_printf(VL_VERBOSE, "\n\nInitialize data structure\n");

    /* Allocate memory for solver */
    *cvode_mem = CVodeCreate(CV_BDF);
    if (*cvode_mem == NULL)
    {
        pihm_printf(VL_ERROR, "Error in allocating memory for solver.\n");
        pihm_exit(EXIT_FAILURE);
    }

    /*
     * Initialize PIHM structure
     */
#if defined(_LUMPEDBGC_)
    pihm->elem = (elem_struct *)malloc((nelem + 1) * sizeof(elem_struct));
#else
    pihm->elem = (elem_struct *)malloc(nelem * sizeof(elem_struct));
#endif
    pihm->river = (river_struct *)malloc(nriver * sizeof(river_struct));

    for (i = 0; i < nelem; i++)
    {
        pihm->elem[i].attrib.soil = pihm->atttbl.soil[i];
#if defined(_DGW_)
        pihm->elem[i].attrib.geol = pihm->atttbl.geol[i];
#endif
        pihm->elem[i].attrib.lc = pihm->atttbl.lc[i];

#if defined(_LUMPEDBGC_)
        soil_counter[pihm->elem[i].attrib.soil]++;
        lc_counter[pihm->elem[i].attrib.lc]++;
#endif

        for (j = 0; j < NUM_EDGE; j++)
        {
            bc = pihm->atttbl.bc[i][j];

            if (bc == NO_FLOW)
            {
                pihm->elem[i].attrib.bc[j] = NO_FLOW;
            }
            else
            {
                /* Adjust bc_type flag so that positive values indicate
                 * Dirichlet type, and negative values indicate Neumann type */
                pihm->elem[i].attrib.bc[j] =
                    (pihm->forc.bc[bc - 1].bc_type == DIRICHLET) ? bc : -bc;
            }

#if defined(_DGW_)
            bc = pihm->atttbl.bc_geol[i][j];

            if (bc == NO_FLOW)
            {
                pihm->elem[i].attrib.bc_geol[j] = NO_FLOW;
            }
            else
            {
                /* Adjust bc_type flag so that positive values indicate
                 * Dirichlet type, and negative values indicate Neumann type */
                pihm->elem[i].attrib.bc_geol[j] =
                    (pihm->forc.bc[bc - 1].bc_type == DIRICHLET) ? bc : -bc;
            }
#endif
        }

        pihm->elem[i].attrib.meteo = pihm->atttbl.meteo[i];
        pihm->elem[i].attrib.lai = pihm->atttbl.lai[i];

#if defined(_RT_)
        pihm->elem[i].attrib.prcp_conc = pihm->atttbl.prcpc[i];
        for (j = 0; j < NCHMVOL; j++)
        {
            pihm->elem[i].attrib.chem_ic[j] = pihm->atttbl.chem_ic[i][j];
        }
#endif
    }

#if defined(_LUMPEDBGC_)
    /* Use the soil type (land cover type) that covers the most number of model
     * grids for the lumped grid */
    pihm->elem[LUMPEDBGC].attrib.soil = 0;
    pihm->elem[LUMPEDBGC].attrib.lc = 0;
    for (i = 0; i < MAX_TYPE; i++)
    {
        pihm->elem[LUMPEDBGC].attrib.soil =
            (soil_counter[i] >
            soil_counter[pihm->elem[LUMPEDBGC].attrib.soil]) ?
            i : pihm->elem[LUMPEDBGC].attrib.soil;
        pihm->elem[LUMPEDBGC].attrib.lc =
            (lc_counter[i] >
            lc_counter[pihm->elem[LUMPEDBGC].attrib.lc]) ?
            i : pihm->elem[LUMPEDBGC].attrib.lc;
    }
#endif

    for (i = 0; i < nriver; i++)
    {
        bc = pihm->rivtbl.bc[i];

        if (bc == NO_FLOW)
        {
            pihm->river[i].attrib.riverbc_type = NO_FLOW;
        }
        else
        {
            /* Adjust bc_type flag so that positive values indicate
             * Dirichlet type, and negative values indicate Neumann type */
            pihm->river[i].attrib.riverbc_type =
                (pihm->forc.riverbc[bc - 1].bc_type == DIRICHLET) ? bc : -bc;
        }
    }

    /* Initialize element mesh structures */
    InitMesh(&pihm->meshtbl, pihm->elem);

    /* Initialize element topography */
    InitTopo(&pihm->meshtbl, pihm->elem);

    /* Calculate average elevation and total area of model domain */
    pihm->siteinfo.zmax = AvgElev(pihm->elem);
    pihm->siteinfo.zmin = AvgZmin(pihm->elem);
    pihm->siteinfo.area = TotalArea(pihm->elem);
#if defined(_LUMPEDBGC_)
    pihm->elem[LUMPEDBGC].topo.zmax = pihm->siteinfo.zmax;
    pihm->elem[LUMPEDBGC].topo.zmin = pihm->siteinfo.zmin;
    pihm->elem[LUMPEDBGC].topo.area = pihm->siteinfo.area;
#endif

    /* Initialize element soil properties */
#if defined(_NOAH_)
    InitSoil(&pihm->soiltbl, &pihm->noahtbl, &pihm->calib, pihm->elem);
#else
    InitSoil(&pihm->soiltbl, &pihm->calib, pihm->elem);
#endif

#if defined(_DGW_)
    /* Initialize element geol properties */
    InitGeol(&pihm->geoltbl, &pihm->calib, pihm->elem);
#endif

    /* Initialize element land cover properties */
    InitLc(&pihm->lctbl, &pihm->calib, pihm->elem);

    /* Initialize element forcing */
#if defined(_RT_)
    InitForcing(&pihm->rttbl, &pihm->calib, &pihm->forc, pihm->elem);
#else
    InitForcing(&pihm->calib, &pihm->forc, pihm->elem);
#endif

    /* Initialize river segment properties */
    InitRiver(&pihm->meshtbl, &pihm->rivtbl, &pihm->shptbl, &pihm->matltbl,
        &pihm->calib, pihm->elem, pihm->river);

    /* Correct element elevations to avoid sinks */
    if (corr_mode)
    {
        CorrectElev(pihm->river, pihm->elem);
    }

    /* Calculate distances between elements */
    InitSurfL(&pihm->meshtbl, pihm->elem);

#if defined(_NOAH_)
    /* Initialize land surface module (Noah) */
    InitLsm(pihm->filename.ice, &pihm->ctrl, &pihm->noahtbl, &pihm->calib,
        pihm->elem);
#endif

#if defined(_CYCLES_)
    InitCycles(&pihm->calib, &pihm->agtbl, pihm->mgmttbl, pihm->croptbl,
        &pihm->soiltbl, pihm->elem);
#endif

#if defined(_BGC_)
    /* Initialize CN (Biome-BGC) module */
    InitBgc(&pihm->epctbl, &pihm->calib, pihm->elem);
#endif

#if defined(_RT_)
    InitChem(pihm->filename.cdbs, &pihm->calib, &pihm->forc, pihm->chemtbl,
        pihm->kintbl, &pihm->rttbl, &pihm->chmictbl, pihm->elem);
#endif

#if defined(_BGC_) || defined(_CYCLES_) || defined(_RT_)
    InitSolute(pihm->elem);
#endif

    /*
     * Create hydrological and land surface initial conditions
     */
    if (pihm->ctrl.init_type == RELAX)
    {
        /* Relaxation mode
         * Noah initialization needs air temperature thus forcing is applied */
#if defined(_RT_)
        ApplyForcing(pihm->ctrl.starttime, pihm->ctrl.rad_mode, &pihm->siteinfo,
            &pihm->rttbl, &pihm->forc, pihm->elem);
#elif defined(_NOAH_)
        ApplyForcing(pihm->ctrl.starttime, pihm->ctrl.rad_mode, &pihm->siteinfo,
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
    /* Initialize Cycles module */
    if (pihm->ctrl.read_cycles_restart)
    {
        ReadCyclesIc(pihm->filename.cyclesic, pihm->elem);
    }
    else
    {
        FirstDay(&pihm->soiltbl, &pihm->ctrl, pihm->elem);
    }

    InitAgVar(pihm->elem, pihm->river, CV_Y);
#endif

#if defined(_BGC_)
    /* Initialize CN variables */
    if (pihm->ctrl.read_bgc_restart)
    {
        ReadBgcIc(pihm->filename.bgcic, pihm->elem, pihm->river);
    }
    else
    {
        FirstDay(&pihm->cninit, pihm->elem, pihm->river);
    }

    InitBgcVar(pihm->elem, pihm->river, CV_Y);
#endif

#if defined(_RT_)
    if (pihm->ctrl.read_rt_restart)
    {
        ReadRtIc(pihm->filename.rtic, pihm->elem);
    }

    InitRTVar(pihm->chemtbl, &pihm->rttbl, pihm->elem, pihm->river, CV_Y);
#endif

    /* Calculate model time steps */
    CalcModelSteps(&pihm->ctrl);

#if defined(_DAILY_)
    InitDailyStruct(pihm->elem);
#endif
}

void CorrectElev(const river_struct river[], elem_struct elem[])
{
    int             i, j;
    int             sink;
    int             river_flag = 0;
    double          nabr_zmax;
    double          new_elevation;

    pihm_printf(VL_VERBOSE, "Correct surface elevation.\n");

    for (i = 0; i < nelem; i++)
    {
        /* Correction of surface elevation (artifacts due to coarse scale
         * discretization). Not needed if there is lake feature. */
        sink = 1;

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].nabr[j] != 0)
            {
                nabr_zmax = (elem[i].nabr_river[j] == 0) ?
                    elem[elem[i].nabr[j] - 1].topo.zmax :
                    river[elem[i].nabr_river[j] - 1].topo.zmax;
                if (elem[i].topo.zmax >= nabr_zmax)
                {
                    sink = 0;
                    break;
                }
            }
        }

        if (sink == 1)
        {
            pihm_printf(VL_NORMAL, "Element %4d is a sink.\n", i + 1);

            /* Note: Following correction is being applied for correction
             * mode only */
            pihm_printf(VL_NORMAL, "  Before correction: surface %7.2lf m, "
                "bedrock %7.2lf m. ", elem[i].topo.zmax, elem[i].topo.zmin);

            new_elevation = 1.0e7;
            for (j = 0; j < NUM_EDGE; j++)
            {
                if (elem[i].nabr[j] != 0)
                {
                    nabr_zmax = (elem[i].nabr_river[j] == 0) ?
                        elem[elem[i].nabr[j] - 1].topo.zmax :
                        river[elem[i].nabr_river[j] - 1].topo.zmax;
                    new_elevation = (nabr_zmax < new_elevation) ?
                        nabr_zmax : new_elevation;
                }
            }

            /* Lift bedrock elevation by the same amount */
            elem[i].topo.zmin += new_elevation - elem[i].topo.zmax;

            /* Apply new surface elevation */
            elem[i].topo.zmax = new_elevation;

            pihm_printf(VL_NORMAL, "Corrected = %7.2lf m, %7.2lf m.\n",
                elem[i].topo.zmax, elem[i].topo.zmin);
        }
    }

    for (i = 0; i < nriver; i++)
    {
        if (river[i].down > 0)
        {
            if (river[i].topo.zbed < river[river[i].down - 1].topo.zbed)
            {
                river_flag = 1;
                pihm_printf(VL_NORMAL,
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
                pihm_printf(VL_NORMAL,
                    "River outlet is higher than the channel (River %d).\n",
                    i + 1);
            }
        }
    }

    if (river_flag == 1)
    {
        pihm_printf(VL_NORMAL, "\nRiver elevation correction needed "
            "but PIHM will continue without fixing river elevation.\n\n");
    }
}

void InitSurfL(const meshtbl_struct *meshtbl, elem_struct elem[])
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

            if (elem[i].nabr[j] == 0)
            {
                elem[i].topo.x_nabr[j] = elem[i].topo.x - 2.0 * distx;
                elem[i].topo.y_nabr[j] = elem[i].topo.y - 2.0 * disty;
                elem[i].topo.dist_nabr[j] =
                    sqrt(pow(elem[i].topo.edge[0] * elem[i].topo.edge[1] *
                    elem[i].topo.edge[2] / (4.0 * elem[i].topo.area), 2) -
                    pow(elem[i].topo.edge[j] / 2.0, 2));
            }
            else
            {
                elem[i].topo.x_nabr[j] = elem[elem[i].nabr[j] - 1].topo.x;
                elem[i].topo.y_nabr[j] = elem[elem[i].nabr[j] - 1].topo.y;
                elem[i].topo.dist_nabr[j] =
                    (elem[i].topo.x - elem[i].topo.x_nabr[j]) *
                    (elem[i].topo.x - elem[i].topo.x_nabr[j]);
                elem[i].topo.dist_nabr[j] +=
                    (elem[i].topo.y - elem[i].topo.y_nabr[j]) *
                    (elem[i].topo.y - elem[i].topo.y_nabr[j]);
                elem[i].topo.dist_nabr[j] = sqrt(elem[i].topo.dist_nabr[j]);
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
                pihm_printf(VL_ERROR,
                    "Error: Return value type %d id not defined.\n", type);
                pihm_exit(EXIT_FAILURE);
        }
    }

    return (type == WS_AREA) ? ans : ans / (double)nelem;
}

void RelaxIc(elem_struct elem[], river_struct river[])
{
    int             i;
    const double    INIT_UNSAT = 0.1;
#if defined(_DGW_)
    const double    INIT_DGW = 5.0;
#endif

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem[i].ic.cmc   = 0.0;
        elem[i].ic.sneqv = 0.0;
        elem[i].ic.surf  = 0.0;
        elem[i].ic.unsat = INIT_UNSAT;
        elem[i].ic.gw    = elem[i].soil.depth - INIT_UNSAT;

#if defined(_DGW_)
        elem[i].ic.gw_geol = MIN(elem[i].geol.depth, INIT_DGW);
        elem[i].ic.unsat_geol = 0.5 * (elem[i].geol.depth - elem[i].ic.gw_geol);
#endif

#if defined(_NOAH_)
        int             j;
        double          sfctmp;

        sfctmp = elem[i].es.sfctmp;

        elem[i].ic.t1 = sfctmp;

        elem[i].ic.stc[0] = sfctmp +
            (sfctmp - elem[i].ps.tbot) / elem[i].ps.zbot *
            elem[i].ps.soil_depth[0] * 0.5;

        for (j = 1; j < MAXLYR; j++)
        {
            elem[i].ic.stc[j] = (elem[i].ps.soil_depth[j] > 0.0) ?
                 elem[i].ic.stc[j - 1] +
                (sfctmp - elem[i].ps.tbot) / elem[i].ps.zbot * 0.5 *
                (elem[i].ps.soil_depth[j - 1] + elem[i].ps.soil_depth[j]) :
                BADVAL;
        }

        for (j = 0; j < MAXLYR; j++)
        {
            elem[i].ic.smc[j] = (elem[i].ps.soil_depth[j] > 0.0) ?
                elem[i].soil.smcmax : BADVAL;
            elem[i].ic.swc[j] = (elem[i].ps.soil_depth[j] > 0.0) ?
                elem[i].soil.smcmax : BADVAL;
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
    }
}

void InitVar(elem_struct elem[], river_struct river[], N_Vector CV_Y)
{
    int             i;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    /* State variables (initial conditions) */
    for (i = 0; i < nelem; i++)
    {
#if defined(_CYCLES_OBSOLETE_)
        elem[i].ws.flatResidueWater = elem[i].ic.cmc;
#else
        elem[i].ws.cmc   = elem[i].ic.cmc;
#endif
        elem[i].ws.sneqv = elem[i].ic.sneqv;

        elem[i].ws.surf  = elem[i].ic.surf;
        elem[i].ws.unsat = elem[i].ic.unsat;
        elem[i].ws.gw    = elem[i].ic.gw;

        NV_Ith(CV_Y, SURF(i))  = elem[i].ic.surf;
        NV_Ith(CV_Y, UNSAT(i)) = elem[i].ic.unsat;
        NV_Ith(CV_Y, GW(i))    = elem[i].ic.gw;

#if defined(_DGW_)
        elem[i].ws.unsat_geol = elem[i].ic.unsat_geol;
        elem[i].ws.gw_geol    = elem[i].ic.gw_geol;

        NV_Ith(CV_Y, UNSAT_GEOL(i)) = elem[i].ic.unsat_geol;
        NV_Ith(CV_Y, GW_GEOL(i))    = elem[i].ic.gw_geol;
#endif

#if defined(_NOAH_)
        int             j;

        elem[i].es.t1    = elem[i].ic.t1;
        elem[i].ps.snowh = elem[i].ic.snowh;

        for (j = 0; j < MAXLYR; j++)
        {
            elem[i].es.stc[j] = elem[i].ic.stc[j];
            elem[i].ws.smc[j] = elem[i].ic.smc[j];
            elem[i].ws.swc[j] = elem[i].ic.swc[j];
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

        NV_Ith(CV_Y, RIVER(i)) = river[i].ic.stage;
    }

    /* Other variables */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        InitWFlux(&elem[i].wf);

#if defined(_NOAH_)
        elem[i].ps.snotime1      = 0.0;
        elem[i].ps.ribb          = 0.0;
        elem[i].ps.fcr           = 1.0;
        elem[i].ps.snoalb        = 0.75;
        elem[i].ps.zlvl          = 3.0;
        elem[i].ps.emissi        = 0.96;
        elem[i].ps.albedo        = 0.18;
        elem[i].ps.z0            = 0.1;
        elem[i].ps.ch            = 1.0e-4;
        elem[i].ps.cm            = 1.0e-4;
        elem[i].ps.beta          = BADVAL;
        elem[i].ps.sncovr        = BADVAL;
        elem[i].ps.rc            = BADVAL;
        elem[i].ps.pc            = BADVAL;
        elem[i].ps.rcs           = BADVAL;
        elem[i].ps.rct           = BADVAL;
        elem[i].ps.rcsoil        = BADVAL;
        elem[i].ps.q1            = BADVAL;
        elem[i].ps.z0brd         = BADVAL;
        elem[i].ps.eta_kinematic = BADVAL;

        elem[i].ef.sheat = BADVAL;
        elem[i].ef.eta   = BADVAL;
        elem[i].ef.fdown = BADVAL;
        elem[i].ef.ec    = BADVAL;
        elem[i].ef.edir  = BADVAL;
        elem[i].ef.ett   = BADVAL;
        elem[i].ef.etp   = BADVAL;
        elem[i].ef.ssoil = BADVAL;
        elem[i].ef.flx1  = BADVAL;
        elem[i].ef.flx2  = BADVAL;
        elem[i].ef.flx3  = BADVAL;
        elem[i].ef.esnow = BADVAL;

        elem[i].wf.runoff2 = BADVAL;
        elem[i].wf.runoff3 = BADVAL;
        elem[i].wf.pcpdrp  = 0.0;
        elem[i].wf.drip    = 0.0;
        elem[i].wf.dew     = BADVAL;
        elem[i].wf.snomlt  = BADVAL;

        elem[i].ws.soilm = BADVAL;
#endif
    }
}

void CalcModelSteps(ctrl_struct *ctrl)
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
        wf->overland[j] = 0.0;
        wf->subsurf[j] = 0.0;
    }
    wf->prcp       = 0.0;
    wf->pcpdrp     = 0.0;
    wf->infil      = 0.0;
    wf->eqv_infil  = 0.0;
    wf->recharge   = 0.0;
    wf->drip       = 0.0;
    wf->edir       = 0.0;
    wf->ett        = 0.0;
    wf->ec         = 0.0;
    wf->etp        = 0.0;
    wf->eta        = 0.0;
    wf->edir_surf  = 0.0;
    wf->edir_unsat = 0.0;
    wf->edir_gw    = 0.0;
    wf->ett_unsat  = 0.0;
    wf->ett_gw     = 0.0;
    wf->esnow      = 0.0;

#if defined(_DGW_)
    wf->infil_geol = 0.0;
    wf->rechg_geol = 0.0;
    for (j = 0; j < NUM_EDGE; j++)
    {
        wf->dgw[j] = 0.0;
    }
# if defined(_LUMPED_)
    wf->dgw_runoff = 0.0;
# endif
#endif

#if defined(_NOAH_)
    int             k;

    for (k = 0; k < MAXLYR; k++)
    {
        wf->et[k]          = 0.0;
        wf->runoff2_lyr[k] = 0.0;
        wf->smflx[k]       = 0.0;
    }
    wf->runoff2 = 0.0;
    wf->runoff3 = 0.0;
    wf->dew     = 0.0;
    wf->snomlt  = 0.0;
    wf->etns    = 0.0;
#endif
#if defined(_CYCLES_)
    wf->irrig   = 0.0;
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

    for (k = 0; k < MAXLYR; k++)
    {
        ef->et[k] = 0.0;
    }
    ef->solnet   = 0.0;
    ef->etp      = 0.0;
    ef->ssoil    = 0.0;
    ef->eta      = 0.0;
    ef->sheat    = 0.0;
    ef->fdown    = 0.0;
    ef->lwdn     = 0.0;
    ef->ec       = 0.0;
    ef->edir     = 0.0;
    ef->ett      = 0.0;
    ef->esnow    = 0.0;
    ef->soldir   = 0.0;
    ef->soldif   = 0.0;
    ef->longwave = 0.0;
    ef->flx1     = 0.0;
    ef->flx2     = 0.0;
    ef->flx3     = 0.0;
#endif
#if defined(_BGC_)
    ef->swabs_per_plaisun   = 0.0;
    ef->swabs_per_plaishade = 0.0;
#endif
}
