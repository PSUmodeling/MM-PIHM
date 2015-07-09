/*****************************************************************************
 * File		: initialize.c							
 * Function	: Initialization of elemental attributes using relational	
 *		  database							
 ****************************************************************************/

#include "pihm.h"

void Initialize (pihm_struct pihm, N_Vector CV_Y)
{
    int             i;

    if (verbose_mode)
    {
        printf ("\n\nInitialize data structure\n");
    }

    pihm->elem = (elem_struct *) malloc (pihm->numele * sizeof (elem_struct));
    pihm->riv = (river_struct *) malloc (pihm->numriv * sizeof (river_struct));

#ifdef _NOAH_
    //pihm->avg_inf = (double *) malloc (pihm->NumEle * sizeof (double));  /* YS */
    //pihm->avg_rech = (double *) malloc (pihm->NumEle * sizeof (double));  /* YS */
    //pihm->avg_subflux = (double **) malloc (pihm->NumEle * sizeof (double *));  /* YS */
    for (i = 0; i < pihm->numele; i++)
    {
        pihm->elem[i].fcr = 1.0;
	//pihm->avg_inf[i] = 0.0;
	//pihm->avg_rech[i] = 0.0;
	//pihm->avg_subflux[i] = (double *) malloc (3 * sizeof (double));
	//pihm->avg_subflux[i][0] = 0.0;
	//pihm->avg_subflux[i][1] = 0.0;
	//pihm->avg_subflux[i][2] = 0.0;
    }
#endif

    InitMeshStruct (pihm->elem, pihm->numele, pihm->mesh_tbl);

    InitTopo (pihm->elem, pihm->numele, pihm->mesh_tbl);

    InitSoil (pihm->elem, pihm->numele, pihm->attrib_tbl, pihm->soil_tbl,
        pihm->geol_tbl, pihm->cal);

    InitLC (pihm->elem, pihm->numele, pihm->attrib_tbl, pihm->lc_tbl,
        pihm->cal);

    InitForcing (pihm->elem, pihm->numele, pihm->riv, pihm->numriv,
        pihm->attrib_tbl, pihm->riv_att_tbl, &pihm->forcing, pihm->cal);

    InitRiver (pihm->riv, pihm->numriv, pihm->elem, pihm->riv_att_tbl, pihm->riv_shp_tbl, pihm->riv_matl_tbl, pihm->mesh_tbl, pihm->cal);

    if (debug_mode)
    {
        CorrectElevation (pihm->elem, pihm->numele, pihm->riv, pihm->numriv);
    }

    InitSurfL (pihm->elem, pihm->numele, pihm->riv, pihm->mesh_tbl);

    if (pihm->ctrl.init_type == 0)
    {
        SaturationIC (pihm->elem, pihm->numele, pihm->riv, pihm->numriv, &pihm->ic);
    }
    else if (pihm->ctrl.init_type == 1)
    {
        for (i = 0; i < pihm->numriv; i++)
        {
            pihm->ic.rivgw[i] = pihm->riv[i].topo.zmax - pihm->riv[i].topo.zmin - 0.1;
        }
    }

    InitStateVrbl (pihm->elem, pihm->numele, pihm->riv, pihm->numriv, CV_Y, pihm->ic);

    CalcModelStep (&pihm->ctrl);
}
    
void InitMeshStruct (elem_struct *elem, int numele, mesh_tbl_struct mesh_tbl)
{
    int             i, j;

    for (i = 0; i < numele; i++)
    {
        for (j = 0; j < 3; j++)
        {
            elem[i].node[j] = mesh_tbl.node[i][j];
            elem[i].nabr[j] = mesh_tbl.nabr[i][j];
        }
    }
}

void InitTopo (elem_struct *elem, int numele, mesh_tbl_struct mesh_tbl)
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
            x[j] = mesh_tbl.x[elem[i].node[j] - 1];
            y[j] = mesh_tbl.y[elem[i].node[j] - 1];
            zmin[j] = mesh_tbl.zmin[elem[i].node[j] - 1];
            zmax[j] = mesh_tbl.zmax[elem[i].node[j] - 1];
        }

        elem[i].topo.area = 0.5 * ((x[1] - x[0]) * (y[2] - y[0]) - (y[1] - y[0]) * (x[2] - x[0]));
        /* Calculate centroid of triangle */
        elem[i].topo.x = (x[0] + x[1] + x[2]) / 3.0;
        elem[i].topo.y = (y[0] + y[1] + y[2]) / 3.0;

        elem[i].topo.zmin = (zmin[0] + zmin[1] + zmin[2]) / 3.0;
        elem[i].topo.zmax = (zmax[0] + zmax[1] + zmax[2]) / 3.0;
        elem[i].topo.edge[0] = sqrt (pow ((x[1] - x[2]), 2) + pow ((y[1] - y[2]), 2));
        elem[i].topo.edge[1] = sqrt (pow ((x[2] - x[0]), 2) + pow ((y[2] - y[0]), 2));
        elem[i].topo.edge[2] = sqrt (pow ((x[0] - x[1]), 2) + pow ((y[0] - y[1]), 2));

    }
}

void InitSoil (elem_struct *elem, int numele, attrib_tbl_struct attrib_tbl, soil_tbl_struct soil_tbl, geol_tbl_struct geol_tbl, calib_struct cal)
{
    int             i;
    int             soil_ind;
    int             geol_ind;
    double          thetaw;
    double          thetaref;

    for (i = 0; i < numele; i++)
    {
        soil_ind = attrib_tbl.soil[i] - 1;
        geol_ind = attrib_tbl.geol[i] - 1;

        elem[i].soil.depth = elem[i].topo.zmax - elem[i].topo.zmin;
        elem[i].soil.ksath = cal.ksath * geol_tbl.ksath[geol_ind];
        elem[i].soil.ksatv = cal.ksatv * geol_tbl.ksatv[geol_ind];
        elem[i].soil.kinfv = cal.kinfv * soil_tbl.ksatv[soil_ind];
        elem[i].soil.porosity = cal.porosity * (soil_tbl.thetas[soil_ind] - soil_tbl.thetar[soil_ind]);
        if (elem[i].soil.porosity > 1.0 || elem[i].soil.porosity <= 0.0)
        {
            printf ("Warning: Porosity value out of bounds for element %d", i + 1);
            exit (1);
        }
        elem[i].soil.dinf = cal.dinf * soil_tbl.dinf[soil_ind];
        elem[i].soil.alpha = cal.alpha * soil_tbl.alpha[soil_ind];
        elem[i].soil.beta = cal.beta * soil_tbl.beta[soil_ind];

        elem[i].soil.thetar = soil_tbl.thetar[soil_ind];
        elem[i].soil.thetas = elem[i].soil.thetar + elem[i].soil.porosity;

        /* Calculate field capacity and wilting point following Chan and Dudhia 2001 MWR, but replacing Campbell with van Genuchten */
        elem[i].soil.thetaw = WiltingPoint (soil_tbl.thetas[soil_ind], soil_tbl.thetar[soil_ind], soil_tbl.alpha[soil_ind], soil_tbl.beta[soil_ind]);
        elem[i].soil.thetaref = FieldCapacity (soil_tbl.alpha[soil_ind], soil_tbl.beta[soil_ind], geol_tbl.ksatv[geol_ind], soil_tbl.thetas[soil_ind], soil_tbl.thetar[soil_ind]);

        elem[i].soil.dmac = cal.dmac * geol_tbl.dmac[geol_ind];
        if (elem[i].soil.dmac > elem[i].soil.depth)
            elem[i].soil.dmac = elem[i].soil.depth;

        elem[i].soil.kmacv = cal.kmacv * soil_tbl.kmacv[soil_ind];
        elem[i].soil.kmach = cal.kmach * geol_tbl.kmach[geol_ind];
        elem[i].soil.areafh = cal.areafh * soil_tbl.areafh[soil_ind];
        elem[i].soil.areafv = cal.areafv * geol_tbl.areafv[geol_ind];
        elem[i].soil.macropore = attrib_tbl.macropore[i];
    }
}

double FieldCapacity (double alpha, double beta, double kv, double thetas, double thetar)
{
    double        satn;
    double        ktemp;
    double        thetaref;
    const double  KFC = 5.79E-9;

    /* Set default value of field capacity, in case a valid value cannot be
     * found using the defination */
    thetaref = 0.75 * thetas;

    for (satn = 0.005; satn < 1.0; satn = satn + 0.001)
    {
        ktemp = kv * KrFunc (alpha, beta, satn);
        if (ktemp >= KFC)
        {
            thetaref = (1.0 / 3.0 + 2.0 / 3.0 * satn) * (thetas - thetar) + thetar;
            break;
        }
    }
    return (thetaref);
}

double WiltingPoint (double thetas, double thetar, double alpha, double beta)
{
    const double    PSIW = 200.0;

    return (0.5 * (thetas - thetar) * pow (1.0 / (1.0 + pow (PSIW * alpha, beta)), 1.0 - 1.0 / beta) + thetar);
}

void InitLC (elem_struct *elem, int numele, attrib_tbl_struct attrib_tbl, lc_tbl_struct lc_tbl, calib_struct cal)
{
    int             i;
    int             lc_ind;

    for (i = 0; i < numele; i++)
    {
        lc_ind = attrib_tbl.lc[i] - 1;
        elem[i].lc.type = lc_ind + 1;
        elem[i].lc.vegfrac = cal.vegfrac * lc_tbl.vegfrac[lc_ind];
        elem[i].lc.rzd = cal.rzd * lc_tbl.rzd[lc_ind];
        elem[i].lc.rsmin = lc_tbl.rsmin[lc_ind];
        elem[i].lc.rgl = lc_tbl.rgl[lc_ind];
        elem[i].lc.hs = lc_tbl.hs[lc_ind];
        elem[i].lc.laimin = lc_tbl.laimin[lc_ind];
        elem[i].lc.laimax = lc_tbl.laimax[lc_ind];
        elem[i].lc.emissmin = lc_tbl.emissmin[lc_ind];
        elem[i].lc.emissmax = lc_tbl.emissmax[lc_ind];
        elem[i].lc.albedomin = cal.albedo * lc_tbl.albedomin[lc_ind];
        elem[i].lc.albedomax = cal.albedo * lc_tbl.albedomax[lc_ind];
        //elem[i].lc.albedo = 0.5 * elem[i].lc.albedo_min + 0.5 * elem[i].lc.albedo_max;
        //if (elem[i].lc.albedo > 1.0 || elem[i].lc.albedo < 0.0)
        //{
        //    printf ("Warning: Albedo out of bounds");
        //    exit (1);
        //}
        elem[i].lc.z0min = lc_tbl.z0min[lc_ind];
        elem[i].lc.z0max = lc_tbl.z0max[lc_ind];
        elem[i].lc.rough = cal.rough * lc_tbl.rough[lc_ind];
        elem[i].lc.intcp_factr = 0.0002;
        elem[i].lc.rsmax = lc_tbl.rsmax;
        elem[i].lc.bare = lc_tbl.bare;
        elem[i].lc.cfactr = lc_tbl.cfactr;
        elem[i].lc.topt = lc_tbl.topt;

//#ifdef _NOAH_
//        elem[i].lc.rsmin *= cal.rsmin;
//        elem[i].lc.intcp_factr *= cal.intcp;
//        elem[i].lc.cfactr *= cal.cfactr;
//        elem[i].lc.rgl *= cal.rgl;
//        elem[i].lc.hs *= cal.hs;
//#endif
    }
}

void InitRiver (river_struct *riv, int numriv, elem_struct *elem, riv_att_tbl_struct riv_att_tbl, riv_shp_tbl_struct riv_shp_tbl, riv_matl_tbl_struct riv_matl_tbl, mesh_tbl_struct mesh_tbl, calib_struct cal)
{
    int             i, j;

    for (i = 0; i < numriv; i++)
    {
        riv[i].leftele = riv_att_tbl.leftele[i];
        riv[i].rightele = riv_att_tbl.rightele[i];
        riv[i].fromnode = riv_att_tbl.fromnode[i];
        riv[i].tonode = riv_att_tbl.tonode[i];
        riv[i].down = riv_att_tbl.down[i];

        for (j = 0; j < 3; j++)
        {
            /* Note: Strategy to use BC < -4 for river identification */
            if (elem[riv[i].leftele - 1].nabr[j] == riv[i].rightele)
                elem[riv[i].leftele - 1].forc.bc_type[j] = -4 * (i + 1);
            if (elem[riv[i].rightele - 1].nabr[j] == riv[i].leftele)
                elem[riv[i].rightele - 1].forc.bc_type[j] = -4 * (i + 1);
        }

        riv[i].topo.x = (mesh_tbl.x[riv[i].fromnode - 1] + mesh_tbl.x[riv[i].tonode - 1]) / 2.0;
        riv[i].topo.y = (mesh_tbl.y[riv[i].fromnode - 1] + mesh_tbl.y[riv[i].tonode - 1]) / 2.0;
        riv[i].topo.zmax = (mesh_tbl.zmax[riv[i].fromnode - 1] + mesh_tbl.zmax[riv[i].tonode - 1]) / 2.0;
        riv[i].topo.zmin = riv[i].topo.zmax - (0.5 * (elem[riv[i].leftele - 1].topo.zmax + elem[riv[i].rightele - 1].topo.zmax)
            - 0.5 * (elem[riv[i].leftele - 1].topo.zmin + elem[riv[i].rightele - 1].topo.zmin));
        riv[i].topo.node_zmax = mesh_tbl.zmax[riv[i].tonode - 1];

        riv[i].shp.depth = cal.rivdepth * riv_shp_tbl.depth[riv_att_tbl.shp[i] - 1];
        riv[i].shp.intrpl_ord = riv_shp_tbl.intrpl_ord[riv_att_tbl.shp[i] - 1];
        riv[i].shp.coeff = cal.rivshpcoeff * riv_shp_tbl.coeff[riv_att_tbl.shp[i] - 1];
        riv[i].shp.length = sqrt (pow (mesh_tbl.x[riv[i].fromnode - 1] - mesh_tbl.x[riv[i].tonode - 1], 2) + pow (mesh_tbl.y[riv[i].fromnode - 1] - mesh_tbl.y[riv[i].tonode - 1], 2));

        riv[i].topo.zbed = riv[i].topo.zmax - riv[i].shp.depth;

        riv[i].matl.rough = cal.rivrough * riv_matl_tbl.rough[riv_att_tbl.matl[i] - 1];
        riv[i].matl.cwr = riv_matl_tbl.cwr[riv_att_tbl.matl[i] - 1];
        riv[i].matl.ksath = cal.rivksath * riv_matl_tbl.ksath[riv_att_tbl.matl[i] - 1];
        riv[i].matl.ksatv = cal.rivksatv * riv_matl_tbl.ksatv[riv_att_tbl.matl[i] - 1];
        riv[i].matl.bedthick = cal.rivbedthick * riv_matl_tbl.bedthick[riv_att_tbl.matl[i] - 1];
        riv[i].matl.porosity = 0.5 * (elem[riv[i].leftele - 1].soil.porosity +  elem[riv[i].rightele - 1].soil.porosity);
    }
}

void InitForcing (elem_struct *elem, int numele, river_struct *riv, int numriv, attrib_tbl_struct attrib_tbl, riv_att_tbl_struct riv_att_tbl, forcing_ts_struct *forcing, calib_struct cal)
{
    int             i, j;

    for (i = 0; i < forcing->nts[METEO_TS]; i++)
    {
        for (j = 0; j < forcing->ts[METEO_TS][i].length; j++)
        {
            forcing->ts[METEO_TS][i].data[j][PRCP_TS] *= cal.prcp;
            forcing->ts[METEO_TS][i].data[j][SFCTMP_TS] *= cal.sfctmp;
        }
    }

    forcing->bc = (double *) malloc (forcing->nts[BC_TS] * sizeof (double));
    for (i = 0; i < NUM_METEO_TS; i++)
    {
        forcing->meteo[i] = (double *) malloc (forcing->nts[METEO_TS] * sizeof (double));
    }
    forcing->lai = (double *) malloc (forcing->nts[LAI_TS] * sizeof (double));
    forcing->riverbc = (double *) malloc (forcing->nts[RIV_TS] * sizeof (double));
    forcing->source = (double *) malloc (forcing->nts[SS_TS] * sizeof (double));

    for (i = 0; i < numele; i++)
    {
        for (j = 0; j < 3; j++)
        {
            elem[i].forc.bc_type[j] = attrib_tbl.bc[i][j];
            if (elem[i].forc.bc_type[j] > 0)
            {
                elem[i].forc.bc[j] = &forcing->bc[elem[i].forc.bc_type[j] - 1];
            }
        }

        for (j = 0; j < NUM_METEO_TS; j++)
        {
            elem[i].forc.meteo[j] = &forcing->meteo[j][attrib_tbl.meteo[i] - 1];
        }
        elem[i].forc.zlvl_wind = forcing->zlvl_wind[attrib_tbl.meteo[i] - 1];

        elem[i].forc.lai_type = attrib_tbl.lai[i];
        if (elem[i].forc.lai_type > 0)
        {
            elem[i].forc.lai = &forcing->lai[attrib_tbl.lai[i] - 1];
        }

        if (attrib_tbl.source[i] > 0)
        {
            elem[i].forc.source = &forcing->source[attrib_tbl.source[i] - 1];
        }

    }

    for (i = 0; i < numriv; i++)
    {
        if (riv_att_tbl.bc[i] > 0)
        {
            riv[i].forc.riverbc = &forcing->riverbc[riv_att_tbl.bc[i] - 1];
        }
    }
}

void CorrectElevation (elem_struct *elem, int numele, river_struct *riv, int numriv)
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
            if (elem[i].nabr[j] > 0)
            {
                new_elevation = elem[i].forc.bc_type[j] > -4 ? elem[elem[i].nabr[j] - 1].topo.zmax : riv[-(elem[i].forc.bc_type[j] / 4) - 1].topo.zmax;
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
                if (elem[i].nabr[j] > 0)
                {
                    elem[i].topo.zmax = (elem[i].forc.bc_type[j] > -4 ? elem[elem[i].nabr[j] - 1].topo.zmax : riv[-(elem[i].forc.bc_type[j] / 4) - 1].topo.zmax);
                    new_elevation = new_elevation > elem[i].topo.zmax ? elem[i].topo.zmax : new_elevation;
                    printf ("(%d)%lf  ", j + 1, (elem[i].forc.bc_type[j] > -4 ? elem[elem[i].nabr[j] - 1].topo.zmax : riv[-(elem[i].forc.bc_type[j] / 4) - 1].topo.zmax));
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
                if (elem[i].nabr[j] > 0)
                {
                    new_elevation = elem[i].forc.bc_type[j] > -4 ? elem[elem[i].nabr[j] - 1].topo.zmin : riv[-(elem[i].forc.bc_type[j] / 4) - 1].topo.zmin;
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
                /* Note: Following correction is being applied for debug==1 case only */
                printf ("\tBefore: %lf Corrected using:", elem[i].topo.zmin);
                new_elevation = 1.0e7;
                for (j = 0; j < 3; j++)
                {
                    if (elem[i].nabr[j] > 0)
                    {
                        elem[i].topo.zmin = (elem[i].forc.bc_type[j] > -4 ? elem[elem[i].nabr[j] - 1].topo.zmin : riv[-(elem[i].forc.bc_type[j] / 4) - 1].topo.zmin);
                        new_elevation = new_elevation > elem[i].topo.zmin ? elem[i].topo.zmin : new_elevation;
                        printf ("(%d)%lf  ", j + 1, (elem[i].forc.bc_type[j] > -4 ? elem[elem[i].nabr[j] - 1].topo.zmin : riv[-(elem[i].forc.bc_type[j] / 4) - 1].topo.zmin));
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

void InitSurfL (elem_struct *ele, int numele, river_struct *riv, mesh_tbl_struct mesh_tbl)
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
            x[j] = mesh_tbl.x[ele[i].node[j] - 1];
            y[j] = mesh_tbl.y[ele[i].node[j] - 1];
            zmin[j] = mesh_tbl.zmin[ele[i].node[j] - 1];
            zmax[j] = mesh_tbl.zmax[ele[i].node[j] - 1];
        }

        for (j = 0; j < 3; j++)
        {
            /* Note: Assumption here is that the forumulation is circumcenter based */
            switch (j)
            {
                case 0:
                    distx = (ele[i].topo.x - 0.5 * (x[1] + x[2]));
                    disty = (ele[i].topo.y - 0.5 * (y[1] + y[2]));
                    break;
                case 1:
                    distx = (ele[i].topo.x - 0.5 * (x[2] + x[0]));
                    disty = (ele[i].topo.y - 0.5 * (y[2] + y[0]));
                    break;
                case 2:
                    distx = (ele[i].topo.x - 0.5 * (x[0] + x[1]));
                    disty = (ele[i].topo.y - 0.5 * (y[0] + y[1]));
                    break;
            }
            ele[i].topo.surfx[j] = ele[i].nabr[j] > 0 ? (ele[i].forc.bc_type[j] > -4 ? ele[ele[i].nabr[j] - 1].topo.x : riv[-(ele[i].forc.bc_type[j] / 4) - 1].topo.x) : (ele[i].topo.x - 2.0 * distx);
            ele[i].topo.surfy[j] = ele[i].nabr[j] > 0 ? (ele[i].forc.bc_type[j] > -4 ? ele[ele[i].nabr[j] - 1].topo.y : riv[-(ele[i].forc.bc_type[j] / 4) - 1].topo.y) : (ele[i].topo.y - 2.0 * disty);
        }
    }
}

void SaturationIC (const elem_struct *ele, int numele, const river_struct *riv, int numriv, ic_struct *ic)
{
    int             i;

    for (i = 0; i < numele; i++)
    {
        ic->intcp[i] = 0.0;
        ic->snow[i] = 0.0;
        ic->surf[i] = 0.0;
        ic->unsat[i] = 0.1;
        ic->gw[i] = ele[i].soil.depth - 0.1;
    }

    for (i = 0; i < numriv; i++)
    {
        ic->stage[i] = 0.0;
        ic->rivgw[i] = riv[i].topo.zbed - riv[i].topo.zmin - 0.1;
    }
}

void InitStateVrbl (elem_struct *elem, int numele, river_struct *riv, int numriv, N_Vector CV_Y, ic_struct ic)
{
    int             i;

    for (i = 0; i < numele; i++)
    {
        elem[i].intcp = ic.intcp[i];
        elem[i].snow = ic.snow[i];

        elem[i].surf0 = ic.surf[i];
        elem[i].unsat0 = ic.unsat[i];
        elem[i].gw0 = ic.gw[i];

        elem[i].surf = ic.surf[i];
        elem[i].unsat = ic.unsat[i];
        elem[i].gw = ic.gw[i];

        NV_Ith_S (CV_Y, i) = ic.surf[i];
        NV_Ith_S (CV_Y, i + numele) = ic.unsat[i];
        NV_Ith_S (CV_Y, i + 2 * numele) = ic.gw[i];

        elem[i].runoff = 0.0;
        elem[i].infil = 0.0;
    }    

    for (i = 0; i < numriv; i++)
    {
        riv[i].stage0 = ic.stage[i];
        riv[i].gw0 = ic.rivgw[i];

        riv[i].stage = ic.stage[i];
        riv[i].gw = ic.rivgw[i];

        NV_Ith_S (CV_Y, i + 3 * numele) = ic.stage[i];
        NV_Ith_S (CV_Y, i + 3 * numele + numriv) = ic.rivgw[i];
    }
}

void CalcModelStep (ctrl_struct *ctrl)
{
    int                 i;

    ctrl->nstep = (ctrl->endtime - ctrl->starttime) / ctrl->stepsize;

    ctrl->tout = (int *) malloc ((ctrl->nstep + 1) * sizeof (int));

    for (i = 0; i < ctrl->nstep + 1; i++)
    {
        if (i == 0)
            ctrl->tout[i] = ctrl->starttime;
        else
            ctrl->tout[i] = ctrl->tout[i - 1] + ctrl->stepsize;
    }

    if (ctrl->tout[ctrl->nstep] < ctrl->endtime)
        ctrl->tout[ctrl->nstep] = ctrl->endtime;
}

void MapOutput (char *simulation, pihm_struct pihm, char *outputdir)
{
    int             i, j, k;
    int             n;

    if (verbose_mode)
        printf ("\nInitializing output files\n");

    n = 0;

    for (i = 0; i < NUM_PRINT; i++)
    {
        if (pihm->ctrl.prtvrbl[i] > 0)
        {
            switch (i)
            {
                case GW_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.gw", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele + pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].gw0;
                    }
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j + pihm->numele] = &pihm->riv[j].gw0;
                    }
                    n++;
                    break;
                case SURF_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.surf", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].surf0;
                    }
                    n++;
                    break;
                case SNOW_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.snow", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].snow;
                    }
                    n++;
                    break;
                case RIVSTG_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.stage", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].stage0;
                    }
                    n++;
                    break;
                case INFIL_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.infil", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].infil;
                    }
                    n++;
                    break;
                case RECHARGE_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.recharge", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].rechg;
                    }
                    n++;
                    break;
                case CMC_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.is", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].intcp;
                    }
                    n++;
                    break;
                case UNSAT_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.unsat", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].unsat0;
                    }
                    n++;
                    break;
                case EC_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.et0", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].et[0];
                    }
                    n++;
                    break;
                case ETT_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.et1", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].et[1];
                    }
                    n++;
                    break;
                case EDIR_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.et2", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].et[2];
                    }
                    n++;
                    break;
                case RIVFLX0_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx0", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].fluxriv[0];
                    }
                    n++;
                    break;
                case RIVFLX1_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx1", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].fluxriv[1];
                    }
                    n++;
                    break;
                case RIVFLX2_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx2", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].fluxriv[2];
                    }
                    n++;
                    break;
                case RIVFLX3_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx3", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].fluxriv[3];
                    }
                    n++;
                    break;
                case RIVFLX4_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx4", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].fluxriv[4];
                    }
                    n++;
                    break;
                case RIVFLX5_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx5", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].fluxriv[5];
                    }
                    n++;
                    break;
                case RIVFLX6_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx6", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].fluxriv[6];
                    }
                    n++;
                    break;
                case RIVFLX7_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx7", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].fluxriv[7];
                    }
                    n++;
                    break;
                case RIVFLX8_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx8", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].fluxriv[8];
                    }
                    n++;
                    break;
                case RIVFLX9_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx9", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].fluxriv[9];
                    }
                    n++;
                    break;
                case RIVFLX10_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx10", outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].fluxriv[10];
                    }
                    n++;
                    break;
                case SUBFLX_CTRL:
                    for (k = 0; k < 3; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.subflx%d", outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvrbl = pihm->numele;
                        pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].fluxsub[k];
                        }
                        n++;
                    }
                    break;
                case TOTALFLX_CTRL:
                    for (k = 0; k < 3; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.totalflx%d", outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvrbl = pihm->numele;
                        pihm->prtctrl[n].vrbl = (double **) malloc (pihm->prtctrl[n].nvrbl * sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].fluxtotal[k];
                        }
                        n++;
                    }
                    break;
                default:
                    break;
            }
        }
    }

    pihm->ctrl.nprint = n;
}

void InitOutputFile (prtctrl_struct *prtctrl, int nprint, int ascii)
{
    FILE           *fid;
    char            ascii_fn[MAXSTRING];
    int             i;

    for (i = 0; i < nprint; i++)
    {
        prtctrl[i].buffer = (double *) calloc (prtctrl[i].nvrbl, sizeof (double));

        fid = fopen (prtctrl[i].name, "w");
        fclose (fid);

        if (ascii)
        {
            sprintf (ascii_fn, "%s.txt", prtctrl[i].name);
            fid = fopen (ascii_fn, "w");
            fclose (fid);
        }
    }
}
