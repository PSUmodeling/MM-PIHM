
/*****************************************************************************
 * File        : read_alloc.c
 * Function    : read parameter files
 *----------------------------------------------------------------------------
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0...........................
 * a) Addition of three new input files: file.calib, file.lc and file.geol
 * b) Declaration and allocation  of new variables for new process, shape
 *    representations  and calibration (e.g. ET, Infiltration, Macropore,
 *    Stormflow, Element beneath river, river shapes, river bed property,
 *    thresholds for root zone, infiltration and macropore depths, land cover
 *    attributes etc)                                                      
 ****************************************************************************/

#include "pihm.h"

void PihmFree (void **ptr)
{
    free (*ptr);
    *ptr = NULL;
}

void ReadAlloc (char *simulation, pihm_struct pihm)
{
    char            project[MAXSTRING];
    char           *token;
    char            tempname[MAXSTRING];

    strcpy (tempname, simulation);
    if (strstr (tempname, ".") != 0)
    {
        token = strtok (tempname, ".");
        strcpy (project, token);
    }
    else
    {
        strcpy (project, simulation);
    }

    if (verbose_mode)
    {
        printf ("\nRead input files:\n");
    }

    ReadRiv (project, &pihm->riv_att_tbl, &pihm->riv_shp_tbl,
        &pihm->riv_matl_tbl, &pihm->riv_ic_tbl, &pihm->ic, &pihm->forcing);
    pihm->numriv = pihm->riv_att_tbl.number;

    ReadMesh (project, &pihm->mesh_tbl);
    pihm->numele = pihm->mesh_tbl.numele;

    ReadAtt (project, &pihm->attrib_tbl, &pihm->ic, pihm->numele);

    ReadSoil (project, &pihm->soil_tbl);

    ReadGeol (project, &pihm->geol_tbl);

    ReadLC (&pihm->lc_tbl);

    ReadForc (project, &pihm->forcing);

    ReadLAI (project, &pihm->forcing, pihm->numele, &pihm->attrib_tbl);

    /* Read source and sink */
    //ReadSS ();
    pihm->forcing.nts[SS_TS] = 0;

    ReadIbc (project, &pihm->forcing);

    ReadPara (project, &pihm->ctrl);

    ReadCalib (project, simulation, &pihm->cal);

    if (pihm->ctrl.init_type == 3)
    {
        ReadInit (project, simulation, &pihm->ic, pihm->numele, pihm->numriv);
    }
}

void ReadRiv (char *project, riv_att_tbl_struct *riv_att_tbl,
    riv_shp_tbl_struct *riv_shp_tbl, riv_matl_tbl_struct *riv_matl_tbl,
    riv_ic_tbl_struct * riv_ic_tbl, ic_struct *ic, forcing_ts_struct *forcing)
{
    int             i, j;
    char            fn[MAXSTRING];
    FILE           *riv_file;   /* Pointer to .riv file */
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    /*
     * Open .riv input file
     */
    sprintf (fn, "input/%s/%s.riv", project, project);
    riv_file = fopen (fn, "r");
    CheckFile (riv_file, fn);

    /*
     * Read river segment block
     */

    /* Read number of river segments */
    FindLine (riv_file, "BOF");
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%d", &riv_att_tbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of river segments!\n");
        printf (".riv file format error!\n");
        exit (1);
    }

    /* Allocate */
    riv_att_tbl->fromnode =
        (int *)malloc (riv_att_tbl->number * sizeof (int));
    riv_att_tbl->tonode = (int *)malloc (riv_att_tbl->number * sizeof (int));
    riv_att_tbl->down = (int *)malloc (riv_att_tbl->number * sizeof (int));
    riv_att_tbl->leftele = (int *)malloc (riv_att_tbl->number * sizeof (int));
    riv_att_tbl->rightele =
        (int *)malloc (riv_att_tbl->number * sizeof (int));
    riv_att_tbl->shp = (int *)malloc (riv_att_tbl->number * sizeof (int));
    riv_att_tbl->matl = (int *)malloc (riv_att_tbl->number * sizeof (int));
    riv_att_tbl->ic = (int *)malloc (riv_att_tbl->number * sizeof (int));
    riv_att_tbl->bc = (int *)malloc (riv_att_tbl->number * sizeof (int));
    riv_att_tbl->rsvr = (int *)malloc (riv_att_tbl->number * sizeof (int));

    /* Read river segment information */
    for (i = 0; i < riv_att_tbl->number; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%d %d %d %d %d %d %d %d %d %d %d",
            &index,
            &riv_att_tbl->fromnode[i], &riv_att_tbl->tonode[i],
            &riv_att_tbl->down[i],
            &riv_att_tbl->leftele[i], &riv_att_tbl->rightele[i],
            &riv_att_tbl->shp[i], &riv_att_tbl->matl[i],
            &riv_att_tbl->ic[i], &riv_att_tbl->bc[i], &riv_att_tbl->bc[i]);
        if (match != 11 || i != index - 1)
        {
            printf
                ("Cannot read river segment information for the %dth segment!\n",
                i + 1);
            printf (".riv file format error!\n");
            exit (1);
        }
    }

    /*
     * Read river shape information
     */
    FindLine (riv_file, "SHAPE");
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%d", &riv_shp_tbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of river shapes!\n");
        printf (".riv file format error!\n");
        exit (1);
    }

    /* Allocate */
    riv_shp_tbl->depth =
        (double *)malloc (riv_shp_tbl->number * sizeof (double));
    riv_shp_tbl->intrpl_ord =
        (int *)malloc (riv_shp_tbl->number * sizeof (int));
    riv_shp_tbl->coeff =
        (double *)malloc (riv_shp_tbl->number * sizeof (double));

    for (i = 0; i < riv_shp_tbl->number; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %d %lf",
            &index, &riv_shp_tbl->depth[i],
            &riv_shp_tbl->intrpl_ord[i], &riv_shp_tbl->coeff[i]);
        if (match != 4 || i != index - 1)
        {
            printf
                ("Cannot read river shape information for the %dth shape!\n",
                i + 1);
            printf (".riv file format error!\n");
            exit (1);
        }
    }

    /*
     * Read river material information
     */
    FindLine (riv_file, "MATERIAL");
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%d", &riv_matl_tbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of river materials!\n");
        printf (".riv file format error!\n");
        exit (1);
    }

    /* Allocate */
    riv_matl_tbl->rough =
        (double *)malloc (riv_matl_tbl->number * sizeof (double));
    riv_matl_tbl->cwr =
        (double *)malloc (riv_matl_tbl->number * sizeof (double));
    riv_matl_tbl->ksath =
        (double *)malloc (riv_matl_tbl->number * sizeof (double));
    riv_matl_tbl->ksatv =
        (double *)malloc (riv_matl_tbl->number * sizeof (double));
    riv_matl_tbl->bedthick =
        (double *)malloc (riv_matl_tbl->number * sizeof (double));

    for (i = 0; i < riv_matl_tbl->number; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf",
            &index,
            &riv_matl_tbl->rough[i], &riv_matl_tbl->cwr[i],
            &riv_matl_tbl->ksath[i], &riv_matl_tbl->ksatv[i],
            &riv_matl_tbl->bedthick[i]);
        if (match != 6 || i != index - 1)
        {
            printf ("Cannot read information for the %dth material!\n",
                i + 1);
            printf (".riv file format error!\n");
            exit (1);
        }
    }

    /*
     * Read river initial condition information
     */
    FindLine (riv_file, "IC");
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%d", &riv_ic_tbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of river materials!\n");
        printf (".riv file format error!\n");
        exit (1);
    }

    /* Allocate */
    riv_ic_tbl->stage =
        (double *)malloc (riv_ic_tbl->number * sizeof (double));

    for (i = 0; i < riv_ic_tbl->number; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf", &index, &riv_ic_tbl->stage[i]);
        if (match != 2 || i != index - 1)
        {
            printf
                ("Cannot read information for the %dth initial condition!\n",
                i + 1);
            printf (".riv file format error!\n");
            exit (1);
        }
    }

    ic->rivgw = (double *)malloc (riv_att_tbl->number * sizeof (double));
    ic->stage = (double *)malloc (riv_att_tbl->number * sizeof (double));

    for (i = 0; i < riv_att_tbl->number; i++)
    {
        ic->stage[i] = riv_ic_tbl->stage[riv_att_tbl->ic[i] - 1];
    }

    /*
     * Read river boundary condition block
     */
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%*s %d", &forcing->nts[RIV_TS]);
    if (match != 1)
    {
        printf ("Cannot read number of river boundary conditions!\n");
        printf (".riv file format error!\n");
        exit (1);
    }

    if (forcing->nts[RIV_TS] > 0)
    {
        forcing->ts[RIV_TS] =
            (ts_struct *)malloc (forcing->nts[RIV_TS] * sizeof (ts_struct));
    }

    for (i = 0; i < forcing->nts[RIV_TS]; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%*s %d", &index);
        if (match != 1 || i != index - 1)
        {
            printf
                ("Cannot read information of the %dth river boudnary condition!\n",
                i);
            printf (".riv file format error!\n");
            exit (1);
        }
        NextLine (riv_file, cmdstr);
        NextLine (riv_file, cmdstr);
        forcing->ts[RIV_TS][i].length =
            CountLine (riv_file, 2, "RIV_TS", "RES");
    }

    FindLine (riv_file, "BC");
    for (i = 0; i < forcing->nts[RIV_TS]; i++)
    {
        NextLine (riv_file, cmdstr);
        NextLine (riv_file, cmdstr);
        NextLine (riv_file, cmdstr);

        forcing->ts[RIV_TS][i].data =
            (double **)malloc ((forcing->ts[RIV_TS][i].length) *
            sizeof (double *));
        forcing->ts[RIV_TS][i].ftime =
            (int *)malloc ((forcing->ts[RIV_TS][i].length) * sizeof (int));
        for (j = 0; j < forcing->ts[RIV_TS][i].length; j++)
        {
            forcing->ts[RIV_TS][i].data[j] =
                (double *)malloc (sizeof (double));
            NextLine (riv_file, cmdstr);
            ReadTS (cmdstr, &forcing->ts[RIV_TS][i].ftime[j],
                &forcing->ts[RIV_TS][i].data[j][0], 1);
        }
    }

    /* Read Reservoir information */
    /* Empty */

    fclose (riv_file);
}

void ReadMesh (char *project, mesh_tbl_struct *mesh_tbl)
{
    FILE           *mesh_file;  /* Pointer to .mesh file */
    char            fn[MAXSTRING];
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    /*
     * Open .mesh input file
     */
    sprintf (fn, "input/%s/%s.mesh", project, project);
    mesh_file = fopen (fn, "r");
    CheckFile (mesh_file, fn);

    /*
     * Read element mesh block
     */
    NextLine (mesh_file, cmdstr);
    match = sscanf (cmdstr, "%d", &mesh_tbl->numele);
    if (match != 1)
    {
        printf ("Cannot read number of elements!\n");
        printf (".mesh file format error!\n");
        exit (1);
    }

    mesh_tbl->node = (int **)malloc (mesh_tbl->numele * sizeof (int *));
    mesh_tbl->nabr = (int **)malloc (mesh_tbl->numele * sizeof (int *));

    for (i = 0; i < mesh_tbl->numele; i++)
    {
        mesh_tbl->node[i] = (int *)malloc (3 * sizeof (int));
        mesh_tbl->nabr[i] = (int *)malloc (3 * sizeof (int));

        NextLine (mesh_file, cmdstr);
        match = sscanf (cmdstr, "%d %d %d %d %d %d %d",
            &index,
            &mesh_tbl->node[i][0], &mesh_tbl->node[i][1],
            &mesh_tbl->node[i][2], &mesh_tbl->nabr[i][0],
            &mesh_tbl->nabr[i][1], &mesh_tbl->nabr[i][2]);
        if (match != 7 || i != index - 1)
        {
            printf ("Cannot read information of the %dth element!\n", i + 1);
            printf (".mesh file format error!\n");
            exit (1);
        }
    }

    /*
     * Read node block
     */
    NextLine (mesh_file, cmdstr);
    match = sscanf (cmdstr, "%d", &mesh_tbl->numnode);
    if (match != 1)
    {
        printf ("Cannot read number of nodes!\n");
        printf (".mesh file format error!\n");
        exit (1);
    }

    mesh_tbl->x = (double *)malloc (mesh_tbl->numnode * sizeof (double));
    mesh_tbl->y = (double *)malloc (mesh_tbl->numnode * sizeof (double));
    mesh_tbl->zmin = (double *)malloc (mesh_tbl->numnode * sizeof (double));
    mesh_tbl->zmax = (double *)malloc (mesh_tbl->numnode * sizeof (double));

    for (i = 0; i < mesh_tbl->numnode; i++)
    {
        NextLine (mesh_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf",
            &index,
            &mesh_tbl->x[i], &mesh_tbl->y[i],
            &mesh_tbl->zmin[i], &mesh_tbl->zmax[i]);
        if (match != 5 || i != index - 1)
        {
            printf ("Cannot read information of the %dth node!\n", i + 1);
            printf (".mesh file format error!\n");
            exit (1);
        }
    }

    /* finish reading mesh_files */
    fclose (mesh_file);
}

void ReadAtt (char *project, attrib_tbl_struct *attrib_tbl, ic_struct *ic,
    int numele)
{
    char            fn[MAXSTRING];
    int             i;
    FILE           *att_file;   /* Pointer to .att file */
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    sprintf (fn, "input/%s/%s.att", project, project);
    att_file = fopen (fn, "r");
    CheckFile (att_file, fn);

    /* start reading att_file */
    ic->intcp = (double *)malloc (numele * sizeof (double));
    ic->snow = (double *)malloc (numele * sizeof (double));
    ic->surf = (double *)malloc (numele * sizeof (double));
    ic->unsat = (double *)malloc (numele * sizeof (double));
    ic->gw = (double *)malloc (numele * sizeof (double));

    attrib_tbl->soil = (int *)malloc (numele * sizeof (int));
    attrib_tbl->geol = (int *)malloc (numele * sizeof (int));
    attrib_tbl->lc = (int *)malloc (numele * sizeof (int));
    attrib_tbl->bc = (int **)malloc (numele * sizeof (int *));
    attrib_tbl->meteo = (int *)malloc (numele * sizeof (int));
    attrib_tbl->lai = (int *)malloc (numele * sizeof (int));
    attrib_tbl->source = (int *)malloc (numele * sizeof (int));
    attrib_tbl->macropore = (int *)malloc (numele * sizeof (int));

    NextLine (att_file, cmdstr);
    for (i = 0; i < numele; i++)
    {
        attrib_tbl->bc[i] = (int *)malloc (3 * sizeof (int));

        NextLine (att_file, cmdstr);
        match =
            sscanf (cmdstr,
            "%d %d %d %d %lf %lf %lf %lf %lf %d %d %d %d %d %d %d", &index,
            &attrib_tbl->soil[i], &attrib_tbl->geol[i], &attrib_tbl->lc[i],
            &ic->intcp[i], &ic->snow[i], &ic->surf[i], &ic->unsat[i],
            &ic->gw[i], &attrib_tbl->meteo[i], &attrib_tbl->lai[i],
            &attrib_tbl->source[i], &attrib_tbl->bc[i][0],
            &attrib_tbl->bc[i][1], &attrib_tbl->bc[i][2],
            &attrib_tbl->macropore[i]);
        if (match != 16)
        {
            printf ("Cannot read information of the %dth element!\n", i + 1);
            printf (".att file format error!\n");
            exit (1);
        }
    }

    /* finish reading att_files */
    fclose (att_file);
}

void ReadSoil (char *project, soil_tbl_struct *soil_tbl)
{
    FILE           *soil_file;  /* Pointer to .soil file */
    char            fn[MAXSTRING];
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    sprintf (fn, "input/%s/%s.soil", project, project);
    soil_file = fopen (fn, "r");
    CheckFile (soil_file, fn);

    /* start reading soil_file */
    NextLine (soil_file, cmdstr);
    match = sscanf (cmdstr, "%d", &soil_tbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of soil types!\n");
        printf (".soil file format error!\n");
        exit (1);
    }

    soil_tbl->ksatv = (double *)malloc (soil_tbl->number * sizeof (double));
    soil_tbl->thetas = (double *)malloc (soil_tbl->number * sizeof (double));
    soil_tbl->thetar = (double *)malloc (soil_tbl->number * sizeof (double));
    soil_tbl->qtz = (double *)malloc (soil_tbl->number * sizeof (double));
    soil_tbl->alpha = (double *)malloc (soil_tbl->number * sizeof (double));
    soil_tbl->beta = (double *)malloc (soil_tbl->number * sizeof (double));
    soil_tbl->areafh = (double *)malloc (soil_tbl->number * sizeof (double));
    soil_tbl->kmacv = (double *)malloc (soil_tbl->number * sizeof (double));
    soil_tbl->dinf = (double *)malloc (soil_tbl->number * sizeof (double));

    for (i = 0; i < soil_tbl->number; i++)
    {
        NextLine (soil_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &index,
            &soil_tbl->ksatv[i],
            &soil_tbl->thetas[i], &soil_tbl->thetar[i],
            &soil_tbl->dinf[i],
            &soil_tbl->alpha[i], &soil_tbl->beta[i],
            &soil_tbl->areafh[i], &soil_tbl->kmacv[i], &soil_tbl->qtz[i]);
        if (match != 10 || i != index - 1)
        {
            printf ("Cannot read information of the %dth soil type!\n",
                i + 1);
            printf (".soil file format error!\n");
            exit (1);
        }
    }

    fclose (soil_file);
}

void ReadGeol (char *project, geol_tbl_struct *geol_tbl)
{
    FILE           *geol_file;  /* Pointer to .geol file */
    char            fn[MAXSTRING];
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    sprintf (fn, "input/%s/%s.geol", project, project);
    geol_file = fopen (fn, "r");
    CheckFile (geol_file, fn);

    /* start reading geol_file */
    NextLine (geol_file, cmdstr);
    match = sscanf (cmdstr, "%d", &geol_tbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of geology types!\n");
        printf (".geol file format error!\n");
        exit (1);
    }

    geol_tbl->ksath = (double *)malloc (geol_tbl->number * sizeof (double));
    geol_tbl->ksatv = (double *)malloc (geol_tbl->number * sizeof (double));
    geol_tbl->thetas = (double *)malloc (geol_tbl->number * sizeof (double));
    geol_tbl->thetar = (double *)malloc (geol_tbl->number * sizeof (double));
    geol_tbl->alpha = (double *)malloc (geol_tbl->number * sizeof (double));
    geol_tbl->beta = (double *)malloc (geol_tbl->number * sizeof (double));
    geol_tbl->areafv = (double *)malloc (geol_tbl->number * sizeof (double));
    geol_tbl->kmach = (double *)malloc (geol_tbl->number * sizeof (double));
    geol_tbl->dmac = (double *)malloc (geol_tbl->number * sizeof (double));

    for (i = 0; i < geol_tbl->number; i++)
    {
        NextLine (geol_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &index,
            &geol_tbl->ksath[i], &geol_tbl->ksatv[i],
            &geol_tbl->thetas[i], &geol_tbl->thetar[i],
            &geol_tbl->alpha[i], &geol_tbl->beta[i],
            &geol_tbl->areafv[i], &geol_tbl->kmach[i], &geol_tbl->dmac[i]);
        if (match != 10 || i != index - 1)
        {
            printf ("Cannot read information of the %dth geology type!\n",
                i + 1);
            printf (".geol file format error!\n");
            exit (1);
        }
    }

    fclose (geol_file);
}

void ReadLC (lc_tbl_struct *lc_tbl)
{
    FILE           *lc_file;    /* Pointer to .lc file */
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    lc_file = fopen ("input/vegprmt.tbl", "r");
    CheckFile (lc_file, "input/vegprmt.tbl");

    /* start reading land cover file */
    NextLine (lc_file, cmdstr);
    match = sscanf (cmdstr, "%d", &lc_tbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of landcover types!\n");
        printf ("Land cover file format error!\n");
        exit (1);
    }

    lc_tbl->laimax = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->laimin = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->vegfrac = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->albedomin = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->albedomax = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->emissmin = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->emissmax = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->z0min = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->z0max = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->hs = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->snup = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->rgl = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->rsmin = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->rough = (double *)malloc (lc_tbl->number * sizeof (double));
    lc_tbl->rzd = (double *)malloc (lc_tbl->number * sizeof (double));

    for (i = 0; i < lc_tbl->number; i++)
    {
        NextLine (lc_file, cmdstr);
        match =
            sscanf (cmdstr,
            "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &index, &lc_tbl->vegfrac[i], &lc_tbl->rzd[i], &lc_tbl->rsmin[i],
            &lc_tbl->rgl[i], &lc_tbl->hs[i], &lc_tbl->snup[i],
            &lc_tbl->laimin[i], &lc_tbl->laimax[i], &lc_tbl->emissmin[i],
            &lc_tbl->emissmax[i], &lc_tbl->albedomin[i],
            &lc_tbl->albedomax[i], &lc_tbl->z0min[i], &lc_tbl->z0max[i],
            &lc_tbl->rough[i]);
        if (match != 16 || i != index - 1)
        {
            printf ("Cannot read information of the %dth landcover type!\n",
                i + 1);
            printf ("Landcover file format error!\n");
            exit (1);
        }
    }

    NextLine (lc_file, cmdstr);
    ReadKeywordDouble (cmdstr, "TOPT_DATA", &lc_tbl->topt);

    NextLine (lc_file, cmdstr);
    ReadKeywordDouble (cmdstr, "CFACTR_DATA", &lc_tbl->cfactr);

    NextLine (lc_file, cmdstr);
    ReadKeywordDouble (cmdstr, "RSMAX_DATA", &lc_tbl->rsmax);

    NextLine (lc_file, cmdstr);
    ReadKeywordInt (cmdstr, "BARE", &lc_tbl->bare);

    NextLine (lc_file, cmdstr);
    ReadKeywordInt (cmdstr, "NATURAL", &lc_tbl->natural);

    fclose (lc_file);
}

void ReadForc (char *project, forcing_ts_struct *forcing)
{
    char            fn[MAXSTRING];
    FILE           *forc_file;  /* Pointer to .forc file */
    char            cmdstr[MAXSTRING];
    int             i, j;
    int             match;
    int             index;

    sprintf (fn, "input/%s/%s.forc", project, project);
    forc_file = fopen (fn, "r");
    CheckFile (forc_file, fn);

    FindLine (forc_file, "BOF");
    NextLine (forc_file, cmdstr);
    match = sscanf (cmdstr, "%*s %d", &forcing->nts[METEO_TS]);
    if (match != 1)
    {
        printf
            ("Cannot read number of meteorological forcing time series!\n");
        printf (".forc file format error!\n");
        exit (1);
    }

    if (forcing->nts[METEO_TS] > 0)
    {
        forcing->ts[METEO_TS] =
            (ts_struct *)malloc (forcing->nts[METEO_TS] * sizeof (ts_struct));
        forcing->zlvl_wind =
            (double *)malloc (forcing->nts[METEO_TS] * sizeof (double));
    }

    for (i = 0; i < forcing->nts[METEO_TS]; i++)
    {
        NextLine (forc_file, cmdstr);
        match =
            sscanf (cmdstr, "%*s %d %*s %lf", &index, &forcing->zlvl_wind[i]);
        if (match != 2 || i != index - 1)
        {
            printf ("Cannot read information of the %dth forcing series!\n",
                i);
            printf (".forc file format error!\n");
            exit (1);
        }
        /* Skip header lines */
        NextLine (forc_file, cmdstr);
        NextLine (forc_file, cmdstr);
        forcing->ts[METEO_TS][i].length =
            CountLine (forc_file, 1, "METEO_TS");
    }

    /* Rewind and read */
    FindLine (forc_file, "NUM_METEO_TS");
    for (i = 0; i < forcing->nts[METEO_TS]; i++)
    {
        /* Skip header lines */
        NextLine (forc_file, cmdstr);
        NextLine (forc_file, cmdstr);
        NextLine (forc_file, cmdstr);

        forcing->ts[METEO_TS][i].ftime =
            (int *)malloc (forcing->ts[METEO_TS][i].length * sizeof (int));
        forcing->ts[METEO_TS][i].data =
            (double **)malloc (forcing->ts[METEO_TS][i].length *
            sizeof (double *));
        for (j = 0; j < forcing->ts[METEO_TS][i].length; j++)
        {
            forcing->ts[METEO_TS][i].data[j] =
                (double *)malloc (NUM_METEO_TS * sizeof (double));
            NextLine (forc_file, cmdstr);
            ReadTS (cmdstr, &forcing->ts[METEO_TS][i].ftime[j],
                &forcing->ts[METEO_TS][i].data[j][0], NUM_METEO_TS);
        }
    }

    fclose (forc_file);

}

void ReadLAI (char *project, forcing_ts_struct *forcing, int numele,
    const attrib_tbl_struct *attrib_tbl)
{
    char            fn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    int             read_lai = 0;
    FILE           *lai_file;
    int             i, j;
    int             index;
    int             match;

    for (i = 0; i < numele; i++)
    {
        if (attrib_tbl->lai[i] > 0)
        {
            read_lai = 1;
            break;
        }
    }

    forcing->nts[LAI_TS] = 0;

    if (read_lai)
    {
        sprintf (fn, "input/%s/%s.lai", project, project);
        lai_file = fopen (fn, "r");
        CheckFile (lai_file, fn);

        /* start reading lai_file */
        FindLine (lai_file, "BOF");
        NextLine (lai_file, cmdstr);
        match = sscanf (cmdstr, "%*s %d", &forcing->nts[LAI_TS]);
        if (match != 1)
        {
            printf ("Cannot read number of LAI time series!\n");
            printf (".lai file format error!\n");
            exit (1);
        }

        if (forcing->nts[LAI_TS] > 0)
        {
            forcing->ts[LAI_TS] =
                (ts_struct *)malloc (forcing->nts[LAI_TS] * sizeof (ts_struct));
        }

        for (i = 0; i < forcing->nts[LAI_TS]; i++)
        {
            NextLine (lai_file, cmdstr);
            match = sscanf (cmdstr, "%*s %d", &index);
            if (match != 1 || i != index - 1)
            {
                printf ("Cannot read information of the %dth LAI series!\n",
                    i);
                printf (".lai file format error!\n");
                exit (1);
            }
            /* Skip header lines */
            NextLine (lai_file, cmdstr);
            NextLine (lai_file, cmdstr);
            forcing->ts[LAI_TS][i].length = CountLine (lai_file, 1, "LAI_TS");
        }

        /* Rewind and read */
        FindLine (lai_file, "NUM_LAI_TS");
        for (i = 0; i < forcing->nts[LAI_TS]; i++)
        {
            /* Skip header lines */
            NextLine (lai_file, cmdstr);
            NextLine (lai_file, cmdstr);
            NextLine (lai_file, cmdstr);

            forcing->ts[LAI_TS][i].ftime =
                (int *)malloc (forcing->ts[LAI_TS][i].length * sizeof (int));
            forcing->ts[LAI_TS][i].data =
                (double **)malloc (forcing->ts[LAI_TS][i].length *
                sizeof (double *));
            for (j = 0; j < forcing->ts[LAI_TS][i].length; j++)
            {
                forcing->ts[LAI_TS][i].data[j] =
                    (double *)malloc (sizeof (double));
                NextLine (lai_file, cmdstr);
                ReadTS (cmdstr, &forcing->ts[LAI_TS][i].ftime[j],
                    &forcing->ts[LAI_TS][i].data[j][0], 1);
            }
        }
        fclose (lai_file);
    }
}

void ReadIbc (char *project, forcing_ts_struct *forcing)
{
    char            fn[MAXSTRING];
    int             i, j;
    FILE           *ibc_file;   /* Pointer to .ibc file */
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    sprintf (fn, "input/%s/%s.ibc", project, project);
    ibc_file = fopen (fn, "r");
    CheckFile (ibc_file, fn);

    /*
     * Start reading ibc_file 
     */
    FindLine (ibc_file, "BOF");
    NextLine (ibc_file, cmdstr);
    match = sscanf (cmdstr, "%*s %d", &forcing->nts[BC_TS]);
    if (match != 1)
    {
        printf ("Cannot read number of boundary condition time series!\n");
        printf (".ibc file format error!\n");
        exit (1);
    }

    if (forcing->nts[BC_TS] > 0)
    {
        forcing->ts[BC_TS] =
            (ts_struct *)malloc (forcing->nts[BC_TS] * sizeof (ts_struct));
    }

    for (i = 0; i < forcing->nts[BC_TS]; i++)
    {
        NextLine (ibc_file, cmdstr);
        match = sscanf (cmdstr, "%*s %d", &index);
        if (match != 1 || i != index - 1)
        {
            printf
                ("Cannot read information of the %dth boundary condition series!\n",
                i);
            printf (".ibc file format error!\n");
            exit (1);
        }
        /* Skip header lines */
        NextLine (ibc_file, cmdstr);
        NextLine (ibc_file, cmdstr);

        forcing->ts[BC_TS][i].length = CountLine (ibc_file, 1, "BC_TS");
    }

    /* Rewind and read */
    FindLine (ibc_file, "NUM_BC_TS");
    for (i = 0; i < forcing->nts[BC_TS]; i++)
    {
        /* Skip header lines */
        NextLine (ibc_file, cmdstr);
        NextLine (ibc_file, cmdstr);
        NextLine (ibc_file, cmdstr);

        forcing->ts[BC_TS][i].ftime =
            (int *)malloc (forcing->ts[BC_TS][i].length * sizeof (int));
        forcing->ts[BC_TS][i].data =
            (double **)malloc (forcing->ts[BC_TS][i].length *
            sizeof (double *));
        for (j = 0; j < forcing->ts[BC_TS][i].length; j++)
        {
            forcing->ts[BC_TS][i].data[j] =
                (double *)malloc (sizeof (double));
            NextLine (ibc_file, cmdstr);
            ReadTS (cmdstr, &forcing->ts[BC_TS][i].ftime[j],
                &forcing->ts[BC_TS][i].data[j][0], 1);
        }
    }

    fclose (ibc_file);
}

void ReadPara (char *project, ctrl_struct *ctrl)
{
    char            fn[MAXSTRING];
    FILE           *para_file;  /* Pointer to .para file */
    char            cmdstr[MAXSTRING];
    int             i;

    for (i = 0; i < NUM_PRINT; i++)
    {
        ctrl->prtvrbl[i] = 0;
    }

    sprintf (fn, "input/%s/%s.para", project, project);
    para_file = fopen (fn, "r");
    CheckFile (para_file, fn);

    /* start reading para_file */
    /* Read through parameter file to find parameters */
    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "INIT_MODE", &ctrl->init_type);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "ASCII_OUTPUT", &ctrl->ascii);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "WRITE_IC", &ctrl->write_ic);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "UNSAT_MODE", &ctrl->unsat_mode);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "SURF_MODE", &ctrl->surf_mode);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIV_MODE", &ctrl->riv_mode);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "SOLVER", &ctrl->solver);

    NextLine (para_file, cmdstr);
    ReadKeywordDouble (cmdstr, "ABSTOL", &ctrl->abstol);

    NextLine (para_file, cmdstr);
    ReadKeywordDouble (cmdstr, "RELTOL", &ctrl->reltol);

    NextLine (para_file, cmdstr);
    ReadKeywordDouble (cmdstr, "INIT_SOLVER_STEP", &ctrl->initstep);

    NextLine (para_file, cmdstr);
    ReadKeywordDouble (cmdstr, "MAX_SOLVER_STEP", &ctrl->maxstep);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "LSM_STEP", &ctrl->etstep);

    NextLine (para_file, cmdstr);
    ReadKeywordTime (cmdstr, "START", &ctrl->starttime);

    NextLine (para_file, cmdstr);
    ReadKeywordTime (cmdstr, "END", &ctrl->endtime);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "MODEL_STEPSIZE", &ctrl->stepsize);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "GW", &ctrl->prtvrbl[GW_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "SURF", &ctrl->prtvrbl[SURF_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "SNOW", &ctrl->prtvrbl[SNOW_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIVSTG", &ctrl->prtvrbl[RIVSTG_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "INFIL", &ctrl->prtvrbl[INFIL_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RECHARGE", &ctrl->prtvrbl[RECHARGE_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "CMC", &ctrl->prtvrbl[CMC_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "UNSAT", &ctrl->prtvrbl[UNSAT_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "EC", &ctrl->prtvrbl[EC_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "ETT", &ctrl->prtvrbl[ETT_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "EDIR", &ctrl->prtvrbl[EDIR_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIVFLX0", &ctrl->prtvrbl[RIVFLX0_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIVFLX1", &ctrl->prtvrbl[RIVFLX1_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIVFLX2", &ctrl->prtvrbl[RIVFLX2_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIVFLX3", &ctrl->prtvrbl[RIVFLX3_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIVFLX4", &ctrl->prtvrbl[RIVFLX4_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIVFLX5", &ctrl->prtvrbl[RIVFLX5_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIVFLX6", &ctrl->prtvrbl[RIVFLX6_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIVFLX7", &ctrl->prtvrbl[RIVFLX7_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIVFLX8", &ctrl->prtvrbl[RIVFLX8_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIVFLX9", &ctrl->prtvrbl[RIVFLX9_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "RIVFLX10", &ctrl->prtvrbl[RIVFLX10_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "SUBFLX", &ctrl->prtvrbl[SUBFLX_CTRL]);

    NextLine (para_file, cmdstr);
    ReadKeywordInt (cmdstr, "TOTALFLX", &ctrl->prtvrbl[TOTALFLX_CTRL]);

    fclose (para_file);

    if (ctrl->etstep < ctrl->stepsize || ctrl->etstep % ctrl->stepsize > 0)
    {
        printf
            ("ERROR: LSM (ET) step size should be an integral multiple of model step size!\n");
        exit (1);
    }
}

void ReadCalib (char *project, char *simulation, calib_struct *cal)
{
    char            fn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    FILE           *global_calib;       /* Pointer to .calib file */

    sprintf (fn, "input/%s/%s.calib", project, simulation);
    global_calib = fopen (fn, "r");
    CheckFile (global_calib, fn);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "KSATH", &cal->ksath);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "KSATV", &cal->ksatv);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "KINF", &cal->kinfv);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "KMACSATH", &cal->kmach);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "KMACSATV", &cal->kmacv);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "DINF", &cal->dinf);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "DROOT", &cal->rzd);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "DMAC", &cal->dmac);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "POROSITY", &cal->porosity);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "ALPHA", &cal->alpha);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "BETA", &cal->beta);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "MACVF", &cal->areafv);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "MACHF", &cal->areafh);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "VEGFRAC", &cal->vegfrac);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "ALBEDO", &cal->albedo);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "ROUGH", &cal->rough);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "PRCP", &cal->prcp);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "SFCTMP", &cal->sfctmp);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "EC", &cal->et[0]);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "ETT", &cal->et[1]);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "EDIR", &cal->et[2]);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "ROUGH_RIV", &cal->rivrough);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "KRIVH", &cal->rivksath);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "KRIVV", &cal->rivksatv);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "BEDTHCK", &cal->rivbedthick);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "RIV_DPTH", &cal->rivdepth);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "RIV_WDTH", &cal->rivshpcoeff);

#ifdef _RT_
    CS->Cal.PCO2 = 1.0;
    CS->Cal.Keq = 1.0;
    CS->Cal.Site_den = 1.0;
    CS->Cal.SSA = 1.0;
    CS->Cal.Prep_conc = 1.0;
#endif

#ifdef _NOAH_
    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "DRIP", &cal->drip);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "CMCMAX", &cal->intcp);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "RS", &cal->rsmin);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "CZIL", &cal->czil);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "FXEXP", &cal->fxexp);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "CFACTR", &cal->cfactr);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "RGL", &cal->rgl);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "HS", &cal->hs);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "REFSMC", &cal->thetaref);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "WLTSMC", &cal->thetaw);
#endif
    /* finish reading calib file */
    fclose (global_calib);
}

void ReadInit (char *project, char *simulation, ic_struct *ic, int numele,
    int numriv)
{
    char            fn[MAXSTRING];
    FILE           *init_file;
    int             i;

    sprintf (fn, "input/%s/%s.init", project, simulation);
    init_file = fopen (fn, "rb");
    CheckFile (init_file, fn);

    for (i = 0; i < numele; i++)
    {
        fread (&ic->intcp[i], sizeof (double), 1, init_file);
        fread (&ic->snow[i], sizeof (double), 1, init_file);
        fread (&ic->surf[i], sizeof (double), 1, init_file);
        fread (&ic->unsat[i], sizeof (double), 1, init_file);
        fread (&ic->gw[i], sizeof (double), 1, init_file);
    }
    for (i = 0; i < numriv; i++)
    {
        fread (&ic->stage[i], sizeof (double), 1, init_file);
        fread (&ic->rivgw[i], sizeof (double), 1, init_file);
    }
}

void FreeData (pihm_struct pihm)
{
    int             i, j, k;

    /* Free river input structure */
    free (pihm->riv_att_tbl.fromnode);
    free (pihm->riv_att_tbl.tonode);
    free (pihm->riv_att_tbl.down);
    free (pihm->riv_att_tbl.leftele);
    free (pihm->riv_att_tbl.rightele);
    free (pihm->riv_att_tbl.shp);
    free (pihm->riv_att_tbl.matl);
    free (pihm->riv_att_tbl.ic);
    free (pihm->riv_att_tbl.bc);
    free (pihm->riv_att_tbl.rsvr);

    free (pihm->riv_shp_tbl.depth);
    free (pihm->riv_shp_tbl.intrpl_ord);
    free (pihm->riv_shp_tbl.coeff);

    free (pihm->riv_matl_tbl.rough);
    free (pihm->riv_matl_tbl.cwr);
    free (pihm->riv_matl_tbl.ksath);
    free (pihm->riv_matl_tbl.ksatv);
    free (pihm->riv_matl_tbl.bedthick);

    free (pihm->riv_ic_tbl.stage);

    /* Free mesh input structure */
    for (i = 0; i < pihm->mesh_tbl.numele; i++)
    {
        free (pihm->mesh_tbl.node[i]);
        free (pihm->mesh_tbl.nabr[i]);
    }
    free (pihm->mesh_tbl.node);
    free (pihm->mesh_tbl.nabr);
    free (pihm->mesh_tbl.x);
    free (pihm->mesh_tbl.y);
    free (pihm->mesh_tbl.zmin);
    free (pihm->mesh_tbl.zmax);

    /* Free attribute input structure */
    for (i = 0; i < pihm->mesh_tbl.numele; i++)
    {
        free (pihm->attrib_tbl.bc[i]);
    }
    free (pihm->attrib_tbl.soil);
    free (pihm->attrib_tbl.geol);
    free (pihm->attrib_tbl.lc);
    free (pihm->attrib_tbl.bc);
    free (pihm->attrib_tbl.meteo);
    free (pihm->attrib_tbl.lai);
    free (pihm->attrib_tbl.source);
    free (pihm->attrib_tbl.macropore);

    /* Free soil input structure */
    free (pihm->soil_tbl.ksatv);
    free (pihm->soil_tbl.thetas);
    free (pihm->soil_tbl.thetar);
    free (pihm->soil_tbl.qtz);
    free (pihm->soil_tbl.alpha);
    free (pihm->soil_tbl.beta);
    free (pihm->soil_tbl.areafh);
    free (pihm->soil_tbl.kmacv);
    free (pihm->soil_tbl.dinf);

    /* Free geol input structure */
    free (pihm->geol_tbl.ksath);
    free (pihm->geol_tbl.ksatv);
    free (pihm->geol_tbl.thetas);
    free (pihm->geol_tbl.thetar);
    free (pihm->geol_tbl.alpha);
    free (pihm->geol_tbl.beta);
    free (pihm->geol_tbl.areafv);
    free (pihm->geol_tbl.kmach);
    free (pihm->geol_tbl.dmac);

    /* Free initial condition */
    free (pihm->ic.intcp);
    free (pihm->ic.snow);
    free (pihm->ic.surf);
    free (pihm->ic.unsat);
    free (pihm->ic.gw);
    free (pihm->ic.rivgw);
    free (pihm->ic.stage);

    /* Free landcover input structure */
    free (pihm->lc_tbl.laimax);
    free (pihm->lc_tbl.laimin);
    free (pihm->lc_tbl.vegfrac);
    free (pihm->lc_tbl.albedomin);
    free (pihm->lc_tbl.albedomax);
    free (pihm->lc_tbl.emissmin);
    free (pihm->lc_tbl.emissmax);
    free (pihm->lc_tbl.z0min);
    free (pihm->lc_tbl.z0max);
    free (pihm->lc_tbl.hs);
    free (pihm->lc_tbl.snup);
    free (pihm->lc_tbl.rgl);
    free (pihm->lc_tbl.rsmin);
    free (pihm->lc_tbl.rough);
    free (pihm->lc_tbl.rzd);

    /* Free forcing input structure */
    for (k = 0; k < NUM_TS; k++)
    {
        if (pihm->forcing.nts[k] > 0)
        {
            for (i = 0; i < pihm->forcing.nts[k]; i++)
            {
                for (j = 0; j < pihm->forcing.ts[k][i].length; j++)
                {
                    free (pihm->forcing.ts[k][i].data[j]);
                }
                free (pihm->forcing.ts[k][i].ftime);
                free (pihm->forcing.ts[k][i].data);
            }
        free (pihm->forcing.ts[k]);
        }
    }

    if (pihm->forcing.nts[BC_TS] > 0)
    {
        free (pihm->forcing.bc);
    }

    if (pihm->forcing.nts[METEO_TS] > 0)
    {
        for (i = 0; i < NUM_METEO_TS; i++)
        {
            free (pihm->forcing.meteo[i]);
        }
        free (pihm->forcing.zlvl_wind);
    }
    if (pihm->forcing.nts[LAI_TS] > 0)
    {
        free (pihm->forcing.lai);
    }
    if (pihm->forcing.nts[RIV_TS] > 0)
    {
        free (pihm->forcing.riverbc);
    }
    if (pihm->forcing.nts[SS_TS] > 0)
    {
        free (pihm->forcing.source);
    }

    free (pihm->ctrl.tout);

    for (i = 0; i < pihm->ctrl.nprint; i++)
    {
        free (pihm->prtctrl[i].vrbl);
        free (pihm->prtctrl[i].buffer);
    }

    free (pihm->elem);
    free (pihm->riv);
}
