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
    if (verbose_mode)
    {
        printf ("\nRead input files:\n");
    }

    /*
     * Set file names of the input files
     */
    sprintf (pihm->filename.riv, "input/%s/%s.riv", project, project);
    sprintf (pihm->filename.mesh, "input/%s/%s.mesh", project, project);
    sprintf (pihm->filename.att, "input/%s/%s.att", project, project);
    sprintf (pihm->filename.soil, "input/%s/%s.soil", project, project);
    sprintf (pihm->filename.geol, "input/%s/%s.geol", project, project);
    sprintf (pihm->filename.lc, "input/vegprmt.tbl");
    sprintf (pihm->filename.forc, "input/%s/%s.forc", project, project);
    sprintf (pihm->filename.lai, "input/%s/%s.lai", project, project);
    sprintf (pihm->filename.ibc, "input/%s/%s.ibc", project, project);
    sprintf (pihm->filename.para, "input/%s/%s.para", project, project);
    sprintf (pihm->filename.calib, "input/%s/%s.calib", project, simulation);
    sprintf (pihm->filename.init, "input/%s/%s.init", project, simulation);
#ifdef _NOAH_
    sprintf (pihm->filename.lsm, "input/%s/%s.lsm", project, project);
    sprintf (pihm->filename.rad, "input/%s/%s.rad", project, project);
    sprintf (pihm->filename.lsminit, "input/%s/%s.lsminit", project, simulation);
#endif

    ReadRiv (pihm->filename.riv, &pihm->rivtbl, &pihm->shptbl, &pihm->matltbl,
        &pihm->forc);
    pihm->numriv = pihm->rivtbl.number;

    ReadMesh (pihm->filename.mesh, &pihm->meshtbl);
    pihm->numele = pihm->meshtbl.numele;

    ReadAtt (pihm->filename.att, &pihm->atttbl, pihm->numele);

    ReadSoil (pihm->filename.soil, &pihm->soiltbl);

    ReadGeol (pihm->filename.geol, &pihm->geoltbl);

    ReadLC (pihm->filename.lc, &pihm->lctbl);

    ReadForc (pihm->filename.forc, &pihm->forc);

    ReadLAI (pihm->filename.lai, &pihm->forc, pihm->numele, &pihm->atttbl);

    /* Read source and sink */
    //ReadSS ();

    ReadIbc (pihm->filename.ibc, &pihm->forc);

    ReadPara (pihm->filename.para, &pihm->ctrl);

    ReadCalib (pihm->filename.calib, &pihm->cal);
#ifdef _NOAH_
    ReadLsm (pihm->filename.lsm, &pihm->latitude, &pihm->longitude,
        &pihm->ctrl, &pihm->noahtbl);

    if (pihm->ctrl.rad_mode == 1)
    {
        ReadRad (pihm->filename.rad, &pihm->forc);
    }
#endif
}

void ReadRiv (char *filename, rivtbl_struct *rivtbl, shptbl_struct *shptbl,
    matltbl_struct *matltbl, forc_struct *forc)
{
    int             i, j;
    FILE           *riv_file;   /* Pointer to .riv file */
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    /*
     * Open .riv input file
     */
    riv_file = fopen (filename, "r");
    CheckFile (riv_file, filename);

    /*
     * Read river segment block
     */
    /* Read number of river segments */
    FindLine (riv_file, "BOF");
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%d", &rivtbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of river segments!\n");
        printf (".riv file format error!\n");
        PihmExit (1);
    }

    /* Allocate */
    rivtbl->fromnode =
        (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->tonode = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->down = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->leftele = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->rightele =
        (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->shp = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->matl = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->bc = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->rsvr = (int *)malloc (rivtbl->number * sizeof (int));

    /* Read river segment information */
    for (i = 0; i < rivtbl->number; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%d %d %d %d %d %d %d %d %d %d",
            &index,
            &rivtbl->fromnode[i], &rivtbl->tonode[i],
            &rivtbl->down[i],
            &rivtbl->leftele[i], &rivtbl->rightele[i],
            &rivtbl->shp[i], &rivtbl->matl[i],
            &rivtbl->bc[i], &rivtbl->rsvr[i]);
        if (match != 10 || i != index - 1)
        {
            printf
                ("Cannot read river segment information for the %dth segment!\n",
                i + 1);
            printf (".riv file format error!\n");
            PihmExit (1);
        }
    }

    /*
     * Read river shape information
     */
    FindLine (riv_file, "SHAPE");
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%d", &shptbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of river shapes!\n");
        printf (".riv file format error!\n");
        PihmExit (1);
    }

    /* Allocate */
    shptbl->depth =
        (double *)malloc (shptbl->number * sizeof (double));
    shptbl->intrpl_ord =
        (int *)malloc (shptbl->number * sizeof (int));
    shptbl->coeff =
        (double *)malloc (shptbl->number * sizeof (double));

    for (i = 0; i < shptbl->number; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %d %lf",
            &index, &shptbl->depth[i],
            &shptbl->intrpl_ord[i], &shptbl->coeff[i]);
        if (match != 4 || i != index - 1)
        {
            printf
                ("Cannot read river shape information for the %dth shape!\n",
                i + 1);
            printf (".riv file format error!\n");
            PihmExit (1);
        }
    }

    /*
     * Read river material information
     */
    FindLine (riv_file, "MATERIAL");
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%d", &matltbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of river materials!\n");
        printf (".riv file format error!\n");
        PihmExit (1);
    }

    /* Allocate */
    matltbl->rough =
        (double *)malloc (matltbl->number * sizeof (double));
    matltbl->cwr =
        (double *)malloc (matltbl->number * sizeof (double));
    matltbl->ksath =
        (double *)malloc (matltbl->number * sizeof (double));
    matltbl->ksatv =
        (double *)malloc (matltbl->number * sizeof (double));
    matltbl->bedthick =
        (double *)malloc (matltbl->number * sizeof (double));

    for (i = 0; i < matltbl->number; i++)
    {
        NextLine (riv_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf",
            &index,
            &matltbl->rough[i], &matltbl->cwr[i],
            &matltbl->ksath[i], &matltbl->ksatv[i],
            &matltbl->bedthick[i]);
        if (match != 6 || i != index - 1)
        {
            printf ("Cannot read information for the %dth material!\n",
                i + 1);
            printf (".riv file format error!\n");
            PihmExit (1);
        }
    }

    /*
     * Read river boundary condition block
     */
    NextLine (riv_file, cmdstr);
    match = sscanf (cmdstr, "%*s %d", &forc->nriverbc);
    if (match != 1)
    {
        printf ("Cannot read number of river boundary conditions!\n");
        printf (".riv file format error!\n");
        PihmExit (1);
    }

    if (forc->nriverbc > 0)
    {
        forc->riverbc =
            (tsdata_struct *)malloc (forc->nriverbc * sizeof (tsdata_struct));

        for (i = 0; i < forc->nriverbc; i++)
        {
            NextLine (riv_file, cmdstr);
            match = sscanf (cmdstr, "%*s %d", &index);
            if (match != 1 || i != index - 1)
            {
                printf
                    ("Cannot read information of the %dth river boudnary condition!\n",
                    i);
                printf (".riv file format error!\n");
                PihmExit (1);
            }
            NextLine (riv_file, cmdstr);
            NextLine (riv_file, cmdstr);
            forc->riverbc[i].length =
                CountLine (riv_file, 2, "RIV_TS", "RES");
        }

        FindLine (riv_file, "BC");
        for (i = 0; i < forc->nriverbc; i++)
        {
            NextLine (riv_file, cmdstr);
            NextLine (riv_file, cmdstr);
            NextLine (riv_file, cmdstr);

            forc->riverbc[i].data =
                (double **)malloc ((forc->riverbc[i].length) *
                sizeof (double *));
            forc->riverbc[i].ftime =
                (int *)malloc ((forc->riverbc[i].length) * sizeof (int));
            for (j = 0; j < forc->riverbc[i].length; j++)
            {
                forc->riverbc[i].data[j] =
                    (double *)malloc (sizeof (double));
                NextLine (riv_file, cmdstr);
                ReadTS (cmdstr, &forc->riverbc[i].ftime[j],
                    &forc->riverbc[i].data[j][0], 1);
            }
        }
    }

    /* Read Reservoir information */
    /* Empty */

    fclose (riv_file);
}

void ReadMesh (char *filename, meshtbl_struct *meshtbl)
{
    FILE           *mesh_file;  /* Pointer to .mesh file */
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    /*
     * Open .mesh input file
     */
    mesh_file = fopen (filename, "r");
    CheckFile (mesh_file, filename);

    /*
     * Read element mesh block
     */
    NextLine (mesh_file, cmdstr);
    match = sscanf (cmdstr, "%d", &meshtbl->numele);
    if (match != 1)
    {
        printf ("Cannot read number of elements!\n");
        printf (".mesh file format error!\n");
        PihmExit (1);
    }

    meshtbl->node = (int **)malloc (meshtbl->numele * sizeof (int *));
    meshtbl->nabr = (int **)malloc (meshtbl->numele * sizeof (int *));

    for (i = 0; i < meshtbl->numele; i++)
    {
        meshtbl->node[i] = (int *)malloc (3 * sizeof (int));
        meshtbl->nabr[i] = (int *)malloc (3 * sizeof (int));

        NextLine (mesh_file, cmdstr);
        match = sscanf (cmdstr, "%d %d %d %d %d %d %d",
            &index,
            &meshtbl->node[i][0], &meshtbl->node[i][1],
            &meshtbl->node[i][2], &meshtbl->nabr[i][0],
            &meshtbl->nabr[i][1], &meshtbl->nabr[i][2]);
        if (match != 7 || i != index - 1)
        {
            printf ("Cannot read information of the %dth element!\n", i + 1);
            printf (".mesh file format error!\n");
            PihmExit (1);
        }
    }

    /*
     * Read node block
     */
    NextLine (mesh_file, cmdstr);
    match = sscanf (cmdstr, "%d", &meshtbl->numnode);
    if (match != 1)
    {
        printf ("Cannot read number of nodes!\n");
        printf (".mesh file format error!\n");
        PihmExit (1);
    }

    meshtbl->x = (double *)malloc (meshtbl->numnode * sizeof (double));
    meshtbl->y = (double *)malloc (meshtbl->numnode * sizeof (double));
    meshtbl->zmin = (double *)malloc (meshtbl->numnode * sizeof (double));
    meshtbl->zmax = (double *)malloc (meshtbl->numnode * sizeof (double));

    for (i = 0; i < meshtbl->numnode; i++)
    {
        NextLine (mesh_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf",
            &index,
            &meshtbl->x[i], &meshtbl->y[i],
            &meshtbl->zmin[i], &meshtbl->zmax[i]);
        if (match != 5 || i != index - 1)
        {
            printf ("Cannot read information of the %dth node!\n", i + 1);
            printf (".mesh file format error!\n");
            PihmExit (1);
        }
    }

    /* finish reading mesh_files */
    fclose (mesh_file);
}

void ReadAtt (char *filename, atttbl_struct *atttbl, int numele)
{
    int             i;
    FILE           *att_file;   /* Pointer to .att file */
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    att_file = fopen (filename, "r");
    CheckFile (att_file, filename);

    atttbl->soil = (int *)malloc (numele * sizeof (int));
    atttbl->geol = (int *)malloc (numele * sizeof (int));
    atttbl->lc = (int *)malloc (numele * sizeof (int));
    atttbl->bc = (int **)malloc (numele * sizeof (int *));
    for (i = 0; i < numele; i++)
    {
        atttbl->bc[i] = (int *)malloc (3 * sizeof (int));
    }
    atttbl->meteo = (int *)malloc (numele * sizeof (int));
    atttbl->lai = (int *)malloc (numele * sizeof (int));
    atttbl->source = (int *)malloc (numele * sizeof (int));
    atttbl->macropore = (int *)malloc (numele * sizeof (int));

    NextLine (att_file, cmdstr);
    for (i = 0; i < numele; i++)
    {
        NextLine (att_file, cmdstr);
        match = sscanf (cmdstr, "%d %d %d %d %d %d %d %d %d %d %d", &index,
            &atttbl->soil[i], &atttbl->geol[i], &atttbl->lc[i],
            &atttbl->meteo[i], &atttbl->lai[i], &atttbl->source[i],
            &atttbl->bc[i][0], &atttbl->bc[i][1], &atttbl->bc[i][2],
            &atttbl->macropore[i]);
        if (match != 11)
        {
            printf ("Cannot read information of the %dth element!\n", i + 1);
            printf (".att file format error!\n");
            PihmExit (1);
        }
    }

    /* finish reading att_files */
    fclose (att_file);
}

void ReadSoil (char *filename, soiltbl_struct *soiltbl)
{
    FILE           *soil_file;  /* Pointer to .soil file */
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    soil_file = fopen (filename, "r");
    CheckFile (soil_file, filename);

    /* start reading soil_file */
    NextLine (soil_file, cmdstr);
    match = sscanf (cmdstr, "%d", &soiltbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of soil types!\n");
        printf (".soil file format error!\n");
        PihmExit (1);
    }

    soiltbl->ksatv = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->smcmax = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->smcmin = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->qtz = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->alpha = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->beta = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->areafh = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->kmacv = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->dinf = (double *)malloc (soiltbl->number * sizeof (double));

    for (i = 0; i < soiltbl->number; i++)
    {
        NextLine (soil_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &index,
            &soiltbl->ksatv[i],
            &soiltbl->smcmax[i], &soiltbl->smcmin[i],
            &soiltbl->dinf[i],
            &soiltbl->alpha[i], &soiltbl->beta[i],
            &soiltbl->areafh[i], &soiltbl->kmacv[i], &soiltbl->qtz[i]);
        if (match != 10 || i != index - 1)
        {
            printf ("Cannot read information of the %dth soil type!\n",
                i + 1);
            printf (".soil file format error!\n");
            PihmExit (1);
        }
    }

    fclose (soil_file);
}

void ReadGeol (char *filename, geoltbl_struct *geoltbl)
{
    FILE           *geol_file;  /* Pointer to .geol file */
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    geol_file = fopen (filename, "r");
    CheckFile (geol_file, filename);

    /* start reading geol_file */
    NextLine (geol_file, cmdstr);
    match = sscanf (cmdstr, "%d", &geoltbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of geology types!\n");
        printf (".geol file format error!\n");
        PihmExit (1);
    }

    geoltbl->ksath = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->ksatv = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->smcmax = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->smcmin = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->alpha = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->beta = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->areafv = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->kmach = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->dmac = (double *)malloc (geoltbl->number * sizeof (double));

    for (i = 0; i < geoltbl->number; i++)
    {
        NextLine (geol_file, cmdstr);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &index,
            &geoltbl->ksath[i], &geoltbl->ksatv[i],
            &geoltbl->smcmax[i], &geoltbl->smcmin[i],
            &geoltbl->alpha[i], &geoltbl->beta[i],
            &geoltbl->areafv[i], &geoltbl->kmach[i], &geoltbl->dmac[i]);
        if (match != 10 || i != index - 1)
        {
            printf ("Cannot read information of the %dth geology type!\n",
                i + 1);
            printf (".geol file format error!\n");
            PihmExit (1);
        }
    }

    fclose (geol_file);
}

void ReadLC (char *filename, lctbl_struct *lctbl)
{
    FILE           *lc_file;    /* Pointer to .lc file */
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    lc_file = fopen (filename, "r");
    CheckFile (lc_file, filename);

    /* start reading land cover file */
    NextLine (lc_file, cmdstr);
    match = sscanf (cmdstr, "%d", &lctbl->number);
    if (match != 1)
    {
        printf ("Cannot read number of landcover types!\n");
        printf ("Land cover file format error!\n");
        PihmExit (1);
    }

    lctbl->laimax = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->laimin = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->vegfrac = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->albedomin = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->albedomax = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->emissmin = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->emissmax = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->z0min = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->z0max = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->hs = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->snup = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->rgl = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->rsmin = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->rough = (double *)malloc (lctbl->number * sizeof (double));
    lctbl->rzd = (double *)malloc (lctbl->number * sizeof (double));

    for (i = 0; i < lctbl->number; i++)
    {
        NextLine (lc_file, cmdstr);
        match =
            sscanf (cmdstr,
            "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &index, &lctbl->vegfrac[i], &lctbl->rzd[i], &lctbl->rsmin[i],
            &lctbl->rgl[i], &lctbl->hs[i], &lctbl->snup[i],
            &lctbl->laimin[i], &lctbl->laimax[i], &lctbl->emissmin[i],
            &lctbl->emissmax[i], &lctbl->albedomin[i],
            &lctbl->albedomax[i], &lctbl->z0min[i], &lctbl->z0max[i],
            &lctbl->rough[i]);
        if (match != 16 || i != index - 1)
        {
            printf ("Cannot read information of the %dth landcover type!\n",
                i + 1);
            printf ("Landcover file format error!\n");
            PihmExit (1);
        }
    }

    NextLine (lc_file, cmdstr);
    ReadKeywordDouble (cmdstr, "TOPT_DATA", &lctbl->topt);

    NextLine (lc_file, cmdstr);
    ReadKeywordDouble (cmdstr, "CFACTR_DATA", &lctbl->cfactr);

    NextLine (lc_file, cmdstr);
    ReadKeywordDouble (cmdstr, "RSMAX_DATA", &lctbl->rsmax);

    NextLine (lc_file, cmdstr);
    ReadKeywordInt (cmdstr, "BARE", &lctbl->bare);

    NextLine (lc_file, cmdstr);
    ReadKeywordInt (cmdstr, "NATURAL", &lctbl->natural);

    fclose (lc_file);
}

void ReadForc (char *filename, forc_struct *forc)
{
    FILE           *forc_file;  /* Pointer to .forc file */
    char            cmdstr[MAXSTRING];
    int             i, j;
    int             match;
    int             index;

    forc_file = fopen (filename, "r");
    CheckFile (forc_file, filename);

    FindLine (forc_file, "BOF");
    NextLine (forc_file, cmdstr);
    match = sscanf (cmdstr, "%*s %d", &forc->nmeteo);
    if (match != 1)
    {
        printf
            ("Cannot read number of meteorological forcing time series!\n");
        printf (".forc file format error!\n");
        PihmExit (1);
    }

    if (forc->nmeteo > 0)
    {
        forc->meteo =
            (tsdata_struct *)malloc (forc->nmeteo * sizeof (tsdata_struct));

        for (i = 0; i < forc->nmeteo; i++)
        {
            NextLine (forc_file, cmdstr);
            match = sscanf (cmdstr, "%*s %d %*s %lf",
                    &index, &forc->meteo[i].zlvl_wind);
            if (match != 2 || i != index - 1)
            {
                printf ("Cannot read information of the %dth forcing series!\n",
                    i);
                printf (".forc file format error!\n");
                PihmExit (1);
            }
            /* Skip header lines */
            NextLine (forc_file, cmdstr);
            NextLine (forc_file, cmdstr);
            forc->meteo[i].length =
                CountLine (forc_file, 1, "METEO_TS");
        }

        /* Rewind and read */
        FindLine (forc_file, "NUM_METEO_TS");
        for (i = 0; i < forc->nmeteo; i++)
        {
            /* Skip header lines */
            NextLine (forc_file, cmdstr);
            NextLine (forc_file, cmdstr);
            NextLine (forc_file, cmdstr);

            forc->meteo[i].ftime =
                (int *)malloc (forc->meteo[i].length * sizeof (int));
            forc->meteo[i].data =
                (double **)malloc (forc->meteo[i].length *
                sizeof (double *));
            for (j = 0; j < forc->meteo[i].length; j++)
            {
                forc->meteo[i].data[j] =
                    (double *)malloc (NUM_METEO_VAR * sizeof (double));
                NextLine (forc_file, cmdstr);
                ReadTS (cmdstr, &forc->meteo[i].ftime[j],
                    &forc->meteo[i].data[j][0], NUM_METEO_VAR);
            }
        }
    }

    fclose (forc_file);
}

void ReadLAI (char *filename, forc_struct *forc, int numele,
    const atttbl_struct *atttbl)
{
    char            cmdstr[MAXSTRING];
    int             read_lai = 0;
    FILE           *lai_file;
    int             i, j;
    int             index;
    int             match;

    for (i = 0; i < numele; i++)
    {
        if (atttbl->lai[i] > 0)
        {
            read_lai = 1;
            break;
        }
    }

    forc->nlai = 0;

    if (read_lai)
    {
        lai_file = fopen (filename, "r");
        CheckFile (lai_file, filename);

        /* start reading lai_file */
        FindLine (lai_file, "BOF");
        NextLine (lai_file, cmdstr);
        match = sscanf (cmdstr, "%*s %d", &forc->nlai);
        if (match != 1)
        {
            printf ("Cannot read number of LAI time series!\n");
            printf (".lai file format error!\n");
            PihmExit (1);
        }

        if (forc->nlai > 0)
        {
            forc->lai =
                (tsdata_struct *)malloc (forc->nlai * sizeof (tsdata_struct));

            for (i = 0; i < forc->nlai; i++)
            {
                NextLine (lai_file, cmdstr);
                match = sscanf (cmdstr, "%*s %d", &index);
                if (match != 1 || i != index - 1)
                {
                    printf ("Cannot read information of the %dth LAI series!\n",
                        i);
                    printf (".lai file format error!\n");
                    PihmExit (1);
                }
                /* Skip header lines */
                NextLine (lai_file, cmdstr);
                NextLine (lai_file, cmdstr);
                forc->lai[i].length = CountLine (lai_file, 1, "LAI_TS");
            }

            /* Rewind and read */
            FindLine (lai_file, "NUM_LAI_TS");
            for (i = 0; i < forc->nlai; i++)
            {
                /* Skip header lines */
                NextLine (lai_file, cmdstr);
                NextLine (lai_file, cmdstr);
                NextLine (lai_file, cmdstr);

                forc->lai[i].ftime =
                    (int *)malloc (forc->lai[i].length * sizeof (int));
                forc->lai[i].data =
                    (double **)malloc (forc->lai[i].length *
                    sizeof (double *));
                for (j = 0; j < forc->lai[i].length; j++)
                {
                    forc->lai[i].data[j] =
                        (double *)malloc (sizeof (double));
                    NextLine (lai_file, cmdstr);
                    ReadTS (cmdstr, &forc->lai[i].ftime[j],
                        &forc->lai[i].data[j][0], 1);
                }
            }
        }

        fclose (lai_file);
    }
}

void ReadIbc (char *filename, forc_struct *forc)
{
    int             i, j;
    FILE           *ibc_file;   /* Pointer to .ibc file */
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;

    ibc_file = fopen (filename, "r");
    CheckFile (ibc_file, filename);

    /*
     * Start reading ibc_file 
     */
    FindLine (ibc_file, "BOF");
    NextLine (ibc_file, cmdstr);
    match = sscanf (cmdstr, "%*s %d", &forc->nbc);
    if (match != 1)
    {
        printf ("Cannot read number of boundary condition time series!\n");
        printf (".ibc file format error!\n");
        PihmExit (1);
    }

    if (forc->nbc > 0)
    {
        forc->bc =
            (tsdata_struct *)malloc (forc->nbc * sizeof (tsdata_struct));

        for (i = 0; i < forc->nbc; i++)
        {
            NextLine (ibc_file, cmdstr);
            match = sscanf (cmdstr, "%*s %d", &index);
            if (match != 1 || i != index - 1)
            {
                printf
                    ("Cannot read information of the %dth boundary condition series!\n",
                    i);
                printf (".ibc file format error!\n");
                PihmExit (1);
            }
            /* Skip header lines */
            NextLine (ibc_file, cmdstr);
            NextLine (ibc_file, cmdstr);

            forc->bc[i].length = CountLine (ibc_file, 1, "BC_TS");
        }

        /* Rewind and read */
        FindLine (ibc_file, "NUM_BC_TS");
        for (i = 0; i < forc->nbc; i++)
        {
            /* Skip header lines */
            NextLine (ibc_file, cmdstr);
            NextLine (ibc_file, cmdstr);
            NextLine (ibc_file, cmdstr);

            forc->bc[i].ftime =
                (int *)malloc (forc->bc[i].length * sizeof (int));
            forc->bc[i].data =
                (double **)malloc (forc->bc[i].length * sizeof (double *));
            for (j = 0; j < forc->bc[i].length; j++)
            {
                forc->bc[i].data[j] = (double *)malloc (sizeof (double));
                NextLine (ibc_file, cmdstr);
                ReadTS (cmdstr, &forc->bc[i].ftime[j],
                    &forc->bc[i].data[j][0], 1);
            }
        }
    }

    fclose (ibc_file);
}

void ReadPara (char *filename, ctrl_struct *ctrl)
{
    FILE           *para_file;  /* Pointer to .para file */
    char            cmdstr[MAXSTRING];
    int             i;

    for (i = 0; i < NUM_PRINT; i++)
    {
        ctrl->prtvrbl[i] = 0;
    }

    para_file = fopen (filename, "r");
    CheckFile (para_file, filename);

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
    ReadKeywordInt (cmdstr, "SURFFLX", &ctrl->prtvrbl[SURFFLX_CTRL]);

    fclose (para_file);

    if (ctrl->etstep < ctrl->stepsize || ctrl->etstep % ctrl->stepsize > 0)
    {
        printf
            ("ERROR: LSM (ET) step size should be an integral multiple of model step size!\n");
        PihmExit (1);
    }
}

void ReadCalib (char *filename, calib_struct *cal)
{
    char            cmdstr[MAXSTRING];
    FILE           *global_calib;       /* Pointer to .calib file */

    global_calib = fopen (filename, "r");
    CheckFile (global_calib, filename);

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
    ReadKeywordDouble (cmdstr, "EC", &cal->ec);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "ETT", &cal->ett);

    NextLine (global_calib, cmdstr);
    ReadKeywordDouble (cmdstr, "EDIR", &cal->edir);

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

void ReadInit (char *filename, elem_struct *elem, int numele,
    river_struct *riv, int numriv)
{
    FILE           *init_file;
    int             i;

    init_file = fopen (filename, "rb");
    CheckFile (init_file, filename);

    for (i = 0; i < numele; i++)
    {
        fread (&elem[i].ic.intcp, sizeof (double), 1, init_file);
        fread (&elem[i].ic.sneqv, sizeof (double), 1, init_file);
        fread (&elem[i].ic.surf, sizeof (double), 1, init_file);
        fread (&elem[i].ic.unsat, sizeof (double), 1, init_file);
        fread (&elem[i].ic.gw, sizeof (double), 1, init_file);
    }
    for (i = 0; i < numriv; i++)
    {
        fread (&riv[i].ic.stage, sizeof (double), 1, init_file);
        fread (&riv[i].ic.gw, sizeof (double), 1, init_file);
    }

    fclose (init_file);
}

//void FreeData (pihm_struct pihm)
//{
//    int             i, j, k;
//
//    /* Free river input structure */
//    free (pihm->rivtbl.fromnode);
//    free (pihm->rivtbl.tonode);
//    free (pihm->rivtbl.down);
//    free (pihm->rivtbl.leftele);
//    free (pihm->rivtbl.rightele);
//    free (pihm->rivtbl.shp);
//    free (pihm->rivtbl.matl);
//    free (pihm->rivtbl.ic);
//    free (pihm->rivtbl.bc);
//    free (pihm->rivtbl.rsvr);
//
//    free (pihm->shptbl.depth);
//    free (pihm->shptbl.intrpl_ord);
//    free (pihm->shptbl.coeff);
//
//    free (pihm->matltbl.rough);
//    free (pihm->matltbl.cwr);
//    free (pihm->matltbl.ksath);
//    free (pihm->matltbl.ksatv);
//    free (pihm->matltbl.bedthick);
//
//    free (pihm->riv_ic_tbl.stage);
//
//    /* Free mesh input structure */
//    for (i = 0; i < pihm->meshtbl.numele; i++)
//    {
//        free (pihm->meshtbl.node[i]);
//        free (pihm->meshtbl.nabr[i]);
//    }
//    free (pihm->meshtbl.node);
//    free (pihm->meshtbl.nabr);
//    free (pihm->meshtbl.x);
//    free (pihm->meshtbl.y);
//    free (pihm->meshtbl.zmin);
//    free (pihm->meshtbl.zmax);
//
//    /* Free attribute input structure */
//    for (i = 0; i < pihm->meshtbl.numele; i++)
//    {
//        free (pihm->attrib_tbl.bc[i]);
//    }
//    free (pihm->attrib_tbl.soil);
//    free (pihm->attrib_tbl.geol);
//    free (pihm->attrib_tbl.lc);
//    free (pihm->attrib_tbl.bc);
//    free (pihm->attrib_tbl.meteo);
//    free (pihm->attrib_tbl.lai);
//    free (pihm->attrib_tbl.source);
//    free (pihm->attrib_tbl.macropore);
//
//    /* Free soil input structure */
//    free (pihm->soil_tbl.ksatv);
//    free (pihm->soil_tbl.smcmax);
//    free (pihm->soil_tbl.smcmin);
//    free (pihm->soil_tbl.qtz);
//    free (pihm->soil_tbl.alpha);
//    free (pihm->soil_tbl.beta);
//    free (pihm->soil_tbl.areafh);
//    free (pihm->soil_tbl.kmacv);
//    free (pihm->soil_tbl.dinf);
//
//    /* Free geol input structure */
//    free (pihm->geol_tbl.ksath);
//    free (pihm->geol_tbl.ksatv);
//    free (pihm->geol_tbl.smcmax);
//    free (pihm->geol_tbl.smcmin);
//    free (pihm->geol_tbl.alpha);
//    free (pihm->geol_tbl.beta);
//    free (pihm->geol_tbl.areafv);
//    free (pihm->geol_tbl.kmach);
//    free (pihm->geol_tbl.dmac);
//
//    /* Free initial condition */
//    free (pihm->ic.intcp);
//    free (pihm->ic.snow);
//    free (pihm->ic.surf);
//    free (pihm->ic.unsat);
//    free (pihm->ic.gw);
//    free (pihm->ic.rivgw);
//    free (pihm->ic.stage);
//
//    /* Free landcover input structure */
//    free (pihm->lc_tbl.laimax);
//    free (pihm->lc_tbl.laimin);
//    free (pihm->lc_tbl.vegfrac);
//    free (pihm->lc_tbl.albedomin);
//    free (pihm->lc_tbl.albedomax);
//    free (pihm->lc_tbl.emissmin);
//    free (pihm->lc_tbl.emissmax);
//    free (pihm->lc_tbl.z0min);
//    free (pihm->lc_tbl.z0max);
//    free (pihm->lc_tbl.hs);
//    free (pihm->lc_tbl.snup);
//    free (pihm->lc_tbl.rgl);
//    free (pihm->lc_tbl.rsmin);
//    free (pihm->lc_tbl.rough);
//    free (pihm->lc_tbl.rzd);
//
//    /* Free forcing input structure */
//    for (k = 0; k < NUM_TS; k++)
//    {
//        if (pihm->forcing.nts[k] > 0)
//        {
//            for (i = 0; i < pihm->forcing.nts[k]; i++)
//            {
//                for (j = 0; j < pihm->forcing.ts[k][i].length; j++)
//                {
//                    free (pihm->forcing.ts[k][i].data[j]);
//                }
//                free (pihm->forcing.ts[k][i].ftime);
//                free (pihm->forcing.ts[k][i].data);
//            }
//        free (pihm->forcing.ts[k]);
//        }
//    }
//
//    if (pihm->forcing.nts[BC_TS] > 0)
//    {
//        free (pihm->forcing.bc);
//    }
//
//    if (pihm->forcing.nts[METEO_TS] > 0)
//    {
//        for (i = 0; i < NUM_METEO_TS; i++)
//        {
//            free (pihm->forcing.meteo[i]);
//        }
//        free (pihm->forcing.zlvl_wind);
//    }
//    if (pihm->forcing.nts[LAI_TS] > 0)
//    {
//        free (pihm->forcing.lai);
//    }
//    if (pihm->forcing.nts[RIV_TS] > 0)
//    {
//        free (pihm->forcing.riverbc);
//    }
//    if (pihm->forcing.nts[SS_TS] > 0)
//    {
//        free (pihm->forcing.source);
//    }
//
//    free (pihm->ctrl.tout);
//
//    for (i = 0; i < pihm->ctrl.nprint; i++)
//    {
//        free (pihm->prtctrl[i].vrbl);
//        free (pihm->prtctrl[i].buffer);
//    }
//
//    free (pihm->elem);
//    free (pihm->riv);
//}
