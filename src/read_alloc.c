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

    /* Set file names of the input files */
    sprintf (pihm->filename.riv, "input/%s/%s.riv", project, project);
    sprintf (pihm->filename.mesh, "input/%s/%s.mesh", project, project);
    sprintf (pihm->filename.att, "input/%s/%s.att", project, project);
    sprintf (pihm->filename.soil, "input/%s/%s.soil", project, project);
    sprintf (pihm->filename.geol, "input/%s/%s.geol", project, project);
    sprintf (pihm->filename.lc, "input/vegprmt.tbl");
    sprintf (pihm->filename.meteo, "input/%s/%s.meteo", project, project);
    sprintf (pihm->filename.lai, "input/%s/%s.lai", project, project);
    sprintf (pihm->filename.bc, "input/%s/%s.bc", project, project);
    sprintf (pihm->filename.para, "input/%s/%s.para", project, project);
    sprintf (pihm->filename.calib, "input/%s/%s.calib", project, simulation);
    sprintf (pihm->filename.ic, "input/%s/%s.ic", project, simulation);
#ifdef _NOAH_
    sprintf (pihm->filename.lsm, "input/%s/%s.lsm", project, project);
    sprintf (pihm->filename.rad, "input/%s/%s.rad", project, project);
#endif
#ifdef _CYCLES_
    sprintf (pihm->filename.cycles, "input/%s/%s.cycles", project, project);
    sprintf (pihm->filename.soilinit, "input/%s/%s.soilinit", project,
        project);
    sprintf (pihm->filename.crop, "input/%s/%s.crop", project, project);
#endif
#ifdef _BGC_
    sprintf (pihm->filename.bgc, "input/%s/%s.bgc", project, project);
    sprintf (pihm->filename.bgcic, "input/%s/%s.bgcic", project, simulation);
#endif

    /* Read river input file */
    ReadRiv (pihm->filename.riv, &pihm->rivtbl, &pihm->shptbl, &pihm->matltbl,
        &pihm->forc);
    pihm->numriv = pihm->rivtbl.number;

    /* Read mesh structure input file */
    ReadMesh (pihm->filename.mesh, &pihm->meshtbl);
    pihm->numele = pihm->meshtbl.numele;

    /* Read attribute table input file */
    ReadAtt (pihm->filename.att, &pihm->atttbl, pihm->numele);

    /* Read soil input file */
    ReadSoil (pihm->filename.soil, &pihm->soiltbl);

    /* Read geology input file */
    //ReadGeol (pihm->filename.geol, &pihm->geoltbl);

    /* Read land cover input file */
    ReadLC (pihm->filename.lc, &pihm->lctbl);

    /* Read meteorological forcing input file */
    ReadForc (pihm->filename.meteo, &pihm->forc);

    /* Read LAI input file */
    ReadLAI (pihm->filename.lai, &pihm->forc, pihm->numele, &pihm->atttbl);

    /* Read source and sink input file */
    pihm->forc.nsource = 0;
    //ReadSS ();

    /* Read boundary condition input file */
    pihm->forc.nbc = 0;
    ReadBC (pihm->filename.bc, &pihm->forc);

    /* Read model control file */
    ReadPara (pihm->filename.para, &pihm->ctrl);

    /* Read calibration input file */
    ReadCalib (pihm->filename.calib, &pihm->cal);

#ifdef _NOAH_
    /* Read LSM input file */
    ReadLsm (pihm->filename.lsm, &pihm->latitude, &pihm->longitude,
        &pihm->ctrl, &pihm->noahtbl);

    if (pihm->ctrl.rad_mode == TOPO_SOL)
    {
        /* Read radiation input file */
        ReadRad (pihm->filename.rad, &pihm->forc);
    }
#endif

#ifdef _CYCLES_
    /* Read Cycles simulation control file */
    ReadCyclesCtrl (pihm->filename.cycles, &pihm->agtbl, &pihm->ctrl,
        pihm->numele);

    /* Read soil initialization file */
    ReadSoilInit (pihm->filename.soilinit, &pihm->soiltbl);

    /* Read crop description file */
    ReadCrop (pihm->filename.crop, &pihm->croptbl);

    /* Read operation file */
    ReadOperation (&pihm->agtbl, &pihm->mgmttbl, &pihm->croptbl);
#endif

#ifdef _BGC_
    ReadBGC (pihm->filename.bgc, &pihm->ctrl, &pihm->co2, &pihm->ndepctrl,
        pihm->filename.co2, pihm->filename.ndep);

    /* Read Biome-BGC epc files */
    ReadEPC (&pihm->epctbl);

    /* Read CO2 and Ndep files */
    if (pihm->co2.varco2)
    {
        pihm->forc.co2 = (tsdata_struct *)malloc (sizeof (tsdata_struct));
        ReadAnnFile (&pihm->forc.co2[0], pihm->filename.co2);
    }

    if (pihm->ndepctrl.varndep)
    {
        pihm->forc.ndep = (tsdata_struct *)malloc (sizeof (tsdata_struct));
        ReadAnnFile (&pihm->forc.ndep[0], pihm->filename.ndep);
    }
#endif

}

void ReadRiv (char *filename, rivtbl_struct *rivtbl, shptbl_struct *shptbl,
    matltbl_struct *matltbl, forc_struct *forc)
{
    int             i, j;
    FILE           *riv_file;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    /** Open .riv input file */
    riv_file = fopen (filename, "r");
    CheckFile (riv_file, filename);
    if (verbose_mode)
    {
        printf (" Reading %s\n", filename);
    }

    /*
     * Read river segment block
     */
    /* Read number of river segments */
    FindLine (riv_file, "BOF", &lno, filename);
    NextLine (riv_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "NUMRIV", &rivtbl->number, 'i', filename, lno);

    /* Allocate */
    rivtbl->fromnode = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->tonode = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->down = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->leftele = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->rightele = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->shp = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->matl = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->bc = (int *)malloc (rivtbl->number * sizeof (int));
    rivtbl->rsvr = (int *)malloc (rivtbl->number * sizeof (int));

    /* Skip header line */
    NextLine (riv_file, cmdstr, &lno);

    /* Read river segment information */
    for (i = 0; i < rivtbl->number; i++)
    {
        NextLine (riv_file, cmdstr, &lno);
        match = sscanf (cmdstr, "%d %d %d %d %d %d %d %d %d %d",
            &index,
            &rivtbl->fromnode[i], &rivtbl->tonode[i],
            &rivtbl->down[i],
            &rivtbl->leftele[i], &rivtbl->rightele[i],
            &rivtbl->shp[i], &rivtbl->matl[i],
            &rivtbl->bc[i], &rivtbl->rsvr[i]);
        if (match != 10 || i != index - 1)
        {
            fprintf (stderr, "Error in .riv file format.\n");
            fprintf (stderr,
                "Cannot read river segment information for the %dth"
                "segment.\n", i + 1);
            PIHMExit (EXIT_FAILURE);
        }
    }

    /*
     * Read river shape information
     */
    NextLine (riv_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SHAPE", &shptbl->number, 'i', filename, lno);

    /* Allocate */
    shptbl->depth = (double *)malloc (shptbl->number * sizeof (double));
    shptbl->intrpl_ord = (int *)malloc (shptbl->number * sizeof (int));
    shptbl->coeff = (double *)malloc (shptbl->number * sizeof (double));

    /* Skip header line */
    NextLine (riv_file, cmdstr, &lno);

    for (i = 0; i < shptbl->number; i++)
    {
        NextLine (riv_file, cmdstr, &lno);
        match = sscanf (cmdstr, "%d %lf %d %lf",
            &index, &shptbl->depth[i],
            &shptbl->intrpl_ord[i], &shptbl->coeff[i]);
        if (match != 4 || i != index - 1)
        {
            fprintf (stderr, "Error in .riv file format.\n");
            fprintf
                (stderr,
                "Cannot read river shape information for the %dth shape.\n",
                i + 1);
            PIHMExit (EXIT_FAILURE);
        }
    }

    /*
     * Read river material information
     */
    NextLine (riv_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "MATERIAL", &matltbl->number, 'i', filename, lno);

    /* Allocate */
    matltbl->rough = (double *)malloc (matltbl->number * sizeof (double));
    matltbl->cwr = (double *)malloc (matltbl->number * sizeof (double));
    matltbl->ksath = (double *)malloc (matltbl->number * sizeof (double));
    matltbl->ksatv = (double *)malloc (matltbl->number * sizeof (double));
    matltbl->bedthick = (double *)malloc (matltbl->number * sizeof (double));

    /* Skip header line */
    NextLine (riv_file, cmdstr, &lno);

    for (i = 0; i < matltbl->number; i++)
    {
        NextLine (riv_file, cmdstr, &lno);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf",
            &index,
            &matltbl->rough[i], &matltbl->cwr[i],
            &matltbl->ksath[i], &matltbl->ksatv[i], &matltbl->bedthick[i]);
        if (match != 6 || i != index - 1)
        {
            fprintf (stderr, "Error in .riv file format.\n");
            fprintf (stderr,
                "Cannot read information for the %dth material.\n", i + 1);
            PIHMExit (EXIT_FAILURE);
        }
    }

    /*
     * Read river boundary condition block
     */
    NextLine (riv_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "BC", &forc->nriverbc, 'i', filename, lno);

    if (forc->nriverbc > 0)
    {
        forc->riverbc =
            (tsdata_struct *)malloc (forc->nriverbc * sizeof (tsdata_struct));

        for (i = 0; i < forc->nriverbc; i++)
        {
            NextLine (riv_file, cmdstr, &lno);
            match = sscanf (cmdstr, "%*s %d", &index);
            if (match != 1 || i != index - 1)
            {
                fprintf (stderr, "Error in .riv file format.\n");
                fprintf (stderr,
                    "Cannot read information of the %dth river boudnary condition!\n",
                    i);
                PIHMExit (EXIT_FAILURE);
            }
            NextLine (riv_file, cmdstr, &lno);
            NextLine (riv_file, cmdstr, &lno);
            forc->riverbc[i].length =
                CountLine (riv_file, cmdstr, 2, "RIV_TS", "RES");
        }

        FindLine (riv_file, "BOF", &lno, filename);
        FindLine (riv_file, "BC", &lno, filename);
        for (i = 0; i < forc->nriverbc; i++)
        {
            NextLine (riv_file, cmdstr, &lno);
            NextLine (riv_file, cmdstr, &lno);
            NextLine (riv_file, cmdstr, &lno);

            forc->riverbc[i].data =
                (double **)malloc ((forc->riverbc[i].length) *
                sizeof (double *));
            forc->riverbc[i].ftime =
                (int *)malloc ((forc->riverbc[i].length) * sizeof (int));
            for (j = 0; j < forc->riverbc[i].length; j++)
            {
                forc->riverbc[i].data[j] = (double *)malloc (sizeof (double));
                NextLine (riv_file, cmdstr, &lno);
                if (!ReadTS (cmdstr, &forc->riverbc[i].ftime[j],
                        &forc->riverbc[i].data[j][0], 1))
                {
                    fprintf (stderr, "Error reading %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }
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
    int             lno = 0;

    /*
     * Open .mesh input file
     */
    mesh_file = fopen (filename, "r");
    CheckFile (mesh_file, filename);
    if (verbose_mode)
    {
        printf (" Reading %s\n", filename);
    }

    /*
     * Read element mesh block
     */
    NextLine (mesh_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "NUMELE", &meshtbl->numele, 'i', filename, lno);

    meshtbl->node = (int **)malloc (meshtbl->numele * sizeof (int *));
    meshtbl->nabr = (int **)malloc (meshtbl->numele * sizeof (int *));

    /* Skip header line */
    NextLine (mesh_file, cmdstr, &lno);

    for (i = 0; i < meshtbl->numele; i++)
    {
        meshtbl->node[i] = (int *)malloc (3 * sizeof (int));
        meshtbl->nabr[i] = (int *)malloc (3 * sizeof (int));

        NextLine (mesh_file, cmdstr, &lno);
        match = sscanf (cmdstr, "%d %d %d %d %d %d %d",
            &index,
            &meshtbl->node[i][0], &meshtbl->node[i][1],
            &meshtbl->node[i][2], &meshtbl->nabr[i][0],
            &meshtbl->nabr[i][1], &meshtbl->nabr[i][2]);
        if (match != 7 || i != index - 1)
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            fprintf (stderr, "Cannot read information of the %dth element.\n",
                i + 1);
            PIHMExit (EXIT_FAILURE);
        }
    }

    /*
     * Read node block
     */
    NextLine (mesh_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "NUMNODE", &meshtbl->numnode, 'i', filename, lno);

    /* Skip header line */
    NextLine (mesh_file, cmdstr, &lno);

    meshtbl->x = (double *)malloc (meshtbl->numnode * sizeof (double));
    meshtbl->y = (double *)malloc (meshtbl->numnode * sizeof (double));
    meshtbl->zmin = (double *)malloc (meshtbl->numnode * sizeof (double));
    meshtbl->zmax = (double *)malloc (meshtbl->numnode * sizeof (double));

    for (i = 0; i < meshtbl->numnode; i++)
    {
        NextLine (mesh_file, cmdstr, &lno);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf",
            &index,
            &meshtbl->x[i], &meshtbl->y[i],
            &meshtbl->zmin[i], &meshtbl->zmax[i]);
        if (match != 5 || i != index - 1)
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            fprintf (stderr, "Cannot read information of the %dth node!\n",
                i + 1);
            PIHMExit (EXIT_FAILURE);
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
    int             lno = 0;

    att_file = fopen (filename, "r");
    CheckFile (att_file, filename);
    if (verbose_mode)
    {
        printf (" Reading %s\n", filename);
    }

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

    NextLine (att_file, cmdstr, &lno);
    for (i = 0; i < numele; i++)
    {
        NextLine (att_file, cmdstr, &lno);
        match = sscanf (cmdstr, "%d %d %d %d %d %d %d %d %d %d", &index,
            &atttbl->soil[i], &atttbl->geol[i], &atttbl->lc[i],
            &atttbl->meteo[i], &atttbl->lai[i], &atttbl->source[i],
            &atttbl->bc[i][0], &atttbl->bc[i][1], &atttbl->bc[i][2]);
        if (match != 10)
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            fprintf (stderr, "Cannot read information of the %dth element.\n",
                i + 1);
            PIHMExit (EXIT_FAILURE);
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
    int             texture;
    const int       TOPSOIL = 1;
    const int       SUBSOIL = 0;
    int             ptf_used = 0;
    int             lno = 0;

    soil_file = fopen (filename, "r");
    CheckFile (soil_file, filename);
    if (verbose_mode)
    {
        printf (" Reading %s\n", filename);
    }

    /* Start reading soil file */
    NextLine (soil_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "NUMSOIL", &soiltbl->number, 'i', filename, lno);

    soiltbl->silt = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->clay = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->om = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->bd = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->kinfv = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->ksatv = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->ksath = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->smcmax = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->smcmin = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->qtz = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->alpha = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->beta = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->areafh = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->areafv = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->dmac = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->smcref = (double *)malloc (soiltbl->number * sizeof (double));
    soiltbl->smcwlt = (double *)malloc (soiltbl->number * sizeof (double));

    /* Skip header line */
    NextLine (soil_file, cmdstr, &lno);

    for (i = 0; i < soiltbl->number; i++)
    {
        NextLine (soil_file, cmdstr, &lno);
        match = sscanf (cmdstr,
            "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &index, &soiltbl->silt[i], &soiltbl->clay[i], &soiltbl->om[i],
            &soiltbl->bd[i],
            &soiltbl->kinfv[i], &soiltbl->ksatv[i], &soiltbl->ksath[i],
            &soiltbl->smcmax[i], &soiltbl->smcmin[i],
            &soiltbl->alpha[i], &soiltbl->beta[i],
            &soiltbl->areafh[i], &soiltbl->areafv[i],
            &soiltbl->dmac[i], &soiltbl->qtz[i]);

        if (match != 16 || i != index - 1)
        {
            fprintf (stderr, "Error reading %s.\n", filename);
            fprintf (stderr,
                "Cannot read information of the %dth soil type!\n", i + 1);
            PIHMExit (EXIT_FAILURE);
        }

        /* Fill in missing organic matter and bulk density values */
        soiltbl->om[i] = (soiltbl->om[i] > 0.0) ? soiltbl->om[i] : 2.5;
        soiltbl->bd[i] = (soiltbl->bd[i] > 0.0) ? soiltbl->bd[i] : 1.3;

        /* Fill missing hydraulic properties using PTFs */
        if (soiltbl->kinfv[i] < 0.0)
        {
            soiltbl->kinfv[i] =
                PtfKV (soiltbl->silt[i], soiltbl->clay[i], soiltbl->om[i],
                soiltbl->bd[i], TOPSOIL);
            ptf_used = 1;
        }
        if (soiltbl->ksatv[i] < 0.0)
        {
            soiltbl->ksatv[i] =
                PtfKV (soiltbl->silt[i], soiltbl->clay[i], soiltbl->om[i],
                soiltbl->bd[i], SUBSOIL);
            ptf_used = 1;
        }
        if (soiltbl->ksath[i] < 0.0)
        {
            soiltbl->ksath[i] = 10.0 * soiltbl->ksatv[i];
            ptf_used = 1;
        }
        if (soiltbl->smcmax[i] < 0.0)
        {
            soiltbl->smcmax[i] =
                PtfThetaS (soiltbl->silt[i], soiltbl->clay[i], soiltbl->om[i],
                soiltbl->bd[i], SUBSOIL);
            ptf_used = 1;
        }
        if (soiltbl->smcmin[i] < 0.0)
        {
            soiltbl->smcmin[i] =
                PtfThetaR (soiltbl->silt[i], soiltbl->clay[i], soiltbl->om[i],
                soiltbl->bd[i], SUBSOIL);
            ptf_used = 1;
        }
        if (soiltbl->alpha[i] < 0.0)
        {
            soiltbl->alpha[i] =
                PtfAlpha (soiltbl->silt[i], soiltbl->clay[i], soiltbl->om[i],
                soiltbl->bd[i], SUBSOIL);
            ptf_used = 1;
        }
        if (soiltbl->beta[i] < 0.0)
        {
            soiltbl->beta[i] =
                PtfBeta (soiltbl->silt[i], soiltbl->clay[i], soiltbl->om[i],
                soiltbl->bd[i], SUBSOIL);
            ptf_used = 1;
        }
        if (soiltbl->qtz[i] < 0.0)
        {
            texture = SoilTex (soiltbl->silt[i], soiltbl->clay[i]);
            soiltbl->qtz[i] = Qtz (texture);
            ptf_used = 1;
        }

        /* Calculate field capacity and wilting point */
        soiltbl->smcref[i] =
            FieldCapacity (soiltbl->alpha[i], soiltbl->beta[i],
            soiltbl->ksatv[i], soiltbl->smcmax[i], soiltbl->smcmin[i]);
        soiltbl->smcwlt[i] =
            WiltingPoint (soiltbl->smcmax[i], soiltbl->smcmin[i],
            soiltbl->alpha[i], soiltbl->beta[i]);
    }

    NextLine (soil_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "DINF", &soiltbl->dinf, 'd', filename, lno);

    NextLine (soil_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "KMACV_RO", &soiltbl->kmacv_ro, 'd', filename, lno);

    NextLine (soil_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "KMACH_RO", &soiltbl->kmach_ro, 'd', filename, lno);

    if (ptf_used)
    {
        printf ("%-7s\t%-15s\t%-15s\t%-15s\t%-7s\t%-7s\t%-7s\t%-7s\t%-7s\n",
            "TYPE", "KINFV", "KSATV", "KSATH", "SMCMAX", "SMCMIN", "ALPHA",
            "BETA", "QTZ");
        for (i = 0; i < soiltbl->number; i++)
        {
            printf
                ("%-7d\t%-15.3le\t%-15.3le\t%-15.3le\t%-7.3lf\t%-7.3lf\t"
                "%-7.3lf\t%-7.3lf\t%-7.3lf\n",
                i + 1, soiltbl->kinfv[i], soiltbl->ksatv[i],
                soiltbl->ksath[i], soiltbl->smcmax[i], soiltbl->smcmin[i],
                soiltbl->alpha[i], soiltbl->beta[i], soiltbl->qtz[i]);
        }
    }

    fclose (soil_file);
}

//void ReadGeol (char *filename, geoltbl_struct *geoltbl)
//{
//    FILE           *geol_file;  /* Pointer to .geol file */
//    int             i;
//    char            cmdstr[MAXSTRING];
//    int             match;
//    int             index;
//
//    geol_file = fopen (filename, "r");
//
//    if (NULL == geol_file)
//    {
//        fprintf (stderr, "Error opening %s.\n", filename);
//        PIHMExit (EXIT_FAILURE);
//    }
//
//    if (verbose_mode)
//    {
//        printf ("Reading %s.\n", filename);
//    }
//
//    /* start reading geol_file */
//    NextLine (geol_file, cmdstr);
//    match = sscanf (cmdstr, "%d", &geoltbl->number);
//    if (match != 1)
//    {
//        printf ("Cannot read number of geology types!\n");
//        printf (".geol file format error!\n");
//        PihmExit (1);
//    }
//
//    geoltbl->silt = (double *)malloc (geoltbl->number * sizeof (double));
//    geoltbl->clay = (double *)malloc (geoltbl->number * sizeof (double));
//    geoltbl->om = (double *)malloc (geoltbl->number * sizeof (double));
//    geoltbl->bd = (double *)malloc (geoltbl->number * sizeof (double));
//    geoltbl->ksath = (double *)malloc (geoltbl->number * sizeof (double));
//    geoltbl->ksatv = (double *)malloc (geoltbl->number * sizeof (double));
//    geoltbl->smcmax = (double *)malloc (geoltbl->number * sizeof (double));
//    geoltbl->smcmin = (double *)malloc (geoltbl->number * sizeof (double));
//    geoltbl->alpha = (double *)malloc (geoltbl->number * sizeof (double));
//    geoltbl->beta = (double *)malloc (geoltbl->number * sizeof (double));
//
//    for (i = 0; i < geoltbl->number; i++)
//    {
//        NextLine (geol_file, cmdstr);
//        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
//            &index,
//            &geoltbl->silt[i], &geoltbl->clay[i], &geoltbl->om[i],
//            &geoltbl->bd[i],
//            &geoltbl->ksatv[i], &geoltbl->ksath[i],
//            &geoltbl->smcmax[i], &geoltbl->smcmin[i],
//            &geoltbl->alpha[i], &geoltbl->beta[i]);
//        if (match != 11 || i != index - 1)
//        {
//            printf ("Cannot read information of the %dth geology type!\n",
//                i + 1);
//            printf (".geol file format error!\n");
//            PihmExit (1);
//        }
//    }
//
//    fclose (geol_file);
//}

void ReadLC (char *filename, lctbl_struct *lctbl)
{
    FILE           *lc_file;    /* Pointer to .lc file */
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    lc_file = fopen (filename, "r");
    CheckFile (lc_file, filename);
    if (verbose_mode)
    {
        printf (" Reading %s\n", filename);
    }

    /* Start reading land cover file */
    NextLine (lc_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "NUMLC", &lctbl->number, 'i', filename, lno);

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

    /* Skip header line */
    NextLine (lc_file, cmdstr, &lno);

    for (i = 0; i < lctbl->number; i++)
    {
        NextLine (lc_file, cmdstr, &lno);
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
            fprintf (stderr, "Error reading %s.\n", filename);
            fprintf (stderr,
                "Cannot read information of the %dth landcover type!\n",
                i + 1);
            PIHMExit (EXIT_FAILURE);
        }
    }

    NextLine (lc_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "TOPT_DATA", &lctbl->topt, 'd', filename, lno);

    NextLine (lc_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "CFACTR_DATA", &lctbl->cfactr, 'd', filename, lno);

    NextLine (lc_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RSMAX_DATA", &lctbl->rsmax, 'd', filename, lno);

    NextLine (lc_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "BARE", &lctbl->bare, 'i', filename, lno);

    NextLine (lc_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "NATURAL", &lctbl->natural, 'i', filename, lno);

    fclose (lc_file);
}

void ReadForc (char *filename, forc_struct *forc)
{
    FILE           *meteo_file; /* Pointer to .forc file */
    char            cmdstr[MAXSTRING];
    int             i, j;
    int             match;
    int             index;
    int             lno = 0;

    meteo_file = fopen (filename, "r");
    CheckFile (meteo_file, filename);
    if (verbose_mode)
    {
        printf (" Reading %s\n", filename);
    }

    FindLine (meteo_file, "BOF", &lno, filename);

    forc->nmeteo = CountOccurance (meteo_file, "METEO_TS");

    FindLine (meteo_file, "BOF", &lno, filename);
    if (forc->nmeteo > 0)
    {
        forc->meteo =
            (tsdata_struct *)malloc (forc->nmeteo * sizeof (tsdata_struct));

        NextLine (meteo_file, cmdstr, &lno);
        for (i = 0; i < forc->nmeteo; i++)
        {
            match = sscanf (cmdstr, "%*s %d %*s %lf",
                &index, &forc->meteo[i].zlvl_wind);
            if (match != 2 || i != index - 1)
            {
                fprintf (stderr, "Error reading %s.\n", filename);
                fprintf (stderr,
                    "Cannot read information of the %dth forcing series!\n",
                    i);
                PIHMExit (EXIT_FAILURE);
            }
            /* Skip header lines */
            NextLine (meteo_file, cmdstr, &lno);
            NextLine (meteo_file, cmdstr, &lno);
            forc->meteo[i].length =
                CountLine (meteo_file, cmdstr, 1, "METEO_TS");
        }

        /* Rewind and read */
        FindLine (meteo_file, "BOF", &lno, filename);
        for (i = 0; i < forc->nmeteo; i++)
        {
            /* Skip header lines */
            NextLine (meteo_file, cmdstr, &lno);
            NextLine (meteo_file, cmdstr, &lno);
            NextLine (meteo_file, cmdstr, &lno);

            forc->meteo[i].ftime =
                (int *)malloc (forc->meteo[i].length * sizeof (int));
            forc->meteo[i].data =
                (double **)malloc (forc->meteo[i].length * sizeof (double *));
            for (j = 0; j < forc->meteo[i].length; j++)
            {
                forc->meteo[i].data[j] =
                    (double *)malloc (NUM_METEO_VAR * sizeof (double));
                NextLine (meteo_file, cmdstr, &lno);
                if (!ReadTS (cmdstr, &forc->meteo[i].ftime[j],
                        &forc->meteo[i].data[j][0], NUM_METEO_VAR))
                {
                    fprintf (stderr, "Error reading %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }
            }
        }
    }

    fclose (meteo_file);
}

void ReadLAI (char *filename, forc_struct *forc, int numele,
    const atttbl_struct *atttbl)
{
    char            cmdstr[MAXSTRING];
    int             read_lai = 0;
    FILE           *lai_file;
    int             i, j;
    int             index;
    int             lno = 0;

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
        if (verbose_mode)
        {
            printf (" Reading %s\n", filename);
        }

        /* start reading lai_file */
        FindLine (lai_file, "BOF", &lno, filename);

        forc->nlai = CountOccurance (lai_file, "LAI_TS");

        FindLine (lai_file, "BOF", &lno, filename);
        if (forc->nlai > 0)
        {
            forc->lai =
                (tsdata_struct *)malloc (forc->nlai * sizeof (tsdata_struct));

            NextLine (lai_file, cmdstr, &lno);
            for (i = 0; i < forc->nlai; i++)
            {
                ReadKeyword (cmdstr, "LAI_TS", &index, 'i', filename, lno);

                if (i != index - 1)
                {
                    fprintf (stderr, "Error reading %s.\n", filename);
                    fprintf (stderr,
                        "Cannot read information of the %dth LAI series.\n",
                        i);
                    PIHMExit (EXIT_FAILURE);
                }
                /* Skip header lines */
                NextLine (lai_file, cmdstr, &lno);
                NextLine (lai_file, cmdstr, &lno);
                forc->lai[i].length =
                    CountLine (lai_file, cmdstr, 1, "LAI_TS");
            }

            /* Rewind and read */
            FindLine (lai_file, "BOF", &lno, filename);
            for (i = 0; i < forc->nlai; i++)
            {
                /* Skip header lines */
                NextLine (lai_file, cmdstr, &lno);
                NextLine (lai_file, cmdstr, &lno);
                NextLine (lai_file, cmdstr, &lno);

                forc->lai[i].ftime =
                    (int *)malloc (forc->lai[i].length * sizeof (int));
                forc->lai[i].data =
                    (double **)malloc (forc->lai[i].length *
                    sizeof (double *));
                for (j = 0; j < forc->lai[i].length; j++)
                {
                    forc->lai[i].data[j] = (double *)malloc (sizeof (double));
                    NextLine (lai_file, cmdstr, &lno);
                    if (!ReadTS (cmdstr, &forc->lai[i].ftime[j],
                            &forc->lai[i].data[j][0], 1))
                    {
                        fprintf (stderr, "Error reading %s.\n", filename);
                        PIHMExit (EXIT_FAILURE);
                    }
                }
            }
        }

        fclose (lai_file);
    }
}

void ReadBC (char *filename, forc_struct *forc)
{
    int             i, j;
    FILE           *bc_file;    /* Pointer to .ibc file */
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    bc_file = fopen (filename, "r");
    CheckFile (bc_file, filename);
    if (verbose_mode)
    {
        printf (" Reading %s\n", filename);
    }

    /*
     * Start reading bc_file 
     */
    FindLine (bc_file, "BOF", &lno, filename);
    NextLine (bc_file, cmdstr, &lno);
    match = sscanf (cmdstr, "%*s %d", &forc->nbc);
    if (match != 1)
    {
        fprintf (stderr, "Error reading %s.\n", filename);
        fprintf (stderr,
            "Cannot read number of boundary condition time series.\n");
        PIHMExit (EXIT_FAILURE);
    }

    if (forc->nbc > 0)
    {
        forc->bc =
            (tsdata_struct *)malloc (forc->nbc * sizeof (tsdata_struct));

        for (i = 0; i < forc->nbc; i++)
        {
            NextLine (bc_file, cmdstr, &lno);
            match = sscanf (cmdstr, "%*s %d", &index);
            if (match != 1 || i != index - 1)
            {
                fprintf (stderr, "Error reading %s.\n", filename);
                fprintf (stderr,
                    "Cannot read information of the %dth boundary condition series!\n",
                    i);
                PIHMExit (EXIT_FAILURE);
            }
            /* Skip header lines */
            NextLine (bc_file, cmdstr, &lno);
            NextLine (bc_file, cmdstr, &lno);

            forc->bc[i].length = CountLine (bc_file, cmdstr, 1, "BC_TS");
        }

        /* Rewind and read */
        FindLine (bc_file, "BOF", &lno, filename);
        FindLine (bc_file, "NUM_BC_TS", &lno, filename);
        for (i = 0; i < forc->nbc; i++)
        {
            /* Skip header lines */
            NextLine (bc_file, cmdstr, &lno);
            NextLine (bc_file, cmdstr, &lno);
            NextLine (bc_file, cmdstr, &lno);

            forc->bc[i].ftime =
                (int *)malloc (forc->bc[i].length * sizeof (int));
            forc->bc[i].data =
                (double **)malloc (forc->bc[i].length * sizeof (double *));
            for (j = 0; j < forc->bc[i].length; j++)
            {
                forc->bc[i].data[j] = (double *)malloc (sizeof (double));
                NextLine (bc_file, cmdstr, &lno);
                if (!ReadTS (cmdstr, &forc->bc[i].ftime[j],
                        &forc->bc[i].data[j][0], 1))
                {
                    fprintf (stderr, "Error reading %s.\n", filename);
                    PIHMExit (EXIT_FAILURE);
                }
            }
        }
    }

    fclose (bc_file);
}

void ReadPara (char *filename, ctrl_struct *ctrl)
{
    FILE           *para_file;  /* Pointer to .para file */
    char            cmdstr[MAXSTRING];
    int             i;
    int             lno = 0;

    for (i = 0; i < MAXPRINT; i++)
    {
        ctrl->prtvrbl[i] = 0;
    }

    para_file = fopen (filename, "r");
    CheckFile (para_file, filename);
    if (verbose_mode)
    {
        printf (" Reading %s\n", filename);
    }

    /* start reading para_file */
    /* Read through parameter file to find parameters */
    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "INIT_MODE", &ctrl->init_type, 'i', filename, lno);
    ctrl->init_type = (ctrl->init_type > 0) ? 1 : 0;

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "ASCII_OUTPUT", &ctrl->ascii, 'i', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "WRITE_IC", &ctrl->write_ic, 'i', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "UNSAT_MODE", &ctrl->unsat_mode, 'i', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SURF_MODE", &ctrl->surf_mode, 'i', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIV_MODE", &ctrl->riv_mode, 'i', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SOLVER", &ctrl->solver, 'i', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "ABSTOL", &ctrl->abstol, 'd', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RELTOL", &ctrl->reltol, 'd', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "INIT_SOLVER_STEP", &ctrl->initstep, 'd', filename,
        lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "MAX_SOLVER_STEP", &ctrl->maxstep, 'd', filename,
        lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "LSM_STEP", &ctrl->etstep, 'i', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "START", &ctrl->starttime, 't', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "END", &ctrl->endtime, 't', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "MODEL_STEPSIZE", &ctrl->stepsize, 'i', filename,
        lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SURF", &ctrl->prtvrbl[SURF_CTRL], 'i', filename,
        lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "UNSAT", &ctrl->prtvrbl[UNSAT_CTRL], 'i', filename,
        lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "GW", &ctrl->prtvrbl[GW_CTRL], 'i', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVSTG", &ctrl->prtvrbl[RIVSTG_CTRL], 'i', filename,
        lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVGW", &ctrl->prtvrbl[RIVGW_CTRL], 'i', filename,
        lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SNOW", &ctrl->prtvrbl[SNOW_CTRL], 'i', filename,
        lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "CMC", &ctrl->prtvrbl[CMC_CTRL], 'i', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "INFIL", &ctrl->prtvrbl[INFIL_CTRL], 'i', filename,
        lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RECHARGE", &ctrl->prtvrbl[RECHARGE_CTRL], 'i',
        filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "EC", &ctrl->prtvrbl[EC_CTRL], 'i', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "ETT", &ctrl->prtvrbl[ETT_CTRL], 'i', filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "EDIR", &ctrl->prtvrbl[EDIR_CTRL], 'i', filename,
        lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVFLX0", &ctrl->prtvrbl[RIVFLX0_CTRL], 'i',
        filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVFLX1", &ctrl->prtvrbl[RIVFLX1_CTRL], 'i',
        filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVFLX2", &ctrl->prtvrbl[RIVFLX2_CTRL], 'i',
        filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVFLX3", &ctrl->prtvrbl[RIVFLX3_CTRL], 'i',
        filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVFLX4", &ctrl->prtvrbl[RIVFLX4_CTRL], 'i',
        filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVFLX5", &ctrl->prtvrbl[RIVFLX5_CTRL], 'i',
        filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVFLX6", &ctrl->prtvrbl[RIVFLX6_CTRL], 'i',
        filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVFLX7", &ctrl->prtvrbl[RIVFLX7_CTRL], 'i',
        filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVFLX8", &ctrl->prtvrbl[RIVFLX8_CTRL], 'i',
        filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVFLX9", &ctrl->prtvrbl[RIVFLX9_CTRL], 'i',
        filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIVFLX10", &ctrl->prtvrbl[RIVFLX10_CTRL], 'i',
        filename, lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SUBFLX", &ctrl->prtvrbl[SUBFLX_CTRL], 'i', filename,
        lno);

    NextLine (para_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "SURFFLX", &ctrl->prtvrbl[SURFFLX_CTRL], 'i',
        filename, lno);

    fclose (para_file);

    if (ctrl->etstep < ctrl->stepsize || ctrl->etstep % ctrl->stepsize > 0)
    {
        fprintf (stderr,
            "Error: LSM (ET) step size should be an integral multiple of"
            "model step size!\n");
        PIHMExit (EXIT_FAILURE);
    }
}

void ReadCalib (char *filename, calib_struct *cal)
{
    char            cmdstr[MAXSTRING];
    FILE           *global_calib;       /* Pointer to .calib file */
    int             lno = 0;

    global_calib = fopen (filename, "r");
    CheckFile (global_calib, filename);
    if (verbose_mode)
    {
        printf (" Reading %s\n", filename);
    }

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "KSATH", &cal->ksath, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "KSATV", &cal->ksatv, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "KINF", &cal->kinfv, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "KMACSATH", &cal->kmach, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "KMACSATV", &cal->kmacv, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "DINF", &cal->dinf, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "DROOT", &cal->rzd, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "DMAC", &cal->dmac, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "POROSITY", &cal->porosity, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "ALPHA", &cal->alpha, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "BETA", &cal->beta, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "MACVF", &cal->areafv, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "MACHF", &cal->areafh, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "VEGFRAC", &cal->vegfrac, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "ALBEDO", &cal->albedo, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "ROUGH", &cal->rough, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "EC", &cal->ec, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "ETT", &cal->ett, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "EDIR", &cal->edir, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "ROUGH_RIV", &cal->rivrough, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "KRIVH", &cal->rivksath, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "KRIVV", &cal->rivksatv, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "BEDTHCK", &cal->rivbedthick, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIV_DPTH", &cal->rivdepth, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "RIV_WDTH", &cal->rivshpcoeff, 'd', filename, lno);

#ifdef _NOAH_
    FindLine (global_calib, "LSM_CALIBRATION", &lno, filename);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "DRIP", &cal->drip, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "CMCMAX", &cal->cmcmax, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "RS", &cal->rsmin, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "CZIL", &cal->czil, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "FXEXP", &cal->fxexp, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "CFACTR", &cal->cfactr, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "RGL", &cal->rgl, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "HS", &cal->hs, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "REFSMC", &cal->smcref, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "WLTSMC", &cal->smcwlt, 'd', filename, lno);
#endif

#ifdef _RT_
    FindLine (global_calib, "RT_CALIBRATION", &lno, filename);

    CS->Cal.PCO2 = 1.0;
    CS->Cal.Keq = 1.0;
    CS->Cal.Site_den = 1.0;
    CS->Cal.SSA = 1.0;
    CS->Cal.Prep_conc = 1.0;
#endif

    /*
     * Scenarios
     */
    FindLine (global_calib, "SCENARIO", &lno, filename);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "PRCP", &cal->prcp, 'd', filename, lno);

    NextLine (global_calib, cmdstr, &lno);
    ReadKeyword (cmdstr, "SFCTMP", &cal->sfctmp, 'd', filename, lno);


    /* Finish reading calib file */
    fclose (global_calib);
}

void ReadIC (char *filename, elem_struct *elem, int numele,
    river_struct *riv, int numriv)
{
    FILE           *ic_file;
    int             i;
    int             size;

    ic_file = fopen (filename, "rb");
    CheckFile (ic_file, filename);
    if (verbose_mode)
    {
        printf (" Reading %s\n", filename);
    }

    fseek (ic_file, 0L, SEEK_END);
    size = ftell (ic_file);

    if (size !=
        sizeof (ic_struct) * numele + sizeof (river_ic_struct) * numriv)
    {
        fprintf (stderr, "Error: %s file size does not match requirement.\n",
            filename);
        fprintf (stderr, "Please use a correct initial condition file.\n");
        PIHMExit (EXIT_FAILURE);
    }

    fseek (ic_file, 0L, SEEK_SET);

    for (i = 0; i < numele; i++)
    {
        fread (&elem[i].ic, sizeof (ic_struct), 1, ic_file);
    }

    for (i = 0; i < numriv; i++)
    {
        fread (&riv[i].ic, sizeof (river_ic_struct), 1, ic_file);
    }

    fclose (ic_file);
}

void FreeData (pihm_struct pihm)
{
    int             i, j;

    /* Free river input structure */
    free (pihm->rivtbl.fromnode);
    free (pihm->rivtbl.tonode);
    free (pihm->rivtbl.down);
    free (pihm->rivtbl.leftele);
    free (pihm->rivtbl.rightele);
    free (pihm->rivtbl.shp);
    free (pihm->rivtbl.matl);
    free (pihm->rivtbl.bc);
    free (pihm->rivtbl.rsvr);

    free (pihm->shptbl.depth);
    free (pihm->shptbl.intrpl_ord);
    free (pihm->shptbl.coeff);

    free (pihm->matltbl.rough);
    free (pihm->matltbl.cwr);
    free (pihm->matltbl.ksath);
    free (pihm->matltbl.ksatv);
    free (pihm->matltbl.bedthick);

    /* Free mesh input structure */
    for (i = 0; i < pihm->meshtbl.numele; i++)
    {
        free (pihm->meshtbl.node[i]);
        free (pihm->meshtbl.nabr[i]);
    }
    free (pihm->meshtbl.node);
    free (pihm->meshtbl.nabr);
    free (pihm->meshtbl.x);
    free (pihm->meshtbl.y);
    free (pihm->meshtbl.zmin);
    free (pihm->meshtbl.zmax);

    /* Free attribute input structure */
    for (i = 0; i < pihm->meshtbl.numele; i++)
    {
        free (pihm->atttbl.bc[i]);
    }
    free (pihm->atttbl.soil);
    free (pihm->atttbl.geol);
    free (pihm->atttbl.lc);
    free (pihm->atttbl.bc);
    free (pihm->atttbl.meteo);
    free (pihm->atttbl.lai);
    free (pihm->atttbl.source);

    /* Free soil input structure */
    free (pihm->soiltbl.silt);
    free (pihm->soiltbl.clay);
    free (pihm->soiltbl.om);
    free (pihm->soiltbl.bd);
    free (pihm->soiltbl.kinfv);
    free (pihm->soiltbl.ksatv);
    free (pihm->soiltbl.ksath);
    free (pihm->soiltbl.smcmax);
    free (pihm->soiltbl.smcmin);
    free (pihm->soiltbl.qtz);
    free (pihm->soiltbl.alpha);
    free (pihm->soiltbl.beta);
    free (pihm->soiltbl.areafh);
    free (pihm->soiltbl.areafv);
    free (pihm->soiltbl.dmac);
    free (pihm->soiltbl.smcref);
    free (pihm->soiltbl.smcwlt);

    /* Free geol input structure */
    //free (pihm->geol_tbl.ksath);
    //free (pihm->geol_tbl.ksatv);
    //free (pihm->geol_tbl.smcmax);
    //free (pihm->geol_tbl.smcmin);
    //free (pihm->geol_tbl.alpha);
    //free (pihm->geol_tbl.beta);
    //free (pihm->geol_tbl.areafv);
    //free (pihm->geol_tbl.kmach);
    //free (pihm->geol_tbl.dmac);

    /* Free landcover input structure */
    free (pihm->lctbl.laimax);
    free (pihm->lctbl.laimin);
    free (pihm->lctbl.vegfrac);
    free (pihm->lctbl.albedomin);
    free (pihm->lctbl.albedomax);
    free (pihm->lctbl.emissmin);
    free (pihm->lctbl.emissmax);
    free (pihm->lctbl.z0min);
    free (pihm->lctbl.z0max);
    free (pihm->lctbl.hs);
    free (pihm->lctbl.snup);
    free (pihm->lctbl.rgl);
    free (pihm->lctbl.rsmin);
    free (pihm->lctbl.rough);
    free (pihm->lctbl.rzd);

    /* Free forcing input structure */
    if (pihm->forc.nriverbc > 0)
    {
        for (i = 0; i < pihm->forc.nriverbc; i++)
        {
            for (j = 0; j < pihm->forc.riverbc[i].length; j++)
            {
                free (pihm->forc.riverbc[i].data[j]);
            }
            free (pihm->forc.riverbc[i].ftime);
            free (pihm->forc.riverbc[i].data);
        }
        free (pihm->forc.riverbc);
    }

    if (pihm->forc.nmeteo > 0)
    {
        for (i = 0; i < pihm->forc.nmeteo; i++)
        {
            for (j = 0; j < pihm->forc.meteo[i].length; j++)
            {
                free (pihm->forc.meteo[i].data[j]);
            }
            free (pihm->forc.meteo[i].ftime);
            free (pihm->forc.meteo[i].data);
            free (pihm->forc.meteo[i].value);
        }
        free (pihm->forc.meteo);
    }

    if (pihm->forc.nlai > 0)
    {
        for (i = 0; i < pihm->forc.nlai; i++)
        {
            for (j = 0; j < pihm->forc.lai[i].length; j++)
            {
                free (pihm->forc.lai[i].data[j]);
            }
            free (pihm->forc.lai[i].ftime);
            free (pihm->forc.lai[i].data);
            free (pihm->forc.lai[i].value);
        }
        free (pihm->forc.lai);
    }

#ifdef _NOAH_
    if (pihm->forc.nrad > 0)
    {
        for (i = 0; i < pihm->forc.nrad; i++)
        {
            for (j = 0; j < pihm->forc.rad[i].length; j++)
            {
                free (pihm->forc.rad[i].data[j]);
            }
            free (pihm->forc.rad[i].ftime);
            free (pihm->forc.rad[i].data);
            free (pihm->forc.rad[i].value);
        }
        free (pihm->forc.rad);
    }
#endif

    free (pihm->ctrl.tout);

    for (i = 0; i < pihm->ctrl.nprint; i++)
    {
        free (pihm->prtctrl[i].var);
        free (pihm->prtctrl[i].buffer);
    }

    free (pihm->elem);
    free (pihm->riv);
}
