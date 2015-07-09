
/*****************************************************************************
 * File		:   lsm_func.c 
 * Function	:   noah related functions
 ****************************************************************************/

#include "pihm.h"
#include "spa.h"
#include "noah.h"

void LsmRead (char *simulation, lsm_struct noah, pihm_struct pihm)
{
    int             i, j;
    char            fn[MAXSTRING];
    struct tm      *timeinfo;
    FILE           *lsm_file;
    FILE           *radn_file;
    FILE           *stream;
    char            cmdstr[MAXSTRING];
    char            buffer[MAXSTRING];
    int             index;
    int             match;
    char            project[MAXSTRING];
    char           *token;
    char            tempname[MAXSTRING];

    timeinfo = (struct tm *)malloc (sizeof (struct tm));

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


    /*
     * Open *.lsm file
     */
    sprintf (fn, "input/%s/%s.lsm", project, project);
    lsm_file = fopen (fn, "r");
    CheckFile (lsm_file, fn);

    /*
     * Start reading lsm_file
     */
    FindLine (lsm_file, "BOF");
    NextLine (lsm_file, cmdstr);
    ReadKeywordDouble (cmdstr, "LATITUDE", &noah->latitude);

    NextLine (lsm_file, cmdstr);
    ReadKeywordDouble (cmdstr, "LONGITUDE", &noah->longitude);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "NSOIL", &noah->std_nsoil);
    if (noah->std_nsoil > MAXLYR - 1)
    {
        printf
            ("Error: the number of soil layers should be smaller than %d\n!",
            MAXLYR - 1);
        exit (1);
    }

    NextLine (lsm_file, cmdstr);
    ReadKeywordStr (cmdstr, "SLDPTH_DATA", buffer);

    stream = fmemopen (buffer, strlen (buffer), "r");
    for (i = 0; i < noah->std_nsoil; i++)
    {
        fscanf (stream, "%lf", &noah->std_sldpth[i]);
    }
    fclose (stream);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "RAD_MODE_DATA", &noah->rad_mode);

    NextLine (lsm_file, cmdstr);
    ReadKeywordDouble (cmdstr, "SBETA_DATA", &noah->genprmt.sbeta_data);

    NextLine (lsm_file, cmdstr);
    ReadKeywordDouble (cmdstr, "FXEXP_DATA", &noah->genprmt.fxexp_data);

    NextLine (lsm_file, cmdstr);
    ReadKeywordDouble (cmdstr, "CSOIL_DATA", &noah->genprmt.csoil_data);

    NextLine (lsm_file, cmdstr);
    ReadKeywordDouble (cmdstr, "SALP_DATA", &noah->genprmt.salp_data);

    NextLine (lsm_file, cmdstr);
    ReadKeywordDouble (cmdstr, "FRZK_DATA", &noah->genprmt.frzk_data);

    NextLine (lsm_file, cmdstr);
    ReadKeywordDouble (cmdstr, "ZBOT_DATA", &noah->genprmt.zbot_data);

    NextLine (lsm_file, cmdstr);
    ReadKeywordDouble (cmdstr, "TBOT_DATA", &noah->genprmt.tbot_data);

    NextLine (lsm_file, cmdstr);
    ReadKeywordDouble (cmdstr, "CZIL_DATA", &noah->genprmt.czil_data);

    NextLine (lsm_file, cmdstr);
    ReadKeywordDouble (cmdstr, "LVCOEF_DATA", &noah->genprmt.lvcoef_data);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "T1", &noah->prtvrbl[T1_CTRL]);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "STC", &noah->prtvrbl[STC_CTRL]);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "SMC", &noah->prtvrbl[SMC_CTRL]);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "SH2O", &noah->prtvrbl[SH2O_CTRL]);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "SNOWH", &noah->prtvrbl[SNOWH_CTRL]);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "ALBEDO", &noah->prtvrbl[ALBEDO_CTRL]);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "LE", &noah->prtvrbl[LE_CTRL]);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "SH", &noah->prtvrbl[SH_CTRL]);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "G", &noah->prtvrbl[G_CTRL]);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "ETP", &noah->prtvrbl[ETP_CTRL]);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "ESNOW", &noah->prtvrbl[ESNOW_CTRL]);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "ROOTW", &noah->prtvrbl[ROOTW_CTRL]);

    NextLine (lsm_file, cmdstr);
    ReadKeywordInt (cmdstr, "SOILM", &noah->prtvrbl[SOILM_CTRL]);

    fclose (lsm_file);

    if (noah->rad_mode == 1)
    {
        sprintf (fn, "input/%s/%s.rad", project, project);
        radn_file = fopen (fn, "r");
        CheckFile (radn_file, fn);

        FindLine (radn_file, "BOF");
        NextLine (radn_file, cmdstr);
        match = sscanf (cmdstr, "%*s %d", &noah->forcing.nts);
        if (match != 1)
        {
            printf ("Cannot read number of radiation forcing time series!\n");
            printf (".rad file format error!\n");
            exit (1);
        }

        noah->forcing.ts =
            (ts_struct *)malloc (noah->forcing.nts * sizeof (ts_struct));

        for (i = 0; i < noah->forcing.nts; i++)
        {
            NextLine (radn_file, cmdstr);
            match = sscanf (cmdstr, "%*s %d", &index);
            if (match != 1 || i != index - 1)
            {
                printf
                    ("Cannot read information of the %dth forcing series!\n",
                    i);
                printf (".forc file format error!\n");
                exit (1);
            }

            /* Skip header lines */
            NextLine (radn_file, cmdstr);
            NextLine (radn_file, cmdstr);
            noah->forcing.ts[i].length = CountLine (radn_file, 1, "RAD_TS");
        }

        /* Rewind and read */
        FindLine (radn_file, "NUM_RAD_TS");
        for (i = 0; i < noah->forcing.nts; i++)
        {
            /* Skip header lines */
            NextLine (radn_file, cmdstr);
            NextLine (radn_file, cmdstr);
            NextLine (radn_file, cmdstr);

            noah->forcing.ts[i].ftime =
                (int *)malloc (noah->forcing.ts[i].length * sizeof (int));
            noah->forcing.ts[i].data =
                (double **)malloc (noah->forcing.ts[i].length *
                sizeof (double *));
            for (j = 0; j < noah->forcing.ts[i].length; j++)
            {
                noah->forcing.ts[i].data[j] =
                    (double *)malloc (2 * sizeof (double));
                NextLine (radn_file, cmdstr);
                ReadTS (cmdstr, &noah->forcing.ts[i].ftime[j],
                    &noah->forcing.ts[i].data[j][0], 2);
            }
        }

        fclose (radn_file);
    }
}

void LsmInitialize (char *simulation, pihm_struct pihm, lsm_struct noah)
{
    grid_struct    *grid;
    int             i, j, kz;
    double          zsoil[MAXLYR];
    double          frzfact;
    char            project[MAXSTRING];
    char           *token;
    char            tempname[MAXSTRING];
    elem_struct    *elem;

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

    if (noah->rad_mode == 1)
    {
        for (i = 0; i < 2; i++)
        {
            noah->forcing.radn[i] =
                (double *)malloc (noah->forcing.nts * sizeof (double));
        }
    }

    noah->grid = (grid_struct *) malloc (pihm->numele * sizeof (grid_struct));

    noah->ic.t1 = (double *)malloc (pihm->numele * sizeof (double));
    noah->ic.snowh = (double *)malloc (pihm->numele * sizeof (double));
    noah->ic.stc = (double **)malloc (pihm->numele * sizeof (double *));
    noah->ic.smc = (double **)malloc (pihm->numele * sizeof (double *));
    noah->ic.sh2o = (double **)malloc (pihm->numele * sizeof (double *));

    for (i = 0; i < pihm->numele; i++)
    {
        noah->ic.stc[i] = (double *)malloc (MAXLYR * sizeof (double));
        noah->ic.smc[i] = (double *)malloc (MAXLYR * sizeof (double));
        noah->ic.sh2o[i] = (double *)malloc (MAXLYR * sizeof (double));
    }

    for (i = 0; i < pihm->numele; i++)
    {
        grid = &noah->grid[i];
        elem = &pihm->elem[i];

        grid->radn[SOLAR_DIR_TS] =
            &noah->forcing.radn[SOLAR_DIR_TS][pihm->attrib_tbl.meteo[i] - 1];
        grid->radn[SOLAR_DIF_TS] =
            &noah->forcing.radn[SOLAR_DIF_TS][pihm->attrib_tbl.meteo[i] - 1];

        grid->avginfil = 0.0;
        for (j = 0; j < 3; j++)
        {
            grid->avgsubflux[j] = 0.0;
        }

        /* Set-up soil layer depths */
        DefSldpth (grid->sldpth, &grid->nsoil, elem->soil.depth,
            noah->std_sldpth, noah->std_nsoil);

        /* Set-up soil parameters */
        grid->csoil = noah->genprmt.csoil_data;
        grid->vgalpha = elem->soil.alpha;
        grid->vgbeta = elem->soil.beta;
        grid->smcmin = elem->soil.thetar;
        grid->macksat = elem->soil.kmacv;
        grid->areaf = elem->soil.areafh;
        grid->nmacd = FindLayer (grid->sldpth, grid->nsoil, elem->soil.dmac);

        grid->dksat = elem->soil.ksatv;
        grid->quartz = pihm->soil_tbl.qtz[pihm->attrib_tbl.soil[i] - 1];
        grid->smcdry = elem->soil.thetaw;
        grid->smcmax = elem->soil.thetas;
        grid->smcref = elem->soil.thetaref;
        grid->smcwlt = elem->soil.thetaw;

        /* Set-up universal parameters (not dependent on soil type or
         * vegetation type */
        grid->zbot = noah->genprmt.zbot_data;
        grid->tbot = noah->genprmt.tbot_data;
        grid->salp = noah->genprmt.salp_data;
        grid->sbeta = noah->genprmt.sbeta_data;
        grid->frzk = noah->genprmt.frzk_data;
        grid->fxexp = noah->genprmt.fxexp_data;
        grid->ptu = 0.0;        /* (not used yet) to satisfy intent (out) */
        grid->czil = noah->genprmt.czil_data;
        grid->lvcoef = noah->genprmt.lvcoef_data;

        /* To adjust frzk parameter to actual soil type */
        frzfact = (grid->smcmax / grid->smcref) * (0.412 / 0.468);
        grid->frzx = grid->frzk * frzfact;

        /* Set-up vegetation parameters */
        grid->vegtyp = elem->lc.type;
        grid->topt = elem->lc.topt;
        grid->cfactr = elem->lc.cfactr;
        grid->rsmax = elem->lc.rsmax;
        //grid->cmcmax = elem->
        grid->nroot = FindLayer (grid->sldpth, grid->nsoil, elem->lc.rzd);
        grid->snup = elem->lc.snup;
        grid->rsmin = elem->lc.rsmin;
        grid->rgl = elem->lc.rgl;
        grid->hs = elem->lc.hs;
        grid->emissmin = elem->lc.emissmin;
        grid->emissmax = elem->lc.emissmax;
        grid->laimin = elem->lc.laimin;
        grid->laimax = elem->lc.laimax;
        grid->z0max = elem->lc.z0max;
        grid->z0min = elem->lc.z0min;
        grid->albedomin = elem->lc.albedomin;
        grid->albedomax = elem->lc.albedomax;
        grid->isurban = 13;
        grid->cmcfactr = elem->lc.intcp_factr;

        if (grid->vegtyp == pihm->lc_tbl.bare)
        {
            grid->shdfac = 0.0;
        }

        /* Calculate root distribution.
         * Present version assumes uniform distribution based on soil layer
         * depths. */

        zsoil[0] = -grid->sldpth[0];
        for (kz = 1; kz < grid->nsoil; kz++)
        {
            zsoil[kz] = -grid->sldpth[kz] + zsoil[kz - 1];
        }

        for (j = 0; j < grid->nroot; j++)
        {
            grid->rtdis[j] = -grid->sldpth[j] / zsoil[grid->nroot - 1];
        }

        /* Apply calibration */
        grid->czil *= pihm->cal.czil;
        grid->fxexp *= pihm->cal.fxexp;
        grid->rsmin *= pihm->cal.rsmin;
        grid->rgl *= pihm->cal.rgl;
        grid->hs *= pihm->cal.hs;
        grid->cfactr *= pihm->cal.cfactr;
        grid->cmcfactr *= pihm->cal.intcp;
        grid->smcref =
            (grid->smcref - grid->smcmin) * pihm->cal.thetaref + grid->smcmin;
        grid->smcwlt =
            (grid->smcwlt - grid->smcmin) * pihm->cal.thetaw + grid->smcmin;

        /* Initialize topographic radiation related parameters */
        if (noah->rad_mode == 1)
        {
            CalcSlopeAspect (grid, *elem, pihm);
        }

        grid->snotime1 = 0.0;
        grid->ribb = 0.0;

        grid->sheat = BADVAL;
        grid->eta_kinematic = BADVAL;
        grid->eta = BADVAL;
        grid->fdown = BADVAL;
        grid->ec = BADVAL;
        grid->edir = BADVAL;
        grid->ett = BADVAL;
        grid->esnow = BADVAL;
        grid->drip = BADVAL;
        grid->dew = BADVAL;
        grid->beta = BADVAL;
        grid->t1 = BADVAL;
        grid->snowh = BADVAL;
        grid->sneqv = BADVAL;
        grid->etp = BADVAL;
        grid->ssoil = BADVAL;
        grid->flx1 = BADVAL;
        grid->flx2 = BADVAL;
        grid->flx3 = BADVAL;
        grid->snomlt = BADVAL;
        grid->sncovr = BADVAL;
        grid->runoff1 = BADVAL;
        grid->runoff2 = BADVAL;
        grid->runoff3 = BADVAL;
        grid->rc = BADVAL;
        grid->pc = BADVAL;
        grid->rcs = BADVAL;
        grid->rct = BADVAL;
        grid->rcsoil = BADVAL;
        grid->soilw = BADVAL;
        grid->soilm = BADVAL;
        grid->q1 = BADVAL;
        //grid->smcwlt = BADVAL;
        //grid->smcdry = BADVAL;
        //grid->smcref = BADVAL;
        //grid->smcmax = BADVAL;
        //grid->rsmin = BADVAL;
        //grid->nroot = BADVAL;

        grid->pcpdrp = 0.0;
        grid->drip = 0.0;

        grid->dt = pihm->ctrl.etstep;

        grid->snoalb = 0.75;
        grid->zlvl = 3.0;
        grid->zlvl_wind = elem->forc.zlvl_wind;
        grid->isurban = 13;
        grid->shdmin = 0.01;
        grid->shdmax = 0.96;

        grid->usemonalb = 0;
        grid->rdlai2d = 0;
        grid->iz0tlnd = 0;

        grid->cosz = BADVAL;
        grid->prcprain = BADVAL;
        grid->solardirect = BADVAL;

        grid->emissi = 0.96;
        grid->albedo = 0.18;
        grid->z0 = 0.1;

        grid->z0brd = BADVAL;

        grid->ch = 1.0e-4;
        grid->cm = 1.0e-4;
    }

    /* Set initial conditions for land surface variables */
    if (pihm->ctrl.init_type < 3)       /* Relaxation */
    {
        ApplyForcing (&pihm->forcing, pihm->ctrl.starttime);
        LsmSaturationIC (&noah->ic, noah->grid, pihm->elem, pihm->numele);
    }
    else
    {
        /* Hot start mode */
        ReadLsmInit (project, simulation, &noah->ic, pihm->numele);
    }

    InitLsmVrbl (noah->grid, pihm->elem, pihm->numele, noah->ic);
}

void MapLsmOutput (char *simulation, lsm_struct noah, int numele,
    char *outputdir)
{
    int             i, j, k;
    int             n;

    if (verbose_mode)
        printf ("\nInitializing output files\n");

    n = 0;

    for (i = 0; i < NUM_PRINT; i++)
    {
        if (noah->prtvrbl[i] > 0)
        {
            switch (i)
            {
                case T1_CTRL:
                    sprintf (noah->prtctrl[n].name, "%s%s.t1", outputdir,
                        simulation);
                    noah->prtctrl[n].intvl = noah->prtvrbl[i];
                    noah->prtctrl[n].nvrbl = numele;
                    noah->prtctrl[n].vrbl =
                        (double **)malloc (noah->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < numele; j++)
                    {
                        noah->prtctrl[n].vrbl[j] = &noah->grid[j].t1;
                    }
                    n++;
                    break;
                case STC_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf (noah->prtctrl[n].name, "%s%s.stc%d",
                            outputdir, simulation, k);
                        noah->prtctrl[n].intvl = noah->prtvrbl[i];
                        noah->prtctrl[n].nvrbl = numele;
                        noah->prtctrl[n].vrbl =
                            (double **)malloc (noah->prtctrl[n].nvrbl *
                            sizeof (double *));
                        for (j = 0; j < numele; j++)
                        {
                            noah->prtctrl[n].vrbl[j] = &noah->grid[j].stc[k];
                        }
                        n++;
                    }
                    break;
                case SMC_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf (noah->prtctrl[n].name, "%s%s.smc%d",
                            outputdir, simulation, k);
                        noah->prtctrl[n].intvl = noah->prtvrbl[i];
                        noah->prtctrl[n].nvrbl = numele;
                        noah->prtctrl[n].vrbl =
                            (double **)malloc (noah->prtctrl[n].nvrbl *
                            sizeof (double *));
                        for (j = 0; j < numele; j++)
                        {
                            noah->prtctrl[n].vrbl[j] = &noah->grid[j].smc[k];
                        }
                        n++;
                    }
                    break;
                case SH2O_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf (noah->prtctrl[n].name, "%s%s.swc%d",
                            outputdir, simulation, k);
                        noah->prtctrl[n].intvl = noah->prtvrbl[i];
                        noah->prtctrl[n].nvrbl = numele;
                        noah->prtctrl[n].vrbl =
                            (double **)malloc (noah->prtctrl[n].nvrbl *
                            sizeof (double *));
                        for (j = 0; j < numele; j++)
                        {
                            noah->prtctrl[n].vrbl[j] = &noah->grid[j].sh2o[k];
                        }
                        n++;
                    }
                    break;
                case SNOWH_CTRL:
                    sprintf (noah->prtctrl[n].name, "%s%s.snowh", outputdir,
                        simulation);
                    noah->prtctrl[n].intvl = noah->prtvrbl[i];
                    noah->prtctrl[n].nvrbl = numele;
                    noah->prtctrl[n].vrbl =
                        (double **)malloc (noah->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < numele; j++)
                    {
                        noah->prtctrl[n].vrbl[j] = &noah->grid[j].snowh;
                    }
                    n++;
                    break;
                case ALBEDO_CTRL:
                    sprintf (noah->prtctrl[n].name, "%s%s.albedo", outputdir,
                        simulation);
                    noah->prtctrl[n].intvl = noah->prtvrbl[i];
                    noah->prtctrl[n].nvrbl = numele;
                    noah->prtctrl[n].vrbl =
                        (double **)malloc (noah->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < numele; j++)
                    {
                        noah->prtctrl[n].vrbl[j] = &noah->grid[j].albedo;
                    }
                    n++;
                    break;
                case LE_CTRL:
                    sprintf (noah->prtctrl[n].name, "%s%s.le", outputdir,
                        simulation);
                    noah->prtctrl[n].intvl = noah->prtvrbl[i];
                    noah->prtctrl[n].nvrbl = numele;
                    noah->prtctrl[n].vrbl =
                        (double **)malloc (noah->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < numele; j++)
                    {
                        noah->prtctrl[n].vrbl[j] = &noah->grid[j].eta;
                    }
                    n++;
                    break;
                case SH_CTRL:
                    sprintf (noah->prtctrl[n].name, "%s%s.sh", outputdir,
                        simulation);
                    noah->prtctrl[n].intvl = noah->prtvrbl[i];
                    noah->prtctrl[n].nvrbl = numele;
                    noah->prtctrl[n].vrbl =
                        (double **)malloc (noah->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < numele; j++)
                    {
                        noah->prtctrl[n].vrbl[j] = &noah->grid[j].sheat;
                    }
                    n++;
                    break;
                case G_CTRL:
                    sprintf (noah->prtctrl[n].name, "%s%s.g", outputdir,
                        simulation);
                    noah->prtctrl[n].intvl = noah->prtvrbl[i];
                    noah->prtctrl[n].nvrbl = numele;
                    noah->prtctrl[n].vrbl =
                        (double **)malloc (noah->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < numele; j++)
                    {
                        noah->prtctrl[n].vrbl[j] = &noah->grid[j].ssoil;
                    }
                    n++;
                    break;
                case ETP_CTRL:
                    sprintf (noah->prtctrl[n].name, "%s%s.etp", outputdir,
                        simulation);
                    noah->prtctrl[n].intvl = noah->prtvrbl[i];
                    noah->prtctrl[n].nvrbl = numele;
                    noah->prtctrl[n].vrbl =
                        (double **)malloc (noah->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < numele; j++)
                    {
                        noah->prtctrl[n].vrbl[j] = &noah->grid[j].etp;
                    }
                    n++;
                    break;
                case ESNOW_CTRL:
                    sprintf (noah->prtctrl[n].name, "%s%s.esnow", outputdir,
                        simulation);
                    noah->prtctrl[n].intvl = noah->prtvrbl[i];
                    noah->prtctrl[n].nvrbl = numele;
                    noah->prtctrl[n].vrbl =
                        (double **)malloc (noah->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < numele; j++)
                    {
                        noah->prtctrl[n].vrbl[j] = &noah->grid[j].esnow;
                    }
                    n++;
                    break;
                case ROOTW_CTRL:
                    sprintf (noah->prtctrl[n].name, "%s%s.rootw", outputdir,
                        simulation);
                    noah->prtctrl[n].intvl = noah->prtvrbl[i];
                    noah->prtctrl[n].nvrbl = numele;
                    noah->prtctrl[n].vrbl =
                        (double **)malloc (noah->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < numele; j++)
                    {
                        noah->prtctrl[n].vrbl[j] = &noah->grid[j].soilw;
                    }
                    n++;
                    break;
                case SOILM_CTRL:
                    sprintf (noah->prtctrl[n].name, "%s%s.soilm", outputdir,
                        simulation);
                    noah->prtctrl[n].intvl = noah->prtvrbl[i];
                    noah->prtctrl[n].nvrbl = numele;
                    noah->prtctrl[n].vrbl =
                        (double **)malloc (noah->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < numele; j++)
                    {
                        noah->prtctrl[n].vrbl[j] = &noah->grid[j].soilm;
                    }
                    n++;
                    break;
            }
        }
    }

    noah->nprint = n;
}

//void LSM_initialize_output (char *filename, Model_Data PIHM, Control_Data CS, LSM_STRUCT noah, char *outputdir)
//{
//    FILE           *Ofile;
//    char           *ascii_name;
//    int             i, j, ensemble_mode, icounter;
//
//    if (strstr (filename, ".") != 0)
//        ensemble_mode = 1;
//    else
//        ensemble_mode = 0;
//
//    if (CS->Verbose)
//        printf ("\nInitializing noah output files ...\n");
//
//    icounter = 0;
//    if (noah->PRINT_T1 > 0)
//    {
//        sprintf (noah->PCtrl[icounter].name, "%s%s.Tsfc", outputdir, filename);
//        noah->PCtrl[icounter].Interval = noah->PRINT_T1;
//        noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//        noah->PCtrl[icounter].PrintVar =
//           (double **)malloc (noah->PCtrl[icounter].NumVar *
//           sizeof (double *));
//        for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//            noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].T1);
//        icounter++;
//    }
//    if (noah->PRINT_STC > 0)
//    {
//        for (j = 0; j < noah->STD_NSOIL + 1; j++)
//        {
//            sprintf (noah->PCtrl[icounter].name, "%s%s.TSOIL%d", outputdir,
//               filename, j);
//            noah->PCtrl[icounter].Interval = noah->PRINT_STC;
//            noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//            noah->PCtrl[icounter].PrintVar =
//               (double **)malloc (noah->PCtrl[icounter].NumVar *
//               sizeof (double *));
//            for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//            noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].STC[j]);
//            icounter++;
//        }
//    }
//    if (noah->PRINT_SMC > 0)
//    {
//        for (j = 0; j < noah->STD_NSOIL + 1; j++)
//        {
//            sprintf (noah->PCtrl[icounter].name, "%s%s.SM%d", outputdir,
//               filename, j);
//            noah->PCtrl[icounter].Interval = noah->PRINT_SMC;
//            noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//            noah->PCtrl[icounter].PrintVar =
//               (double **)malloc (noah->PCtrl[icounter].NumVar *
//               sizeof (double *));
//            for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//                noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].SMC[j]);
//            icounter++;
//        }
//    }
//    if (noah->PRINT_SH2O > 0)
//    {
//        for (j = 0; j < noah->STD_NSOIL + 1; j++)
//        {
//            sprintf (noah->PCtrl[icounter].name, "%s%s.SW%d", outputdir,
//               filename, j);
//            noah->PCtrl[icounter].Interval = noah->PRINT_SH2O;
//            noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//            noah->PCtrl[icounter].PrintVar =
//               (double **)malloc (noah->PCtrl[icounter].NumVar *
//               sizeof (double *));
//            for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//                noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].SH2O[j]);
//            icounter++;
//        }
//    }
//    if (noah->PRINT_SNOWH > 0)
//    {
//        sprintf (noah->PCtrl[icounter].name, "%s%s.snowH", outputdir,
//           filename);
//        noah->PCtrl[icounter].Interval = noah->PRINT_SNOWH;
//        noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//        noah->PCtrl[icounter].PrintVar =
//           (double **)malloc (noah->PCtrl[icounter].NumVar *
//           sizeof (double *));
//        for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//            noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].SNOWH);
//        icounter++;
//    }
//    if (noah->PRINT_ALBEDO > 0)
//    {
//        sprintf (noah->PCtrl[icounter].name, "%s%s.Albedo", outputdir,
//           filename);
//        noah->PCtrl[icounter].Interval = noah->PRINT_ALBEDO;
//        noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//        noah->PCtrl[icounter].PrintVar =
//           (double **)malloc (noah->PCtrl[icounter].NumVar *
//           sizeof (double *));
//        for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//            noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].ALBEDO);
//        icounter++;
//    }
//    if (noah->PRINT_LE > 0)
//    {
//        sprintf (noah->PCtrl[icounter].name, "%s%s.LE", outputdir, filename);
//        noah->PCtrl[icounter].Interval = noah->PRINT_LE;
//        noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//        noah->PCtrl[icounter].PrintVar =
//           (double **)malloc (noah->PCtrl[icounter].NumVar *
//           sizeof (double *));
//        for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//            noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].ETA);
//        icounter++;
//    }
//    if (noah->PRINT_SH > 0)
//    {
//        sprintf (noah->PCtrl[icounter].name, "%s%s.SH", outputdir, filename);
//        noah->PCtrl[icounter].Interval = noah->PRINT_SH;
//        noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//        noah->PCtrl[icounter].PrintVar =
//           (double **)malloc (noah->PCtrl[icounter].NumVar *
//           sizeof (double *));
//        for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//            noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].SHEAT);
//        icounter++;
//    }
//    if (noah->PRINT_G > 0)
//    {
//        sprintf (noah->PCtrl[icounter].name, "%s%s.G", outputdir, filename);
//        noah->PCtrl[icounter].Interval = noah->PRINT_G;
//        noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//        noah->PCtrl[icounter].PrintVar =
//           (double **)malloc (noah->PCtrl[icounter].NumVar *
//           sizeof (double *));
//        for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//            noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].SSOIL);
//        icounter++;
//    }
//    if (noah->PRINT_ETP > 0)
//    {
//        sprintf (noah->PCtrl[icounter].name, "%s%s.ETP", outputdir, filename);
//        noah->PCtrl[icounter].Interval = noah->PRINT_ETP;
//        noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//        noah->PCtrl[icounter].PrintVar =
//           (double **)malloc (noah->PCtrl[icounter].NumVar *
//           sizeof (double *));
//        for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//            noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].ETP);
//        icounter++;
//    }
//    if (noah->PRINT_ESNOW > 0)
//    {
//        sprintf (noah->PCtrl[icounter].name, "%s%s.ESNOW", outputdir, filename);
//        noah->PCtrl[icounter].Interval = noah->PRINT_ESNOW;
//        noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//        noah->PCtrl[icounter].PrintVar =
//           (double **)malloc (noah->PCtrl[icounter].NumVar *
//           sizeof (double *));
//        for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//            noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].ESNOW);
//        icounter++;
//    }
//    if (noah->PRINT_ROOTW > 0)
//    {
//        sprintf (noah->PCtrl[icounter].name, "%s%s.ROOTW", outputdir, filename);
//        noah->PCtrl[icounter].Interval = noah->PRINT_ROOTW;
//        noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//        noah->PCtrl[icounter].PrintVar =
//           (double **)malloc (noah->PCtrl[icounter].NumVar *
//           sizeof (double *));
//        for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//            noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].SOILW);
//        icounter++;
//    }
//    if (noah->PRINT_SOILM > 0)
//    {
//        sprintf (noah->PCtrl[icounter].name, "%s%s.SOILM", outputdir, filename);
//        noah->PCtrl[icounter].Interval = noah->PRINT_SOILM;
//        noah->PCtrl[icounter].NumVar = PIHM->NumEle;
//        noah->PCtrl[icounter].PrintVar =
//           (double **)malloc (noah->PCtrl[icounter].NumVar *
//           sizeof (double *));
//        for (i = 0; i < noah->PCtrl[icounter].NumVar; i++)
//            noah->PCtrl[icounter].PrintVar[i] = &(noah->GRID[i].SOILM);
//        icounter++;
//    }
//    noah->NPRINT = icounter;
//
//    for (i = 0; i < noah->NPRINT; i++)
//    {
//        if (noah->PCtrl[i].Interval < noah->GRID[0].DT)
//        {
//            printf("\nLSM %s print interval must not be smaller than noah step size!\n", noah->PCtrl[i].name);
//            exit(1);
//        }
//        Ofile = fopen (noah->PCtrl[i].name, "w");
//        fclose (Ofile);
//
//      if (CS->Ascii)
//      {
//          ascii_name = (char *)malloc ((strlen (noah->PCtrl[i].name) + 5) * sizeof (char));
//          sprintf (ascii_name, "%s.txt", noah->PCtrl[i].name);
//          Ofile = fopen (ascii_name, "w");
//          fclose (Ofile);
//            free (ascii_name);
//      }
//
//        noah->PCtrl[i].buffer = (double *)calloc (noah->PCtrl[i].NumVar, sizeof (double));
//    }
//}
//
//void LSM_FreeData (Model_Data PIHM, LSM_STRUCT noah)
//{
//    GRID_TYPE      *NOAH;
//
//    int             i;
//    for (i = 0; i < PIHM->NumEle; i++)
//    {
//        NOAH = &(noah->GRID[i]);
//        free (NOAH->STC);
//        free (NOAH->SMC);
//        free (NOAH->SH2O);
//        free (NOAH->SLDPTH);
//        free (NOAH->ET);
//        free (NOAH->SMAV);
//        free (NOAH->RTDIS);
//    }
//    for (i = 0; i < noah->NPRINT; i++)
//    {
//        free (noah->PCtrl[i].PrintVar);
//        free (noah->PCtrl[i].buffer);
//    }
//    free (noah->GRID);
//    free (noah->STD_SLDPTH);
//}
//
//void LSM_PrintInit (Model_Data PIHM, LSM_STRUCT noah, char *filename)
//{
//    FILE           *init_file;
//    char           *init_name;
//    int             i, j;
//    init_name = (char *)malloc ((2 * strlen (filename) + 16) * sizeof (char));
//    sprintf (init_name, "input/%s/%s.lsminit", filename, filename);
//    init_file = fopen (init_name, "w");
//    free (init_name);
//
//    for (i = 0; i < PIHM->NumEle; i++)
//    {
//        fprintf (init_file, "%lf\t%lf", noah->GRID[i].T1, noah->GRID[i].SNOWH);
//        for (j = 0; j < noah->STD_NSOIL + 1; j++)
//            fprintf (init_file, "\t%lf", noah->GRID[i].STC[j]);
//        for (j = 0; j < noah->STD_NSOIL + 1; j++)
//            fprintf (init_file, "\t%lf", noah->GRID[i].SMC[j]);
//        for (j = 0; j < noah->STD_NSOIL + 1; j++)
//            fprintf (init_file, "\t%lf", noah->GRID[i].SH2O[j]);
//        fprintf (init_file, "\n");
//    }
//
//    fclose (init_file);
//}
//
int FindLayer (const double *sldpth, int nsoil, double depth)
{
    int             layer;
    int             j = 0, ind = 0;
    double          dsum = 0.0;

    if (depth <= 0.0)
    {
        layer = 0;
    }
    else
    {
        while (dsum < depth)
        {
            if (sldpth[j] < 0.0)
            {
                break;
            }
            dsum += sldpth[j];
            ind = j;
            j++;
        }
        layer = ind + 1;
        layer = (layer > nsoil) ? nsoil : layer;
    }
    return (layer);
}

double mod (double a, double N)
{
    return (a - N * floor (a / N));
}

double TopoRadiation (double sdir, double sdif, double zenith,
    double azimuth180, double slope, double aspect, double *h_phi, double svf)
{
    double          incidence;
    double          gvf;
    double          soldown;

    if (zenith > h_phi[(int)floor (azimuth180 / 10.0)])
        sdir = 0.0;
    incidence =
        180.0 / PI * acos (cos (zenith * PI / 180.0) * cos (slope * PI /
            180.0) +
        sin (zenith * PI / 180.0) * sin (slope * PI / 180.0) *
        cos ((azimuth180 - aspect) * PI / 180.0));
    incidence = incidence > 90.0 ? 90.0 : incidence;
    gvf = (1.0 + cos (slope * PI / 180.0)) / 2.0 - svf;
    gvf = gvf < 0.0 ? 0.0 : gvf;
    soldown =
        sdir * cos (incidence * PI / 180.0) + svf * sdif +
        0.2 * gvf * (sdir * cos (zenith * PI / 180.0) + sdif);
    soldown = soldown < 0.0 ? 0.0 : soldown;

    return (soldown);
}

void DefSldpth (double *sldpth, int *nsoil, double total_depth,
    double *std_sldpth, int std_nsoil)
{
    int             j, k;
    double          zsoil[MAXLYR];

    zsoil[0] = std_sldpth[0];

    for (j = 1; j < MAXLYR; j++)
    {
        zsoil[j] = zsoil[j - 1] + std_sldpth[j];
    }

    if (total_depth <= zsoil[0])
    {
        sldpth[0] = total_depth;
        *nsoil = 1;
        for (j = 1; j < MAXLYR; j++)
        {
            sldpth[j] = BADVAL;
        }
    }
    else if (total_depth <= zsoil[std_nsoil - 1])
    {
        for (j = 1; j < std_nsoil + 1; j++)
        {
            if (total_depth <= zsoil[j])
            {
                for (k = 0; k < j; k++)
                {
                    sldpth[k] = std_sldpth[k];
                }
                sldpth[j] = total_depth - zsoil[j - 1];
                *nsoil = j + 1;

                /* The following calculations gurantee that each layer is
                 * thicker than the layer on top */
                if (sldpth[j] < sldpth[j - 1])
                {
                    sldpth[j - 1] += sldpth[j];
                    sldpth[j] = BADVAL;
                    *nsoil -= 1;
                }
                for (k = j + 1; k < MAXLYR; k++)
                {
                    sldpth[k] = BADVAL;
                }
                break;
            }
        }
    }
    else
    {
        for (j = 0; j < std_nsoil; j++)
        {
            sldpth[j] = std_sldpth[j];
        }
        sldpth[std_nsoil] = total_depth - zsoil[std_nsoil - 1];
        *nsoil = std_nsoil + 1;
        if (sldpth[std_nsoil] < sldpth[std_nsoil - 1])
        {
            sldpth[std_nsoil - 1] += sldpth[std_nsoil];
            sldpth[std_nsoil] = BADVAL;
            *nsoil -= 1;
        }
    }
}

void CalcSlopeAspect (grid_struct * grid, elem_struct elem, pihm_struct pihm)
{
    double          x[3];
    double          y[3];
    double          zmax[3];
    double          vector1[3], vector2[3], normal_vector[3], vector[3], h, c,
        se, ce;
    int nodes[2];
    double          x1, y1, z1, x2, y2, z2, xc, yc, zc;
    double          c1, c2, ce1, ce2, se1, se2, phi1, phi2;
    int             ind, ind1, ind2;
    int             j, k;

    for (j = 0; j < 3; j++)
    {
        x[j] = pihm->mesh_tbl.x[elem.node[j] - 1];
        y[j] = pihm->mesh_tbl.y[elem.node[j] - 1];
        zmax[j] = pihm->mesh_tbl.zmax[elem.node[j] - 1];
    }

    vector1[0] = x[0] - x[2];
    vector1[1] = y[0] - y[2];
    vector1[2] = zmax[0] - zmax[2];

    vector2[0] = x[1] - x[2];
    vector2[1] = y[1] - y[2];
    vector2[2] = zmax[1] - zmax[2];

    /* Calculate normal vector */
    normal_vector[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1];
    normal_vector[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2];
    normal_vector[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0];

    if (normal_vector[2] < 0.0)
    {
        normal_vector[0] = -normal_vector[0];
        normal_vector[1] = -normal_vector[1];
        normal_vector[2] = -normal_vector[2];
    }

    /* Calculate slope */
    c = sqrt (normal_vector[0] * normal_vector[0] +
        normal_vector[1] * normal_vector[1]);
    grid->slope = atan (c / normal_vector[2]) * 180.0 / PI;

    /* Calculte aspect */
    ce = normal_vector[0] / c;
    se = normal_vector[1] / c;
    grid->aspect = acos (ce) * 180.0 / PI;

    if (se < 0.0)
        grid->aspect = 360.0 - grid->aspect;

    grid->aspect = mod (360.0 - grid->aspect + 270.0, 360.0);

    /* Calculate sky view factor (Dozier and Frew 1990) */
    grid->svf = 0.0;

    for (j = 0; j < 36; j++)
    {
        grid->h_phi[j] = 90.0;
    }

    for (j = 0; j < pihm->numele; j++)
    {
        for (k = 0; k < 3; k++)
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
            x1 = pihm->mesh_tbl.x[pihm->elem[j].node[nodes[0]] - 1];
            y1 = pihm->mesh_tbl.y[pihm->elem[j].node[nodes[0]] - 1];
            z1 = pihm->mesh_tbl.zmax[pihm->elem[j].node[nodes[0]] - 1];
            x2 = pihm->mesh_tbl.x[pihm->elem[j].node[nodes[1]] - 1];
            y2 = pihm->mesh_tbl.y[pihm->elem[j].node[nodes[1]] - 1];
            z2 = pihm->mesh_tbl.zmax[pihm->elem[j].node[nodes[1]] - 1];

            xc = 0.5 * (x1 + x2);
            yc = 0.5 * (y1 + y2);
            zc = 0.5 * (z1 + z2);

            vector[0] = xc - elem.topo.x;
            vector[1] = yc - elem.topo.y;
            vector[2] = zc - elem.topo.zmax;
            c = sqrt (vector[0] * vector[0] + vector[1] * vector[1]);
            h = atan (c / vector[2]) * 180.0 / PI;
            h = (h < 0.0) ? 90.0 : h;

            vector1[0] = x1 - elem.topo.x;
            vector1[1] = y1 - elem.topo.y;
            vector1[2] = z1 - elem.topo.zmax;
            vector2[0] = x2 - elem.topo.x;
            vector2[1] = y2 - elem.topo.y;
            vector2[2] = z2 - elem.topo.zmax;

            c1 = sqrt (vector1[0] * vector1[0] + vector1[1] * vector1[1]);
            c2 = sqrt (vector2[0] * vector2[0] + vector2[1] * vector2[1]);

            ce1 = vector1[0] / c1;
            se1 = vector1[1] / c1;
            phi1 = acos (ce1) * 180.0 / PI;
            if (se1 < 0.0)
            {
                phi1 = 360.0 - phi1;
            }
            phi1 = mod (360.0 - phi1 + 270.0, 360.0);

            ce2 = vector2[0] / c2;
            se2 = vector2[1] / c2;
            phi2 = acos (ce2) * 180.0 / PI;
            if (se2 < 0.0)
            {
                phi2 = 360.0 - phi2;
            }
            phi2 = mod (360.0 - phi2 + 270.0, 360.0);

            if (fabs (phi1 - phi2) > 180.0)
            {
                ind1 = 0;
                ind2 = (int)floor ((phi1 < phi2 ? phi1 : phi2) / 10.0);
                for (ind = ind1; ind <= ind2; ind++)
                {
                    if (h < grid->h_phi[ind])
                    {
                        grid->h_phi[ind] = h;
                    }
                }

                ind1 = (int)floor ((phi1 > phi2 ? phi1 : phi2) / 10.0);
                ind2 = 35;
                for (ind = ind1; ind <= ind2; ind++)
                {
                    if (h < grid->h_phi[ind])
                    {
                        grid->h_phi[ind] = h;
                    }
                }
            }
            else
            {
                ind1 = (int)floor ((phi1 < phi2 ? phi1 : phi2) / 10.0);
                ind2 = (int)floor ((phi1 > phi2 ? phi1 : phi2) / 10.0);
                for (ind = ind1; ind <= ind2; ind++)
                {
                    if (h < grid->h_phi[ind])
                    {
                        grid->h_phi[ind] = h;
                    }
                }
            }
        }
    }

    for (ind = 0; ind < 36; ind++)
    {
        grid->svf +=
            0.5 / PI * (cos (grid->slope * PI / 180.0) *
            (pow (sin (grid->h_phi[ind] * PI / 180.0),
                    2)) + sin (grid->slope * PI / 180.0) * cos ((ind * 10.0 +
                    5.0 -
                    grid->aspect) * PI / 180.0) * grid->h_phi[ind] * PI /
            180.0 -
            sin (grid->h_phi[ind] * PI / 180.0) * cos (grid->h_phi[ind] * PI /
                180.0)) * 10.0 / 180.0 * PI;
    }

    if (verbose_mode)
    {
        printf ("ele: slope = %lf, aspect = %lf, svf = %lf\t", grid->slope,
            grid->aspect, grid->svf);
        for (ind = 0; ind < 36; ind++)
        {
            printf ("%lf\t", grid->h_phi[ind]);
        }
        printf ("\n");
    }
}

void LsmSaturationIC (lsm_ic_struct *ic, const grid_struct * grid,
    const elem_struct *elem, int numele)
{
    double          sfctmp;
    int             i, j;

    for (i = 0; i < numele; i++)
    {
        sfctmp = *elem[i].forc.meteo[SFCTMP_TS];

        ic->t1[i] = sfctmp;

        ic->stc[i][0] =
            sfctmp + (sfctmp -
            grid[i].tbot) / grid[i].zbot * grid[i].sldpth[0] * 0.5;

        for (j = 1; j < MAXLYR; j++)
        {
            if (grid[i].sldpth[j] > 0)
            {
                ic->stc[i][j] =
                    ic->stc[i][j - 1] + (sfctmp -
                    grid[i].tbot) / grid[i].zbot * (grid[i].sldpth[j - 1] +
                    grid[i].sldpth[j]) * 0.5;
            }
            else
            {
                ic->stc[i][j] = BADVAL;
            }
        }

        for (j = 0; j < MAXLYR; j++)
        {
            if (grid[i].sldpth[j] > 0.0)
            {
                ic->smc[i][j] = elem[i].soil.thetas;
                ic->sh2o[i][j] = elem[i].soil.thetas;
            }
            else
            {
                ic->smc[i][j] = BADVAL;
                ic->sh2o[i][j] = BADVAL;
            }
        }
        ic->snowh[i] = 0.0;
    }
}

void ReadLsmInit (char *project, char *simulation, lsm_ic_struct *ic,
    int numele)
{
    char            fn[MAXSTRING];
    FILE           *init_file;
    int             i, j;

    sprintf (fn, "input/%s/%s.lsminit", project, simulation);
    init_file = fopen (fn, "r");
    CheckFile (init_file, fn);

    for (i = 0; i < numele; i++)
    {
        fscanf (init_file, "%lf %lf", &ic->t1[i], &ic->snowh[i]);
        for (j = 0; j < MAXLYR; j++)
        {
            fscanf (init_file, "%lf", &ic->stc[i][j]);
        }
        for (j = 0; j < MAXLYR; j++)
        {
            fscanf (init_file, "%lf", &ic->smc[i][j]);
        }
        for (j = 0; j < MAXLYR; j++)
        {
            fscanf (init_file, "%lf", &ic->sh2o[i][j]);
        }
    }
}

void InitLsmVrbl (grid_struct * grid, elem_struct *elem, int numele,
    lsm_ic_struct ic)
{
    int             i, j;

    for (i = 0; i < numele; i++)
    {
        grid[i].t1 = ic.t1[i];
        grid[i].snowh = ic.snowh[i];

        for (j = 0; j < MAXLYR; j++)
        {
            grid[i].stc[j] = ic.stc[i][j];
            grid[i].smc[j] = ic.smc[i][j];
            grid[i].sh2o[j] = ic.sh2o[i][j];
        }

        grid[i].cmc = elem[i].intcp;
        grid[i].sneqv = elem[i].snow;
    }
}
