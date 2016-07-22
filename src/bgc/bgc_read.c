#include "pihm.h"

void ReadBGC (char *fn, ctrl_struct *ctrl, co2control_struct *co2,
    ndepcontrol_struct *ndepctrl, char *co2_fn, char *ndep_fn)
{
    FILE           *bgc_file;
    struct tm      *timestamp;
    time_t          rawtime;
    char            cmdstr[MAXSTRING];

    timestamp = (struct tm *)malloc (sizeof (struct tm));

    /* Read bgc simulation control file */
    bgc_file = fopen (fn, "r");

    if (NULL == bgc_file)
    {
        fprintf (stderr, "Error opening %s.\n", fn);
        PIHMError (1);
    }

    FindLine (bgc_file, "TIME_DEFINE");
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%d", &ctrl->spinupstartyear);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%d", &ctrl->spinupendyear);

    timestamp->tm_year = ctrl->spinupstartyear - 1900;
    timestamp->tm_mon = 1 - 1;
    timestamp->tm_mday = 1;
    timestamp->tm_hour = 0;
    timestamp->tm_min = 0;
    timestamp->tm_sec = 0;
    rawtime = timegm (timestamp);
    ctrl->spinupstart = (int)rawtime;

    timestamp->tm_year = ctrl->spinupendyear + 1 - 1900;
    rawtime = timegm (timestamp);
    ctrl->spinupend = (int)rawtime;

    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%d", &ctrl->bgc_spinup);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%d", &ctrl->maxspinyears);

    FindLine (bgc_file, "CO2_CONTROL");
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%d", &co2->varco2);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &co2->co2ppm);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%s", co2_fn);

    FindLine (bgc_file, "NDEP_CONTROL");
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%d", &ndepctrl->varndep);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ndepctrl->ndep);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ndepctrl->nfix);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%s", ndep_fn);

    FindLine (bgc_file, "C_STATE");
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->cinit.max_leafc);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->cinit.max_stemc);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->cs.cwdc);
    ctrl->ns.cwdn = BADVAL;
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->cs.litr1c);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->cs.litr2c);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->cs.litr3c);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->cs.litr4c);
    ctrl->ns.litr2n = BADVAL;
    ctrl->ns.litr3n = BADVAL;
    ctrl->ns.litr4n = BADVAL;
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->cs.soil1c);
    ctrl->ns.soil1n = BADVAL;
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->cs.soil2c);
    ctrl->ns.soil2n = BADVAL;
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->cs.soil3c);
    ctrl->ns.soil3n = BADVAL;
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->cs.soil4c);
    ctrl->ns.soil4n = BADVAL;

    FindLine (bgc_file, "N_STATE");
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->ns.litr1n);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ctrl->ns.sminn);

    FindLine (bgc_file, "DAILY_OUTPUT");
    NextLine (bgc_file, cmdstr);
    ReadKeyword (cmdstr, "LAI", &ctrl->prtvrbl[LAI_CTRL], 'i');
    NextLine (bgc_file, cmdstr);
    ReadKeyword (cmdstr, "VEGC", &ctrl->prtvrbl[VEGC_CTRL], 'i');
    NextLine (bgc_file, cmdstr);
    ReadKeyword (cmdstr, "LITRC", &ctrl->prtvrbl[LITRC_CTRL], 'i');
    NextLine (bgc_file, cmdstr);
    ReadKeyword (cmdstr, "SOILC", &ctrl->prtvrbl[SOILC_CTRL], 'i');
    NextLine (bgc_file, cmdstr);
    ReadKeyword (cmdstr, "TOTALC", &ctrl->prtvrbl[TOTALC_CTRL], 'i');
    NextLine (bgc_file, cmdstr);
    ReadKeyword (cmdstr, "NPP", &ctrl->prtvrbl[NPP_CTRL], 'i');
    NextLine (bgc_file, cmdstr);
    ReadKeyword (cmdstr, "NEP", &ctrl->prtvrbl[NEP_CTRL], 'i');
    NextLine (bgc_file, cmdstr);
    ReadKeyword (cmdstr, "NEE", &ctrl->prtvrbl[NEE_CTRL], 'i');
    NextLine (bgc_file, cmdstr);
    ReadKeyword (cmdstr, "GPP", &ctrl->prtvrbl[GPP_CTRL], 'i');
    NextLine (bgc_file, cmdstr);
    ReadKeyword (cmdstr, "SMINN", &ctrl->prtvrbl[SMINN_CTRL], 'i');

    fclose (bgc_file);
}

//    /* Read soil moisture and soil temperature "forcing" */
//    if (ctrl->spinup == 1)
//    {
//        /* Read root zone soil water content forcing */
//        bgc->forcing.ts[SWC_TS] = (ts_struct *) malloc (sizeof (ts_struct));
//        sprintf (fn, "input/%s/%s.rootw.dat", project, project);
//        ReadBinFile (&bgc->forcing.ts[SWC_TS][0], fn, pihm->numele);
//
//        /* Read total soil water storage forcing */
//        bgc->forcing.ts[TOTALW_TS] = (ts_struct *) malloc (sizeof (ts_struct));
//        sprintf (fn, "input/%s/%s.totalw.dat", project, project);
//        ReadBinFile (&bgc->forcing.ts[TOTALW_TS][0], fn, pihm->numele + pihm->numriv);
//
//        /* Read soil temperature forcing */
//        bgc->forcing.ts[STC_TS] = (ts_struct *) malloc (sizeof (ts_struct));
//        sprintf (fn, "input/%s/%s.stc0.dat", project, project);
//        ReadBinFile (&bgc->forcing.ts[STC_TS][0], fn, pihm->numele);
//
//        /* Read subsurface flux forcing */
//        bgc->forcing.ts[SUBFLX_TS] = (ts_struct *) malloc (3 * sizeof (ts_struct));
//        for (k = 0; k < 3; k++)
//        {
//            sprintf (fn, "input/%s/%s.subflx%d.dat", project, project, k);
//            ReadBinFile (&bgc->forcing.ts[SUBFLX_TS][k], fn, pihm->numele);
//        }
//
//        /* Read surface flux forcing */
//        bgc->forcing.ts[SURFFLX_TS] = (ts_struct *) malloc (3 * sizeof (ts_struct));
//        for (k = 0; k < 3; k++)
//        {
//            sprintf (fn, "input/%s/%s.surfflx%d.dat", project, project, k);
//            ReadBinFile (&bgc->forcing.ts[SURFFLX_TS][k], fn, pihm->numele);
//        }
//
//        /* Read river flux forcing */
//        bgc->forcing.ts[RIVFLX_TS] = (ts_struct *) malloc (11 * sizeof (ts_struct));
//        for (k = 0; k < 11; k++)
//        {
//            sprintf (fn, "input/%s/%s.rivflx%d.dat", project, project, k);
//            ReadBinFile (&bgc->forcing.ts[RIVFLX_TS][k], fn, pihm->numriv);
//        }
//    }
//
//    /* Copy initial conditions to every model grid */
//    bgc->grid = (bgc_grid *) malloc (pihm->numele * sizeof (bgc_grid));
//
//    for (i = 0; i < pihm->numele; i++)
//    {
//        bgc->grid[i].ws = ws;
//        bgc->grid[i].cinit = cinit;
//        bgc->grid[i].cs = ctrl->cs;
//        bgc->grid[i].ns = ctrl->ns;
//        if (pihm->attrib_tbl.lc[i] == 4)
//        {
//            bgc->grid[i].epc = bgc->epclist.epc[EPC_DBF];
//        }
//        else if (pihm->attrib_tbl.lc[i] == 1)
//        {
//            bgc->grid[i].epc = bgc->epclist.epc[EPC_ENF];
//        }
//        else if (pihm->attrib_tbl.lc[i] == 5)
//        {
//            bgc->grid[i].epc = bgc->epclist.epc[EPC_MIXED];
//        }
//
//        bgc->grid[i].ns.cwdn = cs.cwdc / bgc->grid[i].epc.deadwood_cn;
//        bgc->grid[i].ns.litr2n = cs.litr2c / bgc->grid[i].epc.leaflitr_cn;
//        bgc->grid[i].ns.litr3n = cs.litr3c / bgc->grid[i].epc.leaflitr_cn;
//        bgc->grid[i].ns.litr4n = cs.litr4c / bgc->grid[i].epc.leaflitr_cn;
//        bgc->grid[i].ns.soil1n = cs.soil1c / SOIL1_CN;
//        bgc->grid[i].ns.soil2n = cs.soil2c / SOIL2_CN;
//        bgc->grid[i].ns.soil3n = cs.soil3c / SOIL3_CN;
//        bgc->grid[i].ns.soil4n = cs.soil4c / SOIL4_CN;
//
//        if (bgc->grid[i].epc.evergreen == 1)
//        {
//            bgc->grid[i].epv.dormant_flag = 0.0;
//        }
//        else
//        {
//            bgc->grid[i].epv.dormant_flag = 1.0;
//        }
//        //bgc->grid[i].epv.days_active = 0.;
//        bgc->grid[i].epv.onset_flag = 0.0;
//        bgc->grid[i].epv.onset_counter = 0.0;
//        bgc->grid[i].epv.onset_gddflag = 0.0;
//        bgc->grid[i].epv.onset_fdd = 0.0;
//        bgc->grid[i].epv.onset_gdd = 0.0;
//        bgc->grid[i].epv.onset_swi = 0.0;
//        bgc->grid[i].epv.offset_flag = 0.0;
//        bgc->grid[i].epv.offset_counter = 0.0;
//        bgc->grid[i].epv.offset_fdd = 0.0;
//        bgc->grid[i].epv.offset_swi = 0.0;
//        bgc->grid[i].epv.lgsf = 0.0;
//        bgc->grid[i].epv.bglfr = 0.0;
//        bgc->grid[i].epv.bgtr = 0.0;
//        bgc->grid[i].epv.annavg_t2m = 280.0;
//        bgc->grid[i].epv.tempavg_t2m = 0.0;
//    }
//
//    bgc->riv = (bgc_river *)malloc (pihm->numriv * sizeof (bgc_river));
//
//    for (i = 0; i < pihm->numriv; i++)
//    {
//        bgc->riv[i].soilw = 0.0;
//        bgc->riv[i].sminn = 0.0;
//        bgc->riv[i].nleached_snk = 0.0;
//        bgc->riv[i].sminn_leached = 0.0;
//    }
//}

void ReadEPC (epctbl_struct *epctbl)
{
    int             i;
    char            fn[MAXSTRING];
    double          t1, t2, t3, t4, r1;
    FILE           *epc_file;
    char            cmdstr[MAXSTRING];

    epctbl->woody = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->evergreen = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->c3_flag = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->phenology_flag = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->onday = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->offday = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->transfer_days = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->litfall_days = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->leaf_turnover = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->froot_turnover = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->livewood_turnover = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->daily_mortality_turnover =
        (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->daily_fire_turnover =
        (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->alloc_frootc_leafc = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->alloc_newstemc_newleafc =
        (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->alloc_newlivewoodc_newwoodc =
        (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->alloc_crootc_stemc = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->alloc_prop_curgrowth =
        (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->avg_proj_sla = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->sla_ratio = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->lai_ratio = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->ext_coef = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->flnr = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->psi_open = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->psi_close = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->vpd_open = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->vpd_close = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->froot_cn = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->leaf_cn = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->livewood_cn = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->deadwood_cn = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->leaflitr_cn = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->leaflitr_flab = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->leaflitr_fucel = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->leaflitr_fscel = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->leaflitr_flig = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->frootlitr_flab = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->frootlitr_fucel = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->frootlitr_fscel = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->frootlitr_flig = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->deadwood_fucel = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->deadwood_fscel = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->deadwood_flig = (double *)malloc (NLCTYPE * sizeof (double));

    /* Read epc files */
    if (verbose_mode)
    {
        printf ("\nRead ecophysiological constant files\n");
    }

    for (i = 0; i < NLCTYPE; i++)
    {
        switch (i + 1)
        {
            case ENF:
                strcpy (fn, "input/epc/enf.epc");
                epc_file = fopen (fn, "r");
                break;
            case EBF:
                strcpy (fn, "input/epc/ebf.epc");
                epc_file = fopen (fn, "r");
                break;
            case DNF:
                strcpy (fn, "input/epc/dnf.epc");
                epc_file = fopen (fn, "r");
                break;
            case DBF:
                strcpy (fn, "input/epc/dbf.epc");
                epc_file = fopen (fn, "r");
                break;
            case GRASS:
                strcpy (fn, "input/epc/c3grass.epc");
                epc_file = fopen (fn, "r");
                break;
            case CLOSE_SHRUB:
                strcpy (fn, "input/epc/shrub.epc");
                epc_file = fopen (fn, "r");
                break;
            case OPEN_SHRUB:
                strcpy (fn, "input/epc/shrub.epc");
                epc_file = fopen (fn, "r");
                break;
            default:
                strcpy (fn, "N/A");
        }

        if (strcasecmp (fn, "N/A") != 0)
        {
            if (NULL == epc_file)
            {
                fprintf (stderr, "Error opening %s.\n", fn);
                PIHMError (1);
            }

            if (verbose_mode)
            {
                printf ("\nReading %s...\n", fn);
            }

            /* Skip header file */
            fgets (cmdstr, MAXSTRING, epc_file);
            /* Read epc */
            /* woody/non-woody flag */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%d", &epctbl->woody[i]);
            /* evergreen/deciduous flag */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%d", &epctbl->evergreen[i]);
            /* C3/C4 flag */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%d", &epctbl->c3_flag[i]);
            /* transfer days */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%d", &epctbl->transfer_days[i]);
            /* litter fall days */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%d", &epctbl->litfall_days[i]);
            /* leaf turnover */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->leaf_turnover[i]);
            /* force leaf turnover fraction to 1.0 if deciduous */
            if (!epctbl->evergreen[i])
            {
                epctbl->leaf_turnover[i] = 1.0;
            }
            epctbl->froot_turnover[i] = epctbl->leaf_turnover[i];
            /* live wood turnover */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->livewood_turnover[i]);
            /* whole-plant mortality */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t1);
            epctbl->daily_mortality_turnover[i] = t1 / 365.0;
            /* fire mortality */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t1);
            epctbl->daily_fire_turnover[i] = t1 / 365.0;
            /* froot C:leaf C */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->alloc_frootc_leafc[i]);
            /* new stem C:new leaf C */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->alloc_newstemc_newleafc[i]);
            /* new livewood C:new wood C */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->alloc_newlivewoodc_newwoodc[i]);
            /* croot C:stem C */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->alloc_crootc_stemc[i]);
            /* new growth:storage growth */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->alloc_prop_curgrowth[i]);
            /* force storage growth to 0.0 if evergreen (following CLM-CN) */
            if (epctbl->evergreen[i])
            {
                epctbl->alloc_prop_curgrowth[i] = 1.0;
            }
            /* average leaf C:N */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->leaf_cn[i]);
            /* leaf litter C:N */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->leaflitr_cn[i]);
            /* test for leaflitter C:N > leaf C:N */
            if (epctbl->leaflitr_cn[i] < epctbl->leaf_cn[i])
            {
                printf ("Error: leaf litter C:N must be >= leaf C:N\n");
                printf
                    ("change the values in ECOPHYS block of initialization file %s\n",
                    fn);
                exit (1);
            }
            /* initial fine root C:N */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->froot_cn[i]);
            /* initial livewood C:N */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->livewood_cn[i]);
            /* initial deadwood C:N */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->deadwood_cn[i]);
            /* test for deadwood C:N > livewood C:N */
            if (epctbl->deadwood_cn[i] < epctbl->livewood_cn[i])
            {
                printf ("Error: livewood C:N must be >= deadwood C:N\n");
                printf
                    ("change the values in ECOPHYS block of initialization file\n");
                exit (1);
            }
            /* leaf litter labile proportion */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t1);
            epctbl->leaflitr_flab[i] = t1;
            /* leaf litter cellulose proportion */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t2);
            /* leaf litter lignin proportion */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t3);
            epctbl->leaflitr_flig[i] = t3;

            /* test for litter fractions sum to 1.0 */
            if (fabs (t1 + t2 + t3 - 1.0) > FLT_COND_TOL)
            {
                printf ("Error:\n");
                printf
                    ("leaf litter proportions of labile, cellulose, and lignin\n");
                printf
                    ("must sum to 1.0. Check initialization file and try again %s.\n",
                    fn);
                exit (1);
            }
            /* calculate shielded and unshielded cellulose fraction */
            r1 = t3 / t2;
            if (r1 <= 0.45)
            {
                epctbl->leaflitr_fscel[i] = 0.0;
                epctbl->leaflitr_fucel[i] = t2;
            }
            else if (r1 > 0.45 && r1 < 0.7)
            {
                t4 = (r1 - 0.45) * 3.2;
                epctbl->leaflitr_fscel[i] = t4 * t2;
                epctbl->leaflitr_fucel[i] = (1.0 - t4) * t2;
            }
            else
            {
                epctbl->leaflitr_fscel[i] = 0.8 * t2;
                epctbl->leaflitr_fucel[i] = 0.2 * t2;
            }
            /* froot litter labile proportion */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t1);
            epctbl->frootlitr_flab[i] = t1;
            /* froot litter cellulose proportion */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t2);
            /* froot litter lignin proportion */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t3);
            epctbl->frootlitr_flig[i] = t3;

            /* test for litter fractions sum to 1.0 */
            if (fabs (t1 + t2 + t3 - 1.0) > FLT_COND_TOL)
            {
                printf ("Error:\n");
                printf
                    ("froot litter proportions of labile, cellulose, and lignin\n");
                printf
                    ("must sum to 1.0. Check initialization file and try again.\n");
                exit (1);
            }
            /* calculate shielded and unshielded cellulose fraction */
            r1 = t3 / t2;
            if (r1 <= 0.45)
            {
                epctbl->frootlitr_fscel[i] = 0.0;
                epctbl->frootlitr_fucel[i] = t2;
            }
            else if (r1 > 0.45 && r1 < 0.7)
            {
                t4 = (r1 - 0.45) * 3.2;
                epctbl->frootlitr_fscel[i] = t4 * t2;
                epctbl->frootlitr_fucel[i] = (1.0 - t4) * t2;
            }
            else
            {
                epctbl->frootlitr_fscel[i] = 0.8 * t2;
                epctbl->frootlitr_fucel[i] = 0.2 * t2;
            }
            /* dead wood cellulose */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t1);
            /* dead wood lignin */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t2);
            epctbl->deadwood_flig[i] = t2;
            /* test for litter fractions sum to 1.0 */
            if (fabs (t1 + t2 - 1.0) > FLT_COND_TOL)
            {
                printf ("Error:\n");
                printf
                    ("deadwood proportions of cellulose and lignin must sum\n");
                printf ("to 1.0. Check initialization file and try again.\n");
                exit (1);
            }
            /* calculate shielded and unshielded cellulose fraction */
            r1 = t2 / t1;
            if (r1 <= 0.45)
            {
                epctbl->deadwood_fscel[i] = 0.0;
                epctbl->deadwood_fucel[i] = t1;
            }
            else if (r1 > 0.45 && r1 < 0.7)
            {
                t4 = (r1 - 0.45) * 3.2;
                epctbl->deadwood_fscel[i] = t4 * t1;
                epctbl->deadwood_fucel[i] = (1.0 - t4) * t1;
            }
            else
            {
                epctbl->deadwood_fscel[i] = 0.8 * t1;
                epctbl->deadwood_fucel[i] = 0.2 * t1;
            }
            /* canopy light ext coef */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->ext_coef[i]);
            /* all to projected LA ratio */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->lai_ratio[i]);
            /* canopy average projected specific leaf area */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->avg_proj_sla[i]);
            /* sunlit SLA ratio */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->sla_ratio[i]);
            /* Rubisco N fraction */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->flnr[i]);
            /* psi_sat */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->psi_open[i]);
            /* psi_close */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->psi_close[i]);
            /* vpd_max */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->vpd_open[i]);
            /* vpd_min */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->vpd_close[i]);

#ifdef _DEBUG_
            printf ("WOODY%d\t", epctbl->woody[i]);
            printf ("EVERGREEN%d\t", epctbl->evergreen[i]);
            printf ("C3%d\t", epctbl->c3_flag[i]);
            printf ("TRANSFERDAY%lf\t", epctbl->transfer_days[i]);
            printf ("LITFALLDAY%lf\t", epctbl->litfall_days[i]);
            printf ("LEAFTURNOVER%lf\t", epctbl->leaf_turnover[i]);
            printf ("FROOTTURNOVER%lf\t", epctbl->froot_turnover[i]);
            printf ("LIVEWOODTURNOVER%lf\t", epctbl->livewood_turnover[i]);
            printf ("DAILYMORTALITYTURNOVER%lf\t",
                epctbl->daily_mortality_turnover[i]);
            printf ("DAILYFIRETURNOVER%lf\t", epctbl->daily_fire_turnover[i]);
            printf ("FROOTC/LEAFC%lf\t", epctbl->alloc_frootc_leafc[i]);
            printf ("NEWSTEMC/NEWLEAFC%lf\t",
                epctbl->alloc_newstemc_newleafc[i]);
            printf ("NEWLIVEWOODC/NEWWOODC%lf\t",
                epctbl->alloc_newlivewoodc_newwoodc[i]);
            printf ("CROOTC/STEMC%lf\t", epctbl->alloc_crootc_stemc[i]);
            printf ("PROP%lf\t", epctbl->alloc_prop_curgrowth[i]);
            printf ("LEAFCN%lf\t", epctbl->leaf_cn[i]);
            printf ("LEAFLITRCN%lf\t", epctbl->leaflitr_cn[i]);
            printf ("FROOTCN%lf\t", epctbl->froot_cn[i]);
            printf ("LIVEWOODCN%lf\t", epctbl->livewood_cn[i]);
            printf ("DEADWOODCN%lf\t", epctbl->deadwood_cn[i]);
            printf ("LEAFLITR %lf %lf %lf\t", epctbl->leaflitr_fscel[i],
                epctbl->leaflitr_fucel[i], epctbl->leaflitr_flig[i]);
            printf ("FROOTLITR %lf %lf %lf\t", epctbl->frootlitr_fscel[i],
                epctbl->frootlitr_fucel[i], epctbl->frootlitr_flig[i]);
            printf ("DEADWOOD %lf %lf %lf\t", epctbl->deadwood_fscel[i],
                epctbl->deadwood_fucel[i], epctbl->deadwood_flig[i]);
            printf ("EXTCOEF %lf\t", epctbl->ext_coef[i]);
            printf ("LAIRATIO %lf\t", epctbl->lai_ratio[i]);
            printf ("PROJSLA %lf\t", epctbl->avg_proj_sla[i]);
            printf ("SLARATIO %lf\t", epctbl->sla_ratio[i]);
            printf ("FLNR %lf\t", epctbl->flnr[i]);
            printf ("GLMAX %lf\t", epctbl->gl_smax[i]);
            printf ("GLC %lf\t", epctbl->gl_c[i]);
            printf ("GL_BL %lf\t", epctbl->gl_bl[i]);
            printf ("PSIOPEN %lf\t", epctbl->psi_open[i]);
            printf ("PSICLOSE %lf\t", epctbl->psi_close[i]);
            printf ("VPDOPEN %lf\t", epctbl->vpd_open[i]);
            printf ("VPDCLOSE%lf\n", epctbl->vpd_close[i]);
#endif
        }
        else
        {
            epctbl->woody[i] = BADVAL;
            epctbl->evergreen[i] = BADVAL;
            epctbl->c3_flag[i] = BADVAL;
            epctbl->phenology_flag[i] = BADVAL;
            epctbl->onday[i] = BADVAL;
            epctbl->offday[i] = BADVAL;
            epctbl->transfer_days[i] = BADVAL;
            epctbl->litfall_days[i] = BADVAL;
            epctbl->leaf_turnover[i] = BADVAL;
            epctbl->froot_turnover[i] = BADVAL;
            epctbl->livewood_turnover[i] = BADVAL;
            epctbl->daily_mortality_turnover[i] = BADVAL;
            epctbl->daily_fire_turnover[i] = BADVAL;
            epctbl->alloc_frootc_leafc[i] = BADVAL;
            epctbl->alloc_newstemc_newleafc[i] = BADVAL;
            epctbl->alloc_newlivewoodc_newwoodc[i] = BADVAL;
            epctbl->alloc_crootc_stemc[i] = BADVAL;
            epctbl->alloc_prop_curgrowth[i] = BADVAL;
            epctbl->avg_proj_sla[i] = BADVAL;
            epctbl->sla_ratio[i] = BADVAL;
            epctbl->lai_ratio[i] = BADVAL;
            epctbl->ext_coef[i] = BADVAL;
            epctbl->flnr[i] = BADVAL;
            epctbl->psi_open[i] = BADVAL;
            epctbl->psi_close[i] = BADVAL;
            epctbl->vpd_open[i] = BADVAL;
            epctbl->vpd_close[i] = BADVAL;
            epctbl->froot_cn[i] = BADVAL;
            epctbl->leaf_cn[i] = BADVAL;
            epctbl->livewood_cn[i] = BADVAL;
            epctbl->deadwood_cn[i] = BADVAL;
            epctbl->leaflitr_cn[i] = BADVAL;
            epctbl->leaflitr_flab[i] = BADVAL;
            epctbl->leaflitr_fucel[i] = BADVAL;
            epctbl->leaflitr_fscel[i] = BADVAL;
            epctbl->leaflitr_flig[i] = BADVAL;
            epctbl->frootlitr_flab[i] = BADVAL;
            epctbl->frootlitr_fucel[i] = BADVAL;
            epctbl->frootlitr_fscel[i] = BADVAL;
            epctbl->frootlitr_flig[i] = BADVAL;
            epctbl->deadwood_fucel[i] = BADVAL;
            epctbl->deadwood_fscel[i] = BADVAL;
            epctbl->deadwood_flig[i] = BADVAL;
        }
    }
}

void ReadAnnFile (tsdata_struct *ts, char *fn)
{
    FILE           *fid;
    time_t          rawtime;
    struct tm      *timeinfo;
    char            cmdstr[MAXSTRING];
    int             i;
    int             match;

    timeinfo = (struct tm *)malloc (sizeof (struct tm));

    fid = fopen (fn, "r");

    if (NULL == fid)
    {
        fprintf (stderr, "Error opening %s.\n", fn);
        PIHMError (1);
    }

    if (verbose_mode)
    {
        printf ("Reading %s...\n", fn);
    }

    ts->length = CountLine (fid, cmdstr, 1, "EOF");
    ts->ftime = (int *)malloc (ts->length * sizeof (int));
    ts->data = (double **)malloc (ts->length * sizeof (double *));

    FindLine (fid, "BOF");
    for (i = 0; i < ts->length; i++)
    {
        ts->data[i] = (double *)malloc (sizeof (double));
        NextLine (fid, cmdstr);
        match =
            sscanf (cmdstr, "%d %lf", &timeinfo->tm_year, &ts->data[i][0]);

        if (match != 2)
        {
            fprintf (stderr, "Error reading %s.\n", fn);
            fprintf (stderr, "Please check file format.\n");
            PIHMError (1);
        }

        timeinfo->tm_year = timeinfo->tm_year - 1900;
        timeinfo->tm_mon = 0;
        timeinfo->tm_mday = 1;
        timeinfo->tm_hour = 0;
        timeinfo->tm_min = 0;
        timeinfo->tm_sec = 0;
        rawtime = timegm (timeinfo);
        ts->ftime[i] = (int)rawtime;
    }

    fclose (fid);
}
