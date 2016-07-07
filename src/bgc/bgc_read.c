#include "pihm.h"

void ReadBGC (char *fn, ctrl_struct *ctrl, co2control_struct *co2, ndepcontrol_struct *ndepctrl, char *co2_fn, char *ndep_fn)
{
    FILE           *bgc_file;
    struct tm      *timestamp;
    time_t          rawtime;
    char            cmdstr[MAXSTRING];

    timestamp = (struct tm *)malloc (sizeof (struct tm));

    /* Read bgc simulation control file */
    bgc_file = fopen (fn, "r");
    CheckFile (bgc_file, fn);

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
    sscanf (cmdstr, "%d", &ctrl->spinup);
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
    ReadKeywordInt (cmdstr, "LAI", &ctrl->prtvrbl[LAI_CTRL]);
    NextLine (bgc_file, cmdstr);
    ReadKeywordInt (cmdstr, "VEGC", &ctrl->prtvrbl[VEGC_CTRL]);
    NextLine (bgc_file, cmdstr);
    ReadKeywordInt (cmdstr, "LITRC", &ctrl->prtvrbl[LITRC_CTRL]);
    NextLine (bgc_file, cmdstr);
    ReadKeywordInt (cmdstr, "SOILC", &ctrl->prtvrbl[SOILC_CTRL]);
    NextLine (bgc_file, cmdstr);
    ReadKeywordInt (cmdstr, "TOTALC", &ctrl->prtvrbl[TOTALC_CTRL]);
    NextLine (bgc_file, cmdstr);
    ReadKeywordInt (cmdstr, "NPP", &ctrl->prtvrbl[NPP_CTRL]);
    NextLine (bgc_file, cmdstr);
    ReadKeywordInt (cmdstr, "NEP", &ctrl->prtvrbl[NEP_CTRL]);
    NextLine (bgc_file, cmdstr);
    ReadKeywordInt (cmdstr, "NEE", &ctrl->prtvrbl[NEE_CTRL]);
    NextLine (bgc_file, cmdstr);
    ReadKeywordInt (cmdstr, "GPP", &ctrl->prtvrbl[GPP_CTRL]);
    NextLine (bgc_file, cmdstr);
    ReadKeywordInt (cmdstr, "SMINN", &ctrl->prtvrbl[SMINN_CTRL]);

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

void ReadEPC (epclist_struct *epclist)
{
    int             i;
    char            fn[MAXSTRING];
    double          t1, t2, t3, t4, r1;
    FILE           *epc_file;
    char            cmdstr[MAXSTRING];

    epconst_struct *epc;

    /* Read epc files */
    epclist->nvegtypes = NVEGTYPES;
    epclist->epc = (epconst_struct *)
        malloc (epclist->nvegtypes * sizeof (epconst_struct));

    if (verbose_mode)
    {
        printf ("\nRead ecophysiological constant files\n");
    }

    for (i = 0; i < epclist->nvegtypes; i++)
    {
        switch (i)
        {
            case EPC_C3GRASS:
                strcpy (fn, "input/epc/c3grass.epc");
                epc_file = fopen (fn, "r");
                break;
            case EPC_C4GRASS:
                strcpy (fn, "input/epc/c4grass.epc");
                epc_file = fopen (fn, "r");
                break;
            case EPC_DBF:
                strcpy (fn, "input/epc/dbf.epc");
                epc_file = fopen (fn, "r");
                break;
            case EPC_DNF:
                strcpy (fn, "input/epc/dnf.epc");
                epc_file = fopen (fn, "r");
                break;
            case EPC_EBF:
                strcpy (fn, "input/epc/ebf.epc");
                epc_file = fopen (fn, "r");
                break;
            case EPC_ENF:
                strcpy (fn, "input/epc/enf.epc");
                epc_file = fopen (fn, "r");
                break;
            case EPC_SHRUB:
                strcpy (fn, "input/epc/shrub.epc");
                epc_file = fopen (fn, "r");
                break;
            //case EPC_MIXED:
            //    strcpy (fn, "input/epc/mix.epc");
            //    epc_file = fopen (fn, "r");
            //    break;
        }

        CheckFile (epc_file, fn);

        epc = &(epclist->epc[i]);

        /* Skip header file */
        fgets (cmdstr, MAXSTRING, epc_file);
        /* Read epc */
        /* woody/non-woody flag */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%d", &epc->woody);
        /* evergreen/deciduous flag */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%d", &epc->evergreen);
        /* C3/C4 flag */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%d", &epc->c3_flag);
        /* transfer days */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->transfer_days);
        /* litter fall days */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->litfall_days);
        /* leaf turnover */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->leaf_turnover);
        /* force leaf turnover fraction to 1.0 if deciduous */
        if (!epc->evergreen)
        {
            epc->leaf_turnover = 1.0;
        }
        epc->froot_turnover = epc->leaf_turnover;
        /* live wood turnover */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->livewood_turnover);
        /* whole-plant mortality */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &t1);
        epc->daily_mortality_turnover = t1 / 365.0;
        /* fire mortality */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &t1);
        epc->daily_fire_turnover = t1 / 365.0;
        /* froot C:leaf C */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->alloc_frootc_leafc);
        /* new stem C:new leaf C */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->alloc_newstemc_newleafc);
        /* new livewood C:new wood C */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->alloc_newlivewoodc_newwoodc);
        /* croot C:stem C */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->alloc_crootc_stemc);
        /* new growth:storage growth */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->alloc_prop_curgrowth);
        /* force storage growth to 0.0 if evergreen (following CLM-CN) */
        if (epc->evergreen)
        {
            epc->alloc_prop_curgrowth = 1.0;
        }
        /* average leaf C:N */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->leaf_cn);
        /* leaf litter C:N */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->leaflitr_cn);
        /* test for leaflitter C:N > leaf C:N */
        if (epc->leaflitr_cn < epc->leaf_cn)
        {
            printf ("Error: leaf litter C:N must be >= leaf C:N\n");
            printf ("change the values in ECOPHYS block of initialization file %s\n", fn);
            exit (1);
        }
        /* initial fine root C:N */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->froot_cn);
        /* initial livewood C:N */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->livewood_cn);
        /* initial deadwood C:N */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->deadwood_cn);
        /* test for deadwood C:N > livewood C:N */
        if (epc->deadwood_cn < epc->livewood_cn)
        {
            printf ("Error: livewood C:N must be >= deadwood C:N\n");
            printf ("change the values in ECOPHYS block of initialization file\n");
            exit (1);
        }
        /* leaf litter labile proportion */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &t1);
        epc->leaflitr_flab = t1;
        /* leaf litter cellulose proportion */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &t2);
        /* leaf litter lignin proportion */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &t3);
        epc->leaflitr_flig = t3;

        /* test for litter fractions sum to 1.0 */
        if (fabs (t1 + t2 + t3 - 1.0) > FLT_COND_TOL)
        {
            printf ("Error:\n");
            printf ("leaf litter proportions of labile, cellulose, and lignin\n");
            printf ("must sum to 1.0. Check initialization file and try again %s.\n", fn);
            exit (1);
        }
        /* calculate shielded and unshielded cellulose fraction */
        r1 = t3 / t2;
        if (r1 <= 0.45)
        {
            epc->leaflitr_fscel = 0.0;
            epc->leaflitr_fucel = t2;
        }
        else if (r1 > 0.45 && r1 < 0.7)
        {
            t4 = (r1 - 0.45) * 3.2;
            epc->leaflitr_fscel = t4 * t2;
            epc->leaflitr_fucel = (1.0 - t4) * t2;
        }
        else
        {
            epc->leaflitr_fscel = 0.8 * t2;
            epc->leaflitr_fucel = 0.2 * t2;
        }
        /* froot litter labile proportion */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &t1);
        epc->frootlitr_flab = t1;
        /* froot litter cellulose proportion */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &t2);
        /* froot litter lignin proportion */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &t3);
        epc->frootlitr_flig = t3;

        /* test for litter fractions sum to 1.0 */
        if (fabs (t1 + t2 + t3 - 1.0) > FLT_COND_TOL)
        {
            printf ("Error:\n");
            printf ("froot litter proportions of labile, cellulose, and lignin\n");
            printf ("must sum to 1.0. Check initialization file and try again.\n");
            exit (1);
        }
        /* calculate shielded and unshielded cellulose fraction */
        r1 = t3 / t2;
        if (r1 <= 0.45)
        {
            epc->frootlitr_fscel = 0.0;
            epc->frootlitr_fucel = t2;
        }
        else if (r1 > 0.45 && r1 < 0.7)
        {
            t4 = (r1 - 0.45) * 3.2;
            epc->frootlitr_fscel = t4 * t2;
            epc->frootlitr_fucel = (1.0 - t4) * t2;
        }
        else
        {
            epc->frootlitr_fscel = 0.8 * t2;
            epc->frootlitr_fucel = 0.2 * t2;
        }
        /* dead wood cellulose */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &t1);
        /* dead wood lignin */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &t2);
        epc->deadwood_flig = t2;
        /* test for litter fractions sum to 1.0 */
        if (fabs (t1 + t2 - 1.0) > FLT_COND_TOL)
        {
            printf ("Error:\n");
            printf ("deadwood proportions of cellulose and lignin must sum\n");
            printf ("to 1.0. Check initialization file and try again.\n");
            exit (1);
        }
        /* calculate shielded and unshielded cellulose fraction */
        r1 = t2 / t1;
        if (r1 <= 0.45)
        {
            epc->deadwood_fscel = 0.0;
            epc->deadwood_fucel = t1;
        }
        else if (r1 > 0.45 && r1 < 0.7)
        {
            t4 = (r1 - 0.45) * 3.2;
            epc->deadwood_fscel = t4 * t1;
            epc->deadwood_fucel = (1.0 - t4) * t1;
        }
        else
        {
            epc->deadwood_fscel = 0.8 * t1;
            epc->deadwood_fucel = 0.2 * t1;
        }
        /* canopy light ext coef */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->ext_coef);
        /* all to projected LA ratio */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->lai_ratio);
        /* canopy average projected specific leaf area */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->avg_proj_sla);
        /* sunlit SLA ratio */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->sla_ratio);
        /* Rubisco N fraction */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->flnr);
        /* psi_sat */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->psi_open);
        /* psi_close */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->psi_close);
        /* vpd_max */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->vpd_open);
        /* vpd_min */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->vpd_close);

#ifdef _DEBUG_
        printf ("WOODY%d\t", epc->woody);
        printf ("EVERGREEN%d\t", epc->evergreen);
        printf ("C3%d\t", epc->c3_flag);
        printf ("TRANSFERDAY%lf\t", epc->transfer_days);
        printf ("LITFALLDAY%lf\t", epc->litfall_days);
        printf ("LEAFTURNOVER%lf\t", epc->leaf_turnover);
        printf ("FROOTTURNOVER%lf\t", epc->froot_turnover);
        printf ("LIVEWOODTURNOVER%lf\t", epc->livewood_turnover);
        printf ("DAILYMORTALITYTURNOVER%lf\t", epc->daily_mortality_turnover);
        printf ("DAILYFIRETURNOVER%lf\t", epc->daily_fire_turnover);
        printf ("FROOTC/LEAFC%lf\t", epc->alloc_frootc_leafc);
        printf ("NEWSTEMC/NEWLEAFC%lf\t", epc->alloc_newstemc_newleafc);
        printf ("NEWLIVEWOODC/NEWWOODC%lf\t", epc->alloc_newlivewoodc_newwoodc);
        printf ("CROOTC/STEMC%lf\t", epc->alloc_crootc_stemc);
        printf ("PROP%lf\t", epc->alloc_prop_curgrowth);
        printf ("LEAFCN%lf\t", epc->leaf_cn);
        printf ("LEAFLITRCN%lf\t", epc->leaflitr_cn);
        printf ("FROOTCN%lf\t", epc->froot_cn);
        printf ("LIVEWOODCN%lf\t", epc->livewood_cn);
        printf ("DEADWOODCN%lf\t", epc->deadwood_cn);
        printf ("LEAFLITR %lf %lf %lf\t", epc->leaflitr_fscel, epc->leaflitr_fucel, epc->leaflitr_flig);
        printf ("FROOTLITR %lf %lf %lf\t", epc->frootlitr_fscel, epc->frootlitr_fucel, epc->frootlitr_flig);
        printf ("DEADWOOD %lf %lf %lf\t", epc->deadwood_fscel, epc->deadwood_fucel, epc->deadwood_flig);
        printf ("EXTCOEF %lf\t", epc->ext_coef);
        printf ("LAIRATIO %lf\t", epc->lai_ratio);
        printf ("PROJSLA %lf\t", epc->avg_proj_sla);
        printf ("SLARATIO %lf\t", epc->sla_ratio);
        printf ("FLNR %lf\t", epc->flnr);
        printf ("GLMAX %lf\t", epc->gl_smax);
        printf ("GLC %lf\t", epc->gl_c);
        printf ("GL_BL %lf\t", epc->gl_bl);
        printf ("PSIOPEN %lf\t", epc->psi_open);
        printf ("PSICLOSE %lf\t", epc->psi_close);
        printf ("VPDOPEN %lf\t", epc->vpd_open);
        printf ("VPDCLOSE%lf\n", epc->vpd_close);
#endif
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
    CheckFile (fid, fn);

    ts->length = CountLine (fid, cmdstr, 1, "EOF");
    ts->ftime = (int *) malloc (ts->length * sizeof (int));
    ts->data = (double **) malloc (ts->length * sizeof (double *));

    FindLine (fid, "BOF");
    for (i = 0; i < ts->length; i++)
    {
        ts->data[i] = (double *) malloc (sizeof (double));
        NextLine (fid, cmdstr);
        match = sscanf (cmdstr, "%d %lf", &timeinfo->tm_year, &ts->data[i][0]);

        if (match != 2)
        {
            printf ("Cannot read annual time series!\n");
            printf ("%s file format error!\n", fn);
            PihmExit (1);
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
