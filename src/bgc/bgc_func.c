#include "bgc.h"

void BgcRead (char *simulation, bgc_struct bgc, pihm_struct pihm)
{
    int             i, k;
    double          t1, t2, t3, t4, r1;
    char            fn[MAXSTRING];
    char            project[MAXSTRING];
    char           *token;
    char            tempname[MAXSTRING];
    FILE           *epc_file;
    FILE           *bgc_file;
    char            co2_fn[MAXSTRING];
    char            ndep_fn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    struct tm      *timeinfo;
    enum epc_vegtype veg_type;
    epconst_struct *epc;
    control_struct *ctrl;
    co2control_struct *co2;
    ndepcontrol_struct *ndepctrl;

    /* Templates for model initial conditions */
    wstate_struct   ws;
    cstate_struct   cs;
    nstate_struct   ns;
    cinit_struct    cinit;

    timeinfo = (struct tm *)malloc (sizeof (struct tm));

    /* Detect if model is running in ensemble mode */
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

    /* Initialize state variables with zeroes */
    presim_state_init (&ws, &cs, &ns, &cinit);

    /* Read epc files */
    bgc->epclist.nvegtypes = NVEGTYPES;
    bgc->epclist.epc = (epconst_struct *) malloc (bgc->epclist.nvegtypes * sizeof (epconst_struct));

    if (verbose_mode)
    {
        printf ("\nRead ecophysiological constant files\n");
    }

    for (veg_type = 0; veg_type < bgc->epclist.nvegtypes; veg_type++)
    {
        switch (veg_type)
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
            case EPC_MIXED:
                strcpy (fn, "input/epc/mix.epc");
                epc_file = fopen (fn, "r");
                break;
        }
        if (epc_file == NULL)
        {
            printf ("\nFatal Error: epc file %s are in use or do not exist!\n", fn);
            exit (1);
        }

        epc = &(bgc->epclist.epc[veg_type]);

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
            epc->leaf_turnover = 1.0;
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
        epc->daily_fire_turnover = t1 / 365;
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
            epc->alloc_prop_curgrowth = 1.0;
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
        /* canopy water int coef */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->int_coef);
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
        /* gl_smax */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->gl_smax);
        /* gl_c */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->gl_c);
        /* gl_bl */
        fgets (cmdstr, MAXSTRING, epc_file);
        sscanf (cmdstr, "%lf", &epc->gl_bl);
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
        printf ("INTCOEF %lf\t", epc->int_coef);
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

    /* Read bgc simulation control file */

    ctrl = &bgc->ctrl;
    co2 = &bgc->co2;
    ndepctrl = &bgc->ndepctrl;

    sprintf (fn, "input/%s/%s.bgc", project, project);
    bgc_file = fopen (fn, "r");
    CheckFile (bgc_file, fn);

    FindLine (bgc_file, "TIME_DEFINE");
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%d", &ctrl->spinupstart);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%d", &ctrl->spinupend);
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
    sscanf (cmdstr, "%lf", &cinit.max_leafc);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &cinit.max_stemc);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &cs.cwdc);
    ns.cwdn = BADVAL;
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &cs.litr1c);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &cs.litr2c);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &cs.litr3c);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &cs.litr4c);
    ns.litr2n = BADVAL;
    ns.litr3n = BADVAL;
    ns.litr4n = BADVAL;
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &cs.soil1c);
    ns.soil1n = BADVAL;
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &cs.soil2c);
    ns.soil2n = BADVAL;
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &cs.soil3c);
    ns.soil3n = BADVAL;
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &cs.soil4c);
    ns.soil4n = BADVAL;

    FindLine (bgc_file, "N_STATE");
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ns.litr1n);
    NextLine (bgc_file, cmdstr);
    sscanf (cmdstr, "%lf", &ns.sminn);

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

    fclose (bgc_file);

    /* Read CO2 and Ndep files */
    if (co2->varco2 == 1)
    {
        bgc->forcing.ts[CO2_TS] = (ts_struct *) malloc (sizeof (ts_struct));
        ReadAnnFile (&bgc->forcing.ts[CO2_TS][0], co2_fn);
    }
    if (ndepctrl->varndep == 1)
    {
        bgc->forcing.ts[NDEP_TS] = (ts_struct *) malloc (sizeof (ts_struct));
        ReadAnnFile (&bgc->forcing.ts[NDEP_TS][0], ndep_fn);
    }

    /* Read soil moisture and soil temperature "forcing" */
    if (ctrl->spinup == 1)
    {
        /* Read root zone soil water content forcing */
        bgc->forcing.ts[SWC_TS] = (ts_struct *) malloc (sizeof (ts_struct));
        sprintf (fn, "input/%s/%s.rootw", project, project);
        ReadBinFile (&bgc->forcing.ts[SWC_TS][0], fn, pihm->numele);

        /* Read total soil water storage forcing */
        bgc->forcing.ts[SOILM_TS] = (ts_struct *) malloc (sizeof (ts_struct));
        sprintf (fn, "input/%s/%s.soilm", project, project);
        ReadBinFile (&bgc->forcing.ts[SOILM_TS][0], fn, pihm->numele);

        /* Read soil temperature forcing */
        bgc->forcing.ts[STC_TS] = (ts_struct *) malloc (sizeof (ts_struct));
        sprintf (fn, "input/%s/%s.stc", project, project);
        ReadBinFile (&bgc->forcing.ts[STC_TS][0], fn, pihm->numele);

        /* Read subsurface flux forcing */
        bgc->forcing.ts[SUBFLX_TS] = (ts_struct *) malloc (3 * sizeof (ts_struct));
        for (k = 0; k < 3; k++)
        {
            sprintf (fn, "input/%s/%s.subflx%d", project, project, k);
            ReadBinFile (&bgc->forcing.ts[SUBFLX_TS][k], fn, pihm->numele);
        }
    }

    /* Copy initial conditions to every model grid */
    bgc->grid = (bgc_grid *) malloc (pihm->numele * sizeof (bgc_grid));

    for (i = 0; i < pihm->numele; i++)
    {
        bgc->grid[i].ws = ws;
        bgc->grid[i].cinit = cinit;
        bgc->grid[i].cs = cs;
        bgc->grid[i].ns = ns;
        if (pihm->attrib_tbl.lc[i] == 4)
        {
            bgc->grid[i].epc = bgc->epclist.epc[EPC_DBF];
        }
        else if (pihm->attrib_tbl.lc[i] == 1)
        {
            bgc->grid[i].epc = bgc->epclist.epc[EPC_ENF];
        }
        else if (pihm->attrib_tbl.lc[i] == 5)
        {
            bgc->grid[i].epc = bgc->epclist.epc[EPC_MIXED];
        }

        bgc->grid[i].ns.cwdn = cs.cwdc / bgc->grid[i].epc.deadwood_cn;
        bgc->grid[i].ns.litr2n = cs.litr2c / bgc->grid[i].epc.leaflitr_cn;
        bgc->grid[i].ns.litr3n = cs.litr3c / bgc->grid[i].epc.leaflitr_cn;
        bgc->grid[i].ns.litr4n = cs.litr4c / bgc->grid[i].epc.leaflitr_cn;
        bgc->grid[i].ns.soil1n = cs.soil1c / SOIL1_CN;
        bgc->grid[i].ns.soil2n = cs.soil2c / SOIL2_CN;
        bgc->grid[i].ns.soil3n = cs.soil3c / SOIL3_CN;
        bgc->grid[i].ns.soil4n = cs.soil4c / SOIL4_CN;

        if (bgc->grid[i].epc.evergreen == 1)
        {
            bgc->grid[i].epv.dormant_flag = 0.0;
        }
        else
        {
            bgc->grid[i].epv.dormant_flag = 1.0;
        }
        //bgc->grid[i].epv.days_active = 0.;
        bgc->grid[i].epv.onset_flag = 0.0;
        bgc->grid[i].epv.onset_counter = 0.0;
        bgc->grid[i].epv.onset_gddflag = 0.0;
        bgc->grid[i].epv.onset_fdd = 0.0;
        bgc->grid[i].epv.onset_gdd = 0.0;
        bgc->grid[i].epv.onset_swi = 0.0;
        bgc->grid[i].epv.offset_flag = 0.0;
        bgc->grid[i].epv.offset_counter = 0.0;
        bgc->grid[i].epv.offset_fdd = 0.0;
        bgc->grid[i].epv.offset_swi = 0.0;
        bgc->grid[i].epv.lgsf = 0.0;
        bgc->grid[i].epv.bglfr = 0.0;
        bgc->grid[i].epv.bgtr = 0.0;
        bgc->grid[i].epv.annavg_t2m = 280.0;
        bgc->grid[i].epv.tempavg_t2m = 0.0;
    }
}

void BgcInit (char *simulation, pihm_struct pihm, lsm_struct noah, bgc_struct bgc)
{
    char            fn[MAXSTRING];
    char            project[MAXSTRING];
    FILE           *init_file;
    int             i, j;
    char           *token;
    char            tempname[MAXSTRING];

    /* Detect if model is running in ensemble mode */
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
        printf ("BGC: Initializing BGC structures\n");
    }

    for (i = 0; i < pihm->numele; i++)
    {
        bgc->grid[i].sitec.soil_alpha = noah->grid[i].vgalpha;
        bgc->grid[i].sitec.soil_beta = noah->grid[i].vgbeta;
        bgc->grid[i].sitec.vwc_sat = noah->grid[i].smcmax;
        bgc->grid[i].sitec.vwc_min = noah->grid[i].smcmin;
        bgc->grid[i].sitec.vwc_fc = noah->grid[i].smcref;
        bgc->grid[i].sitec.lat = noah->latitude;
        bgc->grid[i].sitec.lon = noah->longitude;
        bgc->grid[i].sitec.sw_alb = 0.5 * (noah->grid[i].albedomin + noah->grid[i].albedomax);
        bgc->grid[i].sitec.area = pihm->elem[i].topo.area;
        for (j = 0; j < 3; j++)
        {
            bgc->grid[i].sitec.nabr[j] = pihm->elem[i].nabr[j];
        }

        bgc->grid[i].epv.annavg_t2m = noah->genprmt.tbot_data - 273.15;

        bgc->grid[i].epc.topt = noah->grid[i].topt;
        bgc->grid[i].epc.rgl = noah->grid[i].rgl;
        bgc->grid[i].epc.hs = noah->grid[i].hs;
        bgc->grid[i].epc.gl_smax = 1.0 / noah->grid[i].rsmin;
        bgc->grid[i].epc.gl_c = 1.0 / noah->grid[i].rsmax;
        bgc->grid[i].epc.smcref = noah->grid[i].smcref;
        bgc->grid[i].epc.smcwlt = noah->grid[i].smcwlt;
    }

    /* Read initial conditions */
    if (bgc->ctrl.spinup == 0)
    {
        sprintf (fn, "input/%s/%s.bgcinit", project, simulation);
        init_file = fopen (fn, "rb");
        CheckFile (init_file, fn);

        for (i = 0; i < pihm->numele; i++)
        {
            fread (&(bgc->grid[i].restart_input), sizeof (restart_data_struct), 1, init_file);
            restart_input (&bgc->grid[i].ws, &bgc->grid[i].cs, &bgc->grid[i].ns, &bgc->grid[i].epv, &bgc->grid[i].restart_input);

            noah->grid[i].xlai = bgc->grid[i].cs.leafc * bgc->grid[i].epc.avg_proj_sla;
            noah->grid[i].cmcmax = noah->grid[i].cmcfactr * noah->grid[i].xlai;
        }
        fclose (init_file);
    }
    else
    {
        for (i = 0; i < pihm->numele; i++)
        {
            firstday (&bgc->grid[i].epc, &bgc->grid[i].cinit, &bgc->grid[i].epv, &bgc->grid[i].cs, &bgc->grid[i].ns);
        }
    }

    for (i = 0; i < pihm->numele; i++)
    {
        zero_srcsnk (&bgc->grid[i].cs, &bgc->grid[i].ns, &bgc->grid[i].ws, &bgc->grid[i].summary);
    }
}

void BgcCoupling (int t, int start_time, pihm_struct pihm, lsm_struct noah, bgc_struct bgc)
{
    static int      counter[5000];
    static int      daylight_counter[5000];
    double          dayl, prev_dayl;
    spa_data        spa;
    int             spa_result;
    time_t          rawtime;
    struct tm      *timestamp;
    double          sfctmp;
    double          solar;
    int             i, j, k;
    double          dummy[pihm->numele];
    metvar_struct  *metv;
    static int      first_balance;

    if (t == start_time)
    {
        for (i = 0; i < pihm->numele; i++)
        {
            counter[i] = 0;
            daylight_counter[i] = 0;

            metv = &(bgc->grid[i].metv);
            metv->tmax = -999.0;
            metv->tmin = 999.0;
            metv->tavg = 0.0;
            metv->tsoil = 0.0;
            metv->swc = 0.0;
            metv->soilw = 0.0;
            for (k = 0; k < 3; k++)
            {
                metv->subflux[k] = 0.0;
            }
            metv->tday = 0.0;
            metv->q2d = 0.0;
            metv->pa = 0.0;
            metv->swavgfd = 0.0;
            metv->par = 0.0;
            metv->tnight = 0.0;
        }
        first_balance = 1;
    }

    for (i = 0; i < pihm->numele; i++)
    {
        metv = &(bgc->grid[i].metv);

        sfctmp = noah->grid[i].sfctmp - 273.15;
        metv->tmax = (metv->tmax > sfctmp) ? metv->tmax : sfctmp;
        metv->tmin = (metv->tmin < sfctmp) ? metv->tmin : sfctmp;
        metv->tavg += sfctmp;
        solar = noah->grid[i].soldn;
        metv->tsoil += noah->grid[i].stc[0] - 273.15;
        metv->swc += noah->grid[i].soilw;
        metv->soilw += noah->grid[i].soilm;
        for (k = 0; k < 3; k++)
        {
            //if (pihm->elem[i].bc_type[k] < 0.0 && noah->grid[i].avgsubflux[k] < 0.0)
            //{
            //    metv->subflux[k] = 0.0;
            //}
            //else
            //{
                metv->subflux[k] += noah->grid[i].avgsubflux[k] * 1000.0 * 24.0 * 3600.0 / pihm->elem[i].topo.area;
            //}
        }

        if (solar > 1.0)
        {
            metv->tday += sfctmp;
            metv->q2d += noah->grid[i].q2sat - noah->grid[i].q2;
            metv->pa += noah->grid[i].sfcprs;
            metv->swavgfd += solar;
            metv->par += solar * RAD2PAR;
            daylight_counter[i]++;
        }
        else
        {
            metv->tnight += sfctmp;
        }

        (counter[i])++;
    }

    if ((t - start_time) % 86400 == 0 && t > start_time)
    {
        rawtime = (int) (t - 86400);
        timestamp = gmtime (&rawtime);
        spa.year = timestamp->tm_year + 1900;
        spa.month = timestamp->tm_mon + 1;
        spa.day = timestamp->tm_mday;
        spa.hour = timestamp->tm_hour;
        spa.minute = timestamp->tm_min;
        spa.second = timestamp->tm_sec;

        spa.timezone = 0;
        spa.delta_t = 67;
        spa.delta_ut1 = 0;
        spa.atmos_refract = 0.5667;

        spa.longitude = bgc->grid[0].sitec.lon;
        spa.latitude = bgc->grid[0].sitec.lat;
        spa.elevation = 0.;
        for (i = 0; i < pihm->numele; i++)
        {
            spa.elevation = spa.elevation + (double)pihm->elem[i].topo.zmax;
        }
        spa.elevation = spa.elevation / (double)pihm->numele;
        /*
         * Calculate surface pressure based on FAO 1998 method (Narasimhan 2002) 
         */
        spa.pressure = 1013.25 * pow ((293. - 0.0065 * spa.elevation) / 293., 5.26);
        spa.temperature = noah->genprmt.tbot_data;

        spa.function = SPA_ZA_RTS;
        spa_result = spa_calculate (&spa);

        /* daylength (s) */
        dayl = (spa.sunset - spa.sunrise) * 3600.0;
        dayl = (dayl < 0.0) ? (dayl + 24.0 * 3600.0) : dayl;

        rawtime = rawtime - 24 * 3600;
        timestamp = gmtime (&rawtime);
        spa.year = timestamp->tm_year + 1900;
        spa.month = timestamp->tm_mon + 1;
        spa.day = timestamp->tm_mday;
        spa.hour = timestamp->tm_hour;
        spa.minute = timestamp->tm_min;
        spa.second = timestamp->tm_sec;
        spa_result = spa_calculate (&spa);
        prev_dayl = (spa.sunset - spa.sunrise) * 3600.;
        prev_dayl = (prev_dayl < 0.0) ? (prev_dayl + 12.0 * 3600.0) : prev_dayl;

        for (i = 0; i < pihm->numele; i++)
        {
            metv = &(bgc->grid[i].metv);

            metv->dayl = dayl;
            metv->prev_dayl = prev_dayl;

            metv->tavg /= (double) counter[i];
            metv->tsoil /= (double) counter[i];
            metv->swc /= (double) counter[i];
            metv->soilw /= (double) counter[i];
            for (k = 0; k < 3; k++)
            {
                metv->subflux[k] /= (double) counter[i];
            }

            metv->tday /= (double) daylight_counter[i];
            metv->q2d /= (double) daylight_counter[i];
            metv->pa /= (double) daylight_counter[i];
            metv->swavgfd /= (double) daylight_counter[i];
            metv->par /= (double) daylight_counter[i];

            metv->tnight /= (double) (counter[i] - daylight_counter[i]);
        }

        DailyBgc (bgc, pihm->numele, t, dummy, first_balance);
        first_balance = 0;

        for (j = 0; j < bgc->ctrl.nprint; j++)
        {
            PrintData (&bgc->prtctrl[j], t, 86400, 1);
        }

        
        for (i = 0; i < pihm->numele; i++)
        {
            noah->grid[i].xlai = bgc->grid[i].epv.proj_lai;
            noah->grid[i].cmcmax = noah->grid[i].cmcfactr * noah->grid[i].xlai;

            metv = &(bgc->grid[i].metv);

            counter[i] = 0;
            daylight_counter[i] = 0;

            metv->tmax = -999.0;
            metv->tmin = 999.0;
            metv->tavg = 0.0;
            metv->tsoil = 0.0;
            metv->swc = 0.0;
            metv->soilw = 0.0;
            for (k = 0; k < 3; k++)
                metv->subflux[k] = 0.0;
            metv->tday = 0.0;
            metv->q2d = 0.0;
            metv->pa = 0.0;
            metv->swavgfd = 0.0;
            metv->par = 0.0;
            metv->tnight = 0.0;
        }
    }
}

void MapBgcOutput (char *simulation, bgc_struct bgc, int numele, char *outputdir)
{
    int             i, j;
    int             n;

    n = 0;

    for (i = 0; i < NUM_PRINT; i++)
    {
        if (bgc->ctrl.prtvrbl[i] > 0)
        {
            switch (i)
            {
                case LAI_CTRL:
                    sprintf (bgc->prtctrl[n].name, "%s%s.lai", outputdir, simulation);
                    bgc->prtctrl[n].intvl = 86400;
                    bgc->prtctrl[n].nvrbl = numele;
                    bgc->prtctrl[n].vrbl = (double **) malloc (bgc->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < bgc->prtctrl[n].nvrbl; j++)
                    {
                        bgc->prtctrl[n].vrbl[j] = &bgc->grid[i].epv.proj_lai;
                    }
                    n++;
                    break;
                case VEGC_CTRL:
                    sprintf (bgc->prtctrl[n].name, "%s%s.vegc", outputdir, simulation);
                    bgc->prtctrl[n].intvl = 86400;
                    bgc->prtctrl[n].nvrbl = numele;
                    bgc->prtctrl[n].vrbl = (double **) malloc (bgc->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < bgc->prtctrl[n].nvrbl; j++)
                    {
                        bgc->prtctrl[n].vrbl[j] = &bgc->grid[i].summary.vegc;
                    }
                    n++;
                    break;
                case LITRC_CTRL:
                    sprintf (bgc->prtctrl[n].name, "%s%s.litrc", outputdir, simulation);
                    bgc->prtctrl[n].intvl = 86400;
                    bgc->prtctrl[n].nvrbl = numele;
                    bgc->prtctrl[n].vrbl = (double **) malloc (bgc->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < bgc->prtctrl[n].nvrbl; j++)
                    {
                        bgc->prtctrl[n].vrbl[j] = &bgc->grid[i].summary.litrc;
                    }
                    n++;
                    break;
                case SOILC_CTRL:
                    sprintf (bgc->prtctrl[n].name, "%s%s.soilc", outputdir, simulation);
                    bgc->prtctrl[n].intvl = 86400;
                    bgc->prtctrl[n].nvrbl = numele;
                    bgc->prtctrl[n].vrbl = (double **) malloc (bgc->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < bgc->prtctrl[n].nvrbl; j++)
                    {
                        bgc->prtctrl[n].vrbl[j] = &bgc->grid[i].summary.soilc;
                    }
                    n++;
                    break;
                case TOTALC_CTRL:
                    sprintf (bgc->prtctrl[n].name, "%s%s.totalc", outputdir, simulation);
                    bgc->prtctrl[n].intvl = 86400;
                    bgc->prtctrl[n].nvrbl = numele;
                    bgc->prtctrl[n].vrbl = (double **) malloc (bgc->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < bgc->prtctrl[n].nvrbl; j++)
                    {
                        bgc->prtctrl[n].vrbl[j] = &bgc->grid[i].summary.totalc;
                    }
                    n++;
                    break;
                case NPP_CTRL:
                    sprintf (bgc->prtctrl[n].name, "%s%s.npp", outputdir, simulation);
                    bgc->prtctrl[n].intvl = 86400;
                    bgc->prtctrl[n].nvrbl = numele;
                    bgc->prtctrl[n].vrbl = (double **) malloc (bgc->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < bgc->prtctrl[n].nvrbl; j++)
                    {
                        bgc->prtctrl[n].vrbl[j] = &bgc->grid[i].summary.daily_npp;
                    }
                    n++;
                    break;
                case NEE_CTRL:
                    sprintf (bgc->prtctrl[n].name, "%s%s.nee", outputdir, simulation);
                    bgc->prtctrl[n].intvl = 86400;
                    bgc->prtctrl[n].nvrbl = numele;
                    bgc->prtctrl[n].vrbl = (double **) malloc (bgc->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < bgc->prtctrl[n].nvrbl; j++)
                    {
                        bgc->prtctrl[n].vrbl[j] = &bgc->grid[i].summary.daily_nee;
                    }
                    n++;
                    break;
                case GPP_CTRL:
                    sprintf (bgc->prtctrl[n].name, "%s%s.gpp", outputdir, simulation);
                    bgc->prtctrl[n].intvl = 86400;
                    bgc->prtctrl[n].nvrbl = numele;
                    bgc->prtctrl[n].vrbl = (double **) malloc (bgc->prtctrl[n].nvrbl * sizeof (double *));
                    for (j = 0; j < bgc->prtctrl[n].nvrbl; j++)
                    {
                        bgc->prtctrl[n].vrbl[j] = &bgc->grid[i].summary.daily_gpp;
                    }
                    n++;
                    break;
            }
        }
    }

    bgc->ctrl.nprint = n;
}

void ReadBinFile (ts_struct *ts, char *fn, int numele)
{
    FILE           *fid;
    double          dtime;
    int             i, j;

    fid = fopen (fn, "rb");
    CheckFile (fid, fn);

    fseek (fid, 0L, SEEK_END);
    ts->length = (int) (ftell (fid) / (numele + 1) / 8);
    ts->ftime = (int *) malloc (ts->length * sizeof (int));
    ts->data = (double **) malloc (ts->length * sizeof (double *));

    rewind (fid);
    for (j = 0; j < ts->length; j++)
    {
        ts->data[j] = (double *) malloc (numele * sizeof (double));
        fread (&dtime, sizeof (double), 1, fid);
        ts->ftime[j] = (int) dtime;
        for (i = 0; i < numele; i++)
        {
            fread (&ts->data[j][i], sizeof (double), 1, fid);
        }
    }

    fclose (fid);
}

void ReadAnnFile (ts_struct *ts, char *fn)
{
    FILE           *fid;
    time_t          rawtime;
    struct tm      *timeinfo;
    char            cmdstr[MAXSTRING];
    int             i;

    timeinfo = (struct tm *)malloc (sizeof (struct tm));

    fid = fopen (fn, "r");
    CheckFile (fid, fn);

    ts->length = CountLine (fid, 1, "EOF");
    ts->ftime = (int *) malloc (ts->length * sizeof (int));
    ts->data = (double **) malloc (ts->length * sizeof (double *));

    FindLine (fid, "BOF");
    for (i = 0; i < ts->length; i++)
    {
        ts->data[i] = (double *) malloc (sizeof (double));
        NextLine (fid, cmdstr);
        sscanf (cmdstr, "%d %lf", &timeinfo->tm_year, &ts->data[i][0]);
        timeinfo->tm_year = timeinfo->tm_year - 1900;
        timeinfo->tm_mon = 0;
        timeinfo->tm_mday = 1;
        timeinfo->tm_hour = 0;
        timeinfo->tm_min = 0;
        timeinfo->tm_sec = 0;
        rawtime = timegm (timeinfo);
        ts->ftime[i] = (int) rawtime;
    }
    
    fclose (fid);
}
