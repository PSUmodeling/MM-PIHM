#include "pihm.h"

void ReadBgc(const char fn[], char co2_fn[], char ndep_fn[], ctrl_struct *ctrl,
    co2control_struct *co2, ndepcontrol_struct *ndepctrl,
    cninit_struct * cninit)
{
    FILE           *bgc_file;
    char            cmdstr[MAXSTRING];
    int             lno = 0;
    int             acc_flag = 0;

    /* Read bgc simulation control file */
    bgc_file = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    FindLine(bgc_file, "RESTART", &lno, fn);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%d", &ctrl->read_bgc_restart);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%d", &ctrl->write_bgc_restart);

    FindLine(bgc_file, "TIME_DEFINE", &lno, fn);
    NextLine(bgc_file, cmdstr, &lno);
    if (spinup_mode)
    {
        sscanf(cmdstr, "%d", &acc_flag);
        spinup_mode = (acc_flag == 1) ? ACC_SPINUP_MODE : SPINUP_MODE;
    }

    FindLine(bgc_file, "CO2_CONTROL", &lno, fn);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%d", &co2->varco2);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &co2->co2ppm);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%s", co2_fn);

    FindLine(bgc_file, "NDEP_CONTROL", &lno, fn);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%d", &ndepctrl->varndep);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &ndepctrl->ndep);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &ndepctrl->nfix);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%s", ndep_fn);

    FindLine(bgc_file, "C_STATE", &lno, fn);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->max_leafc);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->max_stemc);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->cwdc);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->litr1c);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->litr2c);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->litr3c);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->litr4c);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->soil1c);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->soil2c);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->soil3c);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->soil4c);

    FindLine(bgc_file, "N_STATE", &lno, fn);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->litr1n);
    NextLine(bgc_file, cmdstr, &lno);
    sscanf(cmdstr, "%lf", &cninit->sminn);

    FindLine(bgc_file, "OUTPUT_CONTROL", &lno, fn);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[LAI_CTRL] = ReadPrintCtrl(cmdstr, "LAI", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[NPP_CTRL] = ReadPrintCtrl(cmdstr, "NPP", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[NEP_CTRL] = ReadPrintCtrl(cmdstr, "NEP", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[NEE_CTRL] = ReadPrintCtrl(cmdstr, "NEE", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[GPP_CTRL] = ReadPrintCtrl(cmdstr, "GPP", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[MR_CTRL] = ReadPrintCtrl(cmdstr, "MR", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[GR_CTRL] = ReadPrintCtrl(cmdstr, "GR", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[HR_CTRL] = ReadPrintCtrl(cmdstr, "HR", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[FIRE_CTRL] = ReadPrintCtrl(cmdstr, "FIRE", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[LITFALLC_CTRL] = ReadPrintCtrl(cmdstr, "LITFALLC", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[VEGC_CTRL] = ReadPrintCtrl(cmdstr, "VEGC", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[AGC_CTRL] = ReadPrintCtrl(cmdstr, "AGC", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[LITRC_CTRL] = ReadPrintCtrl(cmdstr, "LITRC", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[SOILC_CTRL] = ReadPrintCtrl(cmdstr, "SOILC", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[TOTALC_CTRL] = ReadPrintCtrl(cmdstr, "TOTALC", fn, lno);

    NextLine(bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[SMINN_CTRL] = ReadPrintCtrl(cmdstr, "SMINN", fn, lno);

    fclose(bgc_file);
}

void ReadEpc(epctbl_struct *epctbl)
{
    int             i;
    char            fn[MAXSTRING];
    double          t1, t2, t3, t4, r1;
    FILE           *epc_file;
    char            cmdstr[MAXSTRING];

    epctbl->woody = (int *)malloc(NLCTYPE * sizeof(int));
    epctbl->evergreen = (int *)malloc(NLCTYPE * sizeof(int));
    epctbl->c3_flag = (int *)malloc(NLCTYPE * sizeof(int));
    epctbl->phenology_flag = (int *)malloc(NLCTYPE * sizeof(int));
    epctbl->onday = (int *)malloc(NLCTYPE * sizeof(int));
    epctbl->offday = (int *)malloc(NLCTYPE * sizeof(int));
    epctbl->transfer_days = (int *)malloc(NLCTYPE * sizeof(int));
    epctbl->litfall_days = (int *)malloc(NLCTYPE * sizeof(int));
    epctbl->leaf_turnover = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->froot_turnover = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->livewood_turnover = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->daily_mortality_turnover =
        (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->daily_fire_turnover = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->alloc_frootc_leafc = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->alloc_newstemc_newleafc =
        (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->alloc_newlivewoodc_newwoodc =
        (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->alloc_crootc_stemc = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->alloc_prop_curgrowth = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->avg_proj_sla = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->sla_ratio = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->lai_ratio = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->ext_coef = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->flnr = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->psi_open = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->psi_close = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->vpd_open = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->vpd_close = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->froot_cn = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->leaf_cn = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->livewood_cn = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->deadwood_cn = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->leaflitr_cn = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->leaflitr_flab = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->leaflitr_fucel = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->leaflitr_fscel = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->leaflitr_flig = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->frootlitr_flab = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->frootlitr_fucel = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->frootlitr_fscel = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->frootlitr_flig = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->deadwood_fucel = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->deadwood_fscel = (double *)malloc(NLCTYPE * sizeof(double));
    epctbl->deadwood_flig = (double *)malloc(NLCTYPE * sizeof(double));

    /* Read epc files */
    pihm_printf(VL_VERBOSE, "\nRead ecophysiological constant files\n");

    for (i = 0; i < NLCTYPE; i++)
    {
        switch (i + 1)
        {
            case IGBP_ENF:
                strcpy(fn, "input/epc/enf.epc");
                epc_file = pihm_fopen(fn, "r");
                break;
            case IGBP_EBF:
                strcpy(fn, "input/epc/ebf.epc");
                epc_file = pihm_fopen(fn, "r");
                break;
            case IGBP_DNF:
                strcpy(fn, "input/epc/dnf.epc");
                epc_file = pihm_fopen(fn, "r");
                break;
            case IGBP_DBF:
                strcpy(fn, "input/epc/dbf.epc");
                epc_file = pihm_fopen(fn, "r");
                break;
            case IGBP_GRASS:
                strcpy(fn, "input/epc/c3grass.epc");
                epc_file = pihm_fopen(fn, "r");
                break;
            case IGBP_CLOSE_SHRUB:
                strcpy(fn, "input/epc/shrub.epc");
                epc_file = pihm_fopen(fn, "r");
                break;
            case IGBP_OPEN_SHRUB:
                strcpy(fn, "input/epc/shrub.epc");
                epc_file = pihm_fopen(fn, "r");
                break;
            default:
                strcpy(fn, "N/A");
                epc_file = NULL;
        }

        if (strcasecmp(fn, "N/A") != 0)
        {
            pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

            /* Skip header file */
            fgets(cmdstr, MAXSTRING, epc_file);
            /* Read epc */
            /* Woody/non-woody flag */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%d", &epctbl->woody[i]);
            /* Evergreen/deciduous flag */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%d", &epctbl->evergreen[i]);
            /* C3/C4 flag */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%d", &epctbl->c3_flag[i]);
            /* Transfer days */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%d", &epctbl->transfer_days[i]);
            /* Litter fall days */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%d", &epctbl->litfall_days[i]);
            /* Leaf turnover */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->leaf_turnover[i]);
            /* Force leaf turnover fraction to 1.0 if deciduous */
            if (!epctbl->evergreen[i])
            {
                epctbl->leaf_turnover[i] = 1.0;
            }
            epctbl->froot_turnover[i] = epctbl->leaf_turnover[i];
            /* Live wood turnover */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->livewood_turnover[i]);
            /* Whole-plant mortality */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &t1);
            epctbl->daily_mortality_turnover[i] = t1 / 365.0;
            /* Fire mortality */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &t1);
            epctbl->daily_fire_turnover[i] = t1 / 365.0;
            /* Froot C:leaf C */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->alloc_frootc_leafc[i]);
            /* New stem C:new leaf C */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->alloc_newstemc_newleafc[i]);
            /* New livewood C:new wood C */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->alloc_newlivewoodc_newwoodc[i]);
            /* Croot C:stem C */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->alloc_crootc_stemc[i]);
            /* New growth:storage growth */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->alloc_prop_curgrowth[i]);
            /* Force storage growth to 0.0 if evergreen (following CLM-CN) */
            if (epctbl->evergreen[i])
            {
                epctbl->alloc_prop_curgrowth[i] = 1.0;
            }
            /* Average leaf C:N */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->leaf_cn[i]);
            /* Leaf litter C:N */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->leaflitr_cn[i]);
            /* Test for leaflitter C:N > leaf C:N */
            if (epctbl->leaflitr_cn[i] < epctbl->leaf_cn[i])
            {
                pihm_printf(VL_ERROR,
                    "Error: leaf litter C:N must be >= leaf C:N.\n");
                pihm_printf(VL_ERROR, "Change the values in epc file %s.\n", fn);
                pihm_exit(EXIT_FAILURE);
            }
            /* Initial fine root C:N */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->froot_cn[i]);
            /* Initial livewood C:N */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->livewood_cn[i]);
            /* Initial deadwood C:N */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->deadwood_cn[i]);
            /* Test for deadwood C:N > livewood C:N */
            if (epctbl->deadwood_cn[i] < epctbl->livewood_cn[i])
            {
                pihm_printf(VL_ERROR,
                    "Error: livewood C:N must be >= deadwood C:N.\n");
                pihm_printf(VL_ERROR, "Change the values in epc file %s.\n", fn);
                pihm_exit(EXIT_FAILURE);
            }
            /* Leaf litter labile proportion */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &t1);
            epctbl->leaflitr_flab[i] = t1;
            /* Leaf litter cellulose proportion */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &t2);
            /* Leaf litter lignin proportion */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &t3);
            epctbl->leaflitr_flig[i] = t3;

            /* Test for litter fractions sum to 1.0 */
            if (fabs(t1 + t2 + t3 - 1.0) > FLT_COND_TOL)
            {
                pihm_printf(VL_ERROR,
                    "Error: leaf litter proportions of labile, cellulose, "
                    "and lignin\n");
                pihm_printf(VL_ERROR,
                    "must sum to 1.0. Check epc file %s and try again.\n", fn);
                pihm_exit(EXIT_FAILURE);
            }
            /* Calculate shielded and unshielded cellulose fraction */
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
            /* Froot litter labile proportion */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &t1);
            epctbl->frootlitr_flab[i] = t1;
            /* Froot litter cellulose proportion */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &t2);
            /* Froot litter lignin proportion */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &t3);
            epctbl->frootlitr_flig[i] = t3;

            /* Test for litter fractions sum to 1.0 */
            if (fabs(t1 + t2 + t3 - 1.0) > FLT_COND_TOL)
            {
                pihm_printf(VL_ERROR,
                    "Error: froot litter proportions of labile, cellulose, "
                    "and lignin\n");
                pihm_printf(VL_ERROR,
                    "must sum to 1.0. Check epc file %s and try again.\n", fn);
                pihm_exit(EXIT_FAILURE);
            }
            /* Calculate shielded and unshielded cellulose fraction */
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
            /* Dead wood cellulose */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &t1);
            /* Dead wood lignin */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &t2);
            epctbl->deadwood_flig[i] = t2;
            /* Test for litter fractions sum to 1.0 */
            if (fabs(t1 + t2 - 1.0) > FLT_COND_TOL)
            {
                pihm_printf(VL_ERROR,
                    "Error: deadwood proportions of cellulose and lignin "
                    "must sum\n");
                pihm_printf(VL_ERROR,
                    "to 1.0. Check epc file %s and try again.\n", fn);
                pihm_exit(EXIT_FAILURE);
            }
            /* Calculate shielded and unshielded cellulose fraction */
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
            /* Canopy light ext coef */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->ext_coef[i]);
            /* All to projected LA ratio */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->lai_ratio[i]);
            /* Canopy average projected specific leaf area */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->avg_proj_sla[i]);
            /* Sunlit SLA ratio */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->sla_ratio[i]);
            /* Rubisco N fraction */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->flnr[i]);
            /* Psi_sat */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->psi_open[i]);
            /* Psi_close */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->psi_close[i]);
            /* Vpd_max */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->vpd_open[i]);
            /* Vpd_min */
            fgets(cmdstr, MAXSTRING, epc_file);
            sscanf(cmdstr, "%lf", &epctbl->vpd_close[i]);

            fclose(epc_file);
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

void ReadAnnualFile(const char fn[], tsdata_struct *ts)
{
    FILE           *fid;
    char            timestr[MAXSTRING];
    char            cmdstr[MAXSTRING];
    int             k;
    int             match;
    int             lno = 0;

    fid = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    ts->length = CountLine(fid, cmdstr, 1, "EOF");
    ts->ftime = (int *)malloc(ts->length * sizeof(int));
    ts->data = (double **)malloc(ts->length * sizeof(double *));

    FindLine(fid, "BOF", &lno, fn);
    for (k = 0; k < ts->length; k++)
    {
        ts->data[k] = (double *)malloc(sizeof(double));
        NextLine(fid, cmdstr, &lno);
        match = sscanf(cmdstr, "%s %lf", timestr, &ts->data[k][0]);

        if (match != 2)
        {
            pihm_printf(VL_ERROR, "Error reading %s.\n", fn);
            pihm_printf(VL_ERROR,
                "Please check file format near Line %s.\n", lno);
            pihm_exit(EXIT_FAILURE);
        }

        ts->ftime[k] = StrTime(timestr);
    }

    fclose(fid);
}
