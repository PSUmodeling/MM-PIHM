#include "pihm.h"

void InitBGC (elem_struct *elem, int numele, river_struct *riv, int numriv, const epctbl_struct *epctbl, const ctrl_struct *ctrl)
{
    int             i;
    int             epc_ind;

    /* Detect if model is running in ensemble mode */
    if (verbose_mode)
    {
        printf ("BGC: Initializing BGC structures\n");
    }

    for (i = 0; i < numele; i++)
    {
        epc_ind = elem[i].attrib.lc_type - 1;

        if (epc_ind != ENF - 1 &&
            epc_ind != EBF - 1 &&
            epc_ind != DNF - 1 &&
            epc_ind != DBF - 1 &&
            epc_ind != GRASS - 1 &&
            epc_ind != CLOSE_SHRUB - 1 &&
            epc_ind != OPEN_SHRUB - 1)
        {
            fprintf (stderr, "Error: Land cover type %d not been defined in Flux-PIHM-BGC.\n",
                elem[i].attrib.lc_type);
            PIHMError (1);
        }

        elem[i].epc.woody = epctbl->woody[epc_ind];
        elem[i].epc.evergreen = epctbl->evergreen[epc_ind];
        elem[i].epc.c3_flag = epctbl->c3_flag[epc_ind];
        elem[i].epc.phenology_flag = epctbl->phenology_flag[epc_ind];
        elem[i].epc.onday = epctbl->onday[epc_ind];
        elem[i].epc.offday = epctbl->offday[epc_ind];
        elem[i].epc.transfer_days = epctbl->transfer_days[epc_ind];
        elem[i].epc.litfall_days = epctbl->litfall_days[epc_ind];
        elem[i].epc.leaf_turnover = epctbl->leaf_turnover[epc_ind];
        elem[i].epc.froot_turnover = epctbl->froot_turnover[epc_ind];
        elem[i].epc.livewood_turnover = epctbl->livewood_turnover[epc_ind];
        elem[i].epc.daily_mortality_turnover = epctbl->daily_mortality_turnover[epc_ind];
        elem[i].epc.daily_fire_turnover = epctbl->daily_fire_turnover[epc_ind];
        elem[i].epc.alloc_frootc_leafc = epctbl->alloc_frootc_leafc[epc_ind];
        elem[i].epc.alloc_newstemc_newleafc = epctbl->alloc_newstemc_newleafc[epc_ind];
        elem[i].epc.alloc_newlivewoodc_newwoodc = epctbl->alloc_newlivewoodc_newwoodc[epc_ind];
        elem[i].epc.alloc_crootc_stemc = epctbl->alloc_crootc_stemc[epc_ind];
        elem[i].epc.alloc_prop_curgrowth = epctbl->alloc_prop_curgrowth[epc_ind];
        elem[i].epc.avg_proj_sla = epctbl->avg_proj_sla[epc_ind];
        elem[i].epc.sla_ratio = epctbl->sla_ratio[epc_ind];
        elem[i].epc.lai_ratio = epctbl->lai_ratio[epc_ind];
        elem[i].epc.ext_coef = epctbl->ext_coef[epc_ind];
        elem[i].epc.flnr = epctbl->flnr[epc_ind];
        elem[i].epc.psi_open = epctbl->psi_open[epc_ind];
        elem[i].epc.psi_close = epctbl->psi_close[epc_ind];
        elem[i].epc.vpd_open = epctbl->vpd_open[epc_ind];
        elem[i].epc.vpd_close = epctbl->vpd_close[epc_ind];
        elem[i].epc.froot_cn = epctbl->froot_cn[epc_ind];
        elem[i].epc.leaf_cn = epctbl->leaf_cn[epc_ind];
        elem[i].epc.livewood_cn = epctbl->livewood_cn[epc_ind];
        elem[i].epc.deadwood_cn = epctbl->deadwood_cn[epc_ind];
        elem[i].epc.leaflitr_cn = epctbl->leaflitr_cn[epc_ind];
        elem[i].epc.leaflitr_flab = epctbl->leaflitr_flab[epc_ind];
        elem[i].epc.leaflitr_fucel = epctbl->leaflitr_fucel[epc_ind];
        elem[i].epc.leaflitr_fscel = epctbl->leaflitr_fscel[epc_ind];
        elem[i].epc.leaflitr_flig = epctbl->leaflitr_flig[epc_ind];
        elem[i].epc.frootlitr_flab = epctbl->frootlitr_flab[epc_ind];
        elem[i].epc.frootlitr_fucel = epctbl->frootlitr_fucel[epc_ind];
        elem[i].epc.frootlitr_fscel = epctbl->frootlitr_fscel[epc_ind];
        elem[i].epc.frootlitr_flig = epctbl->frootlitr_flig[epc_ind];
        elem[i].epc.deadwood_fucel = epctbl->deadwood_fucel[epc_ind];
        elem[i].epc.deadwood_fscel = epctbl->deadwood_fscel[epc_ind];
        elem[i].epc.deadwood_flig = epctbl->deadwood_flig[epc_ind];

        if (ctrl->bgc_spinup)
        {
            InitElemStor (&elem[i].stor, ctrl->spinupstart, ctrl->spinupend);
        }
    }

    for (i = 0; i < numriv; i++)
    {
        if (ctrl->bgc_spinup)
        {
            InitRiverStor (&riv[i].stor, ctrl->spinupstart, ctrl->spinupend);
        }
    }
}

void InitBGCVar (elem_struct *elem, int numele, river_struct *riv, int numriv, cinit_struct cinit, cstate_struct cs, nstate_struct ns, char *fn, int spinup)
{
    int             i;
    FILE           *init_file;

    /* Read initial conditions */
    if (!spinup)
    {
        init_file = fopen (fn, "rb");

        if (NULL == init_file)
        {
            fprintf (stderr, "Error reading %s.\n", fn);
            PIHMError (1);
        }

        if (verbose_mode)
        {
            printf ("Reading %s...", fn);
        }

        for (i = 0; i < numele; i++)
        {
            fread (&elem[i].restart_input, sizeof (bgc_ic_struct), 1, init_file);

            restart_input (&elem[i].cs, &elem[i].ns, &elem[i].epv,
                &elem[i].restart_input);

            /* Calculate LAI for the coupling with Noah */
            elem[i].ps.proj_lai = elem[i].cs.leafc * elem[i].epc.avg_proj_sla;
            elem[i].epv.annavg_t2m = elem[i].ps.tbot;
        }

        for (i = 0; i < numriv; i++)
        {
            fread (&riv[i].ns.sminn, sizeof (double), 1, init_file);
        }

        fclose (init_file);
    }
    else
    {
        for (i = 0; i < numele; i++)
        {
            elem[i].cinit = cinit;
            elem[i].cs = cs;
            elem[i].ns = ns;

            elem[i].ns.cwdn = cs.cwdc / elem[i].epc.deadwood_cn;
            elem[i].ns.litr2n = cs.litr2c / elem[i].epc.leaflitr_cn;
            elem[i].ns.litr3n = cs.litr3c / elem[i].epc.leaflitr_cn;
            elem[i].ns.litr4n = cs.litr4c / elem[i].epc.leaflitr_cn;
            elem[i].ns.soil1n = cs.soil1c / SOIL1_CN;
            elem[i].ns.soil2n = cs.soil2c / SOIL2_CN;
            elem[i].ns.soil3n = cs.soil3c / SOIL3_CN;
            elem[i].ns.soil4n = cs.soil4c / SOIL4_CN;

            if (elem[i].epc.evergreen)
            {
                elem[i].epv.dormant_flag = 0.0;
            }
            else
            {
                elem[i].epv.dormant_flag = 1.0;
            }

            //elem[i].epv.days_active = 0.;
            elem[i].epv.onset_flag = 0.0;
            elem[i].epv.onset_counter = 0.0;
            elem[i].epv.onset_gddflag = 0.0;
            elem[i].epv.onset_fdd = 0.0;
            elem[i].epv.onset_gdd = 0.0;
            elem[i].epv.onset_swi = 0.0;
            elem[i].epv.offset_flag = 0.0;
            elem[i].epv.offset_counter = 0.0;
            elem[i].epv.offset_fdd = 0.0;
            elem[i].epv.offset_swi = 0.0;
            elem[i].epv.annavg_t2m = elem[i].ps.tbot;

            firstday (&elem[i].epc, &elem[i].cinit, &elem[i].epv, &elem[i].cs, &elem[i].ns);
        }

        for (i = 0; i < numriv; i++)
        {
            riv[i].ns.sminn = 0.0;
        }
    }

    for (i = 0; i < numele; i++)
    {
        zero_srcsnk (&elem[i].cs, &elem[i].ns, &elem[i].summary);
    }

    for (i = 0; i < numriv; i++)
    {
       //riv[i].nleached_snk = 0.0;
       riv[i].nf.sminn_leached = 0.0;
    }
}
