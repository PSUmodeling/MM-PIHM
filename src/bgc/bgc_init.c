#include "pihm.h"

void InitBgc(elem_struct *elem, const epctbl_struct *epctbl,
    const calib_struct *cal)
{
    int             i;
    int             epc_ind;

    PIHMprintf(VL_VERBOSE, "BGC: Initializing BGC structures\n");

#if defined(_LUMPED_)
    i = LUMPED;
#else
    for (i = 0; i < nelem; i++)
#endif
    {
        epc_ind = elem[i].attrib.lc_type - 1;

        if (epc_ind != IGBP_ENF - 1 &&
            epc_ind != IGBP_EBF - 1 &&
            epc_ind != IGBP_DNF - 1 &&
            epc_ind != IGBP_DBF - 1 &&
            epc_ind != IGBP_GRASS - 1 &&
            epc_ind != IGBP_CLOSE_SHRUB - 1 && epc_ind != IGBP_OPEN_SHRUB - 1)
        {
            PIHMprintf(VL_ERROR,
                "Error: Land cover type %d not defined in Flux-PIHM-BGC.\n",
                elem[i].attrib.lc_type);
            PIHMexit(EXIT_FAILURE);
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
        elem[i].epc.daily_mortality_turnover =
            epctbl->daily_mortality_turnover[epc_ind] * cal->mortality;
        elem[i].epc.daily_fire_turnover = epctbl->daily_fire_turnover[epc_ind];
        elem[i].epc.alloc_frootc_leafc = epctbl->alloc_frootc_leafc[epc_ind];
        elem[i].epc.alloc_newstemc_newleafc =
            epctbl->alloc_newstemc_newleafc[epc_ind];
        elem[i].epc.alloc_newlivewoodc_newwoodc =
            epctbl->alloc_newlivewoodc_newwoodc[epc_ind];
        elem[i].epc.alloc_crootc_stemc = epctbl->alloc_crootc_stemc[epc_ind];
        elem[i].epc.alloc_prop_curgrowth =
            epctbl->alloc_prop_curgrowth[epc_ind];
        elem[i].epc.avg_proj_sla = epctbl->avg_proj_sla[epc_ind] * cal->sla;
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
    }
}

void InitBgcVar(elem_struct *elem, river_struct *river, N_Vector CV_Y)
{
    int             i;

#if defined(_LUMPED_)
    i = LUMPED;
#else
    for (i = 0; i < nelem; i++)
#endif
    {
        RestartInput(&elem[i].cs, &elem[i].ns, &elem[i].epv,
            &elem[i].restart_input);

        ZeroSrcSnk(&elem[i].cs, &elem[i].ns, &elem[i].summary, &elem[i].nsol);
        elem[i].epv.annavg_t2m = elem[i].ps.tbot;

#if !defined(_LUMPED_)
        NV_Ith(CV_Y, SURFN(i)) = elem[i].ns.surfn;
        NV_Ith(CV_Y, SMINN(i)) = elem[i].ns.sminn;

        elem[i].nt.surfn0 = elem[i].ns.surfn;
        elem[i].nt.sminn0 = elem[i].ns.sminn;
#endif
    }

#if defined(_LUMPED_)
    NV_Ith(CV_Y, LUMPED_SMINN) = elem[LUMPED].ns.sminn;
    elem[LUMPED].nt.sminn0 = elem[i].ns.sminn;
#endif

#if !defined(_LUMPED_) && !defined(_LEACHING_)
    for (i = 0; i < nriver; i++)
    {
        river[i].ns.streamn = river[i].restart_input.streamn;
        river[i].ns.sminn = river[i].restart_input.sminn;
        river[i].nf.sminn_leached = 0.0;

        NV_Ith(CV_Y, STREAMN(i)) = river[i].ns.streamn;
        NV_Ith(CV_Y, RIVBEDN(i)) = river[i].ns.sminn;
    }
#endif
}
