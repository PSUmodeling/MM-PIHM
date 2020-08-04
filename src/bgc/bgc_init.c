#include "pihm.h"

void InitBgc(const epctbl_struct *epctbl, const calib_struct *calib,
    elem_struct elem[])
{
    int             i;

    pihm_printf(VL_VERBOSE, "BGC: Initializing BGC structures\n");

#if defined(_LUMPEDBGC_)
    i = LUMPEDBGC;
#else
    for (i = 0; i < nelem; i++)
#endif
    {
        int             epc_ind;
        epconst_struct *epc;

        epc_ind = elem[i].attrib.lc;

        if (epc_ind != IGBP_ENF && epc_ind != IGBP_EBF && epc_ind != IGBP_DNF &&
            epc_ind != IGBP_DBF && epc_ind != IGBP_GRASS &&
            epc_ind != IGBP_CLOSE_SHRUB && epc_ind != IGBP_OPEN_SHRUB)
        {
            pihm_printf(VL_ERROR,
                "Error: Land cover type %d not defined in Flux-PIHM-BGC.\n",
                epc_ind);
            pihm_exit(EXIT_FAILURE);
        }

        epc_ind--;
        epc = &elem[i].epc;

        epc->woody = epctbl->woody[epc_ind];
        epc->evergreen = epctbl->evergreen[epc_ind];
        epc->c3_flag = epctbl->c3_flag[epc_ind];
        epc->phenology_flag = epctbl->phenology_flag[epc_ind];
        epc->onday = epctbl->onday[epc_ind];
        epc->offday = epctbl->offday[epc_ind];
        epc->transfer_days = epctbl->transfer_days[epc_ind];
        epc->litfall_days = epctbl->litfall_days[epc_ind];
        epc->leaf_turnover = epctbl->leaf_turnover[epc_ind];
        epc->froot_turnover = epctbl->froot_turnover[epc_ind];
        epc->livewood_turnover = epctbl->livewood_turnover[epc_ind];
        epc->daily_mortality_turnover =
            epctbl->daily_mortality_turnover[epc_ind] * calib->mortality;
        epc->daily_fire_turnover = epctbl->daily_fire_turnover[epc_ind];
        epc->alloc_frootc_leafc = epctbl->alloc_frootc_leafc[epc_ind];
        epc->alloc_newstemc_newleafc = epctbl->alloc_newstemc_newleafc[epc_ind];
        epc->alloc_newlivewoodc_newwoodc =
            epctbl->alloc_newlivewoodc_newwoodc[epc_ind];
        epc->alloc_crootc_stemc = epctbl->alloc_crootc_stemc[epc_ind];
        epc->alloc_prop_curgrowth = epctbl->alloc_prop_curgrowth[epc_ind];
        epc->avg_proj_sla = epctbl->avg_proj_sla[epc_ind] * calib->sla;
        epc->sla_ratio = epctbl->sla_ratio[epc_ind];
        epc->lai_ratio = epctbl->lai_ratio[epc_ind];
        epc->ext_coef = epctbl->ext_coef[epc_ind];
        epc->flnr = epctbl->flnr[epc_ind];
        epc->psi_open = epctbl->psi_open[epc_ind];
        epc->psi_close = epctbl->psi_close[epc_ind];
        epc->vpd_open = epctbl->vpd_open[epc_ind];
        epc->vpd_close = epctbl->vpd_close[epc_ind];
        epc->froot_cn = epctbl->froot_cn[epc_ind];
        epc->leaf_cn = epctbl->leaf_cn[epc_ind];
        epc->livewood_cn = epctbl->livewood_cn[epc_ind];
        epc->deadwood_cn = epctbl->deadwood_cn[epc_ind];
        epc->leaflitr_cn = epctbl->leaflitr_cn[epc_ind];
        epc->leaflitr_flab = epctbl->leaflitr_flab[epc_ind];
        epc->leaflitr_fucel = epctbl->leaflitr_fucel[epc_ind];
        epc->leaflitr_fscel = epctbl->leaflitr_fscel[epc_ind];
        epc->leaflitr_flig = epctbl->leaflitr_flig[epc_ind];
        epc->frootlitr_flab = epctbl->frootlitr_flab[epc_ind];
        epc->frootlitr_fucel = epctbl->frootlitr_fucel[epc_ind];
        epc->frootlitr_fscel = epctbl->frootlitr_fscel[epc_ind];
        epc->frootlitr_flig = epctbl->frootlitr_flig[epc_ind];
        epc->deadwood_fucel = epctbl->deadwood_fucel[epc_ind];
        epc->deadwood_fscel = epctbl->deadwood_fscel[epc_ind];
        epc->deadwood_flig = epctbl->deadwood_flig[epc_ind];
    }
}

void InitBgcVar(elem_struct elem[], river_struct river[], N_Vector CV_Y)
{
    int             i;

#if defined(_LUMPEDBGC_)
    i = LUMPEDBGC;
#else
    for (i = 0; i < nelem; i++)
#endif
    {
        RestartInput(&elem[i].restart_input, &elem[i].epv, &elem[i].cs,
            &elem[i].ns);

        ZeroSrcSnk(&elem[i].cs, &elem[i].ns, &elem[i].summary,
            &elem[i].solute[0]);
        elem[i].epv.annavg_t2m = elem[i].ps.tbot;

#if !defined(_LUMPEDBGC_)
        NV_Ith(CV_Y, SOLUTE_SOIL(i, 0)) = elem[i].ns.sminn;

        elem[i].nt.sminn0 = elem[i].ns.sminn;
#endif
    }

#if defined(_LUMPEDBGC_)
    NV_Ith(CV_Y, LUMPEDBGC_SMINN) = elem[LUMPEDBGC].ns.sminn;
    elem[LUMPEDBGC].nt.sminn0 = elem[i].ns.sminn;
#endif

#if !defined(_LUMPEDBGC_) && !defined(_LEACHING_)
    for (i = 0; i < nriver; i++)
    {
        river[i].ns.streamn = river[i].restart_input.streamn;

        NV_Ith(CV_Y, SOLUTE_RIVER(i, 0)) = river[i].ns.streamn;
    }
#endif
}
