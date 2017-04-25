#include "pihm.h"

void Mortality (const epconst_struct *epc, cstate_struct *cs,
    cflux_struct *cf, nstate_struct *ns, nflux_struct *nf, double dt)
{
    double          mort;
    double          fire;
    /* dead stem combustion proportion */
    const double    DSCP = 0.2;

    mort = epc->daily_mortality_turnover / DAYINSEC;
    fire = epc->daily_fire_turnover / DAYINSEC;

    /*
     * Non-fire mortality: these fluxes all enter litter or CWD pools
     */

    /* carbon fluxes due to mortality */
    /* mortality fluxes out of leaf and fine root pools */
    cf->m_leafc_to_litr1c = mort * cs->leafc * epc->leaflitr_flab;
    cf->m_leafc_to_litr2c = mort * cs->leafc * epc->leaflitr_fucel;
    cf->m_leafc_to_litr3c = mort * cs->leafc * epc->leaflitr_fscel;
    cf->m_leafc_to_litr4c = mort * cs->leafc * epc->leaflitr_flig;
    cf->m_frootc_to_litr1c = mort * cs->frootc * epc->frootlitr_flab;
    cf->m_frootc_to_litr2c = mort * cs->frootc * epc->frootlitr_fucel;
    cf->m_frootc_to_litr3c = mort * cs->frootc * epc->frootlitr_fscel;
    cf->m_frootc_to_litr4c = mort * cs->frootc * epc->frootlitr_flig;

    /* mortality fluxes out of storage and transfer pools */
    cf->m_leafc_storage_to_litr1c = mort * cs->leafc_storage;
    cf->m_frootc_storage_to_litr1c = mort * cs->frootc_storage;
    cf->m_livestemc_storage_to_litr1c = mort * cs->livestemc_storage;
    cf->m_deadstemc_storage_to_litr1c = mort * cs->deadstemc_storage;
    cf->m_livecrootc_storage_to_litr1c = mort * cs->livecrootc_storage;
    cf->m_deadcrootc_storage_to_litr1c = mort * cs->deadcrootc_storage;
    cf->m_leafc_transfer_to_litr1c = mort * cs->leafc_transfer;
    cf->m_frootc_transfer_to_litr1c = mort * cs->frootc_transfer;
    cf->m_livestemc_transfer_to_litr1c = mort * cs->livestemc_transfer;
    cf->m_deadstemc_transfer_to_litr1c = mort * cs->deadstemc_transfer;
    cf->m_livecrootc_transfer_to_litr1c = mort * cs->livecrootc_transfer;
    cf->m_deadcrootc_transfer_to_litr1c = mort * cs->deadcrootc_transfer;
    cf->m_gresp_storage_to_litr1c = mort * cs->gresp_storage;
    cf->m_gresp_transfer_to_litr1c = mort * cs->gresp_transfer;

    /* TREE-specific carbon fluxes */
    if (epc->woody)
    {
        cf->m_livestemc_to_cwdc = mort * cs->livestemc;
        cf->m_deadstemc_to_cwdc =
            (mort + (1.0 - DSCP) * fire) * cs->deadstemc;
        cf->m_livecrootc_to_cwdc = mort * cs->livecrootc;
        cf->m_deadcrootc_to_cwdc =
            (mort + (1.0 - DSCP) * fire) * cs->deadcrootc;
    }

    /* nitrogen fluxes due to mortality */
    /* mortality fluxes out of leaf and fine root pools */
    nf->m_leafn_to_litr1n = cf->m_leafc_to_litr1c / epc->leaf_cn;
    nf->m_leafn_to_litr2n = cf->m_leafc_to_litr2c / epc->leaf_cn;
    nf->m_leafn_to_litr3n = cf->m_leafc_to_litr3c / epc->leaf_cn;
    nf->m_leafn_to_litr4n = cf->m_leafc_to_litr4c / epc->leaf_cn;
    nf->m_frootn_to_litr1n = cf->m_frootc_to_litr1c / epc->froot_cn;
    nf->m_frootn_to_litr2n = cf->m_frootc_to_litr2c / epc->froot_cn;
    nf->m_frootn_to_litr3n = cf->m_frootc_to_litr3c / epc->froot_cn;
    nf->m_frootn_to_litr4n = cf->m_frootc_to_litr4c / epc->froot_cn;

    /* mortality fluxes out of storage and transfer */
    nf->m_leafn_storage_to_litr1n = mort * ns->leafn_storage;
    nf->m_frootn_storage_to_litr1n = mort * ns->frootn_storage;
    nf->m_livestemn_storage_to_litr1n = mort * ns->livestemn_storage;
    nf->m_deadstemn_storage_to_litr1n = mort * ns->deadstemn_storage;
    nf->m_livecrootn_storage_to_litr1n = mort * ns->livecrootn_storage;
    nf->m_deadcrootn_storage_to_litr1n = mort * ns->deadcrootn_storage;
    nf->m_leafn_transfer_to_litr1n = mort * ns->leafn_transfer;
    nf->m_frootn_transfer_to_litr1n = mort * ns->frootn_transfer;
    nf->m_livestemn_transfer_to_litr1n = mort * ns->livestemn_transfer;
    nf->m_deadstemn_transfer_to_litr1n = mort * ns->deadstemn_transfer;
    nf->m_livecrootn_transfer_to_litr1n = mort * ns->livecrootn_transfer;
    nf->m_deadcrootn_transfer_to_litr1n = mort * ns->deadcrootn_transfer;
    nf->m_retransn_to_litr1n = mort * ns->retransn;

    /* woody-specific nitrogen fluxes */
    if (epc->woody)
    {
        nf->m_livestemn_to_cwdn = cf->m_livestemc_to_cwdc / epc->deadwood_cn;
        nf->m_livestemn_to_litr1n =
            (mort * ns->livestemn) - nf->m_livestemn_to_cwdn;
        nf->m_deadstemn_to_cwdn = cf->m_deadstemc_to_cwdc / epc->deadwood_cn;
        nf->m_livecrootn_to_cwdn =
            cf->m_livecrootc_to_cwdc / epc->deadwood_cn;
        nf->m_livecrootn_to_litr1n =
            (mort * ns->livecrootn) - nf->m_livecrootn_to_cwdn;
        nf->m_deadcrootn_to_cwdn =
            cf->m_deadcrootc_to_cwdc / epc->deadwood_cn;
    }

    /* update state variables */
    /* this is the only place other than daily_state_update() routines where
     * state variables get changed.  Mortality is taken care of last and
     * given special treatment for state update so that it doesn't interfere
     * with the other fluxes that are based on proportions of state variables,
     * especially the phenological fluxes */
    /* CARBON mortality state variable update */
    /*   Leaf mortality */
    cs->litr1c += cf->m_leafc_to_litr1c * dt;
    cs->leafc -= cf->m_leafc_to_litr1c * dt;
    cs->litr2c += cf->m_leafc_to_litr2c * dt;
    cs->leafc -= cf->m_leafc_to_litr2c * dt;
    cs->litr3c += cf->m_leafc_to_litr3c * dt;
    cs->leafc -= cf->m_leafc_to_litr3c * dt;
    cs->litr4c += cf->m_leafc_to_litr4c * dt;
    cs->leafc -= cf->m_leafc_to_litr4c * dt;
    /*   Fine root mortality */
    cs->litr1c += cf->m_frootc_to_litr1c * dt;
    cs->frootc -= cf->m_frootc_to_litr1c * dt;
    cs->litr2c += cf->m_frootc_to_litr2c * dt;
    cs->frootc -= cf->m_frootc_to_litr2c * dt;
    cs->litr3c += cf->m_frootc_to_litr3c * dt;
    cs->frootc -= cf->m_frootc_to_litr3c * dt;
    cs->litr4c += cf->m_frootc_to_litr4c * dt;
    cs->frootc -= cf->m_frootc_to_litr4c * dt;
    /*   Storage and transfer mortality */
    cs->litr1c += cf->m_leafc_storage_to_litr1c * dt;
    cs->leafc_storage -= cf->m_leafc_storage_to_litr1c * dt;
    cs->litr1c += cf->m_frootc_storage_to_litr1c * dt;
    cs->frootc_storage -= cf->m_frootc_storage_to_litr1c * dt;
    cs->litr1c += cf->m_livestemc_storage_to_litr1c * dt;
    cs->livestemc_storage -= cf->m_livestemc_storage_to_litr1c * dt;
    cs->litr1c += cf->m_deadstemc_storage_to_litr1c * dt;
    cs->deadstemc_storage -= cf->m_deadstemc_storage_to_litr1c * dt;
    cs->litr1c += cf->m_livecrootc_storage_to_litr1c * dt;
    cs->livecrootc_storage -= cf->m_livecrootc_storage_to_litr1c * dt;
    cs->litr1c += cf->m_deadcrootc_storage_to_litr1c * dt;
    cs->deadcrootc_storage -= cf->m_deadcrootc_storage_to_litr1c * dt;
    cs->litr1c += cf->m_leafc_transfer_to_litr1c * dt;
    cs->leafc_transfer -= cf->m_leafc_transfer_to_litr1c * dt;
    cs->litr1c += cf->m_frootc_transfer_to_litr1c * dt;
    cs->frootc_transfer -= cf->m_frootc_transfer_to_litr1c * dt;
    cs->litr1c += cf->m_livestemc_transfer_to_litr1c * dt;
    cs->livestemc_transfer -= cf->m_livestemc_transfer_to_litr1c * dt;
    cs->litr1c += cf->m_deadstemc_transfer_to_litr1c * dt;
    cs->deadstemc_transfer -= cf->m_deadstemc_transfer_to_litr1c * dt;
    cs->litr1c += cf->m_livecrootc_transfer_to_litr1c * dt;
    cs->livecrootc_transfer -= cf->m_livecrootc_transfer_to_litr1c * dt;
    cs->litr1c += cf->m_deadcrootc_transfer_to_litr1c * dt;
    cs->deadcrootc_transfer -= cf->m_deadcrootc_transfer_to_litr1c * dt;
    cs->litr1c += cf->m_gresp_storage_to_litr1c * dt;
    cs->gresp_storage -= cf->m_gresp_storage_to_litr1c * dt;
    cs->litr1c += cf->m_gresp_transfer_to_litr1c * dt;
    cs->gresp_transfer -= cf->m_gresp_transfer_to_litr1c * dt;
    if (epc->woody)
    {
        /* Stem wood mortality */
        cs->cwdc += cf->m_livestemc_to_cwdc * dt;
        cs->livestemc -= cf->m_livestemc_to_cwdc * dt;
        cs->cwdc += cf->m_deadstemc_to_cwdc * dt;
        cs->deadstemc -= cf->m_deadstemc_to_cwdc * dt;
        /* STEP 1e Coarse root wood mortality */
        cs->cwdc += cf->m_livecrootc_to_cwdc * dt;
        cs->livecrootc -= cf->m_livecrootc_to_cwdc * dt;
        cs->cwdc += cf->m_deadcrootc_to_cwdc * dt;
        cs->deadcrootc -= cf->m_deadcrootc_to_cwdc * dt;
    }
    /* NITROGEN mortality state variable update */
    /*    Leaf mortality */
    ns->litr1n += nf->m_leafn_to_litr1n * dt;
    ns->leafn -= nf->m_leafn_to_litr1n * dt;
    ns->litr2n += nf->m_leafn_to_litr2n * dt;
    ns->leafn -= nf->m_leafn_to_litr2n * dt;
    ns->litr3n += nf->m_leafn_to_litr3n * dt;
    ns->leafn -= nf->m_leafn_to_litr3n * dt;
    ns->litr4n += nf->m_leafn_to_litr4n * dt;
    ns->leafn -= nf->m_leafn_to_litr4n * dt;
    /*    Fine root mortality */
    ns->litr1n += nf->m_frootn_to_litr1n * dt;
    ns->frootn -= nf->m_frootn_to_litr1n * dt;
    ns->litr2n += nf->m_frootn_to_litr2n * dt;
    ns->frootn -= nf->m_frootn_to_litr2n * dt;
    ns->litr3n += nf->m_frootn_to_litr3n * dt;
    ns->frootn -= nf->m_frootn_to_litr3n * dt;
    ns->litr4n += nf->m_frootn_to_litr4n * dt;
    ns->frootn -= nf->m_frootn_to_litr4n * dt;
    /*    Storage, transfer, excess, and npool mortality */
    ns->litr1n += nf->m_leafn_storage_to_litr1n * dt;
    ns->leafn_storage -= nf->m_leafn_storage_to_litr1n * dt;
    ns->litr1n += nf->m_frootn_storage_to_litr1n * dt;
    ns->frootn_storage -= nf->m_frootn_storage_to_litr1n * dt;
    ns->litr1n += nf->m_livestemn_storage_to_litr1n * dt;
    ns->livestemn_storage -= nf->m_livestemn_storage_to_litr1n * dt;
    ns->litr1n += nf->m_deadstemn_storage_to_litr1n * dt;
    ns->deadstemn_storage -= nf->m_deadstemn_storage_to_litr1n * dt;
    ns->litr1n += nf->m_livecrootn_storage_to_litr1n * dt;
    ns->livecrootn_storage -= nf->m_livecrootn_storage_to_litr1n * dt;
    ns->litr1n += nf->m_deadcrootn_storage_to_litr1n * dt;
    ns->deadcrootn_storage -= nf->m_deadcrootn_storage_to_litr1n * dt;
    ns->litr1n += nf->m_leafn_transfer_to_litr1n * dt;
    ns->leafn_transfer -= nf->m_leafn_transfer_to_litr1n * dt;
    ns->litr1n += nf->m_frootn_transfer_to_litr1n * dt;
    ns->frootn_transfer -= nf->m_frootn_transfer_to_litr1n * dt;
    ns->litr1n += nf->m_livestemn_transfer_to_litr1n * dt;
    ns->livestemn_transfer -= nf->m_livestemn_transfer_to_litr1n * dt;
    ns->litr1n += nf->m_deadstemn_transfer_to_litr1n * dt;
    ns->deadstemn_transfer -= nf->m_deadstemn_transfer_to_litr1n * dt;
    ns->litr1n += nf->m_livecrootn_transfer_to_litr1n * dt;
    ns->livecrootn_transfer -= nf->m_livecrootn_transfer_to_litr1n * dt;
    ns->litr1n += nf->m_deadcrootn_transfer_to_litr1n * dt;
    ns->deadcrootn_transfer -= nf->m_deadcrootn_transfer_to_litr1n * dt;
    ns->litr1n += nf->m_retransn_to_litr1n * dt;
    ns->retransn -= nf->m_retransn_to_litr1n * dt;
    if (epc->woody)
    {
        /*    Stem wood mortality */
        ns->litr1n += nf->m_livestemn_to_litr1n * dt;
        ns->livestemn -= nf->m_livestemn_to_litr1n * dt;
        ns->cwdn += nf->m_livestemn_to_cwdn * dt;
        ns->livestemn -= nf->m_livestemn_to_cwdn * dt;
        ns->cwdn += nf->m_deadstemn_to_cwdn * dt;
        ns->deadstemn -= nf->m_deadstemn_to_cwdn * dt;
        ns->litr1n += nf->m_livecrootn_to_litr1n * dt;
        ns->livecrootn -= nf->m_livecrootn_to_litr1n * dt;
        ns->cwdn += nf->m_livecrootn_to_cwdn * dt;
        ns->livecrootn -= nf->m_livecrootn_to_cwdn * dt;
        ns->cwdn += nf->m_deadcrootn_to_cwdn * dt;
        ns->deadcrootn -= nf->m_deadcrootn_to_cwdn * dt;
    }

    /*
     *Fire mortality: these fluxes all enter atmospheric sinks  *
     */
    fire = epc->daily_fire_turnover;

    /* daily carbon fluxes due to fireality */
    /* fireality fluxes out of leaf and fine root pools */
    cf->m_leafc_to_fire = fire * cs->leafc;
    cf->m_frootc_to_fire = fire * cs->frootc;

    /* fireality fluxes out of storage and transfer pools */
    cf->m_leafc_storage_to_fire = fire * cs->leafc_storage;
    cf->m_frootc_storage_to_fire = fire * cs->frootc_storage;
    cf->m_livestemc_storage_to_fire = fire * cs->livestemc_storage;
    cf->m_deadstemc_storage_to_fire = fire * cs->deadstemc_storage;
    cf->m_livecrootc_storage_to_fire = fire * cs->livecrootc_storage;
    cf->m_deadcrootc_storage_to_fire = fire * cs->deadcrootc_storage;
    cf->m_leafc_transfer_to_fire = fire * cs->leafc_transfer;
    cf->m_frootc_transfer_to_fire = fire * cs->frootc_transfer;
    cf->m_livestemc_transfer_to_fire = fire * cs->livestemc_transfer;
    cf->m_deadstemc_transfer_to_fire = fire * cs->deadstemc_transfer;
    cf->m_livecrootc_transfer_to_fire = fire * cs->livecrootc_transfer;
    cf->m_deadcrootc_transfer_to_fire = fire * cs->deadcrootc_transfer;
    cf->m_gresp_storage_to_fire = fire * cs->gresp_storage;
    cf->m_gresp_transfer_to_fire = fire * cs->gresp_transfer;

    /* TREE-specific carbon fluxes */
    if (epc->woody)
    {
        cf->m_livestemc_to_fire = fire * cs->livestemc;
        cf->m_deadstemc_to_fire = DSCP * fire * cs->deadstemc;
        cf->m_livecrootc_to_fire = fire * cs->livecrootc;
        cf->m_deadcrootc_to_fire = DSCP * fire * cs->deadcrootc;
    }

    /* litter and CWD fire fluxes */
    cf->m_litr1c_to_fire = fire * cs->litr1c;
    cf->m_litr2c_to_fire = fire * cs->litr2c;
    cf->m_litr3c_to_fire = fire * cs->litr3c;
    cf->m_litr4c_to_fire = fire * cs->litr4c;
    cf->m_cwdc_to_fire = 0.30 * fire * cs->cwdc;


    /* daily nitrogen fluxes due to fireality */
    /* fireality fluxes out of leaf and fine root pools */
    nf->m_leafn_to_fire = cf->m_leafc_to_fire / epc->leaf_cn;
    nf->m_frootn_to_fire = cf->m_frootc_to_fire / epc->froot_cn;

    /* fireality fluxes out of storage and transfer */
    nf->m_leafn_storage_to_fire = fire * ns->leafn_storage;
    nf->m_frootn_storage_to_fire = fire * ns->frootn_storage;
    nf->m_livestemn_storage_to_fire = fire * ns->livestemn_storage;
    nf->m_deadstemn_storage_to_fire = fire * ns->deadstemn_storage;
    nf->m_livecrootn_storage_to_fire = fire * ns->livecrootn_storage;
    nf->m_deadcrootn_storage_to_fire = fire * ns->deadcrootn_storage;
    nf->m_leafn_transfer_to_fire = fire * ns->leafn_transfer;
    nf->m_frootn_transfer_to_fire = fire * ns->frootn_transfer;
    nf->m_livestemn_transfer_to_fire = fire * ns->livestemn_transfer;
    nf->m_deadstemn_transfer_to_fire = fire * ns->deadstemn_transfer;
    nf->m_livecrootn_transfer_to_fire = fire * ns->livecrootn_transfer;
    nf->m_deadcrootn_transfer_to_fire = fire * ns->deadcrootn_transfer;
    nf->m_retransn_to_fire = fire * ns->retransn;

    /* woody-specific nitrogen fluxes */
    if (epc->woody)
    {
        nf->m_livestemn_to_fire = fire * ns->livestemn;
        nf->m_deadstemn_to_fire = DSCP * fire * ns->deadstemn;
        nf->m_livecrootn_to_fire = fire * ns->livecrootn;
        nf->m_deadcrootn_to_fire = DSCP * fire * ns->deadcrootn;
    }

    /* litter and CWD fire fluxes */
    nf->m_litr1n_to_fire = fire * ns->litr1n;
    nf->m_litr2n_to_fire = fire * ns->litr2n;
    nf->m_litr3n_to_fire = fire * ns->litr3n;
    nf->m_litr4n_to_fire = fire * ns->litr4n;
    nf->m_cwdn_to_fire = 0.30 * fire * ns->cwdn;

    /* update state variables for fire fluxes */
    /* this is the only place other than daily_state_update() routines where
     * state variables get changed.  Mortality is taken care of last and
     * given special treatment for state update so that it doesn't interfere
     * with the other fluxes that are based on proportions of state variables,
     * especially the phenological fluxes */
    /* CARBON fireality state variable update */
    /*   Leaf fireality */
    cs->fire_snk += cf->m_leafc_to_fire * dt;
    cs->leafc -= cf->m_leafc_to_fire * dt;
    /*   Fine root fireality */
    cs->fire_snk += cf->m_frootc_to_fire * dt;
    cs->frootc -= cf->m_frootc_to_fire * dt;
    /*   Storage and transfer fireality */
    cs->fire_snk += cf->m_leafc_storage_to_fire * dt;
    cs->leafc_storage -= cf->m_leafc_storage_to_fire * dt;
    cs->fire_snk += cf->m_frootc_storage_to_fire * dt;
    cs->frootc_storage -= cf->m_frootc_storage_to_fire * dt;
    cs->fire_snk += cf->m_livestemc_storage_to_fire * dt;
    cs->livestemc_storage -= cf->m_livestemc_storage_to_fire * dt;
    cs->fire_snk += cf->m_deadstemc_storage_to_fire * dt;
    cs->deadstemc_storage -= cf->m_deadstemc_storage_to_fire * dt;
    cs->fire_snk += cf->m_livecrootc_storage_to_fire * dt;
    cs->livecrootc_storage -= cf->m_livecrootc_storage_to_fire * dt;
    cs->fire_snk += cf->m_deadcrootc_storage_to_fire * dt;
    cs->deadcrootc_storage -= cf->m_deadcrootc_storage_to_fire * dt;
    cs->fire_snk += cf->m_leafc_transfer_to_fire * dt;
    cs->leafc_transfer -= cf->m_leafc_transfer_to_fire * dt;
    cs->fire_snk += cf->m_frootc_transfer_to_fire * dt;
    cs->frootc_transfer -= cf->m_frootc_transfer_to_fire * dt;
    cs->fire_snk += cf->m_livestemc_transfer_to_fire * dt;
    cs->livestemc_transfer -= cf->m_livestemc_transfer_to_fire * dt;
    cs->fire_snk += cf->m_deadstemc_transfer_to_fire * dt;
    cs->deadstemc_transfer -= cf->m_deadstemc_transfer_to_fire * dt;
    cs->fire_snk += cf->m_livecrootc_transfer_to_fire * dt;
    cs->livecrootc_transfer -= cf->m_livecrootc_transfer_to_fire * dt;
    cs->fire_snk += cf->m_deadcrootc_transfer_to_fire * dt;
    cs->deadcrootc_transfer -= cf->m_deadcrootc_transfer_to_fire * dt;
    cs->fire_snk += cf->m_gresp_storage_to_fire * dt;
    cs->gresp_storage -= cf->m_gresp_storage_to_fire * dt;
    cs->fire_snk += cf->m_gresp_transfer_to_fire * dt;
    cs->gresp_transfer -= cf->m_gresp_transfer_to_fire * dt;
    if (epc->woody)
    {
        /* Stem wood mortality */
        cs->fire_snk += cf->m_livestemc_to_fire * dt;
        cs->livestemc -= cf->m_livestemc_to_fire * dt;
        cs->fire_snk += cf->m_deadstemc_to_fire * dt;
        cs->deadstemc -= cf->m_deadstemc_to_fire * dt;
        /* Coarse root wood mortality */
        cs->fire_snk += cf->m_livecrootc_to_fire * dt;
        cs->livecrootc -= cf->m_livecrootc_to_fire * dt;
        cs->fire_snk += cf->m_deadcrootc_to_fire * dt;
        cs->deadcrootc -= cf->m_deadcrootc_to_fire * dt;
    }
    /* litter and CWD carbon state updates */
    cs->fire_snk += cf->m_litr1c_to_fire * dt;
    cs->litr1c -= cf->m_litr1c_to_fire * dt;
    cs->fire_snk += cf->m_litr2c_to_fire * dt;
    cs->litr2c -= cf->m_litr2c_to_fire * dt;
    cs->fire_snk += cf->m_litr3c_to_fire * dt;
    cs->litr3c -= cf->m_litr3c_to_fire * dt;
    cs->fire_snk += cf->m_litr4c_to_fire * dt;
    cs->litr4c -= cf->m_litr4c_to_fire * dt;
    cs->fire_snk += cf->m_cwdc_to_fire * dt;
    cs->cwdc -= cf->m_cwdc_to_fire * dt;

    /* NITROGEN mortality state variable update */
    /* Leaf mortality */
    ns->fire_snk += nf->m_leafn_to_fire * dt;
    ns->leafn -= nf->m_leafn_to_fire * dt;
    /* Fine root mortality */
    ns->fire_snk += nf->m_frootn_to_fire * dt;
    ns->frootn -= nf->m_frootn_to_fire * dt;
    /* Storage, transfer, excess, and npool mortality */
    ns->fire_snk += nf->m_leafn_storage_to_fire * dt;
    ns->leafn_storage -= nf->m_leafn_storage_to_fire * dt;
    ns->fire_snk += nf->m_frootn_storage_to_fire * dt;
    ns->frootn_storage -= nf->m_frootn_storage_to_fire * dt;
    ns->fire_snk += nf->m_livestemn_storage_to_fire * dt;
    ns->livestemn_storage -= nf->m_livestemn_storage_to_fire * dt;
    ns->fire_snk += nf->m_deadstemn_storage_to_fire * dt;
    ns->deadstemn_storage -= nf->m_deadstemn_storage_to_fire * dt;
    ns->fire_snk += nf->m_livecrootn_storage_to_fire * dt;
    ns->livecrootn_storage -= nf->m_livecrootn_storage_to_fire * dt;
    ns->fire_snk += nf->m_deadcrootn_storage_to_fire * dt;
    ns->deadcrootn_storage -= nf->m_deadcrootn_storage_to_fire * dt;
    ns->fire_snk += nf->m_leafn_transfer_to_fire * dt;
    ns->leafn_transfer -= nf->m_leafn_transfer_to_fire * dt;
    ns->fire_snk += nf->m_frootn_transfer_to_fire * dt;
    ns->frootn_transfer -= nf->m_frootn_transfer_to_fire * dt;
    ns->fire_snk += nf->m_livestemn_transfer_to_fire * dt;
    ns->livestemn_transfer -= nf->m_livestemn_transfer_to_fire * dt;
    ns->fire_snk += nf->m_deadstemn_transfer_to_fire * dt;
    ns->deadstemn_transfer -= nf->m_deadstemn_transfer_to_fire * dt;
    ns->fire_snk += nf->m_livecrootn_transfer_to_fire * dt;
    ns->livecrootn_transfer -= nf->m_livecrootn_transfer_to_fire * dt;
    ns->fire_snk += nf->m_deadcrootn_transfer_to_fire * dt;
    ns->deadcrootn_transfer -= nf->m_deadcrootn_transfer_to_fire * dt;
    ns->fire_snk += nf->m_retransn_to_fire * dt;
    ns->retransn -= nf->m_retransn_to_fire * dt;
    if (epc->woody)
    {
        /* Stem wood mortality */
        ns->fire_snk += nf->m_livestemn_to_fire * dt;
        ns->livestemn -= nf->m_livestemn_to_fire * dt;
        ns->fire_snk += nf->m_deadstemn_to_fire * dt;
        ns->deadstemn -= nf->m_deadstemn_to_fire * dt;
        ns->fire_snk += nf->m_livecrootn_to_fire * dt;
        ns->livecrootn -= nf->m_livecrootn_to_fire * dt;
        ns->fire_snk += nf->m_deadcrootn_to_fire * dt;
        ns->deadcrootn -= nf->m_deadcrootn_to_fire * dt;
    }
    /* litter and CWD nitrogen state updates */
    ns->fire_snk += nf->m_litr1n_to_fire * dt;
    ns->litr1n -= nf->m_litr1n_to_fire * dt;
    ns->fire_snk += nf->m_litr2n_to_fire * dt;
    ns->litr2n -= nf->m_litr2n_to_fire * dt;
    ns->fire_snk += nf->m_litr3n_to_fire * dt;
    ns->litr3n -= nf->m_litr3n_to_fire * dt;
    ns->fire_snk += nf->m_litr4n_to_fire * dt;
    ns->litr4n -= nf->m_litr4n_to_fire * dt;
    ns->fire_snk += nf->m_cwdn_to_fire * dt;
    ns->cwdn -= nf->m_cwdn_to_fire * dt;
}
