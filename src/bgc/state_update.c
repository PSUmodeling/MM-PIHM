
/*
 * state_update.c
 * Resolve the fluxes in bgc() daily loop to update state variables
 *
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 * Biome-BGC version 4.2 (final release)
 * See copyright.txt for Copyright information
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 */

#include "pihm.h"

void CarbonStateUpdate (cflux_struct *cf, cstate_struct *cs, int alloc,
    int woody, int evergreen, double dt)
{
    /* daily update of the carbon state variables */

    /* C state variables are updated below in the order of the relevant fluxes
     * in the daily model loop */

    /* NOTE: Mortality fluxes are all accounted for in a separate routine,
     * which is to be called after this routine.  This is a special case where
     * the updating of state variables is order-sensitive, since otherwise the
     * complications of possibly having mortality fluxes drive the pools
     * negative would create big, unnecessary headaches. */

    /* Phenology fluxes */
    /* leaf and fine root transfer growth */
    cs->leafc += cf->leafc_transfer_to_leafc * dt;
    cs->leafc_transfer -= cf->leafc_transfer_to_leafc * dt;
    cs->frootc += cf->frootc_transfer_to_frootc * dt;
    cs->frootc_transfer -= cf->frootc_transfer_to_frootc * dt;
    if (woody)
    {
        /* Stem and coarse root transfer growth */
        cs->livestemc += cf->livestemc_transfer_to_livestemc * dt;
        cs->livestemc_transfer -= cf->livestemc_transfer_to_livestemc * dt;
        cs->deadstemc += cf->deadstemc_transfer_to_deadstemc * dt;
        cs->deadstemc_transfer -= cf->deadstemc_transfer_to_deadstemc * dt;
        cs->livecrootc += cf->livecrootc_transfer_to_livecrootc * dt;
        cs->livecrootc_transfer -= cf->livecrootc_transfer_to_livecrootc * dt;
        cs->deadcrootc += cf->deadcrootc_transfer_to_deadcrootc * dt;
        cs->deadcrootc_transfer -= cf->deadcrootc_transfer_to_deadcrootc * dt;
    }
    /* Leaf and fine root litterfall */
    cs->litr1c += cf->leafc_to_litr1c * dt;
    cs->leafc -= cf->leafc_to_litr1c * dt;
    cs->litr2c += cf->leafc_to_litr2c * dt;
    cs->leafc -= cf->leafc_to_litr2c * dt;
    cs->litr3c += cf->leafc_to_litr3c * dt;
    cs->leafc -= cf->leafc_to_litr3c * dt;
    cs->litr4c += cf->leafc_to_litr4c * dt;
    cs->leafc -= cf->leafc_to_litr4c * dt;
    cs->litr1c += cf->frootc_to_litr1c * dt;
    cs->frootc -= cf->frootc_to_litr1c * dt;
    cs->litr2c += cf->frootc_to_litr2c * dt;
    cs->frootc -= cf->frootc_to_litr2c * dt;
    cs->litr3c += cf->frootc_to_litr3c * dt;
    cs->frootc -= cf->frootc_to_litr3c * dt;
    cs->litr4c += cf->frootc_to_litr4c * dt;
    cs->frootc -= cf->frootc_to_litr4c * dt;
    /* livewood turnover fluxes */
    if (woody)
    {
        cs->deadstemc += cf->livestemc_to_deadstemc * dt;
        cs->livestemc -= cf->livestemc_to_deadstemc * dt;
        cs->deadcrootc += cf->livecrootc_to_deadcrootc * dt;
        cs->livecrootc -= cf->livecrootc_to_deadcrootc * dt;
    }

    /* Maintenance respiration fluxes */
    cs->leaf_mr_snk += cf->leaf_mr * dt;
    cs->cpool -= cf->leaf_mr * dt;
    cs->froot_mr_snk += cf->froot_mr * dt;
    cs->cpool -= cf->froot_mr * dt;
    if (woody)
    {
        cs->livestem_mr_snk += cf->livestem_mr * dt;
        cs->cpool -= cf->livestem_mr * dt;
        cs->livecroot_mr_snk += cf->livecroot_mr * dt;
        cs->cpool -= cf->livecroot_mr * dt;
    }

    /* Photosynthesis fluxes */
    cs->cpool += cf->psnsun_to_cpool * dt;
    cs->psnsun_src += cf->psnsun_to_cpool * dt;
    cs->cpool += cf->psnshade_to_cpool * dt;
    cs->psnshade_src += cf->psnshade_to_cpool * dt;

    /* Litter decomposition fluxes */
    /* Fluxes out of coarse woody debris into litter pools */
    cs->litr2c += cf->cwdc_to_litr2c * dt;
    cs->cwdc -= cf->cwdc_to_litr2c * dt;
    cs->litr3c += cf->cwdc_to_litr3c * dt;
    cs->cwdc -= cf->cwdc_to_litr3c * dt;
    cs->litr4c += cf->cwdc_to_litr4c * dt;
    cs->cwdc -= cf->cwdc_to_litr4c * dt;
    /* Fluxes out of labile litter pool */
    cs->litr1_hr_snk += cf->litr1_hr * dt;
    cs->litr1c -= cf->litr1_hr * dt;
    cs->soil1c += cf->litr1c_to_soil1c * dt;
    cs->litr1c -= cf->litr1c_to_soil1c * dt;
    /* Fluxes out of cellulose litter pool */
    cs->litr2_hr_snk += cf->litr2_hr * dt;
    cs->litr2c -= cf->litr2_hr * dt;
    cs->soil2c += cf->litr2c_to_soil2c * dt;
    cs->litr2c -= cf->litr2c_to_soil2c * dt;
    /* Fluxes from shielded to unshielded cellulose pools */
    cs->litr2c += cf->litr3c_to_litr2c * dt;
    cs->litr3c -= cf->litr3c_to_litr2c * dt;
    /* Fluxes out of lignin litter pool */
    cs->litr4_hr_snk += cf->litr4_hr * dt;
    cs->litr4c -= cf->litr4_hr * dt;
    cs->soil3c += cf->litr4c_to_soil3c * dt;
    cs->litr4c -= cf->litr4c_to_soil3c * dt;
    /* Fluxes out of fast soil pool */
    cs->soil1_hr_snk += cf->soil1_hr * dt;
    cs->soil1c -= cf->soil1_hr * dt;
    cs->soil2c += cf->soil1c_to_soil2c * dt;
    cs->soil1c -= cf->soil1c_to_soil2c * dt;
    /* Fluxes out of medium soil pool */
    cs->soil2_hr_snk += cf->soil2_hr * dt;
    cs->soil2c -= cf->soil2_hr * dt;
    cs->soil3c += cf->soil2c_to_soil3c * dt;
    cs->soil2c -= cf->soil2c_to_soil3c * dt;
    /* Fluxes out of slow soil pool */
    cs->soil3_hr_snk += cf->soil3_hr * dt;
    cs->soil3c -= cf->soil3_hr * dt;
    cs->soil4c += cf->soil3c_to_soil4c * dt;
    cs->soil3c -= cf->soil3c_to_soil4c * dt;
    /* Fluxes out of recalcitrant SOM pool */
    cs->soil4_hr_snk += cf->soil4_hr * dt;
    cs->soil4c -= cf->soil4_hr * dt;

    /* Daily allocation fluxes */
    /* daily leaf allocation fluxes */
    cs->leafc += cf->cpool_to_leafc * dt;
    cs->cpool -= cf->cpool_to_leafc * dt;
    cs->leafc_storage += cf->cpool_to_leafc_storage * dt;
    cs->cpool -= cf->cpool_to_leafc_storage * dt;
    /* Daily fine root allocation fluxes */
    cs->frootc += cf->cpool_to_frootc * dt;
    cs->cpool -= cf->cpool_to_frootc * dt;
    cs->frootc_storage += cf->cpool_to_frootc_storage * dt;
    cs->cpool -= cf->cpool_to_frootc_storage * dt;
    if (woody)
    {
        /* Daily live stem wood allocation fluxes */
        cs->livestemc += cf->cpool_to_livestemc * dt;
        cs->cpool -= cf->cpool_to_livestemc * dt;
        cs->livestemc_storage += cf->cpool_to_livestemc_storage * dt;
        cs->cpool -= cf->cpool_to_livestemc_storage * dt;
        /* Daily dead stem wood allocation fluxes */
        cs->deadstemc += cf->cpool_to_deadstemc * dt;
        cs->cpool -= cf->cpool_to_deadstemc * dt;
        cs->deadstemc_storage += cf->cpool_to_deadstemc_storage * dt;
        cs->cpool -= cf->cpool_to_deadstemc_storage * dt;
        /* Daily live coarse root wood allocation fluxes */
        cs->livecrootc += cf->cpool_to_livecrootc * dt;
        cs->cpool -= cf->cpool_to_livecrootc * dt;
        cs->livecrootc_storage += cf->cpool_to_livecrootc_storage * dt;
        cs->cpool -= cf->cpool_to_livecrootc_storage * dt;
        /* Daily dead coarse root wood allocation fluxes */
        cs->deadcrootc += cf->cpool_to_deadcrootc * dt;
        cs->cpool -= cf->cpool_to_deadcrootc * dt;
        cs->deadcrootc_storage += cf->cpool_to_deadcrootc_storage * dt;
        cs->cpool -= cf->cpool_to_deadcrootc_storage * dt;
    }
    /* Daily allocation for transfer growth respiration */
    cs->gresp_storage += cf->cpool_to_gresp_storage * dt;
    cs->cpool -= cf->cpool_to_gresp_storage * dt;

    /* Daily growth respiration fluxes */
    /* Leaf growth respiration */
    cs->leaf_gr_snk += cf->cpool_leaf_gr * dt;
    cs->cpool -= cf->cpool_leaf_gr * dt;
    cs->leaf_gr_snk += cf->cpool_leaf_storage_gr * dt;
    cs->cpool -= cf->cpool_leaf_storage_gr * dt;
    cs->leaf_gr_snk += cf->transfer_leaf_gr * dt;
    cs->gresp_transfer -= cf->transfer_leaf_gr * dt;
    /* Fine root growth respiration */
    cs->froot_gr_snk += cf->cpool_froot_gr * dt;
    cs->cpool -= cf->cpool_froot_gr * dt;
    cs->froot_gr_snk += cf->cpool_froot_storage_gr * dt;
    cs->cpool -= cf->cpool_froot_storage_gr * dt;
    cs->froot_gr_snk += cf->transfer_froot_gr * dt;
    cs->gresp_transfer -= cf->transfer_froot_gr * dt;
    if (woody)
    {
        /* Live stem growth respiration */
        cs->livestem_gr_snk += cf->cpool_livestem_gr * dt;
        cs->cpool -= cf->cpool_livestem_gr * dt;
        cs->livestem_gr_snk += cf->cpool_livestem_storage_gr * dt;
        cs->cpool -= cf->cpool_livestem_storage_gr * dt;
        cs->livestem_gr_snk += cf->transfer_livestem_gr * dt;
        cs->gresp_transfer -= cf->transfer_livestem_gr * dt;
        /* Dead stem growth respiration */
        cs->deadstem_gr_snk += cf->cpool_deadstem_gr * dt;
        cs->cpool -= cf->cpool_deadstem_gr * dt;
        cs->deadstem_gr_snk += cf->cpool_deadstem_storage_gr * dt;
        cs->cpool -= cf->cpool_deadstem_storage_gr * dt;
        cs->deadstem_gr_snk += cf->transfer_deadstem_gr * dt;
        cs->gresp_transfer -= cf->transfer_deadstem_gr * dt;
        /* Live coarse root growth respiration */
        cs->livecroot_gr_snk += cf->cpool_livecroot_gr * dt;
        cs->cpool -= cf->cpool_livecroot_gr * dt;
        cs->livecroot_gr_snk += cf->cpool_livecroot_storage_gr * dt;
        cs->cpool -= cf->cpool_livecroot_storage_gr * dt;
        cs->livecroot_gr_snk += cf->transfer_livecroot_gr * dt;
        cs->gresp_transfer -= cf->transfer_livecroot_gr * dt;
        /* Dead coarse root growth respiration */
        cs->deadcroot_gr_snk += cf->cpool_deadcroot_gr * dt;
        cs->cpool -= cf->cpool_deadcroot_gr * dt;
        cs->deadcroot_gr_snk += cf->cpool_deadcroot_storage_gr * dt;
        cs->cpool -= cf->cpool_deadcroot_storage_gr * dt;
        cs->deadcroot_gr_snk += cf->transfer_deadcroot_gr * dt;
        cs->gresp_transfer -= cf->transfer_deadcroot_gr * dt;
    }

    /* Annual allocation fluxes, one day per year */
    if (alloc)
    {
        /* Move storage material into transfer compartments on the annual
         * allocation day. This is a special case, where a flux is assessed
         * in the state_variable update routine.  This is required to have the
         * allocation of excess C and N show up as new growth in the next
         * growing season, instead of two growing seasons from now. */
        cf->leafc_storage_to_leafc_transfer = cs->leafc_storage / dt;
        cf->frootc_storage_to_frootc_transfer = cs->frootc_storage / dt;
        cf->gresp_storage_to_gresp_transfer = cs->gresp_storage / dt;
        if (woody)
        {
            cf->livestemc_storage_to_livestemc_transfer =
                cs->livestemc_storage / dt;
            cf->deadstemc_storage_to_deadstemc_transfer =
                cs->deadstemc_storage / dt;
            cf->livecrootc_storage_to_livecrootc_transfer =
                cs->livecrootc_storage / dt;
            cf->deadcrootc_storage_to_deadcrootc_transfer =
                cs->deadcrootc_storage / dt;
        }
        /* update states variables */
        cs->leafc_transfer += cf->leafc_storage_to_leafc_transfer * dt;
        cs->leafc_storage -= cf->leafc_storage_to_leafc_transfer * dt;
        cs->frootc_transfer += cf->frootc_storage_to_frootc_transfer * dt;
        cs->frootc_storage -= cf->frootc_storage_to_frootc_transfer * dt;
        cs->gresp_transfer += cf->gresp_storage_to_gresp_transfer * dt;
        cs->gresp_storage -= cf->gresp_storage_to_gresp_transfer * dt;
        if (woody)
        {
            cs->livestemc_transfer +=
                cf->livestemc_storage_to_livestemc_transfer * dt;
            cs->livestemc_storage -=
                cf->livestemc_storage_to_livestemc_transfer * dt;
            cs->deadstemc_transfer +=
                cf->deadstemc_storage_to_deadstemc_transfer * dt;
            cs->deadstemc_storage -=
                cf->deadstemc_storage_to_deadstemc_transfer * dt;
            cs->livecrootc_transfer +=
                cf->livecrootc_storage_to_livecrootc_transfer * dt;
            cs->livecrootc_storage -=
                cf->livecrootc_storage_to_livecrootc_transfer * dt;
            cs->deadcrootc_transfer +=
                cf->deadcrootc_storage_to_deadcrootc_transfer * dt;
            cs->deadcrootc_storage -=
                cf->deadcrootc_storage_to_deadcrootc_transfer * dt;
        }

        /* for deciduous system, force leafc and frootc to exactly 0.0 on the
         * last day */
        if (!evergreen)
        {
            if (cs->leafc < 1e-10)
            {
                cs->leafc = 0.0;
            }
            if (cs->frootc < 1e-10)
            {
                cs->frootc = 0.0;
            }
        }
    }                           /* end if allocation day */
}

void NitrogenStateUpdate (nflux_struct *nf, nstate_struct *ns, int alloc,
    int woody, int evergreen, double dt)
{
    /* N state variables are updated below in the order of the relevant fluxes
     * in the daily model loop */

    /* NOTE: Mortality fluxes are all accounted for in a separate routine,
     * which is to be called after this routine.  This is a special case where
     * the updating of state variables is order-sensitive, since otherwise the
     * complications of possibly having mortality fluxes drive the pools
     * negative would create big, unnecessary headaches. */

    /* Phenology fluxes */
    /* Leaf and fine root transfer growth */
    ns->leafn += nf->leafn_transfer_to_leafn * dt;
    ns->leafn_transfer -= nf->leafn_transfer_to_leafn * dt;
    ns->frootn += nf->frootn_transfer_to_frootn * dt;
    ns->frootn_transfer -= nf->frootn_transfer_to_frootn * dt;
    if (woody)
    {
        /* Stem and coarse root transfer growth */
        ns->livestemn += nf->livestemn_transfer_to_livestemn * dt;
        ns->livestemn_transfer -= nf->livestemn_transfer_to_livestemn * dt;
        ns->deadstemn += nf->deadstemn_transfer_to_deadstemn * dt;
        ns->deadstemn_transfer -= nf->deadstemn_transfer_to_deadstemn * dt;
        ns->livecrootn += nf->livecrootn_transfer_to_livecrootn * dt;
        ns->livecrootn_transfer -= nf->livecrootn_transfer_to_livecrootn * dt;
        ns->deadcrootn += nf->deadcrootn_transfer_to_deadcrootn * dt;
        ns->deadcrootn_transfer -= nf->deadcrootn_transfer_to_deadcrootn * dt;
    }
    /* Leaf and fine root litterfall */
    ns->litr1n += nf->leafn_to_litr1n * dt;
    ns->leafn -= nf->leafn_to_litr1n * dt;
    ns->litr2n += nf->leafn_to_litr2n * dt;
    ns->leafn -= nf->leafn_to_litr2n * dt;
    ns->litr3n += nf->leafn_to_litr3n * dt;
    ns->leafn -= nf->leafn_to_litr3n * dt;
    ns->litr4n += nf->leafn_to_litr4n * dt;
    ns->leafn -= nf->leafn_to_litr4n * dt;
    ns->retransn += nf->leafn_to_retransn * dt;      /* N retranslocation */
    ns->leafn -= nf->leafn_to_retransn * dt;
    ns->litr1n += nf->frootn_to_litr1n * dt;
    ns->frootn -= nf->frootn_to_litr1n * dt;
    ns->litr2n += nf->frootn_to_litr2n * dt;
    ns->frootn -= nf->frootn_to_litr2n * dt;
    ns->litr3n += nf->frootn_to_litr3n * dt;
    ns->frootn -= nf->frootn_to_litr3n * dt;
    ns->litr4n += nf->frootn_to_litr4n * dt;
    ns->frootn -= nf->frootn_to_litr4n * dt;
    /* live wood turnover to dead wood */
    ns->deadstemn += nf->livestemn_to_deadstemn * dt;
    ns->livestemn -= nf->livestemn_to_deadstemn * dt;
    ns->retransn += nf->livestemn_to_retransn * dt;  /* N retranslocation */
    ns->livestemn -= nf->livestemn_to_retransn * dt;
    ns->deadcrootn += nf->livecrootn_to_deadcrootn * dt;
    ns->livecrootn -= nf->livecrootn_to_deadcrootn * dt;
    ns->retransn += nf->livecrootn_to_retransn * dt; /* N retranslocation */
    ns->livecrootn -= nf->livecrootn_to_retransn * dt;

    /* Nitrogen deposition */
    ns->sminn += nf->ndep_to_sminn * dt;
    ns->ndep_src += nf->ndep_to_sminn * dt;
    ns->sminn += nf->nfix_to_sminn * dt;
    ns->nfix_src += nf->nfix_to_sminn * dt;

    /* Litter and soil decomposition fluxes */
    /* Fluxes out of coarse woody debris into litter pools */
    ns->litr2n += nf->cwdn_to_litr2n * dt;
    ns->cwdn -= nf->cwdn_to_litr2n * dt;
    ns->litr3n += nf->cwdn_to_litr3n * dt;
    ns->cwdn -= nf->cwdn_to_litr3n * dt;
    ns->litr4n += nf->cwdn_to_litr4n * dt;
    ns->cwdn -= nf->cwdn_to_litr4n * dt;
    /* N fluxes for immobilization and mineralization */
    ns->soil1n += nf->litr1n_to_soil1n * dt;
    ns->litr1n -= nf->litr1n_to_soil1n * dt;
    if (nf->sminn_to_soil1n_l1 < 0.0)
    {
        nf->sminn_to_nvol_l1s1 =
            -DENITRIF_PROPORTION * nf->sminn_to_soil1n_l1;
    }
    else
    {
        nf->sminn_to_nvol_l1s1 = 0.0;
    }
    ns->soil1n += nf->sminn_to_soil1n_l1 * dt;
    ns->sminn -= nf->sminn_to_soil1n_l1 * dt;
    ns->nvol_snk += nf->sminn_to_nvol_l1s1 * dt;
    ns->sminn -= nf->sminn_to_nvol_l1s1 * dt;

    ns->soil2n += nf->litr2n_to_soil2n * dt;
    ns->litr2n -= nf->litr2n_to_soil2n * dt;
    if (nf->sminn_to_soil2n_l2 < 0.0)
    {
        nf->sminn_to_nvol_l2s2 =
            -DENITRIF_PROPORTION * nf->sminn_to_soil2n_l2;
    }
    else
    {
        nf->sminn_to_nvol_l2s2 = 0.0;
    }
    ns->soil2n += nf->sminn_to_soil2n_l2 * dt;
    ns->sminn -= nf->sminn_to_soil2n_l2 * dt;
    ns->nvol_snk += nf->sminn_to_nvol_l2s2 * dt;
    ns->sminn -= nf->sminn_to_nvol_l2s2 * dt;

    ns->litr2n += nf->litr3n_to_litr2n * dt;
    ns->litr3n -= nf->litr3n_to_litr2n * dt;

    ns->soil3n += nf->litr4n_to_soil3n * dt;
    ns->litr4n -= nf->litr4n_to_soil3n * dt;
    if (nf->sminn_to_soil3n_l4 < 0.0)
    {
        nf->sminn_to_nvol_l4s3 =
            -DENITRIF_PROPORTION * nf->sminn_to_soil3n_l4;
    }
    else
    {
        nf->sminn_to_nvol_l4s3 = 0.0;
    }
    ns->soil3n += nf->sminn_to_soil3n_l4 * dt;
    ns->sminn -= nf->sminn_to_soil3n_l4 * dt;
    ns->nvol_snk += nf->sminn_to_nvol_l4s3 * dt;
    ns->sminn -= nf->sminn_to_nvol_l4s3 * dt;

    ns->soil2n += nf->soil1n_to_soil2n * dt;
    ns->soil1n -= nf->soil1n_to_soil2n * dt;
    if (nf->sminn_to_soil2n_s1 < 0.0)
    {
        nf->sminn_to_nvol_s1s2 =
            -DENITRIF_PROPORTION * nf->sminn_to_soil2n_s1;
    }
    else
    {
        nf->sminn_to_nvol_s1s2 = 0.0;
    }
    ns->soil2n += nf->sminn_to_soil2n_s1 * dt;
    ns->sminn -= nf->sminn_to_soil2n_s1 * dt;
    ns->nvol_snk += nf->sminn_to_nvol_s1s2 * dt;
    ns->sminn -= nf->sminn_to_nvol_s1s2 * dt;

    ns->soil3n += nf->soil2n_to_soil3n * dt;
    ns->soil2n -= nf->soil2n_to_soil3n * dt;
    if (nf->sminn_to_soil3n_s2 < 0.0)
    {
        nf->sminn_to_nvol_s2s3 =
            -DENITRIF_PROPORTION * nf->sminn_to_soil3n_s2;
    }
    else
    {
        nf->sminn_to_nvol_s2s3 = 0.0;
    }
    ns->soil3n += nf->sminn_to_soil3n_s2 * dt;
    ns->sminn -= nf->sminn_to_soil3n_s2 * dt;
    ns->nvol_snk += nf->sminn_to_nvol_s2s3 * dt;
    ns->sminn -= nf->sminn_to_nvol_s2s3 * dt;

    ns->soil4n += nf->soil3n_to_soil4n * dt;
    ns->soil3n -= nf->soil3n_to_soil4n * dt;
    if (nf->sminn_to_soil4n_s3 < 0.0)
    {
        nf->sminn_to_nvol_s3s4 =
            -DENITRIF_PROPORTION * nf->sminn_to_soil4n_s3;
    }
    else
    {
        nf->sminn_to_nvol_s3s4 = 0.0;
    }
    ns->soil4n += nf->sminn_to_soil4n_s3 * dt;
    ns->sminn -= nf->sminn_to_soil4n_s3 * dt;
    ns->nvol_snk += nf->sminn_to_nvol_s3s4 * dt;
    ns->sminn -= nf->sminn_to_nvol_s3s4 * dt;

    nf->sminn_to_nvol_s4 = DENITRIF_PROPORTION * nf->soil4n_to_sminn;
    ns->sminn += nf->soil4n_to_sminn * dt;
    ns->soil4n -= nf->soil4n_to_sminn * dt;
    ns->nvol_snk += nf->sminn_to_nvol_s4 * dt;
    ns->sminn -= nf->sminn_to_nvol_s4 * dt;

    /* Bulk denitrification of soil mineral N */
    ns->sminn -= nf->sminn_to_denitrif * dt;
    ns->nvol_snk += nf->sminn_to_denitrif * dt;

    /* Plant allocation flux, from N retrans pool and soil mineral N pool */
    ns->npool += nf->retransn_to_npool * dt;
    ns->retransn -= nf->retransn_to_npool * dt;
    ns->npool += nf->sminn_to_npool * dt;
    ns->sminn -= nf->sminn_to_npool * dt;

    /* Daily allocation fluxes */
    /* Daily leaf allocation fluxes */
    ns->leafn += nf->npool_to_leafn * dt;
    ns->npool -= nf->npool_to_leafn * dt;
    ns->leafn_storage += nf->npool_to_leafn_storage * dt;
    ns->npool -= nf->npool_to_leafn_storage * dt;
    /* Daily fine root allocation fluxes */
    ns->frootn += nf->npool_to_frootn * dt;
    ns->npool -= nf->npool_to_frootn * dt;
    ns->frootn_storage += nf->npool_to_frootn_storage * dt;
    ns->npool -= nf->npool_to_frootn_storage * dt;
    if (woody)
    {
        /* Daily live stem allocation fluxes */
        ns->livestemn += nf->npool_to_livestemn * dt;
        ns->npool -= nf->npool_to_livestemn * dt;
        ns->livestemn_storage += nf->npool_to_livestemn_storage * dt;
        ns->npool -= nf->npool_to_livestemn_storage * dt;
        /* Daily dead stem allocation fluxes */
        ns->deadstemn += nf->npool_to_deadstemn * dt;
        ns->npool -= nf->npool_to_deadstemn * dt;
        ns->deadstemn_storage += nf->npool_to_deadstemn_storage * dt;
        ns->npool -= nf->npool_to_deadstemn_storage * dt;
        /* Daily live coarse root allocation fluxes */
        ns->livecrootn += nf->npool_to_livecrootn * dt;
        ns->npool -= nf->npool_to_livecrootn * dt;
        ns->livecrootn_storage += nf->npool_to_livecrootn_storage * dt;
        ns->npool -= nf->npool_to_livecrootn_storage * dt;
        /* Daily dead coarse root allocation fluxes */
        ns->deadcrootn += nf->npool_to_deadcrootn * dt;
        ns->npool -= nf->npool_to_deadcrootn * dt;
        ns->deadcrootn_storage += nf->npool_to_deadcrootn_storage * dt;
        ns->npool -= nf->npool_to_deadcrootn_storage * dt;
    }

    /* Annual allocation fluxes, one day per year */
    if (alloc)
    {
        /* Move storage material into transfer compartments on the annual
         * allocation day. This is a special case, where a flux is assessed in
         * the state_variable update routine.  This is required to have the
         * allocation of excess C and N show up as new growth in the next
         * growing season, instead of two growing seasons from now. */
        nf->leafn_storage_to_leafn_transfer = ns->leafn_storage / dt;
        nf->frootn_storage_to_frootn_transfer = ns->frootn_storage / dt;
        if (woody)
        {
            nf->livestemn_storage_to_livestemn_transfer =
                ns->livestemn_storage / dt;
            nf->deadstemn_storage_to_deadstemn_transfer =
                ns->deadstemn_storage / dt;
            nf->livecrootn_storage_to_livecrootn_transfer =
                ns->livecrootn_storage / dt;
            nf->deadcrootn_storage_to_deadcrootn_transfer =
                ns->deadcrootn_storage / dt;
        }
        /* update states variables */
        ns->leafn_transfer += nf->leafn_storage_to_leafn_transfer * dt;
        ns->leafn_storage -= nf->leafn_storage_to_leafn_transfer * dt;
        ns->frootn_transfer += nf->frootn_storage_to_frootn_transfer * dt;
        ns->frootn_storage -= nf->frootn_storage_to_frootn_transfer * dt;
        if (woody)
        {
            ns->livestemn_transfer +=
                nf->livestemn_storage_to_livestemn_transfer * dt;
            ns->livestemn_storage -=
                nf->livestemn_storage_to_livestemn_transfer * dt;
            ns->deadstemn_transfer +=
                nf->deadstemn_storage_to_deadstemn_transfer * dt;
            ns->deadstemn_storage -=
                nf->deadstemn_storage_to_deadstemn_transfer * dt;
            ns->livecrootn_transfer +=
                nf->livecrootn_storage_to_livecrootn_transfer * dt;
            ns->livecrootn_storage -=
                nf->livecrootn_storage_to_livecrootn_transfer * dt;
            ns->deadcrootn_transfer +=
                nf->deadcrootn_storage_to_deadcrootn_transfer * dt;
            ns->deadcrootn_storage -=
                nf->deadcrootn_storage_to_deadcrootn_transfer * dt;
        }
        /* for deciduous system, force leafn and frootn to exactly 0.0 on the
         * last day */
        if (!evergreen)
        {
            if (ns->leafn < 1.0e-10)
            {
                ns->leafn = 0.0;
            }
            if (ns->frootn < 1.0e-10)
            {
                ns->frootn = 0.0;
            }
        }
    }                           /* end if annual allocation day */
}
