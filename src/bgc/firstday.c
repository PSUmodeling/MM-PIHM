#include "pihm.h"

void FirstDay(const cninit_struct *cninit, elem_struct elem[], river_struct river[])
{
    int             i;

#if defined(_LUMPEDBGC_)
    i = LUMPEDBGC;
#else
    for (i = 0; i < nelem; i++)
#endif
    {
        bgcic_struct   *restart;
        epconst_struct *epc;
        double          max_leafc, max_frootc, max_stemc;
        double          new_stemc;

        restart = &elem[i].restart_input;
        epc = &elem[i].epc;

        /*
         * Copy from CN initialization structure
         */
        restart->cwdc = cninit->cwdc;
        restart->litr1c = cninit->litr1c;
        restart->litr2c = cninit->litr2c;
        restart->litr3c = cninit->litr3c;
        restart->litr4c = cninit->litr4c;
        restart->soil1c = cninit->soil1c;
        restart->soil2c = cninit->soil2c;
        restart->soil3c = cninit->soil3c;
        restart->soil4c = cninit->soil4c;
        restart->litr1n = cninit->litr1n;
        restart->sminn = cninit->sminn;

        // Calculate N states from C states
        restart->cwdn = restart->cwdc / epc->deadwood_cn;
        restart->litr2n = restart->litr2c / epc->leaflitr_cn;
        restart->litr3n = restart->litr3c / epc->leaflitr_cn;
        restart->litr4n = restart->litr4c / epc->leaflitr_cn;
        restart->soil1n = restart->soil1c / SOIL1_CN;
        restart->soil2n = restart->soil2c / SOIL2_CN;
        restart->soil3n = restart->soil3c / SOIL3_CN;
        restart->soil4n = restart->soil4c / SOIL4_CN;

        // Set phenology flags
        restart->dormant_flag = (epc->evergreen) ? 0 : 1;
        restart->onset_flag = 0;
        restart->onset_counter = 0;
        restart->onset_gddflag = 0;
        restart->onset_fdd = 0.0;
        restart->onset_gdd = 0.0;
        restart->onset_swi = 0.0;
        restart->offset_flag = 0;
        restart->offset_counter = 0;
        restart->offset_fdd = 0.0;
        restart->offset_swi = 0.0;

        // Initialize other C and N storage state variables
        restart->leafc_storage = 0.0;
        restart->frootc_storage = 0.0;
        restart->livestemc_storage = 0.0;
        restart->deadstemc_storage = 0.0;
        restart->livecrootc_storage = 0.0;
        restart->deadcrootc_storage = 0.0;
        restart->gresp_storage = 0.0;
        restart->cpool = 0.0;
        restart->leafn_storage = 0.0;
        restart->frootn_storage = 0.0;
        restart->livestemn_storage = 0.0;
        restart->deadstemn_storage = 0.0;
        restart->livecrootn_storage = 0.0;
        restart->deadcrootn_storage = 0.0;
        restart->retransn = 0.0;
        restart->npool = 0.0;

        // Initialize days-since-rain counter
        restart->dsr = 0.0;

        // Establish the initial partitioning between displayed growth and growth ready for transfer
        max_leafc = cninit->max_leafc;
        restart->leafc_transfer = max_leafc * epc->leaf_turnover;
        restart->leafc = max_leafc - restart->leafc_transfer;
        max_frootc = max_leafc * epc->alloc_frootc_leafc;
        restart->frootc_transfer = cninit->max_leafc * epc->alloc_frootc_leafc * epc->froot_turnover;
        restart->frootc = max_frootc - restart->frootc_transfer;
        if (epc->woody)
        {
            max_stemc = cninit->max_stemc;
            new_stemc = restart->leafc_transfer * epc->alloc_newstemc_newleafc;
            restart->livestemc_transfer = new_stemc * epc->alloc_newlivewoodc_newwoodc;
            restart->livestemc = restart->livestemc_transfer / epc->livewood_turnover;
            restart->deadstemc_transfer = new_stemc - restart->livestemc_transfer;
            restart->deadstemc =
                max_stemc - restart->livestemc_transfer - restart->livestemc - restart->deadstemc_transfer;
            restart->deadstemc = MAX(restart->deadstemc, 0.0);
            restart->livecrootc_transfer = restart->livestemc_transfer * epc->alloc_crootc_stemc;
            restart->livecrootc = restart->livestemc * epc->alloc_crootc_stemc;
            restart->deadcrootc_transfer = restart->deadstemc_transfer * epc->alloc_crootc_stemc;
            restart->deadcrootc = restart->deadstemc * epc->alloc_crootc_stemc;
        }

        // Calculate initial leaf and froot nitrogen pools from carbon pools and user-specified initial C:N for each
        // component
        restart->leafn_transfer = restart->leafc_transfer / epc->leaf_cn;
        restart->leafn = restart->leafc / epc->leaf_cn;
        restart->frootn_transfer = restart->frootc_transfer / epc->froot_cn;
        restart->frootn = restart->frootc / epc->froot_cn;
        if (epc->woody)
        {
            restart->livestemn_transfer = restart->livestemc_transfer / epc->livewood_cn;
            restart->livestemn = restart->livestemc / epc->livewood_cn;
            restart->deadstemn_transfer = restart->deadstemc_transfer / epc->deadwood_cn;
            restart->deadstemn = restart->deadstemc / epc->deadwood_cn;
            restart->livecrootn_transfer = restart->livecrootc_transfer / epc->livewood_cn;
            restart->livecrootn = restart->livecrootc / epc->livewood_cn;
            restart->deadcrootn_transfer = restart->deadcrootc_transfer / epc->deadwood_cn;
            restart->deadcrootn = restart->deadcrootc / epc->deadwood_cn;
        }

        // Add the growth respiration requirement for the first year's leaf and fine root growth from transfer pools to
        // the gresp_transfer pool
        restart->gresp_transfer = 0.0;
        restart->gresp_transfer += (restart->leafc_transfer + restart->frootc_transfer) * GRPERC;
        if (epc->woody)
        {
            restart->gresp_transfer += (restart->livestemc_transfer + restart->deadstemc_transfer +
                restart->livecrootc_transfer + restart->deadcrootc_transfer) * GRPERC;
        }

        // Set the initial rates of litterfall and live wood turnover
        restart->prev_leafc_to_litter = 0.0;
        restart->prev_frootc_to_litter = 0.0;
    }

#if !defined(_LUMPEDBGC_) && !defined(_LEACHING_)
    for (i = 0; i < nriver; i++)
    {
        river[i].restart_input.streamn = 0.0;
    }
#endif
}
