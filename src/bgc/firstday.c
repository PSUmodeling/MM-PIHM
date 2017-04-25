#include "pihm.h"

void FirstDay (const epconst_struct *epc, const cinit_struct *cinit,
    epvar_struct *epv, cstate_struct *cs, nstate_struct *ns)
{
    int             predays, remdays;
    double          max_leafc, max_frootc;
    double          max_stemc, new_stemc;
    double          prop_transfer, transfer;
    double          prop_litfall;

    /* Initialize the c and n storage state variables */
    cs->leafc_storage = 0.0;
    cs->frootc_storage = 0.0;
    cs->livestemc_storage = 0.0;
    cs->deadstemc_storage = 0.0;
    cs->livecrootc_storage = 0.0;
    cs->deadcrootc_storage = 0.0;
    cs->gresp_storage = 0.0;
    cs->cpool = 0.0;
    ns->leafn_storage = 0.0;
    ns->frootn_storage = 0.0;
    ns->livestemn_storage = 0.0;
    ns->deadstemn_storage = 0.0;
    ns->livecrootn_storage = 0.0;
    ns->deadcrootn_storage = 0.0;
    ns->retransn = 0.0;
    ns->npool = 0.0;

    /* Initialize days-since-rain counter */
    epv->dsr = 0.0;

    /* Establish the initial partitioning between displayed growth and growth
     * ready for transfer */
    max_leafc = cinit->max_leafc;
    cs->leafc_transfer = max_leafc * epc->leaf_turnover;
    cs->leafc = max_leafc - cs->leafc_transfer;
    max_frootc = max_leafc * epc->alloc_frootc_leafc;
    cs->frootc_transfer =
        cinit->max_leafc * epc->alloc_frootc_leafc * epc->froot_turnover;
    cs->frootc = max_frootc - cs->frootc_transfer;
    if (epc->woody)
    {
        max_stemc = cinit->max_stemc;
        new_stemc = cs->leafc_transfer * epc->alloc_newstemc_newleafc;
        cs->livestemc_transfer = new_stemc * epc->alloc_newlivewoodc_newwoodc;
        cs->livestemc = cs->livestemc_transfer / epc->livewood_turnover;
        cs->deadstemc_transfer = new_stemc - cs->livestemc_transfer;
        cs->deadstemc = max_stemc - cs->livestemc_transfer - cs->livestemc -
            cs->deadstemc_transfer;
        if (cs->deadstemc < 0.0)
        {
            cs->deadstemc = 0.0;
        }
        cs->livecrootc_transfer =
            cs->livestemc_transfer * epc->alloc_crootc_stemc;
        cs->livecrootc = cs->livestemc * epc->alloc_crootc_stemc;
        cs->deadcrootc_transfer =
            cs->deadstemc_transfer * epc->alloc_crootc_stemc;
        cs->deadcrootc = cs->deadstemc * epc->alloc_crootc_stemc;
    }

    /* Calculate initial leaf and froot nitrogen pools from carbon pools and
     * user-specified initial C:N for each component */
    ns->leafn_transfer = cs->leafc_transfer / epc->leaf_cn;
    ns->leafn = cs->leafc / epc->leaf_cn;
    ns->frootn_transfer = cs->frootc_transfer / epc->froot_cn;
    ns->frootn = cs->frootc / epc->froot_cn;
    if (epc->woody)
    {
        ns->livestemn_transfer = cs->livestemc_transfer / epc->livewood_cn;
        ns->livestemn = cs->livestemc / epc->livewood_cn;
        ns->deadstemn_transfer = cs->deadstemc_transfer / epc->deadwood_cn;
        ns->deadstemn = cs->deadstemc / epc->deadwood_cn;
        ns->livecrootn_transfer = cs->livecrootc_transfer / epc->livewood_cn;
        ns->livecrootn = cs->livecrootc / epc->livewood_cn;
        ns->deadcrootn_transfer = cs->deadcrootc_transfer / epc->deadwood_cn;
        ns->deadcrootn = cs->deadcrootc / epc->deadwood_cn;
    }

    /* Add the growth respiration requirement for the first year's
     * leaf and fine root growth from transfer pools to the
     * gresp_transfer pool */
    cs->gresp_transfer = 0.0;
    cs->gresp_transfer += (cs->leafc_transfer + cs->frootc_transfer) * GRPERC;
    if (epc->woody)
    {
        cs->gresp_transfer +=
            (cs->livestemc_transfer + cs->deadstemc_transfer +
            cs->livecrootc_transfer + cs->deadcrootc_transfer) * GRPERC;
    }

    /* Set the initial rates of litterfall and live wood turnover */
    if (epc->evergreen)
    {
        /* Leaf and fineroot litterfall rates */
        epv->bg_leafc_litfall_rate =
            max_leafc * epc->leaf_turnover / 365.0 / DAYINSEC;
        epv->bg_frootc_litfall_rate =
            max_frootc * epc->froot_turnover / 365.0 / DAYINSEC;
        epv->prev_leafc_to_litter = 0.0;
        epv->prev_frootc_to_litter = 0.0;
    }
    else
    {
        /* Deciduous: reset the litterfall rates to 0.0 for the start of the
         * next litterfall season */
        epv->bg_leafc_litfall_rate = 0.0;
        epv->bg_frootc_litfall_rate = 0.0;
        epv->prev_leafc_to_litter = 0.0;
        epv->prev_frootc_to_litter = 0.0;
    }

    if (epc->woody)
    {
        /* Live wood turnover rates */
        epv->livestemc_turnover_rate =
            cs->livestemc * epc->livewood_turnover / 365.0 / DAYINSEC;
        epv->livecrootc_turnover_rate =
            cs->livecrootc * epc->livewood_turnover / 365.0 / DAYINSEC;
    }
}
