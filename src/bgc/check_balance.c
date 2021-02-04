#include "pihm.h"

void CheckCarbonBalance (const cstate_struct *cs, double *old_balance)
{
    double          in, out, store, balance;

    // DAILY CHECK ON CARBON BALANCE

    // Sum of sources
    in = cs->psnsun_src + cs->psnshade_src;

    // Sum of sinks
    out = cs->leaf_mr_snk + cs->leaf_gr_snk + cs->froot_mr_snk + cs->froot_gr_snk + cs->livestem_mr_snk +
        cs->livestem_gr_snk + cs->deadstem_gr_snk + cs->livecroot_mr_snk + cs->livecroot_gr_snk + cs->deadcroot_gr_snk +
        cs->litr1_hr_snk + cs->litr2_hr_snk + cs->litr4_hr_snk + cs->soil1_hr_snk + cs->soil2_hr_snk +
        cs->soil3_hr_snk + cs->soil4_hr_snk + cs->fire_snk;

    // Sum of current storage
    store = cs->leafc + cs->leafc_storage + cs->leafc_transfer + cs->frootc + cs->frootc_storage + cs->frootc_transfer +
        cs->livestemc + cs->livestemc_storage + cs->livestemc_transfer + cs->deadstemc + cs->deadstemc_storage +
        cs->deadstemc_transfer + cs->livecrootc + cs->livecrootc_storage + cs->livecrootc_transfer + cs->deadcrootc +
        cs->deadcrootc_storage + cs->deadcrootc_transfer + cs->gresp_storage + cs->gresp_transfer + cs->cwdc +
        cs->litr1c + cs->litr2c + cs->litr3c + cs->litr4c + cs->soil1c + cs->soil2c + cs->soil3c + cs->soil4c +
        cs->cpool;

    // Calculate current balance
    balance = in - out - store;

    if (!first_balance)
    {
        if (fabs(*old_balance - balance) > 1E-8)
        {
            pihm_printf (VL_ERROR, "FATAL ERRROR: carbon balance error:\n");
            pihm_printf (VL_ERROR, "Balance from previous day = %lf\n", *old_balance);
            pihm_printf (VL_ERROR, "Balance from current day  = %lf\n", balance);
            pihm_printf (VL_ERROR, "Difference (previous - current) = %lg\n", *old_balance - balance);
            pihm_printf (VL_ERROR, "Components of current balance:\n");
            pihm_printf (VL_ERROR, "Sources (summed over entire run)  = %lf\n", in);
            pihm_printf (VL_ERROR, "Sinks   (summed over entire run)  = %lf\n", out);
            pihm_printf (VL_ERROR, "Storage (current state variables) = %lf\n", store);
            pihm_exit (EXIT_FAILURE);
        }
    }
    *old_balance = balance;
}

void CheckNitrogenBalance (const nstate_struct *ns, double *old_balance)
{
    double          in, out, store, balance;

    // DAILY CHECK ON NITROGEN BALANCE

    // Sum of sources
    in = ns->nfix_src + ns->ndep_src;

    // Sum of sinks
    out = ns->nleached_snk + ns->nvol_snk + ns->fire_snk;

    // Sum of current storage
    store = ns->leafn + ns->leafn_storage + ns->leafn_transfer + ns->frootn + ns->frootn_storage + ns->frootn_transfer +
        ns->livestemn + ns->livestemn_storage + ns->livestemn_transfer + ns->deadstemn + ns->deadstemn_storage +
        ns->deadstemn_transfer + ns->livecrootn + ns->livecrootn_storage + ns->livecrootn_transfer + ns->deadcrootn +
        ns->deadcrootn_storage + ns->deadcrootn_transfer + ns->cwdn + ns->litr1n + ns->litr2n + ns->litr3n +
        ns->litr4n + ns->soil1n + ns->soil2n + ns->soil3n + ns->soil4n + ns->sminn + ns->npool + ns->retransn;

    // Calculate current balance
    balance = in - out - store;

    if (!first_balance)
    {
        if (fabs(*old_balance - balance) > 1E-8)
        {
            pihm_printf (VL_ERROR, "FATAL ERRROR: nitrogen balance error:\n");
            pihm_printf (VL_ERROR, "Balance from previous day = %lg\n", *old_balance);
            pihm_printf (VL_ERROR, "Balance from current day  = %lg\n", balance);
            pihm_printf (VL_ERROR, "Difference (previous - current) = %lg\n", *old_balance - balance);
            pihm_printf (VL_ERROR, "Components of current balance:\n");
            pihm_printf (VL_ERROR, "Sources (summed over entire run)  = %lf\n", in);
            pihm_printf (VL_ERROR, "Sinks   (summed over entire run)  = %lf\n", out);
            pihm_printf (VL_ERROR, "Storage (current state variables) = %lf\n", store);
            pihm_exit (EXIT_FAILURE);
        }
    }
    *old_balance = balance;
}
