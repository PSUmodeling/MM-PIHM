#include "pihm.h"

void DailyAllocation(const epconst_struct *epc, const cstate_struct *cs,
    const nstate_struct *ns, epvar_struct *epv, cflux_struct *cf,
    nflux_struct *nf, ntemp_struct *nt)
{
    double          day_gpp;      // daily gross production
    double          day_mresp;    // daily total maintenance respiration
    double          avail_c;      // total C available for new production
    double          f1;           // RATIO new fine root C : new leaf C
    double          f2;           // RATIO new coarse root C : new stem C
    double          f3;           // RATIO new stem C : new leaf C
    double          f4;           // RATIO new live wood C : new wood C
    double          g1;           // RATIO C respired for growth : C grown
    double          g2;           // proportion of growth resp to release at fixation
    double          cnl;          // RATIO leaf C:N
    double          cnfr;         // RATIO fine root C:N
    double          cnlw;         // RATIO live wood C:N
    double          cndw;         // RATIO dead wood C:N
    double          nlc;          // actual new leaf C, minimum of C and N limits
    double          pnow;         // proportion of growth displayed on current day
    double          gresp_storage;
    double          c_allometry, n_allometry;
    double          plant_ndemand, sum_ndemand;
    double          actual_immob = 0.0;
    double          plant_nalloc, plant_calloc;
    double          fpi = 0.0;
    double          plant_remaining_ndemand;
    double          excess_c;
    int             nlimit;
    double          cn_l1 = 0.0, cn_l2 = 0.0, cn_l4 = 0.0;
    double          cn_s1, cn_s2, cn_s3, cn_s4;
    double          rfl1s1, rfl2s2, rfl4s3, rfs1s2, rfs2s3, rfs3s4;
    double          avail_retransn;
    double          cpool_recovery;
    double          excessn;

    // Assess the carbon availability on the basis of this day's gross production and maintenance respiration costs
    day_gpp = cf->psnsun_to_cpool + cf->psnshade_to_cpool;
    day_mresp = (epc->woody) ?
        cf->leaf_day_mr + cf->leaf_night_mr + cf->froot_mr + cf->livestem_mr + cf->livecroot_mr :
        cf->leaf_day_mr + cf->leaf_night_mr + cf->froot_mr;

    avail_c = day_gpp - day_mresp;

    // No allocation when the daily C balance is negative
    avail_c = MAX(avail_c, 0.0);

    // Test for cpool deficit
    if (cs->cpool < 0.0)
    {
        // Running a deficit in cpool, so the first priority is to let some of today's available C accumulate in cpool.
        // The actual accumulation in the cpool is resolved in day_carbon_state().
        // First determine how much of the deficit should be recovered today
        cpool_recovery = -cs->cpool / DAYSCRECOVER;

        // Potential recovery of cpool deficit is less than the available carbon for the day, so alleviate cpool deficit
        // and use the rest of the available carbon for new growth and storage. Remember that fluxes in and out of the
        // cpool are reconciled at the end of the daily loop, so for now, just keep track of the amount of daily GPP-MR
        // that is not needed to restore a negative cpool.
        if (cpool_recovery < avail_c)
        {
            avail_c -= cpool_recovery;
        }
        // cpool deficit is >= available C, so all of the daily GPP, if any, is used to alleviate negative cpool
        else
        {
            avail_c = 0.0;
        }
    }   // End if negative cpool

    // Assign local values for the allocation control parameters
    f1 = epc->alloc_frootc_leafc;
    f2 = epc->alloc_crootc_stemc;
    f3 = epc->alloc_newstemc_newleafc;
    f4 = epc->alloc_newlivewoodc_newwoodc;
    g1 = GRPERC;
    g2 = GRPNOW;
    cnl = epc->leaf_cn;
    cnfr = epc->froot_cn;
    cnlw = epc->livewood_cn;
    cndw = epc->deadwood_cn;
    pnow = epc->alloc_prop_curgrowth;

    // Given the available C, use constant allometric relationships to determine how much N is required to meet this
    // potential growth demand
    if (epc->woody)
    {
        c_allometry = ((1.0 + g1) * (1.0 + f1 + f3 * (1.0 + f2)));
        n_allometry = 1.0 / cnl + f1 / cnfr + (f3 * f4 * (1.0 + f2)) / cnlw + (f3 * (1.0 - f4) * (1.0 + f2)) / cndw;
    }
    else
    {
        c_allometry = 1.0 + g1 + f1 + f1 * g1;
        n_allometry = 1.0 / cnl + f1 / cnfr;
    }
    plant_ndemand = avail_c * (n_allometry / c_allometry);

    // Now compare the combined decomposition immobilization and plant growth N demands against the available soil
    // mineral N pool.
    avail_retransn = ns->retransn / DAYSNDEPLOY;
    sum_ndemand = plant_ndemand + nt->potential_immob;

    if (sum_ndemand <= ns->sminn)
    {
        // N availability is not limiting immobilization or plant uptake, and both can proceed at their potential rates
        actual_immob = nt->potential_immob;
        nlimit = 0;

        // Determine the split between retranslocation N and soil mineral N to meet the plant demand
        nf->retransn_to_npool = (plant_ndemand > avail_retransn) ? avail_retransn : plant_ndemand;

        nf->sminn_to_npool = plant_ndemand - nf->retransn_to_npool;
        plant_nalloc = nf->retransn_to_npool + nf->sminn_to_npool;
        plant_calloc = avail_c;

        // Under conditions of excess N, some proportion of excess N is assumed to be lost to denitrification, in
        // addition to the constant proportion lost in the decomposition pathways.
        excessn = ns->sminn - sum_ndemand;
        nf->sminn_to_denitrif = excessn * BULK_DENITRIF_PROPORTION;
    }
    else
    {
        // N availability can not satisfy the sum of immobilization and plant growth demands, so these two demands
        // compete for available soil mineral N
        nlimit = 1;
        actual_immob = (sum_ndemand) ? ns->sminn * (nt->potential_immob / sum_ndemand) : actual_immob;
        fpi = (nt->potential_immob) ? actual_immob / nt->potential_immob : 0.0;
        nf->sminn_to_npool = ns->sminn - actual_immob;
        plant_remaining_ndemand = plant_ndemand - nf->sminn_to_npool;
        // The demand not satisfied by uptake from soil mineral N is now sought from the retranslocated N pool
        if (plant_remaining_ndemand <= avail_retransn)
        {
            // There is enough N available from retranslocation pool to satisfy the remaining plant N demand
            nf->retransn_to_npool = plant_remaining_ndemand;
            plant_calloc = avail_c;
        }
        else
        {
            // There is not enough available retranslocation N to satisfy the entire demand. In this case, the remaining
            // unsatisfied N demand is translated back to a C excess, which is deducted proportionally from the sun and
            // shade photosynthesis source terms
            nf->retransn_to_npool = avail_retransn;
            plant_nalloc = nf->retransn_to_npool + nf->sminn_to_npool;
            plant_calloc = plant_nalloc * (c_allometry / n_allometry);
            excess_c = avail_c - plant_calloc;
            cf->psnsun_to_cpool -= excess_c * (cf->psnsun_to_cpool / day_gpp);
            cf->psnshade_to_cpool -= excess_c * (cf->psnshade_to_cpool / day_gpp);
        }
    }

    // Calculate the amount of new leaf C dictated by these allocation decisions, and figure the daily fluxes of C and N
    // to current growth and storage pools
    // Pnow is the proportion of this day's growth that is displayed now, the remainder going into storage for display
    // next year through the transfer pools
    nlc = plant_calloc / c_allometry;
    // Daily C fluxes out of cpool and into new growth or storage
    cf->cpool_to_leafc = nlc * pnow;
    cf->cpool_to_leafc_storage = nlc * (1.0 - pnow);
    cf->cpool_to_frootc = nlc * f1 * pnow;
    cf->cpool_to_frootc_storage = nlc * f1 * (1.0 - pnow);
    if (epc->woody)
    {
        cf->cpool_to_livestemc = nlc * f3 * f4 * pnow;
        cf->cpool_to_livestemc_storage = nlc * f3 * f4 * (1.0 - pnow);
        cf->cpool_to_deadstemc = nlc * f3 * (1.0 - f4) * pnow;
        cf->cpool_to_deadstemc_storage = nlc * f3 * (1.0 - f4) * (1.0 - pnow);
        cf->cpool_to_livecrootc = nlc * f2 * f3 * f4 * pnow;
        cf->cpool_to_livecrootc_storage = nlc * f2 * f3 * f4 * (1.0 - pnow);
        cf->cpool_to_deadcrootc = nlc * f2 * f3 * (1.0 - f4) * pnow;
        cf->cpool_to_deadcrootc_storage =
            nlc * f2 * f3 * (1.0 - f4) * (1.0 - pnow);
    }
    // Daily N fluxes out of npool and into new growth or storage
    nf->npool_to_leafn = (nlc / cnl) * pnow;
    nf->npool_to_leafn_storage = (nlc / cnl) * (1.0 - pnow);
    nf->npool_to_frootn = (nlc * f1 / cnfr) * pnow;
    nf->npool_to_frootn_storage = (nlc * f1 / cnfr) * (1.0 - pnow);
    if (epc->woody)
    {
        nf->npool_to_livestemn = (nlc * f3 * f4 / cnlw) * pnow;
        nf->npool_to_livestemn_storage = (nlc * f3 * f4 / cnlw) * (1.0 - pnow);
        nf->npool_to_deadstemn = (nlc * f3 * (1.0 - f4) / cndw) * pnow;
        nf->npool_to_deadstemn_storage = (nlc * f3 * (1.0 - f4) / cndw) * (1.0 - pnow);
        nf->npool_to_livecrootn = (nlc * f2 * f3 * f4 / cnlw) * pnow;
        nf->npool_to_livecrootn_storage = (nlc * f2 * f3 * f4 / cnlw) * (1.0 - pnow);
        nf->npool_to_deadcrootn = (nlc * f2 * f3 * (1.0 - f4) / cndw) * pnow;
        nf->npool_to_deadcrootn_storage = (nlc * f2 * f3 * (1.0 - f4) / cndw) * (1.0 - pnow);
    }

    // Calculate the amount of carbon that needs to go into growth respiration storage to satisfy all of the storage
    // growth demands. Note that in version 4.1, this function has been changed to allow for the fraction of growth
    // respiration that is released at the time of fixation, versus the remaining fraction that is stored for release at
    // the time of display.
    // Note that all the growth respiration fluxes that get released on a given day are calculated in growth_resp(), but
    // that the storage of C for growth resp during display of transferred growth is assigned here.
    gresp_storage = (epc->woody) ?
        (cf->cpool_to_leafc_storage + cf->cpool_to_frootc_storage + cf->cpool_to_livestemc_storage +
        cf->cpool_to_deadstemc_storage + cf->cpool_to_livecrootc_storage + cf->cpool_to_deadcrootc_storage) *
        g1 * (1.0 - g2) :
        (cf->cpool_to_leafc_storage + cf->cpool_to_frootc_storage) * g1 * (1.0 - g2);
    cf->cpool_to_gresp_storage = gresp_storage;

    // Now use the N limitation information to assess the final decomposition fluxes. Mineralizing fluxes (pmnf 0.0)
    // occur at the potential rate regardless of the competing N demands between microbial processes and plant uptake,
    // but immobilizing fluxes are reduced when soil mineral N is limiting
    // Calculate litter and soil compartment C:N ratios
    cn_l1 = (ns->litr1n > 0.0) ? cs->litr1c / ns->litr1n : cn_l1;
    cn_l2 = (ns->litr2n > 0.0) ? cs->litr2c / ns->litr2n : cn_l2;
    cn_l4 = (ns->litr4n > 0.0) ? cs->litr4c / ns->litr4n : cn_l4;
    cn_s1 = SOIL1_CN;
    cn_s2 = SOIL2_CN;
    cn_s3 = SOIL3_CN;
    cn_s4 = SOIL4_CN;
    // Respiration fractions for fluxes between compartments
    rfl1s1 = RFL1S1;
    rfl2s2 = RFL2S2;
    rfl4s3 = RFL4S3;
    rfs1s2 = RFS1S2;
    rfs2s3 = RFS2S3;
    rfs3s4 = RFS3S4;

    // Labile litter fluxes
    if (cs->litr1c > 0.0)
    {
        if (nlimit && nt->pmnf_l1s1 > 0.0)
        {
            nt->plitr1c_loss *= fpi;
            nt->pmnf_l1s1 *= fpi;
        }
        cf->litr1_hr = rfl1s1 * nt->plitr1c_loss;
        cf->litr1c_to_soil1c = (1.0 - rfl1s1) * nt->plitr1c_loss;
        nf->litr1n_to_soil1n = (ns->litr1n > 0.0) ? nt->plitr1c_loss / cn_l1 : 0.0;
        nf->sminn_to_soil1n_l1 = nt->pmnf_l1s1;
    }

    // Cellulose litter fluxes
    if (cs->litr2c > 0.0)
    {
        if (nlimit && nt->pmnf_l2s2 > 0.0)
        {
            nt->plitr2c_loss *= fpi;
            nt->pmnf_l2s2 *= fpi;
        }
        cf->litr2_hr = rfl2s2 * nt->plitr2c_loss;
        cf->litr2c_to_soil2c = (1.0 - rfl2s2) * nt->plitr2c_loss;
        nf->litr2n_to_soil2n = (ns->litr2n > 0.0) ? nt->plitr2c_loss / cn_l2 : 0.0;
        nf->sminn_to_soil2n_l2 = nt->pmnf_l2s2;
    }

    // Release of shielded cellulose litter, tied to the decay rate of lignin litter
    if (cs->litr3c > 0.0)
    {
        if (nlimit && nt->pmnf_l4s3 > 0.0)
        {
            cf->litr3c_to_litr2c = nt->kl4 * cs->litr3c * fpi;
            nf->litr3n_to_litr2n = nt->kl4 * ns->litr3n * fpi;
        }
        else
        {
            cf->litr3c_to_litr2c = nt->kl4 * cs->litr3c;
            nf->litr3n_to_litr2n = nt->kl4 * ns->litr3n;
        }
    }

    // Lignin litter fluxes
    if (cs->litr4c > 0.0)
    {
        if (nlimit && nt->pmnf_l4s3 > 0.0)
        {
            nt->plitr4c_loss *= fpi;
            nt->pmnf_l4s3 *= fpi;
        }
        cf->litr4_hr = rfl4s3 * nt->plitr4c_loss;
        cf->litr4c_to_soil3c = (1.0 - rfl4s3) * nt->plitr4c_loss;
        nf->litr4n_to_soil3n = (ns->litr4n > 0.0) ? nt->plitr4c_loss / cn_l4 : 0.0;
        nf->sminn_to_soil3n_l4 = nt->pmnf_l4s3;
    }

    // Fast microbial recycling pool
    if (cs->soil1c > 0.0)
    {
        if (nlimit && nt->pmnf_s1s2 > 0.0)
        {
            nt->psoil1c_loss *= fpi;
            nt->pmnf_s1s2 *= fpi;
        }
        cf->soil1_hr = rfs1s2 * nt->psoil1c_loss;
        cf->soil1c_to_soil2c = (1.0 - rfs1s2) * nt->psoil1c_loss;
        nf->soil1n_to_soil2n = nt->psoil1c_loss / cn_s1;
        nf->sminn_to_soil2n_s1 = nt->pmnf_s1s2;
    }

    // Medium microbial recycling pool
    if (cs->soil2c > 0.0)
    {
        if (nlimit && nt->pmnf_s2s3 > 0.0)
        {
            nt->psoil2c_loss *= fpi;
            nt->pmnf_s2s3 *= fpi;
        }
        cf->soil2_hr = rfs2s3 * nt->psoil2c_loss;
        cf->soil2c_to_soil3c = (1.0 - rfs2s3) * nt->psoil2c_loss;
        nf->soil2n_to_soil3n = nt->psoil2c_loss / cn_s2;
        nf->sminn_to_soil3n_s2 = nt->pmnf_s2s3;
    }

    // Slow microbial recycling pool
    if (cs->soil3c > 0.0)
    {
        if (nlimit && nt->pmnf_s3s4 > 0.0)
        {
            nt->psoil3c_loss *= fpi;
            nt->pmnf_s3s4 *= fpi;
        }
        cf->soil3_hr = rfs3s4 * nt->psoil3c_loss;
        cf->soil3c_to_soil4c = (1.0 - rfs3s4) * nt->psoil3c_loss;
        nf->soil3n_to_soil4n = nt->psoil3c_loss / cn_s3;
        nf->sminn_to_soil4n_s3 = nt->pmnf_s3s4;
    }

    // Recalcitrant SOM pool (rf = 1.0, always mineralizing)
    if (cs->soil4c > 0.0)
    {
        cf->soil4_hr = nt->psoil4c_loss;
        nf->soil4n_to_sminn = nt->psoil4c_loss / cn_s4;
    }

    epv->fpi = fpi;
}
