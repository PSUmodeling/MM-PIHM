#include "pihm.h"

void MaintResp (const epconst_struct *epc, epvar_struct *epv,
    const estate_struct *es, const eflux_struct *ef, const cstate_struct *cs,
     cflux_struct *cf, const nstate_struct *ns)
{
    /*
     * Maintenance respiration routine
     * Uses reference values at 20 deg C and an empirical relationship between
     * tissue N content and respiration rate given in:
     * Ryan, M.G., 1991. Effects of climate change on plant respiration
     * Ecological Applications, 1(2):157-167.
     *
     * Uses the same value of Q_10 (2.0) for all compartments, leaf, stem,
     * coarse and fine roots.
     *
     * From Ryan's figures and regressions equations, the maintenance respiration
     * in kgC/day per kg of tissue N is:
     * mrpern = 0.218 (kgC/kgN/d)
     * Convert to per second:
     * mrpern = 2.524E-6
     *
     * Leaf maintenance respiration is calculated separately for day and night,
     * since the PSN routine needs the daylight value.
     * Leaf and fine root respiration are dependent on phenology.
     */

    double          tair;
    double          tsoil;
    double          t1;
    const double    Q10 = 2.0;
    const double    MRPERN = 2.524E-6;
    double          exponent;
    double          n_area_sun, n_area_shade;
    double          dlmr_area_sun, dlmr_area_shade;

    tair = es->sfctmp - TFREEZ;
    tsoil = es->stc[0] - TFREEZ;

    /* Leaf maintenance respiration when leaves on */
    if (cs->leafc)
    {
        t1 = ns->leafn * MRPERN;

        exponent = (tair - 20.0) / 10.0;
        cf->leaf_day_mr = t1 * pow (Q10, exponent);

        /* For day respiration, also determine rates of maintenance respiration
         * per unit of projected leaf area in the sunlit and shaded portions of
         * the canopy, for use in the photosynthesis routine */
        /* First, calculate the mass of N per unit of projected leaf area
         * in each canopy fraction (kg N/m2 projected area) */
        if (ef->soldn > 0.0)
        {
            n_area_sun = 1.0 / (epv->sun_proj_sla * epc->leaf_cn);
            n_area_shade = 1.0 / (epv->shade_proj_sla * epc->leaf_cn);
            /* Convert to respiration flux in kg C/m2 projected area/day, and
             * correct for temperature */
            dlmr_area_sun = n_area_sun * MRPERN * pow (Q10, exponent);
            dlmr_area_shade = n_area_shade * MRPERN * pow (Q10, exponent);
            /* Finally, convert from mass to molar units (umol s-1) */
            epv->dlmr_area_sun = dlmr_area_sun / 12.011e-9;
            epv->dlmr_area_shade = dlmr_area_shade / 12.011e-9;
        }
        else
        {
            epv->dlmr_area_sun = 0.0;
            epv->dlmr_area_shade = 0.0;
        }
    }
    else                        /* No leaves on */
    {
        cf->leaf_day_mr = 0.0;
        epv->dlmr_area_sun = 0.0;
        epv->dlmr_area_shade = 0.0;
    }

    /* Fine root maintenance respiration when fine roots on */
    /* Ammended to consider only the specified n concentration,
     * to avoid excessive MR with n-loading to fine roots */
    if (cs->frootc)
    {
        exponent = (tsoil - 20.0) / 10.0;
        t1 = pow (Q10, exponent);
        cf->froot_mr = ns->frootn * MRPERN * t1;
    }
    else                        /* No fine roots on */
    {
        cf->froot_mr = 0.0;
    }

    /* Tree-specific fluxes */
    if (epc->woody)
    {
        /* Live stem maintenance respiration */
        exponent = (tair - 20.0) / 10.0;
        t1 = pow (Q10, exponent);
        cf->livestem_mr = ns->livestemn * MRPERN * t1;

        /* Live coarse root maintenance respiration */
        exponent = (tair - 20.0) / 10.0;
        t1 = pow (Q10, exponent);
        cf->livecroot_mr = ns->livecrootn * MRPERN * t1;
    }
}
