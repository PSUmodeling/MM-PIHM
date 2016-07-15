#include "pihm.h"

void MaintResp (const cstate_struct *cs, const nstate_struct *ns, const epconst_struct *epc, const daily_struct *daily, cflux_struct *cf, epvar_struct *epv)
{
    /*
     * maintenance respiration routine
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
     *
     * Leaf maintenance respiration is calculated separately for day and night,
     * since the PSN routine needs the daylight value.
     * Leaf and fine root respiration are dependent on phenology.
     */

    double          tday, tnight;
    double          tavg;
    double          tsoil;
    double          t1;
    double          q10 = 2.0;
    double          mrpern = 0.218;
    double          exponent;
    double          n_area_sun, n_area_shade, dlmr_area_sun, dlmr_area_shade;

    tday = daily->tday - TFREEZ;
    tnight = daily->tnight - TFREEZ;
    tavg = daily->avg_sfctmp - TFREEZ;
    tsoil = daily->avg_stc[0] - TFREEZ;

    /* leaf day and night maintenance respiration when leaves on */
    if (cs->leafc)
    {
        t1 = ns->leafn * mrpern;

        /* leaf, day */
        exponent = (tday - 20.0) / 10.0;
        cf->leaf_day_mr = t1 * pow (q10, exponent) * daily->dayl / 86400.0;

        /* for day respiration, also determine rates of maintenance respiration
         * per unit of projected leaf area in the sunlit and shaded portions of
         * the canopy, for use in the photosynthesis routine */
        /* first, calculate the mass of N per unit of projected leaf area
         * in each canopy fraction (kg N/m2 projected area) */
        n_area_sun = 1.0 / (epv->sun_proj_sla * epc->leaf_cn);
        n_area_shade = 1.0 / (epv->shade_proj_sla * epc->leaf_cn);
        /* convert to respiration flux in kg C/m2 projected area/day, and
         * correct for temperature */
        dlmr_area_sun = n_area_sun * mrpern * pow (q10, exponent);
        dlmr_area_shade = n_area_shade * mrpern * pow (q10, exponent);
        /* finally, convert from mass to molar units, and from a daily rate to 
         * a rate per second */
        epv->dlmr_area_sun = dlmr_area_sun / (86400.0 * 12.011e-9);
        epv->dlmr_area_shade = dlmr_area_shade / (86400.0 * 12.011e-9);

        /* leaf, night */
        exponent = (tnight - 20.0) / 10.0;
        cf->leaf_night_mr = t1 * pow (q10, exponent) * (86400.0 - daily->dayl) / 86400.0;
    }
    else                        /* no leaves on */
    {
        cf->leaf_day_mr = 0.0;
        epv->dlmr_area_sun = 0.0;
        epv->dlmr_area_shade = 0.0;
        cf->leaf_night_mr = 0.0;
    }

    /* fine root maintenance respiration when fine roots on */
    /* ammended to consider only the specified n concentration,
     * to avoid excessive MR with n-loading to fine roots */
    if (cs->frootc)
    {
        exponent = (tsoil - 20.0) / 10.0;
        t1 = pow (q10, exponent);
        cf->froot_mr = ns->frootn * mrpern * t1;
    }
    else                        /* no fine roots on */
        cf->froot_mr = 0.0;

    /* TREE-specific fluxes */
    if (epc->woody)
    {
        /* live stem maintenance respiration */
        exponent = (tavg - 20.0) / 10.0;
        t1 = pow (q10, exponent);
        cf->livestem_mr = ns->livestemn * mrpern * t1;

        /* live coarse root maintenance respiration */
        exponent = (tsoil - 20.0) / 10.0;
        t1 = pow (q10, exponent);
        cf->livecroot_mr = ns->livecrootn * mrpern * t1;
    }

}
