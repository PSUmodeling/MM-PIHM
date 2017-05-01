#include "pihm.h"

void CanopyCond (const epconst_struct *epc, epvar_struct *epv,
    const eflux_struct *ef, const pstate_struct *ps,
    const daily_struct *daily)
{
    double          sw, sw_sun, sw_shade;
    double          ff, ff_sun, ff_shade;
    double          rcs, rcs_sun, rcs_shade;
    double          gl_s, gl_s_sun, gl_s_shade;
    double          gl_bl;

    sw = daily->avg_soldn;
    ff = 0.55 * 2.0 * sw / (epc->rgl * ps->proj_lai);
    rcs = (ff + epc->rsmin / epc->rsmax) / (1.0 + ff);
    rcs = (rcs > 0.0001) ? rcs : 0.0001;

    sw_sun = (ef->swabs_per_plaisun * ps->plaisun) /
        (ef->swabs_per_plaisun * ps->plaisun +
        ef->swabs_per_plaishade * ps->plaishade) * sw;
    ff_sun = 0.55 * 2.0 * sw_sun / (epc->rgl * ps->plaisun);
    rcs_sun = (ff_sun + epc->rsmin / epc->rsmax) / (1.0 + ff_sun);
    rcs_sun = (rcs_sun > 0.0001) ? rcs_sun : 0.0001;

    sw_shade = daily->avg_soldn - sw_sun;
    ff_shade = 0.55 * 2.0 * sw_shade / (epc->rgl * ps->plaisun);
    rcs_shade = (ff_shade + epc->rsmin / epc->rsmax) / (1.0 + ff_shade);
    rcs_shade = (rcs_shade > 0.0001) ? rcs_shade : 0.0001;

    /* Stomatal conductance per LAI calculated in Noah */
    gl_s = 1.0 / daily->avg_rc / ps->proj_lai;

    /* Adjust Noah stomatal conductance to calculate sunlit and shaded
     * conductance */
    gl_s_sun = rcs_sun / rcs * gl_s;
    gl_s_shade = rcs_shade / rcs * gl_s;

    gl_bl = 1.0 / daily->avg_ra;

    /* Leaf conductance to transpired water vapor, per unit projected LAI.
     * This formula is derived from stomatal and cuticular conductances in
     * parallel with each other, and both in series with leaf boundary layer
     * conductance. */
    epv->gl_t_wv_sun = gl_bl * gl_s_sun / (gl_bl + gl_s_sun);
    epv->gl_t_wv_shade = gl_bl * gl_s_shade / (gl_bl + gl_s_shade);
}
