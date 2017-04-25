#include "pihm.h"

void CanopyCond (const epconst_struct *epc, epvar_struct *epv,
    eflux_struct *ef, pstate_struct *ps)
{
    double          ff_sun, ff_shade;
    double          rcs_sun, rcs_shade;
    double          gl_s_sun, gl_s_shade;
    double          gl_bl;
    double          ratio;

    ff_sun =
        0.55 * ef->swabs_per_plaisun / (1.0 - ps->albedo) / epc->rgl * 2.0;
    rcs_sun = (ff_sun + epc->rsmin / epc->rsmax) / (1.0 + ff_sun);
    rcs_sun = (rcs_sun > 0.0001) ? rcs_sun : 0.0001;

    ff_shade =
        0.55 * ef->swabs_per_plaishade / (1.0 - ps->albedo) / epc->rgl * 2.0;
    rcs_shade = (ff_shade + epc->rsmin / epc->rsmax) / (1.0 + ff_shade);
    rcs_shade = (rcs_shade > 0.0001) ? rcs_shade : 0.0001;

    gl_s_sun = rcs_sun / ps->rcs / ps->rc;

    gl_s_shade = rcs_shade / ps->rcs / ps->rc;

    gl_bl = ps->ch / ps->proj_lai;

    /* Leaf conductance to transpired water vapor, per unit projected LAI.
     * This formula is derived from stomatal and cuticular conductances in
     * parallel with each other, and both in series with leaf boundary layer
     * conductance. */
    epv->gl_t_wv_sun = gl_bl * gl_s_sun / (gl_bl + gl_s_sun);
    epv->gl_t_wv_shade = gl_bl * gl_s_shade / (gl_bl + gl_s_shade);
}
