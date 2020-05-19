#include "pihm.h"

void CanopyCond(const epconst_struct *epc, epvar_struct *epv,
    const eflux_struct *ef, const phystate_struct *ps, const soil_struct *soil,
    const daily_struct *daily)
{
    int             k;
    double          part[MAXLYR];
    double          gx;
    double          rct = 0.0;
    double          rcq = 0.0;
    double          rcsoil = 0.0;
    double          sw_sun, sw_shade;
    double          ff_sun, ff_shade;
    double          rcs_sun, rcs_shade;
    double          gl_s_sun, gl_s_shade;
    double          gl_bl;

    sw_sun = (ef->swabs_per_plaisun * ps->plaisun) /
        (ef->swabs_per_plaisun * ps->plaisun +
        ef->swabs_per_plaishade * ps->plaishade) * daily->avg_soldn;
    ff_sun = 0.55 * 2.0 * sw_sun / (epc->rgl * ps->plaisun);
    rcs_sun = (ff_sun + epc->rsmin / epc->rsmax) / (1.0 + ff_sun);
    rcs_sun = (rcs_sun > 0.0001) ? rcs_sun : 0.0001;

    sw_shade = daily->avg_soldn - sw_sun;
    ff_shade = 0.55 * 2.0 * sw_shade / (epc->rgl * ps->plaishade);
    rcs_shade = (ff_shade + epc->rsmin / epc->rsmax) / (1.0 + ff_shade);
    rcs_shade = (rcs_shade > 0.0001) ? rcs_shade : 0.0001;

    rct = 1.0 - 0.0016 * pow(epc->topt - daily->tday, 2.0);
    rct = (rct > 0.0001) ? rct : 0.0001;

    rcq = 1.0 / (1.0 + epc->hs * daily->avg_q2d);
    rcq = (rcq > 0.01) ? rcq : 0.01;

    gx = (daily->avg_sh2o[0] - soil->smcwlt) / (soil->smcref - soil->smcwlt);
    gx = (gx > 1.0) ? 1.0 : gx;
    gx = (gx < 0.0) ? 0.0 : gx;

    part[0] = (ps->zsoil[0] / ps->zsoil[ps->nroot - 1]) * gx;
    for (k = 1; k < ps->nroot; k++)
    {
        gx = (daily->avg_sh2o[k] - soil->smcwlt) /
            (soil->smcref - soil->smcwlt);
        gx = (gx > 1.0) ? 1.0 : gx;
        gx = (gx < 0.0) ? 0.0 : gx;

        part[k] = ((ps->zsoil[k] - ps->zsoil[k - 1]) /
            ps->zsoil[ps->nroot - 1]) * gx;
    }

    for (k = 0; k < ps->nroot; k++)
    {
        rcsoil += part[k];
    }
    rcsoil = (rcsoil > 0.0001) ? rcsoil : 0.0001;

    /* Adjust Noah stomatal conductance to calculate sunlit and shaded
     * conductance */
    gl_s_sun = rcs_sun * rct * rcq * rcsoil / epc->rsmin;
    gl_s_shade = rcs_shade * rct * rcq * rcsoil / epc->rsmin;

    gl_bl = daily->avg_ch;

    /* Leaf conductance to transpired water vapor, per unit projected LAI.
     * This formula is derived from stomatal and cuticular conductances in
     * parallel with each other, and both in series with leaf boundary layer
     * conductance. */
    epv->gl_t_wv_sun = gl_bl * gl_s_sun / (gl_bl + gl_s_sun);
    epv->gl_t_wv_shade = gl_bl * gl_s_shade / (gl_bl + gl_s_shade);
}
