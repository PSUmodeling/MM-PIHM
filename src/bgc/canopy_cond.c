#include "pihm.h"

void CanopyCond(const soil_struct *soil, const epconst_struct *epc, const daily_struct *daily, const phystate_struct *ps, const eflux_struct *ef, epvar_struct *epv)
{
    int             kz;
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

    sw_sun = (ef->swabs_per_plaisun * ps->plaisun) / (ef->swabs_per_plaisun * ps->plaisun + ef->swabs_per_plaishade * ps->plaishade) * daily->avg_soldn;
    ff_sun = 0.55 * 2.0 * sw_sun / (epc->rgl * ps->plaisun);
    rcs_sun = (ff_sun + epc->rsmin / epc->rsmax) / (1.0 + ff_sun);
    rcs_sun = MAX(rcs_sun, 1.0E-4);

    sw_shade = daily->avg_soldn - sw_sun;
    ff_shade = 0.55 * 2.0 * sw_shade / (epc->rgl * ps->plaishade);
    rcs_shade = (ff_shade + epc->rsmin / epc->rsmax) / (1.0 + ff_shade);
    rcs_shade = MAX(rcs_shade, 1.0E-4);

    rct = 1.0 - 0.0016 * pow(epc->topt - daily->tday, 2.0);
    rct = MAX(rct, 1.0E-4);

    rcq = 1.0 / (1.0 + epc->hs * daily->avg_q2d);
    rcq = MAX(rcq, 0.01);

    gx = (daily->avg_sh2o[0] - soil->smcwlt) / (soil->smcref - soil->smcwlt);
    gx = MIN(gx, 1.0);
    gx = MAX(gx, 0.0);

    part[0] = (ps->zsoil[0] / ps->zsoil[ps->nroot - 1]) * gx;
    for (kz = 1; kz < ps->nroot; kz++)
    {
        gx = (daily->avg_sh2o[kz] - soil->smcwlt) / (soil->smcref - soil->smcwlt);
        gx = MIN(gx, 1.0);
        gx = MAX(gx, 0.0);

        part[kz] = ((ps->zsoil[kz] - ps->zsoil[kz - 1]) / ps->zsoil[ps->nroot - 1]) * gx;
    }

    for (kz = 0; kz < ps->nroot; kz++)
    {
        rcsoil += part[kz];
    }
    rcsoil = MAX(rcsoil, 1.0E-4);

    // Adjust Noah stomatal conductance to calculate sunlit and shaded conductance
    gl_s_sun = rcs_sun * rct * rcq * rcsoil / epc->rsmin;
    gl_s_shade = rcs_shade * rct * rcq * rcsoil / epc->rsmin;

    gl_bl = daily->avg_ch;

    // Leaf conductance to transpired water vapor, per unit projected LAI. This formula is derived from stomatal and
    // cuticular conductances in parallel with each other, and both in series with leaf boundary layer conductance.
    epv->gl_t_wv_sun = gl_bl * gl_s_sun / (gl_bl + gl_s_sun);
    epv->gl_t_wv_shade = gl_bl * gl_s_shade / (gl_bl + gl_s_shade);
}
