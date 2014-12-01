/* 
 * canopy_et.c
 * A single-function treatment of canopy evaporation and transpiration
 * fluxes.  
 * 
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 * Biome-BGC version 4.2 (final release)
 * See copyright.txt for Copyright information
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 */

#include "bgc.h"

void canopy_et (const metvar_struct *metv, const epconst_struct *epc, epvar_struct * epv, wflux_struct *wf)
{
    double          gl_bl, gl_c, gl_s_sun, gl_s_shade;
    double          gl_e_wv, gl_t_wv_sun, gl_t_wv_shade, gl_sh;
    double          gc_e_wv, gc_sh;
    double          tday;
    double          tmin;
    double          dayl;
    double          vpd, vpd_open, vpd_close;
    double          psi, psi_open, psi_close;
    double          m_ppfd_sun, m_ppfd_shade;
    double          m_tmin, m_psi, m_co2, m_vpd, m_final_sun, m_final_shade;
    double          proj_lai;
    double          canopy_w;
    double          gcorr;

    double          e, cwe, t, trans, trans_sun, trans_shade, e_dayl, t_dayl;
    pmet_struct     pmet_in;

    /* assign variables that are used more than once */
    tday = metv->tday;
    tmin = metv->tmin;
    vpd = metv->vpd;
    dayl = metv->dayl;
    psi = epv->psi;
    proj_lai = epv->proj_lai;
    canopy_w = wf->prcp_to_canopyw;
    psi_open = epc->psi_open;
    psi_close = epc->psi_close;
    vpd_open = epc->vpd_open;
    vpd_close = epc->vpd_close;

    /* temperature and pressure correction factor for conductances */
    gcorr = pow ((metv->tday + 273.15) / 293.15, 1.75) * 101300. / metv->pa;

    /* calculate leaf- and canopy-level conductances to water vapor
     * and sensible heat fluxes */

    /* leaf boundary-layer conductance */
    gl_bl = epc->gl_bl * gcorr;

    /* leaf cuticular conductance */
    gl_c = epc->gl_c * gcorr;

    /* leaf stomatal conductance: first generate multipliers, then apply them
     * to maximum stomatal conductance */
    /* calculate stomatal conductance radiation multiplier:
     * *** NOTE CHANGE FROM BIOME-BGC CODE ***
     * The original Biome-BGC formulation follows the arguments in
     * Rastetter, E.B., A.W. King, B.J. Cosby, G.M. Hornberger, R.V. O'Neill,
     * and J.E. Hobbie, 1992. Aggregating fine-scale ecological knowledge to
     * model coarser-scale attributes of ecosystems. Ecological Applications,
     * 2:55-70.
     *
     * gmult->max = (gsmax/(k*lai))*log((gsmax+rad)/(gsmax+(rad*exp(-k*lai))))
     *
     * I'm using a much simplified form, which doesn't change relative shape
     * as gsmax changes. See Korner, 1995.*/

    /* photosynthetic photon flux density conductance control */
    m_ppfd_sun = metv->ppfd_per_plaisun / (PPFD50 + metv->ppfd_per_plaisun);
    m_ppfd_shade = metv->ppfd_per_plaishade / (PPFD50 + metv->ppfd_per_plaishade);

    /* soil-leaf water potential multiplier */
    //printf ("%lf %lf %lf\n", psi, psi_open, psi_close);
    if (psi > psi_open)         /* no water stress */
        m_psi = 1.0;
    else if (psi <= psi_close)  /* full water stress */
        m_psi = 0.0;
    else                        /* partial water stress */
        m_psi = (psi_close - psi) / (psi_close - psi_open);

    /* CO2 multiplier */
    m_co2 = 1.0;

    /* freezing night minimum temperature multiplier */
    if (tmin > 0.0)             /* no effect */
        m_tmin = 1.0;
    else if (tmin < -8.0)       /* full tmin effect */
        m_tmin = 0.0;
    else                        /* partial reduction (0.0 to -8.0 C) */
        m_tmin = 1.0 + (0.125 * tmin);

    /* vapor pressure deficit multiplier, vpd in Pa */
    if (vpd < vpd_open)         /* no vpd effect */
        m_vpd = 1.0;
    else if (vpd > vpd_close)   /* full vpd effect */
        m_vpd = 0.0;
    else                        /* partial vpd effect */
        m_vpd = (vpd_close - vpd) / (vpd_close - vpd_open);

    //printf ("m_psi %lf m_co2 %lf m_tmin %lf m_vpd %lf\n", m_psi, m_co2, m_tmin, m_vpd);
    /* apply all multipliers to the maximum stomatal conductance */
    m_final_sun = m_ppfd_sun * m_psi * m_co2 * m_tmin * m_vpd;
    if (m_final_sun < 0.00000001)
        m_final_sun = 0.00000001;
    m_final_shade = m_ppfd_shade * m_psi * m_co2 * m_tmin * m_vpd;
    if (m_final_shade < 0.00000001)
        m_final_shade = 0.00000001;
    //printf ("m_final_sun = %lf, m_final_shad = %lf\n", m_final_sun, m_final_shade);
    gl_s_sun = epc->gl_smax * m_final_sun * gcorr;
    gl_s_shade = epc->gl_smax * m_final_shade * gcorr;

    /* calculate leaf-and canopy-level conductances to water vapor and
     * sensible heat fluxes, to be used in Penman-Monteith calculations of
     * canopy evaporation and canopy transpiration. */

    /* Leaf conductance to evaporated water vapor, per unit projected LAI */
    gl_e_wv = gl_bl;

    /* Leaf conductance to transpired water vapor, per unit projected LAI.
     * This formula is derived from stomatal and cuticular conductances in
     * parallel with each other, and both in series with leaf boundary layer
     * conductance. */
    gl_t_wv_sun = (gl_bl * (gl_s_sun + gl_c)) / (gl_bl + gl_s_sun + gl_c);
    gl_t_wv_shade = (gl_bl * (gl_s_shade + gl_c)) / (gl_bl + gl_s_shade + gl_c);

    /* assign leaf-level conductance to transpired water vapor, 
     * for use in calculating co2 conductance for farq_psn() */
    epv->gl_t_wv_sun = gl_t_wv_sun;
    epv->gl_t_wv_shade = gl_t_wv_shade;
}
