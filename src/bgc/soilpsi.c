/* 
 * soilpsi.c
 * soil water potential as a function of volumetric water content and
 * constants related to texture
 * 
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 * Biome-BGC version 4.2 (final release)
 * See copyright.txt for Copyright information
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 */

#include "pihm.h"

void soilpsi (const siteconst_struct * sitec, double vwc, double *psi)
{
    /*
     * Given a list of site constants and the soil water content,
     * this function returns the soil water potential (MPa)
     * Inputs:
     * ws.soilw           (kg/m2) water mass per unit area
     * sitec.soil_depth   (m)     effective soil depth
     * sitec.soil_b       (DIM)   slope of log(psi) vs log(rwc)
     * sitec.vwc_sat      (DIM)   volumetric water content at saturation
     * sitec.psi_sat      (MPa)   soil matric potential at saturation
     * output:
     * psi_s              (MPa)   soil matric potential
     *
     * uses the van Genuchten relation:
     * psi_s = psi_sat * (vwc/vwc_sat)^b
     */

    double          theta;
    double          alpha, beta;

    theta = (vwc - sitec->vwc_min) / (sitec->vwc_sat - sitec->vwc_min);
    theta = theta > 1. ? 1. : theta;
    theta = theta < 0. ? 0.05 : theta;
    alpha = sitec->soil_alpha;
    beta = sitec->soil_beta;

    /* calculate psi */
    *psi = -(pow (pow (1. / theta, sitec->soil_beta / (sitec->soil_beta - 1.)) - 1., 1. / sitec->soil_beta) / sitec->soil_alpha) * 1000. * 9.81 / 1.e6;
    //printf ("theta = %lf, alpha = %lf, beta = %lf, psi = %lf\n", theta, sitec->soil_alpha, sitec->soil_beta, *psi);
}
