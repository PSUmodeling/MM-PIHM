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

void SoilPsi (const soil_struct *soil, double vwc, double *psi)
{
    /*
     * Given a list of site constants and the soil water content,
     * this function returns the soil water potential (MPa)
     * Inputs:
     * vwc                  (kg/m2) water mass per unit area
     * sitec.soil_depth   (m)     effective soil depth
     * sitec.soil_b       (DIM)   slope of log(psi) vs log(rwc)
     * sitec.vwc_sat      (DIM)   volumetric water content at saturation
     * sitec.psi_sat      (MPa)   soil matric potential at saturation
     * output:
     * psi_s              (MPa)   soil matric potential
     *
     * uses the van Genuchten relation
     */

    double          theta;

    theta = (vwc - soil->smcmin) / (soil->smcmax - soil->smcmin);
    theta = (theta > 1.0) ? 1.0 : theta;
    theta = (theta < SATMIN) ? SATMIN : theta;

    /* calculate psi */
    *psi = Psi (theta, soil->alpha, soil->beta) * 1000.0 * GRAV / 1.0e6;
}
