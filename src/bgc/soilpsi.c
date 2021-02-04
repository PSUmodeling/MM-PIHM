#include "pihm.h"

void SoilPsi(const soil_struct *soil, double vwc, double *psi)
{
    // Given a list of site constants and the soil water content,
    // this function returns the soil water potential (MPa)
    // Inputs:
    // vwc                (m3/m3) water mass per unit area
    // psi_s              (MPa)   soil matric potential
    //
    // Uses the van Genuchten relation
    double          theta;

    theta = (vwc - soil->smcmin) / (soil->smcmax - soil->smcmin);
    theta = MIN(theta, 1.0);
    theta = MAX(theta, SATMIN);

    // Calculate psi
    *psi = Psi(theta, soil->alpha, soil->beta) * 1000.0 * GRAV / 1.0E6;
}
