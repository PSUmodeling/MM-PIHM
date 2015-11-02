#include "Cycles.h"

void SoluteTransport (int totalLayers, double Sol_Kd, double WInConc, double *leachate, const double *WFlux, double *soluteMass, const double *BD, const double *thickness, const double *porosity, const double *WCinitial)
{
    /*
     * 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * i                    int         Loop counter
     * soluteConc           double      Solute concentration, current layer
     *                                    [kg solute / m3]
     * soluteFlow           double[]    Solute flow [kg solute/time step]
     * waterInitial         double      Initial layer soil water storage
     *                                    [kg/m2]
     * ratio                double      Ratio of flux out of a layer to the
     *                                    layer porosity (or effective
     *                                    porosity)
     * soluteMassAdsorbed   double      [kg/m2]
     * soluteMassSolution   double      [kg/m2]
     */

    int             i;
    double          soluteConc;
    double          soluteFlow[totalLayers + 1];
    double          waterInitial;
    double          ratio;
    double          soluteMassAdsorbed;
    double          soluteMassSolution;

    soluteFlow[0] = WInConc * WFlux[0];

    for (i = 0; i < totalLayers; i++)
    {
        soluteConc = 0.0;
        soluteMassAdsorbed = 0.0;
        waterInitial = thickness[i] * WATER_DENSITY * WCinitial[i];

        if (soluteMass[i] > 0.0)
            soluteConc = LinearEquilibriumConcentration (Sol_Kd, BD[i], thickness[i], WCinitial[i], soluteMass[i]);

        soluteMassAdsorbed = soluteMass[i] - waterInitial * soluteConc;

        if (WFlux[i + 1] > 0.0)
        {
            /* Note: the 0.67 is temporary, likely not needed, and if needed
             * use a relationship with B */
            ratio = 0.67 * WFlux[i + 1] / (porosity[i] * thickness[i] * WATER_DENSITY);
            soluteMassSolution = (soluteFlow[i] / ratio + (waterInitial * soluteConc - soluteFlow[i] / ratio) * exp (-ratio));
            soluteFlow[i + 1] = waterInitial * soluteConc + soluteFlow[i] - soluteMassSolution;
            if (i == totalLayers - 1)
                *leachate = soluteFlow[i + 1];
        }
        else
        {
            soluteMassSolution = waterInitial * soluteConc + soluteFlow[i];
            soluteFlow[i + 1] = 0.0;
        }

        /* This equilibrium can cause a mass balance error if "mass" is based
         * on concentration at the end ignoring initial mass, not using it */
        //SoluteMass(i) = LinearEquilibriumSoluteMass(Sol_Kd, BD(i), Thickness(i), WCFinal, C(i))
        soluteMass[i] = soluteMassAdsorbed + soluteMassSolution;
    }
}

void SoluteTransportEvaporation (int totalLayers, double Sol_Kd, const double *WFlux, double *soluteMass, const double *BD, const double *thickness, const double *porosity, double *WCinitial)
{
    /*
     * This is a rough attempt to move solutes upward due to evaporation
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * i                    int
     * soluteConc           double      Solute concentration, current layer
     *                                    [kg solute/m3]
     * waterInitial         double      Initial layer soil water storage
     *                                    [kg/m2]
     * soluteFlux           double[]
     */

    int             i;
    double          soluteConc;
    double          waterInitial;
    double          soluteFlux[totalLayers];

    for (i = 1; i < totalLayers; i++)
    {
        waterInitial = thickness[i] * WATER_DENSITY * WCinitial[i];

        if (soluteMass[i] > 0.0)
            soluteConc = LinearEquilibriumConcentration (Sol_Kd, BD[i], thickness[i], WCinitial[i], soluteMass[i]);

        else
            soluteConc = 0.0;

        soluteFlux[i] = WFlux[i] * soluteConc;
    }

    for (i = totalLayers - 2; i >= 0; i--)
    {
        if (soluteFlux[i + 1] > 0)
        {
            soluteMass[i + 1] -= soluteFlux[i + 1];
            soluteMass[i] += soluteFlux[i + 1];
        }
    }
}

/*****************************************************************************
 * FUNCTION NAME:   LinearEquilibriumConcentration
 *
 * ARGUMENT LIST
 *
 * Argument             Type        IO  Description
 * ==========           ==========  ==  ====================
 * Kd                   double      I   Slope of the adsortion isotherm
 *                                        [m3/kg]
 * bulkDensity          double      I   Soil bulk density [kg/m3]
 * layerThickness       double      I
 * waterContent         double      I
 * soluteMass           double      I
 *
 * RETURN VALUE: double (solute equilibrium concentration [kg solute/m3])
 ****************************************************************************/
double LinearEquilibriumConcentration (double Kd, double bulkDensity, double layerThickness, double waterContent, double soluteMass)
{
    /* LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * soilBufferPower      double      [m3/m3]
     */
    double          soilBufferPower;

    bulkDensity *= 1000.0;      /* Convert Mg/m3 to kg/m3 */
    soilBufferPower = Kd * bulkDensity + waterContent;

    return (soluteMass / (soilBufferPower * layerThickness * WATER_DENSITY));
}

/*****************************************************************************
 * FUNCTION NAME:   LinearEquilibriumSoluteMass
 *
 * ARGUMENT LIST
 *
 * Argument             Type        IO  Description
 * ==========           ==========  ==  ====================
 * Kd                   double      I   Slope of the adsortion isotherm
 *                                        [m3/kg]
 * bulkDensity          double      I   Soil bulk density [kg/m3]
 * layerThickness       double      I
 * waterContent         double      I
 * concentration        double      I   Solute equilibrium concentration
 *                                        [kg solute / m3]
 *
 * RETURN VALUE: double (total solute mass based on concentration in the
 *   solution)
 ****************************************************************************/
double LinearEquilibriumSoluteMass (double Kd, double bulkDensity, double layerThickness, double waterContent, double concentration)
{
    /* LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * soilBufferPower      double      [m3/m3]
     */

    double          soilBufferPower;

    bulkDensity *= 1000.0;      /* convert Mg/m3 to kg/m3 */
    soilBufferPower = Kd * bulkDensity + waterContent;

    return (concentration * (soilBufferPower * layerThickness * WATER_DENSITY));
}
