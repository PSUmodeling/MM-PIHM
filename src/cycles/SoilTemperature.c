#include "Cycles.h"

/*****************************************************************************
 * FUNCTION NAME:   Temperature
 *
 * ARGUMENT LIST
 *
 * Argument             Type        IO  Description
 * ==========           ==========  ==  ====================
 * y                    int         I   Simulation year
 * doy                  int         I   Simulation day of year
 * snowCover            double      I
 * cropInterception     double      I
 * Soil                 SoilStruct  IO
 * Weather              WeatherStruct
 *                                  I   Weather structure
 * Residue              ResidueStruct
 *                                  I
 *
 * RETURN VALUE: void
 ****************************************************************************/
void Temperature (int y, int doy, double snowCover, double cropInterception, SoilStruct *Soil, WeatherStruct *Weather, ResidueStruct *Residue)
{
    /* LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * m                    int
     * i                    int         Loop counter
     * CP                   double[]
     * k                    double[]
     * CPsfc                double
     * ksfc                 double
     * a                    double[]
     * b                    double[]
     * c                    double[]
     * d                    double[]
     * T                    double[]
     * Tn                   double[]
     * tsfc                 double      Land surface temperature [C]
     * tAvg                 double      Average air temperature [C]
     * fCover               double      Soil cover factor accouting for flat
     *                                    residue and snow cover
     * soilCover            double      Soil fraction covered by snow and flat
     *                                    residues
     * counter              double      Loop counter
     * f                    double      Constant for setting forward, backward
     *                                    or centered difference method
     * g                    double
     */

    int             m;
    int             i;
    double          CP[Soil->totalLayers];
    double          k[Soil->totalLayers];
    double          CPsfc;
    double          ksfc;
    double          a[Soil->totalLayers + 1], b[Soil->totalLayers + 1], c[Soil->totalLayers + 1], d[Soil->totalLayers + 1];
    double          T[Soil->totalLayers + 1], Tn[Soil->totalLayers + 1];
    double          tsfc;
    double          tAvg;
    double          fCover;
    double          soilCover;
    int             counter;

    double          f = 0.6;
    double          g;

    g = 1.0 - f;

    m = Soil->totalLayers;

    CPsfc = 0.0;            /* Heat capacity of boundary layer  = 0 */
    ksfc = 20.0;            /* Boundary layer conductance (W/m2 K) */

    for (i = 0; i < m; i++)
    {
        /* Calculates heat capacity weighted by solid and water phase only  */
        CP[i] = HeatCapacity (Soil->BD[i], Soil->waterContent[i]) * Soil->layerThickness[i];
        /* Calculates conductance per layer, equation 4.20 for thermal
         * conductivity */
        k[i] = HeatConductivity (Soil->BD[i], Soil->waterContent[i], Soil->Clay[i]) / (Soil->nodeDepth[i + 1] - Soil->nodeDepth[i]);
    }

    /* Recalculates bottom boundary condition */
    T[Soil->totalLayers] = EstimatedSoilTemperature (Soil->nodeDepth[Soil->totalLayers], doy, Weather->annualAverageTemperature[y], Weather->yearlyAmplitude[y], Soil->annualTemperaturePhase, Soil->dampingDepth);
    Tn[Soil->totalLayers] = T[Soil->totalLayers];

    /* Passes previous time step temperatures for all soil layers */
    for (i = 0; i < Soil->totalLayers + 1; i++)
        T[i] = Soil->soilTemperature[i];

    /* Calculates temperature of upper boundary condition
     * uses an empirical factor to weight the effect of air temperature,
     * residue cover, and snow cover on the upper node temperature. This is an
     * empirical approach to allow for residue or snow insulation effect */
    tAvg = 0.5 * (Weather->tMax[y][doy - 1] + Weather->tMin[y][doy - 1]);
    soilCover = 1.0 - (1.0 - cropInterception) * (1.0 - snowCover) * Residue->flatResidueTau;
    fCover = 0.4 * soilCover + 0.3 * snowCover / (soilCover + 0.001);
    tsfc = (1.0 - fCover) * tAvg + fCover * T[0];

    counter = 0;

    do
    {
        /* This loop updates temperatures after the first Thomas loop and is
         * needed due to the inclusion of Tsurface [T(0)] to calculate net
         * radiation at "exchange" surface these two layers are used as
         * criteria to leave the loop (seen end of loop) */
        counter += 1;
        if (counter > 1)
        {
            T[0] = Tn[0];
            T[1] = Tn[1];
        }

        /* Calculates matrix elements */
        a[0] = 0.0;
        for (i = 0; i < Soil->totalLayers; i++)
        {
            c[i] = -k[i] * f;
            a[i + 1] = c[i];
            if (i == 0)
            {
                /* Changed to seconds per day, quite long time step */
                b[i] = f * (k[i] + ksfc) + CP[i] / 86400.0;
                d[i] = g * ksfc * tsfc + (CP[i] / 86400.0 - g * (k[i] + ksfc)) * T[i] + g * k[i] * T[i + 1];
            }
            else
            {
                /* Changed to seconds per day, quite long time step */
                b[i] = f * (k[i] + k[i - 1]) + CP[i] / 86400.0;
                d[i] = g * k[i - 1] * T[i - 1] + (CP[i] / 86400.0 - g * (k[i] + k[i - 1])) * T[i] + g * k[i] * T[i + 1];
            }
        }

        d[0] += ksfc * tsfc * f;
        d[Soil->totalLayers - 1] += k[Soil->totalLayers - 1] * f * Tn[Soil->totalLayers];

        /* Thomas algorithm starts */
        for (i = 0; i < Soil->totalLayers - 1; i++)
        {
            c[i] = c[i] / b[i];
            d[i] = d[i] / b[i];
            b[i + 1] = b[i + 1] - a[i + 1] * c[i];
            d[i + 1] = d[i + 1] - a[i + 1] * d[i];
        }

        Tn[Soil->totalLayers - 1] = d[Soil->totalLayers - 1] / b[Soil->totalLayers - 1];

        for (i = Soil->totalLayers - 2; i >= 0; i--)
            Tn[i] = d[i] - c[i] * Tn[i + 1];
        /* Thomas algorithm ends */
    } while (fabs (T[0] - Tn[0]) > 0.02 || fabs (T[1] - Tn[1]) > 0.02);

    for (i = 0; i < Soil->totalLayers + 1; i++)
        Soil->soilTemperature[i] = Tn[i];
}

double HeatCapacity (double bulkDensity, double volumetricWC)
{
    return (2400000.0 * bulkDensity / 2.65 + 4180000.0 * volumetricWC);
}

double HeatConductivity (double bulkDensity, double volumetricWC, double fractionClay)
{
    double          C1, C2, C3, C4; /* coeff to calculate thermal conductivity
                                     * (page 32+) */
    /* equation 4.27; coeff of 4.20 */
    C1 = 0.65 - 0.78 * bulkDensity + 0.6 * bulkDensity * bulkDensity;

    /* equation 4.25; coeff of 4.20 */
    C2 = 1.06 * bulkDensity;

    /* equation 4.28; coeff of 4.20 */
    C3 = 1.0 + 2.6 / sqrt (fractionClay);

    /* equation 4.22; coeff of 4.20 */
    C4 = 0.03 + 0.1 * bulkDensity * bulkDensity;

    return (C1 + C2 * volumetricWC - (C1 - C4) * exp (-pow (C3 * volumetricWC, 4)));
}

double EstimatedSoilTemperature (double nodeDepth, int doy, double annualAvgTemperature, double yearlyAmplitude, int phase, double dampingDepth)
{
    return (annualAvgTemperature + yearlyAmplitude * exp (-nodeDepth / dampingDepth) * sin (2.0 * PI / 365.0 * (double)(doy - phase) - nodeDepth / dampingDepth));
}
