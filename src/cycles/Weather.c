#include "Cycles.h"

/*****************************************************************************
 * FUNCTION NAME: SatVP
 *
 * ARGUMENT LIST
 *
 * Argument             Type        IO  Description
 * ==========           ==========  ==  ====================
 * T                    double      I   Air temperature
 *
 * RETURN VALUE: double (saturated vapor pressur [kPa])
 ****************************************************************************/

double SatVP (double T)
{
    return 0.6108 * exp (17.27 * T / (T + 237.3));
}

/*****************************************************************************
 * FUNCTION NAME: CalculateDerivedWeather
 *
 * ARGUMENT LIST
 *
 * Argument             Type        IO  Description
 * ==========           ==========  ==  ====================
 * Weather              WeatherStruct
 *                                  IO
 * total_years          int         I
 *
 * RETURN VALUE: void
 ****************************************************************************/
void CalculateDerivedWeather (WeatherStruct *Weather, int total_years)
{
    /* LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * y                    int         Year
     * doy                  int         Day of year
     * annualMaxTemperatureSum
     *                      double
     * annualMinTemperatureSum
     *                      double
     */
    int             y, doy;
    double          annualMaxTemperatureSum;
    double          annualMinTemperatureSum;

    Weather->atmosphericPressure = 101.325 * exp (-Weather->siteAltitude / 8200.0);

    for (y = 0; y < total_years; y++)
    {
        annualMaxTemperatureSum = 0.0;
        annualMinTemperatureSum = 0.0;

        for (doy = 1; doy < Weather->lastDoy[y] + 1; doy++)
        {
            Weather->ETref[y][doy - 1] = CalculatePMET (Weather->siteLatitude, Weather->atmosphericPressure, Weather->screeningHeight, Weather->tMax[y][doy - 1], Weather->tMin[y][doy - 1], Weather->solarRadiation[y][doy - 1], Weather->RHmax[y][doy - 1], Weather->RHmin[y][doy - 1], Weather->wind[y][doy - 1], doy);
            annualMaxTemperatureSum = annualMaxTemperatureSum + Weather->tMax[y][doy - 1];
            annualMinTemperatureSum = annualMinTemperatureSum + Weather->tMin[y][doy - 1];
        }
        Weather->annualAverageTemperature[y] = (annualMaxTemperatureSum + annualMinTemperatureSum) / ((double)(Weather->lastDoy[y])) / 2.0;
        Weather->yearlyAmplitude[y] = (annualMaxTemperatureSum - annualMinTemperatureSum) / ((double)(Weather->lastDoy[y]));
    }
}
