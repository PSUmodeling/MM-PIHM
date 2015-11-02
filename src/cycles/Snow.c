#include "Cycles.h"

void SnowProcesses (SnowStruct *Snow, int y, int doy, WeatherStruct *Weather, double TauStandingRes, double CropInterception)
{
    /*
     * 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * PP		    double	Daily precipitation
     * Tavg		    double	Daily average air temperature
     * Tx		    double	Daily maximum air temperature
     * Tn		    double	Daily minimum air temperature
     */
    double          PP;
    double	    Tavg;
    double	    Tx;
    double	    Tn;

    Snow->snowFall = 0.0;
    Snow->snowMelt = 0.0;
    Snow->snowCover = 0.0;
    Snow->snowEvaporationVol = 0.0;

    PP = Weather->precipitation[y][doy - 1];
    Tx = Weather->tMax[y][doy - 1];
    Tn = Weather->tMin[y][doy - 1];
    Tavg = 0.5 * (Tx + Tn);

    CalculateSnowFall (Snow, Tavg, PP);

    if (Snow->Snow > 0.0)
    {
        CalculateSnowMelt (Snow, Tavg, Tx, Tn);
        if (Snow->Snow > 0.0)
            CalculateSnowEvaporation (Snow, TauStandingRes, CropInterception, Weather->ETref[y][doy - 1]);
        Snow->snowCover = CalculateSnowCover (Snow);
    }
}

void CalculateSnowFall (SnowStruct *Snow, double Tavg, double PP)
{
    if (PP > 0.0 && Tavg < THRESHOLD_TEMPERATURE_SNOWFALL)
    {
        Snow->snowFall = PP;
        Snow->Snow += Snow->snowFall;
    }
}

void CalculateSnowMelt (SnowStruct *Snow, double Tavg, double Tx, double Tn)
{
    /*
     * 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * TTmelt		    double
     */
    double          TTmelt;

    if (Tn > THRESHOLD_TEMPERATURE_SNOWMELT)
        TTmelt = Tavg - THRESHOLD_TEMPERATURE_SNOWMELT;
    else if (Tx < THRESHOLD_TEMPERATURE_SNOWMELT)
        TTmelt = 0.0;
    else
        TTmelt = pow (Tx - THRESHOLD_TEMPERATURE_SNOWMELT, 2.0) / (Tx - Tn);

    Snow->snowMelt = TTmelt * SNOWMELT_RATE;
    if (Snow->snowMelt > Snow->Snow)
        Snow->snowMelt = Snow->Snow;
    Snow->Snow = Snow->Snow - Snow->snowMelt;
}

void CalculateSnowEvaporation (SnowStruct *Snow, double TauStandingRes, double CropInterception, double ETo)
{
    /*
     * 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * evaporativeDemand    double
     */
    double          evaporativeDemand;

    evaporativeDemand = TauStandingRes * (1.0 - CropInterception) * ETo;

    if (evaporativeDemand > Snow->Snow)
        Snow->snowEvaporationVol = Snow->Snow;
    else
        Snow->snowEvaporationVol = evaporativeDemand;

    Snow->Snow = Snow->Snow - Snow->snowEvaporationVol;
}

double CalculateSnowCover (SnowStruct *Snow)
{
    /*
     * 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * snow_cover	    double	[return value]
     */
    double          snow_cover;

    snow_cover = (1.0 - exp (-0.43 * pow (Snow->Snow, 1.14)));

    return (snow_cover);
}
