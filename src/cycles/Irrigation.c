#include "Cycles.h"

double FindIrrigationVolume (int opLayer, double opWaterDepletion, const SoilStruct *Soil)
{
    /*
     * This sub calculates the volume of an automatic irrigation
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * i		    int		Loop counter
     * depletionZonePAW	    double
     * depletionZoneWater   double
     * irrigation_vol	    double	Volume of an automatic irrigation
     *					  [return value]
     */
    int             i;
    double          depletionZonePAW = 0.0;
    double          depletionZoneWater = 0.0;
    double          irrigation_vol = 0.0;

    if (opLayer > Soil->totalLayers)
        opLayer = Soil->totalLayers;

    for (i = 0; i < opLayer; i++)
    {
        if (Soil->waterContent[i] >= Soil->PWP[i])
        {
            depletionZonePAW += Soil->PAW[i] * Soil->layerThickness[i];
            depletionZoneWater += (Soil->waterContent[i] - Soil->PWP[i]) * Soil->layerThickness[i];
        }
        else
	{
	    /* This adds deficit below PWP */
            depletionZonePAW += (Soil->FC[i] - Soil->waterContent[i]) * Soil->layerThickness[i];
	}
    }

    if (depletionZoneWater + Soil->irrigationVol < depletionZonePAW * (1.0 - opWaterDepletion))
    {
        /* Convert m to mm */
        irrigation_vol = (depletionZonePAW - depletionZoneWater) * 1000.0 - Soil->irrigationVol;
    }

    return (irrigation_vol);
}
