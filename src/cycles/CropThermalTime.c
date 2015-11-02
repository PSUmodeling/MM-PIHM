#include "Cycles.h"

void ComputeThermalTime (int total_years, CommunityStruct *Community, WeatherStruct *Weather)
{
    /* 
     * Calculate flowering and Maturity thermal time for each crop
     * in the rotation
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * sumTT		    double
     * sumMTTbyYear	    double
     * sumFTTbyYear	    double
     * cropEvents	    int
     * cropON		    int
     * c		    int
     * y		    int
     * d		    int
     */

    double          sumTT = 0.0;
    double          sumMTTbyYear = 0.0;
    double          sumFTTbyYear = 0.0;

    int             cropEvents;
    int             cropON;
    int             c;
    int		    y;
    int		    d;

    CropStruct *Crop;

    //SelectCropInitialPosition (CropManagement);

    for (c = 0; c < Community->NumCrop; c++)
    {
        Crop = &(Community->Crop[c]);

        if (Crop->calculatedMaturityTT == 0.0)
        {
            /* Computes thermal time for 1st instance of crop in each
             * rotation */
            cropEvents = 0;
            cropON = 0;

            for (y = 0; y < total_years; y++)
            {
                for (d = 1; d <= Weather->lastDoy[y]; d++)
                {
                    if (d == Crop->userSeedingDate)
                        cropON = 1;

                    if (cropON)
                    {
                        sumTT = sumTT + 0.5 * (ThermalTime (Crop->userTemperatureBase, Crop->userTemperatureOptimum, Crop->userTemperatureMaximum, Weather->tMax[y][d - 1]) + ThermalTime (Crop->userTemperatureBase, Crop->userTemperatureOptimum, Crop->userTemperatureMaximum, Weather->tMin[y][d - 1]));

                        if (d == Crop->userFloweringDate)
                            sumFTTbyYear = sumFTTbyYear + sumTT;

                        if (d == Crop->userMaturityDate)
                        {
                            cropEvents = cropEvents + 1;
                            sumMTTbyYear = sumMTTbyYear + sumTT;
                            sumTT = 0.0;
                            cropON = 0;
                        }
                    }
                }
            }

            /* Load flow and mat TT as crop info for the crops in the rotation
             * Order of crops in array are assigned according to order in the
             * rotation */
            if (cropEvents > 0)
            {
                Crop->calculatedFloweringTT = sumFTTbyYear / cropEvents;
                Crop->calculatedMaturityTT = sumMTTbyYear / cropEvents;
            }
            else
            {
                Crop->calculatedFloweringTT = 0.0;
                Crop->calculatedMaturityTT = 0.0;
            }
        }
        sumFTTbyYear = 0.0;
        sumMTTbyYear = 0.0;
    }
}

double ThermalTime (double T_base, double T_op, double T_Max, double Temperature)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * thermal_time	    double	[return value]
     */
    double          thermal_time;

    if (Temperature <= T_base || Temperature >= T_Max)
        thermal_time = 0.0;
    else if (Temperature < T_op)
        thermal_time = Temperature - T_base;
    else
        thermal_time = (T_Max - Temperature) / (T_Max - T_op) * (T_op - T_base);

    return (thermal_time);
}
