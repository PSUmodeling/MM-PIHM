#include "Cycles.h"

void WaterUptake (int y, int doy, CommunityStruct *Community, SoilStruct *Soil, const WeatherStruct *Weather)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * i		    int
     * PTx		    double	Maximum full cover transpiration
     *					  (mm/day)
     * PT		    double	Potential transpiration (mm/day)
     * TE		    double	Expected transpiration, the minimum of
     *					  PT and PWU (mm/day)
     * TA		    double	Attainable transpiration (or TE
     *					  adjusted by LWP) (mm/day)
     * transpirationRatio   double
     * temperatureAvg	    double
     * factorTemperature    double
     * platHC		    double	Plant hydraulic conductance (kg2/J/m2)
     * rootHC		    double	Root hydraulic conductance
     * shootHC		    double	Shoot hydraulic conductance
     * Root_HC_Adjustment   double	Dryness adjustment of hydraulic
     *					  conductance
     * layerPlantHC	    double[]
     * layerRootHC	    double[]
     * layerShootHC	    double[]
     * Layer_Root_HC_Adjustment
     *			    double[]
     * rootActivity	    double[]	Root activity factor to compute plant
     *					  hydraulic conductance
     * rootFraction	    double[]
     * LWP		    double	Leaf water potential (J/kg)
     * LWP_StressOnset	    double	Leaf water potential at the onset of
     *					  stomatal closure
     * LWP_WiltingPoint	    double	Leaf water potential at wilting point
     * SWP_FC		    double	Water potential at field capacity
     *					  (J/kg)
     * SWP_Average	    double	Weighted soil water potential (J/kg)
     * soilWP		    double[]
     * layerSalinityFactor  double[]
     */
    int             i, j;

    //double          PTx = 15.0;
    double          PTx;
    double          PT = 0.00001;
    double          TE;
    double          TA;

    double          transpirationRatio;
    double          temperatureAvg;
    double          factorTemperature;

    double          plantHC;
    double          rootHC;
    double          shootHC;
    double          Root_HC_Adjustment;
    double          layerPlantHC[Soil->totalLayers];
    double          layerRootHC[Soil->totalLayers];
    double          layerShootHC[Soil->totalLayers];
    double          Layer_Root_HC_Adjustment[Soil->totalLayers];
    double          rootActivity[Soil->totalLayers];
    double          rootFraction[Soil->totalLayers];
    double          waterUptake[Soil->totalLayers];

    double          LWP;
    //double          LWP_StressOnset = -1100.0;
    //double          LWP_WiltingPoint = -2000.0;
    double          SWP_FC = -33.0;
    double          SWP_Average;
    double          soilWP[Soil->totalLayers];
    double          layerSalinityFactor[Soil->totalLayers];
    CropStruct     *Crop;

    for (i = 0; i < Soil->totalLayers; i++)
    {
        Soil->waterUptake[i] = 0.0;
        soilWP[i] = SoilWaterPotential (Soil->Porosity[i], Soil->airEntryPotential[i], Soil->B_Value[i], Soil->waterContent[i]);
    }


    for (j = 0; j < Community->NumCrop; j++)
    {
        Crop = &Community->Crop[j];

        for (i = 0; i < Soil->totalLayers; i++)
            waterUptake[i] = 0.0;

        if (Crop->stageGrowth > NO_CROP)
        {
            Crop->svTranspiration = 0.0;
            Crop->svTranspirationPotential = 0.0;
            Crop->svWaterStressFactor = 0.0;

            for (i = 0; i < Soil->totalLayers; i++)
            {
                rootFraction[i] = 0.0;
                layerShootHC[i] = 0.0;
                layerPlantHC[i] = 0.0;
            }

            if (Crop->svTT_Cumulative > Crop->userEmergenceTT)
            {
                temperatureAvg = 0.5 * (Weather->tMax[y][doy - 1] + Weather->tMin[y][doy - 1]);
                factorTemperature = TemperatureLimitation (temperatureAvg, Crop->userTranspirationMinTemperature, Crop->userTranspirationThresholdTemperature);

                /* Calculate potential transpiration rate (kg/m2/d = mm/d) */
                PT = (1.0 + (Crop->userKc - 1.0) * Community->svRadiationInterception) * factorTemperature * Crop->svRadiationInterception * Weather->ETref[y][doy - 1];

                /* Calculate crop maximum water uptake rate (kg/m2/d = mm/d) */
                //PTx = PTx * factorTemperature * Crop->svRadiationInterception;
                PTx = Crop->transpirationMax * factorTemperature * Crop->svRadiationInterception;

                /* Calculate maximum crop transpiration rate (kg/m2/d = mm/d) */
                TE = PT < PTx ? PT : PTx;

                /* Calculate root fraction per soil layer */
                CalcRootFraction (rootFraction, Soil, Crop);

                /* Calculate plant hydraulic conductivity (kg^2)/(m2 J d)
                 * This is the hydraulic conductivity of a canopy fully covering the ground */
                plantHC = PTx / (SWP_FC - Crop->LWP_StressOnset);
                rootHC = plantHC / 0.65;
                shootHC = plantHC / 0.35;

                /* Adjust plant hydraulic conductance based on soil dryness */
                Root_HC_Adjustment = 0.0;
                for (i = 0; i < Soil->totalLayers; i++)
                {
                    rootActivity[i] = 1.0;
                    layerSalinityFactor[i] = 1.0;
                    rootActivity[i] = 1.0 - pow ((soilWP[i] - SWP_FC) / (Crop->LWP_WiltingPoint - SWP_FC), 8.0);
                    rootActivity[i] = rootActivity[i] > 1.0 ? 1.0 : rootActivity[i];
                    rootActivity[i] = rootActivity[i] < 0.0 ? 0.0 : rootActivity[i];

                    Layer_Root_HC_Adjustment[i] = rootActivity[i] * rootFraction[i] * layerSalinityFactor[i];
                    layerRootHC[i] = rootHC * Layer_Root_HC_Adjustment[i];
                    Root_HC_Adjustment += Layer_Root_HC_Adjustment[i];
                }

                for (i = 0; i < Soil->totalLayers; i++)
                {
                    if (Layer_Root_HC_Adjustment[i] > 0.0)
                    {
                        layerShootHC[i] = shootHC * Layer_Root_HC_Adjustment[i] / Root_HC_Adjustment;
                        layerPlantHC[i] = layerRootHC[i] * layerShootHC[i] / (layerRootHC[i] + layerShootHC[i]);
                    }
                    else
                        layerPlantHC[i] = 0.0;
                }

                rootHC *= Root_HC_Adjustment;
                plantHC = (rootHC * shootHC) / (rootHC + shootHC);

                if (plantHC > 0.0)
                {
                    /* Calculate average soil water potential (J/kg) */
                    SWP_Average = 0.0;

                    for (i = 0; i < Soil->totalLayers; i++)
                        SWP_Average += soilWP[i] * Layer_Root_HC_Adjustment[i] / Root_HC_Adjustment;

                    /* Calculate leaf water potential */
                    LWP = SWP_Average - TE / plantHC;

                    if (LWP < Crop->LWP_StressOnset)
                        LWP = (plantHC * SWP_Average * (Crop->LWP_StressOnset - Crop->LWP_WiltingPoint) + Crop->LWP_WiltingPoint * TE) / (plantHC * (Crop->LWP_StressOnset - Crop->LWP_WiltingPoint) + TE);

                    if (LWP < Crop->LWP_WiltingPoint)
                        LWP = Crop->LWP_WiltingPoint;

                    /* Reduce transpiration when LWP < LWP at the onset of stomatal
                     * closure */
                    if (LWP < Crop->LWP_StressOnset)
                    {
                        TA = TE * (LWP - Crop->LWP_WiltingPoint) / (Crop->LWP_StressOnset - Crop->LWP_WiltingPoint);
                        transpirationRatio = TA / TE;
                    }
                    else
                        transpirationRatio = 1.0;

                    /* Calculate crop water uptake (kg/m2/d = mm/d) */
                    for (i = 0; i < Soil->totalLayers; i++)
                    {
                        //Soil->waterUptake[i] += layerPlantHC[i] * (soilWP[i] - LWP) * transpirationRatio;
                        waterUptake[i] = layerPlantHC[i] * (soilWP[i] - LWP);
                        Soil->waterUptake[i] += waterUptake[i];
                    }
                }

                Crop->svTranspiration = 0.0;
                for (i = 0; i < Soil->totalLayers; i++)
                {
                    Crop->svTranspiration += waterUptake[i];
                }
                Crop->svTranspirationPotential = TE;
                Crop->svWaterStressFactor = 1.0 - Crop->svTranspiration / TE;
            }	/* end plant growing */
        }
    }

    for (i = 0; i < Soil->totalLayers; i++)
        Soil->waterContent[i] -= Soil->waterUptake[i] / (Soil->layerThickness[i] * WATER_DENSITY);
}

void CalcRootFraction (double *fractionRootsByLayer, SoilStruct *Soil, CropStruct *Crop)
{
    /* 
     * This function computes root fraction in each layer
     * it assumes that the roots will reach a final distribution
     * Root=a*exp[-b*z], with b~4
     * Root_layer = a / b * (Exp(-b * z1) - Exp(-b * z2)),
     * where z1 and z2 are the top and bottom of layer to compute the
     * progression, the roots are "assumed" to have reached the maximum depth
     * then the root lenght in each layer is computed only the roots within
     * the current rooting depth are considered for the calculation the root
     * fraction is then recalculated for the layers with roots the net effect
     * is that the roots are not sharply exponential during early growth
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * a		    double
     * b		    double	(1/m)
     * rootIntegral	    double
     * rootSum		    double
     * rootDistribution	    double[]
     * cumulativeRootingDepth
     *			    double
     * z1		    double
     * z2		    double
     * i		    int
     * j		    int
     */
    const double    a = 1.0;
    const double    b = 4.0;
    double          rootIntegral;
    double          rootSum;
    double          rootDistribution[Soil->totalLayers];
    double          cumulativeRootingDepth = 0.0;
    double          z1, z2;
    int             i, j;

    rootIntegral = a / b * (exp (-b * 0.0) - exp (-b * Crop->userMaximumRootingDepth));

    j = 0;
    while (cumulativeRootingDepth < Crop->svRootingDepth && j < Soil->totalLayers)
    {
        if (Soil->cumulativeDepth[j] < Crop->svRootingDepth)
            cumulativeRootingDepth = Soil->cumulativeDepth[j];
        else if (Soil->cumulativeDepth[j] >= Crop->svRootingDepth)
            cumulativeRootingDepth = Crop->svRootingDepth;

        if (j == 0)
            z1 = 0.0;
        else
            z1 = Soil->cumulativeDepth[j - 1];
        z2 = cumulativeRootingDepth;
        rootDistribution[j] = ((a / b) * (exp (-b * z1) - exp (-b * z2))) / rootIntegral;
        j++;
    }

    rootSum = 0.0;
    for (i = 0; i < j; i++)
        rootSum = rootSum + rootDistribution[i];

    /* Ensures sum fraction = 1 */
    for (i = 0; i < j; i++)
    {	/* Exits loop on the same layer as the previous loop */
        fractionRootsByLayer[i] = rootDistribution[i] / rootSum;
    }
}

double TemperatureLimitation (double T, double T_Min, double T_Threshold)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * t_limitation	    double	[return value]
     */
    double          t_limitation;

    if (T <= T_Min)
        t_limitation = 0.01;
    else if (T > T_Min && T <= T_Threshold)
        t_limitation = (T - T_Min) / (T_Threshold - T_Min);
    else
        t_limitation = 1.0;

    return (t_limitation);
}
