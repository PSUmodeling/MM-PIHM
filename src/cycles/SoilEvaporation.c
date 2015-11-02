#include "Cycles.h"

void Evaporation (SoilStruct *Soil, const CommunityStruct *Community, ResidueStruct *Residue, double ETo, double SnowCover)
{
    /*
     * 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * i		    int		Loop counter
     * EvaporativeDemand    double	(mm/day)
     * WC_AirDry	    double	(m3/m3)
     * layerTop		    double[]	(m)
     * layerBottom	    double[]	(m)
     * layerMidpoint	    double	(m)
     * WaterAvailable	    double	Layer available water for evaporation
     *					  (mm)
     * WaterSupply	    double	Layer water supply (mm)
     * DepthLimitation	    double	Factor that limits evaporable water
     *					  based on depth
     * WaterContentLimitation
     *			    double	Factor that limits evaporable water
     *					  based on WC
     * Evaporation	    double	Evaporation by layer (mm)
     * WCi		    double[]	Copy initial soil water content
     * EvapFlux		    double[]	Store flux from each layer to compute
     *					  solute transport
     */
    int             i;
    double          EvaporativeDemand;
    double          WC_AirDry;
    double          layerTop[Soil->totalLayers];
    double          layerBottom[Soil->totalLayers];
    double          layerMidpoint;
    double          WaterAvailable;
    double          WaterSupply;
    double          DepthLimitation;
    double          WaterContentLimitation;
    double          Evaporation;
    double          WCi[Soil->totalLayers];
    double          EvapFlux[Soil->totalLayers];

    for (i = 0; i < Soil->totalLayers; i++)
        WCi[i] = Soil->waterContent[i];

    Soil->evaporationVol = 0.0;

    /* It uses the maximum cover of either flat residues or snow to calculate
     * radiation reaching the soil surface */
    if (1.0 - Residue->flatResidueTau >= SnowCover)
        EvaporativeDemand = (1.0 - Residue->residueInterception) * (1.0 - Community->svRadiationInterception) * ETo;
    else
        EvaporativeDemand = Residue->stanResidueTau * (1.0 - SnowCover) * (1.0 - Community->svRadiationInterception) * ETo;

    for (i = 0; i < Soil->totalLayers; i++)
    {
        if (i > 0)
        {
            layerBottom[i] = layerBottom[i - 1] + Soil->layerThickness[i];
            layerTop[i] = layerBottom[i - 1];
        }
        else
        {
            layerBottom[i] = Soil->layerThickness[i];
            layerTop[i] = 0.0;
        }

        layerMidpoint = 0.5 * (layerTop[i] + layerBottom[i]);
        WC_AirDry = Soil->PWP[i] / 3.0;	    /* An approximation to air dry */
        WaterAvailable = (Soil->waterContent[i] - WC_AirDry) * Soil->layerThickness[i] * WATER_DENSITY;

        DepthLimitation = 1.0 / 3.0 * (Depth_Limitation_To_Evaporation (layerTop[i]) + Depth_Limitation_To_Evaporation (layerMidpoint) + Depth_Limitation_To_Evaporation (layerBottom[i]));
        WaterContentLimitation = Water_Content_Limitation_To_Evaporation (Soil->FC[i], WC_AirDry, Soil->waterContent[i]);
        WaterSupply = WaterAvailable * DepthLimitation * WaterContentLimitation;

        if (WaterSupply > EvaporativeDemand)
            Evaporation = EvaporativeDemand;
        else
            Evaporation = WaterSupply;

        EvaporativeDemand -= Evaporation;
        Soil->evaporationVol += Evaporation;
        Soil->waterContent[i] -= Evaporation / (Soil->layerThickness[i] * WATER_DENSITY);
        EvapFlux[i] = Evaporation;

        if (EvaporativeDemand == 0.0)
            break;
    }

    /* Try changing retention constant to slow down nitrate upward movement,
     * try 1 (defaul is 0 for nitrate downward in liquid flow) */
    /* 2015/02/23 The solute transport sub needs to be re-coded to take into
     * account bi-directional solute flow, for now it is removed */
    //if (Soil->evaporationVol > 0.0)
    //{
    //    SoluteTransportEvaporation (Soil->totalLayers, 1.0, EvapFlux, Soil->NO3, Soil->BD, Soil->layerThickness, Soil->Porosity, WCi);
    //    SoluteTransportEvaporation (Soil->totalLayers, 5.6, EvapFlux, Soil->NH4, Soil->BD, Soil->layerThickness, Soil->Porosity, WCi);
    //}
}

double Depth_Limitation_To_Evaporation (double Depth)
{
    return (1.0 / (1.0 + pow (Depth / 0.04, 3.5)));
}

double Water_Content_Limitation_To_Evaporation (double FC, double WC_AirDry, double WC)
{
    return (pow ((WC - WC_AirDry) / (FC - WC_AirDry), 3.0));
}
