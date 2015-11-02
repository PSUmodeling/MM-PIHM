#include "Cycles.h"

void InitializeResidue (ResidueStruct *Residue, int totalYears, int totalLayers)
{
    Residue->residueAbgd = (double *)malloc (totalLayers * sizeof (double));
    Residue->residueRt = (double *)malloc (totalLayers * sizeof (double));
    Residue->residueRz = (double *)malloc (totalLayers * sizeof (double));
    Residue->residueAbgdN = (double *)malloc (totalLayers * sizeof (double));
    Residue->residueRtN = (double *)malloc (totalLayers * sizeof (double));
    Residue->residueRzN = (double *)malloc (totalLayers * sizeof (double));
    Residue->manureC = (double *)malloc (totalLayers * sizeof (double));
    Residue->manureN = (double *)malloc (totalLayers * sizeof (double));

    Residue->residueInterception = 0.0;
    Residue->stanResidueTau = 1.0;
    Residue->flatResidueTau = 1.0;
    Residue->stanResidueMass = 0.0;
    Residue->flatResidueMass = 0.0;
    Residue->stanResidueN = 0.0;
    Residue->flatResidueN = 0.0;
    Residue->manureSurfaceC = 0.0;
    Residue->manureSurfaceN = 0.0;
    Residue->stanResidueWater = 0.0;
    Residue->flatResidueWater = 0.0;

    Residue->yearResidueBiomass = 0.0;
    Residue->yearRootBiomass = 0.0;
    Residue->yearRhizodepositionBiomass = 0.0;
}

void ComputeResidueCover (ResidueStruct *Residue)
{
    /*
     * Compute residue cover
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * stanResidueAI	    double	Standing residue area index
     *					  (m2 residue / m2 ground)
     * flatResidueAI	    double	Flat residue area index
     *					  (m2 residue / m2 ground)
     */
    double          stanResidueAI;
    double          flatResidueAI;

    /* Factor 0.1 converts from Mg/ha to kg/m2 */
    stanResidueAI = STAN_RESIDUE_SA * Residue->stanResidueMass * 0.1;
    flatResidueAI = FLAT_RESIDUE_SA * Residue->flatResidueMass * 0.1;

    Residue->stanResidueTau = exp (-STAN_RESIDUE_K * stanResidueAI);
    Residue->flatResidueTau = exp (-FLAT_RESIDUE_K * flatResidueAI);
    Residue->residueInterception = (1.0 - Residue->stanResidueTau) + Residue->stanResidueTau * (1.0 - Residue->flatResidueTau);
}

void ResidueWetting (ResidueStruct *Residue, SoilStruct *Soil)
{
    /*
     * Compute residue wetting
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * residueMaxWaterConcentration
     *			    double	(kg/kg)
     * flatResidueWaterDeficit
     *			    double	Water need to saturate residue (mm)
     * standingResidueWaterDeficit
     *			    double	Water need to saturate residue (mm)
     * waterWettingResidue  double	Amount of water interceptable by
     *					  residue (mm)
     * waterRetainedResidue double	Water retained in residue and
     *					  discounted for infiltration (mm)
     */
    const double    residueMaxWaterConcentration = 3.3;
    double          flatResidueWaterDeficit;
    double          standingResidueWaterDeficit;
    double          waterWettingResidue;
    double          waterRetainedResidue;

    /* 10 converts residue from Mg/ha to kg/m2 */
    flatResidueWaterDeficit = residueMaxWaterConcentration * Residue->flatResidueMass / 10.0 - Residue->flatResidueWater;
    standingResidueWaterDeficit = residueMaxWaterConcentration * Residue->stanResidueMass / 10.0 - Residue->stanResidueWater;

    waterWettingResidue = Soil->infiltrationVol * Residue->residueInterception;

    waterRetainedResidue = 0.0;

    /* Wet flat residue first */
    if (waterWettingResidue > flatResidueWaterDeficit)
    {
        Residue->flatResidueWater += flatResidueWaterDeficit;
        waterRetainedResidue += flatResidueWaterDeficit;
        waterWettingResidue -= flatResidueWaterDeficit;
    }
    else
    {
        Residue->flatResidueWater += waterWettingResidue;
        waterRetainedResidue += waterWettingResidue;
        waterWettingResidue -= waterWettingResidue;
    }

    if (waterWettingResidue > standingResidueWaterDeficit)
    {
        Residue->stanResidueWater += standingResidueWaterDeficit;
        waterRetainedResidue += standingResidueWaterDeficit;
        waterWettingResidue -= standingResidueWaterDeficit;
    }
    else
    {
        Residue->stanResidueWater += waterWettingResidue;
        waterRetainedResidue += waterWettingResidue;
        waterWettingResidue -= waterWettingResidue;
    }

    Soil->infiltrationVol -= waterRetainedResidue;
}

void ResidueEvaporation (ResidueStruct *Residue, SoilStruct *Soil, const CommunityStruct *Community, double ETo, double snowCover)
{
    /*
     * Compute residue evaporation
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * residueMaxWaterConcentration
     *			    double	(kg water/kg residue)
     * flatEvapFactor	    double
     * standingEvapFactor   double
     * flatEvap		    double
     * standingEvap	    double
     * residueEvapDemand    double
     * xx		    double
     */
    const double    residueMaxWaterConcentration = 3.3;
    double          flatEvapFactor;
    double          standingEvapFactor;
    double          flatEvap;
    double          standingEvap;
    double          residueEvapDemand;
    double          xx;

    Soil->residueEvaporationVol = 0.0;
    if (Residue->stanResidueWater > 0.0 || Residue->flatResidueWater > 0.0)
    {
        Soil->residueEvaporationVol = 0.0;
        residueEvapDemand = Residue->residueInterception * (1.0 - snowCover) * (1.0 - Community->svRadiationInterception) * ETo;
	/* 10 converts residue from Mg/ha to kg/m2 */
        standingEvapFactor = pow (Residue->stanResidueWater / (residueMaxWaterConcentration * Residue->stanResidueMass / 10.0), 2.0);
        flatEvapFactor = pow (Residue->flatResidueWater / (residueMaxWaterConcentration * Residue->flatResidueMass / 10.0), 2.0);

        /* Dry standing residue first */
        xx = residueEvapDemand * standingEvapFactor;
        if (Residue->stanResidueWater >= xx)
            standingEvap = xx;
        else
            standingEvap = Residue->stanResidueWater;
        residueEvapDemand -= standingEvap;

        xx = residueEvapDemand * flatEvapFactor;
        if (Residue->flatResidueWater >= xx)
            flatEvap = xx;
        else
            flatEvap = Residue->flatResidueWater;

        Soil->residueEvaporationVol = standingEvap + flatEvap;
        Residue->stanResidueWater -= standingEvap;
        Residue->flatResidueWater -= flatEvap;
    }
}
