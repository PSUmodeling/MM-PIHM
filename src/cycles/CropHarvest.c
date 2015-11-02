#include "Cycles.h"

void GrainHarvest (int y, int doy, int startYear, CropStruct *Crop, ResidueStruct *Residue, const SoilStruct *Soil, SoilCarbonStruct *SoilCarbon, const WeatherStruct *Weather, const char *project)
{
    /*
     * Update roots and residue biomass at harvest
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * HI		    double	
     * NHI		    double
     * forageMassLoss	    double
     * forageYield	    double
     * retainedResidue	    double
     * grainNitrogenYield   double
     * forageNitrogenYield  double
     * forageNitrogenConcentration
     *                      double
     */
    double          HI, NHI;
    double          forageMassLoss;
    double          forageYield;
    double          retainedResidue;
    double          grainNitrogenYield;
    double          forageNitrogenYield;
    double          forageNitrogenConcentration;

    if (verbose_mode)
        printf ("DOY %3.3d %-20s %s\n", doy, "Grain Harvest", Crop->cropName);

    HI = ComputeHarvestIndex (Crop->userHIx, Crop->userHIo, Crop->userHIk, Crop->svShoot, Crop->svPostFloweringShootBiomass);
    NHI = pow (HI, 1.0 - 2.0 * (1.0 - HI) * HI);
    forageMassLoss = Crop->svShoot * (1.0 - HI);
    forageYield = forageMassLoss * Crop->userFractionResidueRemoved;
    grainNitrogenYield = Crop->svN_Shoot * NHI;
    forageNitrogenYield = (Crop->svN_Shoot - grainNitrogenYield) * Crop->userFractionResidueRemoved;

    if (forageYield > 0.0)
        forageNitrogenConcentration = 100.0 * forageNitrogenYield / forageYield;
    else
        forageNitrogenConcentration = 0.0;

    retainedResidue = forageMassLoss - forageYield;
    Residue->stanResidueMass += retainedResidue * Crop->userFractionResidueStanding;
    Residue->flatResidueMass += retainedResidue * (1.0 - Crop->userFractionResidueStanding);
    Residue->stanResidueN += (Crop->svN_Shoot - grainNitrogenYield - forageNitrogenYield) * Crop->userFractionResidueStanding;
    Residue->flatResidueN += (Crop->svN_Shoot - grainNitrogenYield - forageNitrogenYield) * (1.0 - Crop->userFractionResidueStanding);
    /* Assume 33% residue moisture at harvest */
    Residue->stanResidueWater += retainedResidue * Crop->userFractionResidueStanding / 10.0 * 0.5;
    Residue->flatResidueWater += retainedResidue * (1.0 - Crop->userFractionResidueStanding) / 10.0 * 0.5;

    /* Add roots of harvested annual crop to a root residue pool in each
     * layer */
    DistributeRootDetritus (y, Crop->svRoot, 0.0, Crop->svN_Root, 0.0, Soil, Crop, Residue, SoilCarbon);

    /* Yearly output variables */
    Residue->yearResidueBiomass += retainedResidue;
    Residue->yearRootBiomass += Crop->svRoot;
    Residue->yearRhizodepositionBiomass += Crop->svRizho;

    /* Season outputs */
    Crop->rcBiomass = Crop->svBiomass;
    Crop->rcRoot = Crop->svRoot;
    Crop->rcGrainYield = Crop->svShoot * HI;
    Crop->rcResidueBiomass = retainedResidue;
    Crop->rcForageYield = forageYield;
    Crop->rcHarvestIndex = HI;
    Crop->rcTotalNitrogen = Crop->svN_Shoot + Crop->svN_Root;
    Crop->rcRootNitrogen = Crop->svN_Root;
    Crop->rcGrainNitrogenYield = grainNitrogenYield;
    Crop->rcForageNitrogenYield = forageNitrogenYield;
    Crop->rcNitrogenCumulative = Crop->svN_StressCumulative;
    Crop->rcNitrogenInHarvest = (grainNitrogenYield + forageNitrogenYield) * 1000.0;
    Crop->rcNitrogenInResidue = (Crop->svN_Shoot + Crop->svN_Root - grainNitrogenYield - forageNitrogenYield) * 1000.0;
    Crop->rcNitrogenForageConc = forageNitrogenConcentration;

    PrintSeasonOutput (y, doy, startYear, Weather, Crop, project);

    KillCrop (Crop);
}

void ForageHarvest (int y, int doy, int startYear, CropStruct *Crop, ResidueStruct *Residue, const SoilStruct *Soil, SoilCarbonStruct *SoilCarbon, const WeatherStruct *Weather, const char *project)
{
    /*
     * 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * forageYield	    double	
     * forageMassLoss	    double
     * rootMassDead	    double
     * forageYieldNitrogen  double
     * forageNitrogenLoss	    double
     * rootNitrogenDead	    double
     * NStressCumulative    double	Stored for reporting before resetting
     *					  after harvest
     * fractionalHarvestLoss
     *			    double
     * fractionalNitrogenRemoval
     *			    double
     */
    double          forageYield;
    double          forageMassRetained = 0.0;
    double          forageMassRemoved = 0.0;
    double          forageYieldNitrogen;
    double          forageNitrogenRetained = 0.0;
    double          forageNitrogenRemoved = 0.0;
    double          forageMassLoss;
    double          forageNitrogenLoss;
    double          rootMassDead;
    double          rootNitrogenDead;
    double          manureCarbonMass = 0.0;
    double          manureNitrogenMass = 0.0;
    double          NH4Mass = 0.0;
    double          NStressCumulative;
    const double    fractionalHarvestLoss = 0.05;
    double          fractionalNitrogenRemoval;
    double          clippingFraction;

    if (verbose_mode)
        printf ("DOY %3.3d %-20s %s\n", doy, "Forage Harvest", Crop->cropName);

    /* With this method, forage yield is accumulated over the life of the crop
     * for season harvest (not per year) unless a new variable is created ...
     * Fractional_Nitrogen_Removal allows for N to remain in the crop to avoid
     * harvest-induced N stress. The negative effect is that the removed
     * biomass has a lower N concentration than the remaining biomass, which
     * is likely incorrect to fix this, N should be allocated to different
     * compartments and re-allocated after a clipping event
     *
     * When grazing, 50% of the nitrogen returns as urine, and 50% of the
     * nitrogen returns as manure. For carbon, 50% is respired by the animals
     * and 50% returns as manure */

    NStressCumulative = Crop->svN_StressCumulative;

    clippingFraction = (Crop->svShoot - Crop->userClippingBiomassThresholdLower * (1.0 - exp (-Crop->userPlantingDensity))) / Crop->svShoot / (1.0 - fractionalHarvestLoss);
    clippingFraction = (clippingFraction < Crop->userFractionResidueRemoved) ? clippingFraction : Crop->userFractionResidueRemoved;

    forageYield = Crop->svShoot * clippingFraction * (1.0 - fractionalHarvestLoss);
    fractionalNitrogenRemoval = 1.0 - pow (1.0 - clippingFraction, 0.7);
    forageYieldNitrogen = Crop->svN_Shoot * fractionalNitrogenRemoval * (1.0 - fractionalHarvestLoss);
    forageMassLoss = Crop->svShoot * clippingFraction * fractionalHarvestLoss;
    forageNitrogenLoss = Crop->svN_Shoot * clippingFraction * fractionalHarvestLoss;
    rootMassDead = Crop->svRoot * clippingFraction;
    rootNitrogenDead = Crop->svN_Root * clippingFraction;

    switch (Crop->userClippingDestiny)
    {
        case REMOVE_CLIPPING:
            forageMassRemoved = forageYield;
            forageMassRetained = 0.0;
            forageNitrogenRemoved = forageYieldNitrogen;
            forageNitrogenRetained = 0.0;
            break;
        case RETURN_CLIPPING:
            forageMassRemoved = 0.0;
            forageMassRetained = forageYield;
            forageNitrogenRemoved = 0.0;
            forageNitrogenRetained = forageYieldNitrogen;
            break;
        case GRAZING_CLIPPING:
            forageMassRemoved = 0.0;
            forageMassRetained = 0.0;
            forageNitrogenRemoved = 0.0;
            forageNitrogenRetained = 0.0;
            NH4Mass = 0.5 * forageYieldNitrogen;
            manureNitrogenMass = 0.5 * forageYieldNitrogen;
            manureCarbonMass = 0.5 * forageYield * FRACTION_CARBON_PLANT;
            /* Carbon balance alert:
             * Undocumented carbon respiration losses by animals */
            break;
        default:
            break;
    }

    /* Add roots of clipped crop to a root residue pool in each layer */
    DistributeRootDetritus (y, rootMassDead, 0.0, rootNitrogenDead, 0.0, Soil, Crop, Residue, SoilCarbon);

    Crop->svShoot -= forageYield + forageMassLoss;
    Crop->svRoot -= rootMassDead;
    Crop->svBiomass = Crop->svShoot + Crop->svRoot;
    Crop->svN_Shoot -= (forageYieldNitrogen + forageNitrogenLoss);
    Crop->svN_Root -= rootNitrogenDead;
    /* Resetting of cumulative nitrogen stress after haverst */
    Crop->svN_StressCumulative *= (Crop->svTT_Cumulative - Crop->userEmergenceTT) * (1.0 - pow (clippingFraction, 0.75)) / Crop->calculatedMaturityTT;
    Crop->svTT_Cumulative = Crop->userEmergenceTT + (Crop->svTT_Cumulative - Crop->userEmergenceTT) * (1.0 - pow (clippingFraction, 0.75));
    Crop->svRadiationInterception = Crop->svRadiationInterception * (1.0 - pow (clippingFraction, 0.75));
    Crop->svRadiationInterception_nc = Crop->svRadiationInterception_nc * (1.0 - pow (clippingFraction, 0.75));
    Crop->svShootUnstressed = Crop->svShoot;

    Residue->stanResidueMass += (forageMassLoss + forageMassRetained) * Crop->userFractionResidueStanding;
    Residue->flatResidueMass += (forageMassLoss + forageMassRetained) * (1.0 - Crop->userFractionResidueStanding);
    Residue->stanResidueN += (forageNitrogenLoss + forageNitrogenRetained) * Crop->userFractionResidueStanding;
    Residue->flatResidueN += (forageNitrogenLoss + forageNitrogenRetained) * (1.0 - Crop->userFractionResidueStanding);
    /* Assume 33% residue moisture at harvest */
    Residue->stanResidueWater += (forageMassLoss + forageMassRetained) * Crop->userFractionResidueStanding / 10.0 * 0.5;
    Residue->flatResidueWater += (forageMassLoss + forageMassRetained) * (1.0 - Crop->userFractionResidueStanding) / 10.0 * 0.5;

    Residue->manureSurfaceC += manureCarbonMass;
    Residue->manureSurfaceN += manureNitrogenMass;
    Soil->NH4[0] += NH4Mass;

    /* Yearly output variables */
    Residue->yearResidueBiomass += forageMassLoss + forageMassRetained;
    Residue->yearRootBiomass += rootMassDead;

    /* Season outputs */
    Crop->rcBiomass += forageYield + forageMassLoss + rootMassDead;
    Crop->rcRoot += rootMassDead;
    Crop->rcForageYield += forageYield;
    Crop->rcResidueBiomass += forageMassLoss + forageMassRetained;
    Crop->rcTotalNitrogen = forageYieldNitrogen + forageNitrogenLoss + rootNitrogenDead;
    Crop->rcRootNitrogen = rootNitrogenDead;
    Crop->rcForageNitrogenYield = forageYieldNitrogen;
    //Crop->rcNitrogenCumulative = NStressCumulative / Crop->userClippingTiming;

    PrintSeasonOutput (y, doy, startYear, Weather, Crop, project);
}

void HarvestCrop (int y, int doy, int startYear, CropStruct *Crop, ResidueStruct *Residue, const SoilStruct *Soil, SoilCarbonStruct *SoilCarbon, const WeatherStruct *Weather, const char *project)
{
    /* 
     * Set crop status to Killed
     * Final crop values based on a killed crop performed
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * forageMassLoss	    double
     * rootMassDead	    double
     * forageNitrogenLoss	    double
     * rootNitrogenDead	    double
     */

    double          forageMassLoss;
    double          rootMassDead;
    double          forageNitrogenLoss;
    double          rootNitrogenDead;

    forageMassLoss = Crop->svShoot;
    rootMassDead = Crop->svRoot;
    forageNitrogenLoss = Crop->svN_Shoot;
    rootNitrogenDead = Crop->svN_Root;

    Residue->stanResidueMass += forageMassLoss * Crop->userFractionResidueStanding;
    Residue->flatResidueMass += forageMassLoss * (1.0 - Crop->userFractionResidueStanding);
    Residue->stanResidueN += forageNitrogenLoss * Crop->userFractionResidueStanding;
    Residue->flatResidueN += forageNitrogenLoss * (1.0 - Crop->userFractionResidueStanding);
    /* Assume 33% residue moisture at harvest */
    Residue->stanResidueWater += forageMassLoss * Crop->userFractionResidueStanding / 10.0 * 0.5;
    Residue->flatResidueWater += forageMassLoss * (1.0 - Crop->userFractionResidueStanding) / 10.0 * 0.5;

    /* Add roots of terminated crop to a root residue pool in each layer */
    DistributeRootDetritus (y, rootMassDead, 0.0, rootNitrogenDead, 0.0, Soil, Crop, Residue, SoilCarbon);

    /* Yearly output variables */
    Residue->yearResidueBiomass += forageMassLoss;
    Residue->yearRootBiomass += rootMassDead;
    Residue->yearRhizodepositionBiomass += Crop->svRizho;

    /* Season outputs */
    Crop->rcBiomass += Crop->svBiomass;
    Crop->rcRoot += Crop->svRoot;
    Crop->rcResidueBiomass += forageMassLoss;

    Crop->rcTotalNitrogen = Crop->svN_Shoot + Crop->svN_Root;
    Crop->rcRootNitrogen = Crop->svN_Root;
    Crop->rcNitrogenCumulative = Crop->svN_StressCumulative;
    Crop->rcNitrogenInResidue = (Crop->svN_Shoot + Crop->svN_Root) * 1000.0;

    if (verbose_mode)
        printf ("DOY %3.3d %-20s %s\n", doy, "Harvest", Crop->cropName);

    PrintSeasonOutput (y, doy, startYear, Weather, Crop, project);

    KillCrop (Crop);
}

void DistributeRootDetritus (int y, double rootMass, double rhizoMass, double rootN, double rhizoN, const SoilStruct *Soil, const CropStruct *Crop, ResidueStruct *Residue, SoilCarbonStruct *SoilCarbon)
{
    /*
     * This subroutine distributes dead roots and rizhodeposition in soil
     * layer. The function always computes the fraction of roots based on
     * maximum root depth for the crop then, it corrects the allocation using
     * the actual rooting depth this prevents an over-representation of the
     * roots in the top layer
     * The function is the same than that for root distribution used in water
     * uptake, but is kept in a separate sub in case parameterization
     * changes
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * a		    double
     * b		    double	(1/m)
     * rootIntegral	    double
     * rootSum		    double
     * cumulativeRootingDepth
     *			    double
     * z1		    double
     * z2		    double
     * i		    int
     * j		    int
     * rootDistribution	    double[Soil->totalLayers]
     * fractionRootsByLayer double[Soil->totalLayers]
     */

    /* Function distribution Root=a*exp[-b*z], with b~4 */
    /* Integral to compute mass in each layer
     * Root_layer = a / b * (Exp(-b * z1) - Exp(-b * z2)),
     * where z1 and z2 are the top and bottom of layer */

    const double    a = 1.0;
    const double    b = 4.0;
    double          rootIntegral;
    double          rootSum;
    double          cumulativeRootingDepth = 0.0;
    double          z1, z2;
    int             i, j = 0;
    double          rootDistribution[Soil->totalLayers];
    double          fractionRootsByLayer[Soil->totalLayers];

    rootIntegral = a / b * (exp (-b * 0.0) - exp (-b * Crop->userMaximumRootingDepth));

    for (i = 0; i < Soil->totalLayers; i++)
    {
        rootDistribution[i] = 0.0;
        fractionRootsByLayer[i] = 0.0;
    }

    j = 0;
    while (cumulativeRootingDepth < Crop->svRootingDepth && j < Soil->totalLayers)
    {
        if (Soil->cumulativeDepth[j] < Crop->svRootingDepth)
            cumulativeRootingDepth = Soil->cumulativeDepth[j];
        else
            cumulativeRootingDepth = Crop->svRootingDepth;

        if (j == 0)
            z1 = 0.0;
        else
            z1 = Soil->cumulativeDepth[j - 1];

        z2 = cumulativeRootingDepth;

        rootDistribution[j] = ((a / b) * (exp (-b * z1) - exp (-b * z2))) / rootIntegral;
        j++;
    }

    /* Ensures that cumulative fractional root distribution = 1 */
    rootSum = 0.0;
    for (i = 0; i < j; i++)
        rootSum = rootSum + rootDistribution[i];

    /* Compute input of biomass from roots to each layer */
    for (i = 0; i < j; i++)
    {   /* exits loop on the same layer as the previous loop */
        if (rootMass > 0.0)
        {
            fractionRootsByLayer[i] = rootDistribution[i] / rootSum;
            Residue->residueRt[i] += fractionRootsByLayer[i] * rootMass;
            Residue->residueRtN[i] += fractionRootsByLayer[i] * rootN;
            SoilCarbon->rootBiomassInput[i] += fractionRootsByLayer[i] * rootMass;
        }

        if (rhizoMass > 0.0)
        {
            Residue->residueRz[i] += fractionRootsByLayer[i] * rhizoMass;
            Residue->residueRzN[i] += fractionRootsByLayer[i] * rhizoN;
            SoilCarbon->rhizBiomassInput[i] += fractionRootsByLayer[i] * rhizoMass;
        }
    }
}

double ComputeHarvestIndex (double HIx, double HIo, double HIk, double cumulativeShoot, double cumulativePostFloweringShootBiomass)
{
    /*
     * Update roots and residue biomass at harvest
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * fg		    double	Fractional post-anthesis growth
     * harvest_index	    double	Harvest index [return value]
     */
    double          fg;
    double          harvest_index;

    if (GT (cumulativePostFloweringShootBiomass, 0.0) && GT (cumulativeShoot, 0.0))
    {
        fg = cumulativePostFloweringShootBiomass / cumulativeShoot;
        harvest_index = HIx - (HIx - HIo) * exp (-HIk * fg);
    }
    else
        harvest_index = 0.0;

    return (harvest_index);
}
