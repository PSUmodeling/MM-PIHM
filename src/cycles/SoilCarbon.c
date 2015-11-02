#include "Cycles.h"

void InitializeSoilCarbon (SoilCarbonStruct *SoilCarbon, int totalLayers)
{
    SoilCarbon->factorComposite = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->carbonRespired = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->rootBiomassInput = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->rhizBiomassInput = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->abgdBiomassInput = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->rootCarbonInput = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->rhizCarbonInput = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->manuCarbonInput = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->abgdCarbonInput = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->carbonMassInitial = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->carbonMassFinal = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->annualDecompositionFactor = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->annualSoilCarbonDecompositionRate = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->annualCarbonInputByLayer = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->annualHumifiedCarbonMass = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->annualRespiredCarbonMass = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->annualRespiredResidueCarbonMass = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->annualHumificationCoefficient = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->annualNmineralization = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->annualNImmobilization = (double *)calloc (totalLayers, sizeof (double));
    SoilCarbon->annualNNetMineralization = (double *)calloc (totalLayers, sizeof (double));
}

void ComputeFactorComposite (SoilCarbonStruct *SoilCarbon, int doy, int y, int last_doy, SoilStruct *Soil)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * waterPotential	    double	(J/kg)
     * airContent	    double	(m3/m3)
     * factorMoisture	    double	(-)
     * factorTemperature    double	(-)
     * factorAeration	    double	Cumulative aeration factor accounts
     *					  empirically for air content of
     *					  layers above that considered in the
     *					  calculations (-)
     * i		    int		Loop counter
     */
    double          waterPotential;
    double          airContent;
    double          factorMoisture;
    double          factorTemperature;
    double          factorAeration = 1.0;
    int             i;

    for (i = 0; i < Soil->totalLayers; i++)
    {
        waterPotential = SoilWaterPotential (Soil->Porosity[i], Soil->airEntryPotential[i], Soil->B_Value[i], Soil->waterContent[i]);
        airContent = Soil->Porosity[i] - Soil->waterContent[i];
        factorMoisture = Moisture (waterPotential);
        factorTemperature = TemperatureFunction (Soil->soilTemperature[i]);
        factorAeration *= Aeration (airContent);
        SoilCarbon->factorComposite[i] = factorMoisture * factorAeration * factorTemperature;

        if (doy == 1)
        {
            SoilCarbon->annualDecompositionFactor[i] = 0.0;
        }

        SoilCarbon->annualDecompositionFactor[i] += SoilCarbon->factorComposite[i];

        if (doy == last_doy)
            SoilCarbon->annualDecompositionFactor[i] /= ((double)last_doy);
    }
}

void ComputeSoilCarbonBalanceMB (SoilCarbonStruct *SoilCarbon, int y, ResidueStruct *Residue, SoilStruct *Soil, double *tillageFactor)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * i		    int		Loop counter
     * socDecompositionRate double
     * micrDecompositionRate
     *			    double
     * humifiedCarbon	    double
     * humifiedNitrogen	    double
     * abgdHumificationFactor
     *			    double
     * rootHumificationFactor
     *			    double
     * rhizHumificationFactor
     *			    double
     * manuHumificationFactor
     *			    double
     * micrHumificationFactor
     *			    double
     * soilMass		    double
     * satSOCConc	    double	(g C/kg soil)
     * humificationAdjustmentBySOC
     *			    double
     * decompositionAdjustmentBySOC
     *			    double
     * contactFractionFlat  double	Fraction of surface residues subject
     *					  to decomposition
     * contactFractionStan  double	Fraction of surface residues subject
     *					  to decomposition
     * xx0		    double	Residue mass decomposition (Mg/ha/day)
     * xx1		    double	Residue mass decomposition (Mg/ha/day)
     * xx2		    double	Residue mass decomposition (Mg/ha/day)
     * xx3		    double	Residue mass decomposition (Mg/ha/day)
     * xx4		    double	Residue mass decomposition (Mg/ha/day)
     * xx5		    double	Residue mass decomposition (Mg/ha/day)
     * xx6		    double	Manure carbon decomposition
     *					  (Mg/ha/day)
     * xx7		    double	Manure carbon decomposition
     *					  (Mg/ha/day)
     * xx8		    double	Organic matter decomposition
     *					  (Mg/ha/day)
     * xx9		    double	Microbial pool decomposition
     *					  (Mg/ha/day)
     * nm0		    double	Residue nitrogen net mineralization
     * nm1		    double	Residue nitrogen net mineralization
     * nm2		    double	Residue nitrogen net mineralization
     * nm3		    double	Residue nitrogen net mineralization
     * nm4		    double	Residue nitrogen net mineralization
     * nm5		    double	Residue nitrogen net mineralization
     * nm6		    double	Residue nitrogen net mineralization
     * nm7		    double	Residue nitrogen net mineralization
     * nm8		    double	Residue nitrogen net mineralization
     * nm9		    double	Residue nitrogen net mineralization
     * nr0		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr1		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr2		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr3		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr4		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr5		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr6		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr7		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr8		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr9		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nh0		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh1		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh2		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh3		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh4		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh5		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh6		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh7		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh8		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh9		    double	Residue Nitrogen trasfered from
     *					  microbial biomass to SOM
     * aux1		    double
     * aux2		    double
     * aux3		    double
     * xxPartialSum	    double
     * NMineralization	    double
     * NImmobilization	    double
     * NNetMineralization   double
     * stanCNRatio	    double	CN ratio of standing residues
     * flatCNRatio	    double	CN ratio of flat residues
     * abgdCNRatio	    double	CN ratio of aboveground residues (in
     *					  each soil layer)
     * rootCNRatio	    double	CN ratio of root residues
     * rhizCNRatio	    double	CN ratio of rizhodeposition
     * smanCNRatio	    double	CN ratio of surface manure
     * manuCNRatio	    double	CN ratio of manure
     * micrCNRatio	    double	CN ratio of microbial biomass
     * somCNratio	    double	sman = surface manure, manu = manure
     * CNnew		    double	CN of destiny pool (microbial biomass)
     *					  calculated before each microbial
     *					  attack; for microbial biomass it
     *					  equals CN microbial biomass
     * NMineralConcentration
     *			    double	(g N-NO3/g soil)
     * NMineral		    double	Sum of NO3 and NH4 (be careful with
     *					  units
     * NH4_Fraction	    double	Fraction of NH4 in the sum NO3 + NH4
     * decompReductionFactor
     *			    double
     * NInitial		    double
     * NFinal		    double
     */
    int             i;
    double          socDecompositionRate;
    double          micrDecompositionRate;
    double          humifiedCarbon;
    double          humifiedNitrogen;
    double          abgdHumificationFactor;
    double          rootHumificationFactor;
    double          rhizHumificationFactor;
    double          manuHumificationFactor;
    double          micrHumificationFactor;
    double          soilMass;
    double          satSOCConc;
    double          humificationAdjustmentBySOC;
    double          decompositionAdjustmentBySOC;
    double          contactFractionFlat;
    double          contactFractionStan;

    double          xx0, xx1, xx2, xx3, xx4, xx5, xx6, xx7, xx8, xx9;
    double          nm0, nm1, nm2, nm3, nm4, nm5, nm6, nm7, nm8, nm9;
    double          nr0, nr1, nr2, nr3, nr4, nr5, nr6, nr7, nr8, nr9;
    double          nh0, nh1, nh2, nh3, nh4, nh5, nh6, nh7, nh8, nh9;

    double          aux1, aux2, aux3;
    double          xxPartialSum;
    double          NMineralization;
    double          NImmobilization;
    double          NNetMineralization;
    double          stanCNRatio;
    double          flatCNRatio;
    double          abgdCNRatio;
    double          rootCNRatio; 
    double          rhizCNRatio;
    double          smanCNRatio;
    double          manuCNRatio;
    double          micrCNRatio;
    double          somCNratio;
    double          CNnew;
    double          NMineralConcentration;
    double          NMineral;
    double          NH4_Fraction;
    double          decompReductionFactor;
    double          NInitial, NFinal;

    NInitial = 0.0;
    NFinal = 0.0;
    NInitial = Residue->stanResidueN + Residue->flatResidueN + Residue->manureSurfaceN;
    contactFractionFlat = 0.0;
    contactFractionStan = 0.0;
    Soil->N_Immobilization = 0.0;
    Soil->N_Mineralization = 0.0;
    Soil->N_NetMineralization = 0.0;
    Soil->SOCProfile = 0.0;
    Soil->SONProfile = 0.0;
    Soil->C_Humified = 0.0;
    Soil->C_ResidueRespired = 0.0;
    Soil->C_SoilRespired = 0.0;

    for (i = 0; i < Soil->totalLayers; i++)
    {
        NInitial += Soil->SON_Mass[i] + Soil->MBN_Mass[i] + Soil->NO3[i] + Soil->NH4[i] + Residue->residueAbgdN[i] + Residue->residueRtN[i] + Residue->residueRzN[i] + Residue->manureN[i];

        NMineralConcentration = 0.0;
        NMineral = 0.0;
        NH4_Fraction = 0.0;
        decompReductionFactor = 1.0;

        socDecompositionRate = 0.0;
        micrDecompositionRate = 0.0;
        humifiedCarbon = 0.0;
        humifiedNitrogen = 0.0;
        NMineralization = 0.0;
        NImmobilization = 0.0;
        NNetMineralization = 0.0;
        CNnew = 0.0;
        stanCNRatio = 0.0;
        flatCNRatio = 0.0;
        abgdCNRatio = 0.0;
        rootCNRatio = 0.0;
        rhizCNRatio = 0.0;
        smanCNRatio = 0.0;
        manuCNRatio = 0.0;
        micrCNRatio = 0.0;

        xx0 = xx1 = xx2 = xx3 = xx4 = xx5 = xx6 = xx7 = xx8 = xx9 = 0.0;
        nm0 = nm1 = nm2 = nm3 = nm4 = nm5 = nm6 = nm7 = nm8 = nm9 = 0.0;
        nr0 = nr1 = nr2 = nr3 = nr4 = nr5 = nr6 = nr7 = nr8 = nr9 = 0.0;
        nh0 = nh1 = nh2 = nh3 = nh4 = nh5 = nh6 = nh7 = nh8 = nh9 = 0.0;

        /* Compute auxiliar variables */
	/* 10000 converts from Mg/m2 to Mg/ha */
        soilMass = 10000.0 * Soil->BD[i] * Soil->layerThickness[i];
        NMineral = Soil->NO3[i] + Soil->NH4[i];
        NH4_Fraction = Soil->NH4[i] / NMineral;
        NMineralConcentration = NMineral / soilMass;
        satSOCConc = 21.1 + 0.375 * Soil->Clay[i] * 100.0;

        /* Compute C/N ratios */
        if (i == 0)
        {
            if (Residue->stanResidueMass > 0.0)
                stanCNRatio = Residue->stanResidueMass * FRACTION_CARBON_PLANT / Residue->stanResidueN;
            if (Residue->flatResidueMass > 0.0)
                flatCNRatio = Residue->flatResidueMass * FRACTION_CARBON_PLANT / Residue->flatResidueN;
            if (Residue->manureSurfaceC > 0.0)
                smanCNRatio = Residue->manureSurfaceC / Residue->manureSurfaceN;
        }

        somCNratio = Soil->SOC_Mass[i] / Soil->SON_Mass[i];
        micrCNRatio = Soil->MBC_Mass[i] / Soil->MBN_Mass[i];

        if (Residue->residueAbgd[i] > 0.0)
            abgdCNRatio = Residue->residueAbgd[i] * FRACTION_CARBON_PLANT / Residue->residueAbgdN[i];
        if (Residue->residueRt[i] > 0.0)
            rootCNRatio = Residue->residueRt[i] * FRACTION_CARBON_PLANT / Residue->residueRtN[i];
        if (Residue->residueRz[i] > 0.0)
            rhizCNRatio = Residue->residueRz[i] * FRACTION_CARBON_RIZHO / Residue->residueRzN[i];
        if (Residue->manureC[i] > 0.0)
            manuCNRatio = Residue->manureC[i] / Residue->manureN[i];

        /* Humification */
        /* Humification reduction when C conc approaches saturation */
        humificationAdjustmentBySOC = 1.0 - pow (Soil->SOC_Conc[i] / satSOCConc, SOC_HUMIFICATION_POWER);
        humificationAdjustmentBySOC = humificationAdjustmentBySOC > 0.0 ? humificationAdjustmentBySOC : 0.0;

        abgdHumificationFactor = sqrt (MaximumAbgdHumificationFactor (Soil->Clay[i]) * humificationAdjustmentBySOC);
        rootHumificationFactor = sqrt (MaximumRootHumificationFactor (Soil->Clay[i]) * humificationAdjustmentBySOC);
        rhizHumificationFactor = sqrt (MaximumRhizHumificationFactor (Soil->Clay[i]) * humificationAdjustmentBySOC);
        manuHumificationFactor = sqrt (MaximumManuHumificationFactor (Soil->Clay[i]) * humificationAdjustmentBySOC);

        /* Temporarily assigned abgd humification, then checked if it can be
         * higher if manure is decomposing, but never lower */
        micrHumificationFactor = abgdHumificationFactor;

        /* Residue and manure decomposition */
        if (i == 0)
        {
            contactFractionStan = pow (Residue->stanResidueTau, exp (-1.5 / sqrt (1.0 - Residue->stanResidueTau)));
            contactFractionFlat = pow (Residue->flatResidueTau, exp (-1.5 / sqrt (1.0 - Residue->flatResidueTau)));

            xx2 = SoilCarbon->factorComposite[i] * MAXIMUM_RESIDUE_DECOMPOSITION_RATE * contactFractionStan * Residue->stanResidueMass;
            xx3 = SoilCarbon->factorComposite[i] * MAXIMUM_RESIDUE_DECOMPOSITION_RATE * contactFractionFlat * Residue->flatResidueMass;
            xx6 = SoilCarbon->factorComposite[i] * MAXIMUM_MANURE_DECOMPOSITION_RATE * Residue->manureSurfaceC;
        }

        xx1 = SoilCarbon->factorComposite[i] * MAXIMUM_RESIDUE_DECOMPOSITION_RATE * Residue->residueAbgd[i];
        xx4 = SoilCarbon->factorComposite[i] * MAXIMUM_ROOT_DECOMPOSITION_RATE * Residue->residueRt[i];
        xx5 = SoilCarbon->factorComposite[i] * MAXIMUM_RHIZO_DECOMPOSITION_RATE * Residue->residueRz[i];
        xx7 = SoilCarbon->factorComposite[i] * MAXIMUM_MANURE_DECOMPOSITION_RATE * Residue->manureC[i];

        /* Inorganic N limitation for decomposition */
        /* If decomposition > 0 then compute net N mineralization and
         * accumulate negatives */
        if (i == 0)
        {
            if (xx2 > 0.0)
            {
                CNnew = CNdestiny (NMineralConcentration, stanCNRatio);
                nm2 = NitrogenMineralization (stanCNRatio, CNnew, abgdHumificationFactor, xx2 * FRACTION_CARBON_PLANT);
            }
            if (xx3 > 0.0)
            {
                CNnew = CNdestiny (NMineralConcentration, flatCNRatio);
                nm3 = NitrogenMineralization (flatCNRatio, CNnew, abgdHumificationFactor, xx3 * FRACTION_CARBON_PLANT);
            }
            if (xx6 > 0.0)
            {
                CNnew = CNdestiny (NMineralConcentration, smanCNRatio);
                nm6 = NitrogenMineralization (smanCNRatio, CNnew, manuHumificationFactor, xx6);
            }

            if (nm2 < 0.0)
                nm0 += nm2;
            if (nm3 < 0.0)
                nm0 += nm3;
            if (nm6 < 0.0)
                nm0 += nm6;
        }

        if (xx1 > 0.0)
        {
            CNnew = CNdestiny (NMineralConcentration, abgdCNRatio);
            nm1 = NitrogenMineralization (abgdCNRatio, CNnew, abgdHumificationFactor, xx1 * FRACTION_CARBON_PLANT);
        }
        if (xx4 > 0.0)
        {
            CNnew = CNdestiny (NMineralConcentration, rootCNRatio);
            nm4 = NitrogenMineralization (rootCNRatio, CNnew, rootHumificationFactor, xx4 * FRACTION_CARBON_PLANT);
        }
        if (xx5 > 0.0)
        {
            CNnew = CNdestiny (NMineralConcentration, rhizCNRatio);
            nm5 = NitrogenMineralization (rhizCNRatio, CNnew, rhizHumificationFactor, xx5 * FRACTION_CARBON_RIZHO);
        }
        if (xx7 > 0.0)
        {
            CNnew = CNdestiny (NMineralConcentration, manuCNRatio);
            nm7 = NitrogenMineralization (manuCNRatio, CNnew, manuHumificationFactor, xx7);
        }

        if (nm1 < 0.0)
            nm0 += nm1;
        if (nm4 < 0.0)
            nm0 += nm4;
        if (nm5 < 0.0)
            nm0 += nm5;
        if (nm7 < 0.0)
            nm0 += nm7;

        if (-nm0 > NMineral)
            decompReductionFactor = NMineral / (-nm0);
        else
            decompReductionFactor = 1.0;

        if (decompReductionFactor < 1.0)
        {
            /* Adjust actual decomposition as a function on mineral N
             * availability */
            if (nm1 < 0.0)
                xx1 *= decompReductionFactor;
            if (nm2 < 0.0)
                xx2 *= decompReductionFactor;
            if (nm3 < 0.0)
                xx3 *= decompReductionFactor;
            if (nm4 < 0.0)
                xx4 *= decompReductionFactor;
            if (nm5 < 0.0)
                xx5 *= decompReductionFactor;
            if (nm6 < 0.0)
                xx6 *= decompReductionFactor;
            if (nm7 < 0.0)
                xx7 *= decompReductionFactor;

            /* Recalculate net mineralization only for pools with adjusted
             * decomposition */
            if (nm1 < 0.0)
            {
                CNnew = CNdestiny (NMineralConcentration, abgdCNRatio);
                nm1 = NitrogenMineralization (abgdCNRatio, CNnew, abgdHumificationFactor, xx1 * FRACTION_CARBON_PLANT);
            }
            if (nm2 < 0.0)
            {
                CNnew = CNdestiny (NMineralConcentration, stanCNRatio);
                nm2 = NitrogenMineralization (stanCNRatio, CNnew, abgdHumificationFactor, xx2 * FRACTION_CARBON_PLANT);
            }
            if (nm3 < 0.0)
            {
                CNnew = CNdestiny (NMineralConcentration, flatCNRatio);
                nm3 = NitrogenMineralization (flatCNRatio, CNnew, abgdHumificationFactor, xx3 * FRACTION_CARBON_PLANT);
            }
            if (nm4 < 0.0)
            {
                CNnew = CNdestiny (NMineralConcentration, rootCNRatio);
                nm4 = NitrogenMineralization (rootCNRatio, CNnew, rootHumificationFactor, xx4 * FRACTION_CARBON_PLANT);
            }
            if (nm5 < 0.0)
            {
                CNnew = CNdestiny (NMineralConcentration, rhizCNRatio);
                nm5 = NitrogenMineralization (rhizCNRatio, CNnew, rhizHumificationFactor, xx5 * FRACTION_CARBON_RIZHO);
            }
            if (nm6 < 0.0)
            {
                CNnew = CNdestiny (NMineralConcentration, smanCNRatio);
                nm6 = NitrogenMineralization (smanCNRatio, CNnew, manuHumificationFactor, xx6);
            }
            if (nm7 < 0.0)
            {
                CNnew = CNdestiny (NMineralConcentration, manuCNRatio);
                nm7 = NitrogenMineralization (manuCNRatio, CNnew, manuHumificationFactor, xx7);
            }
        }

        /* Calculate weighted C retention efficiency from microbial biomass to
         * soil organic matter */
        xxPartialSum = xx1 + xx2 + xx3 + xx4 + xx5 + xx6 + xx7;
        if (xxPartialSum > 0.0)
            micrHumificationFactor = micrHumificationFactor > (1.0 / xxPartialSum) * (abgdHumificationFactor * (xx1 + xx2 + xx3) + rootHumificationFactor * xx4 + rhizHumificationFactor * xx5 + manuHumificationFactor * (xx6 + xx7)) ? micrHumificationFactor : (1.0 / xxPartialSum) * (abgdHumificationFactor * (xx1 + xx2 + xx3) + rootHumificationFactor * xx4 + rhizHumificationFactor * xx5 + manuHumificationFactor * (xx6 + xx7));

        /* SOC decomposition and SON mineralization */
        decompositionAdjustmentBySOC = 1.0 - 1.0 / (1.0 + pow ((Soil->SOC_Conc[i] / satSOCConc) / 0.22, 3.0));
        decompositionAdjustmentBySOC = decompositionAdjustmentBySOC < 1.0 ? decompositionAdjustmentBySOC : 1.0;
        socDecompositionRate = SoilCarbon->factorComposite[i] * (1.0 + tillageFactor[i]) * MAXIMUM_UNDISTURBED_SOC_DECOMPOSITION_RATE * decompositionAdjustmentBySOC / (1.0 - pow (micrHumificationFactor, 2.0));

        xx8 = socDecompositionRate * Soil->SOC_Mass[i];
        if (xx8 > 0.0)
        {
            CNnew = CNdestiny (NMineralConcentration, somCNratio);
            nm8 = NitrogenMineralization (somCNratio, CNnew, micrHumificationFactor, xx8);
        }

        /* Microbial biomass decomposition and N mineralization */
	/* Acceleration of microbial * turnover if > 3% of SOC */
        aux1 = (Soil->MBC_Mass[i] / Soil->SOC_Mass[i]) / 0.03;
        aux2 = exp (10.0 * (aux1 - 1.0));

        /* Charlie's steady state km so that Cm = 0.03 of organic carbon
         * km = 0.97ks / 0.03(ex(1-Cs/Cx))1/2.
         * Also notice ks acceleration from apparent to actual turnover
         * ks/(1-e^2) */
        aux3 = aux2 * MAXIMUM_UNDISTURBED_SOC_DECOMPOSITION_RATE * decompositionAdjustmentBySOC / (1.0 - pow (micrHumificationFactor, 2.0)) * (0.97 / 0.03) * (1.0 / micrHumificationFactor);
        micrDecompositionRate = SoilCarbon->factorComposite[i] * (1.0 + tillageFactor[i]) * aux3;
        xx9 = micrDecompositionRate * Soil->MBC_Mass[i];
        if (xx9 > 0.0)
            nm9 = NitrogenMineralization (micrCNRatio, micrCNRatio, micrHumificationFactor, xx9);

        /* Calculate N removal from decomposing pools */
        if (xx1 > 0.0)
            nr1 = xx1 * FRACTION_CARBON_PLANT / abgdCNRatio;
        if (xx2 > 0.0)
            nr2 = xx2 * FRACTION_CARBON_PLANT / stanCNRatio;
        if (xx3 > 0.0)
            nr3 = xx3 * FRACTION_CARBON_PLANT / flatCNRatio;
        if (xx4 > 0.0)
            nr4 = xx4 * FRACTION_CARBON_PLANT / rootCNRatio;
        if (xx5 > 0.0)
            nr5 = xx5 * FRACTION_CARBON_RIZHO / rhizCNRatio;
        if (xx6 > 0.0)
            nr6 = xx6 / smanCNRatio;
        if (xx7 > 0.0)
            nr7 = xx7 / manuCNRatio;
        if (xx8 > 0.0)
            nr8 = xx8 / somCNratio;
        if (xx9 > 0.0)
            nr9 = xx9 / micrCNRatio;

        /* Calculate N contribution (N humification) to microbial pool of each
	 * decomposing pool */
        if (nm1 > 0.0)
            nh1 = nr1 - nm1;
        else
            nh1 = nr1;
        if (nm2 > 0.0)
            nh2 = nr2 - nm2;
        else
            nh2 = nr2;
        if (nm3 > 0.0)
            nh3 = nr3 - nm3;
        else
            nh3 = nr3;
        if (nm4 > 0.0)
            nh4 = nr4 - nm4;
        else
            nh4 = nr4;
        if (nm5 > 0.0)
            nh5 = nr5 - nm5;
        else
            nh5 = nr5;
        if (nm6 > 0.0)
            nh6 = nr6 - nm6;
        else
            nh6 = nr6;
        if (nm7 > 0.0)
            nh7 = nr7 - nm7;
        else
            nh7 = nr7;
        if (nm8 > 0.0)
            nh8 = nr8 - nm8;
        else
            nh8 = nr8;
        if (nm9 > 0.0)
            nh9 = nr9 - nm9;
        else
            nh9 = nr9;

        /* Calculate total residue, manure, and som carbon tansfer to
         * microbial pool */
        humifiedCarbon = abgdHumificationFactor * FRACTION_CARBON_PLANT * (xx1 + xx2 + xx3) + rootHumificationFactor * FRACTION_CARBON_PLANT * xx4 + rhizHumificationFactor * FRACTION_CARBON_RIZHO * xx5 + manuHumificationFactor * (xx6 + xx7) + micrHumificationFactor * xx8;

        /* Calculate total residue, manure, and som nitrogen transfer to
         * microbial pool */
        humifiedNitrogen = nh1 + nh2 + nh3 + nh4 + nh5 + nh6 + nh7 + nh8;

        /* Accumulate N mineralization, immobilization, and net
         * mineralization */
        if (nm1 > 0.0)
            NMineralization += nm1;
        else
            NImmobilization += nm1;
        if (nm2 > 0.0)
            NMineralization += nm2;
        else
            NImmobilization += nm2;
        if (nm3 > 0.0)
            NMineralization += nm3;
        else
            NImmobilization += nm3;
        if (nm4 > 0.0)
            NMineralization += nm4;
        else
            NImmobilization += nm4;
        if (nm5 > 0.0)
            NMineralization += nm5;
        else
            NImmobilization += nm5;
        if (nm6 > 0.0)
            NMineralization += nm6;
        else
            NImmobilization += nm6;
        if (nm7 > 0.0)
            NMineralization += nm7;
        else
            NImmobilization += nm7;
        if (nm8 > 0.0)
            NMineralization += nm8;
        else
            NImmobilization += nm8;
        if (nm9 > 0.0)
            NMineralization += nm9;
        else
            NImmobilization += nm9;

        NNetMineralization = NMineralization + NImmobilization;

        /* Update pools (N immobilization is negative) */
        Soil->NO3[i] += NImmobilization * (1.0 - NH4_Fraction);
        Soil->NH4[i] += NImmobilization * NH4_Fraction + NMineralization;

        Soil->NO3[i] = Soil->NO3[i] < 0.0 ? 0.0 : Soil->NO3[i];
        Soil->NH4[i] = Soil->NH4[i] < 0.0 ? 0.0 : Soil->NH4[i];

        if (i == 0)
        {
            Residue->stanResidueWater *= (1.0 - xx2 / (Residue->stanResidueMass + 1e-10));
            Residue->flatResidueWater *= (1.0 - xx3 / (Residue->flatResidueMass + 1e-10));
            Residue->stanResidueMass -= xx2;
            Residue->flatResidueMass -= xx3;
            Residue->manureSurfaceC -= xx6;
            Residue->stanResidueN -= nr2;
            Residue->flatResidueN -= nr3;
            Residue->manureSurfaceN -= nr6;
        }

        Residue->residueAbgd[i] -= xx1;
        Residue->residueRt[i] -= xx4;
        Residue->residueRz[i] -= xx5;
        Residue->manureC[i] -= xx7;
        Residue->residueAbgdN[i] -= nr1;
        Residue->residueRtN[i] -= nr4;
        Residue->residueRzN[i] -= nr5;
        Residue->manureN[i] -= nr7;

        Soil->SOC_Mass[i] += xx9 * micrHumificationFactor - xx8;
        Soil->SON_Mass[i] += nh9 - nr8;
        Soil->SOC_Conc[i] = Soil->SOC_Mass[i] * 1000.0 / soilMass;
        Soil->MBC_Mass[i] += humifiedCarbon - xx9;
        Soil->MBN_Mass[i] += humifiedNitrogen + (-NImmobilization) - nr9;
        SoilCarbon->carbonRespired[i] = (1.0 - abgdHumificationFactor) * FRACTION_CARBON_PLANT * (xx1 + xx2 + xx3) + (1.0 - rootHumificationFactor) * FRACTION_CARBON_PLANT * xx4 + (1.0 - rhizHumificationFactor) * FRACTION_CARBON_RIZHO * xx5 + (1.0 - manuHumificationFactor) * (xx6 + xx7) + (1.0 - micrHumificationFactor) * (xx8 + xx9);

        /* For output */
        SoilCarbon->annualSoilCarbonDecompositionRate[i] += socDecompositionRate;

        /* Excludes residue and manure */
        SoilCarbon->annualRespiredCarbonMass[i] += (1.0 - micrHumificationFactor) * (xx8 + xx9);

        /* Residue, roots and manure only */
        SoilCarbon->annualRespiredResidueCarbonMass[i] += SoilCarbon->carbonRespired[i] - (1.0 - micrHumificationFactor) * (xx8 + xx9);
        SoilCarbon->abgdCarbonInput[i] += FRACTION_CARBON_PLANT * (xx1 + xx2 + xx3);
        SoilCarbon->rootCarbonInput[i] += FRACTION_CARBON_PLANT * xx4;
        SoilCarbon->rhizCarbonInput[i] += FRACTION_CARBON_RIZHO * xx5;
        SoilCarbon->manuCarbonInput[i] += (xx6 + xx7);
        SoilCarbon->annualCarbonInputByLayer[i] = SoilCarbon->abgdCarbonInput[i] + SoilCarbon->rootCarbonInput[i] + SoilCarbon->rhizCarbonInput[i] + SoilCarbon->manuCarbonInput[i];

        /* C that goes to SOC and C that goes to microbial mass */
        SoilCarbon->annualHumifiedCarbonMass[i] += xx9 * micrHumificationFactor + humifiedCarbon;
        Soil->N_Immobilization += NImmobilization;
        Soil->N_Mineralization += NMineralization;
        Soil->N_NetMineralization += NNetMineralization;
        Soil->SOCProfile += Soil->SOC_Mass[i] + Soil->MBC_Mass[i];
        Soil->SONProfile += Soil->SON_Mass[i] + Soil->MBN_Mass[i];

        /* C that goes to SOC and C that goes to microbial mass */
        Soil->C_Humified += xx9 * micrHumificationFactor + humifiedCarbon;

        /* Residues, roots and manure */
        Soil->C_ResidueRespired += SoilCarbon->carbonRespired[i] - (1.0 - micrHumificationFactor) * (xx8 + xx9);
        Soil->C_SoilRespired += (1.0 - micrHumificationFactor) * (xx8 + xx9);

        SoilCarbon->annualNmineralization[i] += NMineralization * 1000.0;
        SoilCarbon->annualNImmobilization[i] += NImmobilization * 1000.0;
        SoilCarbon->annualNNetMineralization[i] += NNetMineralization * 1000.0;

        NFinal += Soil->SON_Mass[i] + Soil->MBN_Mass[i] + Soil->NO3[i] + Soil->NH4[i] + Residue->residueAbgdN[i] + Residue->residueRtN[i] + Residue->residueRzN[i] + Residue->manureN[i];
    }	/* End soil layer loop */

    NFinal += Residue->stanResidueN + Residue->flatResidueN + Residue->manureSurfaceN;

    SoilCarbon->annualAmmoniumNitrification += Soil->NH4_Nitrification * 1000.0;
    SoilCarbon->annualNitrousOxidefromNitrification += Soil->N2O_Nitrification * 1000.0;
    SoilCarbon->annualAmmoniaVolatilization += Soil->NH4_Volatilization * 1000.0;
    SoilCarbon->annualNO3Denitrification += Soil->NO3_Denitrification * 1000.0;
    SoilCarbon->annualNitrousOxidefromDenitrification += Soil->N2O_Denitrification * 1000.0;
    SoilCarbon->annualNitrateLeaching += Soil->NO3Leaching * 1000.0;
    SoilCarbon->annualAmmoniumLeaching += Soil->NH4Leaching * 1000.0;

    if (fabs (NFinal - NInitial) > 0.00001)
    {
        printf ("ERROR: Soil nitrogen balance error!\n");
        printf ("NInitial = %lf, NFinal = %lf\n", NInitial, NFinal);
	exit (1);
    }
}

void ComputeSoilCarbonBalance (SoilCarbonStruct *SoilCarbon, int y, ResidueStruct *Residue, SoilStruct *Soil, double *tillageFactor)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * i		    int		Loop counter
     * apparentSOCDecompositionRate
     *			    double
     * apparentSOCDecomposition
     *			    double
     * SON_Mineralization   double
     * humifiedCarbon	    double
     * humifiedNitrogen	    double
     * abgdHumificationFactor
     *			    double
     * rootHumificationFactor
     *			    double
     * rhizHumificationFactor
     *			    double
     * manuHumificationFactor
     *			    double
     * soilMass		    double
     * satSOCConc	    double	(g C/kg soil)
     * humificationAdjustmentBySOC
     *			    double
     * decompositionAdjustmentBySOC
     *			    double
     * contactFractionFlat  double	Fraction of surface residues subject
     *					  to decomposition
     * contactFractionStan  double	Fraction of surface residues subject
     *					  to decomposition
     * xx0		    double	Residue mass decomposition (Mg/ha/day)
     * xx1		    double	Residue mass decomposition (Mg/ha/day)
     * xx2		    double	Residue mass decomposition (Mg/ha/day)
     * xx3		    double	Residue mass decomposition (Mg/ha/day)
     * xx4		    double	Residue mass decomposition (Mg/ha/day)
     * xx5		    double	Residue mass decomposition (Mg/ha/day)
     * xx6		    double	Manure carbon decomposition
     *					  (Mg/ha/day)
     * xx7		    double	Manure carbon decomposition
     *					  (Mg/ha/day)
     * nm0		    double	Residue nitrogen net mineralization
     * nm1		    double	Residue nitrogen net mineralization
     * nm2		    double	Residue nitrogen net mineralization
     * nm3		    double	Residue nitrogen net mineralization
     * nm4		    double	Residue nitrogen net mineralization
     * nm5		    double	Residue nitrogen net mineralization
     * nm6		    double	Residue nitrogen net mineralization
     * nm7		    double	Residue nitrogen net mineralization
     * nm8		    double	Residue nitrogen net mineralization
     * nr0		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr1		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr2		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr3		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr4		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr5		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr6		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nr7		    double	Residue nitrogen removed by
     *					  decomposition from each residue pool
     * nh0		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh1		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh2		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh3		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh4		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh5		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh6		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * nh7		    double	Residue nitrogen transfered to
     *					  microbial biomass
     * NMineralization	    double
     * NImmobilization	    double
     * NNetMineralization   double
     * stanCNRatio	    double	CN ratio of standing residues
     * flatCNRatio	    double	CN ratio of flat residues
     * abgdCNRatio	    double	CN ratio of aboveground residues (in
     *					  each soil layer)
     * rootCNRatio	    double	CN ratio of root residues
     * rhizCNRatio	    double	CN ratio of rizhodeposition
     * smanCNRatio	    double	CN ratio of surface manure
     * manuCNRatio	    double	CN ratio of manure
     * CN_SOM		    double	sman = surface manure, manu = manure
     * NO3_Concentration    double	(g N-NO3/g soil)
     * NMineral		    double	Sum of NO3 and NH4 (be careful with
     *					  units
     * NH4_Fraction	    double	Fraction of NH4 in the sum NO3 + NH4
     * decompReductionFactor
     *			    double
     * NInitial		    double
     * NFinal		    double
     */
    int		    i;

    double	    apparentSOCDecompositionRate;
    double	    apparentSOCDecomposition;
    double	    SON_Mineralization;
    double	    humifiedCarbon;
    double	    humifiedNitrogen;

    double	    abgdHumificationFactor;
    double	    rootHumificationFactor;
    double	    rhizHumificationFactor;
    double	    manuHumificationFactor;
    double	    soilMass;
    double	    satSOCConc;

    double	    humificationAdjustmentBySOC;
    double	    decompositionAdjustmentBySOC;

    double	    contactFractionFlat;
    double	    contactFractionStan;

    double	    xx1, xx2, xx3, xx4, xx5, xx6, xx7;
    double	    nm1, nm2, nm3, nm4, nm5, nm6, nm7, nm8;
    double	    nh1, nh2, nh3, nh4, nh5, nh6, nh7;
    double	    nr1, nr2, nr3, nr4, nr5, nr6, nr7;

    double          NMineralization;
    double          NImmobilization;
    double          NNetMineralization;

    double          stanCNRatio;
    double          flatCNRatio;
    double          abgdCNRatio;
    double          rootCNRatio; 
    double          rhizCNRatio;
    double          smanCNRatio;
    double          manuCNRatio;
    double	    CN_SOM;

    double	    NO3_Concentration;
    double          NMineral;
    double          NH4_Fraction;
    double          decompReductionFactor;
    double          NInitial, NFinal;

    NInitial = 0.0;
    NFinal = 0.0;
    NInitial = Residue->stanResidueN + Residue->flatResidueN + Residue->manureSurfaceN;
    contactFractionFlat = 0.0;
    contactFractionStan = 0.0;
    Soil->N_Immobilization = 0.0;
    Soil->N_Mineralization = 0.0;
    Soil->N_NetMineralization = 0.0;
    Soil->SOCProfile = 0.0;
    Soil->SONProfile = 0.0;
    Soil->C_Humified = 0.0;
    Soil->C_ResidueRespired = 0.0;
    Soil->C_SoilRespired = 0.0;

    for (i = 0; i < Soil->totalLayers; i++)
    {
        NInitial += Soil->SON_Mass[i] + Soil->NO3[i] + Soil->NH4[i] + Residue->residueAbgdN[i] + Residue->residueRtN[i] + Residue->residueRzN[i] + Residue->manureN[i];

        apparentSOCDecompositionRate = 0.0;
        apparentSOCDecomposition = 0.0;
        SON_Mineralization = 0.0;
        humifiedCarbon = 0.0;
        humifiedNitrogen = 0.0;

        NNetMineralization = 0.0;
        NMineralization = 0.0;
        NImmobilization = 0.0;

        stanCNRatio = 0.0;
        flatCNRatio = 0.0;
        abgdCNRatio = 0.0;
        rootCNRatio = 0.0;
        rhizCNRatio = 0.0;
        smanCNRatio = 0.0;
        manuCNRatio = 0.0;

        xx1 = xx2 = xx3 = xx4 = xx5 = xx6 = xx7 = 0.0;
        nm1 = nm2 = nm3 = nm4 = nm5 = nm6 = nm7 = nm8 = 0.0;
        nr1 = nr2 = nr3 = nr4 = nr5 = nr6 = nr7 = 0.0;
        nh1 = nh2 = nh3 = nh4 = nh5 = nh6 = nh7 = 0.0;

	/* 10000 converts from Mg/m2 to Mg/ha */
        soilMass = 10000.0 * Soil->BD[i] * Soil->layerThickness[i];
	NO3_Concentration = Soil->NO3[i] / soilMass;
        satSOCConc = 21.1 + 0.375 * Soil->Clay[i] * 100.0;

        /* Compute C/N ratios */
        if (i == 0)
        {
            if (Residue->stanResidueMass > 0.0)
                stanCNRatio = Residue->stanResidueMass * FRACTION_CARBON_PLANT / Residue->stanResidueN;
            if (Residue->flatResidueMass > 0.0)
                flatCNRatio = Residue->flatResidueMass * FRACTION_CARBON_PLANT / Residue->flatResidueN;
            if (Residue->manureSurfaceC > 0.0)
                smanCNRatio = Residue->manureSurfaceC / Residue->manureSurfaceN;
        }

        CN_SOM = Soil->SOC_Mass[i] / Soil->SON_Mass[i];

        if (Residue->residueAbgd[i] > 0.0)
            abgdCNRatio = Residue->residueAbgd[i] * FRACTION_CARBON_PLANT / Residue->residueAbgdN[i];
        if (Residue->residueRt[i] > 0.0)
            rootCNRatio = Residue->residueRt[i] * FRACTION_CARBON_PLANT / Residue->residueRtN[i];
        if (Residue->residueRz[i] > 0.0)
            rhizCNRatio = Residue->residueRz[i] * FRACTION_CARBON_RIZHO / Residue->residueRzN[i];
        if (Residue->manureC[i] > 0.0)
            manuCNRatio = Residue->manureC[i] / Residue->manureN[i];

        /* Humification */
        /* Humification reduction when C conc approaches saturation */
        humificationAdjustmentBySOC = 1.0 - pow (Soil->SOC_Conc[i] / satSOCConc, SOC_HUMIFICATION_POWER);
        humificationAdjustmentBySOC = humificationAdjustmentBySOC > 0.0 ? humificationAdjustmentBySOC : 0.0;

        abgdHumificationFactor = MaximumAbgdHumificationFactor (Soil->Clay[i]) * humificationAdjustmentBySOC;
        rootHumificationFactor = MaximumRootHumificationFactor (Soil->Clay[i]) * humificationAdjustmentBySOC;
        rhizHumificationFactor = MaximumRhizHumificationFactor (Soil->Clay[i]) * humificationAdjustmentBySOC;
        manuHumificationFactor = MaximumManuHumificationFactor (Soil->Clay[i]) * humificationAdjustmentBySOC;

        /* Residue and manure decomposition */
        if (i == 0)
        {
            contactFractionStan = pow (Residue->stanResidueTau, exp (-1.5 / sqrt (1.0 - Residue->stanResidueTau)));
            contactFractionFlat = pow (Residue->flatResidueTau, exp (-1.5 / sqrt (1.0 - Residue->flatResidueTau)));

            xx2 = SoilCarbon->factorComposite[i] * MAXIMUM_RESIDUE_DECOMPOSITION_RATE * contactFractionStan * Residue->stanResidueMass;
            xx3 = SoilCarbon->factorComposite[i] * MAXIMUM_RESIDUE_DECOMPOSITION_RATE * contactFractionFlat * Residue->flatResidueMass;
            xx6 = SoilCarbon->factorComposite[i] * MAXIMUM_MANURE_DECOMPOSITION_RATE * Residue->manureSurfaceC;
        }

        xx1 = SoilCarbon->factorComposite[i] * MAXIMUM_RESIDUE_DECOMPOSITION_RATE * Residue->residueAbgd[i];
        xx4 = SoilCarbon->factorComposite[i] * MAXIMUM_ROOT_DECOMPOSITION_RATE * Residue->residueRt[i];
        xx5 = SoilCarbon->factorComposite[i] * MAXIMUM_RHIZO_DECOMPOSITION_RATE * Residue->residueRz[i];
        xx7 = SoilCarbon->factorComposite[i] * MAXIMUM_MANURE_DECOMPOSITION_RATE * Residue->manureC[i];

        /* Inorganic N limitation for decomposition */
        /* If decomposition > 0 then compute net N mineralization and
	 * accumulate negatives */
	if (i == 0)
	{
	    if (xx2 > 0.0)
		nm2 = PoolNitrogenMineralization (NO3_Concentration, stanCNRatio, abgdHumificationFactor, xx2, FRACTION_CARBON_PLANT);
            if (xx3 > 0.0)
		nm3 = PoolNitrogenMineralization (NO3_Concentration, flatCNRatio, abgdHumificationFactor, xx3, FRACTION_CARBON_PLANT);
            if (xx6 > 0.0)
		nm6 = PoolNitrogenMineralization (NO3_Concentration, smanCNRatio, manuHumificationFactor, xx6, 1.0);

	    if (nm2 < 0.0)
		nm8 += nm2;
            if (nm3 < 0.0)
		nm8 += nm3;
            if (nm6 < 0.0)
		nm8 += nm6;
	}

        if (xx1 > 0.0)
	    nm1 = PoolNitrogenMineralization (NO3_Concentration, abgdCNRatio, abgdHumificationFactor, xx1, FRACTION_CARBON_PLANT);
        if (xx4 > 0.0)
	    nm4 = PoolNitrogenMineralization (NO3_Concentration, rootCNRatio, rootHumificationFactor, xx4, FRACTION_CARBON_PLANT);
        if (xx5 > 0.0)
	    nm5 = PoolNitrogenMineralization (NO3_Concentration, rhizCNRatio, rhizHumificationFactor, xx5, FRACTION_CARBON_RIZHO);
        if (xx7 > 0.0)
	    nm7 = PoolNitrogenMineralization (NO3_Concentration, manuCNRatio, manuHumificationFactor, xx7, 1.0);

        if (nm1 < 0.0)
	    nm8 += nm1;
        if (nm4 < 0.0)
	    nm8 += nm4;
        if (nm5 < 0.0)
	    nm8 += nm5;
        if (nm7 < 0.0)
	    nm8 += nm7;

        NMineral = Soil->NO3[i] + Soil->NH4[i];
        NH4_Fraction = Soil->NH4[i] / NMineral;
        if (-nm8 > NMineral)
	    decompReductionFactor = NMineral / (-nm8);
        else
	    decompReductionFactor = 1.0;

        if (decompReductionFactor < 1.0)
	{
	    /* Adjust actual decomposition as a function on mineral N
	     * availability */
            if (nm1 < 0.0)
		xx1 *= decompReductionFactor;
            if (nm2 < 0.0)
		xx2 *= decompReductionFactor;
            if (nm3 < 0.0)
		xx3 *= decompReductionFactor;
            if (nm4 < 0.0)
		xx4 *= decompReductionFactor;
            if (nm5 < 0.0)
		xx5 *= decompReductionFactor;
            if (nm6 < 0.0)
		xx6 *= decompReductionFactor;
            if (nm7 < 0.0)
		xx7 *= decompReductionFactor;

            /* Recalculate net mineralization only for pools with adjusted
	     * decomposition */
            if (nm1 < 0.0)
		nm1 = PoolNitrogenMineralization (NO3_Concentration, abgdCNRatio, abgdHumificationFactor, xx1, FRACTION_CARBON_PLANT);
            if (nm2 < 0.0)
		nm2 = PoolNitrogenMineralization (NO3_Concentration, stanCNRatio, abgdHumificationFactor, xx2, FRACTION_CARBON_PLANT);
            if (nm3 < 0.0)
		nm3 = PoolNitrogenMineralization (NO3_Concentration, flatCNRatio, abgdHumificationFactor, xx3, FRACTION_CARBON_PLANT);
            if (nm4 < 0.0)
		nm4 = PoolNitrogenMineralization (NO3_Concentration, rootCNRatio, rootHumificationFactor, xx4, FRACTION_CARBON_PLANT);
            if (nm5 < 0.0)
		nm5 = PoolNitrogenMineralization (NO3_Concentration, rhizCNRatio, rhizHumificationFactor, xx5, FRACTION_CARBON_RIZHO);
            if (nm6 < 0.0)
		nm6 = PoolNitrogenMineralization (NO3_Concentration, smanCNRatio, manuHumificationFactor, xx6, 1.0);
            if (nm7 < 0.0)
		nm7 = PoolNitrogenMineralization (NO3_Concentration, manuCNRatio, manuHumificationFactor, xx7, 1.0);
	}

        /* Calculate N removal from decomposing pools */
        if (xx1 > 0.0)
            nr1 = xx1 * FRACTION_CARBON_PLANT / abgdCNRatio;
        if (xx2 > 0.0)
            nr2 = xx2 * FRACTION_CARBON_PLANT / stanCNRatio;
        if (xx3 > 0.0)
            nr3 = xx3 * FRACTION_CARBON_PLANT / flatCNRatio;
        if (xx4 > 0.0)
            nr4 = xx4 * FRACTION_CARBON_PLANT / rootCNRatio;
        if (xx5 > 0.0)
            nr5 = xx5 * FRACTION_CARBON_RIZHO / rhizCNRatio;
        if (xx6 > 0.0)
            nr6 = xx6 / smanCNRatio;
        if (xx7 > 0.0)
            nr7 = xx7 / manuCNRatio;

        /* Calculate N contribution (N humification) to microbial pool of each
	 * decomposing pool */
        if (nm1 > 0.0)
            nh1 = nr1 - nm1;
        else
            nh1 = nr1;
        if (nm2 > 0.0)
            nh2 = nr2 - nm2;
        else
            nh2 = nr2;
        if (nm3 > 0.0)
            nh3 = nr3 - nm3;
        else
            nh3 = nr3;
        if (nm4 > 0.0)
            nh4 = nr4 - nm4;
        else
            nh4 = nr4;
        if (nm5 > 0.0)
            nh5 = nr5 - nm5;
        else
            nh5 = nr5;
        if (nm6 > 0.0)
            nh6 = nr6 - nm6;
        else
            nh6 = nr6;
        if (nm7 > 0.0)
            nh7 = nr7 - nm7;
        else
            nh7 = nr7;

        /* Calculate total residue and manure carbon humification */
        humifiedCarbon = abgdHumificationFactor * FRACTION_CARBON_PLANT * (xx1 + xx2 + xx3) + rootHumificationFactor * FRACTION_CARBON_PLANT * xx4 + rhizHumificationFactor * FRACTION_CARBON_RIZHO * xx5 + manuHumificationFactor * (xx6 + xx7);

        /* Calculate total residue and manure nitrogen humification */
        humifiedNitrogen = nh1 + nh2 + nh3 + nh4 + nh5 + nh6 + nh7;

        /* SOC DECOMPOSITION */
        decompositionAdjustmentBySOC = 1.0 - 1.0 / (1.0 + pow (Soil->SOC_Conc[i] / satSOCConc / 0.22, 3.0));
        if (decompositionAdjustmentBySOC > 1.0)
	    decompositionAdjustmentBySOC = 1.0;
        apparentSOCDecompositionRate = SoilCarbon->factorComposite[i] * (1.0 + tillageFactor[i]) * MAXIMUM_UNDISTURBED_SOC_DECOMPOSITION_RATE * decompositionAdjustmentBySOC;
        apparentSOCDecomposition = apparentSOCDecompositionRate * Soil->SOC_Mass[i];
        SON_Mineralization = apparentSOCDecomposition / CN_SOM;

        /* Accumulate N mineralization, immobilization, and net
	 * mineralization */
        if (nm1 > 0.0)
	    NMineralization += nm1;
	else
	    NImmobilization += nm1;
	if (nm2 > 0.0)
	    NMineralization += nm2;
	else
	    NImmobilization += nm2;
        if (nm3 > 0.0)
	    NMineralization += nm3;
	else
	    NImmobilization += nm3;
        if (nm4 > 0.0)
	    NMineralization += nm4;
	else
	    NImmobilization += nm4;
        if (nm5 > 0.0)
	    NMineralization += nm5;
	else
	    NImmobilization += nm5;
        if (nm6 > 0.0)
	    NMineralization += nm6;
	else
	    NImmobilization += nm6;
	if (nm7 > 0.0)
	    NMineralization += nm7;
	else
	    NImmobilization += nm7;

        NMineralization += SON_Mineralization;
        NNetMineralization = NMineralization + NImmobilization;

        /* Update pools (N immbilization is negative) */
        Soil->NO3[i] += NImmobilization * (1.0 - NH4_Fraction);
        Soil->NH4[i] += NImmobilization * NH4_Fraction + NMineralization;

        if (i == 0)
        {
            Residue->stanResidueWater *= (1.0 - xx2 / (Residue->stanResidueMass + 1e-10));
            Residue->flatResidueWater *= (1.0 - xx3 / (Residue->flatResidueMass + 1e-10));
            Residue->stanResidueMass -= xx2;
            Residue->flatResidueMass -= xx3;
            Residue->manureSurfaceC -= xx6;
            Residue->stanResidueN -= nr2;
            Residue->flatResidueN -= nr3;
            Residue->manureSurfaceN -= nr6;
        }

        Residue->residueAbgd[i] -= xx1;
        Residue->residueRt[i] -= xx4;
        Residue->residueRz[i] -= xx5;
        Residue->manureC[i] -= xx7;
        Residue->residueAbgdN[i] -= nr1;
        Residue->residueRtN[i] -= nr4;
        Residue->residueRzN[i] -= nr5;
        Residue->manureN[i] -= nr7;

        Soil->SOC_Mass[i] += humifiedCarbon - apparentSOCDecomposition;
        Soil->SON_Mass[i] += humifiedNitrogen + (-NImmobilization) - SON_Mineralization;
        Soil->SOC_Conc[i] = Soil->SOC_Mass[i] * 1000.0 / soilMass;
        SoilCarbon->carbonRespired[i] = apparentSOCDecomposition + (1.0 - abgdHumificationFactor) * FRACTION_CARBON_PLANT * (xx1 + xx2 + xx3) + (1.0 - rootHumificationFactor) * FRACTION_CARBON_PLANT * xx4 + (1.0 - rhizHumificationFactor) * FRACTION_CARBON_RIZHO * xx5 + (1.0 - manuHumificationFactor) * (xx6 + xx7);

        /* For output */
        SoilCarbon->annualSoilCarbonDecompositionRate[i] += apparentSOCDecompositionRate;
	/* Excludes residue and manure */
        SoilCarbon->annualRespiredCarbonMass[i] += apparentSOCDecomposition;
	/* Residue and manure only */
        SoilCarbon->annualRespiredResidueCarbonMass[i] += SoilCarbon->carbonRespired[i] - apparentSOCDecomposition;
        SoilCarbon->abgdCarbonInput[i] += FRACTION_CARBON_PLANT * (xx1 + xx2 + xx3);
        SoilCarbon->rootCarbonInput[i] += FRACTION_CARBON_PLANT * xx4;
        SoilCarbon->rhizCarbonInput[i] += FRACTION_CARBON_RIZHO * xx5;
        SoilCarbon->manuCarbonInput[i] += (xx6 + xx7);

	SoilCarbon->annualCarbonInputByLayer[i] = SoilCarbon->abgdCarbonInput[i] + SoilCarbon->rootCarbonInput[i] + SoilCarbon->rhizCarbonInput[i] + SoilCarbon->manuCarbonInput[i];
	SoilCarbon->annualHumifiedCarbonMass[i] += humifiedCarbon;

	Soil->N_Immobilization += NImmobilization;
	Soil->N_Mineralization += NMineralization;
	Soil->N_NetMineralization += NNetMineralization;
	Soil->SOCProfile += Soil->SOC_Mass[i];
	Soil->SONProfile += Soil->SON_Mass[i];
	Soil->C_Humified += humifiedCarbon;
	/* Includes residue, roots, and manure */
	Soil->C_ResidueRespired += (SoilCarbon->carbonRespired[i] - apparentSOCDecomposition);
	Soil->C_SoilRespired += apparentSOCDecomposition;

	NFinal += Soil->SON_Mass[i] + Soil->NO3[i] + Soil->NH4[i] + Residue->residueAbgdN[i] + Residue->residueRtN[i] + Residue->residueRzN[i] + Residue->manureN[i];
    }	/* End soil layer loop */

    NFinal += Residue->stanResidueN + Residue->flatResidueN + Residue->manureSurfaceN;

    if (fabs(NFinal - NInitial) > 0.00001)
    {
        printf ("ERROR: Soil nitrogen balance error!\n");
        printf ("NInitial = %lf, NFinal = %lf\n", NInitial, NFinal);
	exit (1);
    }
}

void StoreOutput (SoilCarbonStruct *SoilCarbon, int y, int totalLayers, double *SOCMass)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * i		    int		Loop counter
     */
    int             i;

    for (i = 0; i < totalLayers; i++)
        SoilCarbon->carbonMassFinal[i] = SOCMass[i];
}

double Aeration (double AC)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * A1		    double
     * A2		    double
     */
    /* AC = soil air content */
    const double    A1 = 0.05;
    const double    A2 = 4.0;

    return (1.0 - 0.6 / (1.0 + pow (AC / A1, A2)));
}

double Moisture (double wp)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * M1		    double
     * M2		    double
     */
    /* WP = soil water potential */
    const double    M1 = -600.0;
    const double    M2 = 3.0;

    return (1.0 / (1.0 + pow (wp / M1, M2)));
}

double TemperatureFunction (double T)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * temp		    double
     * Q		    double
     * tMin		    double
     * tOpt		    double
     * tMax		    double
     */
    double          temp;
    double          Q;
    const double    tMin = -5.0;
    const double    tOpt = 35.0;
    const double    tMax = 50.0;

    if (T < 0.0 || T > tMax)
        temp = 0.0;
    else
    {
        Q = (tMin - tOpt) / (tOpt - tMax);
        temp = (pow (T - tMin, Q) * (tMax - T)) / (pow (tOpt - tMin, Q) * (tMax - tOpt));
        if (temp > 1.0)
            temp = 1.0;
    }

    return (temp);
}

double MaximumAbgdHumificationFactor (double clayFraction)
{
    return (0.092 + 0.104 * (1.0 - exp (-5.5 * clayFraction)));
}

double MaximumRootHumificationFactor (double clayFraction)
{
    return (0.092 + 0.104 * (1.0 - exp (-5.5 * clayFraction)));
}

double MaximumRhizHumificationFactor (double clayFraction)
{
    return (0.0 + 0.084 * (1.0 - exp (-5.5 * clayFraction)));
}

double MaximumManuHumificationFactor (double clayFraction)
{
    return (0.15 + 0.25 * (1.0 - exp (-5.5 * clayFraction)));
}

double NitrogenMineralization (double CNDecomposing, double CNnew, double humRate, double decomposedMass)
{
    return (decomposedMass * (1.0 / CNDecomposing - humRate / CNnew));
}

double CNdestiny (double NmineralConc, double CNdecomposing)
{

    /* Returns CN ratio of newly formed microbial biomass based on CN or
     * decomposing residue and N mineral in soil. Same function that for one-
     * pool model, but applied to microbial biomass
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * YY2		    double
     */
    double          YY2;

    NmineralConc = NmineralConc + 0.0000001;
    YY2 = 5.5 * (1.0 - 1.0 / (1.0 + pow (CNdecomposing / 110.0, 3.0)));

    return (8.5 + YY2 * (0.5 + 0.5 / (1.0 + pow (NmineralConc / 0.000008, 3.0))));
}

double PoolNitrogenMineralization (double NmineralConc, double CNRatioDecomposing, double humRate, double decomposedMass, double carbonConc)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * newCN		    double      CN of new organic matter (humified
     *					  residue)
     */
    double          newCN;

    decomposedMass *= carbonConc;
    newCN = Function_CNnew (NmineralConc, CNRatioDecomposing);

    return (decomposedMass * (1.0 / CNRatioDecomposing - humRate / newCN));
}

double Function_CNnew (double NmineralConc, double CNDecomposingPool)
{

    /* Returns CN ratio of newly formed organic matter based on CN or
     * decomposing residue and N mineral in soil
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * YY2		    double
     */
    double          YY2;

    NmineralConc = NmineralConc + 0.0000001;
    YY2 = 5.5 * (1.0 - 1.0 / (1.0 + pow (CNDecomposingPool / 110.0, 3.0)));

    return (8.5 + YY2 * (0.5 + 0.5 / (1.0 + pow (NmineralConc / 0.000008, 3.0))));
}
