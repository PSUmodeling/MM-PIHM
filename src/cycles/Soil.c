#include "Cycles.h"

void InitializeSoil (SoilStruct *Soil, WeatherStruct *Weather, SimControlStruct *SimControl)
{
    /*
     * 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * sb		    double	Saxton's b coefficient
     * sAP		    double	Saxton's air entry potential
     * sBD		    double	Saxton's bulk density
     * s33		    double	Saxton's volumetric WC at 33 J/kg
     * s1500		    double	Saxton's volumetric WC at 1500 J/kg
     * rr		    double
     * i		    int		Loop counter
     */
    double          sb;
    double          sAP;
    double          sBD;
    double          s33;
    double          s1500;
    double          rr;

    int             i;

    /* Initialize soil vairables */
    Soil->nodeDepth = (double *)malloc ((Soil->totalLayers + 1) * sizeof (double));
    Soil->cumulativeDepth = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->waterContent = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->Porosity = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->PAW = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->FC_WaterPotential = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->airEntryPotential = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->B_Value = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->M_Value = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->SOC_Conc = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->SOC_Mass = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->SON_Mass = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->MBC_Mass = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->MBN_Mass = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->waterUptake = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->pH = (double *)malloc ((Soil->totalLayers) * sizeof (double));
    Soil->n2o = (double *)malloc ((Soil->totalLayers) * sizeof (double));

    for (i = 0; i < Soil->totalLayers; i++)
        Soil->cumulativeDepth[i] = 0.0;

    Soil->cumulativeDepth[0] = Soil->layerThickness[0];
    Soil->nodeDepth[0] = 0.5 * Soil->layerThickness[0];

    for (i = 0; i < Soil->totalLayers; i++)
    {
        Soil->Clay[i] = Soil->Clay[i] / 100.0;
        Soil->Sand[i] = Soil->Sand[i] / 100.0;
        Soil->IOM[i] = Soil->IOM[i];
        Soil->NO3[i] = Soil->NO3[i] / 1000.0;
        Soil->NH4[i] = Soil->NH4[i] / 1000.0;
        if (i > 0)
        {
            Soil->cumulativeDepth[i] = Soil->cumulativeDepth[i - 1] + Soil->layerThickness[i];
            Soil->nodeDepth[i] = Soil->cumulativeDepth[i - 1] + Soil->layerThickness[i] / 2.0;
        }
    }

    Soil->nodeDepth[Soil->totalLayers] = Soil->cumulativeDepth[Soil->totalLayers - 1] + Soil->layerThickness[Soil->totalLayers - 1] / 2.0;

    /* Compute hydraulic properties */
    for (i = 0; i < Soil->totalLayers; i++)
    {
        sBD = BulkDensity (Soil->Clay[i], Soil->Sand[i], Soil->IOM[i]);
        s33 = VolumetricWCAt33Jkg (Soil->Clay[i], Soil->Sand[i], Soil->IOM[i]);
        s1500 = VolumetricWCAt1500Jkg (Soil->Clay[i], Soil->Sand[i], Soil->IOM[i]);
        sb = (log (1500.0) - log (33.0)) / (log (s33) - log (s1500));
        sAP = -33.0 * pow (s33 / (1.0 - sBD / 2.65), sb);

	/* Buld Density switch */
        if ((int)Soil->BD[i] != (int)BADVAL)  
            Soil->airEntryPotential[i] = sAP * pow (Soil->BD[i] / sBD, 0.67 * sb);
        else
        {
            Soil->BD[i] = sBD;
            Soil->airEntryPotential[i] = sAP;
        }

        Soil->Porosity[i] = 1.0 - Soil->BD[i] / 2.65;
        Soil->B_Value[i] = sb;
        Soil->M_Value[i] = 2.0 * Soil->B_Value[i] + 3.0;

	/* Field Capacity switch */
        if ((int)Soil->FC[i] == (int)BADVAL)
        {
            Soil->FC_WaterPotential[i] = -0.35088 * Soil->Clay[i] * 100.0 - 28.947;
            Soil->FC[i] = SoilWaterContent (Soil->Porosity[i], Soil->airEntryPotential[i], Soil->B_Value[i], Soil->FC_WaterPotential[i]);
        }
        else
            Soil->FC_WaterPotential[i] = SoilWaterPotential (Soil->Porosity[i], Soil->airEntryPotential[i], Soil->B_Value[i], Soil->FC[i]);

	/* Permanent Wilting Point switch */
        if ((int)Soil->PWP[i] == (int)BADVAL)
            Soil->PWP[i] = SoilWaterContent (Soil->Porosity[i], Soil->airEntryPotential[i], Soil->B_Value[i], -1500.0);

        if (Soil->PWP[i] >= Soil->FC[i])
        {
            printf ("ERROR: Permanent Wilting Point must be less than Field Capacity in Layer %d\n", i + 1);
            exit (1);
        }

        if (Soil->Porosity[i] < Soil->FC[i])
        {
            /* In case the air entry potential adjustment fails to deliver
             * FC < Porosity, this function estimates the ratio of volume
             * (sat - fc) / (fc - pwp) used to recalculate FC if bulk density
             * is too high (source: guess based on Saxton and Rawls data */
            rr = 0.63 + 0.59 * exp (0.6 * Soil->Sand[i] * 100.0 - 2.84);
            Soil->FC[i] = (rr * Soil->PWP[i] + Soil->Porosity[i]) / (1.0 + rr);
        }
    }

    /* Initialize variables depending on previous loop */
    for (i = 0; i < Soil->totalLayers; i++)
    {
        Soil->SOC_Conc[i] = Soil->IOM[i] * 10.0 * 0.58;
        Soil->SOC_Mass[i] = Soil->IOM[i] / 100.0 * 0.58 * Soil->layerThickness[i] * Soil->BD[i] * 10000.0;
	/* Initializes with CN ratio = 10 */
        Soil->SON_Mass[i] = Soil->SOC_Mass[i] / 10.0;
	/* Initializes as 3% of SOC_Mass but "added" C */
        Soil->MBC_Mass[i] = 0.03 * Soil->SOC_Mass[i];
	/* Initializes with CN ratio = 10 */
        Soil->MBN_Mass[i] = Soil->MBC_Mass[i] / 10.0;
        Soil->PAW[i] = Soil->FC[i] - Soil->PWP[i];
        Soil->waterContent[i] = (Soil->FC[i] + Soil->PWP[i]) / 2.0;
    }

    /* Initializes soil temperature in first day of simulation */
    Soil->dampingDepth = 2.0;

    if (Weather->siteLatitude >= 0.0)
        Soil->annualTemperaturePhase = 100.0;
    else
        Soil->annualTemperaturePhase = 280.0;

    Soil->soilTemperature = (double *)malloc ((Soil->totalLayers + 1) * sizeof (double));

    for (i = 0; i < Soil->totalLayers + 1; i++)
        Soil->soilTemperature[i] = EstimatedSoilTemperature (Soil->nodeDepth[i], 1, Weather->annualAverageTemperature[0], Weather->yearlyAmplitude[0], Soil->annualTemperaturePhase, Soil->dampingDepth);

    if (debug_mode)
    {
        printf ("\n*Soil properties after initialization:\n");
        printf ("\n*%-30s\t", "Clay");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->Clay[i]);
        printf ("\n*%-30s\t", "Sand");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->Sand[i]);
        printf ("\n*%-30s\t", "IOM");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->IOM[i]);
        printf ("\n*%-30s\t", "NO3");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->NO3[i]);
        printf ("\n*%-30s\t", "NH4");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->NH4[i]);
        printf ("\n*%-30s\t", "Node Depth");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->nodeDepth[i]);
        printf ("\n*%-30s\t", "Cumulative Depth");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->cumulativeDepth[i]);
        printf ("\n*%-30s\t", "Bulk density");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->BD[i]);
        printf ("\n*%-30s\t", "Porosity");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->Porosity[i]);
        printf ("\n*%-30s\t", "B");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->B_Value[i]);
        printf ("\n*%-30s\t", "Air Entry Potential");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->airEntryPotential[i]);
        printf ("\n*%-30s\t", "M");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->M_Value[i]);
        printf ("\n*%-30s\t", "Field capacity potential");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->FC_WaterPotential[i]);
        printf ("\n*%-30s\t", "Field capacity");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->FC[i]);
        printf ("\n*%-30s\t", "Wilting Point");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->PWP[i]);
        printf ("\n*%-30s\t", "SOC concentration");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->SOC_Conc[i]);
        printf ("\n*%-30s\t", "SOC mass");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->SOC_Mass[i]);
        printf ("\n*%-30s\t", "SON mass");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->SON_Mass[i]);
        printf ("\n*%-30s\t", "MBC mass");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->MBC_Mass[i]);
        printf ("\n*%-30s\t", "MBN mass");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->MBN_Mass[i]);
        printf ("\n*%-30s\t", "PAW");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->PAW[i]);
        printf ("\n*%-30s\t", "Water content");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->waterContent[i]);
        printf ("\n*%-30s\t", "Initial soil temp");
        for (i = 0; i < Soil->totalLayers; i++)
            printf ("%-10.3lf\t", Soil->soilTemperature[i]);
	printf ("\n");
    }
}

double SoilWaterPotential (double SaturationWC, double AirEntryPot, double Campbell_b, double WC)
{
    /*
     * 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * swp		    double	Soil water potential [return value]
     */
    double          swp;

    swp = AirEntryPot * pow (WC / SaturationWC, -Campbell_b);

    return (swp);
}

double VolumetricWCAt33Jkg (double Clay, double Sand, double OM)
{
    /*
     * Calculate volumetric water content at 33 Jkg
     * Saxton and Rawls 2006 SSSAJ 70:1569-1578 Eq 2 (r2 = 0.63)
     * Clay and sand fractional, OM as %
     * (original paper says % for everything, results make no sense)
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * x1		    double
     * vwc		    double	Volumetric soil water content [return
     *					  value]
     */
    double          x1;
    double          vwc;

    x1 = 0.299 - 0.251 * Sand + 0.195 * Clay + 0.011 * OM + 0.006 * Sand * OM - 0.027 * Clay * OM + 0.452 * Sand * Clay;
    vwc = -0.015 + 0.636 * x1 + 1.283 * pow (x1, 2.0);

    return (vwc);
}

double VolumetricWCAt1500Jkg (double Clay, double Sand, double OM)
{
    /*
     * Calculate volumetric water content at 33 Jkg
     * Saxton and Rawls 2006 SSSAJ 70:1569-1578 Eq 1 (r2 = 0.86)
     * Clay and sand fractional, OM as %
     * (original paper says % for everything, results make no sense)
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * x1		    double
     * vwc		    double	Volumetric soil water content [return
     *					  value]
     */
    double          x1;
    double          vwc;

    x1 = 0.031 - 0.024 * Sand + 0.487 * Clay + 0.006 * OM + 0.005 * Sand * OM - 0.013 * Clay * OM + 0.068 * Sand * Clay;
    vwc = -0.02 + 1.14 * x1;

    return (vwc);
}

double SoilWaterContent (double SaturationWC, double AirEntryPot, double Campbell_b, double Water_Potential)
{
    /* 
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * swc		    double	Soil water content [return value]
     */
    double          swc;

    swc = SaturationWC * pow (Water_Potential / AirEntryPot, -1.0 / Campbell_b);

    return (swc);
}

double BulkDensity (double Clay, double Sand, double OM)
{
    /*
     * Saxton and Rawls 2006 SSSAJ 70:1569-1578 Eq 6,5 (r2 = 0.29)
     * really poor fit
     * Clay and sand fractional, OM as %
     * (original paper says % for everything, results make no sense)
     * Note: X2 is Eq 3, supposedly representing moisture from FC to
     * saturation;
     * however, porosity is further adjusted by sand, an inconsistency
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * x1		    double
     * x2		    double
     * FC		    double
     * Porosity		    double
     * bd		    double	Bulk density [return value]
     */
    double          x1;
    double          x2;
    double          FC;
    double          Porosity;
    double          bd;

    x1 = 0.078 + 0.278 * Sand + 0.034 * Clay + 0.022 * OM - 0.018 * Sand * OM - 0.027 * Clay * OM - 0.584 * Sand * Clay;
    x2 = -0.107 + 1.636 * x1;
    FC = VolumetricWCAt33Jkg (Clay, Sand, OM);
    Porosity = 0.043 + FC + x2 - 0.097 * Sand;
    bd = (1.0 - Porosity) * 2.65;

    return (bd);
}
