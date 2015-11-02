#include "Cycles.h"

int             verbose_mode;
int             debug_mode;

int main (int argc, char *argv[])
{
    /*
     * Cycles main function
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * rotationYear	    int		Rotation year
     * nextSeedingDate	    int
     * nextSeedingYear	    int
     * y		    int
     * doy		    int
     * i		    int
     * c		    int
     * begin_t		    time_t	Time Cycles simulation begins
     * end_t		    time_t	Time Cycles simulation ends
     * Cycles		    CyclesStruct
     * project		    char*	Name of project
     */
    int             rotationYear;
    //int             nextSeedingDate;
    //int             nextSeedingYear;
    int             y;
    int             doy;
    int             i;
    int             c;
    time_t          begin_t, end_t;

    CyclesStruct    Cycles;     /* Model structure */
    char           *project;    /* Name of simulation */

    time (&begin_t);

    Cycles = (CyclesStruct)malloc (sizeof (*Cycles));

    printf ("\n\n");
    printf ("\t\t ######  ##    ##  ######  ##       ########  ######\n");
    printf ("\t\t##    ##  ##  ##  ##    ## ##       ##       ##    ##\n");
    printf ("\t\t##         ####   ##       ##       ##       ##\n");
    printf ("\t\t##          ##    ##       ##       ######    ######\n");
    printf ("\t\t##          ##    ##       ##       ##             ##\n");
    printf ("\t\t##    ##    ##    ##    ## ##       ##       ##    ##\n");
    printf ("\t\t ######     ##     ######  ######## ########  ######\n\n\n");

    verbose_mode = 0;

    while ((c = getopt (argc, argv, "vd")) != -1)
    {
        if (optind >= argc)
        {
            printf ("\nUsage: ./Cycles [-v] [-d] <project name>\n");
            printf ("\t-v Verbose mode\n");
            printf ("\t-d Debug mode\n");
            exit (1);
        }
        switch (c)
        {
            case 'v':
                verbose_mode = 1;
                printf ("Verbose mode turned on.\n");
                break;
            case 'd':
                debug_mode = 1;
                printf ("Debug mode turned on.\n");
                break;
            case '?':
                printf ("Option not recognisable %s\n", argv[optind]);
                break;
            default:
                break;
        }
    }

    if (optind >= argc)
    {
        printf ("ERROR: Please specify the name of project!\n");
        printf ("\nUsage: ./Cycles [-v] [-d] <project name>\n");
        printf ("\t-v Verbose mode\n");
        printf ("\t-d Ddebug mode\n");
        exit (1);
    }
    else
    {
        project = (char *)malloc ((strlen (argv[optind]) + 1) * sizeof (char));
        strcpy (project, argv[optind]);
    }

    printf ("Now running the %s simulation.\n\n", project);

    /* Read simulation control input file */
    ReadSimControl (project, &Cycles->SimControl);
    if (debug_mode)
        PrintSimContrl (Cycles->SimControl);

    /* Read soil description file */
    ReadSoil (Cycles->SimControl.soil_filename, &Cycles->Soil);
    if (debug_mode)
        PrintSoil (Cycles->Soil);

    /* Read crop description file */
    ReadCrop (Cycles->SimControl.crop_filename, &Cycles->Community);
    if (debug_mode)
        PrintCrop (Cycles->Community);

    /* Read field operation file */
    ReadOperation (Cycles->SimControl.operation_filename, &Cycles->CropManagement, &Cycles->Community, Cycles->SimControl.yearsInRotation);
    if (debug_mode)
        PrintOperation (Cycles->CropManagement.plantingOrder, Cycles->CropManagement.totalCropsPerRotation, Cycles->CropManagement.ForcedHarvest, Cycles->CropManagement.numHarvest, Cycles->CropManagement.Tillage, Cycles->CropManagement.numTillage, Cycles->CropManagement.FixedIrrigation, Cycles->CropManagement.numIrrigation, Cycles->CropManagement.FixedFertilization, Cycles->CropManagement.numFertilization);

    /* Read meteorological driver */
    ReadWeather (Cycles->SimControl.weather_filename, &Cycles->Weather, Cycles->SimControl.simStartYear, Cycles->SimControl.totalYears);
    if (debug_mode)
        PrintWeather (Cycles->Weather);

    InitializeOutput (project, &Cycles->Community, Cycles->Soil.totalLayers);

    /* Initialize model variables and parameters */
    Initialize (&Cycles->SimControl, &Cycles->Weather, &Cycles->Soil, &Cycles->Residue, &Cycles->SoilCarbon, &Cycles->Community, &Cycles->CropManagement, &Cycles->Snow);

    /* Compute crop thermal time */
    if (verbose_mode)
	printf ("Compute crop thermal time.\n");
    ComputeThermalTime (Cycles->SimControl.totalYears, &Cycles->Community, &Cycles->Weather);

    rotationYear = 0;

    printf ("\nSimulation running ...\n");

    Cycles->Summary = (SummaryStruct) {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (y = 0; y < Cycles->SimControl.totalYears; y++)
    {
        printf ("Year %4d (%4d)\n", y + 1, Cycles->SimControl.simStartYear + y);

        if (rotationYear < Cycles->SimControl.yearsInRotation)
            rotationYear++;
        else
            rotationYear = 1;
        if (debug_mode)
            printf ("*%-15s = %-d\n", "Rotation year", rotationYear);

	/* Initialize annual variables */
        for (i = 0; i < Cycles->Soil.totalLayers; i++)
        {
            Cycles->SoilCarbon.carbonMassInitial[i] = Cycles->Soil.SOC_Mass[i];
            Cycles->SoilCarbon.carbonMassFinal[i] = 0.0;
            Cycles->SoilCarbon.annualHumifiedCarbonMass[i] = 0.0;
            Cycles->SoilCarbon.annualRespiredCarbonMass[i] = 0.0;
            Cycles->SoilCarbon.annualRespiredResidueCarbonMass[i] = 0.0;
            Cycles->SoilCarbon.annualSoilCarbonDecompositionRate[i] = 0.0;
            Cycles->SoilCarbon.abgdBiomassInput[i] = 0.0;
            Cycles->SoilCarbon.rootBiomassInput[i] = 0.0;
            Cycles->SoilCarbon.rhizBiomassInput[i] = 0.0;
            Cycles->SoilCarbon.abgdCarbonInput[i] = 0.0;
            Cycles->SoilCarbon.rootCarbonInput[i] = 0.0;
            Cycles->SoilCarbon.annualNmineralization[i] = 0.0;
            Cycles->SoilCarbon.annualNImmobilization[i] = 0.0;
            Cycles->SoilCarbon.annualNNetMineralization[i] = 0.0;
            Cycles->SoilCarbon.annualAmmoniumNitrification = 0.0;
            Cycles->SoilCarbon.annualNitrousOxidefromNitrification = 0.0;
            Cycles->SoilCarbon.annualAmmoniaVolatilization = 0.0;
            Cycles->SoilCarbon.annualNO3Denitrification = 0.0;
            Cycles->SoilCarbon.annualNitrousOxidefromDenitrification = 0.0;
            Cycles->SoilCarbon.annualNitrateLeaching = 0.0;
            Cycles->SoilCarbon.annualAmmoniumLeaching = 0.0;

            Cycles->Residue.yearResidueBiomass = 0.0;
            Cycles->Residue.yearRootBiomass = 0.0;
            Cycles->Residue.yearRhizodepositionBiomass = 0.0;
        }

        /* Daily operations */
        for (doy = 1; doy < Cycles->Weather.lastDoy[y] + 1; doy++)
        {
            if (debug_mode)
                printf ("DOY %3.3d\n", doy);
            DailyOperations (rotationYear, y, doy, &Cycles->CropManagement, &Cycles->Community, &Cycles->Residue, &Cycles->SimControl, &Cycles->Snow, &Cycles->Soil, &Cycles->SoilCarbon, &Cycles->Weather, project);
            PrintDailyOutput (y, doy, Cycles->SimControl.simStartYear, &Cycles->Weather, &Cycles->Community, &Cycles->Soil, &Cycles->Snow, &Cycles->Residue, project);
        }

        for (i = 0; i < Cycles->Soil.totalLayers; i++)
            Cycles->SoilCarbon.carbonMassFinal[i] = Cycles->Soil.SOC_Mass[i];

        PrintAnnualOutput (y, Cycles->SimControl.simStartYear, &Cycles->Soil, &Cycles->SoilCarbon, project);
        PrintCarbonEvolution (y, Cycles->SimControl.simStartYear, Cycles->Soil.totalLayers, &Cycles->Soil, &Cycles->SoilCarbon, &Cycles->Residue, project);
        StoreSummary (&Cycles->Summary, &Cycles->SoilCarbon, &Cycles->Residue, Cycles->Soil.totalLayers, y);
    }

    PrintSummary (&Cycles->Summary, Cycles->SimControl.totalYears, project);

    FreeCyclesStruct (Cycles, Cycles->SimControl.totalYears);
    free (project);
    free (Cycles);

    time (&end_t);

    printf ("\nSimulation time: %-d seconds.\n", (int)(end_t - begin_t));

    return (0);
}
