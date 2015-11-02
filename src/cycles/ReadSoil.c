#include "Cycles.h"

void ReadSoil (char *filename, SoilStruct *Soil)
{
    /*
     * Read soil description file
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * soil_file	    FILE*	File pointer of soil description file
     * fullname		    char*	Full file name of the soil description
     *					  file
     * cmdstr		    char[MAXSTRING]
     *					Command string
     * i		    int		Loop counter
     */
    FILE           *soil_file;
    char           *fullname;
    char            cmdstr[MAXSTRING];
    int             i;

    /* Open soil file */
    fullname = (char *)malloc ((strlen (filename) + 7) * sizeof (char));
    sprintf (fullname, "input/%s", filename);
    soil_file = fopen (fullname, "r");

    if (soil_file == NULL)
    {
        printf ("\nERROR: Cannot find the soil file %s!\n", filename);
        exit (1);
    }
    else
        printf ("%-30s input/%s.\n", "Read soil initialization file:", filename);

    free (fullname);

    /* Read soil file */
    FindLine (soil_file, "BOF");

    NextLine (soil_file, cmdstr);
    ReadKeywordDouble (cmdstr, "CURVE_NUMBER", &Soil->Curve_Number);

    NextLine (soil_file, cmdstr);
    ReadKeywordDouble (cmdstr, "SLOPE", &Soil->Percent_Slope);
    Soil->Percent_Slope /= 100.0;

    NextLine (soil_file, cmdstr);
    ReadKeywordInt (cmdstr, "TOTAL_LAYERS", &Soil->totalLayers);

    /* Allocate memories for soil class */
    Soil->layerThickness = (double *)malloc (Soil->totalLayers * sizeof (double));
    Soil->Clay = (double *)malloc (Soil->totalLayers * sizeof (double));
    Soil->Sand = (double *)malloc (Soil->totalLayers * sizeof (double));
    Soil->IOM = (double *)malloc (Soil->totalLayers * sizeof (double));
    Soil->BD = (double *)malloc (Soil->totalLayers * sizeof (double));
    Soil->FC = (double *)malloc (Soil->totalLayers * sizeof (double));
    Soil->PWP = (double *)malloc (Soil->totalLayers * sizeof (double));
    Soil->NO3 = (double *)malloc (Soil->totalLayers * sizeof (double));
    Soil->NH4 = (double *)malloc (Soil->totalLayers * sizeof (double));

    /* Skip header */
    NextLine (soil_file, cmdstr);

    for (i = 0; i < Soil->totalLayers; i++)
    {
        NextLine (soil_file, cmdstr);
        sscanf (cmdstr, "%*d %lg %lf %lf %lf %lf %lf %lf %lf %lf", &Soil->layerThickness[i], &Soil->Clay[i], &Soil->Sand[i], &Soil->IOM[i], &Soil->BD[i], &Soil->FC[i], &Soil->PWP[i], &Soil->NO3[i], &Soil->NH4[i]);
    }

    fclose (soil_file);
}
