#include "Cycles.h"

void ReadWeather (char *filename, WeatherStruct *Weather, int start_year, int total_years)
{
    /*
     * Read weather file
     * -----------------------------------------------------------------------
     * LOCAL VARIABLES
     *
     * Variable             Type        Description
     * ==========           ==========  ====================
     * weather_file	    FILE*	File pointer of weather file
     * fullname		    char*	Full file name of the weather file
     * cmdstr		    char[MAXSTRING]
     *					Command string
     * optstr		    char[MAXSTRING]
     *					Option argument string
     * y		    int		year counter
     * doy		    int		Day of year counter
     * temp_year	    int
     * temp_doy		    int
     * start_year_str	    char[5]
     */
    FILE           *weather_file;
    char           *fullname;
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];
    int             y, doy;
    int             temp_year, temp_doy;
    char            start_year_str[5];

    /* Open simulation control file */
    fullname = (char *)malloc ((strlen (filename) + 7) * sizeof (char));
    sprintf (fullname, "input/%s", filename);
    weather_file = fopen (fullname, "r");

    if (weather_file == NULL)
    {
        printf ("\nError: Cannot find the weather file %s!\n", filename);
        exit (1);
    }
    else
        printf ("%-30s input/%s.\n", "Read weather file:", filename);

    free (fullname);

    /* Allocate */
    Weather->wind = (double **)malloc (total_years * sizeof (double *));
    Weather->ETref = (double **)malloc (total_years * sizeof (double *));
    Weather->precipitation = (double **)malloc (total_years * sizeof (double *));
    Weather->RHmax = (double **)malloc (total_years * sizeof (double *));
    Weather->RHmin = (double **)malloc (total_years * sizeof (double *));
    Weather->solarRadiation = (double **)malloc (total_years * sizeof (double *));
    Weather->tMax = (double **)malloc (total_years * sizeof (double *));
    Weather->tMin = (double **)malloc (total_years * sizeof (double *));
    Weather->yearlyAmplitude = (double *)malloc (total_years * sizeof (double));
    Weather->annualAverageTemperature = (double *)malloc (total_years * sizeof (double));
    Weather->lastDoy = (int *)malloc (total_years * sizeof (int));
    for (y = 0; y < total_years; y++)
    {
        Weather->wind[y] = (double *)malloc (366 * sizeof (double));
        Weather->ETref[y] = (double *)malloc (366 * sizeof (double));
        Weather->precipitation[y] = (double *)malloc (366 * sizeof (double));
        Weather->RHmax[y] = (double *)malloc (366 * sizeof (double));
        Weather->RHmin[y] = (double *)malloc (366 * sizeof (double));
        Weather->solarRadiation[y] = (double *)malloc (366 * sizeof (double));
        Weather->tMax[y] = (double *)malloc (366 * sizeof (double));
        Weather->tMin[y] = (double *)malloc (366 * sizeof (double));
    }

    sprintf (start_year_str, "%4.4d", start_year);

    /* Read weather file */
    FindLine (weather_file, "BOF");

    /* Read in site information and count number of weather records */
    NextLine (weather_file, cmdstr);
    ReadKeywordDouble (cmdstr, "LATITUDE", &Weather->siteLatitude);

    NextLine (weather_file, cmdstr);
    ReadKeywordDouble (cmdstr, "ALTITUDE", &Weather->siteAltitude);

    NextLine (weather_file, cmdstr);
    ReadKeywordDouble (cmdstr, "SCREENING_HEIGHT", &Weather->screeningHeight);

    /* Skip header */
    NextLine (weather_file, cmdstr);

    strcpy (cmdstr, "\0");

    NextLine (weather_file, cmdstr);
    while (strcasecmp (cmdstr, "EOF") != 0)
    {
        sscanf (cmdstr, "%s", optstr);
        if (strcasecmp (start_year_str, optstr) == 0)
        {
            for (y = 0; y < total_years; y++)
            {
                for (doy = 1; doy < 367; doy++)
                {
                    sscanf (cmdstr, "%d %d %*f %*f %*f %*f %*f %*f %*f", &temp_year, &temp_doy);
                    if (temp_year == y + start_year && temp_doy == doy)
                    {
                        sscanf (cmdstr, "%*d %*d %lf %lf %lf %lf %lf %lf %lf", &Weather->precipitation[y][doy - 1], &Weather->tMax[y][doy - 1], &Weather->tMin[y][doy - 1], &Weather->solarRadiation[y][doy - 1], &Weather->RHmax[y][doy - 1], &Weather->RHmin[y][doy - 1], &Weather->wind[y][doy - 1]);
                        if (doy == 366)
                            Weather->lastDoy[y] = 366;
                        NextLine (weather_file, cmdstr);
                    }
                    else if (doy == 366 && temp_year == y + start_year + 1 && temp_doy == 1)
                        Weather->lastDoy[y] = 365;
                    else if (doy == 366 && feof (weather_file))
                        Weather->lastDoy[y] = 365;
                    else
                    {
                        printf ("ERROR: Please check your weather input file near YEAR: %4.4d, DOY: %-d, expecting %4.4d-%-d, eof status %d\n", temp_year, temp_doy, y, doy, feof (weather_file));
                        exit (1);
                    }
                }
            }
            break;
        }

        NextLine (weather_file, cmdstr);
    }

    fclose (weather_file);
}
