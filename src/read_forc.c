#include "pihm.h"

void ReadMeteo(const char fn[], forc_struct *forc)
{
    FILE           *fp;
    char            cmdstr[MAXSTRING];
    char            tempstr[2][MAXSTRING];
    int             i, j;
    int             match;
    int             index;
    int             lno = 0;

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    FindLine(fp, "BOF", &lno, fn);

    forc->nmeteo = CountOccurr(fp, "METEO_TS");

    FindLine(fp, "BOF", &lno, fn);
    if (forc->nmeteo > 0)
    {
        forc->meteo =
            (tsdata_struct *)malloc(forc->nmeteo * sizeof(tsdata_struct));

        NextLine(fp, cmdstr, &lno);
        for (i = 0; i < forc->nmeteo; i++)
        {
            match = sscanf(cmdstr, "%s %d %s %lf", tempstr[0], &index,
                tempstr[1], &forc->meteo[i].zlvl_wind);
            if (match != 4 || i != index - 1 ||
                strcasecmp(tempstr[0], "METEO_TS") != 0 ||
                strcasecmp(tempstr[1], "WIND_LVL") != 0)
            {
                pihm_printf(VL_ERROR,
                    "Error reading the %dth meteorological forcing"
                    " time series.\n", i + 1);
                pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", fn, lno);
                pihm_exit(EXIT_FAILURE);
            }
            /* Check header lines */
            NextLine(fp, cmdstr, &lno);
            if (!CheckHeader(cmdstr, 8, "TIME", "PRCP", "SFCTMP", "RH",
                "SFCSPD", "SOLAR", "LONGWV", "PRES"))
            {
                pihm_printf(VL_ERROR, "Meteorological forcing file header error.\n");
                pihm_exit(EXIT_FAILURE);
            }
            forc->meteo[i].length = CountLine(fp, cmdstr, 1, "METEO_TS");
        }

        /* Rewind and read */
        FindLine(fp, "BOF", &lno, fn);
        for (i = 0; i < forc->nmeteo; i++)
        {
            /* Skip header lines */
            NextLine(fp, cmdstr, &lno);
            NextLine(fp, cmdstr, &lno);

            forc->meteo[i].ftime =
                (int *)malloc(forc->meteo[i].length * sizeof(int));
            forc->meteo[i].data =
                (double **)malloc(forc->meteo[i].length * sizeof(double *));
            for (j = 0; j < forc->meteo[i].length; j++)
            {
                forc->meteo[i].data[j] =
                    (double *)malloc(NUM_METEO_VAR * sizeof(double));
                NextLine(fp, cmdstr, &lno);
                if (!ReadTs(cmdstr, NUM_METEO_VAR, &forc->meteo[i].ftime[j],
                    &forc->meteo[i].data[j][0]))
                {
                    pihm_printf(VL_ERROR,
                        "Error reading meteorological forcing.");
                    pihm_printf(VL_ERROR, "Error in %s near Line %d.\n",
                        fn, lno);
                    pihm_exit(EXIT_FAILURE);
                }
            }
        }
    }

    fclose(fp);
}
