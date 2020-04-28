#include "pihm.h"

void ReadForc(const char *filename, forc_struct *forc)
{
    FILE           *meteo_file;
    char            cmdstr[MAXSTRING];
    int             i, j;
    int             match;
    int             index;
    int             lno = 0;

    meteo_file = fopen(filename, "r");
    CheckFile(meteo_file, filename);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filename);

    FindLine(meteo_file, "BOF", &lno, filename);

    forc->nmeteo = CountOccurr(meteo_file, "METEO_TS");

    FindLine(meteo_file, "BOF", &lno, filename);
    if (forc->nmeteo > 0)
    {
        forc->meteo =
            (tsdata_struct *)malloc(forc->nmeteo * sizeof(tsdata_struct));

        NextLine(meteo_file, cmdstr, &lno);
        for (i = 0; i < forc->nmeteo; i++)
        {
            match = sscanf(cmdstr, "%*s %d %*s %lf",
                &index, &forc->meteo[i].zlvl_wind);
            if (match != 2 || i != index - 1)
            {
                PIHMprintf(VL_ERROR,
                    "Error reading the %dth meteorological forcing"
                    " time series.\n", i + 1);
                PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n",
                    filename, lno);
                PIHMexit(EXIT_FAILURE);
            }
            /* Skip header lines */
            NextLine(meteo_file, cmdstr, &lno);
            NextLine(meteo_file, cmdstr, &lno);
            forc->meteo[i].length =
                CountLine(meteo_file, cmdstr, 1, "METEO_TS");
        }

        /* Rewind and read */
        FindLine(meteo_file, "BOF", &lno, filename);
        for (i = 0; i < forc->nmeteo; i++)
        {
            /* Skip header lines */
            NextLine(meteo_file, cmdstr, &lno);
            NextLine(meteo_file, cmdstr, &lno);
            NextLine(meteo_file, cmdstr, &lno);

            forc->meteo[i].ftime =
                (int *)malloc(forc->meteo[i].length * sizeof(int));
            forc->meteo[i].data =
                (double **)malloc(forc->meteo[i].length * sizeof(double *));
            for (j = 0; j < forc->meteo[i].length; j++)
            {
                forc->meteo[i].data[j] =
                    (double *)malloc(NUM_METEO_VAR * sizeof(double));
                NextLine(meteo_file, cmdstr, &lno);
                if (!ReadTS(cmdstr, &forc->meteo[i].ftime[j],
                    &forc->meteo[i].data[j][0], NUM_METEO_VAR))
                {
                    PIHMprintf(VL_ERROR,
                        "Error reading meteorological forcing.");
                    PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n",
                        filename, lno);
                    PIHMexit(EXIT_FAILURE);
                }
            }
        }
    }

    fclose(meteo_file);
}
