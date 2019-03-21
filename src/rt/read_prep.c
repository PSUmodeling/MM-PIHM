#include "pihm.h"

void ReadPrep(const char filen[], const chemtbl_struct chemtbl[],
    const double prcp_conc[], forc_struct *forc)
{
    FILE           *fp;
    int             lno = 0;
    int             i;
    int             bytes_now;
    int             bytes_consumed = 0;
    int             index[MAXSPS];
    int             nsps;
    double          temp_conc[MAXSPS];
    char            cmdstr[MAXSTRING];

    fp = fopen(filen, "r");
    CheckFile(fp, filen);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filen);

    /* Currently, only one precipitation concentration time series is
     * supported */
    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "PCONC", &nsps, 'i', filen, lno);

    NextLine(fp, cmdstr, &lno);
    for (i = 0; i < nsps; i++)
    {
        if (sscanf(cmdstr + bytes_consumed, "%d%n", &index[i], &bytes_now) != 1)
        {
            PIHMprintf(VL_ERROR,
                "Error reading precipitation conc. in %s near Line %d.\n",
                filen, lno);
            PIHMexit(EXIT_FAILURE);
        }
        bytes_consumed += bytes_now;

        if (index[i] > 0)
        {
            if (index[i] <= NumSpc)
            {
                PIHMprintf(VL_VERBOSE,
                    "  Precipitation concentration of '%s' is a time series.\n",
                    chemtbl[index[i] - 1].ChemName);
            }
            else
            {
                PIHMprintf(VL_VERBOSE,
                    "Error: Precipitation species index is larger than "
                    "number of primary species.\n");
                PIHMexit(EXIT_FAILURE);
            }
        }
    }

    forc->TSD_prepconc.length = CountLine(fp, cmdstr, 1, "PCONC");

    forc->TSD_prepconc.ftime =
        (int *)malloc((forc->TSD_prepconc.length) * sizeof(int));
    forc->TSD_prepconc.data =
        (double **)malloc((forc->TSD_prepconc.length) * sizeof(double *));
    forc->TSD_prepconc.value =
        (double *)malloc(NumSpc * sizeof(double));

    FindLine(fp, "BOF", &lno, filen);
    NextLine(fp, cmdstr, &lno);
    NextLine(fp, cmdstr, &lno);

    for (i = 0; i < forc->TSD_prepconc.length; i++)
    {
        int             k, kk;

        forc->TSD_prepconc.data[i] =
            (double *)malloc(NumSpc * sizeof(double));

        NextLine(fp, cmdstr, &lno);
        ReadTS(cmdstr, &forc->TSD_prepconc.ftime[i], temp_conc, nsps);

        for (k = 0; k < NumSpc; k++)
        {
            /* Species not described in the forcing file will be filled with
             * the concentrations in .chem file */
            forc->TSD_prepconc.data[i][k] = prcp_conc[k];

            for (kk = 0; kk < nsps; kk++)
            {
                if (index[kk] - 1 == k)
                {
                    if (strcmp(chemtbl[k].ChemName, "pH") == 0)
                    {
                        /* Convert pH to H+ concentration */
                        forc->TSD_prepconc.data[i][k] = (temp_conc[kk] < 7.0) ?
                            pow(10, -temp_conc[kk]) :
                            -pow(10, -temp_conc[kk] - 14);
                    }
                    else
                    {
                        forc->TSD_prepconc.data[i][k] = temp_conc[kk];
                    }
                    break;
                }
            }
        }
    }

    fclose(fp);
}
