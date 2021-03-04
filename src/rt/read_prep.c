#include "pihm.h"

void ReadPrep(const char fn[], const chemtbl_struct chemtbl[], const rttbl_struct *rttbl, forc_struct *forc)
{
    FILE           *fp;
    int             lno = 0;
    int             i, j;
    int             match;
    int             bytes_now;
    int             bytes_consumed = 0;
    int             index[MAXSPS];
    int             tsind;
    int            *nsps;
    double          temp_conc[MAXSPS];
    char            chemn[MAXSTRING];
    char            cmdstr[MAXSTRING];
    char            tempstr[2][MAXSTRING];

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    forc->nprcpc = CountOccurr(fp, "PRCP_CONC_TS");

    FindLine(fp, "BOF", &lno, fn);

    if (forc->nprcpc > 0)
    {
        forc->prcpc = (tsdata_struct *)malloc(forc->nprcpc * sizeof(tsdata_struct));
        nsps = (int *)malloc(forc->nprcpc * sizeof(int));

        NextLine(fp, cmdstr, &lno);
        for (i = 0; i < forc->nprcpc; i++)
        {
            match = sscanf(cmdstr, "%s %d %s %d", tempstr[0], &tsind, tempstr[1], &nsps[i]);
            if (match != 4 || i != tsind - 1 || strcasecmp(tempstr[0], "PRCP_CONC_TS") != 0 ||
                strcasecmp(tempstr[1], "PCONC") != 0)
            {
                pihm_error(ERR_WRONG_FORMAT, fn, lno);
            }
            // Skip header lines
            NextLine(fp, cmdstr, &lno);
            forc->prcpc[i].length = CountLines(fp, cmdstr, 1, "PRCP_CONC_TS");
        }

        // Rewind and read
        FindLine(fp, "BOF", &lno, fn);
        for (i = 0; i < forc->nprcpc; i++)
        {
            // Skip header lines
            NextLine(fp, cmdstr, &lno);

            NextLine(fp, cmdstr, &lno);
            bytes_consumed = 0;
            // Check header line with TIME keyword
            if (sscanf(cmdstr + bytes_consumed, "%s%n", tempstr[0], &bytes_now) != 1)
            {
                pihm_error(ERR_WRONG_FORMAT, fn, lno);
            }
            bytes_consumed += bytes_now;
            if (strcasecmp(tempstr[0], "TIME") != 0)
            {
                pihm_printf(VL_ERROR, "Expect header line starting with \"TIME\".\n");
                pihm_error(ERR_WRONG_FORMAT, fn, lno);
            }
            for (j = 0; j < nsps[i]; j++)
            {
                if (sscanf(cmdstr + bytes_consumed, "%s%n", chemn, &bytes_now) != 1)
                {
                    pihm_printf(VL_ERROR, "Error reading precipitation conc. in %s near Line %d.\n", fn, lno);
                    pihm_exit(EXIT_FAILURE);
                }
                bytes_consumed += bytes_now;

                index[j] = FindChem(chemn, chemtbl, rttbl->num_stc);

                if (index[j] < 0 || index[j] >= rttbl->num_spc)
                {
                    pihm_printf(VL_ERROR, "Error: Precipitation species %s is not defined in the simulation.\n", chemn);
                    pihm_exit(EXIT_FAILURE);
                }
                else
                {
                    pihm_printf(VL_VERBOSE, "  Precipitation concentration of '%s' is a time series.\n", chemn);
                }
            }

            forc->prcpc[i].ftime = (int *)malloc((forc->prcpc[i].length) * sizeof(int));
            forc->prcpc[i].data = (double **)malloc((forc->prcpc[i].length) * sizeof(double *));
            forc->prcpc[i].value = (double *)malloc(rttbl->num_spc * sizeof(double));

            for (j = 0; j < forc->prcpc[i].length; j++)
            {
                int             k, kk;

                forc->prcpc[i].data[j] = (double *)malloc(rttbl->num_spc * sizeof(double));

                NextLine(fp, cmdstr, &lno);
                if (!ReadTs(cmdstr, nsps[i], &forc->prcpc[i].ftime[j], temp_conc))
                {
                    pihm_error(ERR_WRONG_FORMAT, fn, lno);
                }

                for (k = 0; k < rttbl->num_spc; k++)
                {
                    // Species not described in the forcing file will be filled with the concentrations in .chem file
                    forc->prcpc[i].data[j][k] = rttbl->prcp_conc[k];

                    for (kk = 0; kk < nsps[i]; kk++)
                    {
                        if (index[kk] == k)
                        {
                            if (strcmp(chemtbl[k].name, "pH") == 0)
                            {
                                // Convert pH to H+ concentration
                                forc->prcpc[i].data[j][k] = (temp_conc[kk] < 7.0) ?
                                    pow(10, -temp_conc[kk]) : -pow(10, -temp_conc[kk] - 14);
                            }
                            else
                            {
                                forc->prcpc[i].data[j][k] = temp_conc[kk];
                            }
                            break;
                        }
                    }
                }
            }
        }
    }

    fclose(fp);
}
