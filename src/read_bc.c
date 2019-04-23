#include "pihm.h"

#if defined(_RT_)
void ReadBc(const char *filename, forc_struct *forc,
    const atttbl_struct *atttbl, const rttbl_struct *rttbl,
    const chemtbl_struct chemtbl[])
#else
void ReadBc(const char *filename, forc_struct *forc,
    const atttbl_struct *atttbl)
#endif
{
    int             i, j;
    FILE           *bc_file;    /* Pointer to .ibc file */
    int             read_bc = 0;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;
#if defined(_RT_)
    int             bytes_now;
    int             bytes_consumed = 0;
    int             ind[MAXSPS];
    char            chemn[MAXSTRING];
    double          bcval[MAXSPS + 1];
#endif

#if defined(_RT_)
    for (j = 0; j < MAXSPS; j++)
    {
        ind[j] = BADVAL;
    }
#endif

    for (i = 0; i < nelem; i++)
    {
        for (j = 0; j < NUM_EDGE; j++)
        {
            if (atttbl->bc[i][j] != 0)
            {
                read_bc = 1;
                break;
            }
#if defined(_FBR_)
             if (atttbl->fbr_bc[i][j] != 0)
             {
                read_bc = 1;
                break;
            }
#endif
        }
    }

    forc->nbc = 0;

    if (read_bc)
    {
        bc_file = fopen(filename, "r");
        CheckFile(bc_file, filename);
        PIHMprintf(VL_VERBOSE, " Reading %s\n", filename);

        FindLine(bc_file, "BOF", &lno, filename);

        forc->nbc = CountOccurr(bc_file, "BC_TS");

        FindLine(bc_file, "BOF", &lno, filename);
        if (forc->nbc > 0)
        {
            forc->bc =
                (tsdata_struct *)malloc(forc->nbc * sizeof(tsdata_struct));

            NextLine(bc_file, cmdstr, &lno);
            for (i = 0; i < forc->nbc; i++)
            {
                match = sscanf(cmdstr, "%*s %d", &index);
                if (match != 1 || i != index - 1)
                {
                    PIHMprintf(VL_ERROR,
                        "Error reading the %dth boundary condition "
                        "time series.\n", i + 1);
                    PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n",
                        filename, lno);
                    PIHMexit(EXIT_FAILURE);
                }
#if defined(_RT_)
                int             k;

                /* When reactive transport is turned on, the header line
                 * contains the names of species that need to be read */
                NextLine(bc_file, cmdstr, &lno);

                /* Skip the first two columns (TIME and HEAD/FLUX) */
                sscanf(cmdstr + bytes_consumed, "%*s %*s%n", &bytes_now);
                bytes_consumed += bytes_now;

                for (k = 0; k < rttbl->NumStc; k++)
                {
                    if (sscanf(cmdstr + bytes_consumed, "%s%n", chemn,
                        &bytes_now) == 1)
                    {
                        bytes_consumed += bytes_now;

                        ind[k] = FindChem(chemn, chemtbl, rttbl->NumStc);

                        if (ind[k] < 0)
                        {
                            PIHMprintf(VL_ERROR, "Error finding chemical %s.\n",
                                chemn);
                            PIHMexit(EXIT_FAILURE);
                        }
                    }
                    else
                    {
                        PIHMprintf(VL_ERROR,
                            "Error reading primary species concentrations.\n");
                        PIHMexit(EXIT_FAILURE);
                    }
                }

                /* Skip unit header line */
                NextLine(bc_file, cmdstr, &lno);
#else
                /* Skip header lines */
                NextLine(bc_file, cmdstr, &lno);
                NextLine(bc_file, cmdstr, &lno);
#endif
                forc->bc[i].length = CountLine(bc_file, cmdstr, 1, "BC_TS");
            }

            /* Rewind and read */
            FindLine(bc_file, "BOF", &lno, filename);
            for (i = 0; i < forc->nbc; i++)
            {
                /* Skip header lines */
                NextLine(bc_file, cmdstr, &lno);
                NextLine(bc_file, cmdstr, &lno);
                NextLine(bc_file, cmdstr, &lno);

                forc->bc[i].ftime =
                    (int *)malloc(forc->bc[i].length * sizeof(int));
                forc->bc[i].data =
                    (double **)malloc(forc->bc[i].length * sizeof(double *));
                for (j = 0; j < forc->bc[i].length; j++)
                {
#if defined(_RT_)
                    forc->bc[i].data[j] =
                        (double *)malloc((rttbl->NumStc + 1) * sizeof(double));
#else
                    forc->bc[i].data[j] = (double *)malloc(sizeof(double));
#endif
                    NextLine(bc_file, cmdstr, &lno);
#if defined(_RT_)
                    if (!ReadTS(cmdstr, &forc->bc[i].ftime[j],
                        bcval, rttbl->NumStc + 1))
#else
                    if (!ReadTS(cmdstr, &forc->bc[i].ftime[j],
                        &forc->bc[i].data[j][0], 1))
#endif
                    {
                        PIHMprintf(VL_ERROR,
                            "Error reading boundary condition.");
                        PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n",
                            filename, lno);
                        PIHMexit(EXIT_FAILURE);
                    }
#if defined(_RT_)
                    else
                    {
                        int             k;

                        forc->bc[i].data[j][0] = bcval[0];

                        for (k = 0; k < rttbl->NumStc; k++)
                        {
                            if (strcmp(chemtbl[k].ChemName, "pH") == 0)
                            {
                                /* Convert pH to H+ concentration */
                                forc->bc[i].data[j][k + 1] =
                                    (bcval[1 + ind[k]] < 7.0) ?
                                    pow(10, -bcval[1 + ind[k]]) :
                                    -pow(10, -bcval[1 + ind[k]] - 14);
                            }
                            else
                            {
                                forc->bc[i].data[j][k + 1] = bcval[1 + ind[k]];
                            }
                        }
                    }
#endif
                }
            }
        }

        fclose(bc_file);
    }
}
