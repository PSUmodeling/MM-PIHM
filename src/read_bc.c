#include "pihm.h"

void ReadBc(const char *filename, forc_struct *forc,
    const atttbl_struct *atttbl)
{
    int             i, j;
    FILE           *bc_file;    /* Pointer to .ibc file */
    int             read_bc = 0;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

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
                /* Skip header lines */
                NextLine(bc_file, cmdstr, &lno);
                NextLine(bc_file, cmdstr, &lno);
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
                    forc->bc[i].data[j] = (double *)malloc(sizeof(double));
                    NextLine(bc_file, cmdstr, &lno);
                    if (!ReadTS(cmdstr, &forc->bc[i].ftime[j],
                        &forc->bc[i].data[j][0], 1))
                    {
                        PIHMprintf(VL_ERROR,
                            "Error reading boundary condition.");
                        PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n",
                            filename, lno);
                        PIHMexit(EXIT_FAILURE);
                    }
                }
            }
        }

        fclose(bc_file);
    }
}
