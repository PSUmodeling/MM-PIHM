#include "pihm.h"

void ReadLai(const char *filename, forc_struct *forc,
    const atttbl_struct *atttbl)
{
    char            cmdstr[MAXSTRING];
    int             read_lai = 0;
    FILE           *lai_file;
    int             i, j;
    int             index;
    int             lno = 0;

    for (i = 0; i < nelem; i++)
    {
        if (atttbl->lai[i] > 0)
        {
            read_lai = 1;
            break;
        }
    }

    forc->nlai = 0;

    if (read_lai)
    {
        lai_file = pihm_fopen(filename, "r");
        pihm_printf(VL_VERBOSE, " Reading %s\n", filename);

        /* Start reading lai_file */
        FindLine(lai_file, "BOF", &lno, filename);

        forc->nlai = CountOccurr(lai_file, "LAI_TS");

        FindLine(lai_file, "BOF", &lno, filename);
        if (forc->nlai > 0)
        {
            forc->lai =
                (tsdata_struct *)malloc(forc->nlai * sizeof(tsdata_struct));

            NextLine(lai_file, cmdstr, &lno);
            for (i = 0; i < forc->nlai; i++)
            {
                ReadKeyword(cmdstr, "LAI_TS", &index, 'i', filename, lno);

                if (i != index - 1)
                {
                    pihm_printf(VL_ERROR,
                        "Error reading the %dth LAI time series.\n", i + 1);
                    pihm_printf(VL_ERROR, "Error in %s near Line %d.\n",
                        filename, lno);
                    pihm_exit(EXIT_FAILURE);
                }
                /* Skip header lines */
                NextLine(lai_file, cmdstr, &lno);
                NextLine(lai_file, cmdstr, &lno);
                forc->lai[i].length = CountLine(lai_file, cmdstr, 1, "LAI_TS");
            }

            /* Rewind and read */
            FindLine(lai_file, "BOF", &lno, filename);
            for (i = 0; i < forc->nlai; i++)
            {
                /* Skip header lines */
                NextLine(lai_file, cmdstr, &lno);
                NextLine(lai_file, cmdstr, &lno);
                NextLine(lai_file, cmdstr, &lno);

                forc->lai[i].ftime =
                    (int *)malloc(forc->lai[i].length * sizeof(int));
                forc->lai[i].data =
                    (double **)malloc(forc->lai[i].length * sizeof(double *));
                for (j = 0; j < forc->lai[i].length; j++)
                {
                    forc->lai[i].data[j] = (double *)malloc(sizeof(double));
                    NextLine(lai_file, cmdstr, &lno);
                    if (!ReadTs(cmdstr, 1, &forc->lai[i].ftime[j],
                        &forc->lai[i].data[j][0]))
                    {
                        pihm_printf(VL_ERROR, "Error reading LAI forcing.");
                        pihm_printf(VL_ERROR, "Error in %s near Line %d.\n",
                            filename, lno);
                        pihm_exit(EXIT_FAILURE);
                    }
                }
            }
        }

        fclose(lai_file);
    }
}
