#include "pihm.h"

void ReadLai(const char fn[], const atttbl_struct *atttbl,
    forc_struct *forc)
{
    char            cmdstr[MAXSTRING];
    int             read_lai = 0;
    FILE           *fp;
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
        fp = pihm_fopen(fn, "r");
        pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

        /* Start reading lai_file */
        FindLine(fp, "BOF", &lno, fn);

        forc->nlai = CountOccurr(fp, "LAI_TS");

        FindLine(fp, "BOF", &lno, fn);
        if (forc->nlai > 0)
        {
            forc->lai =
                (tsdata_struct *)malloc(forc->nlai * sizeof(tsdata_struct));

            NextLine(fp, cmdstr, &lno);
            for (i = 0; i < forc->nlai; i++)
            {
                ReadKeyword(cmdstr, "LAI_TS", 'i', fn, lno, &index);

                if (i != index - 1)
                {
                    pihm_printf(VL_ERROR,
                        "Error reading the %dth LAI time series.\n", i + 1);
                    pihm_printf(VL_ERROR, "Error in %s near Line %d.\n",
                        fn, lno);
                    pihm_exit(EXIT_FAILURE);
                }
                /* Check header lines */
                NextLine(fp, cmdstr, &lno);
                if (!CheckHeader(cmdstr, 2, "TIME", "LAI"))
                {
                    pihm_printf(VL_ERROR, "LAI forcing file header error.\n");
                    pihm_exit(EXIT_FAILURE);
                }
                forc->lai[i].length = CountLine(fp, cmdstr, 1, "LAI_TS");
            }

            /* Rewind and read */
            FindLine(fp, "BOF", &lno, fn);
            for (i = 0; i < forc->nlai; i++)
            {
                /* Skip header lines */
                NextLine(fp, cmdstr, &lno);
                NextLine(fp, cmdstr, &lno);

                forc->lai[i].ftime =
                    (int *)malloc(forc->lai[i].length * sizeof(int));
                forc->lai[i].data =
                    (double **)malloc(forc->lai[i].length * sizeof(double *));
                for (j = 0; j < forc->lai[i].length; j++)
                {
                    forc->lai[i].data[j] = (double *)malloc(sizeof(double));
                    NextLine(fp, cmdstr, &lno);
                    if (!ReadTs(cmdstr, 1, &forc->lai[i].ftime[j],
                        &forc->lai[i].data[j][0]))
                    {
                        pihm_printf(VL_ERROR, "Error reading LAI forcing."
                            "Error in %s near Line %d.\n", fn, lno);
                        pihm_exit(EXIT_FAILURE);
                    }
                }
            }
        }

        fclose(fp);
    }
}
