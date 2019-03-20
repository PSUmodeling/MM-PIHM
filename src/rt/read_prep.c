#include "pihm.h"

void ReadPrep(const char filen[], const chemtbl_struct chemtbl[],
    Chem_Data rt)
{
    FILE           *fp;
    int             lno = 0;
    int             i, k;
    int             match;
    int             bytes_now;
    int             bytes_consumed = 0;
    char            cmdstr[MAXSTRING];

    fp = fopen(filen, "r");
    CheckFile(fp, filen);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filen);

    rt->TSD_prepconc = (tsdata_struct *)malloc(sizeof(tsdata_struct));

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "PCONC", &rt->TSD_prepconc[0].nspec, 'i', filen, lno);

    rt->prepconcindex = (int *)malloc(rt->TSD_prepconc[0].nspec * sizeof(int));
    /* The number of primary species must be equal to the number of primary
     * species specified before. */

    NextLine(fp, cmdstr, &lno);
    for (i = 0; i < rt->TSD_prepconc[0].nspec; i++)
    {
        match = sscanf(cmdstr + bytes_consumed, "%d%n", &rt->prepconcindex[i],
            &bytes_now);
        if (match != 1)
        {
            PIHMprintf(VL_ERROR,
                "Error reading precipitation conc. in %s near Line %d.\n",
                filen, lno);
            PIHMexit(EXIT_FAILURE);
        }
        bytes_consumed += bytes_now;

        if (rt->prepconcindex[i] > 0)
        {
            assert(rt->prepconcindex[i] <= NumSpc);
            PIHMprintf(VL_NORMAL,
                "  Precipitation conc of '%s' is a time series. \n",
                chemtbl[rt->prepconcindex[i] - 1].ChemName);
        }
    }

    rt->TSD_prepconc[0].length = CountLine(fp, cmdstr, 1, "PCONC");

    rt->TSD_prepconc[0].ftime =
        (int *)malloc((rt->TSD_prepconc[0].length) * sizeof(int));
    rt->TSD_prepconc[0].data =
        (double **)malloc((rt->TSD_prepconc[0].length) * sizeof(double *));
    rt->TSD_prepconc[0].value =
        (double *)malloc(rt->TSD_prepconc[0].nspec * sizeof(double));

    FindLine(fp, "BOF", &lno, filen);
    NextLine(fp, cmdstr, &lno);
    NextLine(fp, cmdstr, &lno);

    for (i = 0; i < rt->TSD_prepconc[0].length; i++)
    {
        rt->TSD_prepconc[0].data[i] =
            (double *)malloc(rt->TSD_prepconc[0].nspec * sizeof(double));

        NextLine(fp, cmdstr, &lno);
        ReadTS(cmdstr, &rt->TSD_prepconc[0].ftime[i],
            &rt->TSD_prepconc[0].data[i][0], rt->TSD_prepconc[0].nspec);
    }

    /* Convert pH to H+ concentration */
    for (i = 0; i < rt->TSD_prepconc[0].nspec; i++)
    {
        if (rt->prepconcindex[i] > 0 &&
            !strcmp(chemtbl[rt->prepconcindex[i] - 1].ChemName,
            "pH"))
        {
            for (k = 0; k < rt->TSD_prepconc[0].length; k++)
            {
                rt->TSD_prepconc[0].data[k][i] =
                    (rt->TSD_prepconc[0].data[k][i] < 7.0) ?
                    pow(10, -rt->TSD_prepconc[0].data[k][i]) :
                    -pow(10, -rt->TSD_prepconc[0].data[k][i] - 14);
            }
            break;
        }
    }

    fclose(fp);
}
