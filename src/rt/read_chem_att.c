#include "pihm.h"

void ReadChemAtt(const char *filename, atttbl_struct *atttbl)
{
    int             i;
    FILE           *fp;   /* Pointer to .att file */
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;
#if defined(_FBR_)
    const int       NVOL = 4;
#else
    const int       NVOL = 2;
#endif

    fp = fopen(filename, "r");
    CheckFile(fp, filename);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filename);

    atttbl->prcpc = (int *)malloc(nelem * sizeof(int));
    atttbl->chem_ic = (int **)malloc(nelem * sizeof(int *));

    /* Skip header line */
    NextLine(fp, cmdstr, &lno);
    for (i = 0; i < nelem; i++)
    {
        atttbl->chem_ic[i] = (int *)malloc(4 * sizeof(int));

        NextLine(fp, cmdstr, &lno);
#if defined(_FBR_)
        match = sscanf(cmdstr, "%d %d %d %d %d %d %d %d %d %d",
            &index, &atttbl->prcpc[i],
            &atttbl->chem_ic[i][0], &atttbl->chem_ic[i][1],
            &atttbl->chem_ic[i][2], &atttbl->chem_ic[i][3]);
#else
        match = sscanf(cmdstr, "%d %d %d %d",
            &index, &atttbl->prcpc[i],
            &atttbl->chem_ic[i][0], &atttbl->chem_ic[i][1]);
#endif

        if (match != NVOL + 2)
        {
            PIHMprintf(VL_ERROR,
                "Error reading chemistry attribute of the %dth element.\n",
                i + 1);
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            PIHMexit(EXIT_FAILURE);
        }
    }

    fclose(fp);
}
