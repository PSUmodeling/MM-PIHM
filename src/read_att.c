#include "pihm.h"

void ReadAtt(const char fn[], atttbl_struct *atttbl)
{
    int             i;
    FILE           *fp;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    atttbl->soil   = (int *)malloc(nelem * sizeof(int));
    atttbl->geol   = (int *)malloc(nelem * sizeof(int));
    atttbl->lc     = (int *)malloc(nelem * sizeof(int));
    atttbl->bc     = (int **)malloc(nelem * sizeof(int *));
    atttbl->meteo  = (int *)malloc(nelem * sizeof(int));
    atttbl->lai    = (int *)malloc(nelem * sizeof(int));
    for (i = 0; i < nelem; i++)
    {
        atttbl->bc[i] = (int *)malloc(NUM_EDGE * sizeof(int));
    }

    /* Check header line */
    NextLine(fp, cmdstr, &lno);
    if (!CheckHeader(cmdstr, 9, "INDEX", "SOIL", "GEOL", "LC", "METEO", "LAI",
        "BC1", "BC2", "BC3"))
    {
        pihm_printf(VL_ERROR, "Attribute file header error.\n");
        pihm_exit(EXIT_FAILURE);
    }
    for (i = 0; i < nelem; i++)
    {
        NextLine(fp, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %d %d %d %d %d %d %d %d",
            &index,
            &atttbl->soil[i], &atttbl->geol[i], &atttbl->lc[i],
            &atttbl->meteo[i], &atttbl->lai[i],
            &atttbl->bc[i][0], &atttbl->bc[i][1], &atttbl->bc[i][2]);
        if (match != 9 || index != i + 1)
        {
            pihm_printf(VL_ERROR,
                "Error reading attribute of the %dth element.\n", i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", fn, lno);
            pihm_exit(EXIT_FAILURE);
        }
    }

    fclose(fp);
}
