#include "pihm.h"

void ReadAtt(const char filename[], atttbl_struct *atttbl)
{
    int             i;
    FILE           *att_file;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    att_file = pihm_fopen(filename, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filename);

    atttbl->soil   = (int *)malloc(nelem * sizeof(int));
    atttbl->geol   = (int *)malloc(nelem * sizeof(int));
    atttbl->lc     = (int *)malloc(nelem * sizeof(int));
    atttbl->bc     = (int **)malloc(nelem * sizeof(int *));
    atttbl->meteo  = (int *)malloc(nelem * sizeof(int));
    atttbl->lai    = (int *)malloc(nelem * sizeof(int));
    atttbl->source = (int *)malloc(nelem * sizeof(int));
    for (i = 0; i < nelem; i++)
    {
        atttbl->bc[i] = (int *)malloc(NUM_EDGE * sizeof(int));
    }

    NextLine(att_file, cmdstr, &lno);
    for (i = 0; i < nelem; i++)
    {
        NextLine(att_file, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %d %d %d %d %d %d %d %d %d",
            &index,
            &atttbl->soil[i], &atttbl->geol[i], &atttbl->lc[i],
            &atttbl->meteo[i], &atttbl->lai[i], &atttbl->source[i],
            &atttbl->bc[i][0], &atttbl->bc[i][1], &atttbl->bc[i][2]);
        if (match != 10)
        {
            pihm_printf(VL_ERROR,
                "Error reading attribute of the %dth element.\n", i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            pihm_exit(EXIT_FAILURE);
        }
    }

    fclose(att_file);
}
