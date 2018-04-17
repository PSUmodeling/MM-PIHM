#include "pihm.h"

void ReadAtt(const char *filename, atttbl_struct *atttbl)
{
    int             i;
    FILE           *att_file;   /* Pointer to .att file */
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    att_file = fopen(filename, "r");
    CheckFile(att_file, filename);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filename);

    atttbl->soil = (int *)malloc(nelem * sizeof(int));
    atttbl->geol = (int *)malloc(nelem * sizeof(int));
    atttbl->lc = (int *)malloc(nelem * sizeof(int));
    atttbl->bc = (int **)malloc(nelem * sizeof(int *));
    for (i = 0; i < nelem; i++)
    {
        atttbl->bc[i] = (int *)malloc(NUM_EDGE * sizeof(int));
    }
    atttbl->meteo = (int *)malloc(nelem * sizeof(int));
    atttbl->lai = (int *)malloc(nelem * sizeof(int));
    atttbl->source = (int *)malloc(nelem * sizeof(int));

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
            PIHMprintf(VL_ERROR,
                "Error reading attribute of the %dth element.\n", i + 1);
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            PIHMexit(EXIT_FAILURE);
        }
    }

    fclose(att_file);
}
