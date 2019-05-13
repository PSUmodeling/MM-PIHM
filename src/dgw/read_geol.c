#include "pihm.h"

void ReadGeol(const char *filename, geoltbl_struct *geoltbl)
{
    FILE           *geol_file;
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    geol_file = fopen (filename, "r");
    CheckFile (geol_file, filename);
    PIHMprintf (VL_VERBOSE, " Reading %s\n", filename);

    /* Start reading soil file */
    NextLine (geol_file, cmdstr, &lno);
    ReadKeyword (cmdstr, "NUMGEOL", &geoltbl->number, 'i', filename, lno);

    geoltbl->ksatv = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->ksath = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->smcmax = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->smcmin = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->alpha = (double *)malloc (geoltbl->number * sizeof (double));
    geoltbl->beta = (double *)malloc (geoltbl->number * sizeof (double));

    /* Skip header line */
    NextLine (geol_file, cmdstr, &lno);

    for (i = 0; i < geoltbl->number; i++)
    {
        NextLine (geol_file, cmdstr, &lno);
        match = sscanf (cmdstr, "%d %lf %lf %lf %lf %lf %lf",
            &index, &geoltbl->ksatv[i], &geoltbl->ksath[i],
            &geoltbl->smcmax[i], &geoltbl->smcmin[i],
            &geoltbl->alpha[i], &geoltbl->beta[i]);

        if (match != 7 || i != index - 1)
        {
            PIHMprintf (VL_ERROR,
                "Error reading properties of the %dth geology type.\n", i + 1);
            PIHMprintf (VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            PIHMexit (EXIT_FAILURE);
        }
    }

    fclose (geol_file);
}
