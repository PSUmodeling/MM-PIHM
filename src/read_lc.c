#include "pihm.h"

void ReadLc(const char filename[], lctbl_struct *lctbl)
{
    FILE           *lc_file;    /* Pointer to .lc file */
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    lc_file = pihm_fopen(filename, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filename);

    /* Start reading land cover file */
    NextLine(lc_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUMLC", 'i', filename, lno, &lctbl->number);

    lctbl->laimax    = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->laimin    = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->vegfrac   = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->albedomin = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->albedomax = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->emissmin  = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->emissmax  = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->z0min     = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->z0max     = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->hs        = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->snup      = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->rgl       = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->rsmin     = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->rough     = (double *)malloc(lctbl->number * sizeof(double));
    lctbl->rzd       = (double *)malloc(lctbl->number * sizeof(double));

    /* Skip header line */
    NextLine(lc_file, cmdstr, &lno);

    for (i = 0; i < lctbl->number; i++)
    {
        NextLine(lc_file, cmdstr, &lno);
        match =
            sscanf(cmdstr,
            "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &index,
            &lctbl->vegfrac[i], &lctbl->rzd[i], &lctbl->rsmin[i],
            &lctbl->rgl[i], &lctbl->hs[i], &lctbl->snup[i],
            &lctbl->laimin[i], &lctbl->laimax[i],
            &lctbl->emissmin[i], &lctbl->emissmax[i],
            &lctbl->albedomin[i], &lctbl->albedomax[i],
            &lctbl->z0min[i], &lctbl->z0max[i], &lctbl->rough[i]);
        if (match != 16 || i != index - 1)
        {
            pihm_printf(VL_ERROR,
                "Error reading properties of the %dth landcover type.\n",
                i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            pihm_exit(EXIT_FAILURE);
        }
    }

    NextLine(lc_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "TOPT_DATA", 'd', filename, lno, &lctbl->topt);

    NextLine(lc_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "CFACTR_DATA", 'd', filename, lno, &lctbl->cfactr);

    NextLine(lc_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "RSMAX_DATA", 'd', filename, lno, &lctbl->rsmax);

    fclose(lc_file);
}
