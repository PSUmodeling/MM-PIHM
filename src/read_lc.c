#include "pihm.h"

void ReadLc(const char fn[], lctbl_struct *lctbl)
{
    FILE           *fp;
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    /* Start reading land cover file */
    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUMLC", 'i', fn, lno, &lctbl->number);

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

    /* Check header line */
    NextLine(fp, cmdstr, &lno);
    if (!CheckHeader(cmdstr, 16, "INDEX", "SHDFAC", "DROOT", "RS", "RGL", "HS",
        "SNUP", "LAIMIN", "LAIMAX", "EMISMIN", "EMISMAX", "ALBMIN", "ALBMAX",
        "Z0MIN", "Z0MAX", "ROUGH"))
    {
        pihm_error(ERR_WRONG_FORMAT, fn, lno);
    }

    for (i = 0; i < lctbl->number; i++)
    {
        NextLine(fp, cmdstr, &lno);
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
            pihm_error(ERR_WRONG_FORMAT, fn, lno);
        }
    }

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "TOPT_DATA", 'd', fn, lno, &lctbl->topt);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "CFACTR_DATA", 'd', fn, lno, &lctbl->cfactr);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "RSMAX_DATA", 'd', fn, lno, &lctbl->rsmax);

    fclose(fp);
}
