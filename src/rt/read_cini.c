#include "pihm.h"

void ReadCini(const char fn[], const chemtbl_struct *chemtbl, int num_stc,
    atttbl_struct *atttbl, chmictbl_struct *chmictbl)
{
    FILE           *fp;
    int             i, k;
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    int             match;
    int             index;
    int             ind;
    int             lno = 0;
    int             convert = 0;

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    atttbl->prcpc = (int *)malloc(nelem * sizeof(int));
    atttbl->chem_ic = (int **)malloc(nelem * sizeof(int *));

    /* Check header line */
    NextLine(fp, cmdstr, &lno);
#if defined(_DGW_)
    if (!CheckHeader(cmdstr, 4, "ELEM", "PRCPC", "SOIL", "GEOL"))
#else
    if (!CheckHeader(cmdstr, 3, "ELEM", "PRCPC", "SOIL", "GEOL"))
#endif
    {
        pihm_printf(VL_ERROR, "Cini file header error.\n");
        pihm_exit(EXIT_FAILURE);
    }
    for (i = 0; i < nelem; i++)
    {
        atttbl->chem_ic[i] = (int *)malloc(NCHMVOL * sizeof(int));

        NextLine(fp, cmdstr, &lno);
#if defined(_DGW_)
        match = sscanf(cmdstr, "%d %d %d %d",
            &index, &atttbl->prcpc[i],
            &atttbl->chem_ic[i][SOIL_CHMVOL], &atttbl->chem_ic[i][GEOL_CHMVOL]);
#else
        match = sscanf(cmdstr, "%d %d %d",
            &index, &atttbl->prcpc[i],
            &atttbl->chem_ic[i][SOIL_CHMVOL]);
#endif

        if (match != NCHMVOL + 2)
        {
            pihm_printf(VL_ERROR,
                "Error reading chemistry attribute of the %dth element.\n",
                i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", fn, lno);
            pihm_exit(EXIT_FAILURE);
        }
    }

    chmictbl->nic = CountOccurr(fp, "CONDITION");
    FindLine(fp, "BOF", &lno, fn);

    chmictbl->conc = (double **)malloc(chmictbl->nic * sizeof(double *));
    chmictbl->ssa = (double **)malloc(chmictbl->nic * sizeof(double *));

    /*
     * Read initial conditions
     */
    for (i = 0; i < chmictbl->nic; i++)
    {
        chmictbl->conc[i] = (double *)calloc(num_stc, sizeof(double));
        chmictbl->ssa[i] = (double *)calloc(num_stc, sizeof(double));

        FindLine(fp, "CONDITION", &lno, fn);

        for (k = 0; k < num_stc; k++)
        {
            NextLine(fp, cmdstr, &lno);
            sscanf(cmdstr, "%s", temp_str);
            if (strcmp("pH", temp_str) == 0)
            {
                convert = 1;
            }

            ind = FindChem(temp_str, chemtbl, num_stc);
            if (ind < 0)
            {
                pihm_printf(VL_ERROR, "Error finding chemical %s.\n", temp_str);
                pihm_exit(EXIT_FAILURE);
            }

            if (chemtbl[ind].itype == MINERAL)
            {
                if (sscanf(cmdstr, "%*s %lf %*s %lf",
                    &chmictbl->conc[i][ind], &chmictbl->ssa[i][ind]) !=2)
                {
                    pihm_printf(VL_ERROR,
                        "Error reading initial condition in %s at Line %d.\n",
                        fn, lno);
                }
            }
            else
            {
                if (sscanf(cmdstr, "%*s %lf", &chmictbl->conc[i][ind]) != 1)
                {
                    pihm_printf(VL_ERROR,
                        "Error reading initial condition in %s at Line %d.\n",
                        fn, lno);
                }
            }

            chmictbl->conc[i][ind] =
                (strcmp(chemtbl[ind].name, "pH") == 0 && convert == 1) ?
                pow(10, -chmictbl->conc[i][ind]) : chmictbl->conc[i][ind];
        }
    }

    fclose(fp);
}
