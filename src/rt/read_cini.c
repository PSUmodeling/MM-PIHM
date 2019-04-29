#include "pihm.h"

void ReadCini(const char filen[], const chemtbl_struct *chemtbl, int NumStc,
    const calib_struct *cal, atttbl_struct *atttbl, chmictbl_struct *chmictbl,
    elem_struct elem[])
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
#if defined(_FBR_)
    const int       NVOL = 4;
#else
    const int       NVOL = 2;
#endif

    fp = fopen(filen, "r");
    CheckFile(fp, filen);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filen);

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
            PIHMprintf(VL_ERROR, "Error in %s near Line %d.\n", filen, lno);
            PIHMexit(EXIT_FAILURE);
        }
    }

    chmictbl->nic = CountOccurr(fp, "CONDITION");
    FindLine(fp, "BOF", &lno, filen);

    chmictbl->conc = (double **)malloc(chmictbl->nic * sizeof(double *));
    chmictbl->ssa = (double **)malloc(chmictbl->nic * sizeof(double *));

    /*
     * Read intial conditions
     */
    for (i = 0; i < chmictbl->nic; i++)
    {
        chmictbl->conc[i] = (double *)calloc(NumStc, sizeof(double));
        chmictbl->ssa[i] = (double *)calloc(NumStc, sizeof(double));

        FindLine(fp, "CONDITION", &lno, filen);

        for (k = 0; k < NumStc; k++)
        {
            NextLine(fp, cmdstr, &lno);
            sscanf(cmdstr, "%s", temp_str);
            if (strcmp("pH", temp_str) == 0)
            {
                convert = 1;
            }

            ind = FindChem(temp_str, chemtbl, NumStc);
            if (ind < 0)
            {
                PIHMprintf(VL_ERROR, "Error finding chemical %s.\n", temp_str);
                PIHMexit(EXIT_FAILURE);
            }

            if (chemtbl[ind].itype == MINERAL)
            {
                if (sscanf(cmdstr, "%*s %lf %*s %lf",
                    &chmictbl->conc[i][ind], &chmictbl->ssa[i][ind]) !=2)
                {
                    PIHMprintf(VL_ERROR,
                        "Error reading initial condition in %s at Line %d.\n",
                        filen, lno);
                }
            }
            else
            {
                if (sscanf(cmdstr, "%*s %lf", &chmictbl->conc[i][ind]) != 1)
                {
                    PIHMprintf(VL_ERROR,
                        "Error reading initial condition in %s at Line %d.\n",
                        filen, lno);
                }
            }

            chmictbl->conc[i][ind] =
                (strcmp(chemtbl[ind].ChemName, "pH") == 0 && convert == 1) ?
                pow(10, -chmictbl->conc[i][ind]) : chmictbl->conc[i][ind];
        }
    }

    fclose(fp);
}
