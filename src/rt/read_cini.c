#include "pihm.h"

void ReadCini(const char filen[], const chemtbl_struct *chemtbl, int NumStc,
    vol_conc *Vcele)
{
    FILE           *fp;
    int             i, k;
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    int             nvol;
    int             match;
    int             elem_ind;
    int             river_ind;
    int             nic;
    int             ind;
    int             lno = 0;
    int            *ic_ind;
    double        **conc;
    double        **ssa;

    fp = fopen(filen, "r");
    CheckFile(fp, filen);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filen);

    FindLine(fp, "ELEM", &lno, filen);

#if defined(_FBR_)
    nvol = 4 * nelem + nriver;
#else
    nvol = 2 * nelem + nriver;
#endif

    ic_ind = (int *)malloc(nvol * sizeof(int));

    for (i = 0; i < nelem; i++)
    {
        NextLine(fp, cmdstr, &lno);
#if defined(_FBR_)
        match = sscanf(cmdstr, "%d %d %d %d %d",
            &elem_ind, &ic_ind[RT_UNSAT(i)], &ic_ind[RT_GW(i)],
            &ic_ind[RT_FBR_UNSAT(i)], &ic_ind[RT_FBR_GW(i)]);
        if (match != 5 || elem_ind != i + 1)
        {
            PIHMprintf(VL_ERROR, "Error reading %s at Line %d.\n", filen, lno);
            PIHMexit(EXIT_FAILURE);
        }
#else
        match = sscanf(cmdstr, "%d %d %d",
            &elem_ind, &ic_ind[RT_UNSAT(i)], &ic_ind[RT_GW(i)]);
        if (match != 3 || elem_ind != i + 1)
        {
            PIHMprintf(VL_ERROR, "Error reading %s at Line %d.\n", filen, lno);
            PIHMexit(EXIT_FAILURE);
        }
#endif
    }

    FindLine(fp, "RIVER", &lno, filen);
    for (i = 0; i < nriver; i++)
    {
        NextLine(fp, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %d",
            &river_ind, &ic_ind[RT_RIVER(i)]);
        if (match != 2 || river_ind != i + 1)
        {
            PIHMprintf(VL_ERROR, "Error reading %s at Line %d.\n", filen, lno);
            PIHMexit(EXIT_FAILURE);
        }
    }

    nic = CountOccurr(fp, "CONDITION");
    FindLine(fp, "BOF", &lno, filen);

    conc = (double **)malloc(nic * sizeof(double *));
    ssa = (double **)malloc(nic * sizeof(double *));

    /*
     * Read intial conditions
     */
    for (i = 0; i < nic; i++)
    {
        conc[i] = (double *)calloc(NumStc, sizeof(double));
        ssa[i] = (double *)calloc(NumStc, sizeof(double));

        FindLine(fp, "CONDITION", &lno, filen);

        for (k = 0; k < NumStc; k++)
        {
            NextLine(fp, cmdstr, &lno);
            sscanf(cmdstr, "%s", temp_str);

            ind = FindChem(temp_str, chemtbl, NumStc);
            if (ind < 0)
            {
                PIHMprintf(VL_ERROR, "Error finding chemical %s.\n", temp_str);
                PIHMexit(EXIT_FAILURE);
            }

            if (chemtbl[ind].itype == MINERAL)
            {
                if (sscanf(cmdstr, "%*s %lf %*s %lf", &conc[i][ind], &ssa[i][ind]) !=2 )
                {
                    PIHMprintf(VL_ERROR,
                        "Error reading initial condition in %s at Line %d.\n",
                        filen, lno);
                }
            }
            else
            {
                if (sscanf(cmdstr, "%*s %lf", &conc[i][ind]) != 1)
                {
                    PIHMprintf(VL_ERROR,
                        "Error reading initial condition in %s at Line %d.\n",
                        filen, lno);
                }

                conc[i][ind] = (strcmp(chemtbl[ind].ChemName, "pH") == 0) ?
                    pow(10, -conc[i][ind]) : conc[i][ind];
            }
        }
    }

    /*
     * Assign initial conditions to different volumes
     */
    for (i = 0; i < nvol; i++)
    {
        for (k = 0; k < NumStc; k++)
        {
            Vcele[i].ic.t_conc[k] = conc[ic_ind[i] - 1][k];
            Vcele[i].ic.p_para[k] = ssa[ic_ind[i] - 1][k];
        }
    }

    free(ic_ind);
}
