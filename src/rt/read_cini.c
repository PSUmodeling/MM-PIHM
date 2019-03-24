#include "pihm.h"

void ReadCini(const char filen[], const chemtbl_struct *chemtbl, int NumStc,
    const calib_struct *cal, elem_struct elem[])
{
    FILE           *fp;
    int             i, k;
    char            cmdstr[MAXSTRING];
    char            temp_str[MAXSTRING];
    int             match;
    int             elem_ind;
    int             nic;
    int             ind;
    int             lno = 0;
    int           **ic_ind;
    int             convert = 0;
    const int       UNSAT_IND = 0;
    const int       GW_IND = 1;
    double        **conc;
    double        **ssa;

    fp = fopen(filen, "r");
    CheckFile(fp, filen);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", filen);

    FindLine(fp, "ELEM", &lno, filen);

    ic_ind = (int **)malloc(nelem * sizeof(int *));

    for (i = 0; i < nelem; i++)
    {
        ic_ind[i] = (int *)malloc(2 * sizeof(int));

        NextLine(fp, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %d %d",
            &elem_ind, &ic_ind[i][UNSAT_IND], &ic_ind[i][GW_IND]);
        if (match != 3 || elem_ind != i + 1)
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
            if (strcmp("pH", temp_str) == 0)
            {
                strcpy(temp_str, "H+");
                convert = 1;
            }
            wrap(temp_str);

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
                ssa[i][ind] *= cal->ssa;
            }
            else
            {
                if (sscanf(cmdstr, "%*s %lf", &conc[i][ind]) != 1)
                {
                    PIHMprintf(VL_ERROR,
                        "Error reading initial condition in %s at Line %d.\n",
                        filen, lno);
                }
            }

            conc[i][ind] *= (strcmp(chemtbl[ind].ChemName, "'DOC'") == 0) ?
                cal->initconc : 1.0;

            conc[i][ind] =
                (strcmp(chemtbl[ind].ChemName, "'H+'") == 0 && convert == 1) ?
                pow(10, -conc[i][ind]) : conc[i][ind];
        }
    }

    /*
     * Assign initial conditions to different volumes
     */
    for (i = 0; i < nelem; i++)
    {
        for (k = 0; k < NumStc; k++)
        {
            elem[i].restart_input.tconc_unsat[k] =
                conc[ic_ind[i][UNSAT_IND] - 1][k];
            elem[i].restart_input.tconc_gw[k] = conc[ic_ind[i][GW_IND] - 1][k];

            elem[i].restart_input.ssa_unsat[k] =
                ssa[ic_ind[i][UNSAT_IND] - 1][k];
            elem[i].restart_input.ssa_gw[k] =
                ssa[ic_ind[i][GW_IND] - 1][k];
        }
    }

    for (i = 0; i < nic; i++)
    {
        free(conc[i]);
        free(ssa[i]);
    }

    for (i = 0; i < nelem; i++)
    {
        free(ic_ind[i]);
    }

    free(conc);
    free(ssa);
    free(ic_ind);
}
