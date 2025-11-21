#include "pihm.h"

void ReadBc(const char fn[], const atttbl_struct *atttbl, forc_struct *forc)
{
    int             i, j;
    FILE           *fp;
    int             read_bc = 0;
    char            cmdstr[MAXSTRING];
    char            tempstr[2][MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    for (i = 0; i < nelem; i++)
    {
        for (j = 0; j < NUM_EDGE; j++)
        {
            if (atttbl->bc[i][j] != 0)
            {
                read_bc = 1;
                break;
            }
        }
    }

    forc->nbc = 0;

    if (read_bc)
    {
        fp = pihm_fopen(fn, "r");
        pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

        FindLine(fp, "BOF", &lno, fn);

        forc->nbc = CountOccurr(fp, "BC_TS");

        FindLine(fp, "BOF", &lno, fn);
        if (forc->nbc > 0)
        {
            forc->bc = (tsdata_struct *)malloc(forc->nbc * sizeof(tsdata_struct));

            NextLine(fp, cmdstr, &lno);
            for (i = 0; i < forc->nbc; i++)
            {
                match = sscanf(cmdstr, "%s %d %s %d", tempstr[0], &index, tempstr[1], &forc->bc[i].bc_type);
                if (match != 4 || i != index - 1 || strcasecmp(tempstr[0], "BC_TS") != 0 ||
                    strcasecmp(tempstr[1], "TYPE") != 0)
                {
                    pihm_error(ERROR, ERR_WRONG_FORMAT, fn, lno);
                }
                if (forc->bc[i].bc_type != DIRICHLET && forc->bc[i].bc_type != NEUMANN)
                {
                    pihm_printf(VL_ERROR, "Boundary condition type should be either Dirichlet (1) or Neumann (2).\n");
                    pihm_error(ERROR, ERR_WRONG_FORMAT, fn, lno);
                }
                // Check header lines
                NextLine(fp, cmdstr, &lno);
                if (!CheckHeader(cmdstr, 2, "TIME", (forc->bc[i].bc_type == DIRICHLET) ? "HEAD" : "FLUX"))
                {
                    pihm_error(ERROR, ERR_WRONG_FORMAT, fn, lno);
                }
                forc->bc[i].length = CountLines(fp, cmdstr, 1, "BC_TS");
            }

            // Rewind and read
            FindLine(fp, "BOF", &lno, fn);
            for (i = 0; i < forc->nbc; i++)
            {
                // Skip header lines
                NextLine(fp, cmdstr, &lno);
                NextLine(fp, cmdstr, &lno);

                forc->bc[i].ftime = (int *)malloc(forc->bc[i].length * sizeof(int));
                forc->bc[i].data = (double **)malloc(forc->bc[i].length * sizeof(double *));
                for (j = 0; j < forc->bc[i].length; j++)
                {
                    forc->bc[i].data[j] = (double *)malloc(sizeof(double));
                    NextLine(fp, cmdstr, &lno);
                    if (!ReadTs(cmdstr, 1, &forc->bc[i].ftime[j], &forc->bc[i].data[j][0]))
                    {
                        pihm_error(ERROR, ERR_WRONG_FORMAT, fn, lno);
                    }
                }
            }
        }

        fclose(fp);
    }
}
