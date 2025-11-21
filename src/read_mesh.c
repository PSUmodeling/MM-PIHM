#include "pihm.h"

void ReadMesh(const char *fn, meshtbl_struct *meshtbl)
{
    FILE           *fp;
    int             i;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    // Read element mesh block
    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUMELE", 'i', fn, lno, &nelem);

    meshtbl->node = (int **)malloc(nelem * sizeof(int *));
    meshtbl->nabr = (int **)malloc(nelem * sizeof(int *));

    // Check header line
    NextLine(fp, cmdstr, &lno);
    if (!CheckHeader(cmdstr, 7, "INDEX", "NODE1", "NODE2", "NODE3", "NABR1", "NABR2", "NABR3"))
    {
        pihm_error(ERROR, ERR_WRONG_FORMAT, fn, lno);
    }

    for (i = 0; i < nelem; i++)
    {
        meshtbl->node[i] = (int *)malloc(NUM_EDGE * sizeof(int));
        meshtbl->nabr[i] = (int *)malloc(NUM_EDGE * sizeof(int));

        NextLine(fp, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %d %d %d %d %d %d",
            &index, &meshtbl->node[i][0], &meshtbl->node[i][1], &meshtbl->node[i][2], &meshtbl->nabr[i][0], &meshtbl->nabr[i][1], &meshtbl->nabr[i][2]);
        if (match != 7 || i != index - 1)
        {
            pihm_error(ERROR, ERR_WRONG_FORMAT, fn, lno);
        }
    }

    // Read node block
    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUMNODE", 'i', fn, lno, &meshtbl->numnodes);

    // Check header line
    NextLine(fp, cmdstr, &lno);
    if (!CheckHeader(cmdstr, 5, "INDEX", "X", "Y", "ZMIN", "ZMAX"))
    {
        pihm_error(ERROR, ERR_WRONG_FORMAT, fn, lno);
    }

    meshtbl->x = (double *)malloc(meshtbl->numnodes * sizeof(double));
    meshtbl->y = (double *)malloc(meshtbl->numnodes * sizeof(double));
    meshtbl->zmin = (double *)malloc(meshtbl->numnodes * sizeof(double));
    meshtbl->zmax = (double *)malloc(meshtbl->numnodes * sizeof(double));

    for (i = 0; i < meshtbl->numnodes; i++)
    {
        NextLine(fp, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %lf %lf %lf %lf", &index, &meshtbl->x[i], &meshtbl->y[i], &meshtbl->zmin[i], &meshtbl->zmax[i]);
        if (match != 5 || i != index - 1)
        {
            pihm_error(ERROR, ERR_WRONG_FORMAT, fn, lno);
        }
    }

    fclose(fp);
}
