#include "pihm.h"

void ReadRiver(const char fn[], rivtbl_struct *rivtbl,
    shptbl_struct *shptbl, matltbl_struct *matltbl, forc_struct *forc)
{
    int             i, j;
    FILE           *fp;
    char            cmdstr[MAXSTRING];
    char            tempstr[2][MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    /*
     * Read river segment block
     */
    /* Read number of river segments */
    FindLine(fp, "BOF", &lno, fn);
    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUMRIV", 'i', fn, lno, &nriver);

#if defined(_LUMPED_)
    if (nriver != 1)
    {
        pihm_printf(VL_ERROR,
            "Error: Number of river segments should be 1 in lumped mode.\n");
        pihm_exit(EXIT_FAILURE);
    }
#endif

    /* Allocate */
    rivtbl->from  = (int *)malloc(nriver * sizeof(int));
    rivtbl->to    = (int *)malloc(nriver * sizeof(int));
    rivtbl->down  = (int *)malloc(nriver * sizeof(int));
    rivtbl->left  = (int *)malloc(nriver * sizeof(int));
    rivtbl->right = (int *)malloc(nriver * sizeof(int));
    rivtbl->shp   = (int *)malloc(nriver * sizeof(int));
    rivtbl->matl  = (int *)malloc(nriver * sizeof(int));
    rivtbl->bc    = (int *)malloc(nriver * sizeof(int));
    rivtbl->rsvr  = (int *)malloc(nriver * sizeof(int));

    /* Check header line */
    NextLine(fp, cmdstr, &lno);
    if (!CheckHeader(cmdstr, 10, "INDEX", "FROM", "TO", "DOWN", "LEFT", "RIGHT",
        "SHAPE", "MATL", "BC", "RES"))
    {
        pihm_error(ERR_WRONG_FORMAT, fn, lno);
    }

    /* Read river segment information */
    for (i = 0; i < nriver; i++)
    {
        NextLine(fp, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %d %d %d %d %d %d %d %d %d", &index,
            &rivtbl->from[i], &rivtbl->to[i], &rivtbl->down[i],
            &rivtbl->left[i], &rivtbl->right[i], &rivtbl->shp[i],
            &rivtbl->matl[i], &rivtbl->bc[i], &rivtbl->rsvr[i]);
        if (match != 10 || i != index - 1)
        {
            pihm_error(ERR_WRONG_FORMAT, fn, lno);
        }
    }

    /*
     * Read river shape information
     */
    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "SHAPE", 'i', fn, lno, &shptbl->number);

    /* Allocate */
    shptbl->depth = (double *)malloc(shptbl->number * sizeof(double));
    shptbl->intrpl_ord = (int *)malloc(shptbl->number * sizeof(int));
    shptbl->coeff = (double *)malloc(shptbl->number * sizeof(double));

    /* Check header line */
    NextLine(fp, cmdstr, &lno);
    if (!CheckHeader(cmdstr, 4, "INDEX", "DPTH", "OINT", "CWID"))
    {
        pihm_error(ERR_WRONG_FORMAT, fn, lno);
    }

    for (i = 0; i < shptbl->number; i++)
    {
        NextLine(fp, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %lf %d %lf", &index,
            &shptbl->depth[i], &shptbl->intrpl_ord[i], &shptbl->coeff[i]);
        if (match != 4 || i != index - 1)
        {
            pihm_error(ERR_WRONG_FORMAT, fn, lno);
        }
    }

    /*
     * Read river material information
     */
    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "MATERIAL", 'i', fn, lno, &matltbl->number);

    /* Allocate */
    matltbl->rough = (double *)malloc(matltbl->number * sizeof(double));
    matltbl->cwr   = (double *)malloc(matltbl->number * sizeof(double));
    matltbl->ksath = (double *)malloc(matltbl->number * sizeof(double));

    /* Check header line */
    NextLine(fp, cmdstr, &lno);
    if (!CheckHeader(cmdstr, 4, "INDEX", "ROUGH", "CWR", "KH"))
    {
        pihm_error(ERR_WRONG_FORMAT, fn, lno);
    }

    for (i = 0; i < matltbl->number; i++)
    {
        NextLine(fp, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %lf %lf %lf", &index,
            &matltbl->rough[i], &matltbl->cwr[i], &matltbl->ksath[i]);
        if (match != 4 || i != index - 1)
        {
            pihm_error(ERR_WRONG_FORMAT, fn, lno);
        }
    }

    /*
     * Read river boundary condition block
     */
    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "BC", 'i', fn, lno, &forc->nriverbc);

    if (forc->nriverbc > 0)
    {
        forc->riverbc =
            (tsdata_struct *)malloc(forc->nriverbc * sizeof(tsdata_struct));

        NextLine(fp, cmdstr, &lno);
        for (i = 0; i < forc->nriverbc; i++)
        {
            match = sscanf(cmdstr, "%s %d %s %d", tempstr[0], &index,
                tempstr[1], &forc->riverbc[i].bc_type);
            if (match != 4 || i != index - 1 ||
                strcasecmp(tempstr[0], "RIV_TS") != 0 ||
                strcasecmp(tempstr[1], "TYPE") != 0)
            {
                pihm_error(ERR_WRONG_FORMAT, fn, lno);
            }
            if (forc->riverbc[i].bc_type != DIRICHLET &&
                forc->riverbc[i].bc_type != NEUMANN)
            {
                pihm_printf(VL_ERROR, "Boundary condition type should be "
                    "either Dirichlet (1) or Neumann (2).\n");
                pihm_error(ERR_WRONG_FORMAT, fn, lno);
            }
            /* Check header */
            NextLine(fp, cmdstr, &lno);
            if (!CheckHeader(cmdstr, 2, "TIME",
                (forc->riverbc[i].bc_type == DIRICHLET) ? "HEAD" : "FLUX"))
            {
                pihm_error(ERR_WRONG_FORMAT, fn, lno);
            }

            forc->riverbc[i].length =
                CountLine(fp, cmdstr, 2, "RIV_TS", "RES");
        }

        FindLine(fp, "BOF", &lno, fn);
        FindLine(fp, "BC", &lno, fn);

        for (i = 0; i < forc->nriverbc; i++)
        {
            NextLine(fp, cmdstr, &lno);
            NextLine(fp, cmdstr, &lno);

            forc->riverbc[i].data =
                (double **)malloc((forc->riverbc[i].length) * sizeof(double *));
            forc->riverbc[i].ftime =
                (int *)malloc((forc->riverbc[i].length) * sizeof(int));
            for (j = 0; j < forc->riverbc[i].length; j++)
            {
                forc->riverbc[i].data[j] = (double *)malloc(sizeof(double));
                NextLine(fp, cmdstr, &lno);
                if (!ReadTs(cmdstr, 1, &forc->riverbc[i].ftime[j],
                    &forc->riverbc[i].data[j][0]))
                {
                    pihm_error(ERR_WRONG_FORMAT, fn, lno);
                }
            }
        }
    }

#if NOT_YET_IMPLEMENTED
    /* Read Reservoir information */
#endif

    fclose(fp);
}
