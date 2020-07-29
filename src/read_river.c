#include "pihm.h"

void ReadRiver(const char *filename, rivtbl_struct *rivtbl,
    shptbl_struct *shptbl, matltbl_struct *matltbl, forc_struct *forc)
{
    int             i, j;
    FILE           *riv_file;
    char            cmdstr[MAXSTRING];
    int             match;
    int             index;
    int             lno = 0;

    riv_file = pihm_fopen(filename, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", filename);

    /*
     * Read river segment block
     */
    /* Read number of river segments */
    FindLine(riv_file, "BOF", &lno, filename);
    NextLine(riv_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUMRIV", 'i', filename, lno, &nriver);

#if defined(_LUMPED_)
    if (nriver != 1)
    {
        pihm_printf(VL_ERROR,
            "Error: Number of river segments should be 1 in lumped mode.\n");
        pihm_exit(EXIT_FAILURE);
    }
#endif

    /* Allocate */
    rivtbl->from = (int *)malloc(nriver * sizeof(int));
    rivtbl->to = (int *)malloc(nriver * sizeof(int));
    rivtbl->down = (int *)malloc(nriver * sizeof(int));
    rivtbl->left = (int *)malloc(nriver * sizeof(int));
    rivtbl->right = (int *)malloc(nriver * sizeof(int));
    rivtbl->shp = (int *)malloc(nriver * sizeof(int));
    rivtbl->matl = (int *)malloc(nriver * sizeof(int));
    rivtbl->bc = (int *)malloc(nriver * sizeof(int));
    rivtbl->rsvr = (int *)malloc(nriver * sizeof(int));

    /* Skip header line */
    NextLine(riv_file, cmdstr, &lno);

    /* Read river segment information */
    for (i = 0; i < nriver; i++)
    {
        NextLine(riv_file, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %d %d %d %d %d %d %d %d %d",
            &index,
            &rivtbl->from[i], &rivtbl->to[i], &rivtbl->down[i],
            &rivtbl->left[i], &rivtbl->right[i], &rivtbl->shp[i],
            &rivtbl->matl[i], &rivtbl->bc[i], &rivtbl->rsvr[i]);
        if (match != 10 || i != index - 1)
        {
            pihm_printf(VL_ERROR,
                "Error reading river attribute for the %dth segment.\n", i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            pihm_exit(EXIT_FAILURE);
        }
    }

    /*
     * Read river shape information
     */
    NextLine(riv_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "SHAPE", 'i', filename, lno, &shptbl->number);

    /* Allocate */
    shptbl->depth = (double *)malloc(shptbl->number * sizeof(double));
    shptbl->intrpl_ord = (int *)malloc(shptbl->number * sizeof(int));
    shptbl->coeff = (double *)malloc(shptbl->number * sizeof(double));

    /* Skip header line */
    NextLine(riv_file, cmdstr, &lno);

    for (i = 0; i < shptbl->number; i++)
    {
        NextLine(riv_file, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %lf %d %lf",
            &index,
            &shptbl->depth[i], &shptbl->intrpl_ord[i], &shptbl->coeff[i]);
        if (match != 4 || i != index - 1)
        {
            pihm_printf(VL_ERROR,
                "Error reading river shape description for the %dth shape.\n",
                i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            pihm_exit(EXIT_FAILURE);
        }
    }

    /*
     * Read river material information
     */
    NextLine(riv_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "MATERIAL", 'i', filename, lno, &matltbl->number);

    /* Allocate */
    matltbl->rough = (double *)malloc(matltbl->number * sizeof(double));
    matltbl->cwr = (double *)malloc(matltbl->number * sizeof(double));
    matltbl->ksath = (double *)malloc(matltbl->number * sizeof(double));

    /* Skip header line */
    NextLine(riv_file, cmdstr, &lno);

    for (i = 0; i < matltbl->number; i++)
    {
        NextLine(riv_file, cmdstr, &lno);
        match = sscanf(cmdstr, "%d %lf %lf %lf",
            &index, &matltbl->rough[i], &matltbl->cwr[i],
            &matltbl->ksath[i]);
        if (match != 4 || i != index - 1)
        {
            pihm_printf(VL_ERROR,
                "Error reading description of the %dth material.\n", i + 1);
            pihm_printf(VL_ERROR, "Error in %s near Line %d.\n", filename, lno);
            pihm_exit(EXIT_FAILURE);
        }
    }

    /*
     * Read river boundary condition block
     */
    NextLine(riv_file, cmdstr, &lno);
    ReadKeyword(cmdstr, "BC", 'i', filename, lno, &forc->nriverbc);

    if (forc->nriverbc > 0)
    {
        forc->riverbc =
            (tsdata_struct *)malloc(forc->nriverbc * sizeof(tsdata_struct));

        NextLine(riv_file, cmdstr, &lno);
        for (i = 0; i < forc->nriverbc; i++)
        {
            match = sscanf(cmdstr, "%*s %d %*s %d",
                &index, &forc->riverbc[i].bc_type);
            if (match != 2 || i != index - 1)
            {
                pihm_printf(VL_ERROR,
                    "Error reading description "
                    "of the %dth river boundary condition.\n", i);
                pihm_printf(VL_ERROR, "Error in %s near Line %d.\n",
                    filename, lno);
                pihm_exit(EXIT_FAILURE);
            }
            if (forc->riverbc[i].bc_type != DIRICHLET &&
                forc->riverbc[i].bc_type != NEUMANN)
            {
                pihm_printf(VL_ERROR,
                    "Error reading the %dth river boundary condition "
                    "time series.\n", i + 1);
                pihm_printf(VL_ERROR,
                    "Boundary condition type should be "
                    "either Dirichlet (1) or Neumann (2).\n");
                pihm_printf(VL_ERROR, "Error in %s near Line %d.\n",
                    filename, lno);
                pihm_exit(EXIT_FAILURE);
            }
            NextLine(riv_file, cmdstr, &lno);
            NextLine(riv_file, cmdstr, &lno);
            forc->riverbc[i].length =
                CountLine(riv_file, cmdstr, 2, "RIV_TS", "RES");
        }

        FindLine(riv_file, "BOF", &lno, filename);
        FindLine(riv_file, "BC", &lno, filename);

        for (i = 0; i < forc->nriverbc; i++)
        {
            NextLine(riv_file, cmdstr, &lno);
            NextLine(riv_file, cmdstr, &lno);
            NextLine(riv_file, cmdstr, &lno);

            forc->riverbc[i].data =
                (double **)malloc((forc->riverbc[i].length) * sizeof(double *));
            forc->riverbc[i].ftime =
                (int *)malloc((forc->riverbc[i].length) * sizeof(int));
            for (j = 0; j < forc->riverbc[i].length; j++)
            {
                forc->riverbc[i].data[j] = (double *)malloc(sizeof(double));
                NextLine(riv_file, cmdstr, &lno);
                if (!ReadTs(cmdstr, 1, &forc->riverbc[i].ftime[j],
                        &forc->riverbc[i].data[j][0]))
                {
                    pihm_printf(VL_ERROR,
                        "Error reading river boundary condition.\n");
                    pihm_printf(VL_ERROR, "Error in %s near Line %d.\n",
                        filename, lno);
                    pihm_exit(EXIT_FAILURE);
                }
            }
        }
    }

#if NOT_YET_IMPLEMENTED
    /* Read Reservoir information */
#endif

    fclose(riv_file);
}
