#include "pihm.h"

void ReadTracer(const char fn[], ctrl_struct *ctrl, solutetbl_struct *solutetbl)
{
    int             kelem, ksolute;
    FILE           *fp;
    char            cmdstr[MAXSTRING];
    char            buffer[MAXSTRING];
    int             bytes_now;
    int             bytes_consumed = 0;
    int             match;
    int             index;
    int             lno = 0;

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    // Start reading lsm_file
    FindLine(fp, "BOF", &lno, fn);

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "NUMTRACER", 'i', fn, lno, &nsolute);

    solutetbl->amount = (double **)malloc(nsolute * sizeof(double *));
    for (ksolute = 0; ksolute < nsolute; ksolute++)
    {
        solutetbl->amount[ksolute] = (double *)malloc(nelem * sizeof(double));
    }

    NextLine(fp, cmdstr, &lno);
    ReadKeyword(cmdstr, "INDEX", 's', fn, lno, buffer);

    for (ksolute = 0; ksolute < nsolute; ksolute++)
    {
        if (sscanf(buffer + bytes_consumed, "%s%n", solutetbl->name[ksolute], &bytes_now) != 1)
        {
            pihm_error(ERROR, ERR_WRONG_FORMAT, fn, lno);
        }
        bytes_consumed += bytes_now;
    }

    for (kelem = 0; kelem < nelem; kelem++)
    {
        bytes_consumed = 0;
        NextLine(fp, cmdstr, &lno);

        if (sscanf(cmdstr + bytes_consumed, "%d%n", &index, &bytes_now) != 1 || index != kelem + 1)
        {
            pihm_error(ERROR, ERR_WRONG_FORMAT, fn, lno);
        }
        bytes_consumed += bytes_now;

        for (ksolute = 0; ksolute < nsolute; ksolute++)
        {
            if (sscanf(cmdstr + bytes_consumed, "%lf%n", &solutetbl->amount[ksolute][kelem], &bytes_now) != 1)
            {
                pihm_error(ERROR, ERR_WRONG_FORMAT, fn, lno);
            }
            bytes_consumed += bytes_now;
        }
    }

    NextLine(fp, cmdstr, &lno);
    ctrl->prtvrbl[SOLUTE_CTRL] = ReadPrintCtrl(cmdstr, "OUTPUT_INTERVAL", fn, lno);

    fclose(fp);
}
