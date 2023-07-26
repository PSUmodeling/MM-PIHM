#include "pihm.h"

void ReadAnnualFile(const char fn[], tsdata_struct *ts)
{
    FILE           *fp;
    char            timestr[MAXSTRING];
    char            cmdstr[MAXSTRING];
    int             k;
    int             lno = 0;

    fp = pihm_fopen(fn, "r");
    pihm_printf(VL_VERBOSE, " Reading %s\n", fn);

    ts->length = CountLines(fp, cmdstr, 1, "EOF");
    ts->ftime = (int *)malloc(ts->length * sizeof(int));
    ts->data = (double **)malloc(ts->length * sizeof(double *));

    FindLine(fp, "BOF", &lno, fn);
    for (k = 0; k < ts->length; k++)
    {
        ts->data[k] = (double *)malloc(sizeof(double));
        NextLine(fp, cmdstr, &lno);
        if (sscanf(cmdstr, "%s %lf", timestr, &ts->data[k][0]) != 2)
        {
            pihm_error(ERR_WRONG_FORMAT, fn, lno);
        }

        ts->ftime[k] = StrTime(timestr);
        if (ts->ftime[k] == BADVAL)
        {
            pihm_error(ERR_WRONG_FORMAT, fn, lno);
        }
    }

    fclose(fp);
}
