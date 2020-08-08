#include "pihm.h"

int ReadTs(const char cmdstr[], int nvrbl, int *ftime, double *data)
{
    int             match;
    char            timestr[MAXSTRING], ts1[MAXSTRING], ts2[MAXSTRING];
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i;
    int             success = 1;

    match = sscanf(cmdstr + bytes_consumed, "%s %s%n", ts1, ts2, &bytes_now);
    bytes_consumed += bytes_now;

    if (match != 2)
    {
        success = 0;
    }
    else
    {
        for (i = 0; i < nvrbl; i++)
        {
            match =
                sscanf(cmdstr + bytes_consumed, "%lf%n", &data[i], &bytes_now);
            if (match != 1)
            {
                success = 0;
            }
            bytes_consumed += bytes_now;
        }

        sprintf(timestr, "%s %s", ts1, ts2);
        *ftime = StrTime(timestr);
    }

    return success;
}

int ReadKeyword(const char buffer[], const char keyword[], char type,
    const char fn[], int lno, void *value)
{
    char            timestr[MAXSTRING], ts1[MAXSTRING], ts2[MAXSTRING];
    char            optstr[MAXSTRING];
    int             success = 1;

    if (NULL == value)
    {
        if (sscanf(buffer, "%s", optstr) != 0)
        {
            success = 0;
        }
        else if (strcasecmp(keyword, optstr) != 0)
        {
            pihm_printf(VL_ERROR, "Expected keyword \"%s\", "
                "detected keyword \"%s\".\n", keyword, optstr);
            success = 0;
        }
    }
    else
    {
        switch (type)
        {
            case 'd':
                if (sscanf(buffer, "%s %lf", optstr, (double *)value) != 2)
                {
                    success = 0;
                }
                else if (strcasecmp(keyword, optstr) != 0)
                {
                    pihm_printf(VL_ERROR, "Expected keyword \"%s\", "
                        "detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                break;
            case 'i':
                if (sscanf(buffer, "%s %d", optstr, (int *)value) != 2)
                {
                    success = 0;
                }
                else if (strcasecmp(keyword, optstr) != 0)
                {
                    pihm_printf(VL_ERROR, "Expected keyword \"%s\", "
                        "detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                break;
            case 's':
                if (sscanf(buffer, "%s %[^\n]", optstr, (char *)value) != 2)
                {
                    success = 0;
                }
                else if (strcasecmp(keyword, optstr) != 0)
                {
                    pihm_printf(VL_ERROR, "Expected keyword \"%s\", "
                        "detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                break;
            case 'w':
                if (sscanf(buffer, "%s %s", optstr, (char *)value) != 2)
                {
                    success = 0;
                }
                else if (strcasecmp(keyword, optstr) != 0)
                {
                    pihm_printf(VL_ERROR, "Expected keyword \"%s\", "
                        "detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                break;
            case 't':
                if (sscanf(buffer, "%s %s %s", optstr, ts1, ts2) != 3)
                {
                    success = 0;
                }
                else if (strcasecmp(keyword, optstr) != 0)
                {
                    pihm_printf(VL_ERROR, "Expected keyword \"%s\", "
                        "detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                else
                {
                    sprintf(timestr, "%s %s", ts1, ts2);
                    *((int *)value) = StrTime(timestr);
                }
                break;
            default:
                pihm_printf(VL_ERROR,
                    "Error: Keyword type \'%c\' is not defined.\n", type);
                pihm_exit(EXIT_FAILURE);
        }
    }

    if (!success)
    {
        pihm_error(ERR_WRONG_FORMAT, fn, lno);
    }

    return success;
}

int ReadPrintCtrl(const char buffer[], const char keyword[],
    const char fn[], int lno)
{
    int             prtvrbl;
    char            ctrlstr[MAXSTRING];
    char            optstr[MAXSTRING];

    if (sscanf(buffer, "%s %s", optstr, ctrlstr) != 2)
    {
        pihm_error(ERR_WRONG_FORMAT, fn, lno);
    }
    else if (strcasecmp(keyword, optstr) != 0)
    {
        pihm_printf(VL_ERROR, "Expected keyword \"%s\", "
            "detected keyword \"%s\".\n", keyword, optstr);
        pihm_error(ERR_WRONG_FORMAT, fn, lno);
    }

    if (strcasecmp(ctrlstr, "YEARLY") == 0)
    {
        prtvrbl = YEARLY_OUTPUT;
    }
    else if (strcasecmp(ctrlstr, "MONTHLY") == 0)
    {
        prtvrbl = MONTHLY_OUTPUT;
    }
    else if (strcasecmp(ctrlstr, "DAILY") == 0)
    {
        prtvrbl = DAILY_OUTPUT;
    }
    else if (strcasecmp(ctrlstr, "HOURLY") == 0)
    {
        prtvrbl = HOURLY_OUTPUT;
    }
    else
    {
        if (sscanf(ctrlstr, "%d", &prtvrbl) != 1)
        {
            pihm_printf(VL_ERROR, "Unknown output control option %s.\n",
                ctrlstr);
            pihm_error(ERR_WRONG_FORMAT, fn, lno);
        }
    }

    return prtvrbl;
}

int CheckHeader(const char buffer[], int nvar, ...)
{
    va_list         valist;
    char            token[MAXSTRING];
    char            var[MAXSTRING];
    char            expected[MAXSTRING];
    int             k;
    int             bytes_now;
    int             bytes_consumed = 0;
    int             success = 1;

    /* Initialize valist for num number of arguments */
    va_start(valist, nvar);

    expected[0] = '\0';

    for (k = 0; k < nvar; k++)
    {
        /* Expected header from input */
        strcpy(token, va_arg(valist, char *));

        strcat(expected, token);
        strcat(expected, " ");

        /* Read header from line */
        if (sscanf(buffer + bytes_consumed, "%s%n", var, &bytes_now) != 1)
        {
            success = 0;
        }

        if (strcasecmp(token, var) != 0)
        {
            success = 0;
        }

        bytes_consumed += bytes_now;
    }

    /* Clean memory reserved for valist */
    va_end(valist);

    if (success == 0)
    {
        pihm_printf(VL_ERROR, "Expected header: %s\n", expected);
    }

    return success;
}
