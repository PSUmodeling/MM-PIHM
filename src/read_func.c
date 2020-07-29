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

int ReadKeyword(const char *buffer, const char *keyword, void *value, char type,
    const char *filename, int lno)
{
    int             match;
    char            timestr[MAXSTRING], ts1[MAXSTRING], ts2[MAXSTRING];
    char            optstr[MAXSTRING];
    int             success = 1;

    if (NULL == value)
    {
        match = sscanf(buffer, "%s", optstr);
        if (match != 1 || strcasecmp(keyword, optstr) != 0)
        {
            success = 0;
        }
    }
    else
    {
        switch (type)
        {
            case 'd':
                match = sscanf(buffer, "%s %lf", optstr, (double *)value);
                if (match != 2 || strcasecmp(keyword, optstr) != 0)
                {
                    pihm_printf(VL_ERROR, "Expected keyword \"%s\", "
                        "detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                break;
            case 'i':
                match = sscanf(buffer, "%s %d", optstr, (int *)value);
                if (match != 2 || strcasecmp(keyword, optstr) != 0)
                {
                    pihm_printf(VL_ERROR, "Expected keyword \"%s\", "
                        "detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                break;
            case 's':
                match = sscanf(buffer, "%s %[^\n]", optstr, (char *)value);
                if (match != 2 || strcasecmp(keyword, optstr) != 0)
                {
                    pihm_printf(VL_ERROR, "Expected keyword \"%s\", "
                        "detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                break;
            case 'w':
                match = sscanf(buffer, "%s %s", optstr, (char *)value);
                if (match != 2 || strcasecmp(keyword, optstr) != 0)
                {
                    pihm_printf(VL_ERROR, "Expected keyword \"%s\", "
                        "detected keyword \"%s\".\n", keyword, optstr);
                    success = 0;
                }
                break;
            case 't':
                match = sscanf(buffer, "%s %s %s", optstr, ts1, ts2);
                if (match != 3 || strcasecmp(keyword, optstr) != 0)
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
        pihm_printf(VL_ERROR, "Error reading %s near Line %d.\n", filename, lno);
        pihm_exit(EXIT_FAILURE);
    }

    return success;
}

int ReadPrtCtrl(const char *buffer, const char *keyword, const char *filename,
    int lno)
{
    int             match;
    int             prtvrbl;
    char            ctrlstr[MAXSTRING];
    char            optstr[MAXSTRING];

    match = sscanf(buffer, "%s %s", optstr, ctrlstr);
    if (match != 2 || strcasecmp(keyword, optstr) != 0)
    {
        pihm_printf(VL_ERROR, "Expected keyword \"%s\", "
            "detected keyword \"%s\".\n", keyword, optstr);
        pihm_exit(EXIT_FAILURE);
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
        match = sscanf(ctrlstr, "%d", &prtvrbl);
        if (match != 1)
        {
            pihm_printf(VL_ERROR, "Unknown output control option %s "
                "in %s near Line %d.\n", ctrlstr, filename, lno);
            pihm_exit(EXIT_FAILURE);
        }
    }

    return prtvrbl;
}
