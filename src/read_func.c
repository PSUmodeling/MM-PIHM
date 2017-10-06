#include "pihm.h"

int Readable(const char *cmdstr)
{
    int             readable;
    int             i;
    char            ch;

    for (i = 0; i < (int)(strlen(cmdstr)); i++)
    {
        if (cmdstr[i] == 32 || cmdstr[i] == '\t' || cmdstr[i] == ' ')
        {
            continue;
        }
        else
        {
            break;
        }
    }

    if (i >= (int)(strlen(cmdstr)))
    {
        readable = 0;
    }
    else
    {
        ch = cmdstr[i];

        if (ch != '#' && ch != '\n' && ch != '\0' && ch != '\r')
        {
            readable = 1;
        }
        else
        {
            readable = 0;
        }
    }

    return readable;
}

void FindLine(FILE *fid, const char *token, int *lno, const char *filename)
{
    int             success = 0;
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];

    if (strcasecmp("BOF", token) == 0)
    {
        rewind(fid);
        *lno = 0;
        success = 1;
    }
    else
    {
        /* Initialize cmdstr */
        strcpy(cmdstr, "\0");

        while (!feof(fid))
        {
            if (Readable(cmdstr))
            {
                sscanf(cmdstr, "%s", optstr);
                if (strcasecmp(token, optstr) == 0)
                {
                    success = 1;
                    break;
                }
            }

            fgets(cmdstr, MAXSTRING, fid);
            (*lno)++;
        }
    }

    if (!success)
    {
        PIHMprintf(VL_ERROR, "Cannot find required keyword %s.\n", token);
        PIHMprintf(VL_ERROR, "Error reading %s.\n", filename);
        PIHMexit(EXIT_FAILURE);
    }
}

void NextLine(FILE *fid, char *cmdstr, int *lno)
{
    /*
     * Read a non-blank line into cmdstr
     */
    strcpy(cmdstr, "\0");

    while (!Readable(cmdstr))
    {
        if (fgets(cmdstr, MAXSTRING, fid) == NULL)
        {
            strcpy(cmdstr, "EOF");
            break;
        }
        else
        {
            (*lno)++;
        }
    }
}

int CountLine(FILE *fid, char *cmdstr, int num_arg, ...)
{
    /*
     * Count number of non-blank lines between current location to where
     * token occurs
     */
    va_list         valist;
    char            optstr[MAXSTRING];
    char            token[MAXSTRING];
    int             count;
    int             success = 0;
    int             i;

    /* access all the arguments assigned to valist */
    /* Initialize cmdstr */
    strcpy(cmdstr, "\0");
    count = 0;

    while (!feof(fid))
    {
        if (Readable(cmdstr))
        {
            sscanf(cmdstr, "%s", optstr);

            /* initialize valist for num number of arguments */
            va_start(valist, num_arg);
            for (i = 0; i < num_arg; i++)
            {
                strcpy(token, va_arg(valist, char *));
                if (strcasecmp(token, optstr) == 0)
                {
                    success = 1;
                }
            }
            /* clean memory reserved for valist */
            va_end(valist);

            if (success)
            {
                break;
            }
            else
            {
                count++;
            }
        }

        fgets(cmdstr, MAXSTRING, fid);
    }

    return count;
}

void CheckFile(const FILE *fid, const char *fn)
{
    if (NULL == fid)
    {
        PIHMprintf(VL_ERROR, "Error opening %s.\n", fn);
        PIHMexit(EXIT_FAILURE);
    }
}

int ReadTS(const char *cmdstr, int *ftime, double *data, int nvrbl)
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

    switch (type)
    {
        case 'd':
            match = sscanf(buffer, "%s %lf", optstr, (double *)value);
            if (match != 2 || strcasecmp(keyword, optstr) != 0)
            {
                PIHMprintf(VL_ERROR, "Expected keyword \"%s\", "
                    "detected keyword \"%s\".\n", keyword, optstr);
                success = 0;
            }
            break;
        case 'i':
            match = sscanf(buffer, "%s %d", optstr, (int *)value);
            if (match != 2 || strcasecmp(keyword, optstr) != 0)
            {
                PIHMprintf(VL_ERROR, "Expected keyword \"%s\", "
                    "detected keyword \"%s\".\n", keyword, optstr);
                success = 0;
            }
            break;
        case 's':
            match = sscanf(buffer, "%s %[^\n]", optstr, (char *)value);
            if (match != 2 || strcasecmp(keyword, optstr) != 0)
            {
                PIHMprintf(VL_ERROR, "Expected keyword \"%s\", "
                    "detected keyword \"%s\".\n", keyword, optstr);
                success = 0;
            }
            break;
        case 't':
            match = sscanf(buffer, "%s %s %s", optstr, ts1, ts2);
            if (match != 3 || strcasecmp(keyword, optstr) != 0)
            {
                PIHMprintf(VL_ERROR, "Expected keyword \"%s\", "
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
            PIHMprintf(VL_ERROR,
                "Error: Keyword type \'%c\' is not defined.\n", type);
            PIHMexit(EXIT_FAILURE);
    }

    if (!success)
    {
        PIHMprintf(VL_ERROR, "Error reading %s near Line %d.\n", filename, lno);
        PIHMexit(EXIT_FAILURE);
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

    match = sscanf(buffer, "%s %[^\n]", optstr, ctrlstr);
    if (match != 2 || strcasecmp(keyword, optstr) != 0)
    {
        PIHMprintf(VL_ERROR, "Expected keyword \"%s\", "
            "detected keyword \"%s\".\n", keyword, optstr);
        PIHMexit(EXIT_FAILURE);
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
            PIHMprintf(VL_ERROR, "Unknown output control option %s "
                "in %s near Line %d.\n", ctrlstr, filename, lno);
            PIHMexit(EXIT_FAILURE);
        }
    }

    return prtvrbl;
}

int CountOccurr(FILE *fid, const char *token)
{
    /*
     * Count number of occurrence of keyword from the current line to the end
     * of file
     */
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];
    int             count;

    /* Initialize cmdstr */
    strcpy(cmdstr, "\0");
    count = 0;

    while (!feof(fid))
    {
        if (Readable(cmdstr))
        {
            sscanf(cmdstr, "%s", optstr);
            if (strcasecmp(token, optstr) == 0)
            {
                count++;
            }
        }

        fgets(cmdstr, MAXSTRING, fid);
    }

    return count;
}
