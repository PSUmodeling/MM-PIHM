#include "pihm.h"

int Readable (char *cmdstr)
{
    int             readable;
    int             i;
    char            ch;

    for (i = 0; i < strlen (cmdstr); i++)
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

    if (i >= strlen (cmdstr))
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

    return (readable);
}

int FindLine (FILE * fid, char *token)
{
    int             success = 0;
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];

    rewind (fid);

    if (strcasecmp ("BOF", token) == 0)
        success = 1;
    else
    {
        /* Initialize cmdstr */
        strcpy (cmdstr, "\0");

        while (!feof (fid))
        {
            if (Readable (cmdstr))
            {
                sscanf (cmdstr, "%s", optstr);
                if (strcasecmp (token, optstr) == 0)
                {
                    success = 1;
                    break;
                }
            }

            fgets (cmdstr, MAXSTRING, fid);
        }
    }

    return (success);
}

void NextLine (FILE * fid, char *cmdstr)
{
    /*
     * Read a non-blank line into cmdstr
     */
    strcpy (cmdstr, "\0");

    while (!Readable (cmdstr))
    {
        if (fgets (cmdstr, MAXSTRING, fid) == NULL)
        {
            strcpy (cmdstr, "EOF");
            break;
        }

    }
}

int CountLine (FILE * fid, char *cmdstr, int num_arg, ...)
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
    strcpy (cmdstr, "\0");
    count = 0;

    while (!feof (fid))
    {
        if (Readable (cmdstr))
        {
            sscanf (cmdstr, "%s", optstr);

            /* initialize valist for num number of arguments */
            va_start (valist, num_arg);
            for (i = 0; i < num_arg; i++)
            {
                strcpy (token, va_arg (valist, char *));
                if (strcasecmp (token, optstr) == 0)
                    success = 1;
            }
            /* clean memory reserved for valist */
            va_end (valist);

            if (success)
                break;
            else
                count++;
        }

        fgets (cmdstr, MAXSTRING, fid);
    }

    return (count);
}

int ReadTS (char *cmdstr, int *ftime, double *data, int nvrbl)
{
    int             match;
    struct tm      *timeinfo;
    int             bytes_now;
    int             bytes_consumed = 0;
    int             i;
    int             success = 1;

    timeinfo = (struct tm *)malloc (sizeof (struct tm));

    if (1 == nvrbl)
    {
        match = sscanf (cmdstr, "%d-%d-%d %d:%d %lf",
            &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday,
            &timeinfo->tm_hour, &timeinfo->tm_min, &data[0]);
        timeinfo->tm_sec = 0;
        if (match != nvrbl + 5)
        {
            success = 0;
        }
        else
        {
            timeinfo->tm_year = timeinfo->tm_year - 1900;
            timeinfo->tm_mon = timeinfo->tm_mon - 1;
            *ftime = timegm (timeinfo);
        }
    }
    else
    {
        match = sscanf (cmdstr + bytes_consumed, "%d-%d-%d %d:%d%n",
            &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday,
            &timeinfo->tm_hour, &timeinfo->tm_min, &bytes_now);
        bytes_consumed += bytes_now;

        if (match != 5)
        {
            success = 0;
        }
        else
        {
            for (i = 0; i < nvrbl; i++)
            {
                match =
                    sscanf (cmdstr + bytes_consumed, "%lf%n", &data[i],
                    &bytes_now);
                if (match != 1)
                {
                    success = 0;
                }
                bytes_consumed += bytes_now;
            }

            timeinfo->tm_year = timeinfo->tm_year - 1900;
            timeinfo->tm_mon = timeinfo->tm_mon - 1;
            timeinfo->tm_sec = 0;
            *ftime = timegm (timeinfo);
        }
    }

    free (timeinfo);

    return (success);
}

int ReadKeyword (char *buffer, char *keyword, void *value, char type)
{
    int             match;
    char            optstr[MAXSTRING];
    int             success = 1;
    struct tm      *timeinfo;

    switch (type)
    {
        case 'd':
            match = sscanf (buffer, "%s %lf", optstr, (double *)value);
            if (match != 2 || strcasecmp (keyword, optstr) != 0)
            {
                printf ("Expected keyword \"%s\", detected keyword \"%s\".\n",
                    keyword, optstr);
                success = 0;
            }
            break;
        case 'i':
            match = sscanf (buffer, "%s %d", optstr, (int *)value);
            if (match != 2 || strcasecmp (keyword, optstr) != 0)
            {
                printf ("Expected keyword \"%s\", detected keyword \"%s\".\n",
                    keyword, optstr);
                success = 0;
            }
            break;
        case 's':
            match = sscanf (buffer, "%s %[^\n]", optstr, (char *)value);
            if (match != 2 || strcasecmp (keyword, optstr) != 0)
            {
                printf ("Expected keyword \"%s\", detected keyword \"%s\".\n",
                    keyword, optstr);
                success = 0;
            }
            break;
        case 't':
            timeinfo = (struct tm *)malloc (sizeof (struct tm));

            match = sscanf (buffer, "%s %d-%d-%d %d:%d", optstr,
                &timeinfo->tm_year, &timeinfo->tm_mon, &timeinfo->tm_mday,
                &timeinfo->tm_hour, &timeinfo->tm_min);
            timeinfo->tm_sec = 0;
            if (match != 6 || strcasecmp (keyword, optstr) != 0)
            {
                printf ("Expected keyword \"%s\", detected keyword \"%s\".\n",
                    keyword, optstr);
                success = 0;
            }
            else
            {
                timeinfo->tm_year = timeinfo->tm_year - 1900;
                timeinfo->tm_mon = timeinfo->tm_mon - 1;
                *((int *)value) = timegm (timeinfo);
            }

            free (timeinfo);
            break;
        default:
            fprintf (stderr, "Error: Keyword type \'%c\' is not defined.\n",
                type);
            PIHMExit (EXIT_FAILURE);
    }

    return (success);
}

int CountOccurance (FILE * fid, char *token)
{
    /*
     * Count number of occurance of keyword from the current line to the end
     * of file
     */

    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];
    int             count;

    /* Initialize cmdstr */
    strcpy (cmdstr, "\0");
    count = 0;

    while (!feof (fid))
    {
        if (Readable (cmdstr))
        {
            sscanf (cmdstr, "%s", optstr);
            if (strcasecmp (token, optstr) == 0)
                count++;
        }

        fgets (cmdstr, MAXSTRING, fid);
    }

    return (count);
}
