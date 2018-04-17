#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "custom_io.h"

void _custom_exit(const char *fn, int lineno, const char *func, int debug,
    int error)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Exiting from %s", func);
    if (debug)
    {
        fprintf(stderr, " (%s, Line %d)", fn, lineno);
    }
    fprintf(stderr, "...\n\n");
    fflush(stderr);

    exit(error);
}

void _custom_printf(const char *fn, int lineno, const char *func, int debug,
    int model_verbosity, int verbosity, const char *fmt, ...)
{
    va_list         va;

    va_start(va, fmt);

    if (VL_ERROR == verbosity)
    {
        vfprintf(stderr, fmt, va);
        if (debug)
        {
            fprintf(stderr, "Printed from %s", func);
            fprintf(stderr, " (%s, Line %d).\n", fn, lineno);
        }
        fflush(stderr);
    }
    else if (verbosity <= model_verbosity)
    {
        vfprintf(stdout, fmt, va);
        if (debug)
        {
            printf("Printed from %s", func);
            printf(" (%s, Line %d).\n", fn, lineno);
        }
        fflush(stderr);
    }

    va_end(va);
}

void CheckFile(const FILE *fp, const char *fn)
{
    if (NULL == fp)
    {
        fprintf(stderr, "Error opening %s.\n", fn);
        fflush(stderr);
        exit(EXIT_FAILURE);
    }
}

int CountLine(FILE *fp, char *cmdstr, int num_arg, ...)
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

    /* Access all the arguments assigned to valist */
    /* Initialize cmdstr */
    strcpy(cmdstr, "\0");
    count = 0;

    while (!feof(fp))
    {
        if (Readable(cmdstr))
        {
            sscanf(cmdstr, "%s", optstr);

            /* Initialize valist for num number of arguments */
            va_start(valist, num_arg);
            for (i = 0; i < num_arg; i++)
            {
                strcpy(token, va_arg(valist, char *));
                if (strcasecmp(token, optstr) == 0)
                {
                    success = 1;
                }
            }
            /* Clean memory reserved for valist */
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

        fgets(cmdstr, MAXSTRING, fp);
    }

    return count;
}

int CountOccurr(FILE *fp, const char *token)
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

    while (!feof(fp))
    {
        if (Readable(cmdstr))
        {
            sscanf(cmdstr, "%s", optstr);
            if (strcasecmp(token, optstr) == 0)
            {
                count++;
            }
        }

        fgets(cmdstr, MAXSTRING, fp);
    }

    return count;
}

void FindLine(FILE *fp, const char *token, int *lno, const char *filename)
{
    int             success = 0;
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];

    if (strcasecmp("BOF", token) == 0)
    {
        rewind(fp);
        *lno = 0;
        success = 1;
    }
    else
    {
        /* Initialize cmdstr */
        strcpy(cmdstr, "\0");

        while (!feof(fp))
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

            fgets(cmdstr, MAXSTRING, fp);
            (*lno)++;
        }
    }

    if (!success)
    {
        fprintf(stderr, "Cannot find required keyword %s.\n", token);
        fprintf(stderr, "Error reading %s.\n", filename);
        fflush(stderr);
        exit(EXIT_FAILURE);
    }
}

void NextLine(FILE *fp, char *cmdstr, int *lno)
{
    /*
     * Read a non-blank line into cmdstr
     */
    strcpy(cmdstr, "\0");

    while (!Readable(cmdstr))
    {
        if (fgets(cmdstr, MAXSTRING, fp) == NULL)
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
