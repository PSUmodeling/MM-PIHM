#include "custom_io.h"

void _custom_exit(const char *fn, int lno, const char *func, int debug, int error)
{
    if (error != EXIT_SUCCESS)
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Exiting from %s", func);
        if (debug)
        {
            fprintf(stderr, " (%s, Line %d)", fn, lno);
        }
        fprintf(stderr, "...\n\n");
        fflush(stderr);
    }

    exit(error);
}

void _custom_printf(int model_verbosity, int verbosity, const char *fmt, ...)
{
    va_list         va;

    va_start(va, fmt);

    if (VL_ERROR == verbosity)
    {
        vfprintf(stderr, fmt, va);
        fflush(stderr);
    }
    else if (verbosity <= model_verbosity)
    {
        vfprintf(stdout, fmt, va);
        fflush(stderr);
    }

    va_end(va);
}

FILE* _custom_fopen(const char *fn, const char *mode)
{
    FILE           *fp;

#if defined(_MSC_VER)
    if (fopen_s(&fp, fn, mode) != 0)
#else
    fp = fopen(fn, mode);

    if (fp == NULL)
#endif
    {
        fprintf(stderr, "Error opening %s.\n", fn);
        fflush(stderr);
        exit(EXIT_FAILURE);
    }

    return fp;
}

int NonBlank(char *cmdstr)
{
    int             k;
    char            ch;

    // Remove UTF-8 BOM
    if (strncasecmp("\357\273\277", cmdstr, 3) == 0)
    {
        for (k = 0; k < (int)strlen(cmdstr) - 3 + 1; k++)
        {
            cmdstr[k] = cmdstr[k + 3];
        }
    }

    // Go to the first non-whitespace character
    for (k = 0; k < (int)(strlen(cmdstr)); k++)
    {
        if (cmdstr[k] != 32 && cmdstr[k] != '\t')   // ASCII code 32 = space
        {
            break;
        }
    }

    if (k >= (int)(strlen(cmdstr)))
    {
        return 0;
    }
    else
    {
        ch = cmdstr[k];

        if (ch != '#' && ch != '\n' && ch != '\0' && ch != '\r')
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}

// Read a non-blank line into cmdstr
int NextLine(FILE *fp, char *cmdstr, int *lno)
{
    while (fgets(cmdstr, MAXSTRING, fp) != NULL)
    {
        (*lno)++;

        if (NonBlank(cmdstr))
        {
            return 1;
        }
    }

    return 0;    // Return 0 if reaches end of file
}

// Count number of non-blank lines between current location to where token occurs
int CountLines(FILE *fp, char *cmdstr, int num_arg, ...)
{
    va_list         valist;
    char            optstr[MAXSTRING];
    char            token[MAXSTRING];
    int             count = 0;
    int             success = 0;
    int             dummy = 0;
    int             k;

    // Access all the arguments assigned to valist
    while (NextLine(fp, cmdstr, &dummy) != 0)
    {
#if defined(_MSC_VER)
        sscanf_s(cmdstr, "%s", optstr, (unsigned)_countof(optstr));
#else
        sscanf(cmdstr, "%s", optstr);
#endif

        // Initialize valist for num number of arguments
        va_start(valist, num_arg);
        for (k = 0; k < num_arg; k++)
        {
#if defined(_MSC_VER)
            int err = strcpy_s(token, _countof(token), va_arg(valist, char*));

            if (err < 0)
            {
                fprintf(stderr, "\n Error CountLine");
                fflush(stderr);
                exit(EXIT_FAILURE);
            }
#else
            strcpy(token, va_arg(valist, char *));
#endif
            if (strcasecmp(token, optstr) == 0)
            {
                success = 1;
            }
        }
        // Clean memory reserved for valist
        va_end(valist);

        if (success)
        {
            break;
        }

        count++;
    }

    return count;
}

// Count number of occurrence of keyword from the current line to the end of file
int CountOccurr(FILE *fp, const char *token)
{
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];
    int             count = 0;
    int             dummy = 0;

    while (NextLine(fp, cmdstr, &dummy) != 0)
    {
#if defined(_MSC_VER)
        sscanf_s(cmdstr, "%s", optstr, MAXSTRING);
#else
        sscanf(cmdstr, "%s", optstr);
#endif

        if (strcasecmp(token, optstr) == 0)
        {
            count++;
        }
    }

    return count;
}

void FindLine(FILE *fp, const char *token, int *lno, const char *fn)
{
    char            cmdstr[MAXSTRING];
    char            optstr[MAXSTRING];

    if (strcasecmp("BOF", token) == 0)
    {
        rewind(fp);
        *lno = 0;

        return;
    }
    else
    {
        while (NextLine(fp, cmdstr, lno) != 0)
        {
#if defined(_MSC_VER)
            sscanf_s(cmdstr, "%s", optstr, MAXSTRING);
#else
            sscanf(cmdstr, "%s", optstr);
#endif
            if (strcasecmp(token, optstr) == 0)
            {
                return;
            }
        }
    }

    fprintf(stderr, "Cannot find required keyword %s.\n", token);
    fprintf(stderr, "Error reading %s.\n", fn);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

#if !defined(_PIHM_)
void _error(const char *fn, int lno, const char *func, const char *name, int multi, int debug, int error_level, const char *fmt, ...)
#else
void _error(const char *fn, int lno, const char *func, int debug, int error_level, const char *fmt, ...)
#endif
{
    va_list         va;
    char            error_msg[MAXSTRING];

    va_start(va, fmt);

#if defined (_MSC_VER)
    vsprintf_s(error_msg, MAXSTRING, fmt, va);
#else
    vsprintf(error_msg, fmt, va);
#endif

    va_end(va);

    fprintf(stderr, "\n");
#if !defined(_PIHM_)
    if (strcasecmp(fn, "src/main.c") != 0)
    {
        if (multi == 1)
        {
            fprintf(stderr, "Simulation %s:\n", name);
        }
    }
#endif
    if (error_level == ERROR)
    {
        fprintf(stderr, "Error: %s\n", error_msg);
        fflush(stderr);
        _custom_exit(fn, lno, func, debug, EXIT_FAILURE);
    }
    else if (error_level == WARNING)
    {
        fprintf(stderr, "\nWarning: %s\n\n", error_msg);
        fflush(stderr);
    }
}
