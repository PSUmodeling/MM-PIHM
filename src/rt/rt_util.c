#include "pihm.h"

#define LINE_WIDTH 512
#define WORDS_LINE 40
#define WORD_WIDTH 80

int keymatch(const char *line, const char *keyword, double *value,
    char **strval)
{
    /*
     * A very general and convenient way of reading datafile and input file
     * Find keyword in line, assign the value after keyword to value array if
     * there is any store both numbers and strings in order for later use,
     * buffer required.
     * If is keyword not found return 0. If comments, return 2. Otherwise return
     * 1
     */
    int             i;

    for (i = 0; i < WORDS_LINE; i++)
        value[i] = 0.0;

    if ((line[0] == '!') || (line[0] == '#'))
    {
        /* assign a special flag for comments */
        return (2);
    }

    int             j, k;
    int             words_line = WORDS_LINE;
    int             keyfoundflag = 0;

    char          **words;
    words = (char **)malloc(WORDS_LINE * sizeof(char *));

    for (i = 0; i < WORDS_LINE; i++)
    {
        words[i] = (char *)malloc(WORD_WIDTH * sizeof(char));
        memset(words[i], 0, WORD_WIDTH);
    }
    i = j = k = 0;

    /* Partition the line into words */
    while (i < (int)strlen(line))
    {
        if (line[i] != 39)
        {
            while (line[i] != 9 && line[i] != 0 && line[i] != 10
                && line[i] != 32 && line[i] != 13)
            {
                words[k][j++] = line[i++];
                if (line[i] == 9 || line[i] == 32 || line[i] == 13)
                {
                    k++;
                    j = 0;
                }
            }
        }
        else
        {
            words[k][j++] = line[i++];
            while (line[i] != 39)
            {
                words[k][j++] = line[i++];
            }
            words[k++][j] = line[i++];
            j = 0;
        }
        i++;
    }

    words_line = k + 1;

    for (i = 0; i < words_line; i++)
        if (strcmp(words[i], keyword) == 0)
            keyfoundflag = 1;

    j = k = 0;
    for (i = 0; i < words_line; i++)
    {
        strcpy(strval[k++], words[i]);
        if (realcheck(words[i]) == 1)
            value[j++] = atof(words[i]);
    }

    for (i = 0; i < WORDS_LINE; i++)
        free(words[i]);
    free(words);
    return (keyfoundflag);

}

int realcheck(const char *words)
{
    int             flg = 1, i;
    if (((words[0] >= '0') && (words[0] <= '9')) ||
        (words[0] == '.') || (words[0] == '-') || (words[0] == '+'))
    {
        for (i = 0; i < (int)strlen(words); i++)
        {
            /* Ascii 10 is new line and 13 is carriage return */
            if ((words[i] > '9' || words[i] < '+') && (words[i] != 'E')
                && (words[i] != 'e') && (words[i] != 10) && (words[i] != 13))
            {
                flg = 0;
            }
        }
    }
    else
    {
        flg = 0;
    }
    return (flg);
}

int FindChem(const char chemn[MAXSTRING], const chemtbl_struct  chemtbl[],
    int nsps)
{
    int             i;
    int             ind = -1;

    for (i = 0; i < nsps; i++)
    {
        if (strcmp(chemn, chemtbl[i].ChemName) == 0)
        {
            ind = i;
            break;
        }
    }

    return ind;
}
