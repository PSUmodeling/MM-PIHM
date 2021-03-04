#include "pihm.h"

int FindChem(const char chemn[MAXSTRING], const chemtbl_struct  chemtbl[],
    int nsps)
{
    int             i;
    int             ind = BADVAL;

    for (i = 0; i < nsps; i++)
    {
        if (strcmp(chemn, chemtbl[i].name) == 0)
        {
            ind = i;
            break;
        }
    }

    return ind;
}

int MatchWrappedKey(const char cmdstr[], const char key[])
{
    char            optstr[MAXSTRING];

    if (sscanf(cmdstr, "'%[^']'", optstr) != 1)
    {
        return 1;
    }
    else
    {
        Wrap(optstr);
        return (strcmp(optstr, key) == 0) ? 0 : 1;
    }
}
