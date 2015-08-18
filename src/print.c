#include "pihm.h"

void PrintData (prtctrl_struct *prtctrl, int nprint, int t, int lapse, int dt,
    int ascii)
{
    int             i, j;
    struct tm      *timestamp;
    time_t          rawtime;
    char            ascii_fn[MAXSTRING];
    char            dat_fn[MAXSTRING];
    FILE           *fid;
    double          outval;
    double          outtime;

    for (i = 0; i < nprint; i++)
    {
        for (j = 0; j < prtctrl[i].nvrbl; j++)
        {
            prtctrl[i].buffer[j] += *prtctrl[i].vrbl[j];
        }

        if (lapse % prtctrl[i].intvl == 0 && lapse > 0)
        {
            rawtime = t;
            timestamp = gmtime (&rawtime);

            if (ascii)
            {
                sprintf (ascii_fn, "%s.txt", prtctrl[i].name);
                fid = fopen (ascii_fn, "a");
                if (NULL == fid)
                {
                    printf ("ERROR: opening output files (%s)!\n", ascii_fn);
                    PihmExit (1);
                }
                fprintf (fid, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
                    timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                    timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
                for (j = 0; j < prtctrl[i].nvrbl; j++)
                {
                    if (prtctrl[i].intvl > dt)
                    {
                        fprintf (fid, "\t%lf",
                            prtctrl[i].buffer[j] / ((double)(prtctrl[i].intvl / dt)));
                    }
                    else
                    {
                        fprintf (fid, "\t%lf", prtctrl[i].buffer[j]);
                    }
                }
                fprintf (fid, "\n");
                fflush (fid);
                fclose (fid);
            }

            sprintf (dat_fn, "%s.dat", prtctrl[i].name);
            fid = fopen (dat_fn, "ab");
            if (NULL == fid)
            {
                printf ("ERROR: opening output files (.%s)!\n", prtctrl[i].name);
                PihmExit (1);
            }

            outtime = (double) t;
            fwrite (&outtime, sizeof (double), 1, fid);
            for (j = 0; j < prtctrl[i].nvrbl; j++)
            {
                if (prtctrl[i].intvl > dt)
                {
                    outval =
                        prtctrl[i].buffer[j] / ((double)(prtctrl[i].intvl / dt));
                }
                else
                {
                    outval = prtctrl[i].buffer[j];
                }
                fwrite (&outval, sizeof (double), 1, fid);

                prtctrl[i].buffer[j] = 0.0;
            }
            fflush (fid);
            fclose (fid);
        }
    }
}

void PrtInit (pihm_struct pihm, char *simulation)
{
    FILE           *init_file;
    char            fn[MAXSTRING];
    char            project[MAXSTRING];
    char           *token;
    char            tempname[MAXSTRING];
    int             i;

    strcpy (tempname, simulation);
    if (strstr (tempname, ".") != 0)
    {
        token = strtok (tempname, ".");
        strcpy (project, token);
    }
    else
    {
        strcpy (project, simulation);
    }

    sprintf (fn, "input/%s/%s.init", project, simulation);
    init_file = fopen (fn, "wb");

    for (i = 0; i < pihm->numele; i++)
    {
        fwrite (&pihm->elem[i].intcp, sizeof (double), 1, init_file);
        fwrite (&pihm->elem[i].snow, sizeof (double), 1, init_file);
        fwrite (&pihm->elem[i].surf0, sizeof (double), 1, init_file);
        fwrite (&pihm->elem[i].unsat0, sizeof (double), 1, init_file);
        fwrite (&pihm->elem[i].gw0, sizeof (double), 1, init_file);
    }

    for (i = 0; i < pihm->numriv; i++)
    {
        fwrite (&pihm->riv[i].stage0, sizeof (double), 1, init_file);
        fwrite (&pihm->riv[i].gw0, sizeof (double), 1, init_file);
    }

    fclose (init_file);
}
