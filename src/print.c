#include "pihm.h"

void PrintData (prtctrl_struct *prtctrl, int t, int lapse, int dt, int ascii)
{
    int             j;
    struct tm      *timestamp;
    time_t          rawtime;
    char            ascii_fn[MAXSTRING];
    FILE           *fid;
    double          outval;
    double          outtime;

    for (j = 0; j < prtctrl->nvrbl; j++)
    {
        prtctrl->buffer[j] += *prtctrl->vrbl[j];
    }

    if (lapse % prtctrl->intvl == 0 && lapse > 0)
    {
        rawtime = t;
        timestamp = gmtime (&rawtime);

        if (ascii)
        {
            sprintf (ascii_fn, "%s.txt", prtctrl->name);
            fid = fopen (ascii_fn, "a");
            if (NULL == fid)
            {
                printf ("ERROR: opening output files (%s)!\n", ascii_fn);
                exit (1);
            }
            fprintf (fid, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
                timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
            for (j = 0; j < prtctrl->nvrbl; j++)
            {
                if (prtctrl->intvl > dt)
                {
                    fprintf (fid, "\t%lf",
                        prtctrl->buffer[j] / ((double)(prtctrl->intvl / dt)));
                }
                else
                {
                    fprintf (fid, "\t%lf", prtctrl->buffer[j]);
                }
            }
            fprintf (fid, "\n");
            fflush (fid);
            fclose (fid);
        }

        fid = fopen (prtctrl->name, "ab");
        if (NULL == fid)
        {
            printf ("ERROR: opening output files (.%s)!\n", prtctrl->name);
            exit (1);
        }

        outtime = (double) t;
        fwrite (&outtime, sizeof (double), 1, fid);
        for (j = 0; j < prtctrl->nvrbl; j++)
        {
            if (prtctrl->intvl > dt)
            {
                outval =
                    prtctrl->buffer[j] / ((realtype)(prtctrl->intvl / dt));
            }
            else
            {
                outval = prtctrl->buffer[j];
            }
            fwrite (&outval, sizeof (double), 1, fid);
            prtctrl->buffer[j] = 0.0;
        }
        fflush (fid);
        fclose (fid);
    }
}

//void PrintInit (Model_Data DS, char *filename)
//{
//    FILE           *init_file;
//    char           *init_name;
//    int             i;
//
//    init_name = (char *)malloc ((2 * strlen (filename) + 13) * sizeof (char));
//    sprintf (init_name, "input/%s/%s.init", filename, filename);
//    init_file = fopen (init_name, "w");
//    free (init_name);
//
//    for (i = 0; i < DS->NumEle; i++)
//        fprintf (init_file, "%lf\t%lf\t%lf\t%lf\t%lf\n", DS->EleIS[i], DS->EleSnow[i], DS->EleSurf[i], DS->EleUnsat[i], DS->EleGW[i]);
//    for (i = 0; i < DS->NumRiv; i++)
//        fprintf (init_file, "%lf\t%lf\n", DS->RivStg[i], DS->EleGW[i + DS->NumEle]);
//    fclose (init_file);
//}
