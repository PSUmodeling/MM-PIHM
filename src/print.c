/*****************************************************************************
 * File        : print.c
 * Version     : June, 2014
 * Function    : print out model results output files
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0...........................
 * a) This file is downgraded from Version 1.0, as no ancillary results are 
 *    being output
 * b) Only state variables and flux to/in/accross river and its bed are being 
 *    output
 * c) Addition of Average Function to output average variables at regular time
 *    intervals
 ********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "pihm.h"

/*
 * Temporal average of State vectors 
 */
void PrintData (Print_Ctrl PCtrl, realtype tmpt, realtype dt, int Ascii)
{
    int             j;
    struct tm      *timestamp;
    time_t         *rawtime;
    char           *ascii_name;
    FILE           *fpin;
    realtype        outval, outtime;

    rawtime = (time_t *) malloc (sizeof (time_t));

    for (j = 0; j < PCtrl.NumVar; j++)
        PCtrl.buffer[j] = PCtrl.buffer[j] + *PCtrl.PrintVar[j];
    if (((int)tmpt % PCtrl.Interval) == 0)
    {
        *rawtime = (int)tmpt;
        timestamp = gmtime (rawtime);
        outtime = (realtype) (*rawtime);

        if (Ascii)
        {
            ascii_name = (char *)malloc (strlen (PCtrl.name) + 5 * sizeof (char));
            sprintf (ascii_name, "%s.txt", PCtrl.name);
            fpin = fopen (ascii_name, "a");
            if (NULL == fpin)
            {
                printf ("\t ERROR: opening output files (%s)!", ascii_name);
                exit (1);
            }
            fprintf (fpin, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"\t", timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);
            for (j = 0; j < PCtrl.NumVar; j++)
            {
                if ((realtype) PCtrl.Interval > dt)
                    fprintf (fpin, "%lf\t", PCtrl.buffer[j] / ( (realtype) PCtrl.Interval / dt));
                else
                    fprintf (fpin, "%lf\t", PCtrl.buffer[j]);
            }
            fprintf (fpin, "\n");
            fflush (fpin);
            fclose (fpin);
            free (ascii_name);
        }
        fpin = fopen (PCtrl.name, "ab");
        if (NULL == fpin)
        {
            printf ("\t ERROR: opening output files (.%s)!", PCtrl.name);
            exit (1);
        }

        fwrite (&outtime, sizeof (realtype), 1, fpin);
        for (j = 0; j < PCtrl.NumVar; j++)
        {
            if ((realtype) PCtrl.Interval > dt)
                outval = PCtrl.buffer[j] / ( (realtype) PCtrl.Interval / dt);
            else
                outval = PCtrl.buffer[j];
            fwrite (&outval, sizeof (realtype), 1, fpin);
            PCtrl.buffer[j] = 0;
        }
        fflush (fpin);
        fclose (fpin);
    }
    free (rawtime);
}

void PrintInit (Model_Data DS, char *filename)
{
    FILE           *init_file;
    char           *init_name;
    int             i;
    init_name = (char *)malloc ((2 * strlen (filename) + 13) * sizeof (char));
    sprintf (init_name, "input/%s/%s.init", filename);
    init_file = fopen (init_name, "w");
    free (init_name);

    for (i = 0; i < DS->NumEle; i++)
        fprintf (init_file, "%lf\t%lf\t%lf\t%lf\t%lf\n", DS->EleIS[i], DS->EleSnow[i], DS->EleSurf[i], DS->EleUnsat[i], DS->EleGW[i]);
    for (i = 0; i < DS->NumRiv; i++)
        fprintf (init_file, "%lf\t%lf\n", DS->RivStg[i], DS->EleGW[i + DS->NumEle]);
    fclose (init_file);
}
