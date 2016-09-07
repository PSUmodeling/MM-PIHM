#include "pihm.h"

void PrintEnKFStatus (int starttime, int endtime)
{
    time_t          rawtime;
    struct tm      *timestamp;

    PIHMprintf (VL_NORMAL, "\nRunning ensemble members from ");

    rawtime = (int)starttime;
    timestamp = gmtime (&rawtime);
    PIHMprintf (VL_NORMAL, "%4.4d-%2.2d-%2.2d %2.2d:%2.2d to ",
        timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday,
        timestamp->tm_hour, timestamp->tm_min);

    rawtime = (int)endtime;
    timestamp = gmtime (&rawtime);
    PIHMprintf (VL_NORMAL, "%4.4d-%2.2d-%2.2d %2.2d:%2.2d\n",
        timestamp->tm_year + 1900, timestamp->tm_mon + 1, timestamp->tm_mday,
        timestamp->tm_hour, timestamp->tm_min);
}

void WriteEnKFOut (char *project, enkf_struct ens, char *outputdir, int t)
{
    int             i, j, k;
    int             id;
    time_t          rawtime;
    struct tm      *timestamp;
    FILE           *fid;
    FILE           *fid1;
    char            varn[MAXSTRING];
    char            fn[MAXSTRING];
    double         *x;
    int             ne;
    double          outtime;

    rawtime = (time_t) t;
    timestamp = gmtime (&rawtime);

    ne = ens->ne;

    x = (double *)malloc (ne * sizeof (double));

    /*
     * Write parameter output
     */
    for (i = 0; i < MAXPARAM; i++)
    {
        if (ens->param[i].update == 1 && ens->update_param == 1)
        {
            for (j = 0; j < ne; j++)
            {
                x[j] = ens->member[j].param[i];
            }

            sprintf (fn, "%s%s.txt", outputdir, ens->param[i].name);
            fid = fopen (fn, "a");
            fprintf (fid, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
                timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                timestamp->tm_mday, timestamp->tm_hour, timestamp->tm_min);

            for (j = 0; j < ne; j++)
            {
                fprintf (fid, "\t%lf", x[j]);
            }
            fprintf (fid, "\n");
            fflush (fid);
            fclose (fid);
        }
    }

    for (i = 0; i < ne; i++)
    {
        /*
         * Write variable/flux output
         */
        for (k = 0; k < MAXVAR; k++)
        {
            if (ens->var[k].dim > 0)
            {
                sprintf (fn, "%s%s.%3.3d.%s.dat",
                    outputdir, project, i + 1, ens->var[k].name);
                fid = fopen (fn, "ab");
                outtime = (double)t;
                fwrite (&outtime, sizeof (double), 1, fid);
                for (j = 0; j < ens->var[k].dim; j++)
                {
                    fwrite (&ens->member[i].var[k][j], sizeof (double), 1,
                        fid);
                }
                fflush (fid);
                fclose (fid);

                if (ens->ascii)
                {
                    sprintf (fn, "%s%s.%3.3d.%s.txt",
                        outputdir, project, i + 1, ens->var[k].name);
                    fid = fopen (fn, "a");
                    fprintf (fid, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
                        timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                        timestamp->tm_mday, timestamp->tm_hour,
                        timestamp->tm_min);
                    for (j = 0; j < ens->var[k].dim; j++)
                    {
                        fprintf (fid, "\t%lf", ens->member[i].var[k][j]);
                    }
                    fprintf (fid, "\n");
                    fflush (fid);
                    fclose (fid);
                }
            }
        }

        /*
         * Write PIHM initial condition
         */
        sprintf (fn, "input/%s/%s.%3.3d.ic", project, project, i + 1);
        fid = fopen (fn, "wb");
        CheckFile (fid, fn);

        sprintf (fn, "input/%s/%s.%3.3d.ic.txt", project, project, i + 1);
        fid1 = fopen (fn, "w");
        CheckFile (fid1, fn);

        for (j = 0; j < ens->numele; j++)
        {
            id = FindVar (ens->var, "is");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "snow");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "surf");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "unsat");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "gw");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "t1");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "snowh");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            for (k = 0; k < MAXLYR; k++)
            {
                sprintf (varn, "stc%d", k);
                id = FindVar (ens->var, varn);
                fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
                fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);
            }
            for (k = 0; k < MAXLYR; k++)
            {
                sprintf (varn, "smc%d", k);
                id = FindVar (ens->var, varn);
                fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
                fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);
            }
            for (k = 0; k < MAXLYR; k++)
            {
                sprintf (varn, "swc%d", k);
                id = FindVar (ens->var, varn);
                fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
                fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);
            }

            fprintf (fid1, "\n");
        }

        for (j = 0; j < ens->numriv; j++)
        {
            id = FindVar (ens->var, "stage");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\t", ens->member[i].var[id][j]);

            id = FindVar (ens->var, "rivgw");
            fwrite (&ens->member[i].var[id][j], sizeof (double), 1, fid);
            fprintf (fid1, "%lf\n", ens->member[i].var[id][j]);
        }
        fflush (fid);
        fclose (fid);
        fflush (fid1);
        fclose (fid1);
    }

    free (x);
}

void WriteParamOutput (int rawtime, enkf_struct ens, int ind, char *outputdir)
{
    char            fn[MAXSTRING];
    time_t          timevar;
    struct tm      *timeinfo;
    FILE           *fid;
    int             i;

    timevar = (time_t) rawtime;
    timeinfo = gmtime (&timevar);

    sprintf (fn, "%s/%s.txt", outputdir, ens->param[ind].name);
    fid = fopen (fn, "w");

    fprintf (fid, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
        timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday,
        timeinfo->tm_hour, timeinfo->tm_min);
    for (i = 0; i < ens->ne; i++)
    {
        fprintf (fid, "\t%lf", ens->member[i].param[ind]);
    }
    fprintf (fid, "\n");
    fflush (fid);
    fclose (fid);
}

void WriteCalFile (enkf_struct ens, char *project)
{
    char            fn[MAXSTRING];
    FILE           *fid;
    int             i, j;

    for (i = 0; i < ens->ne; i++)
    {
        sprintf (fn, "input/%s/%s.%3.3d.calib", project, project, i + 1);
        fid = fopen (fn, "w");

        for (j = 0; j <= RIVSHPCOEFF; j++)
        {
            fprintf (fid, "%-16s%-6lf\n", ens->param[j].name,
                ens->member[i].param[j]);
        }

        fprintf (fid, "\nLSM_CALIBRATION\n");

        for (j = DRIP; j <= THETAW; j++)
        {
            fprintf (fid, "%-16s%-6lf\n", ens->param[j].name,
                ens->member[i].param[j]);
        }

        fprintf (fid, "\nSCENARIO\n");

        for (j = PRCP; j <= SFCTMP; j++)
        {
            fprintf (fid, "%-16s%-6lf\n", ens->param[j].name,
                ens->member[i].param[j]);
        }

        fflush (fid);
        fclose (fid);
    }

}
