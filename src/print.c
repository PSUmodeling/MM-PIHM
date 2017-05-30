#include "pihm.h"

void AsciiArt ()
{
    PIHMprintf (VL_NORMAL, "\n");
    PIHMprintf (VL_NORMAL, "\t\t########  #### ##     ## ##     ##\n");
    PIHMprintf (VL_NORMAL, "\t\t##     ##  ##  ##     ## ###   ###\n");
    PIHMprintf (VL_NORMAL, "\t\t##     ##  ##  ##     ## #### ####\n");
    PIHMprintf (VL_NORMAL, "\t\t########   ##  ######### ## ### ##\n");
    PIHMprintf (VL_NORMAL, "\t\t##         ##  ##     ## ##     ##\n");
    PIHMprintf (VL_NORMAL, "\t\t##         ##  ##     ## ##     ##\n");
    PIHMprintf (VL_NORMAL, "\t\t##        #### ##     ## ##     ##\n");
    PIHMprintf (VL_NORMAL,
        "\n\t    The Penn State Integrated Hydrologic Model\n\n");

#ifdef _NOAH_
    PIHMprintf (VL_NORMAL, "\t    * Land surface module turned on.\n");
#endif
#ifdef _RT_
    PIHMprintf (VL_NORMAL, "\t    * Reactive transport module turned on.\n");
#endif
#ifdef _BGC_
    PIHMprintf (VL_NORMAL, "\t    * Biogeochemistry module turned on.\n");
#endif
#ifdef _CYCLES_
    PIHMprintf (VL_NORMAL, "\t    * Crop module turned on.\n");
#endif
#ifdef _OPENMP
    PIHMprintf (VL_NORMAL, "\t    * OpenMP (# of threads = %d).\n", nthreads);
#endif

    PIHMprintf (VL_NORMAL, "\n");
}

void _PIHMprintf (const char *fn, int lineno, const char *func, int verbosity,
    const char *fmt, ...)
{
    va_list         va;

    va_start (va, fmt);

    if (VL_ERROR == verbosity)
    {
        vfprintf (stderr, fmt, va);
        if (debug_mode)
        {
            fprintf (stderr, "Printed from %s", func);
            fprintf (stderr, " (%s, Line %d.)\n", fn, lineno);
        }
        fflush (stderr);
    }
    else if (verbosity <= verbose_mode)
    {
        vfprintf (stdout, fmt, va);
        if (debug_mode)
        {
            printf ("Printed from %s", func);
            printf (" (%s, Line %d.)\n", fn, lineno);
        }
        fflush (stderr);
    }

    va_end (va);
}

void InitOutputFile (prtctrl_struct *prtctrl, int nprint, int ascii)
{
    char            ascii_fn[MAXSTRING];
    char            dat_fn[MAXSTRING];
    int             i;

    for (i = 0; i < nprint; i++)
    {
        sprintf (dat_fn, "%s.dat", prtctrl[i].name);
        prtctrl[i].datfile = fopen (dat_fn, "w");

        if (ascii)
        {
            sprintf (ascii_fn, "%s.txt", prtctrl[i].name);
            prtctrl[i].txtfile = fopen (ascii_fn, "w");
        }
    }
}

void UpdPrintVar (prtctrl_struct *prtctrl, int nprint, int module_step)
{
    int             i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nprint; i++)
    {
        int         j;

        if (prtctrl[i].upd_intvl == module_step)
        {
            for (j = 0; j < prtctrl[i].nvar; j++)
            {
                prtctrl[i].buffer[j] += *prtctrl[i].var[j];
            }

            prtctrl[i].counter++;
        }
    }
}

void PrintData (prtctrl_struct *prtctrl, int nprint, int t, int lapse,
    int ascii)
{
    int             i;
    pihm_t_struct   pihm_time;

    pihm_time = PIHMTime (t);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nprint; i++)
    {
        int         j;
        int         print = 0;
        double      outval;
        double      outtime;

        switch (prtctrl[i].intvl)
        {
            case YEARLY_OUTPUT:
                if (pihm_time.month == 1 && pihm_time.day == 1 &&
                    pihm_time.hour == 0 && pihm_time.minute == 0)
                {
                    print = 1;
                }
                break;
            case MONTHLY_OUTPUT:
                if (pihm_time.day == 1 && pihm_time.hour == 0 &&
                    pihm_time.minute == 0)
                {
                    print = 1;
                }
                break;
            case DAILY_OUTPUT:
                if (pihm_time.hour == 0 && pihm_time.minute == 0)
                {
                    print = 1;
                }
                break;
            case HOURLY_OUTPUT:
                if (pihm_time.minute == 0)
                {
                    print = 1;
                }
                break;
            default:
                if (lapse % prtctrl[i].intvl == 0 && lapse > 0)
                {
                    print = 1;
                }
        }

        if (print)
        {
            if (ascii)
            {
                fprintf (prtctrl[i].txtfile, "\"%s\"", pihm_time.str);
                for (j = 0; j < prtctrl[i].nvar; j++)
                {
                    if (prtctrl[i].counter > 0)
                    {
                        fprintf (prtctrl[i].txtfile, "\t%lf",
                            prtctrl[i].buffer[j] /
                            (double)prtctrl[i].counter);
                    }
                    else
                    {
                        fprintf (prtctrl[i].txtfile, "\t%lf",
                            prtctrl[i].buffer[j]);
                    }
                }
                fprintf (prtctrl[i].txtfile, "\n");
                fflush (prtctrl[i].txtfile);
            }

            outtime = (double)t;
            fwrite (&outtime, sizeof (double), 1, prtctrl[i].datfile);
            for (j = 0; j < prtctrl[i].nvar; j++)
            {
                if (prtctrl[i].counter > 0)
                {
                    outval = prtctrl[i].buffer[j] /
                        (double)prtctrl[i].counter;
                }
                else
                {
                    outval = prtctrl[i].buffer[j];
                }
                fwrite (&outval, sizeof (double), 1, prtctrl[i].datfile);

                prtctrl[i].buffer[j] = 0.0;
            }
            prtctrl[i].counter = 0;
            fflush (prtctrl[i].datfile);
        }
    }
}

void PrtInit (elem_struct *elem, river_struct *river, char *simulation)
{
    FILE           *init_file;
    char            fn[MAXSTRING];
    int             i;
#ifdef _NOAH_
    int             j;
#endif

    sprintf (fn, "input/%s/%s.ic", project, simulation);
    init_file = fopen (fn, "wb");
    CheckFile (init_file, fn);
    PIHMprintf (VL_ERROR, "Writing initial conditions.\n");

    for (i = 0; i < nelem; i++)
    {
        fwrite (&elem[i].ws.cmc, sizeof (double), 1, init_file);
        fwrite (&elem[i].ws.sneqv, sizeof (double), 1, init_file);
        fwrite (&elem[i].ws.surf, sizeof (double), 1, init_file);
        fwrite (&elem[i].ws.unsat, sizeof (double), 1, init_file);
        fwrite (&elem[i].ws.gw, sizeof (double), 1, init_file);
#ifdef _NOAH_
        fwrite (&elem[i].es.t1, sizeof (double), 1, init_file);
        fwrite (&elem[i].ps.snowh, sizeof (double), 1, init_file);
        for (j = 0; j < MAXLYR; j++)
        {
            fwrite (&elem[i].es.stc[j], sizeof (double), 1, init_file);
        }
        for (j = 0; j < MAXLYR; j++)
        {
            fwrite (&elem[i].ws.smc[j], sizeof (double), 1, init_file);
        }
        for (j = 0; j < MAXLYR; j++)
        {
            fwrite (&elem[i].ws.sh2o[j], sizeof (double), 1, init_file);
        }
#endif
    }

    for (i = 0; i < nriver; i++)
    {
        fwrite (&river[i].ws.stage, sizeof (double), 1, init_file);
        fwrite (&river[i].ws.gw, sizeof (double), 1, init_file);
    }

    fclose (init_file);
}
