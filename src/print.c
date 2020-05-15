#include "pihm.h"

#define WB_HEADER           "VARIABLES = \"TIME (s)\" \"SRC (m)\" \"SNK (m)\" \"STRG (m)\""

void StartupScreen(void)
{
    PIHMprintf(VL_NORMAL, "\n");
    PIHMprintf(VL_NORMAL, "    ########    ####   ##     ##   ##     ##\n");
    PIHMprintf(VL_NORMAL, "    ##     ##    ##    ##     ##   ###   ###\n");
    PIHMprintf(VL_NORMAL, "    ##     ##    ##    ##     ##   #### ####\n");
    PIHMprintf(VL_NORMAL, "    ########     ##    #########   ## ### ##\n");
    PIHMprintf(VL_NORMAL, "    ##           ##    ##     ##   ##     ##\n");
    PIHMprintf(VL_NORMAL, "    ##           ##    ##     ##   ##     ##\n");
    PIHMprintf(VL_NORMAL, "    ##          ####   ##     ##   ##     ##\n");

    PIHMprintf(VL_BRIEF, "\n");
    PIHMprintf(VL_BRIEF, "    The Penn State Integrated Hydrologic Model\n");
    PIHMprintf(VL_BRIEF, "          Version %s\n\n", VERSION);
#if defined(_NOAH_)
    PIHMprintf(VL_BRIEF, "    * Land surface module turned on.\n");
#endif
#if defined(_RT_)
    PIHMprintf(VL_BRIEF, "    * Reactive transport module turned on.\n");
#endif
#if defined(_BGC_)
    PIHMprintf(VL_BRIEF, "    * Biogeochemistry module turned on.\n");
#endif
#if defined(_CYCLES_OBSOLETE_)
    PIHMprintf(VL_BRIEF, "    * Crop module turned on.\n");
#endif
#if defined(_FBR_)
    PIHMprintf(VL_BRIEF, "    * Deep groundwater module turned on.\n");
#endif
#if defined(_OPENMP)
    PIHMprintf(VL_BRIEF, "    * OpenMP (# of threads = %d).\n", nthreads);
#endif
    PIHMprintf(VL_BRIEF, "\n");

    if (1 == corr_mode)
    {
        PIHMprintf(VL_NORMAL,
            "    Surface elevation correction mode turned on.\n");
    }
    if (1 == debug_mode)
    {
        PIHMprintf(VL_NORMAL,
            "    Debug mode turned on.\n");
    }
    if (VL_BRIEF == verbose_mode)
    {
        PIHMprintf(VL_NORMAL,
            "    Brief mode turned on.\n");
    }
    if (VL_VERBOSE == verbose_mode)
    {
        PIHMprintf(VL_NORMAL,
            "    Verbose mode turned on.\n");
    }
    if (1 == append_mode)
    {
        PIHMprintf(VL_NORMAL,
            "    Append mode turned on.\n");
    }
}

#if defined(_TGM_) && defined(_RT_)
void InitOutputFile(const char *outputdir, int watbal, int ascii,
    const chemtbl_struct chemtbl[], const rttbl_struct *rttbl,
    print_struct *print)
#else
void InitOutputFile(const char *outputdir, int watbal, int ascii,
    print_struct *print)
#endif
{
    char            ascii_fn[MAXSTRING];
    char            dat_fn[MAXSTRING];
    char            watbal_fn[MAXSTRING];
    char            perf_fn[MAXSTRING];
    int             i;
    char            mode[2];
    char            bin_mode[3];

    if (append_mode)
    {
        strcpy(mode, "a");
        strcpy(bin_mode, "ab");
    }
    else
    {
        strcpy(mode, "w");
        strcpy(bin_mode, "wb");
    }

    /* Initialize water balance file*/
    if (watbal)
    {
        sprintf(watbal_fn, "%s%s.watbal.plt", outputdir, project);
        print->watbal_file = PIHMfopen(watbal_fn, mode);
    }

    /* Initialize cvode output files */
    if (debug_mode)
    {
        sprintf(perf_fn, "%s%s.cvode.log", outputdir, project);
        print->cvodeperf_file = PIHMfopen(perf_fn, mode);
        /* Print header lines */
        fprintf(print->cvodeperf_file,
            "%-8s%-8s%-16s%-8s%-8s%-8s%-8s%-8s%-8s\n",
            "step", "cpu_dt", "cputime", "maxstep",
            "nsteps", "niters", "nevals", "nefails", "ncfails");
    }

    /*
     * Initialize model variable output files
     */
    for (i = 0; i < print->nprint; i++)
    {
        sprintf(dat_fn, "%s.dat", print->varctrl[i].name);
        print->varctrl[i].datfile = PIHMfopen(dat_fn, bin_mode);

        if (ascii)
        {
            sprintf(ascii_fn, "%s.txt", print->varctrl[i].name);
            print->varctrl[i].txtfile = PIHMfopen(ascii_fn, mode);
        }
    }

#if defined(_TGM_)
    int             n = 0;
    int             k;

    for (i = 0; i < 2; i++)
    {
        /* wflux file header */
        fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "INFIL");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "RECHG");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "EC");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "ETT");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "EDIR");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "SUBSURF1");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "SUBSURF2");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "SUBSURF3");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "OVLFLOW1");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "OVLFLOW2");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "OVLFLOW3");
# if defined(_FBR_)
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "DEEPINFIL");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "DEEPRECHG");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "DEEPLAT1");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "DEEPLAT2");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "DEEPLAT3");
#endif
        fprintf(print->varctrl[n].txtfile, "\n");

        fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m3/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m3/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m3/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m3/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m3/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m3/s");
# if defined(_FBR_)
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m3/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m3/s");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m3/s");
#endif
        fprintf(print->varctrl[n].txtfile, "\n");
        fflush(print->varctrl[n].txtfile);
        n++;

        /* wstate file header */
        fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "IS");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "SNEQV");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "SURF");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "UNSAT");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "GW");
# if defined(_FBR_)
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "DEEPUNSAT");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "DEEPGW");
#endif
        fprintf(print->varctrl[n].txtfile, "\n");

        fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m");
# if defined(_FBR_)
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m");
#endif
        fprintf(print->varctrl[n].txtfile, "\n");
        fflush(print->varctrl[n].txtfile);
        n++;

        /* Soil moisture content file header */
        fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
        for (k = 0; k < MAXLYR; k++)
        {
            fprintf(print->varctrl[n].txtfile, "\tSMC_L%-4d", k + 1);
        }
        fprintf(print->varctrl[n].txtfile, "\n");

        fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
        for (k = 0; k < MAXLYR; k++)
        {
            fprintf(print->varctrl[n].txtfile, "\t%-9s", "m3/m3");
        }
        fprintf(print->varctrl[n].txtfile, "\n");
        fflush(print->varctrl[n].txtfile);
        n++;

        /* Soil water content file header */
        fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
        for (k = 0; k < MAXLYR; k++)
        {
            fprintf(print->varctrl[n].txtfile, "\tSWC_L%-4d", k + 1);
        }
        fprintf(print->varctrl[n].txtfile, "\n");

        fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
        for (k = 0; k < MAXLYR; k++)
        {
            fprintf(print->varctrl[n].txtfile, "\t%-9s", "m3/m3");
        }
        fprintf(print->varctrl[n].txtfile, "\n");
        fflush(print->varctrl[n].txtfile);
        n++;

        /* Soil temperature file header */
        fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
        for (k = 0; k < MAXLYR; k++)
        {
            fprintf(print->varctrl[n].txtfile, "\tSTC_L%-4d", k + 1);
        }
        fprintf(print->varctrl[n].txtfile, "\n");

        fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
        for (k = 0; k < MAXLYR; k++)
        {
            fprintf(print->varctrl[n].txtfile, "\t%-9s", "K");
        }
        fprintf(print->varctrl[n].txtfile, "\n");
        fflush(print->varctrl[n].txtfile);
        n++;

        /* Land surface file header */
        fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "T1");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "SNOWH");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "ALBEDO");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "LE");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "SHEAT");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "GHEAT");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "ETP");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "ESNOW");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "ROOTW");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "SOILM");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "SOLDN");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "CH");
        fprintf(print->varctrl[n].txtfile, "\n");

        fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "K");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "-");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "W/m2");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "W/m2");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "W/m2");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "W/m2");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "W/m2");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "-");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "W/m2");
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m/s");
        fprintf(print->varctrl[n].txtfile, "\n");
        fflush(print->varctrl[n].txtfile);
        n++;

# if defined(_RT_)
        /* Soil concentration file header */
        fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
        for (k = 0; k < rttbl->num_stc + rttbl->num_ssc; k++)
        {
            char                chemn[MAXSTRING];
            Unwrap(chemn, chemtbl[k].name);

            fprintf(print->varctrl[n].txtfile, "\t%-9s", chemn);
        }
        fprintf(print->varctrl[n].txtfile, "\n");

        fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
        for (k = 0; k < rttbl->num_stc + rttbl->num_ssc; k++)
        {
            fprintf(print->varctrl[n].txtfile, "\t%-9s", "mole/L");
        }
        fprintf(print->varctrl[n].txtfile, "\n");
        fflush(print->varctrl[n].txtfile);
        n++;

#  if defined(_FBR_)
        /* Deep aquifer concentration file header */
        fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
        for (k = 0; k < rttbl->num_stc + rttbl->num_ssc; k++)
        {
            char                chemn[MAXSTRING];
            Unwrap(chemn, chemtbl[k].name);

            fprintf(print->varctrl[n].txtfile, "\t%-9s", chemn);
        }
        fprintf(print->varctrl[n].txtfile, "\n");

        fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
        for (k = 0; k < rttbl->num_stc + rttbl->num_ssc; k++)
        {
            fprintf(print->varctrl[n].txtfile, "\t%-9s", "mole/L");
        }
        fprintf(print->varctrl[n].txtfile, "\n");
        fflush(print->varctrl[n].txtfile);
        n++;
#  endif
# endif
    }

    /* River flux file header */
    fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
    for (k = 0; k < NUM_RIVFLX; k++)
    {
        fprintf(print->varctrl[n].txtfile, "\tFLUX%-5d", k);
    }
    fprintf(print->varctrl[n].txtfile, "\n");

    fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
    for (k = 0; k < NUM_RIVFLX; k++)
    {
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "m3/s");
    }
    fprintf(print->varctrl[n].txtfile, "\n");
    fflush(print->varctrl[n].txtfile);
    n++;

    /* River state header */
    fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
    fprintf(print->varctrl[n].txtfile, "\t%-9s", "STAGE");
    fprintf(print->varctrl[n].txtfile, "\n");

    fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
    fprintf(print->varctrl[n].txtfile, "\t%-9s", "m");
    fprintf(print->varctrl[n].txtfile, "\n");
    fflush(print->varctrl[n].txtfile);
    n++;

# if defined(_RT_)
    /* River concentration file header */
    fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
    for (k = 0; k < rttbl->num_stc; k++)
    {
        char                chemn[MAXSTRING];
        Unwrap(chemn, chemtbl[k].name);

        fprintf(print->varctrl[n].txtfile, "\t%-9s", chemn);
    }
    fprintf(print->varctrl[n].txtfile, "\n");

    fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
    for (k = 0; k < rttbl->num_stc; k++)
    {
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "mole/L");
    }
    fprintf(print->varctrl[n].txtfile, "\n");
    fflush(print->varctrl[n].txtfile);
    n++;

    /* River leaching file header */
    fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
    for (k = 0; k < rttbl->num_stc; k++)
    {
        char                chemn[MAXSTRING];
        Unwrap(chemn, chemtbl[k].name);

        fprintf(print->varctrl[n].txtfile, "\t%-9s", chemn);
    }
    fprintf(print->varctrl[n].txtfile, "\n");

    fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
    for (k = 0; k < rttbl->num_stc; k++)
    {
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "kmole/s");
    }
    fprintf(print->varctrl[n].txtfile, "\n");
    fflush(print->varctrl[n].txtfile);
    n++;

#  if defined(_FBR_)
    /* River deep leaching file header */
    fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
    for (k = 0; k < rttbl->num_stc; k++)
    {
        char                chemn[MAXSTRING];
        Unwrap(chemn, chemtbl[k].name);

        fprintf(print->varctrl[n].txtfile, "\t%-9s", chemn);
    }
    fprintf(print->varctrl[n].txtfile, "\n");

    fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
    for (k = 0; k < rttbl->num_stc; k++)
    {
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "kmole/s");
    }
    fprintf(print->varctrl[n].txtfile, "\n");
    fflush(print->varctrl[n].txtfile);
    n++;

    /* River deep leaching file header */
    fprintf(print->varctrl[n].txtfile, "%-18s",  "TIME");
    for (k = 0; k < rttbl->num_stc; k++)
    {
        char                chemn[MAXSTRING];
        Unwrap(chemn, chemtbl[k].name);

        fprintf(print->varctrl[n].txtfile, "\t%-9s", chemn);
    }
    fprintf(print->varctrl[n].txtfile, "\n");

    fprintf(print->varctrl[n].txtfile, "%-18s",  "\"YYYY-MM-DD hh:mm\"");
    for (k = 0; k < rttbl->num_stc; k++)
    {
        fprintf(print->varctrl[n].txtfile, "\t%-9s", "kmole/s");
    }
    fprintf(print->varctrl[n].txtfile, "\n");
    fflush(print->varctrl[n].txtfile);
    n++;
#  endif
# endif
#endif
}

void UpdPrintVar(varctrl_struct *varctrl, int nprint, int module_step)
{
    int             i;
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nprint; i++)
    {
        int             j;

        if (varctrl[i].upd_intvl == module_step)
        {
            for (j = 0; j < varctrl[i].nvar; j++)
            {
                varctrl[i].buffer[j] += *varctrl[i].var[j];
            }

            varctrl[i].counter++;
        }
    }
}

void PrintData(varctrl_struct *varctrl, int nprint, int t, int lapse, int ascii)
{
    int             i;
    pihm_t_struct   pihm_time;

    pihm_time = PIHMTime(t);

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nprint; i++)
    {
        int             j;
        double          outval;
        double          outtime;

        if(PrintNow(varctrl[i].intvl, lapse, &pihm_time))
        {
            if (ascii)
            {
                fprintf(varctrl[i].txtfile, "\"%s\"", pihm_time.str);
                for (j = 0; j < varctrl[i].nvar; j++)
                {
                    outval = (varctrl[i].counter > 0) ?
                        varctrl[i].buffer[j] / (double)varctrl[i].counter :
                        varctrl[i].buffer[j];

                    fprintf(varctrl[i].txtfile,
                        (roundi(outval) == BADVAL) ? "\t%-8.0lf" :
                        ((outval == 0.0 || fabs(outval) > 1.0E-3) ?
                        "\t%lf" : "\t%.2le"), outval);
                }
                fprintf(varctrl[i].txtfile, "\n");
                fflush(varctrl[i].txtfile);
            }

            outtime = (double)t;
            fwrite(&outtime, sizeof(double), 1, varctrl[i].datfile);
            for (j = 0; j < varctrl[i].nvar; j++)
            {
                outval = (varctrl[i].counter > 0) ?
                    varctrl[i].buffer[j] / (double)varctrl[i].counter :
                    varctrl[i].buffer[j];

                fwrite(&outval, sizeof(double), 1, varctrl[i].datfile);

                varctrl[i].buffer[j] = 0.0;
            }
            varctrl[i].counter = 0;
            fflush(varctrl[i].datfile);
        }
    }
}

void PrintInit(const elem_struct *elem, const river_struct *river,
    const char *outputdir, int t, int starttime, int endtime, int intvl)
{
    pihm_t_struct   pihm_time;

    pihm_time = PIHMTime(t);

    if(PrintNow(intvl, t - starttime, &pihm_time) || t == endtime)
    {
        FILE           *init_file;
        char            fn[MAXSTRING];
        int             i;

        sprintf(fn, "%s/restart/%s.%s.ic", outputdir, project,
            pihm_time.strshort);

        init_file = PIHMfopen(fn, "wb");

        for (i = 0; i < nelem; i++)
        {
#if defined(_CYCLES_OBSOLETE_)
            fwrite(&elem[i].ws.flatResidueWater, sizeof(double), 1, init_file);
#else
            fwrite(&elem[i].ws.cmc,       sizeof(double), 1, init_file);
#endif
            fwrite(&elem[i].ws.sneqv,     sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.surf,      sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.unsat,     sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.gw,        sizeof(double), 1, init_file);
#if defined(_FBR_)
            fwrite(&elem[i].ws.fbr_unsat, sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.fbr_gw,    sizeof(double), 1, init_file);
#endif
#if defined(_NOAH_)
            fwrite(&elem[i].es.t1,        sizeof(double), 1, init_file);
            fwrite(&elem[i].ps.snowh,     sizeof(double), 1, init_file);

            int             j;

            for (j = 0; j < MAXLYR; j++)
            {
                fwrite(&elem[i].es.stc[j], sizeof(double), 1, init_file);
            }
            for (j = 0; j < MAXLYR; j++)
            {
                fwrite(&elem[i].ws.smc[j], sizeof(double), 1, init_file);
            }
            for (j = 0; j < MAXLYR; j++)
            {
                fwrite(&elem[i].ws.swc[j], sizeof(double), 1, init_file);
            }
#endif
        }

        for (i = 0; i < nriver; i++)
        {
            fwrite(&river[i].ws.stage, sizeof(double), 1, init_file);
        }

        fflush(init_file);
        fclose(init_file);
    }
}

void PrintPerf(void *cvode_mem, int t, int starttime, double cputime_dt,
    double cputime, double maxstep, FILE *perf_file)
{
    static double   dt;
    static long int nst0, nfe0, nni0, ncfn0, netf0;
    long int        nst, nfe, nni, ncfn, netf;
    int             cv_flag;

    /* Gets the cumulative number of internal steps taken by the solver (total
     * so far) */
    cv_flag = CVodeGetNumSteps(cvode_mem, &nst);
    CheckCVodeFlag(cv_flag);

    /* Gets the number of calls to the user's right-hand side function */
    cv_flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    CheckCVodeFlag(cv_flag);

    /* Gets the number of nonlinear iterations performed */
    cv_flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    CheckCVodeFlag(cv_flag);

    /* Gets the number of nonlinear convergence failures that have occurred */
    cv_flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    CheckCVodeFlag(cv_flag);

    /* Gets the number of local error test failures that have occurred */
    cv_flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
    CheckCVodeFlag(cv_flag);

    fprintf(perf_file, "%-8d%-8.3f%-16.3f%-8.2f",
        t - starttime, cputime_dt, cputime, maxstep);
    fprintf(perf_file, "%-8ld%-8ld%-8ld%-8ld%-8ld\n",
        nst - nst0, nni - nni0, nfe - nfe0, netf - netf0, ncfn - ncfn0);
    fflush(perf_file);

    dt = 0.0;

    nst0 = nst;
    nni0 = nni;
    nfe0 = nfe;
    netf0 = netf;
    ncfn0 = ncfn;

    dt += cputime_dt;
}

void PrintWaterBal(FILE *watbal_file, int t, int tstart, int dt,
    const elem_struct *elem, const river_struct *river)
{
    int             i;
    double          tot_src = 0.0, tot_snk = 0.0, tot_strg = 0.0;
    static double   tot_strg_prev = 0.0;
    static double   error = 0.0;

    if (t == tstart + dt)
    {
        fprintf(watbal_file, "%s\n", WB_HEADER);
    }

    for (i = 0; i < nelem; i++)
    {
        tot_src += elem[i].wf.prcp * elem[i].topo.area * dt;
#if defined(_NOAH_)
        tot_src += (elem[i].wf.dew + elem[i].wf.snomlt) *
            elem[i].topo.area * dt;
#endif

        tot_snk += (elem[i].wf.edir + elem[i].wf.ett + elem[i].wf.ec) *
            elem[i].topo.area * dt;
#if defined(_NOAH_)
        tot_snk += elem[i].wf.esnow * elem[i].topo.area * dt;
#endif

#if defined(_CYCLES_OBSOLETE_)
        tot_strg += (elem[i].ws.flatResidueWater + elem[i].ws.stanResidueWater +
            elem[i].ws.sneqv + elem[i].ws.surf +
            (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity) *
            elem[i].topo.area;
#else
        tot_strg += (elem[i].ws.cmc + elem[i].ws.sneqv + elem[i].ws.surf +
            (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity) *
            elem[i].topo.area;
#endif
    }

    for (i = 0; i < nriver; i++)
    {
        tot_strg += river[i].ws.stage * river[i].topo.area;

        if (river[i].down < 0)
        {
            tot_snk += river[i].wf.rivflow[DOWN_CHANL2CHANL] * dt;
        }
    }

    if (tot_strg_prev != 0.0)
    {
        error += tot_src - tot_snk - (tot_strg - tot_strg_prev);

        fprintf(watbal_file, "%d %lg %lg %lg %lg %lg\n", t - tstart,
            tot_src, tot_snk, tot_strg - tot_strg_prev,
            tot_src - tot_snk - (tot_strg - tot_strg_prev), error);
        fflush(watbal_file);
    }

    tot_strg_prev = tot_strg;
}

void PrintCVodeFinalStats(void *cvode_mem)
{
    int             cv_flag;
    long int        nst;
    long int        nfe;
    long int        netf;
    long int        nni;
    long int        ncfn;

    cv_flag = CVodeGetNumSteps(cvode_mem, &nst);
    CheckCVodeFlag(cv_flag);

    cv_flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    CheckCVodeFlag(cv_flag);

    cv_flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
    CheckCVodeFlag(cv_flag);

    cv_flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    CheckCVodeFlag(cv_flag);

    cv_flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    CheckCVodeFlag(cv_flag);

    PIHMprintf(VL_NORMAL, "\n");
    PIHMprintf(VL_NORMAL,
        "num of steps = %-6ld num of rhs evals = %-6ld\n", nst, nfe);
    PIHMprintf(VL_NORMAL,
        "num of nonlin solv iters = %-6ld "
        "num of nonlin solv conv fails = %-6ld "
        "num of err test fails = %-6ld\n",
        nni, ncfn, netf);
}

int PrintNow(int intvl, int lapse, const pihm_t_struct *pihm_time)
{
    int             print = 0;

    if (intvl != 0)
    {
        switch (intvl)
        {
            case YEARLY_OUTPUT:
                if (pihm_time->month == 1 && pihm_time->day == 1 &&
                    pihm_time->hour == 0 && pihm_time->minute == 0)
                {
                    print = 1;
                }
                break;
            case MONTHLY_OUTPUT:
                if (pihm_time->day == 1 && pihm_time->hour == 0 &&
                    pihm_time->minute == 0)
                {
                    print = 1;
                }
                break;
            case DAILY_OUTPUT:
                if (pihm_time->hour == 0 && pihm_time->minute == 0)
                {
                    print = 1;
                }
                break;
            case HOURLY_OUTPUT:
                if (pihm_time->minute == 0)
                {
                    print = 1;
                }
                break;
            default:
                if (lapse % intvl == 0)
                {
                    print = 1;
                }
        }
    }

    return print;
}

void ProgressBar(double progress)
{
    int             i;
    const int       BAR_LENGTH = 50;
    int             length;

    length = (int)(progress * (double)BAR_LENGTH);

    PIHMprintf(VL_NORMAL, "[");

    for (i = 0; i < length; i++)
    {
        PIHMprintf(VL_NORMAL, "=");
    }
    for (i = length; i < BAR_LENGTH; i++)
    {
        PIHMprintf(VL_NORMAL, " ");
    }

    PIHMprintf(VL_NORMAL, "] %d%%", (int)(progress * 100.0));

    if (length == BAR_LENGTH)
    {
        PIHMprintf(VL_NORMAL, "\n");
    }
}
