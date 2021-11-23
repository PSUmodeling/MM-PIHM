#include "pihm.h"

#define WB_HEADER           "VARIABLES = \"TIME (s)\" \"SRC (m)\" \"SNK (m)\" \"STRG (m)\""

void StartupScreen(void)
{
    pihm_printf(VL_NORMAL, "\n"
        "    ########    ####   ##     ##   ##     ##\n"
        "    ##     ##    ##    ##     ##   ###   ###\n"
        "    ##     ##    ##    ##     ##   #### ####\n"
        "    ########     ##    #########   ## ### ##\n"
        "    ##           ##    ##     ##   ##     ##\n"
        "    ##           ##    ##     ##   ##     ##\n"
        "    ##          ####   ##     ##   ##     ##\n");

    pihm_printf(VL_BRIEF, "\n"
        "    The Penn State Integrated Hydrologic Model\n"
        "          Version %s\n\n", VERSION);
#if defined(_NOAH_)
    pihm_printf(VL_BRIEF, "    * Land surface module turned on.\n");
#endif
#if defined(_RT_)
    pihm_printf(VL_BRIEF, "    * Reactive transport module turned on.\n");
#endif
#if defined(_BGC_)
    pihm_printf(VL_BRIEF, "    * Biogeochemistry module turned on.\n");
#endif
#if defined(_CYCLES_)
    pihm_printf(VL_BRIEF, "    * Agroecosystem module (Cycles Version %s) turned on.\n", CYCLES_VERSION);
# if defined(_AVGN_)
    pihm_printf(VL_BRIEF, "    * Average profile N mode turned on.\n");
# endif
#endif
#if defined(_DGW_)
    pihm_printf(VL_BRIEF, "    * Deep groundwater module turned on.\n");
#endif
#if defined(_OPENMP)
    pihm_printf(VL_BRIEF, "    * OpenMP (# of threads = %d).\n", nthreads);
#endif
    pihm_printf(VL_BRIEF, "\n");

    if (1 == corr_mode)
    {
        pihm_printf(VL_NORMAL, "    Surface elevation correction mode turned on.\n");
    }
    if (1 == debug_mode)
    {
        pihm_printf(VL_NORMAL, "    Debug mode turned on.\n");
    }
    if (VL_BRIEF == verbose_mode)
    {
        pihm_printf(VL_NORMAL, "    Brief mode turned on.\n");
    }
    if (VL_VERBOSE == verbose_mode)
    {
        pihm_printf(VL_NORMAL, "    Verbose mode turned on.\n");
    }
    if (1 == append_mode)
    {
        pihm_printf(VL_NORMAL, "    Append mode turned on.\n");
    }
}

void InitOutputFiles(const char outputdir[], int watbal, int ascii, print_struct *print)
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

    // Initialize water balance file
    if (watbal)
    {
        sprintf(watbal_fn, "%s%s.watbal.plt", outputdir, project);
        print->watbal_file = pihm_fopen(watbal_fn, mode);
    }

    // Initialize CVODE output files
    if (debug_mode)
    {
        sprintf(perf_fn, "%s%s.cvode.log", outputdir, project);
        print->cvodeperf_file = pihm_fopen(perf_fn, mode);
        // Print header lines
        fprintf(print->cvodeperf_file, "%-8s%-8s%-16s%-8s%-8s%-8s%-8s%-8s%-8s\n",
            "step", "cpu_dt", "cputime", "maxstep", "nsteps", "niters", "nevals", "nefails", "ncfails");
    }

    // Initialize model variable output files
    for (i = 0; i < print->nprint; i++)
    {
        sprintf(dat_fn, "%s.dat", print->varctrl[i].name);
        print->varctrl[i].datfile = pihm_fopen(dat_fn, bin_mode);

        if (ascii)
        {
            sprintf(ascii_fn, "%s.txt", print->varctrl[i].name);
            print->varctrl[i].txtfile = pihm_fopen(ascii_fn, mode);
        }
    }
}

void UpdatePrintVar(int nprint, int module_step, varctrl_struct *varctrl)
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

void PrintData(int nprint, int t, int lapse, int ascii, varctrl_struct *varctrl)
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

        if(PrintNow(varctrl[i].intvl, lapse, pihm_time))
        {
            if (ascii)
            {
                fprintf(varctrl[i].txtfile, "\"%s\"", pihm_time.str);
                for (j = 0; j < varctrl[i].nvar; j++)
                {
                    outval = (varctrl[i].counter > 0) ?
                        varctrl[i].buffer[j] / (double)varctrl[i].counter : varctrl[i].buffer[j];

                    fprintf(varctrl[i].txtfile, (roundi(outval) == BADVAL) ?
                        "\t%-8.0lf" : ((outval == 0.0 || fabs(outval) > 1.0E-3) ? "\t%lf" : "\t%.2le"), outval);
                }
                fprintf(varctrl[i].txtfile, "\n");
                fflush(varctrl[i].txtfile);
            }

            outtime = (double)t;
            fwrite(&outtime, sizeof(double), 1, varctrl[i].datfile);
            for (j = 0; j < varctrl[i].nvar; j++)
            {
                outval = (varctrl[i].counter > 0) ?
                    varctrl[i].buffer[j] / (double)varctrl[i].counter : varctrl[i].buffer[j];

                fwrite(&outval, sizeof(double), 1, varctrl[i].datfile);

                varctrl[i].buffer[j] = 0.0;
            }
            varctrl[i].counter = 0;
            fflush(varctrl[i].datfile);
        }
    }
}

void PrintInit(const char outputdir[], int t, int starttime, int endtime, int intvl, const elem_struct elem[],
    const river_struct river[])
{
    pihm_t_struct   pihm_time;

    pihm_time = PIHMTime(t);

    if(PrintNow(intvl, t - starttime, pihm_time) || t == endtime)
    {
        FILE           *init_file;
        char            fn[MAXSTRING];
        int             i;

        sprintf(fn, "%s/restart/%s.%s.ic", outputdir, project, pihm_time.strshort);

        init_file = pihm_fopen(fn, "wb");

        for (i = 0; i < nelem; i++)
        {
            fwrite(&elem[i].ws.cmc, sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.sneqv, sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.surf, sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.unsat, sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.gw, sizeof(double), 1, init_file);
#if defined(_DGW_)
            fwrite(&elem[i].ws.unsat_geol, sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.gw_geol, sizeof(double), 1, init_file);
#endif
#if defined(_NOAH_)
            fwrite(&elem[i].es.t1, sizeof(double), 1, init_file);
            fwrite(&elem[i].ps.snowh, sizeof(double), 1, init_file);

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

void PrintPerf(int t, int starttime, double cputime_dt, double cputime, double maxstep, FILE *perf_file,
    void *cvode_mem)
{
    static double   dt;
    static long int nst0, nfe0, nni0, ncfn0, netf0;
    long int        nst, nfe, nni, ncfn, netf;
    int             cv_flag;

    // Get the cumulative number of internal steps taken by the solver (total so far)
    cv_flag = CVodeGetNumSteps(cvode_mem, &nst);
    CheckCVodeFlag(cv_flag);

    // Get the number of calls to the user's right-hand side function
    cv_flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    CheckCVodeFlag(cv_flag);

    // Get the number of nonlinear iterations performed
    cv_flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    CheckCVodeFlag(cv_flag);

    // Get the number of nonlinear convergence failures that have occurred
    cv_flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    CheckCVodeFlag(cv_flag);

    // Get the number of local error test failures that have occurred
    cv_flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
    CheckCVodeFlag(cv_flag);

    fprintf(perf_file, "%-8d%-8.3f%-16.3f%-8.2f", t - starttime, cputime_dt, cputime, maxstep);
    fprintf(perf_file, "%-8ld%-8ld%-8ld%-8ld%-8ld\n", nst - nst0, nni - nni0, nfe - nfe0, netf - netf0, ncfn - ncfn0);
    fflush(perf_file);

    dt = 0.0;

    nst0 = nst;
    nni0 = nni;
    nfe0 = nfe;
    netf0 = netf;
    ncfn0 = ncfn;

    dt += cputime_dt;
}

void PrintWaterBalance(int t, int tstart, int dt, const elem_struct elem[], const river_struct river[],
    FILE *watbal_file)
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
        tot_src += (elem[i].wf.dew + elem[i].wf.snomlt) * elem[i].topo.area * dt;
#endif

        tot_snk += (elem[i].wf.edir + elem[i].wf.ett + elem[i].wf.ec) * elem[i].topo.area * dt;
#if defined(_NOAH_)
        tot_snk += elem[i].wf.esnow * elem[i].topo.area * dt;
#endif

        tot_strg += (elem[i].ws.cmc + elem[i].ws.sneqv + elem[i].ws.surf +
            (elem[i].ws.unsat + elem[i].ws.gw) * elem[i].soil.porosity) * elem[i].topo.area;
    }

    for (i = 0; i < nriver; i++)
    {
        tot_strg += river[i].ws.stage * river[i].topo.area;

        if (river[i].down < 0)
        {
            tot_snk += river[i].wf.rivflow[DOWNSTREAM] * dt;
        }
    }

    if (tot_strg_prev != 0.0)
    {
        error += tot_src - tot_snk - (tot_strg - tot_strg_prev);

        fprintf(watbal_file, "%d %lg %lg %lg %lg %lg\n", t - tstart, tot_src, tot_snk, tot_strg - tot_strg_prev,
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

    pihm_printf(VL_NORMAL, "\n");
    pihm_printf(VL_NORMAL, "num of steps = %-6ld num of rhs evals = %-6ld\n", nst, nfe);
    pihm_printf(VL_NORMAL, "num of nonlin solv iters = %-6ld "
        "num of nonlin solv conv fails = %-6ld "
        "num of err test fails = %-6ld\n", nni, ncfn, netf);
}

int PrintNow(int intvl, int lapse, pihm_t_struct pihm_time)
{
    int             print = 0;

    if (intvl != 0)
    {
        switch (intvl)
        {
            case YEARLY_OUTPUT:
                print = (pihm_time.month == 1 && pihm_time.day == 1 && pihm_time.hour == 0 && pihm_time.minute == 0) ?
                    1 : 0;
                break;
            case MONTHLY_OUTPUT:
                print = (pihm_time.day == 1 && pihm_time.hour == 0 && pihm_time.minute == 0) ? 1 : 0;
                break;
            case DAILY_OUTPUT:
                print = (pihm_time.hour == 0 && pihm_time.minute == 0) ? 1 : 0;
                break;
            case HOURLY_OUTPUT:
                print = (pihm_time.minute == 0) ? 1 : 0;
                break;
            default:
                print = (lapse % intvl == 0) ? 1 : 0;
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

    pihm_printf(VL_NORMAL, "[");

    for (i = 0; i < length; i++)
    {
        pihm_printf(VL_NORMAL, "=");
    }
    for (i = length; i < BAR_LENGTH; i++)
    {
        pihm_printf(VL_NORMAL, " ");
    }

    pihm_printf(VL_NORMAL, "] %d%%", (int)(progress * 100.0));

    if (length == BAR_LENGTH)
    {
        pihm_printf(VL_NORMAL, "\n");
    }
}
