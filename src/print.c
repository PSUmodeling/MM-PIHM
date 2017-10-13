#include "pihm.h"

#define TEC_HEADER          "VARIABLES = \"X\" \"Y\" \"Zmin\" \"Zmax\" \"h\""
#define WB_HEADER           "VARIABLES = \"TIME (s)\" \"Outflow (cms)\" \"Surf2Chan (cms)\" \"AqF2Chan (cms)\" \"Chan_LKG (cms)\" \"Precipitation (cms)\" \"NetPrec (cms)\" \"Infiltration (cms)\" \"Recharge (cms)\" \"E_soil (cms)\" \"ET_plant (cms)\" \"E_canopy (cms)\" \"PET (cms)\" \"ET (cms)\" \"E_Surface (cms)\" \"E_Unsat (cms)\" \"E_GW (cms)\" \"T_Unsat (cms)\" \"T_GW (cms)\""
#define RIVER_TEC_HEADER2   "ZONE T = \"Water Depth River\""
#define RIVER_TEC_HEADER3   "StrandID=1, SolutionTime="
#define ELEM_TEC_HEADER3    "VARSHARELIST = ([1, 2, 3, 4]=1), CONNECTIVITYSHAREZONE = 1"

void AsciiArt()
{
    PIHMprintf(VL_NORMAL, "\n");
    PIHMprintf(VL_NORMAL, "\t########    ####   ##     ##   ##     ##\n");
    PIHMprintf(VL_NORMAL, "\t##     ##    ##    ##     ##   ###   ###\n");
    PIHMprintf(VL_NORMAL, "\t##     ##    ##    ##     ##   #### ####\n");
    PIHMprintf(VL_NORMAL, "\t########     ##    #########   ## ### ##\n");
    PIHMprintf(VL_NORMAL, "\t##           ##    ##     ##   ##     ##\n");
    PIHMprintf(VL_NORMAL, "\t##           ##    ##     ##   ##     ##\n");
    PIHMprintf(VL_NORMAL, "\t##          ####   ##     ##   ##     ##\n");
    PIHMprintf(VL_NORMAL, "\n\tThe Penn State Integrated Hydrologic Model\n\n");
#ifdef _NOAH_
    PIHMprintf(VL_NORMAL, "\t* Land surface module turned on.\n");
#endif
#ifdef _RT_
    PIHMprintf(VL_NORMAL, "\t* Reactive transport module turned on.\n");
#endif
#ifdef _BGC_
    PIHMprintf(VL_NORMAL, "\t* Biogeochemistry module turned on.\n");
#endif
#ifdef _CYCLES_
    PIHMprintf(VL_NORMAL, "\t* Crop module turned on.\n");
#endif
#ifdef _OPENMP
    PIHMprintf(VL_NORMAL, "\t* OpenMP (# of threads = %d).\n", nthreads);
#endif
    PIHMprintf(VL_NORMAL, "\n");
}

void _PIHMprintf(const char *fn, int lineno, const char *func, int verbosity,
    const char *fmt, ...)
{
    va_list         va;

    va_start(va, fmt);

    if (VL_ERROR == verbosity)
    {
        vfprintf(stderr, fmt, va);
        if (debug_mode)
        {
            fprintf(stderr, "Printed from %s", func);
            fprintf(stderr, " (%s, Line %d).\n", fn, lineno);
        }
        fflush(stderr);
    }
    else if (verbosity <= verbose_mode)
    {
        vfprintf(stdout, fmt, va);
        if (debug_mode)
        {
            printf("Printed from %s", func);
            printf(" (%s, Line %d).\n", fn, lineno);
        }
        fflush(stderr);
    }

    va_end(va);
}

void InitOutputFile(print_struct *print, const char *outputdir, int watbal,
    int ascii)
{
    char            ascii_fn[MAXSTRING];
    char            dat_fn[MAXSTRING];
    char            watbal_fn[MAXSTRING];
    char            perf_fn[MAXSTRING];
    int             i;

    /* Initialize water balance file*/
    if (watbal)
    {
        sprintf(watbal_fn, "%s%s.watbal.plt", outputdir, project);
        print->watbal_file = fopen(watbal_fn, "w");
        CheckFile(print->watbal_file, watbal_fn);
    }

    /* Initialize cvode output files */
    if (debug_mode)
    {
        sprintf(perf_fn, "%s%s.cvode.log", outputdir, project);
        print->cvodeperf_file = fopen(perf_fn, "w");
        CheckFile(print->cvodeperf_file, perf_fn);
        /* Print header lines */
        fprintf(print->cvodeperf_file,
            "%-8s%-8s%-16s%-8s%-8s%-8s%-8s%-8s%-8s\n",
            "step", "cpu_dt", "cputime", "maxstep",
            "nsteps", "nevals", "niters", "ncfails", "nefails");
    }

    /*
     * Initialize model variable output files
     */
    for (i = 0; i < print->nprint; i++)
    {
        sprintf(dat_fn, "%s.dat", print->varctrl[i].name);
        print->varctrl[i].datfile = fopen(dat_fn, "w");
        CheckFile(print->varctrl[i].datfile, dat_fn);

        if (ascii)
        {
            sprintf(ascii_fn, "%s.txt", print->varctrl[i].name);
            print->varctrl[i].txtfile = fopen(ascii_fn, "w");
            CheckFile(print->varctrl[i].txtfile, ascii_fn);
        }
    }

    /* Tecplot files */
    if (tecplot)
    {
        for (i = 0; i < print->ntpprint; i++)
        {
            int             j;

            sprintf(dat_fn, "%s.plt", print->tp_varctrl[i].name);
            print->tp_varctrl[i].datfile = fopen(dat_fn, "w");
            CheckFile(print->tp_varctrl[i].datfile, dat_fn);

            if (print->tp_varctrl[i].intr == 0)
            {
                fprintf(print->tp_varctrl[i].datfile, "%s \n", TEC_HEADER);
                fprintf(print->tp_varctrl[i].datfile,
                    "ZONE T=\"%s\", N=%d, E=%d, DATAPACKING=%s, "
                    "SOLUTIONTIME=%lf, ZONETYPE=%s\n",
                    print->tp_varctrl[i].name, print->tp_varctrl[i].nnodes,
                    print->tp_varctrl[i].nvar, "POINT", 0.0000, "FETRIANGLE");

                for (j = 0; j < print->tp_varctrl[i].nnodes; j++)
                {
                    fprintf(print->tp_varctrl[i].datfile,
                        "%lf %lf %lf %lf %lf\n",
                        print->tp_varctrl[i].x[j], print->tp_varctrl[i].y[j],
                        print->tp_varctrl[i].zmin[j],
                        print->tp_varctrl[i].zmax[j],
                        0.000001);
                }
                for (j = 0; j < print->tp_varctrl[i].nvar; j++)
                {
                    fprintf(print->tp_varctrl[i].datfile, "%d %d %d\n",
                        print->tp_varctrl[i].node0[j],
                        print->tp_varctrl[i].node1[j],
                        print->tp_varctrl[i].node2[j]);
                }
            }
        }
    }
}

void UpdPrintVar(varctrl_struct *varctrl, int nprint, int module_step)
{
    int             i;
#ifdef _OPENMP
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

#ifdef _OPENMP
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
                    if (varctrl[i].counter > 0)
                    {
                        fprintf(varctrl[i].txtfile, "\t%lf",
                            varctrl[i].buffer[j] / (double)varctrl[i].counter);
                    }
                    else
                    {
                        fprintf(varctrl[i].txtfile, "\t%lf",
                            varctrl[i].buffer[j]);
                    }
                }
                fprintf(varctrl[i].txtfile, "\n");
                fflush(varctrl[i].txtfile);
            }

            outtime = (double)t;
            fwrite(&outtime, sizeof(double), 1, varctrl[i].datfile);
            for (j = 0; j < varctrl[i].nvar; j++)
            {
                if (varctrl[i].counter > 0)
                {
                    outval = varctrl[i].buffer[j] / (double)varctrl[i].counter;
                }
                else
                {
                    outval = varctrl[i].buffer[j];
                }
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

        init_file = fopen(fn, "wb");
        CheckFile(init_file, fn);

        for (i = 0; i < nelem; i++)
        {
            fwrite(&elem[i].ws.cmc, sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.sneqv, sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.surf, sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.unsat, sizeof(double), 1, init_file);
            fwrite(&elem[i].ws.gw, sizeof(double), 1, init_file);
#ifdef _NOAH_
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
                fwrite(&elem[i].ws.sh2o[j], sizeof(double), 1, init_file);
            }
#endif
        }

        for (i = 0; i < nriver; i++)
        {
            fwrite(&river[i].ws.stage, sizeof(double), 1, init_file);
            fwrite(&river[i].ws.gw, sizeof(double), 1, init_file);
        }

        fflush(init_file);
        fclose(init_file);
    }
}

void PrintDataTecplot(varctrl_struct *varctrl, int nprint, int t, int lapse)
{
    int             i;
    pihm_t_struct   pihm_time;
    double         *hnodes;    /* h at nodes */
    int            *inodes;

    pihm_time = PIHMTime(t);

    for (i = 0; i < nprint; i++)
    {
        int             j;
        double          outval;
        double          outtime;

        if(PrintNow(varctrl[i].intvl, lapse, &pihm_time))
        {
            outtime = (double)t;

            if (varctrl[i].intr == RIVERVAR)
            {
                /*Print river files */
                fprintf(varctrl[i].datfile, "%s\n", TEC_HEADER);
                fprintf(varctrl[i].datfile, "%s\n", RIVER_TEC_HEADER2);
                fprintf(varctrl[i].datfile, "%s %d\n", RIVER_TEC_HEADER3, t);
                for (j = 0; j < varctrl[i].nvar; j++)
                {
                    if (varctrl[i].counter > 0)
                    {
                        outval = varctrl[i].buffer[j] /
                            (double)varctrl[i].counter;
                    }
                    else
                    {
                        outval = varctrl[i].buffer[j];
                    }

                    fprintf(varctrl[i].datfile, "%lf %lf %lf %lf %lf\n",
                        varctrl[i].x[j], varctrl[i].y[j],
                        varctrl[i].zmin[j], varctrl[i].zmax[j], outval);
                    varctrl[i].buffer[j] = 0.0;
                }
            }
            else
            {
                /*Print element files */
                hnodes = (double *)calloc(varctrl[i].nnodes, sizeof(double));
                inodes = (int *)calloc(varctrl[i].nnodes, sizeof(int));
                for (j = 0; j < varctrl[i].nnodes; j++)
                {
                    hnodes[j] = 0.0;
                    inodes[j] = 0;
                }
                fprintf(varctrl[i].datfile,
                    "ZONE T=\"%s\", N=%d, E=%d, DATAPACKING=%s, "
                    "SOLUTIONTIME=%lf, ZONETYPE=%s\n",
                    varctrl[i].name, varctrl[i].nnodes, varctrl[i].nvar,
                    "POINT", outtime, "FETRIANGLE");
                fprintf(varctrl[i].datfile, "%s\n", ELEM_TEC_HEADER3);
                for (j = 0; j < varctrl[i].nvar; j++)
                {
                    if (varctrl[i].counter > 0)
                    {
                        outval = varctrl[i].buffer[j] /
                            (double)varctrl[i].counter;
                    }
                    else
                    {
                        outval = varctrl[i].buffer[j];
                    }

                    hnodes[varctrl[i].node0[j] - 1] += outval;
                    hnodes[varctrl[i].node1[j] - 1] += outval;
                    hnodes[varctrl[i].node2[j] - 1] += outval;
                    inodes[varctrl[i].node0[j] - 1] += 1;
                    inodes[varctrl[i].node1[j] - 1] += 1;
                    inodes[varctrl[i].node2[j] - 1] += 1;
                    varctrl[i].buffer[j] = 0.0;
                }
                for (j = 0; j < varctrl[i].nnodes; j++)
                {
                    if (inodes[j] == 0)
                    {
                        fprintf(varctrl[i].datfile, "%8.6f\n", 0.0);
                    }
                    else
                    {
                        fprintf(varctrl[i].datfile, "%8.6f\n",
                            hnodes[j] / inodes[j]);
                    }
                }

                free(hnodes);
                free(inodes);
            }

            varctrl[i].counter = 0;
            fflush(varctrl[i].datfile);
        }
    }
}

void PrintPerf(void *cvode_mem, int t, int starttime, double cputime_dt,
    double cputime, double maxstep, FILE *perf_file)
{
    static double   dt;
    static long int nst0, nfe0, nni0, ncfn0, netf0;
    long int        nst, nfe, nni, ncfn, netf;
    int             flag;

    flag = CVodeGetNumSteps(cvode_mem, &nst);
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    flag = CVodeGetNumErrTestFails(cvode_mem, &netf);

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
    int             i, j;
    double          totarea = 0.0;
    double          totlength = 0.0;
    double          tot_prep = 0.0;
    double          tot_net_prep = 0.0;
    double          tot_inf = 0.0;
    double          tot_rechg = 0.0;
    double          tot_esoil = 0.0;
    double          tot_ett = 0.0;
    double          tot_ec = 0.0;
    double          tot_pet = 0.0;
    double          tot_et = 0.0;
    double          tot_edir_surf = 0.0;
    double          tot_edir_unsat = 0.0;
    double          tot_edir_gw = 0.0;
    double          tot_ett_unsat = 0.0;
    double          tot_ett_gw = 0.0;
    double          outflow = 0.0;
    double          riv2elem_olf = 0.0;
    double          riv_base = 0.0;
    double          riv_lkg = 0.0;

    if (t == tstart + dt)
    {
        fprintf(watbal_file, "%s\n", WB_HEADER);
    }
    for (i = 0; i < nelem; i++)
    {
        totarea += elem[i].topo.area;
        tot_prep += elem[i].wf.prcp * elem[i].topo.area;
        tot_net_prep += elem[i].wf.pcpdrp * elem[i].topo.area;
        tot_inf += elem[i].wf.infil * elem[i].topo.area;
        tot_rechg += elem[i].wf.rechg * elem[i].topo.area;
        tot_esoil += elem[i].wf.edir * elem[i].topo.area;
        tot_ett += elem[i].wf.ett * elem[i].topo.area;
        tot_ec += elem[i].wf.ec * elem[i].topo.area;
        tot_pet += elem[i].wf.etp * elem[i].topo.area;
        tot_et += elem[i].wf.eta * elem[i].topo.area;
        tot_edir_surf += elem[i].wf.edir_surf * elem[i].topo.area;
        tot_edir_unsat += elem[i].wf.edir_unsat * elem[i].topo.area;
        tot_edir_gw += elem[i].wf.edir_gw * elem[i].topo.area;
        tot_ett_unsat += elem[i].wf.ett_unsat * elem[i].topo.area;
        tot_ett_gw += elem[i].wf.ett_gw * elem[i].topo.area;
    }

    for (j = 0; j < nriver; j++)
    {
        totlength += river[j].shp.length;
        if (river[j].down < 0)
        {
            outflow = river[j].wf.rivflow[DOWN_CHANL2CHANL];
        }
        riv2elem_olf +=
            river[j].wf.rivflow[LEFT_SURF2CHANL] +
            river[j].wf.rivflow[RIGHT_SURF2CHANL];
        riv_base +=
            river[j].wf.rivflow[LEFT_AQUIF2CHANL] +
            river[j].wf.rivflow[RIGHT_AQUIF2CHANL];
        riv_lkg += river[j].wf.rivflow[CHANL_LKG];
    }

    fprintf(watbal_file,
        "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
        "%lf %lf\n",
        t - tstart, outflow, riv2elem_olf, riv_base, riv_lkg, tot_prep,
        tot_net_prep, tot_inf, tot_rechg, tot_esoil, tot_ett, tot_ec, tot_pet,
        tot_et, tot_edir_surf, tot_edir_unsat, tot_edir_gw, tot_ett_unsat,
        tot_ett_gw);
}

void PrintCVodeFinalStats(void *cvode_mem)
{
    int             flag;
    long int        nst;
    long int        nfe;
    long int        netf;
    long int        nni;
    long int        ncfn;

    flag = CVodeGetNumSteps(cvode_mem, &nst);
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
    flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);

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

    return print;
}
