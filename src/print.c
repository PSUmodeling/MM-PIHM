#include "pihm.h"

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
    PIHMprintf(VL_NORMAL,
        "\n\tThe Penn State Integrated Hydrologic Model\n\n");

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
            fprintf(stderr, " (%s, Line %d.)\n", fn, lineno);
        }
        fflush(stderr);
    }
    else if (verbosity <= verbose_mode)
    {
        vfprintf(stdout, fmt, va);
        if (debug_mode)
        {
            printf("Printed from %s", func);
            printf(" (%s, Line %d.)\n", fn, lineno);
        }
        fflush(stderr);
    }

    va_end(va);
}

void InitOutputFile(print_struct *print, int ascii)
{
    char            ascii_fn[MAXSTRING];
    char            dat_fn[MAXSTRING];
    char            tec_fn[MAXSTRING];
    int             i;

    for (i = 0; i < print->nprint; i++)
    {
        sprintf(dat_fn, "%s.dat", print->varctrl[i].name);
        print->varctrl[i].datfile = fopen(dat_fn, "w");

        if (ascii)
        {
            sprintf(ascii_fn, "%s.txt", print->varctrl[i].name);
            print->varctrl[i].txtfile = fopen(ascii_fn, "w");
        }
    }
    if (tecplot)
    {
        for (i = 0; i < print->ntpprint; i++)
        {
            sprintf(tec_fn, "%s.plt", print->tp_varctrl[i].name);
            print->tp_varctrl[i].datfile = fopen(tec_fn, "w");
        }
    }
}

void UpdPrintVar(varctrl_struct *prtctrl, int nprint, int module_step)
{
    int             i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nprint; i++)
    {
        int             j;

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

void UpdPrintVarT(varctrl_struct *prtctrlT, int nprintT)
{
    int             i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nprintT; i++)
    {
        int             j;
        for (j = 0; j < prtctrlT[i].nvar; j++)
        {
            prtctrlT[i].buffer[j] += *prtctrlT[i].var[j];
        }

        prtctrlT[i].counter++;
    }
}

void PrintData(varctrl_struct *prtctrl, int nprint, int t, int lapse, int ascii)
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
        int             print = 0;
        double          outval;
        double          outtime;

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
                fprintf(prtctrl[i].txtfile, "\"%s\"", pihm_time.str);
                for (j = 0; j < prtctrl[i].nvar; j++)
                {
                    if (prtctrl[i].counter > 0)
                    {
                        fprintf(prtctrl[i].txtfile, "\t%lf",
                            prtctrl[i].buffer[j] / (double)prtctrl[i].counter);
                    }
                    else
                    {
                        fprintf(prtctrl[i].txtfile, "\t%lf",
                            prtctrl[i].buffer[j]);
                    }
                }
                fprintf(prtctrl[i].txtfile, "\n");
                fflush(prtctrl[i].txtfile);
            }

            outtime = (double)t;
            fwrite(&outtime, sizeof(double), 1, prtctrl[i].datfile);
            for (j = 0; j < prtctrl[i].nvar; j++)
            {
                if (prtctrl[i].counter > 0)
                {
                    outval = prtctrl[i].buffer[j] / (double)prtctrl[i].counter;
                }
                else
                {
                    outval = prtctrl[i].buffer[j];
                }
                fwrite(&outval, sizeof(double), 1, prtctrl[i].datfile);

                prtctrl[i].buffer[j] = 0.0;
            }
            prtctrl[i].counter = 0;
            fflush(prtctrl[i].datfile);
        }
    }
}

void PrtInit(elem_struct *elem, river_struct *river, char *simulation, int t)
{
    FILE           *init_file;
    char            fn[MAXSTRING];
    int             i;
    char            name[20];

    pihm_t_struct   pihm_time;

    pihm_time = PIHMTime(t);
    strcpy(name, pihm_time.str);
    name[13] = 0;

#ifdef _NOAH_
    int             j;
#endif

    sprintf(fn, "input/%s/ic/%s %s.ic", project, simulation, name);

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

    fclose(init_file);
}

void PrintDataTecplot(varctrl_struct *varctrl, int nprint, int t, int lapse)
{
    int             i;
    pihm_t_struct   pihm_time;
    const char     *str1_Tec, *str2_Tec, *str3_Tec;
    pihm_time = PIHMTime(t);
    realtype       *hnodes;    /* h at nodes */
    int            *inodes;

    str1_Tec = "VARIABLES = \"X\" \"Y\" \"Zmin\" \"Zmax\" \"h\"";
    for (i = 0; i < nprint; i++)
    {
        int             j;
        int             print = 0;
        double          outval;
        double          outtime;

        switch (varctrl[i].intvl)
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
                if (lapse % varctrl[i].intvl == 0 && lapse > 0)
                {
                    print = 1;
                }
        }
        if (print)
        {
            outtime = (double)t;

            if (varctrl[i].intr == 1)
            {
                /*Print river files */
                str2_Tec = "ZONE T = \"Water Depth River\" ";
                str3_Tec = "StrandID=1, SolutionTime=";

                fprintf(varctrl[i].datfile, "%s \n", str1_Tec);
                fprintf(varctrl[i].datfile, "%s \n", str2_Tec);
                fprintf(varctrl[i].datfile, "%s %d \n", str3_Tec, t);
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

                    fprintf(varctrl[i].datfile, "%lf %lf %lf %lf %lf \n",
                        *varctrl[i].x[j], *varctrl[i].y[j],
                        *varctrl[i].zmin[j], *varctrl[i].zmax[j], outval);
                    varctrl[i].buffer[j] = 0.0;
                }
                varctrl[i].counter = 0;
                fflush(varctrl[i].datfile);
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
                if (varctrl[i].first)
                {
                    fprintf(varctrl[i].datfile, "%s \n", str1_Tec);
                    fprintf(varctrl[i].datfile,
                        "%s %s %s %d %s %d %s %lf %s\n", "ZONE T=\"",
                        varctrl[i].name, "\", N=", varctrl[i].nnodes, ", E=",
                        varctrl[i].nvar, "DATAPACKING=POINT, SOLUTIONTIME = ",
                        0.0000, ", ZONETYPE=FETRIANGLE");

                    for (j = 0; j < varctrl[i].nnodes; j++)
                    {
                        fprintf(varctrl[i].datfile, "%lf %lf %lf %lf %lf\n",
                            *varctrl[i].x[j], *varctrl[i].y[j],
                            *varctrl[i].zmin[j], *varctrl[i].zmax[j],
                            0.000001);
                    }
                    for (j = 0; j < varctrl[i].nvar; j++)
                    {
                        fprintf(varctrl[i].datfile, "%d %d %d \n",
                            *varctrl[i].node0[j], *varctrl[i].node1[j],
                            *varctrl[i].node2[j]);
                    }
                    varctrl[i].first = 0;
                }
                str3_Tec =
                    "VARSHARELIST = ([1, 2, 3, 4]=1), "
                    "CONNECTIVITYSHAREZONE = 1";
                fprintf(varctrl[i].datfile, "%s %s %s %d %s %d %s %lf %s\n",
                    "ZONE T=\"", varctrl[i].name, "\", N=", varctrl[i].nnodes,
                    ", E=", varctrl[i].nvar,
                    "DATAPACKING=POINT, SOLUTIONTIME = ", outtime,
                    ", ZONETYPE=FETRIANGLE,");
                fprintf(varctrl[i].datfile, "%s \n", str3_Tec);
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

                    hnodes[*varctrl[i].node0[j] - 1] =
                        hnodes[*varctrl[i].node0[j] - 1] + outval;
                    hnodes[*varctrl[i].node1[j] - 1] =
                        hnodes[*varctrl[i].node1[j] - 1] + outval;
                    hnodes[*varctrl[i].node2[j] - 1] =
                        hnodes[*varctrl[i].node2[j] - 1] + outval;
                    inodes[*varctrl[i].node0[j] - 1] =
                        inodes[*varctrl[i].node0[j] - 1] + 1;
                    inodes[*varctrl[i].node1[j] - 1] =
                        inodes[*varctrl[i].node1[j] - 1] + 1;
                    inodes[*varctrl[i].node2[j] - 1] =
                        inodes[*varctrl[i].node2[j] - 1] + 1;
                    varctrl[i].buffer[j] = 0.0;
                }
                for (j = 0; j < varctrl[i].nnodes; j++)
                {
                    if (inodes[j] == 0)
                    {
                        fprintf(varctrl[i].datfile, "%8.6f \n", 0.0);
                    }
                    else
                    {
                        fprintf(varctrl[i].datfile, "%8.6f \n",
                            hnodes[j] / inodes[j]);
                    }
                }
                varctrl[i].counter = 0;
                fflush(varctrl[i].datfile);
            }
        }
    }
}

void PrintStats(void *cvode_mem, FILE *Conv)
{
    long int        nst, nfe, nni, ncfn, netf;
    int             flag;

    flag = CVodeGetNumSteps(cvode_mem, &nst);
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    flag = CVodeGetNumErrTestFails(cvode_mem, &netf);

    fprintf(Conv,
        "nst = %-6ld nni = %-6ld nfe = %-6ld netf = %-6ld ncfn = %-6ld\n", nst,
        nni, nfe, netf, ncfn);
}


void PrintWaterBalance(FILE *WaterBalance, int t, int tstart, int dt,
    elem_struct *elem, int numele, river_struct *riv, int numriv)
{
    long int        i, j;
    const char     *str1_Tec;
    realtype        totarea = 0., totlenght = 0.;
    realtype        totPrep = 0., totNetPrep = 0., totInf = 0., totRecharge =
        0., totEsoil = 0., totETplant = 0., totEcan = 0., totPET = 0., totET =
        0., totES = 0., totEU = 0., totEGW = 0., totTU = 0., totTGW = 0.;
    realtype        outflow = 0.0, RE_OLF = 0., R_Exf = 0., R_LKG = 0.;


    str1_Tec =
        "VARIABLES = \"TIME (s)\" \"Outflow (cms)\" \"Surf2Chan (cms)\""
        "\"AqF2Chan (cms)\" \"Chan_LKG (cms)\" \"Precipitation (cms)\""
        "\"NetPrec (cms)\" \"Infiltration (cms)\" \"Recharge (cms)\""
        "\"E_soil (cms)\" \"ET_plant (cms)\" \"E_canopy (cms)\" \"PET (cms)\""
        "\"ET (cms)\" \"E_Surface (cms)\" \"E_Unsat (cms)\" \"E_GW (cms)\""
        "\"T_Unsat (cms)\" \"T_GW (cms)\"";

    if (t == tstart + dt)
    {
        fprintf(WaterBalance, "%s\n", str1_Tec);
    }
    for (i = 0; i < numele; i++)
    {
        totarea = totarea + elem[i].topo.area;
        totPrep = totPrep + elem[i].wf.prcp * elem[i].topo.area;
        totNetPrep = totNetPrep + elem[i].wf.pcpdrp * elem[i].topo.area;
        totInf = totInf + elem[i].wf.infil * elem[i].topo.area;
        totRecharge = totRecharge + elem[i].wf.rechg * elem[i].topo.area;
        totEsoil = totEsoil + elem[i].wf.edir * elem[i].topo.area;
        totETplant = totETplant + elem[i].wf.ett * elem[i].topo.area;
        totEcan = totEcan + elem[i].wf.ec * elem[i].topo.area;
        totPET = totPET + elem[i].wf.etp * elem[i].topo.area;
        totET = totET + elem[i].wf.eta * elem[i].topo.area;
        totES = totES + elem[i].wf.edir_surf * elem[i].topo.area;
        totEU = totEU + elem[i].wf.edir_unsat * elem[i].topo.area;
        totEGW = totEGW + elem[i].wf.edir_gw * elem[i].topo.area;
        totTU = totTU + elem[i].wf.ett_unsat * elem[i].topo.area;
        totTGW = totTGW + elem[i].wf.ett_gw * elem[i].topo.area;
    }

    for (j = 0; j < numriv; j++)
    {
        totlenght = totlenght + riv[j].shp.length;
        if (riv[j].down < 0)
        {
            outflow = riv[j].wf.rivflow[DOWN_CHANL2CHANL];
        }
        RE_OLF =
            RE_OLF + (riv[j].wf.rivflow[LEFT_SURF2CHANL] +
            riv[j].wf.rivflow[RIGHT_SURF2CHANL]);
        R_Exf =
            R_Exf + (riv[j].wf.rivflow[LEFT_AQUIF2CHANL] +
            riv[j].wf.rivflow[RIGHT_AQUIF2CHANL]);
        R_LKG = R_LKG + riv[j].wf.rivflow[CHANL_LKG];

    }

    fprintf(WaterBalance,
        "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
        "%lf %lf\n",
        (t - tstart), outflow, RE_OLF, R_Exf, R_LKG, totPrep, totNetPrep,
        totInf, totRecharge, totEsoil, totETplant, totEcan, totPET, totET,
        totES, totEU, totEGW, totTU, totTGW);
}
