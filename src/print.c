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
    PIHMprintf (VL_NORMAL, "\n\t    The Penn State Integrated Hydrologic Model\n\n");

#ifdef _NOAH_
    PIHMprintf (VL_NORMAL, "\t    * Land surface module turned on.\n");
#endif
#ifdef _RT_
    PIHMprintf (VL_NORMAL, "\t       * Reactive transport module turned on.\n");
#endif
#ifdef _BGC_
    PIHMprintf (VL_NORMAL, "\t    * Biogeochemistry module turned on.\n");
#endif
#ifdef _ENKF_
    PIHMprintf (VL_NORMAL, "\t    * Ensemble Kalman filter turned on.\n");
#endif
#ifdef _CYCLES_
    PIHMprintf (VL_NORMAL, "\t    * Crop module turned on.\n");
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
            fprintf (stderr, " (%s, Line %d.\n)", fn, lineno);
        }
        fflush (stderr);
    }
    else if (verbosity <= verbose_mode)
    {
        vfprintf (stdout, fmt, va);
        if (debug_mode)
        {
            printf ("Printed from %s", func);
            printf (" (%s, Line %d.\n)", fn, lineno);
        }
        fflush (stderr);
    }

    va_end (va);
}

void MapOutput (char *simulation, pihm_struct pihm, char *outputdir)
{
    int             i, j, k;
    int             n, nT = 0;

    PIHMprintf (VL_VERBOSE, "\nInitializing PIHM output files\n");

    n = 0;

    for (i = 0; i < MAXPRINT; i++)
    {
        if (pihm->ctrl.prtvrbl[i] > 0)
        {
            switch (i)
            {
                case SURF_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.surf", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
					pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ws.surf;
                    }
                    n++;
                    break;
                case UNSAT_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.unsat", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ws.unsat;
                    }
                    n++;
                    break;
                case GW_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.gw", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ws.gw;
                    }
                    n++;
                    break;
                case RIVSTG_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.stage", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->riv[j].ws.stage;
                    }
                    n++;
                    break;
                case RIVGW_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivgw", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->riv[j].ws.gw;
                    }
                    n++;
                    break;
                case SNOW_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.snow", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ws.sneqv;
                    }
                    n++;
                    break;
                case CMC_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.is", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ws.cmc;
                    }
                    n++;
                    break;
                case INFIL_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.infil", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].wf.infil;
                    }
                    n++;
                    break;
                case RECHARGE_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.recharge",
                        outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].wf.rechg;
                    }
                    n++;
                    break;
                case EC_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.ec", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].wf.ec;
                    }
                    n++;
                    break;
                case ETT_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.ett", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].wf.ett;
                    }
                    n++;
                    break;
                case EDIR_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.edir", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].wf.edir;
                    }
                    n++;
                    break;
                case RIVFLX0_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx0", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].wf.rivflow[0];
                    }
                    n++;
                    break;
                case RIVFLX1_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx1", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].wf.rivflow[1];
                    }
                    n++;
                    break;
                case RIVFLX2_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx2", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].wf.rivflow[2];
                    }
                    n++;
                    break;
                case RIVFLX3_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx3", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].wf.rivflow[3];
                    }
                    n++;
                    break;
                case RIVFLX4_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx4", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].wf.rivflow[4];
                    }
                    n++;
                    break;
                case RIVFLX5_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx5", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].wf.rivflow[5];
                    }
                    n++;
                    break;
                case RIVFLX6_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx6", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].wf.rivflow[6];
                    }
                    n++;
                    break;
                case RIVFLX7_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx7", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].wf.rivflow[7];
                    }
                    n++;
                    break;
                case RIVFLX8_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx8", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].wf.rivflow[8];
                    }
                    n++;
                    break;
                case RIVFLX9_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx9", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].wf.rivflow[9];
                    }
                    n++;
                    break;
                case RIVFLX10_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx10",
                        outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].wf.rivflow[10];
                    }
                    n++;
                    break;
                case SUBFLX_CTRL:
                    for (k = 0; k < 3; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.subflx%d",
                            outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->elem[j].wf.subsurf[k];
                        }
                        n++;
                    }
                    break;
                case SURFFLX_CTRL:
                    for (k = 0; k < 3; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.surfflx%d",
                            outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->elem[j].wf.ovlflow[k];
                        }
                        n++;
                    }
                    break;
#ifdef _NOAH_
                case T1_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.t1", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].es.t1;
                    }
                    n++;
                    break;
                case STC_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.stc%d",
                            outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->elem[j].es.stc[k];
                        }
                        n++;
                    }
                    break;
                case SMC_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.smc%d",
                            outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->elem[j].ws.smc[k];
                        }
                        n++;
                    }
                    break;
                case SH2O_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.swc%d",
                            outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->elem[j].ws.sh2o[k];
                        }
                        n++;
                    }
                    break;
                case SNOWH_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.snowh", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ps.snowh;
                    }
                    n++;
                    break;
                case ALBEDO_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.albedo", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ps.albedo;
                    }
                    n++;
                    break;
                case LE_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.le", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ef.eta;
                    }
                    n++;
                    break;
                case SH_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.sh", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ef.sheat;
                    }
                    n++;
                    break;
                case G_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.g", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ef.ssoil;
                    }
                    n++;
                    break;
                case ETP_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.etp", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ef.etp;
                    }
                    n++;
                    break;
                case ESNOW_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.esnow", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ef.esnow;
                    }
                    n++;
                    break;
                case ROOTW_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rootw", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ps.soilw;
                    }
                    n++;
                    break;
                case SOILM_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.soilm", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ws.soilm;
                    }
                    n++;
                    break;
                case SOLAR_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.solar", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ef.soldn;
                    }
                    n++;
                    break;
#endif

#ifdef _BGC_
                case LAI_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.lai", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ps.proj_lai;
                    }
                    n++;
                    break;
                case VEGC_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.vegc", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].summary.vegc;
                    }
                    n++;
                    break;
                case LITRC_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.litrc", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].summary.litrc;
                    }
                    n++;
                    break;
                case SOILC_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.soilc", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].summary.soilc;
                    }
                    n++;
                    break;
                case TOTALC_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.totalc", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].summary.totalc;
                    }
                    n++;
                    break;
                case NPP_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.npp", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].summary.daily_npp;
                    }
                    n++;
                    break;
                case NEE_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.nee", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].summary.daily_nee;
                    }
                    n++;
                    break;
                case GPP_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.gpp", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].summary.daily_gpp;
                    }
                    n++;
                    break;
                case SMINN_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.sminn", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ns.sminn;
                    }
                    n++;
                    break;
#endif
#ifdef _CYCLES_
                case BIOMASS_CTRL:
                    for (k = 0; k < pihm->elem[0].comm.NumCrop; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.%s.biomass",
                            outputdir, simulation,
                            pihm->elem[0].comm.Crop[k].cropName);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->elem[j].comm.Crop[k].svBiomass;
                        }
                        n++;
                    }

                    sprintf (pihm->prtctrl[n].name, "%s%s.comm.biomass",
                        outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].comm.svBiomass;
                    }
                    n++;
                    break;
                case RADNINTCP_CTRL:
                    for (k = 0; k < pihm->elem[0].comm.NumCrop; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.%s.radintcp",
                            outputdir, simulation,
                            pihm->elem[0].comm.Crop[k].cropName);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->
                                elem[j].comm.Crop[k].svRadiationInterception;
                        }
                        n++;
                    }

                    sprintf (pihm->prtctrl[n].name, "%s%s.comm.radintcp",
                        outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].comm.svRadiationInterception;
                    }
                    n++;
                    break;
                case WATER_STS_CTRL:
                    for (k = 0; k < pihm->elem[0].comm.NumCrop; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.%s.waterstress",
                            outputdir, simulation,
                            pihm->elem[0].comm.Crop[k].cropName);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->
                                elem[j].comm.Crop[k].svWaterStressFactor;
                        }
                        n++;
                    }

                    sprintf (pihm->prtctrl[n].name, "%s%s.comm.waterstress",
                        outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].comm.svWaterStressFactor;
                    }
                    n++;
                    break;
                case N_STS_CTRL:
                    for (k = 0; k < pihm->elem[0].comm.NumCrop; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.%s.nstress",
                            outputdir, simulation,
                            pihm->elem[0].comm.Crop[k].cropName);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->elem[j].comm.Crop[k].svN_StressFactor;
                        }
                        n++;
                    }

                    sprintf (pihm->prtctrl[n].name, "%s%s.comm.nstress",
                        outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].comm.svN_StressFactor;
                    }
                    n++;
                    break;
                case CROP_TR_CTRL:
                    for (k = 0; k < pihm->elem[0].comm.NumCrop; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.%s.transp",
                            outputdir, simulation,
                            pihm->elem[0].comm.Crop[k].cropName);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->elem[j].comm.Crop[k].svTranspiration;
                        }
                        n++;
                    }

                    sprintf (pihm->prtctrl[n].name, "%s%s.comm.transp",
                        outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].comm.svTranspiration;
                    }
                    n++;
                    break;
                case CROP_POTTR_CTRL:
                    for (k = 0; k < pihm->elem[0].comm.NumCrop; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.%s.pottransp",
                            outputdir, simulation,
                            pihm->elem[0].comm.Crop[k].cropName);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->
                                elem[j].comm.Crop[k].svTranspirationPotential;
                        }
                        n++;
                    }

                    sprintf (pihm->prtctrl[n].name, "%s%s.comm.pottransp",
                        outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].comm.svTranspirationPotential;
                    }
                    n++;
                    break;
                case RES_EVAP_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.eres", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].wf.eres;
                    }
                    n++;
                    break;
                case NO3_PROF_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.NO3", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].soil.NO3Profile;
                    }
                    n++;
                    break;
                case NO3_RIVER_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivNO3", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].NO3sol.soluteMass;
                    }
                    n++;
                    break;
                case NH4_PROF_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.NH4", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].soil.NH4Profile;
                    }
                    n++;
                    break;
                case NH4_RIVER_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivNH4", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numriv;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->riv[j].NH4sol.soluteMass;
                    }
                    n++;
                    break;
                case NO3_DENIT_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.NO3denitrif",
                        outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] =
                            &pihm->elem[j].soil.NO3_Denitrification;
                    }
                    n++;
                    break;
                case NO3_LEACH_CTRL:
                    for (k = 0; k < NUM_EDGE; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.NO3leaching%d",
                            outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->elem[j].soil.NO3Leaching[k];
                        }
                        n++;
                    }
                    break;
                case NH4_LEACH_CTRL:
                    for (k = 0; k < NUM_EDGE; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.NH4leaching%d",
                            outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numele;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->elem[j].soil.NH4Leaching[k];
                        }
                        n++;
                    }
                    break;
                case NO3_LEACH_RIVER_CTRL:
                    for (k = 0; k < 4; k++)
                    {
                        sprintf (pihm->prtctrl[n].name,
                            "%s%s.riverNO3leaching%d",
                            outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numriv;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numriv; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->riv[j].NO3Leaching[k];
                        }
                        n++;
                    }
                    break;
                case NH4_LEACH_RIVER_CTRL:
                    for (k = 0; k < 4; k++)
                    {
                        sprintf (pihm->prtctrl[n].name,
                            "%s%s.riverNH4leaching%d",
                            outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvar = pihm->numriv;
                        pihm->prtctrl[n].var =
                            (double **)malloc (pihm->prtctrl[n].nvar *
                            sizeof (double *));
                        for (j = 0; j < pihm->numriv; j++)
                        {
                            pihm->prtctrl[n].var[j] =
                                &pihm->riv[j].NH4Leaching[k];
                        }
                        n++;
                    }
                    break;
                case LAI_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.lai", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvar = pihm->numele;
                    pihm->prtctrl[n].var =
                        (double **)malloc (pihm->prtctrl[n].nvar *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].var[j] = &pihm->elem[j].ps.proj_lai;
                    }
                    n++;
                    break;

#endif
				case SURFTEC_CTRL:
					sprintf(pihm->prtctrlT[nT].name, "%s%s_surf", outputdir,
						simulation);
					pihm->prtctrlT[nT].intvl = pihm->ctrl.prtvrbl[i];
					pihm->prtctrlT[nT].intr = 0;
					pihm->prtctrlT[nT].nvar = pihm->numele;
					pihm->prtctrlT[nT].nnodes = pihm->meshtbl.numnode;

					pihm->prtctrlT[nT].x =
						(double **)malloc(pihm->prtctrlT[nT].nnodes *
							sizeof(double *));
					pihm->prtctrlT[nT].y =
						(double **)malloc(pihm->prtctrlT[nT].nnodes *
							sizeof(double *));
					pihm->prtctrlT[nT].zmin =
						(double **)malloc(pihm->prtctrlT[nT].nnodes *
							sizeof(double *));
					pihm->prtctrlT[nT].zmax =
						(double **)malloc(pihm->prtctrlT[nT].nnodes *
							sizeof(double *));
					for (j = 0; j < pihm->prtctrlT[nT].nnodes; j++)
					{
						pihm->prtctrlT[nT].x[j] = &pihm->meshtbl.x[j];
						pihm->prtctrlT[nT].y[j] = &pihm->meshtbl.y[j];
						pihm->prtctrlT[nT].zmax[j] = &pihm->meshtbl.zmax[j];
						pihm->prtctrlT[nT].zmin[j] = &pihm->meshtbl.zmin[j];

					}
					pihm->prtctrlT[nT].node0 =
						(int **) malloc(pihm->prtctrlT[nT].nvar *
							sizeof(int *));
					pihm->prtctrlT[nT].node1 =
						(int **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(int *));
					pihm->prtctrlT[nT].node2 =
						(int **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(int *));
					pihm->prtctrlT[nT].var =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));

					for (j = 0; j < pihm->numele; j++)
					{
						pihm->prtctrlT[nT].var[j] = &pihm->elem[j].ws.surf;
						pihm->prtctrlT[nT].node0[j] = &pihm->elem[j].node[0];
						pihm->prtctrlT[nT].node1[j] = &pihm->elem[j].node[1];
						pihm->prtctrlT[nT].node2[j] = &pihm->elem[j].node[2];
					}
					nT++;
					break;
				case UNSATTEC_CTRL:
					sprintf(pihm->prtctrlT[nT].name, "%s%s_unsat", outputdir,
						simulation);
					pihm->prtctrlT[nT].intvl = pihm->ctrl.prtvrbl[i];
					pihm->prtctrlT[nT].intr = 0;
					pihm->prtctrlT[nT].nvar = pihm->numele;
					pihm->prtctrlT[nT].nnodes = pihm->meshtbl.numnode;

					pihm->prtctrlT[nT].x =
						(double **)malloc(pihm->prtctrlT[nT].nnodes *
							sizeof(double *));
					pihm->prtctrlT[nT].y =
						(double **)malloc(pihm->prtctrlT[nT].nnodes *
							sizeof(double *));
					pihm->prtctrlT[nT].zmin =
						(double **)malloc(pihm->prtctrlT[nT].nnodes *
							sizeof(double *));
					pihm->prtctrlT[nT].zmax =
						(double **)malloc(pihm->prtctrlT[nT].nnodes *
							sizeof(double *));
					for (j = 0; j < pihm->prtctrlT[nT].nnodes; j++)
					{
						pihm->prtctrlT[nT].x[j] = &pihm->meshtbl.x[j];
						pihm->prtctrlT[nT].y[j] = &pihm->meshtbl.y[j];
						pihm->prtctrlT[nT].zmax[j] = &pihm->meshtbl.zmax[j];
						pihm->prtctrlT[nT].zmin[j] = &pihm->meshtbl.zmin[j];

					}
					pihm->prtctrlT[nT].node0 =
						(int **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(int *));
					pihm->prtctrlT[nT].node1 =
						(int **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(int *));
					pihm->prtctrlT[nT].node2 =
						(int **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(int *));
					pihm->prtctrlT[nT].var =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));

					for (j = 0; j < pihm->numele; j++)
					{
						pihm->prtctrlT[nT].var[j] = &pihm->elem[j].ws.unsat;
						pihm->prtctrlT[nT].node0[j] = &pihm->elem[j].node[0];
						pihm->prtctrlT[nT].node1[j] = &pihm->elem[j].node[1];
						pihm->prtctrlT[nT].node2[j] = &pihm->elem[j].node[2];
					}
					nT++;
					break;
				case GWTEC_CTRL:
					sprintf(pihm->prtctrlT[nT].name, "%s%s_gw", outputdir,
						simulation);
					pihm->prtctrlT[nT].intvl = pihm->ctrl.prtvrbl[i];
					pihm->prtctrlT[nT].intr = 0;
					pihm->prtctrlT[nT].nvar = pihm->numele;
					pihm->prtctrlT[nT].nnodes = pihm->meshtbl.numnode;

					pihm->prtctrlT[nT].x =
						(double **)malloc(pihm->prtctrlT[nT].nnodes *
							sizeof(double *));
					pihm->prtctrlT[nT].y =
						(double **)malloc(pihm->prtctrlT[nT].nnodes *
							sizeof(double *));
					pihm->prtctrlT[nT].zmin =
						(double **)malloc(pihm->prtctrlT[nT].nnodes *
							sizeof(double *));
					pihm->prtctrlT[nT].zmax =
						(double **)malloc(pihm->prtctrlT[nT].nnodes *
							sizeof(double *));
					for (j = 0; j < pihm->prtctrlT[nT].nnodes; j++)
					{
						pihm->prtctrlT[nT].x[j] = &pihm->meshtbl.x[j];
						pihm->prtctrlT[nT].y[j] = &pihm->meshtbl.y[j];
						pihm->prtctrlT[nT].zmax[j] = &pihm->meshtbl.zmax[j];
						pihm->prtctrlT[nT].zmin[j] = &pihm->meshtbl.zmin[j];

					}
					pihm->prtctrlT[nT].node0 =
						(int **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(int *));
					pihm->prtctrlT[nT].node1 =
						(int **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(int *));
					pihm->prtctrlT[nT].node2 =
						(int **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(int *));
					pihm->prtctrlT[nT].var =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));

					for (j = 0; j < pihm->numele; j++)
					{
						pihm->prtctrlT[nT].var[j] = &pihm->elem[j].ws.gw;
						pihm->prtctrlT[nT].node0[j] = &pihm->elem[j].node[0];
						pihm->prtctrlT[nT].node1[j] = &pihm->elem[j].node[1];
						pihm->prtctrlT[nT].node2[j] = &pihm->elem[j].node[2];
					}
					nT++;
					break;
				case RIVSTGTEC_CTRL:
					sprintf(pihm->prtctrlT[nT].name, "%s%s_stage", outputdir,
						simulation);
					pihm->prtctrlT[nT].intvl = pihm->ctrl.prtvrbl[i];
					pihm->prtctrlT[nT].intr = 1;
					pihm->prtctrlT[nT].nvar = pihm->numriv;
					pihm->prtctrlT[nT].var =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));
					pihm->prtctrlT[nT].x =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));
					pihm->prtctrlT[nT].y =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));
					pihm->prtctrlT[nT].zmin =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));
					pihm->prtctrlT[nT].zmax =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));
					for (j = 0; j < pihm->numriv; j++)
					{
						pihm->prtctrlT[nT].var[j] = &pihm->riv[j].ws.stage;
						pihm->prtctrlT[nT].x[j] = &pihm->riv[j].topo.x;
						pihm->prtctrlT[nT].y[j] = &pihm->riv[j].topo.y;
						pihm->prtctrlT[nT].zmax[j] = &pihm->riv[j].topo.zmax;
						pihm->prtctrlT[nT].zmin[j] = &pihm->riv[j].topo.zmin;
					}
					nT++;
					break;
				case RIVGWTEC_CTRL:
					sprintf(pihm->prtctrlT[nT].name, "%s%s_rivgw", outputdir,
						simulation);
					pihm->prtctrlT[nT].intvl = pihm->ctrl.prtvrbl[i];
					pihm->prtctrlT[nT].intr = 1;
					pihm->prtctrlT[nT].nvar = pihm->numriv;
					pihm->prtctrlT[nT].var =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));
					pihm->prtctrlT[nT].x =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));
					pihm->prtctrlT[nT].y =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));
					pihm->prtctrlT[nT].zmin =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));
					pihm->prtctrlT[nT].zmax =
						(double **)malloc(pihm->prtctrlT[nT].nvar *
							sizeof(double *));
					for (j = 0; j < pihm->numriv; j++)
					{
						pihm->prtctrlT[nT].var[j] = &pihm->riv[j].ws.gw;
						pihm->prtctrlT[nT].x[j] = &pihm->riv[j].topo.x;
						pihm->prtctrlT[nT].y[j] = &pihm->riv[j].topo.y;
						pihm->prtctrlT[nT].zmax[j] = &pihm->riv[j].topo.zmax;
						pihm->prtctrlT[nT].zmin[j] = &pihm->riv[j].topo.zmin;
					}
					nT++;
					break;

                default:
                    break;
            }
        }
    }

    if (n > MAXPRINT)
    {
        PIHMprintf (VL_ERROR, "Error: Too many output files. ");
        PIHMprintf (VL_ERROR, "The maximum number of output files is %d.\n",
            MAXPRINT);
        PIHMexit (EXIT_FAILURE);
    }

    pihm->ctrl.nprint = n;
	pihm->ctrl.nprintT = nT;

    for (i = 0; i < pihm->ctrl.nprint; i++)
    {
        pihm->prtctrl[i].buffer =
            (double *)calloc (pihm->prtctrl[i].nvar, sizeof (double));
    }
	for (i = 0; i < pihm->ctrl.nprintT; i++)
	{
		pihm->prtctrlT[i].buffer =
			(double *)calloc(pihm->prtctrlT[i].nvar, sizeof(double));
	}
}

void InitOutputFile (prtctrl_struct *prtctrl, int nprint, int ascii)
{
    FILE           *fid;
    char            ascii_fn[MAXSTRING];
    char            dat_fn[MAXSTRING];
    int             i;

    for (i = 0; i < nprint; i++)
    {
        sprintf (dat_fn, "%s.dat", prtctrl[i].name);
        fid = fopen (dat_fn, "w");
        fclose (fid);

        if (ascii)
        {
            sprintf (ascii_fn, "%s.txt", prtctrl[i].name);
            fid = fopen (ascii_fn, "w");
            fclose (fid);
        }
    }
}

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

        for (j = 0; j < prtctrl[i].nvar; j++)
        {
            prtctrl[i].buffer[j] += *prtctrl[i].var[j];
/* 		if (j = 55) {
			printf("test");
		}*/       
		}

        if (lapse % prtctrl[i].intvl == 0 && lapse > 0)
        {
            rawtime = t;
            timestamp = gmtime (&rawtime);

            if (ascii)
            {
                sprintf (ascii_fn, "%s.txt", prtctrl[i].name);
                fid = fopen (ascii_fn, "a");
                CheckFile (fid, ascii_fn);

                fprintf (fid, "\"%4.4d-%2.2d-%2.2d %2.2d:%2.2d\"",
                    timestamp->tm_year + 1900, timestamp->tm_mon + 1,
                    timestamp->tm_mday, timestamp->tm_hour,
                    timestamp->tm_min);
                for (j = 0; j < prtctrl[i].nvar; j++)
                {
                    if (prtctrl[i].intvl > dt)
                    {
                        fprintf (fid, "\t%lf",
                            prtctrl[i].buffer[j] /
                            ((double)(prtctrl[i].intvl / dt)));
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
            CheckFile (fid, dat_fn);

            outtime = (double)t;
            fwrite (&outtime, sizeof (double), 1, fid);
            for (j = 0; j < prtctrl[i].nvar; j++)
            {
                if (prtctrl[i].intvl > dt)
                {
                    outval =
                        prtctrl[i].buffer[j] / ((double)(prtctrl[i].intvl /
                            dt));
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

void PrintDataTecplot(prtctrlT_struct *prtctrlT, int nprintT, int t, int lapse, int dt, int count)
{
	int             i, j;
	char            dat_fn[MAXSTRING];
	FILE           *fid;
	double          outval;
	double          outtime;
	const char		*str1_Tec, *str2_Tec, *str3_Tec;
	realtype        *hnodes;   /* h at nodes */
	int				*inodes;


	outtime = (double)t;	
    str1_Tec = "VARIABLES = \"X\" \"Y\" \"Zmin\" \"Zmax\" \"h\"";
	for (i = 0; i < nprintT; i++)
	{
		
		for (j = 0; j < prtctrlT[i].nvar; j++)
		{
			prtctrlT[i].buffer[j] += *prtctrlT[i].var[j];
		}

		  sprintf(dat_fn, "%s.plt", prtctrlT[i].name);
		  fid = fopen(dat_fn, "a");
		  CheckFile(fid, dat_fn);	
		if (lapse % prtctrlT[i].intvl == 0 && lapse > 0)
		{
	
  		  if (prtctrlT[i].intr == 1)
			{
				/*Print river files */
				str2_Tec = "ZONE T = \"Water Depth River\" ";
				str3_Tec = "StrandID=1, SolutionTime=";

				fprintf(fid, "%s \n", str1_Tec);
				fprintf(fid, "%s \n", str2_Tec);
				fprintf(fid, "%s %d \n", str3_Tec, t);
				for (j = 0; j < prtctrlT[i].nvar; j++)
				{
					if (prtctrlT[i].intvl > dt)
					{
						outval =
							prtctrlT[i].buffer[j] / ((double)(prtctrlT[i].intvl /
								dt));
					}
					else
					{
						outval = prtctrlT[i].buffer[j];
					}

					fprintf(fid, "%lf %lf %lf %lf %lf \n", *prtctrlT[i].x[j], *prtctrlT[i].y[j], *prtctrlT[i].zmin[j], *prtctrlT[i].zmax[j], outval);
					prtctrlT[i].buffer[j] = 0.0;
				}
			}
			else
			{

				/*Print element files */
				hnodes = (double *)calloc(prtctrlT[i].nnodes, sizeof(double));
				inodes = (int *)calloc(prtctrlT[i].nnodes, sizeof(int));
				for (j = 0; j < prtctrlT[i].nnodes; j++)
				{
					hnodes[j] = 0.0;
					inodes[j] = 0;
				}
				if ( (int)(count/dt) <= 0)
				{
					fprintf(fid, "%s \n", str1_Tec);
					fprintf(fid, "%s %s %s %d %s %d %s %lf %s\n", "ZONE T=\"", prtctrlT[i].name, "\", N=", prtctrlT[i].nnodes, ", E=", prtctrlT[i].nvar, "DATAPACKING=POINT, SOLUTIONTIME = ", 0.0000, ", ZONETYPE=FETRIANGLE");

					for (j = 0; j < prtctrlT[i].nnodes; j++)
					{
						fprintf(fid, "%lf %lf %lf %lf %lf\n", *prtctrlT[i].x[j], *prtctrlT[i].y[j], *prtctrlT[i].zmin[j], *prtctrlT[i].zmax[j], 0.000001);
					}
					for (j = 0; j < prtctrlT[i].nvar; j++)
					{
						fprintf(fid, "%d %d %d \n", *prtctrlT[i].node0[j], *prtctrlT[i].node1[j], *prtctrlT[i].node2[j]);
					}
				}
				str3_Tec = "VARSHARELIST = ([1, 2, 3, 4]=1), CONNECTIVITYSHAREZONE = 1";
				fprintf(fid, "%s %s %s %d %s %d %s %lf %s\n", "ZONE T=\"", prtctrlT[i].name, "\", N=", prtctrlT[i].nnodes, ", E=", prtctrlT[i].nvar, "DATAPACKING=POINT, SOLUTIONTIME = ", outtime, ", ZONETYPE=FETRIANGLE,");
				fprintf(fid, "%s \n", str3_Tec);
				for (j = 0; j < prtctrlT[i].nvar; j++)
				{
					if (prtctrlT[i].intvl > dt)
					{
						outval =
							prtctrlT[i].buffer[j] / ((double)(prtctrlT[i].intvl /
								dt));
					}
					else
					{
						outval = prtctrlT[i].buffer[j];
					}

					hnodes[*prtctrlT[i].node0[j] - 1] = hnodes[*prtctrlT[i].node0[j] - 1] + outval;
					hnodes[*prtctrlT[i].node1[j] - 1] = hnodes[*prtctrlT[i].node1[j] - 1] + outval;
					hnodes[*prtctrlT[i].node2[j] - 1] = hnodes[*prtctrlT[i].node2[j] - 1] + outval;
					inodes[*prtctrlT[i].node0[j] - 1] = inodes[*prtctrlT[i].node0[j] - 1] + 1;
					inodes[*prtctrlT[i].node1[j] - 1] = inodes[*prtctrlT[i].node1[j] - 1] + 1;
					inodes[*prtctrlT[i].node2[j] - 1] = inodes[*prtctrlT[i].node2[j] - 1] + 1;
					prtctrlT[i].buffer[j] = 0.0;
				}
				for (j = 0; j < prtctrlT[i].nnodes; j++) 
				{
					if (inodes[j] == 0) {
						fprintf(fid, "%8.6f \n", 0.0);
					}
					else {
						fprintf(fid, "%8.6f \n", hnodes[j] / inodes[j]);
					}			
				}
			}
			
		}
		fclose(fid);
	}

}
void PrtInit (pihm_struct pihm, char *simulation, int t)
{
    FILE           *init_file;
    char            fn[MAXSTRING];
    int             i;
	struct tm      *timestamp;
	time_t          rawtime;
#ifdef _NOAH_
    int             j;
#endif

	rawtime = t;
	timestamp = gmtime(&rawtime);



    //sprintf (fn, "input/%s/%s%d.ic", project, simulation, t);
	sprintf(fn, "input/%s/%s%4.4d-%2.2d-%2.2d %2.2d:%2.2d.ic", project, simulation, timestamp->tm_year + 1900, timestamp->tm_mon + 1,
		timestamp->tm_mday, timestamp->tm_hour,
		timestamp->tm_min);

    init_file = fopen (fn, "wb");

    for (i = 0; i < pihm->numele; i++)
    {
        fwrite (&pihm->elem[i].ws.cmc, sizeof (double), 1, init_file);
        fwrite (&pihm->elem[i].ws.sneqv, sizeof (double), 1, init_file);
        fwrite (&pihm->elem[i].ws.surf, sizeof (double), 1, init_file);
        fwrite (&pihm->elem[i].ws.unsat, sizeof (double), 1, init_file);
        fwrite (&pihm->elem[i].ws.gw, sizeof (double), 1, init_file);
#ifdef _NOAH_
        fwrite (&pihm->elem[i].es.t1, sizeof (double), 1, init_file);
        fwrite (&pihm->elem[i].ps.snowh, sizeof (double), 1, init_file);
        for (j = 0; j < MAXLYR; j++)
        {
            fwrite (&pihm->elem[i].es.stc[j], sizeof (double), 1, init_file);
        }
        for (j = 0; j < MAXLYR; j++)
        {
            fwrite (&pihm->elem[i].ws.smc[j], sizeof (double), 1, init_file);
        }
        for (j = 0; j < MAXLYR; j++)
        {
            fwrite (&pihm->elem[i].ws.sh2o[j], sizeof (double), 1, init_file);
        }
#endif
    }

    for (i = 0; i < pihm->numriv; i++)
    {
        fwrite (&pihm->riv[i].ws.stage, sizeof (double), 1, init_file);
        fwrite (&pihm->riv[i].ws.gw, sizeof (double), 1, init_file);
    }

    fclose (init_file);
}

void PrintStats(void *cvode_mem)
{
	long int nst, nfe, nfeLS, nni, ncfn, netf;
	int flag;
	FILE *Converg; /* Convergence file */
	Converg = fopen("CVODE.log", "a");

	flag = CVodeGetNumSteps(cvode_mem, &nst);
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
	flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
	flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
	flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
	
	//fprintf(Converg, "\nCVODE Statistics:\n");
	fprintf(Converg, "nst = %-6ld nni = %-6ld nfe = %-6ld netf = %-6ld ncfn = %-6ld nfeLS = %-6ld\n",
		nst, nni, nfe, netf, ncfn, nfeLS);
	fclose(Converg);
}

void PrintWaterBalance(FILE *WaterBalance, int dt, elem_struct *elem, int numele, river_struct *riv, int numriv, int count)
{
	long int i,j;
	const char		*str1_Tec;
	realtype        totarea = 0., totlenght = 0.;
	realtype        totPrep = 0., totNetPrep = 0., totInf = 0., totRecharge = 0., totEsoil = 0., \
		            totETplant = 0., totEcan = 0., totPET = 0., totET = 0., totES = 0., totEU = 0., \
		            totEGW = 0., totTU = 0., totTGW = 0.;
	realtype        outflow, RE_OLF = 0., R_Exf = 0.;


	str1_Tec = "VARIABLES = \"TIME (s)\" \"Outflow (cms)\" \"Surf2Chan (cms)\" \"AqF2Chan (cms)\" \
          \"Precipitation (cms)\" \"NetPrec (cms)\" \"Infiltration (cms)\" \"Recharge (cms)\" \
          \"E_soil (cms)\" \"ET_plant (cms)\" \"E_canopy (cms)\" \"PET (cms)\" \"ET (cms)\" \
          \"E_Surface (cms)\" \"E_Unsat (cms)\" \"E_GW (cms)\" \"T_Unsat (cms)\" \"T_GW (cms)\"";
	
	if (count == 0) {
		fprintf(WaterBalance, "%s\n", str1_Tec);
	}

	for (i = 0; i < numele; i++) {
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

	for (j = 0; j < numriv; j++) {
		totlenght = totlenght + riv[j].shp.length;
		if (riv[j].down < 0) {
			outflow = riv[j].wf.rivflow[DOWN_CHANL2CHANL];
		}
		RE_OLF = RE_OLF + (riv[j].wf.rivflow[LEFT_SURF2CHANL] + riv[j].wf.rivflow[RIGHT_SURF2CHANL]);
		R_Exf = R_Exf = +(riv[j].wf.rivflow[LEFT_AQUIF2CHANL] + riv[j].wf.rivflow[RIGHT_AQUIF2CHANL]);

	}

	fprintf(WaterBalance, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
		dt, outflow, RE_OLF, R_Exf, totPrep, totNetPrep, totInf, totRecharge, totEsoil, totETplant, totEcan, totPET, totET, totES, totEU, totEGW, totTU, totTGW);
}
