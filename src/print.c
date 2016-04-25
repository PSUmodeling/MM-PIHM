#include "pihm.h"

void MapOutput (char *simulation, pihm_struct pihm, char *outputdir)
{
    int             i, j, k;
    int             n;

    if (verbose_mode)
        printf ("\nInitializing PIHM output files\n");

    n = 0;

    for (i = 0; i < NUM_PRINT; i++)
    {
        if (pihm->ctrl.prtvrbl[i] > 0)
        {
            switch (i)
            {
                case SURF_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.surf", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ws.surf;
                    }
                    n++;
                    break;
                case UNSAT_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.unsat", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ws.unsat;
                    }
                    n++;
                    break;
                case GW_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.gw", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ws.gw;
                    }
                    n++;
                    break;
                case RIVSTG_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.stage", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].ws.stage;
                    }
                    n++;
                    break;
                case RIVGW_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivgw", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->riv[j].ws.gw;
                    }
                    n++;
                    break;
                case SNOW_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.snow", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ws.sneqv;
                    }
                    n++;
                    break;
                case CMC_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.is", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ws.cmc;
                    }
                    n++;
                    break;
                case INFIL_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.infil", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].wf.infil;
                    }
                    n++;
                    break;
                case RECHARGE_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.recharge",
                        outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].wf.rechg;
                    }
                    n++;
                    break;
                case EC_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.ec", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].wf.ec;
                    }
                    n++;
                    break;
                case ETT_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.ett", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].wf.ett;
                    }
                    n++;
                    break;
                case EDIR_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.edir", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].wf.edir;
                    }
                    n++;
                    break;
                case RIVFLX0_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx0", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] =
                            &pihm->riv[j].wf.fluxriv[0];
                    }
                    n++;
                    break;
                case RIVFLX1_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx1", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] =
                            &pihm->riv[j].wf.fluxriv[1];
                    }
                    n++;
                    break;
                case RIVFLX2_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx2", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] =
                            &pihm->riv[j].wf.fluxriv[2];
                    }
                    n++;
                    break;
                case RIVFLX3_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx3", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] =
                            &pihm->riv[j].wf.fluxriv[3];
                    }
                    n++;
                    break;
                case RIVFLX4_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx4", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] =
                            &pihm->riv[j].wf.fluxriv[4];
                    }
                    n++;
                    break;
                case RIVFLX5_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx5", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] =
                            &pihm->riv[j].wf.fluxriv[5];
                    }
                    n++;
                    break;
                case RIVFLX6_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx6", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] =
                            &pihm->riv[j].wf.fluxriv[6];
                    }
                    n++;
                    break;
                case RIVFLX7_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx7", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] =
                            &pihm->riv[j].wf.fluxriv[7];
                    }
                    n++;
                    break;
                case RIVFLX8_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx8", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] =
                            &pihm->riv[j].wf.fluxriv[8];
                    }
                    n++;
                    break;
                case RIVFLX9_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx9", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] =
                            &pihm->riv[j].wf.fluxriv[9];
                    }
                    n++;
                    break;
                case RIVFLX10_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rivflx10",
                        outputdir, simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numriv;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numriv; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] =
                            &pihm->riv[j].wf.fluxriv[10];
                    }
                    n++;
                    break;
                case SUBFLX_CTRL:
                    for (k = 0; k < 3; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.subflx%d",
                            outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvrbl = pihm->numele;
                        pihm->prtctrl[n].vrbl =
                            (double **)malloc (pihm->prtctrl[n].nvrbl *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].vrbl[j] =
                                &pihm->elem[j].wf.fluxsub[k];
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
                        pihm->prtctrl[n].nvrbl = pihm->numele;
                        pihm->prtctrl[n].vrbl =
                            (double **)malloc (pihm->prtctrl[n].nvrbl *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].vrbl[j] =
                                &pihm->elem[j].wf.fluxsurf[k];
                        }
                        n++;
                    }
                    break;
#ifdef _NOAH_
                case T1_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.t1", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].es.t1;
                    }
                    n++;
                    break;
                case STC_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf (pihm->prtctrl[n].name, "%s%s.stc%d",
                            outputdir, simulation, k);
                        pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                        pihm->prtctrl[n].nvrbl = pihm->numele;
                        pihm->prtctrl[n].vrbl =
                            (double **)malloc (pihm->prtctrl[n].nvrbl *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].vrbl[j] =
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
                        pihm->prtctrl[n].nvrbl = pihm->numele;
                        pihm->prtctrl[n].vrbl =
                            (double **)malloc (pihm->prtctrl[n].nvrbl *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].vrbl[j] =
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
                        pihm->prtctrl[n].nvrbl = pihm->numele;
                        pihm->prtctrl[n].vrbl =
                            (double **)malloc (pihm->prtctrl[n].nvrbl *
                            sizeof (double *));
                        for (j = 0; j < pihm->numele; j++)
                        {
                            pihm->prtctrl[n].vrbl[j] =
                                &pihm->elem[j].ws.sh2o[k];
                        }
                        n++;
                    }
                    break;
                case SNOWH_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.snowh", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ps.snowh;
                    }
                    n++;
                    break;
                case ALBEDO_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.albedo", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ps.albedo;
                    }
                    n++;
                    break;
                case LE_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.le", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ef.eta;
                    }
                    n++;
                    break;
                case SH_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.sh", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ef.sheat;
                    }
                    n++;
                    break;
                case G_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.g", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ef.ssoil;
                    }
                    n++;
                    break;
                case ETP_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.etp", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ef.etp;
                    }
                    n++;
                    break;
                case ESNOW_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.esnow", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ef.esnow;
                    }
                    n++;
                    break;
                case ROOTW_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.rootw", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ws.soilw;
                    }
                    n++;
                    break;
                case SOILM_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.soilm", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ws.soilm;
                    }
                    n++;
                    break;
                case SOLAR_CTRL:
                    sprintf (pihm->prtctrl[n].name, "%s%s.solar", outputdir,
                        simulation);
                    pihm->prtctrl[n].intvl = pihm->ctrl.prtvrbl[i];
                    pihm->prtctrl[n].nvrbl = pihm->numele;
                    pihm->prtctrl[n].vrbl =
                        (double **)malloc (pihm->prtctrl[n].nvrbl *
                        sizeof (double *));
                    for (j = 0; j < pihm->numele; j++)
                    {
                        pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].ef.soldn;
                    }
                    n++;
                    break;
#endif
                default:
                    break;
            }
        }
    }

#ifdef _CYCLES_
                sprintf (pihm->prtctrl[n].name, "%s%s.biomass", outputdir,
                    simulation);
                pihm->prtctrl[n].intvl = 86400;
                pihm->prtctrl[n].nvrbl = pihm->numele;
                pihm->prtctrl[n].vrbl =
                    (double **)malloc (pihm->prtctrl[n].nvrbl *
                    sizeof (double *));
                for (j = 0; j < pihm->numele; j++)
                {
                    pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].comm.svBiomass;
                }
                n++;

                sprintf (pihm->prtctrl[n].name, "%s%s.eres", outputdir,
                    simulation);
                pihm->prtctrl[n].intvl = 3600;
                pihm->prtctrl[n].nvrbl = pihm->numele;
                pihm->prtctrl[n].vrbl =
                    (double **)malloc (pihm->prtctrl[n].nvrbl *
                    sizeof (double *));
                for (j = 0; j < pihm->numele; j++)
                {
                    pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].wf.eres;
                }
                n++;

                sprintf (pihm->prtctrl[n].name, "%s%s.radintcp", outputdir,
                    simulation);
                pihm->prtctrl[n].intvl = 86400;
                pihm->prtctrl[n].nvrbl = pihm->numele;
                pihm->prtctrl[n].vrbl =
                    (double **)malloc (pihm->prtctrl[n].nvrbl *
                    sizeof (double *));
                for (j = 0; j < pihm->numele; j++)
                {
                    pihm->prtctrl[n].vrbl[j] = &pihm->elem[j].comm.svRadiationInterception;
                }
                n++;
#endif

    pihm->ctrl.nprint = n;

    for (i = 0; i < pihm->ctrl.nprint; i++)
    {
        pihm->prtctrl[i].buffer =
            (double *)calloc (pihm->prtctrl[i].nvrbl, sizeof (double));
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
                    timestamp->tm_mday, timestamp->tm_hour,
                    timestamp->tm_min);
                for (j = 0; j < prtctrl[i].nvrbl; j++)
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
            if (NULL == fid)
            {
                printf ("ERROR: opening output files (.%s)!\n",
                    prtctrl[i].name);
                PihmExit (1);
            }

            outtime = (double)t;
            fwrite (&outtime, sizeof (double), 1, fid);
            for (j = 0; j < prtctrl[i].nvrbl; j++)
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

void PrtInit (pihm_struct pihm, char *simulation)
{
    FILE           *init_file;
    char            fn[MAXSTRING];
    int             i;
#ifdef _NOAH_
    int             j;
#endif

    sprintf (fn, "input/%s/%s.init", project, simulation);
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
