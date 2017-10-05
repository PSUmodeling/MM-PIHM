#include "pihm.h"

void MapOutput(pihm_struct pihm, const char *outputdir)
{
    int             i, j, k;
    int             n;
    char            ext[MAXSTRING];

    PIHMprintf(VL_VERBOSE, "\nInitializing PIHM output files\n");

    n = 0;

    for (i = 0; i < MAXPRINT; i++)
    {
        if (pihm->ctrl.prtvrbl[i] != 0)
        {
            switch (i)
            {
                case SURF_CTRL:
                    InitPrtVarCtrl(outputdir, "surf", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ws.surf;
                    }
                    n++;
                    break;
                case UNSAT_CTRL:
                    InitPrtVarCtrl(outputdir, "unsat", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ws.unsat;
                    }
                    n++;
                    break;
                case GW_CTRL:
                    InitPrtVarCtrl(outputdir, "gw", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ws.gw;
                    }
                    n++;
                    break;
                case RIVSTG_CTRL:
                    InitPrtVarCtrl(outputdir, "stage", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].ws.stage;
                    }
                    n++;
                    break;
                case RIVGW_CTRL:
                    InitPrtVarCtrl(outputdir, "rivgw", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->river[j].ws.gw;
                    }
                    n++;
                    break;
                case SNOW_CTRL:
                    InitPrtVarCtrl(outputdir, "snow", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ws.sneqv;
                    }
                    n++;
                    break;
                case CMC_CTRL:
                    InitPrtVarCtrl(outputdir, "is", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ws.cmc;
                    }
                    n++;
                    break;
                case INFIL_CTRL:
                    InitPrtVarCtrl(outputdir, "infil", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].wf.infil;
                    }
                    n++;
                    break;
                case RECHARGE_CTRL:
                    InitPrtVarCtrl(outputdir, "recharge", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].wf.rechg;
                    }
                    n++;
                    break;
                case EC_CTRL:
                    InitPrtVarCtrl(outputdir, "ec", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].wf.ec;
                    }
                    n++;
                    break;
                case ETT_CTRL:
                    InitPrtVarCtrl(outputdir, "ett", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].wf.ett;
                    }
                    n++;
                    break;
                case EDIR_CTRL:
                    InitPrtVarCtrl(outputdir, "edir", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].wf.edir;
                    }
                    n++;
                    break;
                case RIVFLX0_CTRL:
                    InitPrtVarCtrl(outputdir, "rivflx0", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].wf.rivflow[0];
                    }
                    n++;
                    break;
                case RIVFLX1_CTRL:
                    InitPrtVarCtrl(outputdir, "rivflx1", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].wf.rivflow[1];
                    }
                    n++;
                    break;
                case RIVFLX2_CTRL:
                    InitPrtVarCtrl(outputdir, "rivflx2", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].wf.rivflow[2];
                    }
                    n++;
                    break;
                case RIVFLX3_CTRL:
                    InitPrtVarCtrl(outputdir, "rivflx3", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].wf.rivflow[3];
                    }
                    n++;
                    break;
                case RIVFLX4_CTRL:
                    InitPrtVarCtrl(outputdir, "rivflx4", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].wf.rivflow[4];
                    }
                    n++;
                    break;
                case RIVFLX5_CTRL:
                    InitPrtVarCtrl(outputdir, "rivflx5", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].wf.rivflow[5];
                    }
                    n++;
                    break;
                case RIVFLX6_CTRL:
                    InitPrtVarCtrl(outputdir, "rivflx6", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].wf.rivflow[6];
                    }
                    n++;
                    break;
                case RIVFLX7_CTRL:
                    InitPrtVarCtrl(outputdir, "rivflx7", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].wf.rivflow[7];
                    }
                    n++;
                    break;
                case RIVFLX8_CTRL:
                    InitPrtVarCtrl(outputdir, "rivflx8", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].wf.rivflow[8];
                    }
                    n++;
                    break;
                case RIVFLX9_CTRL:
                    InitPrtVarCtrl(outputdir, "rivflx9", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].wf.rivflow[9];
                    }
                    n++;
                    break;
                case RIVFLX10_CTRL:
                    InitPrtVarCtrl(outputdir, "rivflx10", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].wf.rivflow[10];
                    }
                    n++;
                    break;
                case SUBFLX_CTRL:
                    for (k = 0; k < NUM_EDGE; k++)
                    {
                        sprintf(ext, "subflx%d", k);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].wf.subsurf[k];
                        }
                        n++;
                    }
                    break;
                case SURFFLX_CTRL:
                    for (k = 0; k < NUM_EDGE; k++)
                    {
                        sprintf(ext, "surfflx%d", k);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].wf.ovlflow[k];
                        }
                        n++;
                    }
                    break;
#ifdef _NOAH_
                case T1_CTRL:
                    InitPrtVarCtrl(outputdir, "t1", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].es.t1;
                    }
                    n++;
                    break;
                case STC_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf(ext, "stc%d", k);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            LS_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].es.stc[k];
                        }
                        n++;
                    }
                    break;
                case SMC_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf(ext, "smc%d", k);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].ws.smc[k];
                        }
                        n++;
                    }
                    break;
                case SH2O_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf(ext, "swc%d", k);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].ws.sh2o[k];
                        }
                        n++;
                    }
                    break;
                case SNOWH_CTRL:
                    InitPrtVarCtrl(outputdir, "snowh", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ps.snowh;
                    }
                    n++;
                    break;
                case ALBEDO_CTRL:
                    InitPrtVarCtrl(outputdir, "albedo", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].ps.albedo;
                    }
                    n++;
                    break;
                case LE_CTRL:
                    InitPrtVarCtrl(outputdir, "le", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ef.eta;
                    }
                    n++;
                    break;
                case SH_CTRL:
                    InitPrtVarCtrl(outputdir, "sh", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ef.sheat;
                    }
                    n++;
                    break;
                case G_CTRL:
                    InitPrtVarCtrl(outputdir, "g", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ef.ssoil;
                    }
                    n++;
                    break;
                case ETP_CTRL:
                    InitPrtVarCtrl(outputdir, "etp", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ef.etp;
                    }
                    n++;
                    break;
                case ESNOW_CTRL:
                    InitPrtVarCtrl(outputdir, "esnow", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ef.esnow;
                    }
                    n++;
                    break;
                case ROOTW_CTRL:
                    InitPrtVarCtrl(outputdir, "rootw", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ps.soilw;
                    }
                    n++;
                    break;
                case SOILM_CTRL:
                    InitPrtVarCtrl(outputdir, "soilm", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ws.soilm;
                    }
                    n++;
                    break;
                case SOLAR_CTRL:
                    InitPrtVarCtrl(outputdir, "solar", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ef.soldn;
                    }
                    n++;
                    break;
                case CH_CTRL:
                    InitPrtVarCtrl(outputdir, "ch", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ps.ch;
                    }
                    n++;
                    break;
#endif
#ifdef _BGC_
                case LAI_CTRL:
                    InitPrtVarCtrl(outputdir, "lai", pihm->ctrl.prtvrbl[i],
                        CN_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].ps.proj_lai;
                    }
                    n++;
                    break;
                case VEGC_CTRL:
                    InitPrtVarCtrl(outputdir, "vegc", pihm->ctrl.prtvrbl[i],
                        CN_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].summary.vegc;
                    }
                    n++;
                    break;
                case LITRC_CTRL:
                    InitPrtVarCtrl(outputdir, "litrc", pihm->ctrl.prtvrbl[i],
                        CN_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].summary.litrc;
                    }
                    n++;
                    break;
                case SOILC_CTRL:
                    InitPrtVarCtrl(outputdir, "soilc", pihm->ctrl.prtvrbl[i],
                        CN_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].summary.soilc;
                    }
                    n++;
                    break;
                case TOTALC_CTRL:
                    InitPrtVarCtrl(outputdir, "totalc", pihm->ctrl.prtvrbl[i],
                        CN_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].summary.totalc;
                    }
                    n++;
                    break;
                case NPP_CTRL:
                    InitPrtVarCtrl(outputdir, "npp", pihm->ctrl.prtvrbl[i],
                        CN_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].summary.daily_npp;
                    }
                    n++;
                    break;
                case NEP_CTRL:
                    InitPrtVarCtrl(outputdir, "nep", pihm->ctrl.prtvrbl[i],
                        CN_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].summary.daily_nep;
                    }
                    n++;
                    break;
                case NEE_CTRL:
                    InitPrtVarCtrl(outputdir, "nee", pihm->ctrl.prtvrbl[i],
                        CN_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].summary.daily_nee;
                    }
                    n++;
                    break;
                case GPP_CTRL:
                    InitPrtVarCtrl(outputdir, "gpp", pihm->ctrl.prtvrbl[i],
                        CN_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].summary.daily_gpp;
                    }
                    n++;
                    break;
                case SMINN_CTRL:
                    InitPrtVarCtrl(outputdir, "sminn", pihm->ctrl.prtvrbl[i],
                        CN_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].ns.sminn;
                    }
                    n++;
                    break;
                case LEAFC_CTRL:
                    InitPrtVarCtrl(outputdir, "leafc", pihm->ctrl.prtvrbl[i],
                        CN_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].cs.leafc;
                    }
                    n++;
                    break;
                case LIVESTEMC_CTRL:
                    InitPrtVarCtrl(outputdir, "livestemc",
                        pihm->ctrl.prtvrbl[i], CN_STEP, nelem,
                        &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].cs.livestemc;
                    }
                    n++;
                    break;
                case DEADSTEMC_CTRL:
                    InitPrtVarCtrl(outputdir, "deadstemc",
                        pihm->ctrl.prtvrbl[i], CN_STEP, nelem,
                        &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].cs.deadstemc;
                    }
                    n++;
                    break;
#endif
#ifdef _CYCLES_
                case BIOMASS_CTRL:
                    for (k = 0; k < pihm->elem[0].comm.NumCrop; k++)
                    {
                        sprintf(ext, "%s.biomass",
                            pihm->elem[0].comm.Crop[k].cropName);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            CN_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].comm.Crop[k].svBiomass;
                        }
                        n++;
                    }

                    InitPrtVarCtrl(outputdir, "comm.biomass",
                        pihm->ctrl.prtvrbl[i], CN_STEP, nelem,
                        &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].comm.svBiomass;
                    }
                    n++;
                    break;
                case RADNINTCP_CTRL:
                    for (k = 0; k < pihm->elem[0].comm.NumCrop; k++)
                    {
                        sprintf(ext, "%s.radintcp",
                            pihm->elem[0].comm.Crop[k].cropName);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            CN_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].comm.Crop[k].
                                svRadiationInterception;
                        }
                        n++;
                    }

                    InitPrtVarCtrl(outputdir, "comm.radintcp",
                        pihm->ctrl.prtvrbl[i], CN_STEP, nelem,
                        &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].comm.svRadiationInterception;
                    }
                    n++;
                    break;
                case WATER_STS_CTRL:
                    for (k = 0; k < pihm->elem[0].comm.NumCrop; k++)
                    {
                        sprintf(ext, "%s.waterstress",
                            pihm->elem[0].comm.Crop[k].cropName);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            CN_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].comm.Crop[k].svWaterStressFactor;
                        }
                        n++;
                    }

                    InitPrtVarCtrl(outputdir, "comm.waterstress",
                        pihm->ctrl.prtvrbl[i], CN_STEP, nelem,
                        &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].comm.svWaterStressFactor;
                    }
                    n++;
                    break;
                case N_STS_CTRL:
                    for (k = 0; k < pihm->elem[0].comm.NumCrop; k++)
                    {
                        sprintf(ext, "%s.nstress",
                            pihm->elem[0].comm.Crop[k].cropName);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            CN_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].comm.Crop[k].svN_StressFactor;
                        }
                        n++;
                    }

                    InitPrtVarCtrl(outputdir, "comm.nstress",
                        pihm->ctrl.prtvrbl[i], CN_STEP, nelem,
                        &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].comm.svN_StressFactor;
                    }
                    n++;
                    break;
                case CROP_TR_CTRL:
                    for (k = 0; k < pihm->elem[0].comm.NumCrop; k++)
                    {
                        sprintf(ext, "%s.transp",
                            pihm->elem[0].comm.Crop[k].cropName);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            LS_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].comm.Crop[k].svTranspiration;
                        }
                        n++;
                    }

                    InitPrtVarCtrl(outputdir, "comm.transp",
                        pihm->ctrl.prtvrbl[i], LS_STEP, nelem,
                        &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].comm.svTranspiration;
                    }
                    n++;
                    break;
                case CROP_POTTR_CTRL:
                    for (k = 0; k < pihm->elem[0].comm.NumCrop; k++)
                    {
                        sprintf("%s.pottransp",
                            pihm->elem[0].comm.Crop[k].cropName);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            LS_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].comm.Crop[k].
                                svTranspirationPotential;
                        }
                        n++;
                    }

                    InitPrtVarCtrl(outputdir, "comm.pottransp",
                        pihm->ctrl.prtvrbl[i], LS_STEP, nelem,
                        &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].comm.svTranspirationPotential;
                    }
                    n++;
                    break;
                case RES_EVAP_CTRL:
                    InitPrtVarCtrl(outputdir, "eres", pihm->ctrl.prtvrbl[i],
                        LS_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] = &pihm->elem[j].wf.eres;
                    }
                    n++;
                    break;
                case NO3_PROF_CTRL:
                    InitPrtVarCtrl(outputdir, "NO3", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].soil.NO3Profile;
                    }
                    n++;
                    break;
                case NO3_RIVER_CTRL:
                    InitPrtVarCtrl(outputdir, "rivNO3", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].NO3sol.soluteMass;
                    }
                    n++;
                    break;
                case NH4_PROF_CTRL:
                    InitPrtVarCtrl(outputdir, "NH4", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].soil.NH4Profile;
                    }
                    n++;
                    break;
                case NH4_RIVER_CTRL:
                    InitPrtVarCtrl(outputdir, "rivNH4", pihm->ctrl.prtvrbl[i],
                        HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->river[j].NH4sol.soluteMass;
                    }
                    n++;
                    break;
                case NO3_DENIT_CTRL:
                    InitPrtVarCtrl(outputdir, "NO3denitrif",
                        pihm->ctrl.prtvrbl[i], CN_STEP, nelem,
                        &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].soil.NO3_Denitrification;
                    }
                    n++;
                    break;
                case NO3_LEACH_CTRL:
                    for (k = 0; k < NUM_EDGE; k++)
                    {
                        sprintf(ext, "NO3leach%d", k);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].soil.NO3Leaching[k];
                        }
                        n++;
                    }
                    break;
                case NH4_LEACH_CTRL:
                    for (k = 0; k < NUM_EDGE; k++)
                    {
                        sprintf(ext, "NH4leach%d", k);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            HYDROL_STEP, nelem, &pihm->print.varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->elem[j].soil.NH4Leaching[k];
                        }
                        n++;
                    }
                    break;
                case NO3_LEACH_RIVER_CTRL:
                    for (k = 0; k < 4; k++)
                    {
                        sprintf(ext, "rivNO3leach%d", k);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                        for (j = 0; j < nriver; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->river[j].NO3Leaching[k];
                        }
                        n++;
                    }
                    break;
                case NH4_LEACH_RIVER_CTRL:
                    for (k = 0; k < 4; k++)
                    {
                        sprintf(ext, "rivNH4leach%d", k);
                        InitPrtVarCtrl(outputdir, ext, pihm->ctrl.prtvrbl[i],
                            HYDROL_STEP, nriver, &pihm->print.varctrl[n]);
                        for (j = 0; j < nriver; j++)
                        {
                            pihm->print.varctrl[n].var[j] =
                                &pihm->river[j].NH4Leaching[k];
                        }
                        n++;
                    }
                    break;
                case LAI_CTRL:
                    InitPrtVarCtrl(outputdir, "lai", pihm->ctrl.prtvrbl[i],
                        CN_STEP, nelem, &pihm->print.varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.varctrl[n].var[j] =
                            &pihm->elem[j].ps.proj_lai;
                    }
                    n++;
                    break;
#endif
                default:
                    break;
            }
        }
    }

    if (n > MAXPRINT)
    {
        PIHMprintf(VL_ERROR, "Error: Too many output files. ");
        PIHMprintf(VL_ERROR, "The maximum number of output files is %d.\n",
            MAXPRINT);
        PIHMexit(EXIT_FAILURE);
    }

    pihm->print.nprint = n;

    for (i = 0; i < pihm->print.nprint; i++)
    {
        if (spinup_mode)
        {
            pihm->print.varctrl[i].intvl = MONTHLY_OUTPUT;
        }

    }

    /*
     * Tecplot output
     */
    n = 0;

    for (i = 0; i < MAXPRINT; i++)
    {
        if (pihm->ctrl.tpprtvrbl[i] != 0)
        {
            switch (i)
            {
                case SURF_CTRL:
                    InitTecPrtVarCtrl(outputdir, "surf",
                        pihm->ctrl.tpprtvrbl[i], 0, HYDROL_STEP, nelem,
                        pihm->meshtbl.numnode, &pihm->print.tp_varctrl[n]);

                    for (j = 0; j < pihm->print.tp_varctrl[n].nnodes; j++)
                    {
                        pihm->print.tp_varctrl[n].x[j] = pihm->meshtbl.x[j];
                        pihm->print.tp_varctrl[n].y[j] = pihm->meshtbl.y[j];
                        pihm->print.tp_varctrl[n].zmax[j] =
                            pihm->meshtbl.zmax[j];
                        pihm->print.tp_varctrl[n].zmin[j] =
                            pihm->meshtbl.zmin[j];

                    }
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.tp_varctrl[n].var[j] =
                            &pihm->elem[j].ws.surf;
                        pihm->print.tp_varctrl[n].node0[j] =
                            pihm->elem[j].node[0];
                        pihm->print.tp_varctrl[n].node1[j] =
                            pihm->elem[j].node[1];
                        pihm->print.tp_varctrl[n].node2[j] =
                            pihm->elem[j].node[2];
                    }
                    n++;
                    break;
                case UNSAT_CTRL:
                    InitTecPrtVarCtrl(outputdir, "unsat",
                        pihm->ctrl.tpprtvrbl[i], 0, HYDROL_STEP, nelem,
                        pihm->meshtbl.numnode, &pihm->print.tp_varctrl[n]);
                    for (j = 0; j < pihm->print.tp_varctrl[n].nnodes; j++)
                    {
                        pihm->print.tp_varctrl[n].x[j] = pihm->meshtbl.x[j];
                        pihm->print.tp_varctrl[n].y[j] = pihm->meshtbl.y[j];
                        pihm->print.tp_varctrl[n].zmax[j] =
                            pihm->meshtbl.zmax[j];
                        pihm->print.tp_varctrl[n].zmin[j] =
                            pihm->meshtbl.zmin[j];
                    }
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.tp_varctrl[n].var[j] =
                            &pihm->elem[j].ws.unsat;
                        pihm->print.tp_varctrl[n].node0[j] =
                            pihm->elem[j].node[0];
                        pihm->print.tp_varctrl[n].node1[j] =
                            pihm->elem[j].node[1];
                        pihm->print.tp_varctrl[n].node2[j] =
                            pihm->elem[j].node[2];
                    }
                    n++;
                    break;
                case GW_CTRL:
                    InitTecPrtVarCtrl(outputdir, "gw", pihm->ctrl.tpprtvrbl[i],
                        0, HYDROL_STEP, nelem, pihm->meshtbl.numnode,
                        &pihm->print.tp_varctrl[n]);
                    for (j = 0; j < pihm->print.tp_varctrl[n].nnodes; j++)
                    {
                        pihm->print.tp_varctrl[n].x[j] = pihm->meshtbl.x[j];
                        pihm->print.tp_varctrl[n].y[j] = pihm->meshtbl.y[j];
                        pihm->print.tp_varctrl[n].zmax[j] =
                            pihm->meshtbl.zmax[j];
                        pihm->print.tp_varctrl[n].zmin[j] =
                            pihm->meshtbl.zmin[j];

                    }
                    for (j = 0; j < nelem; j++)
                    {
                        pihm->print.tp_varctrl[n].var[j] = &pihm->elem[j].ws.gw;
                        pihm->print.tp_varctrl[n].node0[j] =
                            pihm->elem[j].node[0];
                        pihm->print.tp_varctrl[n].node1[j] =
                            pihm->elem[j].node[1];
                        pihm->print.tp_varctrl[n].node2[j] =
                            pihm->elem[j].node[2];
                    }
                    n++;
                    break;
                case RIVSTG_CTRL:
                    InitTecPrtVarCtrl(outputdir, "stage",
                        pihm->ctrl.tpprtvrbl[i], 1, HYDROL_STEP, nriver, nriver,
                        &pihm->print.tp_varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.tp_varctrl[n].var[j] =
                            &pihm->river[j].ws.stage;
                        pihm->print.tp_varctrl[n].x[j] = pihm->river[j].topo.x;
                        pihm->print.tp_varctrl[n].y[j] = pihm->river[j].topo.y;
                        pihm->print.tp_varctrl[n].zmax[j] =
                            pihm->river[j].topo.zmax;
                        pihm->print.tp_varctrl[n].zmin[j] =
                            pihm->river[j].topo.zmin;
                    }
                    n++;
                    break;
                case RIVGW_CTRL:
                    InitTecPrtVarCtrl(outputdir, "rivgw",
                        pihm->ctrl.tpprtvrbl[i], 1, HYDROL_STEP, nriver, nriver,
                        &pihm->print.tp_varctrl[n]);

                    for (j = 0; j < nriver; j++)
                    {
                        pihm->print.tp_varctrl[n].var[j] =
                            &pihm->river[j].ws.gw;
                        pihm->print.tp_varctrl[n].x[j] = pihm->river[j].topo.x;
                        pihm->print.tp_varctrl[n].y[j] = pihm->river[j].topo.y;
                        pihm->print.tp_varctrl[n].zmax[j] =
                            pihm->river[j].topo.zmax;
                        pihm->print.tp_varctrl[n].zmin[j] =
                            pihm->river[j].topo.zmin;
                    }
                    n++;
                    break;
                default:
                    break;
            }
        }
    }

    if (n > MAXPRINT)
    {
        PIHMprintf(VL_ERROR, "Error: Too many Tecplot output files. ");
        PIHMprintf(VL_ERROR,
            "The maximum number of Tecplot output files is %d.\n", MAXPRINT);
        PIHMexit(EXIT_FAILURE);
    }

    pihm->print.ntpprint = n;
}

void InitPrtVarCtrl(const char *outputdir, const char *ext, int intvl,
    int upd_intvl, int nvar, varctrl_struct *varctrl)
{
    sprintf(varctrl->name, "%s%s.%s", outputdir, project, ext);
    varctrl->intvl = intvl;
    varctrl->upd_intvl = upd_intvl;
    varctrl->nvar = nvar;
    varctrl->var = (double **)malloc(nvar * sizeof(double *));
    varctrl->buffer = (double *)calloc(nvar, sizeof(double));
    varctrl->counter = 0;
}

void InitTecPrtVarCtrl(const char *outputdir, const char *ext, int intvl,
    int intr, int upd_intvl, int nvar, int nnode, varctrl_struct *varctrl)
{
    sprintf(varctrl->name, "%s%s.%s", outputdir, project, ext);
    varctrl->intvl = intvl;
    varctrl->intr = intr;
    varctrl->upd_intvl = upd_intvl;
    varctrl->nvar = nvar;
    varctrl->nnodes = nnode;

    varctrl->x = (double *)malloc(nnode * sizeof(double));
    varctrl->y = (double *)malloc(nnode * sizeof(double));
    varctrl->zmin = (double *)malloc(nnode * sizeof(double));
    varctrl->zmax = (double *)malloc(nnode * sizeof(double));
    if (intr == 0)
    {
        varctrl->node0 = (int *)malloc(nvar * sizeof(int));
        varctrl->node1 = (int *)malloc(nvar * sizeof(int));
        varctrl->node2 = (int *)malloc(nvar * sizeof(int));
    }
    varctrl->var = (double **)malloc(nvar * sizeof(double *));
    varctrl->buffer = (double *)calloc(nvar, sizeof(double));
    varctrl->counter = 0;
}
