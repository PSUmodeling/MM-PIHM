#include "pihm.h"

#if defined(_CYCLES_)
void MapOutput(const char outputdir[], const int prtvrbl[], const crop_struct crop[], const elem_struct elem[],
    const river_struct river[], print_struct *print)
#elif defined(_RT_)
void MapOutput(const char outputdir[], const int prtvrbl[], const chemtbl_struct chemtbl[], const rttbl_struct *rttbl,
    const elem_struct elem[], const river_struct river[], print_struct *print)
#else
void MapOutput(const char outputdir[], const int prtvrbl[], const elem_struct elem[], const river_struct river[],
    print_struct *print)
#endif
{
    int             i, j, k;
    int             n;
    char            ext[MAXSTRING];

    pihm_printf(VL_VERBOSE, "\nInitializing PIHM output files\n");

    n = 0;

    for (i = 0; i < MAXPRINT; i++)
    {
        if (prtvrbl[i] != 0)
        {
            switch (i)
            {
                case SURF_CTRL:
                    InitPrintCtrl(outputdir, "surf", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.surf;
                    }
                    n++;
                    break;
                case UNSAT_CTRL:
                    InitPrintCtrl(outputdir, "unsat", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.unsat;
                    }
                    n++;
                    break;
                case GW_CTRL:
                    InitPrintCtrl(outputdir, "gw", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.gw;
                    }
                    n++;
                    break;
                case STAGE_CTRL:
                    InitPrintCtrl(outputdir, "river.stage", prtvrbl[i], HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].ws.stage;
                    }
                    n++;
                    break;
                case SNOW_CTRL:
                    InitPrintCtrl(outputdir, "snow", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.sneqv;
                    }
                    n++;
                    break;
                case CMC_CTRL:
                    InitPrintCtrl(outputdir, "is", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.cmc;
                    }
                    n++;
                    break;
                case INFIL_CTRL:
                    InitPrintCtrl(outputdir, "infil", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.eqv_infil;
                    }
                    n++;
                    break;
                case RECHARGE_CTRL:
                    InitPrintCtrl(outputdir, "recharge", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.recharge;
                    }
                    n++;
                    break;
                case EC_CTRL:
#if defined(_CYCLES_)
                    InitPrintCtrl(outputdir, "eres", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
#else
                    InitPrintCtrl(outputdir, "ec", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
#endif
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.ec;
                    }
                    n++;
                    break;
                case ETT_CTRL:
                    InitPrintCtrl(outputdir, "ett", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.ett;
                    }
                    n++;
                    break;
                case EDIR_CTRL:
                    InitPrintCtrl(outputdir, "edir", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.edir;
                    }
                    n++;
                    break;
                case RIVFLX0_CTRL:
                    InitPrintCtrl(outputdir, "river.flx0", prtvrbl[i], HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].wf.rivflow[0];
                    }
                    n++;
                    break;
                case RIVFLX1_CTRL:
                    InitPrintCtrl(outputdir, "river.flx1", prtvrbl[i], HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].wf.rivflow[1];
                    }
                    n++;
                    break;
                case RIVFLX2_CTRL:
                    InitPrintCtrl(outputdir, "river.flx2", prtvrbl[i], HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].wf.rivflow[2];
                    }
                    n++;
                    break;
                case RIVFLX3_CTRL:
                    InitPrintCtrl(outputdir, "river.flx3", prtvrbl[i], HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].wf.rivflow[3];
                    }
                    n++;
                    break;
                case RIVFLX4_CTRL:
                    InitPrintCtrl(outputdir, "river.flx4", prtvrbl[i], HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].wf.rivflow[4];
                    }
                    n++;
                    break;
                case RIVFLX5_CTRL:
                    InitPrintCtrl(outputdir, "river.flx5", prtvrbl[i], HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].wf.rivflow[5];
                    }
                    n++;
                    break;
                case SUBFLX_CTRL:
                    for (k = 0; k < NUM_EDGE; k++)
                    {
                        sprintf(ext, "subflx%d", k);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].wf.subsurf[k];
                        }
                        n++;
                    }
                    break;
                case SURFFLX_CTRL:
                    for (k = 0; k < NUM_EDGE; k++)
                    {
                        sprintf(ext, "surfflx%d", k);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].wf.overland[k];
                        }
                        n++;
                    }
                    break;
#if defined(_NOAH_)
                case T1_CTRL:
                    InitPrintCtrl(outputdir, "t1", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].es.t1;
                    }
                    n++;
                    break;
                case STC_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf(ext, "stc%d", k);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].es.stc[k];
                        }
                        n++;
                    }
                    break;
                case SMC_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf(ext, "smc%d", k);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].ws.smc[k];
                        }
                        n++;
                    }
                    break;
                case SH2O_CTRL:
                    for (k = 0; k < MAXLYR; k++)
                    {
                        sprintf(ext, "swc%d", k);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].ws.swc[k];
                        }
                        n++;
                    }
                    break;
                case SNOWH_CTRL:
                    InitPrintCtrl(outputdir, "snowh", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.snowh;
                    }
                    n++;

                    InitPrintCtrl(outputdir, "iceh", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.iceh;
                    }
                    n++;
                    break;
                case ALBEDO_CTRL:
                    InitPrintCtrl(outputdir, "albedo", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.albedo;
                    }
                    n++;
                    break;
                case LE_CTRL:
                    InitPrintCtrl(outputdir, "le", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ef.eta;
                    }
                    n++;
                    break;
                case SH_CTRL:
                    InitPrintCtrl(outputdir, "sh", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ef.sheat;
                    }
                    n++;
                    break;
                case G_CTRL:
                    InitPrintCtrl(outputdir, "g", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ef.ssoil;
                    }
                    n++;
                    break;
                case ETP_CTRL:
                    InitPrintCtrl(outputdir, "etp", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ef.etp;
                    }
                    n++;
                    break;
                case ESNOW_CTRL:
                    InitPrintCtrl(outputdir, "esnow", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ef.esnow;
                    }
                    n++;
                    break;
                case ROOTW_CTRL:
                    InitPrintCtrl(outputdir, "rootw", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.soilw;
                    }
                    n++;
                    break;
                case SOILM_CTRL:
                    InitPrintCtrl(outputdir, "soilm", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.soilm;
                    }
                    n++;
                    break;
                case SOLAR_CTRL:
                    InitPrintCtrl(outputdir, "solar", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ef.soldn;
                    }
                    n++;
                    break;
                case CH_CTRL:
                    InitPrintCtrl(outputdir, "ch", prtvrbl[i], LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.ch;
                    }
                    n++;
                    break;
#endif
#if defined(_BGC_)
                case LAI_CTRL:
                    InitPrintCtrl(outputdir, "lai", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.proj_lai;
                    }
                    n++;
                    break;
                case NPP_CTRL:
                    InitPrintCtrl(outputdir, "npp", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_npp;
                    }
                    n++;
                    break;
                case NEP_CTRL:
                    InitPrintCtrl(outputdir, "nep", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_nep;
                    }
                    n++;
                    break;
                case NEE_CTRL:
                    InitPrintCtrl(outputdir, "nee", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_nee;
                    }
                    n++;
                    break;
                case GPP_CTRL:
                    InitPrintCtrl(outputdir, "gpp", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_gpp;
                    }
                    n++;
                    break;
                case MR_CTRL:
                    InitPrintCtrl(outputdir, "mr", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_mr;
                    }
                    n++;
                    break;
                case GR_CTRL:
                    InitPrintCtrl(outputdir, "gr", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_gr;
                    }
                    n++;
                    break;
                case HR_CTRL:
                    InitPrintCtrl(outputdir, "hr", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_hr;
                    }
                    n++;
                    break;
                case FIRE_CTRL:
                    InitPrintCtrl(outputdir, "fire", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_fire;
                    }
                    n++;
                    break;
                case LITFALLC_CTRL:
                    InitPrintCtrl(outputdir, "litfallc", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_litfallc;
                    }
                    n++;
                    break;
                case VEGC_CTRL:
                    InitPrintCtrl(outputdir, "vegc", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.vegc;
                    }
                    n++;
                    break;
                case AGC_CTRL:
                    InitPrintCtrl(outputdir, "agc", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.agc;
                    }
                    n++;
                    break;
                case LITRC_CTRL:
                    InitPrintCtrl(outputdir, "litrc", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.litrc;
                    }
                    n++;
                    break;
                case SOILC_CTRL:
                    InitPrintCtrl(outputdir, "soilc", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.soilc;
                    }
                    n++;
                    break;
                case TOTALC_CTRL:
                    InitPrintCtrl(outputdir, "totalc", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.totalc;
                    }
                    n++;
                    break;
                case SMINN_CTRL:
                    InitPrintCtrl(outputdir, "sminn", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ns.sminn;
                    }
                    n++;
                    break;
#endif
#if defined(_CYCLES_)
                case YIELD_CTRL:
                    for (k = 0; k < MAXCROP && '\0' != crop[k].epc.name[0]; k++)
                    {
                        if (crop[k].stage_growth == NOT_USED)
                        {
                            continue;
                        }

                        sprintf(ext, "grain_yield.%s", crop[k].epc.name);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].crop[k].grain_yield;
                        }
                        n++;

                        sprintf(ext, "forage_yield.%s", crop[k].epc.name);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].crop[k].forage_yield;
                        }
                        n++;
                    }
                    break;
                case BIOMASS_CTRL:
                    for (k = 0; k < MAXCROP && '\0' != crop[k].epc.name[0]; k++)
                    {
                        if (crop[k].stage_growth == NOT_USED)
                        {
                            continue;
                        }

                        sprintf(ext, "shoot.%s", crop[k].epc.name);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].crop[k].shoot;
                        }
                        n++;

                        sprintf(ext, "root.%s", crop[k].epc.name);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].crop[k].root;
                        }
                        n++;
                    }
                    break;
                case RADNINTCP_CTRL:
                    for (k = 0; k < MAXCROP && '\0' != crop[k].epc.name[0]; k++)
                    {
                        if (crop[k].stage_growth == NOT_USED)
                        {
                            continue;
                        }

                        sprintf(ext, "radintcp.%s", crop[k].epc.name);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].crop[k].rad_intcp;
                        }
                        n++;
                    }
                    break;
                case WATER_STS_CTRL:
                    for (k = 0; k < MAXCROP && '\0' != crop[k].epc.name[0]; k++)
                    {
                        if (crop[k].stage_growth == NOT_USED)
                        {
                            continue;
                        }

                        sprintf(ext, "wstress.%s", crop[k].epc.name);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].crop[k].water_stress;
                        }
                        n++;
                    }
                    break;
                case N_STS_CTRL:
                    for (k = 0; k < MAXCROP && '\0' != crop[k].epc.name[0]; k++)
                    {
                        if (crop[k].stage_growth == NOT_USED)
                        {
                            continue;
                        }

                        sprintf(ext, "nstress.%s", crop[k].epc.name);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].crop[k].n_stress;
                        }
                        n++;
                    }
                    break;
                case CROP_TR_CTRL:
                    for (k = 0; k < MAXCROP && '\0' != crop[k].epc.name[0]; k++)
                    {
                        if (crop[k].stage_growth == NOT_USED)
                        {
                            continue;
                        }

                        sprintf(ext, "transp.%s", crop[k].epc.name);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].crop[k].transp;
                        }
                        n++;
                    }
                    break;
                case CROP_POTTR_CTRL:
                    for (k = 0; k < MAXCROP && '\0' != crop[k].epc.name[0]; k++)
                    {
                        if (crop[k].stage_growth == NOT_USED)
                        {
                            continue;
                        }

                        sprintf(ext, "pottransp.%s", crop[k].epc.name);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].crop[k].transp_pot;
                        }
                        n++;
                    }
                    break;
                case N_PROFILE_CTRL:
                    InitPrintCtrl(outputdir, "NO3", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.no3;
                    }
                    n++;

                    InitPrintCtrl(outputdir, "NH4", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.nh4;
                    }
                    n++;

                    InitPrintCtrl(outputdir, "SON", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.son;
                    }
                    n++;
                    break;
                case N_RIVER_CTRL:
                    InitPrintCtrl(outputdir, "river.NO3", prtvrbl[i], HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].ns.no3;
                    }
                    n++;

                    InitPrintCtrl(outputdir, "river.NH4", prtvrbl[i], HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].ns.nh4;
                    }
                    n++;
                    break;
                case DENITRIF_CTRL:
                    InitPrintCtrl(outputdir, "denitrif", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.denitrif;
                    }
                    n++;
                    break;
                case NITRIF_CTRL:
                    InitPrintCtrl(outputdir, "nitrif", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.nitrif;
                    }
                    n++;
                    break;
                case IMMOBIL_CTRL:
                    InitPrintCtrl(outputdir, "immobil", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.immobil;
                    }
                    n++;
                    break;
                case MINERAL_CTRL:
                    InitPrintCtrl(outputdir, "mineral", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.mineral;
                    }
                    n++;
                    break;
                case VOLATIL_CTRL:
                    InitPrintCtrl(outputdir, "volatil", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.volatil;
                    }
                    n++;
                    break;
                case LEACHING_CTRL:
                    InitPrintCtrl(outputdir, "NO3leaching", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].nf.no3_leach;
                    }
                    n++;

                    InitPrintCtrl(outputdir, "NH4leaching", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].nf.nh4_leach;
                    }
                    n++;
                    break;
                case SOC_CTRL:
                    InitPrintCtrl(outputdir, "soc", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.soc;
                    }
                    n++;
                    break;
                case N2O_CTRL:
                    InitPrintCtrl(outputdir, "n2o", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.n2o_emis;
                    }
                    n++;
                    break;
                case N_HARVEST_CTRL:
                    InitPrintCtrl(outputdir, "n_harvest", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].nf.harvest;
                    }
                    n++;
                    break;
                case N_INPUT_CTRL:
                    InitPrintCtrl(outputdir, "n_fert", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.n_fert;
                    }
                    n++;

                    InitPrintCtrl(outputdir, "n_auto", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].nf.auto_added;
                    }
                    n++;

                    InitPrintCtrl(outputdir, "n_fix", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].nf.fixation;
                    }
                    n++;
                    break;
                case NEP_CTRL:
                    InitPrintCtrl(outputdir, "nep", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.nep;
                    }
                    n++;
                    break;
                case LAI_CTRL:
                    InitPrintCtrl(outputdir, "lai", prtvrbl[i], CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.proj_lai;
                    }
                    n++;
                    break;
#endif
#if defined(_DGW_)
                case GEOLUNSAT_CTRL:
                    InitPrintCtrl(outputdir, "deep.unsat", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.unsat_geol;
                    }
                    n++;
                    break;
                case GEOLGW_CTRL:
                    InitPrintCtrl(outputdir, "deep.gw", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.gw_geol;
                    }
                    n++;
                    break;
                case GEOLINFIL_CTRL:
                    InitPrintCtrl(outputdir, "deep.infil", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.infil_geol;
                    }
                    n++;
                    break;
                case GEOLRECHG_CTRL:
                    InitPrintCtrl(outputdir, "deep.recharge", prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.rechg_geol;
                    }
                    n++;
                    break;
                case DGWFLOW_CTRL:
                    for (k = 0; k < NUM_EDGE; k++)
                    {
                        sprintf(ext, "deep.flow%d", k);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], HYDROL_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].wf.dgw[k];
                        }
                        n++;
                    }
                    break;
#endif
#if defined(_RT_)
                case CHEM_CTRL:
                    // Primary species
                    for (k = 0; k < rttbl->num_stc; k++)
                    {
                        char            chemn[MAXSTRING];
                        Unwrap(chemtbl[k].name, chemn);

                        // Unsaturated zone concentration
                        sprintf(ext, "conc.%s", chemn);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], RT_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].chms.prim_conc[k];
                        }
                        n++;
# if defined(_DGW_)
                        // Deep zone concentration
                        sprintf(ext, "deep.conc.%s", chemn);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], RT_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].chms_geol.prim_conc[k];
                        }
                        n++;
# endif

                        // River concentration
                        sprintf(ext, "river.conc.%s", chemn);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], RT_STEP, nriver, &print->varctrl[n]);
                        for (j = 0; j < nriver; j++)
                        {
                            print->varctrl[n].var[j] = &river[j].chms.prim_conc[k];
                        }
                        n++;

                        // River fluxes
                        sprintf(ext, "river.chflx.%s", chemn);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], RT_STEP, nriver, &print->varctrl[n]);
                        for (j = 0; j < nriver; j++)
                        {
                            print->varctrl[n].var[j] = &river[j].solute[k].flux[DOWNSTREAM];
                        }
                        n++;
                    }

                    // Secondary species
                    for (k = 0; k < rttbl->num_ssc; k++)
                    {
                        char            chemn[MAXSTRING];
                        Unwrap(chemtbl[k + rttbl->num_stc].name, chemn);

                        sprintf(ext, "conc.%s", chemn);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], RT_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].chms.sec_conc[k];
                        }
                        n++;

# if defined(_DGW_)
                        sprintf(ext, "deep.conc.%s", chemn);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], RT_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].chms_geol.sec_conc[k];
                        }
                        n++;
# endif

                        // River concentration
                        sprintf(ext, "river.conc.%s", chemn);
                        InitPrintCtrl(outputdir, ext, prtvrbl[i], RT_STEP, nriver, &print->varctrl[n]);
                        for (j = 0; j < nriver; j++)
                        {
                            print->varctrl[n].var[j] = &river[j].chms.sec_conc[k];
                        }
                        n++;
                    }
#endif
                default:
                    break;
            }
        }
    }


    if (n > MAXPRINT)
    {
        pihm_printf(VL_ERROR, "Error: Too many output files. The maximum number of output files is %d.\n", MAXPRINT);
        pihm_exit(EXIT_FAILURE);
    }

    print->nprint = n;
}

void InitPrintCtrl(const char outputdir[], const char ext[], int intvl, int upd_intvl, int nvar,
    varctrl_struct *varctrl)
{
    sprintf(varctrl->name, "%s%s.%s", outputdir, project, ext);

    // When spinning-up, print interval is set to monthly
    varctrl->intvl = (spinup_mode) ? MONTHLY_OUTPUT: intvl;
    varctrl->upd_intvl = upd_intvl;
    varctrl->nvar = nvar;
    varctrl->var = (const double **)malloc(nvar * sizeof(double *));
    varctrl->buffer = (double *)calloc(nvar, sizeof(double));
    varctrl->counter = 0;
}

#if defined(_RT_)
void Unwrap(const char wrapped_str[], char str[])
{
    int             i, j = 0;

    for (i = 0; i < (int)strlen(wrapped_str); i++)
    {
        if (wrapped_str[i] != '\'')
        {
            str[j] = wrapped_str[i];
            j++;
        }
    }

    str[j] = '\0';
}

#endif
