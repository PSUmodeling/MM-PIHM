#include "pihm.h"

#if defined(_CYCLES_)
void MapOutput(const int *prtvrbl, const crop_struct crop[],
    const elem_struct *elem, const river_struct *river, const char *outputdir,
    print_struct *print)
#elif defined(_RT_)
void MapOutput(const int *prtvrbl, const chemtbl_struct chemtbl[],
    const rttbl_struct *rttbl, const elem_struct *elem,
    const river_struct *river, const char *outputdir, print_struct *print)
#else
void MapOutput(const int *prtvrbl, const elem_struct *elem,
    const river_struct *river, const char *outputdir, print_struct *print)
#endif
#if !defined(_LUMPED_)
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
                    InitPrtVarCtrl(outputdir, "surf", prtvrbl[i],
                        HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.surf;
                    }
                    n++;
                    break;
                case UNSAT_CTRL:
                    InitPrtVarCtrl(outputdir, "unsat", prtvrbl[i],
                        HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.unsat;
                    }
                    n++;
                    break;
                case GW_CTRL:
                    InitPrtVarCtrl(outputdir, "gw", prtvrbl[i],
                        HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.gw;
                    }
                    n++;
                    break;
                case RIVSTG_CTRL:
                    InitPrtVarCtrl(outputdir, "river.stage", prtvrbl[i],
                        HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].ws.stage;
                    }
                    n++;
                    break;
                case SNOW_CTRL:
                    InitPrtVarCtrl(outputdir, "snow", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.sneqv;
                    }
                    n++;
                    break;
                case CMC_CTRL:
#if defined(_CYCLES_OBSOLETE_)
                    InitPrtVarCtrl(outputdir, "stanresw", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.stanResidueWater;
                    }
                    n++;

                    InitPrtVarCtrl(outputdir, "flatresw", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.flatResidueWater;
                    }
                    n++;
#else
                    InitPrtVarCtrl(outputdir, "is", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.cmc;
                    }
                    n++;
#endif
                    break;
                case INFIL_CTRL:
                    InitPrtVarCtrl(outputdir, "infil", prtvrbl[i],
                        HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.eqv_infil;
                    }
                    n++;
                    break;
                case RECHARGE_CTRL:
                    InitPrtVarCtrl(outputdir, "recharge", prtvrbl[i],
                        HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.rechg;
                    }
                    n++;
                    break;
                case EC_CTRL:
#if defined(_CYCLES_)
                    InitPrtVarCtrl(outputdir, "eres", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
#else
                    InitPrtVarCtrl(outputdir, "ec", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
#endif
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.ec;
                    }
                    n++;
                    break;
                case ETT_CTRL:
                    InitPrtVarCtrl(outputdir, "ett", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.ett;
                    }
                    n++;
                    break;
                case EDIR_CTRL:
                    InitPrtVarCtrl(outputdir, "edir", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.edir;
                    }
                    n++;
                    break;
                case RIVFLX0_CTRL:
                    InitPrtVarCtrl(outputdir, "river.flx0", prtvrbl[i],
                        HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].wf.rivflow[0];
                    }
                    n++;
                    break;
                case RIVFLX1_CTRL:
                    InitPrtVarCtrl(outputdir, "river.flx1", prtvrbl[i],
                        HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].wf.rivflow[1];
                    }
                    n++;
                    break;
                case RIVFLX2_CTRL:
                    InitPrtVarCtrl(outputdir, "river.flx2", prtvrbl[i],
                        HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].wf.rivflow[2];
                    }
                    n++;
                    break;
                case RIVFLX3_CTRL:
                    InitPrtVarCtrl(outputdir, "river.flx3", prtvrbl[i],
                        HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].wf.rivflow[3];
                    }
                    n++;
                    break;
                case RIVFLX4_CTRL:
                    InitPrtVarCtrl(outputdir, "river.flx4", prtvrbl[i],
                        HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].wf.rivflow[4];
                    }
                    n++;
                    break;
                case RIVFLX5_CTRL:
                    InitPrtVarCtrl(outputdir, "river.flx5", prtvrbl[i],
                        HYDROL_STEP, nriver, &print->varctrl[n]);
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
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            HYDROL_STEP, nelem, &print->varctrl[n]);
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
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            HYDROL_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].wf.ovlflow[k];
                        }
                        n++;
                    }
                    break;
#if defined(_NOAH_)
                case T1_CTRL:
                    InitPrtVarCtrl(outputdir, "t1", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
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
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            LS_STEP, nelem, &print->varctrl[n]);
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
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            HYDROL_STEP, nelem, &print->varctrl[n]);
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
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            HYDROL_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].ws.swc[k];
                        }
                        n++;
                    }
                    break;
                case SNOWH_CTRL:
                    InitPrtVarCtrl(outputdir, "snowh", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.snowh;
                    }
                    n++;

                    InitPrtVarCtrl(outputdir, "iceh", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.iceh;
                    }
                    n++;
                    break;
                case ALBEDO_CTRL:
                    InitPrtVarCtrl(outputdir, "albedo", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.albedo;
                    }
                    n++;
                    break;
                case LE_CTRL:
                    InitPrtVarCtrl(outputdir, "le", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ef.eta;
                    }
                    n++;
                    break;
                case SH_CTRL:
                    InitPrtVarCtrl(outputdir, "sh", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ef.sheat;
                    }
                    n++;
                    break;
                case G_CTRL:
                    InitPrtVarCtrl(outputdir, "g", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ef.ssoil;
                    }
                    n++;
                    break;
                case ETP_CTRL:
                    InitPrtVarCtrl(outputdir, "etp", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ef.etp;
                    }
                    n++;
                    break;
                case ESNOW_CTRL:
                    InitPrtVarCtrl(outputdir, "esnow", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ef.esnow;
                    }
                    n++;
                    break;
                case ROOTW_CTRL:
                    InitPrtVarCtrl(outputdir, "rootw", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.soilw;
                    }
                    n++;
                    break;
                case SOILM_CTRL:
                    InitPrtVarCtrl(outputdir, "soilm", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.soilm;
                    }
                    n++;
                    break;
                case SOLAR_CTRL:
                    InitPrtVarCtrl(outputdir, "solar", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ef.soldn;
                    }
                    n++;
                    break;
                case CH_CTRL:
                    InitPrtVarCtrl(outputdir, "ch", prtvrbl[i],
                        LS_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.ch;
                    }
                    n++;
                    break;
#endif
#if defined(_BGC_)
# if defined(_LUMPEDBGC_)
                case LAI_CTRL:
                    InitPrtVarCtrl(outputdir, "lai", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].ps.proj_lai;
                    n++;
                    break;
                case NPP_CTRL:
                    InitPrtVarCtrl(outputdir, "npp", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.daily_npp;
                    n++;
                    break;
                case NEP_CTRL:
                    InitPrtVarCtrl(outputdir, "nep", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.daily_nep;
                    n++;
                    break;
                case NEE_CTRL:
                    InitPrtVarCtrl(outputdir, "nee", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.daily_nee;
                    n++;
                    break;
                case GPP_CTRL:
                    InitPrtVarCtrl(outputdir, "gpp", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.daily_gpp;
                    n++;
                    break;
                case MR_CTRL:
                    InitPrtVarCtrl(outputdir, "mr", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.daily_mr;
                    n++;
                    break;
                case GR_CTRL:
                    InitPrtVarCtrl(outputdir, "gr", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.daily_gr;
                    n++;
                    break;
                case HR_CTRL:
                    InitPrtVarCtrl(outputdir, "hr", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.daily_hr;
                    n++;
                    break;
                case FIRE_CTRL:
                    InitPrtVarCtrl(outputdir, "fire", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.daily_fire;
                    n++;
                    break;
                case LITFALLC_CTRL:
                    InitPrtVarCtrl(outputdir, "litfallc", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] =
                        &elem[LUMPEDBGC].summary.daily_litfallc;
                    n++;
                    break;
                case VEGC_CTRL:
                    InitPrtVarCtrl(outputdir, "vegc", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.vegc;
                    n++;
                    break;
                case AGC_CTRL:
                    InitPrtVarCtrl(outputdir, "agc", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.agc;
                    n++;
                    break;
                case LITRC_CTRL:
                    InitPrtVarCtrl(outputdir, "litrc", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.litrc;
                    n++;
                    break;
                case SOILC_CTRL:
                    InitPrtVarCtrl(outputdir, "soilc", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.soilc;
                    n++;
                    break;
                case TOTALC_CTRL:
                    InitPrtVarCtrl(outputdir, "totalc", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].summary.totalc;
                    n++;
                    break;
                case SMINN_CTRL:
                    InitPrtVarCtrl(outputdir, "sminn", prtvrbl[i],
                        CN_STEP, 1, &print->varctrl[n]);
                    print->varctrl[n].var[0] = &elem[LUMPEDBGC].ns.sminn;
                    n++;
                    break;
# else
                case LAI_CTRL:
                    InitPrtVarCtrl(outputdir, "lai", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.proj_lai;
                    }
                    n++;
                    break;
                case NPP_CTRL:
                    InitPrtVarCtrl(outputdir, "npp", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_npp;
                    }
                    n++;
                    break;
                case NEP_CTRL:
                    InitPrtVarCtrl(outputdir, "nep", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_nep;
                    }
                    n++;
                    break;
                case NEE_CTRL:
                    InitPrtVarCtrl(outputdir, "nee", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_nee;
                    }
                    n++;
                    break;
                case GPP_CTRL:
                    InitPrtVarCtrl(outputdir, "gpp", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_gpp;
                    }
                    n++;
                    break;
                case MR_CTRL:
                    InitPrtVarCtrl(outputdir, "mr", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_mr;
                    }
                    n++;
                    break;
                case GR_CTRL:
                    InitPrtVarCtrl(outputdir, "gr", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_gr;
                    }
                    n++;
                    break;
                case HR_CTRL:
                    InitPrtVarCtrl(outputdir, "hr", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_hr;
                    }
                    n++;
                    break;
                case FIRE_CTRL:
                    InitPrtVarCtrl(outputdir, "fire", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.daily_fire;
                    }
                    n++;
                    break;
                case LITFALLC_CTRL:
                    InitPrtVarCtrl(outputdir, "litfallc", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] =
                            &elem[j].summary.daily_litfallc;
                    }
                    n++;
                    break;
                case VEGC_CTRL:
                    InitPrtVarCtrl(outputdir, "vegc", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.vegc;
                    }
                    n++;
                    break;
                case AGC_CTRL:
                    InitPrtVarCtrl(outputdir, "agc", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.agc;
                    }
                    n++;
                    break;
                case LITRC_CTRL:
                    InitPrtVarCtrl(outputdir, "litrc", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.litrc;
                    }
                    n++;
                    break;
                case SOILC_CTRL:
                    InitPrtVarCtrl(outputdir, "soilc", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.soilc;
                    }
                    n++;
                    break;
                case TOTALC_CTRL:
                    InitPrtVarCtrl(outputdir, "totalc", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].summary.totalc;
                    }
                    n++;
                    break;
                case SMINN_CTRL:
                    InitPrtVarCtrl(outputdir, "sminn", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ns.sminn;
                    }
                    n++;
                    break;
# endif
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
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] =
                                &elem[j].crop[k].grain_yield;
                        }
                        n++;

                        sprintf(ext, "forage_yield.%s", crop[k].epc.name);
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] =
                                &elem[j].crop[k].forage_yield;
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
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] =
                                &elem[j].crop[k].shoot;
                        }
                        n++;

                        sprintf(ext, "root.%s", crop[k].epc.name);
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] =
                                &elem[j].crop[k].root;
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
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] =
                                &elem[j].crop[k].rad_intcp;
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
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] =
                                &elem[j].crop[k].water_stress;
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
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] =
                                &elem[j].crop[k].n_stress;
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
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] =
                                &elem[j].crop[k].transp;
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
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            CN_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] =
                                &elem[j].crop[k].transp_pot;
                        }
                        n++;
                    }
                    break;
                case N_PROFILE_CTRL:
                    InitPrtVarCtrl(outputdir, "NO3", prtvrbl[i],
                        HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.no3;
                    }
                    n++;

                    InitPrtVarCtrl(outputdir, "NH4", prtvrbl[i],
                        HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.nh4;
                    }
                    n++;
                    break;
                case N_RIVER_CTRL:
                    InitPrtVarCtrl(outputdir, "river.NO3", prtvrbl[i],
                        HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].ns.no3;
                    }
                    n++;

                    InitPrtVarCtrl(outputdir, "river.NH4", prtvrbl[i],
                        HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] = &river[j].ns.nh4;
                    }
                    n++;
                    break;
                case DENITRIF_CTRL:
                    InitPrtVarCtrl(outputdir, "denitrif", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.denitrif;
                    }
                    n++;
                    break;
                case LEACHING_CTRL:
                    InitPrtVarCtrl(outputdir, "river.NO3leaching", prtvrbl[i],
                        HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] =
                            &river[j].solute[NO3].flux[DOWN_CHANL2CHANL];
                    }
                    n++;

                    InitPrtVarCtrl(outputdir, "river.NH4leaching", prtvrbl[i],
                        HYDROL_STEP, nriver, &print->varctrl[n]);
                    for (j = 0; j < nriver; j++)
                    {
                        print->varctrl[n].var[j] =
                            &river[j].solute[NH4].flux[DOWN_CHANL2CHANL];
                    }
                    n++;
                    break;
                case LAI_CTRL:
                    InitPrtVarCtrl(outputdir, "lai", prtvrbl[i],
                        CN_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ps.proj_lai;
                    }
                    n++;
                    break;
#endif
#if defined(_FBR_)
                case FBRUNSAT_CTRL:
                    InitPrtVarCtrl(outputdir, "deep.unsat", prtvrbl[i],
                        HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.fbr_unsat;
                    }
                    n++;
                    break;
                case FBRGW_CTRL:
                    InitPrtVarCtrl(outputdir, "deep.gw", prtvrbl[i],
                        HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].ws.fbr_gw;
                    }
                    n++;
                    break;
                case FBRINFIL_CTRL:
                    InitPrtVarCtrl(outputdir, "deep.infil", prtvrbl[i],
                        HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.fbr_infil;
                    }
                    n++;
                    break;
                case FBRRECHG_CTRL:
                    InitPrtVarCtrl(outputdir, "deep.rechg", prtvrbl[i],
                        HYDROL_STEP, nelem, &print->varctrl[n]);
                    for (j = 0; j < nelem; j++)
                    {
                        print->varctrl[n].var[j] = &elem[j].wf.fbr_rechg;
                    }
                    n++;
                    break;
                case FBRFLOW_CTRL:
                    for (k = 0; k < NUM_EDGE; k++)
                    {
                        sprintf(ext, "deep.flow%d", k);
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            HYDROL_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].wf.fbrflow[k];
                        }
                        n++;
                    }
                    break;
#endif
#if defined(_RT_)
                case CHEM_CTRL:
                    /* Primary species */
                    for (k = 0; k < rttbl->num_stc; k++)
                    {
                        char            chemn[MAXSTRING];
                        Unwrap(chemn, chemtbl[k].name);

                        /* Unsaturated zone concentration */
                        sprintf(ext, "conc.%s", chemn);
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            RT_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].chms.prim_conc[k];
                        }
                        n++;
# if defined(_FBR_)
                        /* Fractured unsaturated bedrock layer concentration */
                        sprintf(ext, "deep.conc.%s", chemn);
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            RT_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] =
                                &elem[j].chms_geol.prim_conc[k];
                        }
                        n++;
# endif

                        /* River concentration */
                        sprintf(ext, "river.conc.%s", chemn);
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            RT_STEP, nriver, &print->varctrl[n]);
                        for (j = 0; j < nriver; j++)
                        {
                            print->varctrl[n].var[j] =
                                &river[j].chms.prim_conc[k];
                        }
                        n++;

                        /* River fluxes */
                        sprintf(ext, "river.chflx.%s", chemn);
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            RT_STEP, nriver, &print->varctrl[n]);
                        for (j = 0; j < nriver; j++)
                        {
                            print->varctrl[n].var[j] =
                                &river[j].solute[k].flux[DOWN_CHANL2CHANL];
                        }
                        n++;
                    }

                    /* Secondary species */
                    for (k = 0; k < rttbl->num_ssc; k++)
                    {
                        char            chemn[MAXSTRING];
                        Unwrap(chemn, chemtbl[k + rttbl->num_stc].name);

                        sprintf(ext, "conc.%s", chemn);
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            RT_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] = &elem[j].chms.sec_conc[k];
                        }
                        n++;

# if defined(_FBR_)
                        sprintf(ext, "deep.conc.%s", chemn);
                        InitPrtVarCtrl(outputdir, ext, prtvrbl[i],
                            RT_STEP, nelem, &print->varctrl[n]);
                        for (j = 0; j < nelem; j++)
                        {
                            print->varctrl[n].var[j] =
                                &elem[j].chms_geol.sec_conc[k];
                        }
                        n++;
# endif
                    }
#endif
                default:
                    break;
            }
        }
    }


    if (n > MAXPRINT)
    {
        pihm_printf(VL_ERROR, "Error: Too many output files. ");
        pihm_printf(VL_ERROR, "The maximum number of output files is %d.\n",
            MAXPRINT);
        pihm_exit(EXIT_FAILURE);
    }

    print->nprint = n;
}
#else
{
    int             i, k;
    int             n;
    char            ext[MAXSTRING];

    pihm_printf(VL_VERBOSE, "\nInitializing PIHM output files\n");

    n = 0;

    for (i = 0; i < 2; i++)
    {
        sprintf(ext, "elem%d.wflux", i + 1);
# if defined(_FBR_)
        InitPrtVarCtrl(outputdir, ext, DAILY_OUTPUT, HYDROL_STEP, 16,
            &print->varctrl[n]);
# else
        InitPrtVarCtrl(outputdir, ext, DAILY_OUTPUT, HYDROL_STEP, 11,
            &print->varctrl[n]);
# endif
        print->varctrl[n].var[0] = &elem[i].wf.infil;
        print->varctrl[n].var[1] = &elem[i].wf.rechg;
        print->varctrl[n].var[2] = &elem[i].wf.ec;
        print->varctrl[n].var[3] = &elem[i].wf.ett;
        print->varctrl[n].var[4] = &elem[i].wf.edir;
        for (k = 0; k < NUM_EDGE; k++)
        {
            print->varctrl[n].var[5 + k] = &elem[i].wf.subsurf[k];
            print->varctrl[n].var[8 + k] = &elem[i].wf.ovlflow[k];
        }
# if defined(_FBR_)
        print->varctrl[n].var[11] = &elem[i].wf.fbr_infil;
        print->varctrl[n].var[12] = &elem[i].wf.fbr_rechg;
        for (k = 0; k < NUM_EDGE; k++)
        {
            print->varctrl[n].var[13 + k] = &elem[i].wf.fbrflow[k];
        }
# endif
        n++;

        sprintf(ext, "elem%d.wstate", i + 1);
# if defined(_FBR_)
        InitPrtVarCtrl(outputdir, ext, DAILY_OUTPUT, HYDROL_STEP, 7,
            &print->varctrl[n]);
# else
        InitPrtVarCtrl(outputdir, ext, DAILY_OUTPUT, HYDROL_STEP, 5,
            &print->varctrl[n]);
# endif
        print->varctrl[n].var[0] = &elem[i].ws.cmc;
        print->varctrl[n].var[1] = &elem[i].ws.sneqv;
        print->varctrl[n].var[2] = &elem[i].ws.surf;
        print->varctrl[n].var[3] = &elem[i].ws.unsat;
        print->varctrl[n].var[4] = &elem[i].ws.gw;
# if defined(_FBR_)
        print->varctrl[n].var[5] = &elem[i].ws.fbr_unsat;
        print->varctrl[n].var[6] = &elem[i].ws.fbr_gw;
# endif
        n++;

# if defined(_NOAH_)
        sprintf(ext, "elem%d.smc", i + 1);
        InitPrtVarCtrl(outputdir, ext, DAILY_OUTPUT, HYDROL_STEP, MAXLYR,
            &print->varctrl[n]);
        for (k = 0; k < MAXLYR; k++)
        {
            print->varctrl[n].var[k] = &elem[i].ws.smc[k];
        }
        n++;

        sprintf(ext, "elem%d.swc", i + 1);
        InitPrtVarCtrl(outputdir, ext, DAILY_OUTPUT, HYDROL_STEP, MAXLYR,
            &print->varctrl[n]);
        for (k = 0; k < MAXLYR; k++)
        {
            print->varctrl[n].var[k] = &elem[i].ws.swc[k];
        }
        n++;

        sprintf(ext, "elem%d.stc", i + 1);
        InitPrtVarCtrl(outputdir, ext, DAILY_OUTPUT, HYDROL_STEP, MAXLYR,
            &print->varctrl[n]);
        for (k = 0; k < MAXLYR; k++)
        {
            print->varctrl[n].var[k] = &elem[i].es.stc[k];
        }
        n++;

        sprintf(ext, "elem%d.ls", i + 1);
        InitPrtVarCtrl(outputdir, ext, DAILY_OUTPUT, LS_STEP, 12,
            &print->varctrl[n]);
        print->varctrl[n].var[0] = &elem[i].es.t1;
        print->varctrl[n].var[1] = &elem[i].ps.snowh;
        print->varctrl[n].var[2] = &elem[i].ps.albedo;
        print->varctrl[n].var[3] = &elem[i].ef.eta;
        print->varctrl[n].var[4] = &elem[i].ef.sheat;
        print->varctrl[n].var[5] = &elem[i].ef.ssoil;
        print->varctrl[n].var[6] = &elem[i].ef.etp;
        print->varctrl[n].var[7] = &elem[i].ef.esnow;
        print->varctrl[n].var[8] = &elem[i].ps.soilw;
        print->varctrl[n].var[9] = &elem[i].ws.soilm;
        print->varctrl[n].var[10] = &elem[i].ef.soldn;
        print->varctrl[n].var[11] = &elem[i].ps.ch;
        n++;
# endif

# if defined(_RT_)
        sprintf(ext, "elem%d.conc", i + 1);
        InitPrtVarCtrl(outputdir, ext, DAILY_OUTPUT, RT_STEP,
            rttbl->num_stc + rttbl->num_ssc, &print->varctrl[n]);
        for (k = 0; k < rttbl->num_stc; k++)
        {
            print->varctrl[n].var[k] = &elem[i].chms.prim_conc[k];
        }
        for (k = 0; k < rttbl->num_ssc; k++)
        {
            print->varctrl[n].var[rttbl->num_stc + k] =
                &elem[i].chms.sec_conc[k];
        }
        n++;

#  if defined(_FBR_)
        sprintf(ext, "elem%d.deep.conc", i + 1);
        InitPrtVarCtrl(outputdir, ext, DAILY_OUTPUT, RT_STEP,
            rttbl->num_stc + rttbl->num_ssc, &print->varctrl[n]);
        for (k = 0; k < rttbl->num_stc; k++)
        {
            print->varctrl[n].var[k] = &elem[i].chms_geol.prim_conc[k];
        }
        for (k = 0; k < rttbl->num_ssc; k++)
        {
            print->varctrl[n].var[rttbl->num_stc + k] =
                &elem[i].chms_geol.sec_conc[k];
        }
        n++;
#  endif
# endif
    }

    InitPrtVarCtrl(outputdir, "river.wflux", DAILY_OUTPUT, HYDROL_STEP,
        NUM_RIVFLX, &print->varctrl[n]);
    for (k = 0; k < NUM_RIVFLX; k++)
    {
        print->varctrl[n].var[k] = &river[0].wf.rivflow[k];
    }
    n++;

    InitPrtVarCtrl(outputdir, "river.wstate", DAILY_OUTPUT, HYDROL_STEP, 1,
        &print->varctrl[n]);
    print->varctrl[n].var[0] = &river[0].ws.stage;
    n++;

# if defined(_RT_)
    InitPrtVarCtrl(outputdir, "river.conc", DAILY_OUTPUT, RT_STEP,
        rttbl->num_stc, &print->varctrl[n]);
    for (k = 0; k < rttbl->num_stc; k++)
    {
        print->varctrl[n].var[k] = &river[0].chms.prim_conc[k];
    }
    n++;

    InitPrtVarCtrl(outputdir, "leach", DAILY_OUTPUT, RT_STEP, rttbl->num_stc,
        &print->varctrl[n]);
    for (k = 0; k < rttbl->num_stc; k++)
    {
        print->varctrl[n].var[k] = &river[0].solute[k].flux[DOWN_CHANL2CHANL];
    }
    n++;

#  if defined(_FBR_)
    InitPrtVarCtrl(outputdir, "left_leach", DAILY_OUTPUT, RT_STEP,
        rttbl->num_stc, &print->varctrl[n]);
    for (k = 0; k < rttbl->num_stc; k++)
    {
        print->varctrl[n].var[k] = &river[0].solute[k].flux[LEFT_FBR2CHANL];
    }
    n++;

    InitPrtVarCtrl(outputdir, "right_leach", DAILY_OUTPUT, RT_STEP,
        rttbl->num_stc, &print->varctrl[n]);
    for (k = 0; k < rttbl->num_stc; k++)
    {
        print->varctrl[n].var[k] = &river[0].solute[k].flux[RIGHT_FBR2CHANL];
    }
    n++;
#  endif
# endif

    if (n > MAXPRINT)
    {
        pihm_printf(VL_ERROR, "Error: Too many output files. ");
        pihm_printf(VL_ERROR, "The maximum number of output files is %d.\n",
            MAXPRINT);
        pihm_exit(EXIT_FAILURE);
    }

    print->nprint = n;
}
#endif

void InitPrtVarCtrl(const char *outputdir, const char *ext, int intvl,
    int upd_intvl, int nvar, varctrl_struct *varctrl)
{
    sprintf(varctrl->name, "%s%s.%s", outputdir, project, ext);
    if (spinup_mode)
    {
        /* When spinning-up, print interval is set to monthly */
        varctrl->intvl = MONTHLY_OUTPUT;
    }
    else
    {
        varctrl->intvl = intvl;
    }
    varctrl->upd_intvl = upd_intvl;
    varctrl->nvar = nvar;
    varctrl->var = (const double **)malloc(nvar * sizeof(double *));
    varctrl->buffer = (double *)calloc(nvar, sizeof(double));
    varctrl->counter = 0;
}

#if defined(_RT_)
void Unwrap(char *str, const char *str0)
{
    int             i, j = 0;

    for (i = 0; i < (int)strlen(str0); i++)
    {
        if (str0[i] != '\'')
        {
            str[j] = str0[i];
            j++;
        }
    }

    str[j] = '\0';
}
#endif
