/*******************************************************************************
* RT-Flux-PIHM is a finite volume based, reactive transport module that operates
* on top of the hydrological land surface processes described by Flux-PIHM.
* RT-Flux-PIHM tracks the transportation and reaction in a given watershed. It
* uses operator splitting technique to couple transport and reaction.
*****************************************************************************/
#include "pihm.h"

/* Begin global variable definition (MACRO) */
#define ZERO   1E-20
#define LINE_WIDTH 512
#define WORDS_LINE 40
#define WORD_WIDTH 80
#define INFTYSMALL  1E-6

void InitChem(const char cdbs_filen[], const char cini_filen[],
    const ctrl_struct *ctrl, const calib_struct *cal, forc_struct *forc,
    chemtbl_struct chemtbl[], kintbl_struct kintbl[], rttbl_struct *rttbl,
    elem_struct elem[], N_Vector CV_Y)
{
    int             i, j, k;
    int             PRCP_VOL;
    int             BOUND_VOL;
    int             chem_ind;
    FILE           *fp;

    fp = fopen(cdbs_filen, "r");
    CheckFile(fp, cdbs_filen);

    ReadCini(cini_filen, chemtbl, rttbl->NumStc, elem);

    /*
     * Look up database to find required parameters and dependencies for
     * chemical species
     */
    Lookup(fp, cal, chemtbl, kintbl, rttbl);

    /*
     * Apply calibration
     */
    rttbl->pumps[0].Injection_rate *= cal->gwinflux;
    rttbl->pumps[0].flow_rate *= cal->gwinflux;

    for (i = 0; i < nelem; i++)
    {
        for (k = 0; k < rttbl->NumStc; k++)
        {
            elem[i].restart_input.ssa_unsat[k] *= cal->ssa;
            elem[i].restart_input.ssa_gw[k] *= cal->ssa;
        }
    }

    chem_ind = FindChem("'DOC'", chemtbl, rttbl->NumStc);
    if (chem_ind >= 0)
    {
        for (i = 0; i < nelem; i++)
        {
            elem[i].restart_input.tconc_unsat[chem_ind] *= cal->initconc;
            elem[i].restart_input.tconc_gw[chem_ind] *= cal->initconc;
        }

        rttbl->prcp_conc[chem_ind] *= cal->prcpconc;
        if (ctrl->PrpFlg == 2)
        {
            for (i = 0; i < forc->TSD_prepconc.length; i++)
            {
                forc->TSD_prepconc.data[i][chem_ind] *= cal->prcpconc;
            }
        }
    }

#if TEMP_DISABLED
    /* Initializing volumetric parameters, inherit from PIHM
     * That is, if PIHM is started from a hot start, rt is also
     * initialized with the hot data */
    PRCP_VOL = CD->NumVol - 1;
    BOUND_VOL = CD->NumVol;

    for (i = 0; i < nelem; i++)
    {
        double          heqv;
        double          satn;

        /* Initializing volumetrics for groundwater (GW) cells */
        InitVcele(pihm->elem[i].ws.gw, pihm->elem[i].topo.area,
            pihm->elem[i].soil.smcmax, 1.0, LAND_VOL, &CD->Vcele[RT_GW(i)]);

        /* Initializing volumetrics for unsaturated cells */
        /* Porosity in PIHM is
         * Effective Porosity = Porosity - Residue Water Porosity
         * Porosity in RT is total Porosity, therefore, the water height in the
         * unsaturated zone needs be converted as well */
        heqv = EqvUnsatH(pihm->elem[i].soil.smcmax,
            pihm->elem[i].soil.smcmin, pihm->elem[i].soil.depth,
            pihm->elem[i].ws.unsat, pihm->elem[i].ws.gw);
        satn = UnsatSatRatio(pihm->elem[i].soil.depth,
            pihm->elem[i].ws.unsat, pihm->elem[i].ws.gw);

        InitVcele(heqv, pihm->elem[i].topo.area, pihm->elem[i].soil.smcmax,
            satn, LAND_VOL, &CD->Vcele[RT_UNSAT(i)]);

#if defined(_FBR_)
        /* Initializing volumetrics for deep groundwater (FBR GW) cells */
        InitVcele(pihm->elem[i].ws.fbr_gw, pihm->elem[i].topo.area,
            pihm->elem[i].geol.smcmax, 1.0, LAND_VOL,
            &CD->Vcele[RT_FBR_GW(i)]);

        /* Initializing volumetrics for bedrock unsaturated cells */
        heqv = EqvUnsatH(pihm->elem[i].geol.smcmax,
            pihm->elem[i].geol.smcmin, pihm->elem[i].geol.depth,
            pihm->elem[i].ws.fbr_unsat, pihm->elem[i].ws.fbr_gw);
        satn = UnsatSatRatio(pihm->elem[i].geol.depth,
            pihm->elem[i].ws.fbr_unsat, pihm->elem[i].ws.fbr_gw);

        InitVcele(heqv, pihm->elem[i].topo.area, pihm->elem[i].geol.smcmax,
            satn, LAND_VOL, &CD->Vcele[RT_FBR_UNSAT(i)]);
#endif
    }

    /* Initializing volumetrics for river cells */
    for (i = 0; i < nriver; i++)
    {
        InitVcele(pihm->river[i].ws.gw, pihm->river[i].topo.area, 1.0, 1.0,
            RIVER_VOL, &CD->Vcele[RT_RIVER(i)]);
    }

    /* Initialize virtual cell */
    InitVcele(0.0, 0.0, 0.0, 0.0, VIRTUAL_VOL, &CD->Vcele[PRCP_VOL - 1]);
    InitVcele(1.0, 1.0, 1.0, 1.0, VIRTUAL_VOL, &CD->Vcele[BOUND_VOL - 1]);

    /* Initializing concentration distributions */
    PIHMprintf(VL_NORMAL,
        "\n Initializing concentration, Vcele [i, 0 ~ NumVol]... \n");

    for (i = 0; i < CD->NumVol; i++)
    {
        CD->Vcele[i].index = i + 1;

        for (j = 0; j < pihm->rttbl.NumStc; j++)
        {
            if (strcmp(pihm->chemtbl[j].ChemName, "'H+'") == 0)
            {
                CD->Vcele[i].chms.t_conc[j] = CD->Vcele[i].ic.t_conc[j];
                CD->Vcele[i].chms.p_actv[j] = CD->Vcele[i].chms.t_conc[j];
                CD->Vcele[i].chms.p_conc[j] = CD->Vcele[i].chms.t_conc[j];
            }
            else if (pihm->chemtbl[j].itype == MINERAL)
            {
                CD->Vcele[i].chms.t_conc[j] = CD->Vcele[i].ic.t_conc[j];
                /* Update the concentration of mineral after get the molar volume of
                 * mineral */
                CD->Vcele[i].chms.t_conc[j] *= (pihm->rttbl.RelMin == 0) ?
                    /* Absolute mineral volume fraction */
                    1000.0 / pihm->chemtbl[j].MolarVolume / CD->Vcele[i].porosity :
                    /* Relative mineral volume fraction */
                    (1.0 - CD->Vcele[i].porosity + INFTYSMALL) * 1000.0 /
                    pihm->chemtbl[j].MolarVolume / CD->Vcele[i].porosity;
                CD->Vcele[i].chms.p_actv[j] = 1.0;
                CD->Vcele[i].chms.p_conc[j] = CD->Vcele[i].chms.t_conc[j];
                CD->Vcele[i].chms.ssa[j] = CD->Vcele[i].ic.p_para[j];
            }
            else if ((pihm->chemtbl[j].itype == CATION_ECHG) ||
                (pihm->chemtbl[j].itype == ADSORPTION))
            {
                CD->Vcele[i].chms.t_conc[j] = CD->Vcele[i].ic.t_conc[j];
                CD->Vcele[i].chms.p_actv[j] = CD->Vcele[i].chms.t_conc[j] * 0.5;
                /* Change the unit of CEC (eq/g) into C(ion site)
                 * (eq/L porous space), assuming density of solid is always
                 * 2650 g/L solid */
                CD->Vcele[i].chms.t_conc[j] *= (1.0 - CD->Vcele[i].porosity) * 2650.0;
                CD->Vcele[i].chms.p_conc[j] = CD->Vcele[i].chms.t_conc[j];
            }
            else
            {
                CD->Vcele[i].chms.t_conc[j] = CD->Vcele[i].ic.t_conc[j];
                CD->Vcele[i].chms.p_actv[j] = CD->Vcele[i].chms.t_conc[j] * 0.5;
                CD->Vcele[i].chms.p_conc[j] = CD->Vcele[i].chms.t_conc[j] * 0.5;
                CD->Vcele[i].chms.ssa[j] = CD->Vcele[i].ic.p_para[j];
            }
        }

        for (j = 0; j < pihm->rttbl.NumSsc; j++)
        {
            CD->Vcele[i].chms.s_conc[j] = ZERO;
        }
    }

    /*
     * Beginning configuring the connectivity for flux
     */
#if defined(_FBR_)
    CD->NumFac = NUM_EDGE * nelem * 4 + 7 * nelem + 6 * nriver;
#else
    CD->NumFac = NUM_EDGE * nelem * 2 + 3 * nelem + 6 * nriver;
#endif

    /* Configuring the lateral connectivity of GW grid blocks */
    PIHMprintf(VL_NORMAL,
        "\n Configuring the lateral connectivity of GW grid blocks... \n");

    CD->Flux = (face *) malloc(CD->NumFac * sizeof(face));

    for (i = 0; i < nelem; i++)
    {
        int             elemlo;
        int             elemuu;
        int             elemll;
        double          distance;

        for (j = 0; j < NUM_EDGE; j++)
        {
            distance = pihm->elem[i].topo.nabrdist[j];

            if (pihm->elem[i].nabr[j] > 0)
            {
                elemlo = pihm->elem[i].nabr[j];
                elemuu = 0;
                elemll = 0;

                /* Initialize GW fluxes */
                InitFlux(CD->Vcele[RT_GW(i)].index,
                    CD->Vcele[RT_GW(elemlo - 1)].index, 0,
                    (elemuu > 0) ?  CD->Vcele[RT_GW(elemuu - 1)].index : 0,
                    (elemll > 0) ?  CD->Vcele[RT_GW(elemll - 1)].index : 0,
                    DISPERSION, distance, &CD->Flux[RT_LAT_GW(i, j)]);

                /* Initialize unsat zone fluxes */
                InitFlux(CD->Vcele[RT_UNSAT(i)].index,
                    CD->Vcele[RT_UNSAT(elemlo - 1)].index, 0,
                    (elemuu > 0) ?  CD->Vcele[RT_UNSAT(elemuu - 1)].index :0,
                    (elemll > 0) ?  CD->Vcele[RT_UNSAT(elemll - 1)].index : 0,
                    DISPERSION, distance, &CD->Flux[RT_LAT_UNSAT(i, j)]);
            }
            else if (pihm->elem[i].nabr[j] < 0)
            {
                elemlo = -pihm->elem[i].nabr[j];
                elemuu = 0;
                elemll = 0;

                /* Initialize GW fluxes */
                InitFlux(CD->Vcele[RT_GW(i)].index,
                    CD->Vcele[RT_RIVER(elemlo - 1)].index,
                    0, 0, 0,
                    DISPERSION, distance, &CD->Flux[RT_LAT_GW(i, j)]);

                /* Initialize unsat zone fluxes */
                InitFlux(CD->Vcele[RT_UNSAT(i)].index,
                    CD->Vcele[RT_RIVER(elemlo - 1)].index,
                    0, 0, 0,
                    DISPERSION, distance, &CD->Flux[RT_LAT_UNSAT(i, j)]);
            }
            else
            {
                if (pihm->elem[i].attrib.bc_type[j] == NO_FLOW)
                {
                    InitFlux(CD->Vcele[RT_GW(i)].index, BOUND_VOL, 0, 0, 0,
                        NO_FLOW, distance, &CD->Flux[RT_LAT_GW(i, j)]);
                }
                else
                {
                    InitFlux(CD->Vcele[RT_GW(i)].index, PRCP_VOL, 0, 0, 0,
                        NO_DISP, distance, &CD->Flux[RT_LAT_GW(i, j)]);
                }

                InitFlux(CD->Vcele[RT_UNSAT(i)].index, BOUND_VOL, 0, 0, 0,
                    NO_FLOW, distance, &CD->Flux[RT_LAT_UNSAT(i, j)]);
            }
        }

#if defined(_FBR_)
        for (j = 0; j < NUM_EDGE; j++)
        {
            if (pihm->elem[i].nabr[j] == 0)
            {
                distance = elem[i].topo.nabrdist[j];

                if (pihm->elem[i].attrib.fbrbc_type[j] == NO_FLOW)
                {
                    InitFlux(CD->Vcele[RT_FBR_GW(i)].index, BOUND_VOL, 0, 0,
                        0, NO_FLOW, distance, &CD->Flux[RT_LAT_FBR_GW(i, j)]);
                }
                else
                {
                    InitFlux(CD->Vcele[RT_FBR_GW(i)].index, PRCP_VOL, 0, 0, 0,
                        NO_DISP, distance, &CD->Flux[RT_LAT_FBR_GW(i, j)]);
                }

                InitFlux(CD->Vcele[RT_FBR_UNSAT(i)].index, BOUND_VOL, 0, 0, 0,
                    NO_FLOW, distance, &CD->Flux[RT_LAT_FBR_UNSAT(i, j)]);
            }
            else
            {
                if (pihm->elem[i].nabr[j] > 0)
                {
                    elemlo = pihm->elem[i].nabr[j];
                    distance = elem[i].topo.nabrdist[j];
                }
                else
                {
                    elemlo =
                        (pihm->river[-pihm->elem[i].nabr[j] - 1].leftele ==
                            pihm->elem[i].ind) ?
                        &pihm->river[-elem[i].nabr[j] - 1].rightele :
                        &pihm->river[-elem[i].nabr[j] - 1].leftele;
                    int             jj;

                    for (jj = 0; jj < NUM_EDGE; jj++)
                    {
                        if (pihm->elem[elemlo - 1].nabr[jj] == pihm->elem[i].nabr[j])
                        {
                            distance = pihm->elem[i].topo.nabrdist[j] +
                                pihm->elem[elemlo - 1].topo.nabrdist[jj];
                            break;
                        }
                    }
                }

                elemuu = 0;
                elemll = 0;

                /* Initialize GW fluxes */
                InitFlux(CD->Vcele[RT_FBR_GW(i)].index,
                    CD->Vcele[RT_FBR_GW(elemlo - 1)].index, 0,
                    (elemuu > 0) ?  CD->Vcele[RT_FBR_GW(elemuu - 1)].index : 0,
                    (elemll > 0) ?  CD->Vcele[RT_FBR_GW(elemll - 1)].index : 0,
                    DISPERSION, distance, &CD->Flux[RT_LAT_FBR_GW(i, j)]);

                /* Initialize unsat zone fluxes */
                InitFlux(CD->Vcele[RT_FBR_UNSAT(i)].index,
                    CD->Vcele[RT_FBR_UNSAT(elemlo - 1)].index, 0,
                    (elemuu > 0) ?  CD->Vcele[RT_FBR_UNSAT(elemuu - 1)].index :0,
                    (elemll > 0) ?  CD->Vcele[RT_FBR_UNSAT(elemll - 1)].index : 0,
                    DISPERSION, distance, &CD->Flux[RT_LAT_FBR_UNSAT(i, j)]);
            }
        }
#endif

        /* Infiltration */
        InitFlux(CD->Vcele[RT_UNSAT(i)].index, PRCP_VOL, 0, 0, 0, NO_DISP, 0.0,
            &CD->Flux[RT_INFIL(i)]);

        /* Rechage centered at unsat blocks */
        InitFlux(CD->Vcele[RT_UNSAT(i)].index, CD->Vcele[RT_GW(i)].index,
            0, 0, 0, DISPERSION, 0.5 * pihm->elem[i].soil.depth,
            &CD->Flux[RT_RECHG_UNSAT(i)]);

        /* Recharge centered at gw blocks */
        InitFlux(CD->Vcele[RT_GW(i)].index, CD->Vcele[RT_UNSAT(i)].index,
            0, 0, 0, DISPERSION, 0.5 * pihm->elem[i].soil.depth,
            &CD->Flux[RT_RECHG_GW(i)]);

#if defined(_FBR_)
        /* Bedrock leakage */
        InitFlux(CD->Vcele[RT_GW(i)].index, CD->Vcele[RT_FBR_UNSAT(i)].index,
            0, 0, 0, DISPERSION,
            0.1 * (pihm->elem[i].soil.depth + pihm->elem[i].geol.depth),
            &CD->Flux[RT_FBR_LKG(i)]);

        /* Bedrock infiltration */
        InitFlux(CD->Vcele[RT_FBR_UNSAT(i)].index, CD->Vcele[RT_GW(i)].index,
            0, 0, 0, DISPERSION,
            0.1 * (pihm->elem[i].soil.depth + pihm->elem[i].geol.depth),
            &CD->Flux[RT_FBR_INFIL(i)]);

        /* Bedrock recharge centered at unsat blocks */
        InitFlux(CD->Vcele[RT_FBR_UNSAT(i)].index,
            CD->Vcele[RT_FBR_GW(i)].index, 0, 0, 0, DISPERSION,
            0.5 * pihm->elem[i].geol.depth, &CD->Flux[RT_FBR_RECHG_UNSAT(i)]);

        /* Bedrock recharge centered at groundwater blocks */
        InitFlux(CD->Vcele[RT_FBR_GW(i)].index,
            CD->Vcele[RT_FBR_UNSAT(i)].index, 0, 0, 0, DISPERSION,
            0.5 * pihm->elem[i].geol.depth, &CD->Flux[RT_FBR_RECHG_GW(i)]);
#endif
    }

    /* Configuring the vertical connectivity of UNSAT - GW blocks */
    PIHMprintf(VL_NORMAL,
        "\n Configuring the vertical connectivity of UNSAT - GW grid blocks... \n");

    /* Configuring the connectivity of RIVER and EBR blocks */
    PIHMprintf(VL_NORMAL,
        "\n Configuring the connectivity of RIVER & EBR grid blocks... \n");

    for (i = 0; i < nriver; i++)
    {
        double          distance;

        /* Between River and Left */
        /* River to left OFL 2 */
        InitFlux(CD->Vcele[RT_RIVER(i)].index, BOUND_VOL, 0, 0, 0, NO_DISP,
            0.0, &CD->Flux[RT_LEFT_SURF2RIVER(i)]);

        /* Between River and Right */
        /* River to right OFL 3 */
        InitFlux(CD->Vcele[RT_RIVER(i)].index, BOUND_VOL, 0, 0, 0, NO_DISP,
            0.0, &CD->Flux[RT_RIGHT_SURF2RIVER(i)]);

        /* Between Left and EBR */
        /* EBR to left  7 + 4 */
        for (j = 0; j < NUM_EDGE; j++)
        {
            if (-pihm->elem[pihm->river[i].leftele - 1].nabr[j] == i + 1)
            {
                distance =
                    CD->Flux[RT_LAT_GW(pihm->river[i].leftele - 1, j)].distance;
                break;
            }
        }
        InitFlux(CD->Vcele[RT_RIVER(i)].index,
            CD->Vcele[RT_GW(pihm->river[i].leftele - 1)].index,
            0, 0, 0, DISPERSION, distance, &CD->Flux[RT_LEFT_AQIF2RIVER(i)]);

        /* Between Right and EBR */
        /* EBR to right 8 + 5 */
        for (j = 0; j < NUM_EDGE; j++)
        {
            if (-pihm->elem[pihm->river[i].rightele - 1].nabr[j] == i + 1)
            {
                distance =
                    CD->Flux[RT_LAT_GW(pihm->river[i].rightele - 1, j)].distance;
                break;
            }
        }
        InitFlux(CD->Vcele[RT_RIVER(i)].index,
            CD->Vcele[RT_GW(pihm->river[i].rightele - 1)].index,
            0, 0, 0, DISPERSION, distance, &CD->Flux[RT_RIGHT_AQIF2RIVER(i)]);

        /* Between EBR */
        /* To downstream EBR 9 */
        InitFlux(CD->Vcele[RT_RIVER(i)].index,
            (pihm->river[i].down < 0) ?
            BOUND_VOL : CD->Vcele[RT_RIVER(pihm->river[i].down - 1)].index,
            0, 0, 0, NO_DISP, 0.0, &CD->Flux[RT_DOWN_RIVER2RIVER(i)]);

        /* From upstream EBR 10 */
        InitFlux(CD->Vcele[RT_RIVER(i)].index,
            (pihm->river[i].up[0] < 0) ?
            BOUND_VOL : CD->Vcele[RT_RIVER(pihm->river[i].up[0] - 1)].index,
            (pihm->river[i].up[1] < 0) ?
            pihm->river[i].up[1] :
            CD->Vcele[RT_RIVER(pihm->river[i].up[1] - 1)].index,
            0, 0, NO_DISP, 0.0, &CD->Flux[RT_UP_RIVER2RIVER(i)]);
    }

    for (k = 0; k < CD->NumFac; k++)
    {
        CD->Flux[k].velocity = 0.0;
        CD->Flux[k].flux = 0.0; /* Initialize 0.0 for above sections of GW-GW,
                                 * UNSAT-UNSAT, GW-UNSAT, UNSAT-GW */
        CD->Flux[k].flux_trib = 0.0;
        CD->Flux[k].s_area = 0.0;
    }

    if (!pihm->rttbl.RecFlg)
    {
        for (i = 0; i < nelem; i++)
        {
            Speciation(pihm->chemtbl, &pihm->rttbl, 1,
                &CD->Vcele[RT_GW(i)].chms);
            Speciation(pihm->chemtbl, &pihm->rttbl, 1,
                &CD->Vcele[RT_UNSAT(i)].chms);
#if defined(_FBR_)
            Speciation(CD, RT_FBR_GW(i), 1);
#endif
        }
    }

    /* Initialize river concentrations */
    for (i = 0; i < nriver; i++)
    {
        for (k = 0; k < pihm->rttbl.NumStc; k++)
        {
            if (pihm->chemtbl[k].itype != AQUEOUS)
            {
                CD->Vcele[RT_RIVER(i)].chms.t_conc[k] = 1.0E-20;
                CD->Vcele[RT_RIVER(i)].chms.p_conc[k] = 1.0E-20;
                CD->Vcele[RT_RIVER(i)].chms.p_actv[k] = 1.0E-20;
            }
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        for (k = 0; k < NumSpc; k++)
        {
            CD->Vcele[RT_UNSAT(i)].chms.t_mole[k] =
                CD->Vcele[RT_UNSAT(i)].chms.t_conc[k] * CD->Vcele[RT_UNSAT(i)].vol * CD->Vcele[RT_UNSAT(i)].porosity;

            NV_Ith(CV_Y, UNSAT_MOLE(i, k)) = CD->Vcele[RT_UNSAT(i)].chms.t_mole[k];

            CD->Vcele[RT_GW(i)].chms.t_mole[k] =
                CD->Vcele[RT_GW(i)].chms.t_conc[k] * CD->Vcele[RT_GW(i)].vol * CD->Vcele[RT_GW(i)].porosity;

            NV_Ith(CV_Y, GW_MOLE(i, k)) = CD->Vcele[RT_GW(i)].chms.t_mole[k];
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        for (k = 0; k < NumSpc; k++)
        {
            CD->Vcele[RT_RIVER(i)].chms.t_mole[k] =
                CD->Vcele[RT_RIVER(i)].chms.t_conc[k] * CD->Vcele[RT_RIVER(i)].vol * CD->Vcele[RT_RIVER(i)].porosity;

            NV_Ith(CV_Y, RIVER_MOLE(i, k)) = CD->Vcele[RT_RIVER(i)].chms.t_mole[k];
        }
    }

    fclose(fp);
#endif
}

void FreeChem(Chem_Data CD)
{
//    int             i;
//
//    free(CD->BTC_loc);
//    free(CD->prepconcindex);
//
//    // CD->Vcele
//    for (i = 0; i < CD->NumVol; i++)
//    {
//        free(CD->Vcele[i].log10_pconc);
//        free(CD->Vcele[i].log10_sconc);
//        free(CD->Vcele[i].p_para);
//        free(CD->Vcele[i].btcv_pconc);
//    }
//    free(CD->Vcele);
//
//    free(CD->Flux);
//
//    if (CD->NumPUMP > 0)
//    {
//        free(CD->pumps);
//    }
//
//    // CD->TSD_prepconc
//    for (i = 0; i < CD->TSD_prepconc[0].length; i++)
//    {
//        free(CD->TSD_prepconc[0].data[i]);
//    }
//    free(CD->TSD_prepconc[0].data);
//    free(CD->TSD_prepconc[0].ftime);
//    free(CD->TSD_prepconc[0].value);
//    free(CD->TSD_prepconc);
//
//    free(CD->Precipitation.chms.t_conc);
//    free(CD->Precipitation.chms.p_conc);
//    free(CD->Precipitation.p_para);
//
}

void InitVcele(double height, double area, double porosity, double sat,
    int type, vol_conc *Vcele)
{
    Vcele->height_o = height;
    Vcele->height_t = height;
    Vcele->area = area;
    Vcele->porosity = porosity;
    Vcele->vol = height * area;
    Vcele->sat = sat;
    Vcele->type = type;

    if (sat > 1.0)
    {
        PIHMprintf(VL_ERROR,
            "Fatal Error, Unsaturated Zone Initialization For RT Failed!\n");
    }
}

void InitFlux(int nodeup, int nodelo, int node_trib, int nodeuu, int nodell,
    int bc, double distance, face *flux)
{
    flux->nodeup = nodeup;
    flux->nodelo = nodelo;
    flux->node_trib = node_trib;
    flux->nodeuu = nodeuu;
    flux->nodell = nodell;
    flux->BC = bc;
    flux->distance = distance;
}

void UpdateVcele(double height, double sat, vol_conc *Vcele)
{
    Vcele->height_o = Vcele->height_t;
    Vcele->height_t = height;
    Vcele->vol = Vcele->area * height;
    Vcele->sat = sat;
}

double EqvUnsatH(double smcmax, double smcmin, double depth, double unsat,
    double gw)
{
    double          deficit;

    deficit = depth - gw;
    deficit = (deficit < depth) ? deficit : depth;
    deficit = (deficit > 0.0) ? deficit : 0.0;

    return (unsat * (smcmax - smcmin) + deficit * smcmin) / smcmax;
}

double UnsatSatRatio(double depth, double unsat, double gw)
{
    return ((unsat < 0.0) ? 0.0 : ((gw > depth) ? 1.0 : unsat / (depth - gw)));
}

