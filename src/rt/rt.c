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
    pihm_struct pihm, Chem_Data CD, N_Vector CV_Y)
{
    int             i, j, k;
    int             speciation_flg = 0;
    int             PRCP_VOL;
    int             BOUND_VOL;
    FILE           *fp;

    fp = fopen(cdbs_filen, "r");
    CheckFile(fp, cdbs_filen);

    /*
     * Begin updating variables
     */
#if defined(_FBR_)
    CD->NumVol = 4 * nelem + nriver + 2;
#else
    CD->NumVol = 2 * nelem + nriver + 2;
#endif
    CD->Vcele = (vol_conc *) malloc(CD->NumVol * sizeof(vol_conc));

    ReadCini(pihm->filename.cini, pihm->chemtbl, pihm->rttbl.NumStc, CD->Vcele);

    /*
     * Initialize chemical parameters
     */
    for (i = 0; i < pihm->rttbl.NumStc + pihm->rttbl.NumSsc; i++)
    {
        pihm->chemtbl[i].DiffCoe = pihm->rttbl.DiffCoe;
        pihm->chemtbl[i].DispCoe = pihm->rttbl.DispCoe;
        pihm->chemtbl[i].Charge = 0.0;
        pihm->chemtbl[i].SizeF = 1.0;
    }

    for (i = 0; i < MAXSPS; i++)
    {
        for (j = 0; j < MAXSPS; j++)
        {
            pihm->rttbl.Dependency[i][j] = 0.0;      /* NumSsc x NumSdc */
            pihm->rttbl.Dep_kinetic[i][j] = 0.0;     /* (NumMkr + NumAkr) x NumStc */
            pihm->rttbl.Dep_kinetic_all[i][j] = 0.0; /* NumMin x NumStc */
            pihm->rttbl.Totalconc[i][j] = 0.0;       /* NumStc x (NumStc + NumSsc) */
#if NOT_YET_IMPLEMENTED
            pihm->rttbl.Totalconck[i][j] = 0.0;      /* NumStc x (NumStc + NumSsc) */
#endif
        }

        /* Keqs of equilibrium/ kinetic and kinetic all */
        pihm->rttbl.Keq[i] = 0.0;                    /* NumSsc */
        pihm->rttbl.KeqKinect[i] = 0.0;              /* NumMkr + NumAkr */
        pihm->rttbl.KeqKinect_all[i] = 0.0;          /* NumMin */
    }

    /* Primary species table */
    PIHMprintf(VL_NORMAL,
        "\n Primary species and their types: [1], aqueous; [2], adsorption; [3], cation exchange; [4], mineral. \n");
    /* Number of total species in the rt simulator */
    for (i = 0; i < pihm->rttbl.NumStc; i++)
    {
        PIHMprintf(VL_NORMAL, "  %-20s %10d\n", pihm->chemtbl[i].ChemName,
            pihm->chemtbl[i].itype);
    }

    CD->CalPorosity = pihm->cal.porosity;
    CD->CalRate = pihm->cal.rate;
    CD->CalSSA = pihm->cal.ssa;
    CD->CalPrcpconc = pihm->cal.prcpconc;
    CD->CalInitconc = pihm->cal.initconc;
    CD->CalXsorption = pihm->cal.Xsorption;

    for (i = 0; i < NumSpc; i++)
    {
        if (strcmp(pihm->chemtbl[i].ChemName, "pH") == 0)
        {
            strcpy(pihm->chemtbl[i].ChemName, "H+");
            speciation_flg = 1;
        }
    }

    Lookup(fp, pihm->chemtbl, pihm->kintbl, &pihm->rttbl, CD);

    /*
     * Apply calibration
     */
    pihm->rttbl.pumps[0].Injection_rate *= pihm->cal.gwinflux;
    pihm->rttbl.pumps[0].flow_rate *= pihm->cal.gwinflux;

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
        CD->Vcele[i].t_conc = (double *)calloc(pihm->rttbl.NumStc, sizeof(double));
        CD->Vcele[i].t_mole = (double *)calloc(NumSpc, sizeof(double));
        CD->Vcele[i].transp_flux = (double *)calloc(NumSpc, sizeof(double));
        CD->Vcele[i].react_flux = (double *)calloc(NumSpc, sizeof(double));
        CD->Vcele[i].p_conc = (double *)calloc(pihm->rttbl.NumStc, sizeof(double));
        CD->Vcele[i].s_conc = (double *)calloc(pihm->rttbl.NumSsc, sizeof(double));
        CD->Vcele[i].p_actv = (double *)calloc(pihm->rttbl.NumStc, sizeof(double));
        CD->Vcele[i].p_para = (double *)calloc(pihm->rttbl.NumStc, sizeof(double));
        CD->Vcele[i].log10_pconc = (double *)calloc(pihm->rttbl.NumStc, sizeof(double));
        CD->Vcele[i].log10_sconc = (double *)calloc(pihm->rttbl.NumSsc, sizeof(double));
        CD->Vcele[i].btcv_pconc = (double *)calloc(pihm->rttbl.NumStc, sizeof(double));

        CD->Vcele[i].illness = 0;

        for (j = 0; j < pihm->rttbl.NumStc; j++)
        {
            if (strcmp(pihm->chemtbl[j].ChemName, "'H+'") == 0)
            {
                CD->Vcele[i].t_conc[j] = CD->Vcele[i].ic.t_conc[j];
                CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                CD->Vcele[i].p_actv[j] = CD->Vcele[i].p_conc[j];
            }
            else if (pihm->chemtbl[j].itype == MINERAL)
            {
                CD->Vcele[i].t_conc[j] = CD->Vcele[i].ic.t_conc[j];
                CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                CD->Vcele[i].p_actv[j] = 1.0;
                CD->Vcele[i].p_para[j] = CD->Vcele[i].ic.p_para[j];
            }
            else
            {
                CD->Vcele[i].t_conc[j] = CD->Vcele[i].ic.t_conc[j];
                CD->Vcele[i].t_conc[j] *=
                    (strcmp(pihm->chemtbl[j].ChemName, "'DOC'") == 0) ?
                    CD->CalInitconc : 1.0;
                CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j] * 0.5;
                CD->Vcele[i].p_actv[j] = CD->Vcele[i].p_conc[j];
                CD->Vcele[i].p_para[j] = CD->Vcele[i].ic.p_para[j];
            }
        }
        for (j = 0; j < pihm->rttbl.NumSsc; j++)
        {
            CD->Vcele[i].s_conc[j] = ZERO;
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
                elemuu = upstream(pihm->elem[i],
                    pihm->elem[pihm->elem[i].nabr[j] - 1], pihm);
                elemll = upstream(pihm->elem[pihm->elem[i].nabr[j] - 1],
                    pihm->elem[i], pihm);

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

                elemuu = upstream(pihm->elem[i],
                    pihm->elem[elemlo - 1], pihm);
                elemll = upstream(pihm->elem[elemlo - 1],
                    pihm->elem[i], pihm);

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

    /* Update the concentration of mineral after get the molar volume of
     * mineral */

    double          Cal_PCO2 = 1.0;
    double          Cal_Keq = 1.0;
    for (i = 0; i < pihm->rttbl.NumAkr + pihm->rttbl.NumMkr; i++)
    {
        pihm->rttbl.KeqKinect[i] += (!strcmp(
            pihm->chemtbl[i + NumSpc + pihm->rttbl.NumAds + pihm->rttbl.NumCex].ChemName,
            "'CO2(*g)'")) ?
            log10(Cal_PCO2) : log10(Cal_Keq);
    }

    PIHMprintf(VL_NORMAL, "\n Kinetic Mass Matrx (calibrated Keq)! \n");
    PIHMprintf(VL_NORMAL, "%-15s", " ");
    for (i = 0; i < pihm->rttbl.NumStc; i++)
        PIHMprintf(VL_NORMAL, "%-14s", pihm->chemtbl[i].ChemName);
    PIHMprintf(VL_NORMAL, "\n");
    for (j = 0; j < pihm->rttbl.NumMkr + pihm->rttbl.NumAkr; j++)
    {
        PIHMprintf(VL_NORMAL, " %-14s",
            pihm->chemtbl[j + NumSpc + pihm->rttbl.NumAds + pihm->rttbl.NumCex].ChemName);
        for (i = 0; i < pihm->rttbl.NumStc; i++)
        {
            PIHMprintf(VL_NORMAL, "%-14.2f", pihm->rttbl.Dep_kinetic[j][i]);
        }
        PIHMprintf(VL_NORMAL, " Keq = %-6.2f\n", pihm->rttbl.KeqKinect[j]);
    }
    PIHMprintf(VL_NORMAL, "\n");
    /* Use calibration coefficient to produce new Keq values for
     * 1) CO2, 2) other kinetic reaction */

    PIHMprintf(VL_NORMAL,
        " \n Mass action species type determination (0: immobile, 1: mobile, 2: Mixed) \n");
    for (i = 0; i < NumSpc; i++)
    {
        pihm->chemtbl[i].mtype = (pihm->chemtbl[i].itype == AQUEOUS) ?
             MOBILE_MA : IMMOBILE_MA;

        for (j = 0; j < pihm->rttbl.NumStc + pihm->rttbl.NumSsc; j++)
        {
            if (pihm->rttbl.Totalconc[i][j] != 0 &&
                pihm->chemtbl[j].itype != pihm->chemtbl[i].mtype)
            {
                pihm->chemtbl[i].mtype = MIXED_MA;
            }
        }
        PIHMprintf(VL_NORMAL, " %12s\t%10d\n", pihm->chemtbl[i].ChemName,
            pihm->chemtbl[i].mtype);
    }

    PIHMprintf(VL_NORMAL,
        " \n Individual species type determination (1: aqueous, 2: adsorption, 3: ion exchange, 4: solid) \n");
    for (i = 0; i < pihm->rttbl.NumStc + pihm->rttbl.NumSsc; i++)
    {
        PIHMprintf(VL_NORMAL, " %12s\t%10d\n", pihm->chemtbl[i].ChemName,
            pihm->chemtbl[i].itype);
    }

    for (i = 0; i < CD->NumVol; i++)
    {
        for (j = 0; j < pihm->rttbl.NumStc; j++)
        {
            if (pihm->chemtbl[j].itype == MINERAL)
            {
                if (pihm->rttbl.RelMin == 0)
                {
                    /* Absolute mineral volume fraction */
                    CD->Vcele[i].t_conc[j] =
                        CD->Vcele[i].t_conc[j] * 1000 /
                        pihm->chemtbl[j].MolarVolume / CD->Vcele[i].porosity;
                    CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                }
                if (pihm->rttbl.RelMin == 1)
                {
                    /* Relative mineral volume fraction */
                    /* Porosity can be 1.0 so the relative fraction option needs
                     * a small modification */
                    CD->Vcele[i].t_conc[j] = CD->Vcele[i].t_conc[j] *
                        (1 - CD->Vcele[i].porosity + INFTYSMALL) * 1000 /
                        pihm->chemtbl[j].MolarVolume / CD->Vcele[i].porosity;
                    CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                }
            }
            if ((pihm->chemtbl[j].itype == CATION_ECHG) ||
                (pihm->chemtbl[j].itype == ADSORPTION))
            {
                /* Change the unit of CEC (eq/g) into C(ion site)
                 * (eq/L porous space), assuming density of solid is always
                 * 2650 g/L solid */
                CD->Vcele[i].t_conc[j] =
                    CD->Vcele[i].t_conc[j] * (1 - CD->Vcele[i].porosity) * 2650;
                CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
            }
        }
    }

    if (!pihm->rttbl.RecFlg)
    {
        for (i = 0; i < nelem; i++)
        {
            Speciation(pihm->chemtbl, &pihm->rttbl, CD, RT_GW(i), speciation_flg);
#if defined(_FBR_)
            Speciation(CD, RT_FBR_GW(i), speciation_flg);
#endif
        }
    }

    /* Initialize unsaturated zone concentrations to be the same as in saturated
     * zone */
    for (i = 0; i < nelem; i++)
    {
        for (k = 0; k < pihm->rttbl.NumStc; k++)
        {
            CD->Vcele[RT_UNSAT(i)].t_conc[k] = CD->Vcele[RT_GW(i)].t_conc[k];
            CD->Vcele[RT_UNSAT(i)].p_conc[k] = CD->Vcele[RT_GW(i)].p_conc[k];
            CD->Vcele[RT_UNSAT(i)].p_actv[k] = CD->Vcele[RT_GW(i)].p_actv[k];
#if defined(_FBR_)
            CD->Vcele[RT_FBR_UNSAT(i)].t_conc[k] = CD->Vcele[RT_GW(i)].t_conc[k];
            CD->Vcele[RT_FBR_UNSAT(i)].p_conc[k] = CD->Vcele[RT_GW(i)].p_conc[k];
            CD->Vcele[RT_FBR_UNSAT(i)].p_actv[k] = CD->Vcele[RT_GW(i)].p_actv[k];
            CD->Vcele[RT_FBR_GW(i)].t_conc[k] = CD->Vcele[RT_GW(i)].t_conc[k];
            CD->Vcele[RT_FBR_GW(i)].p_conc[k] = CD->Vcele[RT_GW(i)].p_conc[k];
            CD->Vcele[RT_FBR_GW(i)].p_actv[k] = CD->Vcele[RT_GW(i)].p_actv[k];
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
                CD->Vcele[RT_RIVER(i)].t_conc[k] = 1.0E-20;
                CD->Vcele[RT_RIVER(i)].p_conc[k] = 1.0E-20;
                CD->Vcele[RT_RIVER(i)].p_actv[k] = 1.0E-20;
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
            CD->Vcele[RT_UNSAT(i)].t_mole[k] =
                CD->Vcele[RT_UNSAT(i)].t_conc[k] * CD->Vcele[RT_UNSAT(i)].vol * CD->Vcele[RT_UNSAT(i)].porosity;

            NV_Ith(CV_Y, UNSAT_MOLE(i, k)) = CD->Vcele[RT_UNSAT(i)].t_mole[k];

            CD->Vcele[RT_GW(i)].t_mole[k] =
                CD->Vcele[RT_GW(i)].t_conc[k] * CD->Vcele[RT_GW(i)].vol * CD->Vcele[RT_GW(i)].porosity;

            NV_Ith(CV_Y, GW_MOLE(i, k)) = CD->Vcele[RT_GW(i)].t_mole[k];
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        for (k = 0; k < NumSpc; k++)
        {
            CD->Vcele[RT_RIVER(i)].t_mole[k] =
                CD->Vcele[RT_RIVER(i)].t_conc[k] * CD->Vcele[RT_RIVER(i)].vol * CD->Vcele[RT_RIVER(i)].porosity;

            NV_Ith(CV_Y, RIVER_MOLE(i, k)) = CD->Vcele[RT_RIVER(i)].t_mole[k];
        }
    }

    fclose(fp);
}

int upstream(elem_struct up, elem_struct lo, const pihm_struct pihm)
{
    /* Locate the upstream grid of up -> lo flow */
    /* Require verification                      */
    /* only determines points in triangular elements */
    double          x_, y_;
    int             i;

    x_ = 2 * up.topo.x - lo.topo.x;
    y_ = 2 * up.topo.y - lo.topo.y;

    for (i = 0; i < nelem; i++)
    {
        double          x_a, x_b, x_c;
        double          y_a, y_b, y_c;
        double          dot00, dot01, dot02, dot11, dot12, u, v, invDenom;

        /* Find point lies in which triangular element, a very interesting
         * method */
        if ((i != (up.ind - 1)) && (i != (lo.ind - 1)))
        {
            x_a = pihm->meshtbl.x[pihm->elem[i].node[0] - 1];
            x_b = pihm->meshtbl.x[pihm->elem[i].node[1] - 1];
            x_c = pihm->meshtbl.x[pihm->elem[i].node[2] - 1];
            y_a = pihm->meshtbl.y[pihm->elem[i].node[0] - 1];
            y_b = pihm->meshtbl.y[pihm->elem[i].node[1] - 1];
            y_c = pihm->meshtbl.y[pihm->elem[i].node[2] - 1];
            dot00 = (x_c - x_a) * (x_c - x_a) + (y_c - y_a) * (y_c - y_a);
            dot01 = (x_c - x_a) * (x_b - x_a) + (y_c - y_a) * (y_b - y_a);
            dot02 = (x_c - x_a) * (x_ - x_a) + (y_c - y_a) * (y_ - y_a);
            dot11 = (x_b - x_a) * (x_b - x_a) + (y_b - y_a) * (y_b - y_a);
            dot12 = (x_b - x_a) * (x_ - x_a) + (y_b - y_a) * (y_ - y_a);
            invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
            u = (dot11 * dot02 - dot01 * dot12) * invDenom;
            v = (dot00 * dot12 - dot01 * dot02) * invDenom;
            if ((u > 0.0) && (v > 0.0) && (u + v < 1.0))
            {
                return pihm->elem[i].ind;
            }
        }
    }

    return 0;
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
//        free(CD->Vcele[i].t_conc);
//        free(CD->Vcele[i].p_conc);
//        free(CD->Vcele[i].s_conc);
//        free(CD->Vcele[i].log10_pconc);
//        free(CD->Vcele[i].log10_sconc);
//        free(CD->Vcele[i].p_actv);
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
//    free(CD->Precipitation.t_conc);
//    free(CD->Precipitation.p_conc);
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
    Vcele->vol_o = height * area;
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
    Vcele->height_int = height;
    Vcele->vol_o = Vcele->area * Vcele->height_o;
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

void SortChem(char chemn[MAXSPS][MAXSTRING], const int p_type[MAXSPS], int nsps,
    chemtbl_struct chemtbl[])
{
    int             i, j;
    int             temp;
    int             rank[MAXSPS];
    int             ranked_type[MAXSPS];

    for (i = 0; i < nsps; i++)
    {
        rank[i] = i;
        ranked_type[i] = p_type[i];
    }

    for (i = 0; i < nsps - 1; i++)
    {
        for (j = 0; j < nsps - i - 1; j++)
        {
            if (ranked_type[j] > ranked_type[j + 1])
            {
                temp = rank[j];
                rank[j] = rank[j + 1];
                rank[j + 1] = temp;

                temp = ranked_type[j];
                ranked_type[j] = ranked_type[j + 1];
                ranked_type[j + 1] = temp;
            }
        }
    }

    for (i = 0; i < nsps; i++)
    {
        strcpy(chemtbl[i].ChemName, chemn[rank[i]]);
        chemtbl[i].itype = p_type[rank[i]];
    }
}

int FindChem(const char chemn[MAXSTRING], const chemtbl_struct  chemtbl[], int nsps)
{
    int             i;
    int             ind = -1;

    for (i = 0; i < nsps; i++)
    {
        if (strcmp(chemn, chemtbl[i].ChemName) == 0)
        {
            ind = i;
            break;
        }
    }

    return ind;
}
