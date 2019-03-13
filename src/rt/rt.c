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

void InitChem(char *filename, const char cini_filen[], const pihm_struct pihm,
    Chem_Data CD)
{
    int             i, j, k;
    int             speciation_flg = 0;
    int             PRCP_VOL;
    int             BOUND_VOL;

    assert(pihm != NULL);

    FILE           *database = fopen(pihm->filename.cdbs, "r");

    if (database == NULL)
    {
        fprintf(stderr, "\n  Fatal Error: %s.cdbs does not exist! \n",
            filename);
        exit(1);
    }

    /*
     * Begin updating variables
     */
#if defined(_FBR_)
    CD->NumVol = 4 * nelem + nriver + 2;
#else
    CD->NumVol = 2 * nelem + nriver + 2;
#endif
    CD->Vcele = (vol_conc *) malloc(CD->NumVol * sizeof(vol_conc));

    ReadCini(pihm->filename.cini, CD->chemtype, CD->NumStc, CD->Vcele);

    for (i = 0; i < CD->NumStc + CD->NumSsc; i++)
    {
        CD->chemtype[i].DiffCoe = CD->DiffCoe;

        CD->chemtype[i].DispCoe = CD->DispCoe;

        CD->chemtype[i].Charge = 0.0;
        CD->chemtype[i].SizeF = 1.0;
    }

    PRCP_VOL = CD->NumVol - 1;
    BOUND_VOL = CD->NumVol;

    CD->Dependency = (double **)malloc(CD->NumSsc * sizeof(double *));
    for (i = 0; i < CD->NumSsc; i++)
    {
        CD->Dependency[i] = (double *)malloc(CD->NumSdc * sizeof(double));
        /* Convert secondary species as an expression of primary species */
        for (j = 0; j < CD->NumSdc; j++)
            CD->Dependency[i][j] = 0.0;
    }

    CD->Dep_kinetic =
        (double **)malloc((CD->NumMkr + CD->NumAkr) * sizeof(double *));
    for (i = 0; i < CD->NumMkr + CD->NumAkr; i++)
    {
        CD->Dep_kinetic[i] = (double *)malloc(CD->NumStc * sizeof(double));
        /* Express kinetic species as function of primary species */
        for (j = 0; j < CD->NumStc; j++)
            CD->Dep_kinetic[i][j] = 0.0;
    }

    CD->Dep_kinetic_all = (double **)malloc((CD->NumMin) * sizeof(double *));
    for (i = 0; i < CD->NumMin; i++)
    {
        CD->Dep_kinetic_all[i] = (double *)malloc(CD->NumStc * sizeof(double));
        /* Dependencies of minearls, all */
        for (j = 0; j < CD->NumStc; j++)
            CD->Dep_kinetic_all[i][j] = 0.0;
    }

    /* Keqs of equilibrium/ kinetic and kinetic all */
    CD->Keq = (double *)malloc(CD->NumSsc * sizeof(double));
    CD->KeqKinect =
        (double *)malloc((CD->NumMkr + CD->NumAkr) * sizeof(double));
    CD->KeqKinect_all = (double *)malloc(CD->NumMin * sizeof(double));

    /* Convert total concentration as an expression of all species */
    CD->Totalconc = (double **)malloc(CD->NumStc * sizeof(double *));
    for (i = 0; i < CD->NumStc; i++)
        CD->Totalconc[i] =
            (double *)malloc((CD->NumStc + CD->NumSsc) * sizeof(double));

#if NOT_YET_IMPLEMENTED
    /* Convert total concentration as an expression of all species */
    CD->Totalconck = (double **)malloc(CD->NumStc * sizeof(double *));
    for (i = 0; i < CD->NumStc; i++)
        CD->Totalconck[i] =
            (double *)malloc((CD->NumStc + CD->NumSsc) * sizeof(double));
#endif

    for (i = 0; i < CD->NumStc; i++)
        for (j = 0; j < CD->NumStc + CD->NumSsc; j++)
        {
            CD->Totalconc[i][j] = 0.0;
#if NOT_YET_IMPLEMENTED
            CD->Totalconck[i][j] = 0.0;
#endif
        }

    /* Primary species table */
    PIHMprintf(VL_NORMAL,
        "\n Primary species and their types: [1], aqueous; [2], adsorption; [3], cation exchange; [4], mineral. \n");
    /* Number of total species in the rt simulator */
    for (i = 0; i < CD->NumStc; i++)
    {
        PIHMprintf(VL_NORMAL, "  %-20s %10d\n", CD->chemtype[i].ChemName,
            CD->chemtype[i].itype);
    }


    for (i = 0; i < CD->NumPUMP; i++)
    {
        CD->pumps[i].Position_Species = -1;
        for (j = 0; j < CD->NumStc; j++)
        {
            if (!strcmp(CD->pumps[i].Name_Species,
                    CD->chemtype[j].ChemName))
            {
                CD->pumps[i].Position_Species = j;
            }
        }
    }

    /* End of reading input files */


    /* Initializing volumetric parameters, inherit from PIHM
     * That is, if PIHM is started from a hot start, rt is also
     * initialized with the hot data */
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

    CD->CalPorosity = pihm->cal.porosity;
    CD->CalRate = pihm->cal.rate;
    CD->CalSSA = pihm->cal.ssa;
    CD->CalPrcpconc = pihm->cal.prcpconc;
    CD->CalInitconc = pihm->cal.initconc;
    CD->CalXsorption = pihm->cal.Xsorption;

    /* Initializing volumetrics for river cells */
    for (i = 0; i < nriver; i++)
    {
        InitVcele(pihm->river[i].ws.gw, pihm->river[i].topo.area, 1.0, 1.0,
            RIVER_VOL, &CD->Vcele[RT_RIVER(i)]);
    }

    /* Initialize virtual cell */
    InitVcele(0.0, 0.0, 0.0, 0.0, VIRTUAL_VOL, &CD->Vcele[PRCP_VOL - 1]);
    InitVcele(1.0, 1.0, 1.0, 1.0, VIRTUAL_VOL, &CD->Vcele[BOUND_VOL - 1]);

    for (i = 0; i < NumSpc; i++)
    {
        if (strcmp(CD->chemtype[i].ChemName, "pH") == 0)
        {
            strcpy(CD->chemtype[i].ChemName, "H+");
            speciation_flg = 1;
        }
    }

    /* Initializing concentration distributions */
    PIHMprintf(VL_NORMAL,
        "\n Initializing concentration, Vcele [i, 0 ~ NumVol]... \n");

    for (i = 0; i < CD->NumVol; i++)
    {
        CD->Vcele[i].index = i + 1;
        CD->Vcele[i].t_conc = (double *)calloc(CD->NumStc, sizeof(double));
        CD->Vcele[i].p_conc = (double *)calloc(CD->NumStc, sizeof(double));
        CD->Vcele[i].s_conc = (double *)calloc(CD->NumSsc, sizeof(double));
        CD->Vcele[i].p_actv = (double *)calloc(CD->NumStc, sizeof(double));
        CD->Vcele[i].p_para = (double *)calloc(CD->NumStc, sizeof(double));
        CD->Vcele[i].log10_pconc = (double *)calloc(CD->NumStc, sizeof(double));
        CD->Vcele[i].log10_sconc = (double *)calloc(CD->NumSsc, sizeof(double));
        CD->Vcele[i].btcv_pconc = (double *)calloc(CD->NumStc, sizeof(double));

        CD->Vcele[i].illness = 0;

        for (j = 0; j < CD->NumStc; j++)
        {
            if (speciation_flg == 1 &&
                strcmp(CD->chemtype[j].ChemName, "H+") == 0)
            {
                CD->Vcele[i].p_conc[j] = pow(10,
                    -(CD->Vcele[i].ic.t_conc[j]));
                CD->Vcele[i].t_conc[j] = CD->Vcele[i].p_conc[j];
                CD->Vcele[i].p_actv[j] = CD->Vcele[i].p_conc[j];
            }
            else if (CD->chemtype[j].itype == MINERAL)
            {
                CD->Vcele[i].t_conc[j] =
                    CD->Vcele[i].ic.t_conc[j];
                CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                CD->Vcele[i].p_actv[j] = 1.0;
                CD->Vcele[i].p_para[j] =
                    CD->Vcele[i].ic.p_para[j];
            }
            else
            {
                CD->Vcele[i].t_conc[j] =
                    CD->Vcele[i].ic.t_conc[j];
                CD->Vcele[i].t_conc[j] *=
                    (strcmp(CD->chemtype[j].ChemName, "DOC") == 0) ?
                    CD->CalInitconc : 1.0;
                CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j] * 0.5;
                CD->Vcele[i].p_actv[j] = CD->Vcele[i].p_conc[j];
                CD->Vcele[i].p_para[j] =
                    CD->Vcele[i].ic.p_para[j];
            }
        }
        for (j = 0; j < CD->NumSsc; j++)
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
            0, 0, 0, DISPERSION, 0.1 * pihm->elem[i].soil.depth,
            &CD->Flux[RT_RECHG_UNSAT(i)]);

        /* Recharge centered at gw blocks */
        InitFlux(CD->Vcele[RT_GW(i)].index, CD->Vcele[RT_UNSAT(i)].index,
            0, 0, 0, DISPERSION, 0.1 * pihm->elem[i].soil.depth,
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
            0.1 * pihm->elem[i].geol.depth, &CD->Flux[RT_FBR_RECHG_UNSAT(i)]);

        /* Bedrock recharge centered at groundwater blocks */
        InitFlux(CD->Vcele[RT_FBR_GW(i)].index,
            CD->Vcele[RT_FBR_UNSAT(i)].index, 0, 0, 0, DISPERSION,
            0.1 * pihm->elem[i].geol.depth, &CD->Flux[RT_FBR_RECHG_GW(i)]);
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

    CD->SPCFlg = speciation_flg;
    Lookup(database, CD);
    /* Update the concentration of mineral after get the molar volume of
     * mineral */

    double          Cal_PCO2 = 1.0;
    double          Cal_Keq = 1.0;
    for (i = 0; i < CD->NumAkr + CD->NumMkr; i++)
    {
        CD->KeqKinect[i] += (!strcmp(
            CD->chemtype[i + NumSpc + CD->NumAds + CD->NumCex].ChemName,
            "'CO2(*g)'")) ?
            log10(Cal_PCO2) : log10(Cal_Keq);
    }

    PIHMprintf(VL_NORMAL, "\n Kinetic Mass Matrx (calibrated Keq)! \n");
    PIHMprintf(VL_NORMAL, "%-15s", " ");
    for (i = 0; i < CD->NumStc; i++)
        PIHMprintf(VL_NORMAL, "%-14s", CD->chemtype[i].ChemName);
    PIHMprintf(VL_NORMAL, "\n");
    for (j = 0; j < CD->NumMkr + CD->NumAkr; j++)
    {
        PIHMprintf(VL_NORMAL, " %-14s",
            CD->chemtype[j + NumSpc + CD->NumAds + CD->NumCex].ChemName);
        for (i = 0; i < CD->NumStc; i++)
        {
            PIHMprintf(VL_NORMAL, "%-14.2f", CD->Dep_kinetic[j][i]);
        }
        PIHMprintf(VL_NORMAL, " Keq = %-6.2f\n", CD->KeqKinect[j]);
    }
    PIHMprintf(VL_NORMAL, "\n");
    /* Use calibration coefficient to produce new Keq values for
     * 1) CO2, 2) other kinetic reaction */

    PIHMprintf(VL_NORMAL,
        " \n Mass action species type determination (0: immobile, 1: mobile, 2: Mixed) \n");
    for (i = 0; i < NumSpc; i++)
    {
        CD->chemtype[i].mtype = (CD->chemtype[i].itype == AQUEOUS) ?
             MOBILE_MA : IMMOBILE_MA;

        for (j = 0; j < CD->NumStc + CD->NumSsc; j++)
        {
            if (CD->Totalconc[i][j] != 0 &&
                CD->chemtype[j].itype != CD->chemtype[i].mtype)
            {
                CD->chemtype[i].mtype = MIXED_MA;
            }
        }
        PIHMprintf(VL_NORMAL, " %12s\t%10d\n", CD->chemtype[i].ChemName,
            CD->chemtype[i].mtype);
    }

    PIHMprintf(VL_NORMAL,
        " \n Individual species type determination (1: aqueous, 2: adsorption, 3: ion exchange, 4: solid) \n");
    for (i = 0; i < CD->NumStc + CD->NumSsc; i++)
    {
        PIHMprintf(VL_NORMAL, " %12s\t%10d\n", CD->chemtype[i].ChemName,
            CD->chemtype[i].itype);
    }

    for (i = 0; i < CD->NumVol; i++)
    {
        for (j = 0; j < CD->NumStc; j++)
        {
            if (CD->chemtype[j].itype == MINERAL)
            {
                if (CD->RelMin == 0)
                {
                    /* Absolute mineral volume fraction */
                    CD->Vcele[i].t_conc[j] =
                        CD->Vcele[i].t_conc[j] * 1000 /
                        CD->chemtype[j].MolarVolume / CD->Vcele[i].porosity;
                    CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                }
                if (CD->RelMin == 1)
                {
                    /* Relative mineral volume fraction */
                    /* Porosity can be 1.0 so the relative fraction option needs
                     * a small modification */
                    CD->Vcele[i].t_conc[j] = CD->Vcele[i].t_conc[j] *
                        (1 - CD->Vcele[i].porosity + INFTYSMALL) * 1000 /
                        CD->chemtype[j].MolarVolume / CD->Vcele[i].porosity;
                    CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                }
            }
            if ((CD->chemtype[j].itype == CATION_ECHG) ||
                (CD->chemtype[j].itype == ADSORPTION))
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

    CD->SPCFlg = 1;
    if (!CD->RecFlg)
    {
        for (i = 0; i < nelem; i++)
        {
            Speciation(CD, RT_GW(i));
#if defined(_FBR_)
            Speciation(CD, RT_FBR_GW(i));
#endif
        }
    }
    CD->SPCFlg = 0;

    /* Initialize unsaturated zone concentrations to be the same as in saturated
     * zone */
    for (i = 0; i < nelem; i++)
    {
        for (k = 0; k < CD->NumStc; k++)
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
        for (k = 0; k < CD->NumStc; k++)
        {
            if (CD->chemtype[k].itype != AQUEOUS)
            {
                CD->Vcele[RT_RIVER(i)].t_conc[k] = 1.0E-20;
                CD->Vcele[RT_RIVER(i)].p_conc[k] = 1.0E-20;
                CD->Vcele[RT_RIVER(i)].p_actv[k] = 1.0E-20;
            }
        }
    }

    fclose(database);
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
    int             i;

    free(CD->BTC_loc);
    free(CD->prepconcindex);

    for (i = 0; i < CD->NumSsc; i++)
    {
        free(CD->Dependency[i]);
    }
    free(CD->Dependency);

    for (i = 0; i < CD->NumMkr + CD->NumAkr; i++)
    {
        free(CD->Dep_kinetic[i]);
    }
    free(CD->Dep_kinetic);

    for (i = 0; i < CD->NumMin; i++)
    {
        free(CD->Dep_kinetic_all[i]);
    }
    free(CD->Dep_kinetic_all);

    for (i = 0; i < CD->NumStc; i++)
    {
        free(CD->Totalconc[i]);
#if NOT_YET_IMPLEMENTED
        free(CD->Totalconck[i]);
#endif
    }
    free(CD->Totalconc);
#if NOT_YET_IMPLEMENTED
    free(CD->Totalconck);
#endif

    free(CD->Keq);
    free(CD->KeqKinect);
    free(CD->KeqKinect_all);

    // CD->Vcele
    for (i = 0; i < CD->NumVol; i++)
    {
        free(CD->Vcele[i].t_conc);
        free(CD->Vcele[i].p_conc);
        free(CD->Vcele[i].s_conc);
        free(CD->Vcele[i].log10_pconc);
        free(CD->Vcele[i].log10_sconc);
        free(CD->Vcele[i].p_actv);
        free(CD->Vcele[i].p_para);
        free(CD->Vcele[i].btcv_pconc);
    }
    free(CD->Vcele);

    free(CD->Flux);

    if (CD->NumPUMP > 0)
    {
        free(CD->pumps);
    }

    // CD->TSD_prepconc
    for (i = 0; i < CD->TSD_prepconc[0].length; i++)
    {
        free(CD->TSD_prepconc[0].data[i]);
    }
    free(CD->TSD_prepconc[0].data);
    free(CD->TSD_prepconc[0].ftime);
    free(CD->TSD_prepconc[0].value);
    free(CD->TSD_prepconc);

    free(CD->Precipitation.t_conc);
    free(CD->Precipitation.p_conc);
    free(CD->Precipitation.p_para);

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
    species chem[])
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
        //strcpy(chem[rank[i]].ChemName, chemn[i]);
        //chem[rank[i]].itype = p_type[i];
        strcpy(chem[i].ChemName, chemn[rank[i]]);
        chem[i].itype = p_type[rank[i]];
    }
}

int FindChem(const char chemn[MAXSTRING], const species chemtype[], int nsps)
{
    int             i;
    int             ind = -1;

    for (i = 0; i < nsps; i++)
    {
        if (strcasecmp(chemn, chemtype[i].ChemName) == 0)
        {
            ind = i;
            break;
        }
    }

    if (ind < 0)
    {
        PIHMprintf(VL_ERROR, "Error finding chemical %s.\n");
        PIHMexit(EXIT_FAILURE);
    }

    return ind;
}
