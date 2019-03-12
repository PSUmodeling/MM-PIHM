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
#define MIN(a,b) (((a)<(b))? (a):(b))
#define MAX(a,b) (((a)>(b))? (a):(b))

void InitChem(char *filename, const char chem_filen[], const pihm_struct pihm,
    Chem_Data CD)
{
    int             i, j, k;
    int             species_counter, min_counter, ads_counter, cex_counter, num_other,
        num_conditions = 0;
    int             line_width = LINE_WIDTH, words_line =
        WORDS_LINE, word_width = WORD_WIDTH;
    int             speciation_flg = 0, specflg;
    double          tmpval[WORDS_LINE];
    char            cmdstr[MAXSTRING];
    int             lno = 0;
    int             PRCP_VOL;
    int             BOUND_VOL;

    assert(pihm != NULL);

    char            line[256];
    char          **tmpstr = (char **)malloc(WORDS_LINE * sizeof(char *));

    for (i = 0; i < words_line; i++)
        tmpstr[i] = (char *)malloc(WORD_WIDTH * sizeof(char));

    FILE           *chem_fp;

    char           *datafn =
        (char *)malloc((strlen(filename) * 2 + 100) * sizeof(char));
    sprintf(datafn, "input/%s/%s.cdbs", filename, filename);
    FILE           *database = fopen(datafn, "r");

    char           *forcfn =
        (char *)malloc((strlen(filename) * 2 + 100) * sizeof(char));
    sprintf(forcfn, "input/%s/%s.prep", filename, filename);

    chem_fp = fopen(chem_filen, "r");
    CheckFile(chem_fp, chem_filen);
    PIHMprintf(VL_VERBOSE, " Reading %s\n", chem_filen);

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
    CD->NumEle = nelem;
    CD->NumRiv = nriver;

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

    /* INITIAL_CONDITIONS block */
    PIHMprintf(VL_NORMAL, " Reading '%s.chem' INITIAL_CONDITIONS: \n", filename);
    rewind(chem_fp);
    fgets(line, line_width, chem_fp);
    while (keymatch(line, "INITIAL_CONDITIONS", tmpval, tmpstr) != 1)
        fgets(line, line_width, chem_fp);
    fgets(line, line_width, chem_fp);
    while (keymatch(line, "END", tmpval, tmpstr) != 1)
    {
        if (keymatch(line, " ", tmpval, tmpstr) != 2)
        {
            num_conditions++;
        }
        fgets(line, line_width, chem_fp);
    }
    PIHMprintf(VL_NORMAL, "  %d conditions assigned. \n", num_conditions);

    char          **chemcon = (char **)malloc(num_conditions * sizeof(char *));
    for (i = 0; i < num_conditions; i++)
        chemcon[i] = (char *)malloc(word_width * sizeof(char));

    int            *condition_index = (int *)malloc(CD->NumVol * sizeof(int));
    /* When user assign conditions to blocks, they start from 1 */

    for (i = 0; i < CD->NumVol; i++)
    {
        condition_index[i] = 0;
    }

    vol_conc       *Condition_vcele =
        (vol_conc *) malloc(num_conditions * sizeof(vol_conc));
    for (i = 0; i < num_conditions; i++)
    {
        Condition_vcele[i].index = i + 1;
        Condition_vcele[i].t_conc =
            (double *)malloc(CD->NumStc * sizeof(double));
        Condition_vcele[i].p_conc =
            (double *)malloc(CD->NumStc * sizeof(double));
        Condition_vcele[i].p_para =
            (double *)malloc(CD->NumStc * sizeof(double));
        Condition_vcele[i].s_conc = NULL;
        /* We do not input cocentration for secondary speices in rt */
        for (j = 0; j < CD->NumStc; j++)
        {
            Condition_vcele[i].t_conc[j] = ZERO;
            Condition_vcele[i].p_conc[j] = ZERO;
        }
    }

    for (i = 0; i < CD->NumStc + CD->NumSsc; i++)
    {
        CD->chemtype[i].DiffCoe = CD->DiffCoe;

        CD->chemtype[i].DispCoe = CD->DispCoe;

        CD->chemtype[i].Charge = 0.0;
        CD->chemtype[i].SizeF = 1.0;
    }

    k = 0;
    int             initfile = 0;
    FILE           *cheminitfile = NULL;
    rewind(chem_fp);
    fgets(line, line_width, chem_fp);
    while (keymatch(line, "INITIAL_CONDITIONS", tmpval, tmpstr) != 1)
        fgets(line, line_width, chem_fp);
    if (strcmp(tmpstr[1], "FILE") == 0)
    {
        /* Initialize chemical distribution from file evoked. This will nullify
         * all the condition assignment given in the next lines.
         * But for now, please keep those lines to let the code work. */

        initfile = 1;
        PIHMprintf(VL_NORMAL, "  Specifiying the initial chemical distribution from file '%s.cini'. \n", filename);

        char           *cheminit =
            (char *)malloc((strlen(filename) * 2 + 100) * sizeof(char));
        sprintf(cheminit, "input/%s/%s.cini", filename, filename);
        cheminitfile = fopen(cheminit, "r");

        if (cheminitfile == NULL)
        {
            PIHMprintf(VL_NORMAL, "  Fatal Error: %s.cini does not exist! \n",
                filename);
            exit(1);
        }
        else
        {
            PIHMprintf(VL_NORMAL, "  Reading the '%s.cini'!! \n", filename);
        }

        free(cheminit);         // 10.02
    }

    fgets(line, line_width, chem_fp);
    while (keymatch(line, "END", tmpval, tmpstr) != 1)
    {
        if (keymatch(line, " ", tmpval, tmpstr) != 2)
        {
            strcpy(chemcon[k++], tmpstr[0]);
            if (initfile == 0)
            {
                PIHMprintf(VL_ERROR,
                    "Assigning initial conditions in .chem file is temporarily"
                    " disabled. Please use a .cini file.\n");
                PIHMexit(EXIT_FAILURE);
            }
        }
        fgets(line, line_width, chem_fp);
    }
    if (initfile == 1)
    {
        for (i = 0; i < CD->NumVol; i++)
        {
            fscanf(cheminitfile, "%d %d", &k, &condition_index[i]);
        }
    }

    if (cheminitfile != NULL)
        fclose(cheminitfile);

    /* CONDITIONS block */
    PIHMprintf(VL_NORMAL, "\n Reading '%s.chem' CONDITIONS: ", filename);
    for (i = 0; i < num_conditions; i++)
    {
        rewind(chem_fp);
        species_counter = 0;
        min_counter = 0;
        ads_counter = 0;
        cex_counter = 0;
        num_other = 0;
        fgets(line, line_width, chem_fp);
        while ((keymatch(line, "Condition", tmpval, tmpstr) != 1) ||
            (keymatch(line, chemcon[i], tmpval, tmpstr) != 1))
            fgets(line, line_width, chem_fp);
        if (strcmp(tmpstr[1], chemcon[i]) == 0)
            PIHMprintf(VL_NORMAL, "\n  %s", line);
        fgets(line, line_width, chem_fp);
        while (keymatch(line, "END", tmpval, tmpstr) != 1)
        {
            if (keymatch(line, "NULL", tmpval, tmpstr) != 2)
            {
                specflg = SpeciationType(database, tmpstr[0]);

                if (specflg == AQUEOUS)
                {
                    /* Arrange the concentration of the primary species in such a
                     * way that all the mobile species are at the beginning. */
                    num_other = min_counter + ads_counter + cex_counter;
                    Condition_vcele[i].t_conc[species_counter - num_other] =
                        tmpval[0];
                    PIHMprintf(VL_NORMAL, "  %-28s %g \n",
                        tmpstr[0], tmpval[0]);
                }
                if (specflg == MINERAL)
                {
                    Condition_vcele[i].t_conc[CD->NumSpc + CD->NumAds +
                        CD->NumCex + min_counter] = tmpval[0];
                    if (strcmp(tmpstr[2], "-ssa") == 0)
                        Condition_vcele[i].p_para[CD->NumSpc + CD->NumAds +
                            CD->NumCex + min_counter] = tmpval[1] * 1.0;
                    PIHMprintf(VL_NORMAL,
                        "  mineral %-20s %6.4f \t specific surface area \t%6.4f \n",
                        tmpstr[0], tmpval[0], tmpval[1]);
                    min_counter++;
                }
                if ((tmpstr[0][0] == '>') || (specflg == ADSORPTION))
                {
                    /* Adsorptive sites and species start with > */
                    /* Condition_vcele[i].t_conc[CD->NumSpc + ads_counter] = tmpval[0] * CS->Cal.Site_den;  09.25 temporal comment-out */
                    Condition_vcele[i].t_conc[CD->NumSpc + ads_counter] =
                        tmpval[0] * 1.0;
                    Condition_vcele[i].p_para[CD->NumSpc + ads_counter] = 0;
                    /* Update when fill in the parameters for adsorption */
                    PIHMprintf(VL_NORMAL, "  surface complex %s\t\t%6.4f \n",
                        tmpstr[0], tmpval[0]);
                    ads_counter++;
                    /* under construction */
                }
                if (specflg == CATION_ECHG)
                {
                    Condition_vcele[i].t_conc[CD->NumSpc + CD->NumAds +
                        cex_counter] = tmpval[0];
                    Condition_vcele[i].p_para[CD->NumSpc + CD->NumAds +
                        cex_counter] = 0;
                    /* update when fill in the parameters for cation exchange. */
                    PIHMprintf(VL_NORMAL, "  cation exchange %s\t\t%6.4f \n",
                        tmpstr[0], tmpval[0]);
                    cex_counter++;
                    /* under construction */
                }
                species_counter++;
            }
            fgets(line, line_width, chem_fp);
        }
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

    CD->Vcele = (vol_conc *) malloc(CD->NumVol * sizeof(vol_conc));

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

    for (i = 0; i < CD->NumSpc; i++)
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
                    -(Condition_vcele[condition_index[i] - 1].t_conc[j]));
                CD->Vcele[i].t_conc[j] = CD->Vcele[i].p_conc[j];
                CD->Vcele[i].p_actv[j] = CD->Vcele[i].p_conc[j];
            }
            else if (CD->chemtype[j].itype == MINERAL)
            {
                CD->Vcele[i].t_conc[j] =
                    Condition_vcele[condition_index[i] - 1].t_conc[j];
                CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j];
                CD->Vcele[i].p_actv[j] = 1.0;
                CD->Vcele[i].p_para[j] =
                    Condition_vcele[condition_index[i] - 1].p_para[j];
            }
            else
            {
                CD->Vcele[i].t_conc[j] =
                    Condition_vcele[condition_index[i] - 1].t_conc[j];
                CD->Vcele[i].t_conc[j] *=
                    (strcmp(CD->chemtype[j].ChemName, "DOC") == 0) ?
                    CD->CalInitconc : 1.0;
                CD->Vcele[i].p_conc[j] = CD->Vcele[i].t_conc[j] * 0.5;
                CD->Vcele[i].p_actv[j] = CD->Vcele[i].p_conc[j];
                CD->Vcele[i].p_para[j] =
                    Condition_vcele[condition_index[i] - 1].p_para[j];
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
            CD->chemtype[i + CD->NumSpc + CD->NumAds + CD->NumCex].ChemName,
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
            CD->chemtype[j + CD->NumSpc + CD->NumAds + CD->NumCex].ChemName);
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
    for (i = 0; i < CD->NumSpc; i++)
    {
        CD->chemtype[i].mtype = (CD->chemtype[i].itype == AQUEOUS) ?
             1 : 0;

        for (j = 0; j < CD->NumStc + CD->NumSsc; j++)
        {
            if (CD->Totalconc[i][j] != 0 &&
                CD->chemtype[j].itype != CD->chemtype[i].mtype)
            {
                CD->chemtype[i].mtype = 2;
            }
        }
        /*
         * if (strcmp( CD->chemtype[i].ChemName, "'H+'") == 0)
         * CD->chemtype[i].mtype = 1;
         */
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

    for (i = 0; i < num_conditions; i++)
    {
        free(Condition_vcele[i].t_conc);
        free(Condition_vcele[i].p_conc);
        free(Condition_vcele[i].p_para);
    }
    free(Condition_vcele);

    for (i = 0; i < num_conditions; i++)
        free(chemcon[i]);
    free(chemcon);

    free(datafn);
    free(forcfn);
    free(condition_index);

    for (i = 0; i < words_line; i++)
    {
        free(tmpstr[i]);
    }
    free(tmpstr);

    fclose(chem_fp);
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

int realcheck(const char *words)
{
    int             flg = 1, i;
    if (((words[0] >= '0') && (words[0] <= '9')) ||
        (words[0] == '.') || (words[0] == '-') || (words[0] == '+'))
    {
        for (i = 0; i < (int)strlen(words); i++)
        {
            /* Ascii 10 is new line and 13 is carriage return */
            if ((words[i] > '9' || words[i] < '+') && (words[i] != 'E')
                && (words[i] != 'e') && (words[i] != 10) && (words[i] != 13))
            {
                flg = 0;
            }
        }
    }
    else
    {
        flg = 0;
    }
    return (flg);
}

int keymatch(const char *line, const char *keyword, double *value, char **strval)
{
    /* A very general and convinient way of reading datafile and input file */
    /* find keyword in line, assign the value after keyword to value array if there is any */
    /* store both numbers and strings in order for later use, buffer required */
    /* if is keyword not found return 0. If comments, return 2. Otherwise return 1 */
    int             i;

    for (i = 0; i < WORDS_LINE; i++)
        value[i] = 0.0;

    if ((line[0] == '!') || (line[0] == '#'))
    {
        /* assign a special flag for comments */
        return (2);
    }

    int             j, k;
    int             words_line = WORDS_LINE;
    int             keyfoundflag = 0;

    char          **words;
    words = (char **)malloc(WORDS_LINE * sizeof(char *));

    for (i = 0; i < WORDS_LINE; i++)
    {
        words[i] = (char *)malloc(WORD_WIDTH * sizeof(char));
        memset(words[i], 0, WORD_WIDTH);
    }
    i = j = k = 0;

    /* Partition the line into words */
    while (i < (int)strlen(line))
    {
        if (line[i] != 39)
        {
            while (line[i] != 9 && line[i] != 0 && line[i] != 10
                && line[i] != 32 && line[i] != 13)
            {
                words[k][j++] = line[i++];
                if (line[i] == 9 || line[i] == 32 || line[i] == 13)
                {
                    k++;
                    j = 0;
                }
            }
        }
        else
        {
            words[k][j++] = line[i++];
            while (line[i] != 39)
            {
                words[k][j++] = line[i++];
            }
            words[k++][j] = line[i++];
            j = 0;
        }
        i++;
    }

    words_line = k + 1;

    for (i = 0; i < words_line; i++)
        if (strcmp(words[i], keyword) == 0)
            keyfoundflag = 1;

    j = k = 0;
    for (i = 0; i < words_line; i++)
    {
        strcpy(strval[k++], words[i]);
        if (realcheck(words[i]) == 1)
            value[j++] = atof(words[i]);
    }

    for (i = 0; i < WORDS_LINE; i++)
        free(words[i]);
    free(words);
    return (keyfoundflag);

}


void fluxtrans(int t, int stepsize, const pihm_struct pihm, Chem_Data CD,
    double *t_duration_transp, double *t_duration_react)
{
    /* unit of t: second
     * unit of stepsize: second
     * swi irreducible water saturation
     * hn  non mobile water height
     * ht  transient zone height
     */
    int             i, k = 0;
    int             BOUND_VOL = CD->NumVol;
    int             PRCP_VOL = CD->NumVol - 1;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        for (j = 0; j < 3; j++)
        {
            /* Flux for GW lateral flow */
            CD->Flux[RT_LAT_GW(i, j)].flux = pihm->elem[i].wf.subsurf[j];

            /* Flux for UNSAT lateral flow */
            CD->Flux[RT_LAT_UNSAT(i, j)].s_area = 0.5 *
                pihm->elem[i].topo.edge[j] *
                (CD->Vcele[CD->Flux[RT_LAT_UNSAT(i, j)].nodeup - 1].height_t +
                CD->Vcele[CD->Flux[RT_LAT_UNSAT(i, j)].nodelo - 1].height_t);

#if defined(_FBR_)
            /* Flux for deep lateral flow */
            CD->Flux[RT_LAT_FBR_GW(i, j)].flux = pihm->elem[i].wf.fbrflow[j];

            /* Flux for bedrock unsat lateral flow */
            CD->Flux[RT_LAT_FBR_UNSAT(i, j)].s_area = 0.5 *
                pihm->elem[i].topo.edge[j] *
                (CD->Vcele[CD->Flux[RT_LAT_FBR_UNSAT(i, j)].nodeup - 1].height_t +
                CD->Vcele[CD->Flux[RT_LAT_FBR_UNSAT(i, j)].nodelo - 1].height_t);

#endif
        }

        /* Flux for UNSAT - GW vertical flow */
        CD->Flux[RT_RECHG_UNSAT(i)].flux = pihm->elem[i].wf.rechg *
            CD->Vcele[RT_UNSAT(i)].area;

        CD->Flux[RT_RECHG_GW(i)].flux = -pihm->elem[i].wf.rechg *
            CD->Vcele[RT_GW(i)].area;

        CD->Flux[RT_INFIL(i)].flux =
            -((pihm->elem[i].wf.infil > 0.0) ? pihm->elem[i].wf.infil : 0.0) *
            pihm->elem[i].topo.area;

#if defined(_FBR_)
        CD->Flux[RT_FBR_RECHG_UNSAT(i)].flux = pihm->elem[i].wf.fbr_rechg * CD->Vcele[RT_FBR_RECHG_UNSAT(i)].area;
        CD->Flux[RT_FBR_RECHG_GW(i)].flux = -pihm->elem[i].wf.fbr_rechg * CD->Vcele[RT_FBR_RECHG_GW(i)].area;
#endif
    }

    /* Flux for RIVER flow */
    for (i = 0; i < nriver; i++)
    {
        if (pihm->river[i].down < 0)
        {
            CD->riv = pihm->river[i].wf.rivflow[1] * 86400;
        }
    }

    if ((t - pihm->ctrl.starttime) % DAYINSEC == 0)
    {
        CD->rivd = CD->riv / 1440;  /* Averaging the sum of 1440 mins for a
                                     * daily discharge, rivFlx1 */
        CD->riv = 0;
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        CD->Flux[RT_LEFT_SURF2RIVER(i)].flux = pihm->river[i].wf.rivflow[2];
        CD->Flux[RT_RIGHT_SURF2RIVER(i)].flux = pihm->river[i].wf.rivflow[3];
        CD->Flux[RT_LEFT_AQIF2RIVER(i)].flux = pihm->river[i].wf.rivflow[7] +
            pihm->river[i].wf.rivflow[4];
        CD->Flux[RT_RIGHT_AQIF2RIVER(i)].flux = pihm->river[i].wf.rivflow[8] +
            pihm->river[i].wf.rivflow[5];
        CD->Flux[RT_DOWN_RIVER2RIVER(i)].flux = pihm->river[i].wf.rivflow[9] +
            pihm->river[i].wf.rivflow[1];
        CD->Flux[RT_UP_RIVER2RIVER(i)].flux = pihm->river[i].wf.rivflow[10] +
            pihm->river[i].wf.rivflow[0];

        if (CD->Flux[RT_UP_RIVER2RIVER(i)].node_trib > 0)
        {
            CD->Flux[RT_UP_RIVER2RIVER(i)].flux_trib =
                -(pihm->river[pihm->river[i].up[1] - 1].wf.rivflow[9] +
                pihm->river[pihm->river[i].up[1] - 1].wf.rivflow[1]);
        }
    }

    /* Update the concentration in precipitation here. */
    if (CD->PrpFlg == 2)
    {
        IntrplForc(&CD->TSD_prepconc[0], t, CD->TSD_prepconc[0].nspec,
            NO_INTRPL);

#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (i = 0; i < CD->TSD_prepconc[0].nspec; i++)
        {
            if (CD->prepconcindex[i] > 0)
            {
                int             ind;

                ind = CD->prepconcindex[i] - 1;
                if (CD->Precipitation.t_conc[ind] !=
                    CD->TSD_prepconc[0].value[i])
                {
                    CD->Precipitation.t_conc[ind] =
                        CD->TSD_prepconc[0].value[i];
                    PIHMprintf(VL_NORMAL,
                        "  %s in precipitation is changed to %6.4g\n",
                        CD->chemtype[ind].ChemName,
                        CD->Precipitation.t_conc[ind]);
                }
            }
        }
    }

#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        double          heqv;
        double          satn;

        UpdateVcele(MAX(pihm->elem[i].ws.gw, 1.0E-5), 1.0,
            &CD->Vcele[RT_GW(i)]);

        heqv = EqvUnsatH(pihm->elem[i].soil.smcmax,
            pihm->elem[i].soil.smcmin, pihm->elem[i].soil.depth,
            pihm->elem[i].ws.unsat, pihm->elem[i].ws.gw);

        satn = UnsatSatRatio(pihm->elem[i].soil.depth,
            pihm->elem[i].ws.unsat, pihm->elem[i].ws.gw);

        /* Update the unsaturated zone (vadoze) */
        UpdateVcele(MAX(heqv, 1.0E-5), satn, &CD->Vcele[RT_UNSAT(i)]);

#if defined(_FBR_)
        UpdateVcele(MAX(pihm->elem[i].ws.fbr_gw, 1.0E-5), 1.0,
            &CD->Vcele[RT_FBR_GW(i)]);

        heqv = EqvUnsatH(pihm->elem[i].geol.smcmax,
            pihm->elem[i].geol.smcmin, pihm->elem[i].geol.depth,
            pihm->elem[i].ws.fbr_unsat, pihm->elem[i].ws.fbr_gw);

        satn = UnsatSatRatio(pihm->elem[i].geol.depth,
            pihm->elem[i].ws.fbr_unsat, pihm->elem[i].ws.fbr_gw);

        /* Update the unsaturated zone (vadoze) */
        UpdateVcele(MAX(heqv, 1.0E-5), satn, &CD->Vcele[RT_FBR_UNSAT(i)]);
#endif
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    /* Update river cells */
    for (i = 0; i < nriver; i++)
    {
        UpdateVcele(MAX(pihm->river[i].ws.gw, 1.0E-5) +
            MAX(pihm->river[i].ws.stage, 1.0E-5) /
            CD->Vcele[RT_RIVER(i)].porosity, 1.0, &CD->Vcele[RT_RIVER(i)]);
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        /* For gw cells, contact area is needed for dispersion; */
        for (j = 0; j < 3; j++)
        {
            double              h1, h2;

            if (CD->Flux[RT_LAT_GW(i, j)].BC == NO_FLOW)
            {
                continue;
            }

            if (pihm->elem[i].nabr[j] > 0)
            {
                h1 = 0.5 *
                    (CD->Vcele[RT_GW(i)].height_o +
                    CD->Vcele[RT_GW(i)].height_t);
                h2 = 0.5 *
                    (CD->Vcele[RT_GW(pihm->elem[i].nabr[j] - 1)].height_o +
                    CD->Vcele[RT_GW(pihm->elem[i].nabr[j] - 1)].height_t);

                CD->Flux[RT_LAT_GW(i, j)].s_area =
                    (CD->Flux[RT_LAT_GW(i, j)].flux > 0.0) ?
                    pihm->elem[i].topo.edge[j] * h1 :
                    pihm->elem[i].topo.edge[j] * h2;
            }
            else if (pihm->elem[i].nabr[j] < 0)
            {
                h1 = 0.5 *
                    (CD->Vcele[RT_GW(i)].height_o +
                    CD->Vcele[RT_GW(i)].height_t);
                h2 = 0.5 *
                    (CD->Vcele[RT_RIVER(-pihm->elem[i].nabr[j] - 1)].height_o +
                    CD->Vcele[RT_RIVER(-pihm->elem[i].nabr[j] - 1)].height_t);

                CD->Flux[RT_LAT_GW(i, j)].s_area =
                    (CD->Flux[RT_LAT_GW(i, j)].flux > 0.0) ?
                    pihm->elem[i].topo.edge[j] * h1 :
                    pihm->elem[i].topo.edge[j] * h2;
            }

            /* Calculate velocity according to flux and area */
            CD->Flux[RT_LAT_GW(i, j)].velocity =
                (CD->Flux[RT_LAT_GW(i, j)].s_area > 1.0E-4) ?
                CD->Flux[RT_LAT_GW(i, j)].flux /
                CD->Flux[RT_LAT_GW(i, j)].s_area : 1.0E-10 / 86400;
        }

        CD->Flux[RT_RECHG_UNSAT(i)].s_area = pihm->elem[i].topo.area;
        CD->Flux[RT_RECHG_UNSAT(i)].velocity =
            CD->Flux[RT_RECHG_UNSAT(i)].flux / pihm->elem[i].topo.area;

        CD->Flux[RT_RECHG_GW(i)].s_area = pihm->elem[i].topo.area;
        CD->Flux[RT_RECHG_GW(i)].velocity =
            CD->Flux[RT_RECHG_GW(i)].flux / pihm->elem[i].topo.area;
    }

    /* Correct river flux area and velocity */
#ifdef _OPENMP
# pragma omp parallel for
#endif
    for (i = 0; i < nriver; i++)
    {
        int             j;

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (-pihm->elem[pihm->river[i].leftele - 1].nabr[j] == i + 1)
            {
                CD->Flux[RT_LEFT_AQIF2RIVER(i)].s_area =
                CD->Flux[RT_LAT_GW(pihm->river[i].leftele - 1, j)].s_area;
                CD->Flux[RT_LEFT_AQIF2RIVER(i)].velocity =
                -CD->Flux[RT_LAT_GW(pihm->river[i].leftele - 1, j)].velocity;
                break;
            }
        }

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (-pihm->elem[pihm->river[i].rightele - 1].nabr[j] == i + 1)
            {
                CD->Flux[RT_RIGHT_AQIF2RIVER(i)].s_area =
                CD->Flux[RT_LAT_GW(pihm->river[i].rightele - 1, j)].s_area;
                CD->Flux[RT_RIGHT_AQIF2RIVER(i)].velocity =
                -CD->Flux[RT_LAT_GW(pihm->river[i].rightele - 1, j)].velocity;
                break;
            }
        }
    }

    /* Update virtual volume */

    if (CD->PrpFlg)
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (k = 0; k < CD->NumSpc; k++)
        {
            CD->Vcele[PRCP_VOL - 1].t_conc[k] =
                (strcmp(CD->chemtype[k].ChemName, "'DOC'") == 0) ?
                CD->Precipitation.t_conc[k] * CD->Condensation *
                CD->CalPrcpconc :
                CD->Precipitation.t_conc[k] * CD->Condensation;
        }
    }
    else
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (k = 0; k < CD->NumSpc; k++)
        {
            CD->Vcele[PRCP_VOL - 1].t_conc[k] = 0.0;
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (k = 0; k < CD->NumStc; k++)
    {
        CD->Vcele[BOUND_VOL - 1].t_conc[k] =
            CD->Precipitation.t_conc[k] * CD->Condensation;
        CD->Vcele[BOUND_VOL - 1].p_conc[k] =
            CD->Precipitation.t_conc[k] * CD->Condensation;
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < CD->NumVol; i++)
    {
        CD->Vcele[i].rt_step = 0.0;
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < CD->NumFac; i++)
    {
        int             j;
        double          peclet;

        if (CD->Flux[i].BC == DISPERSION)
        {
            for (j = 0; j < CD->NumSpc; j++)
            {
                peclet = fabs(CD->Flux[i].velocity * CD->Flux[i].distance /
                    (CD->chemtype[j].DiffCoe +
                    CD->chemtype[j].DispCoe * CD->Flux[i].velocity));
                peclet = MAX(peclet, 1.0E-8);
            }

            CD->Vcele[CD->Flux[i].nodeup - 1].rt_step +=
                fabs(CD->Flux[i].flux / CD->Vcele[CD->Flux[i].nodeup - 1].vol) *
                (1 + peclet) / peclet;
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < CD->NumVol; i++)
    {
        if (CD->Vcele[i].type != VIRTUAL_VOL)
        {
            CD->Vcele[i].rt_step = 0.6 / CD->Vcele[i].rt_step;
            CD->Vcele[i].rt_step = (CD->Vcele[i].rt_step >= stepsize) ?
                stepsize : CD->Vcele[i].rt_step;
        }
    }

    /*
     * RT step control begins
     */
    if (t - pihm->ctrl.starttime >= CD->RT_delay)
    {
        /*
         * Transport
         */
        AdptTime(CD, (double)stepsize, t_duration_transp, t_duration_react);

        /*
         * Reaction
         */
        if ((!CD->RecFlg) &&
            (t - pihm->ctrl.starttime) % (CD->AvgScl * stepsize) == 0)
        {
#ifdef _OPENMP
# pragma omp parallel for
#endif
            for (i = 0; i < nelem; i++)
            {
                React((double)stepsize, CD, &CD->Vcele[RT_GW(i)]);
                React((double)stepsize, CD, &CD->Vcele[RT_UNSAT(i)]);
            }
        }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (i = 0; i < CD->NumEle; i++)
        {
            int             j;

            for (j = 0; j < CD->NumStc; j++)
            {
                if (CD->chemtype[j].itype == MINERAL)
                {
                    /* Averaging mineral concentration to ensure mass
                     * conservation !! */
                    CD->Vcele[RT_GW(i)].t_conc[j] =
                        (CD->Vcele[RT_GW(i)].t_conc[j] *
                        CD->Vcele[RT_GW(i)].height_t +
                        CD->Vcele[RT_UNSAT(i)].t_conc[j] *
                        (pihm->elem[i].soil.depth -
                        CD->Vcele[RT_GW(i)].height_t)) /
                        pihm->elem[i].soil.depth;
                    CD->Vcele[RT_UNSAT(i)].t_conc[j] =
                        CD->Vcele[RT_GW(i)].t_conc[j];
                    CD->Vcele[RT_GW(i)].p_conc[j] =
                        CD->Vcele[RT_GW(i)].t_conc[j];
                    CD->Vcele[RT_UNSAT(i)].p_conc[j] =
                        CD->Vcele[RT_GW(i)].t_conc[j];
                }
            }
        }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (i = 0; i < CD->NumVol; i++)
        {
            int             j;

            if (CD->Vcele[i].type == VIRTUAL_VOL)
            {
                continue;
            }

            /* Make sure intrapolation worked well */
            if (fabs(CD->Vcele[i].height_t - CD->Vcele[i].height_int) >
                1.0E-6)
                PIHMprintf(VL_NORMAL, "%d %6.4f\t%6.4f\n", i,
                    CD->Vcele[i].height_t, CD->Vcele[i].height_int);
            assert(fabs(CD->Vcele[i].height_t - CD->Vcele[i].height_int) <
                1.0E-6);
            if (CD->Vcele[i].illness >= 20)
            {
                for (j = 0; j < CD->NumStc; j++)
                    CD->Vcele[i].t_conc[j] = 1.0E-10;
                PIHMprintf(VL_NORMAL,
                    " Cell %d isolated due to proneness to err!\n",
                    CD->Vcele[i].index);
            }
        }
    } /* RT step control ends */

    /* Reset fluxes for next averaging stage */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (k = 0; k < CD->NumFac; k++)
    {
        CD->Flux[k].velocity = 0.0;
        CD->Flux[k].flux = 0.0;
        CD->Flux[k].flux_trib = 0.0;
        CD->Flux[k].s_area = 0.0;
    }

    /* Every hour */
    if ((t - pihm->ctrl.starttime) % 3600 == 0)
    {
        CD->SPCFlg = 0;

        if (!CD->RecFlg)
        {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
            for (i = 0; i < CD->NumStc; i++)
            {
                int             j;

                for (j = 0; j < nriver; j++)
                {
                    CD->Vcele[RT_RIVER(j)].p_conc[i] =
                        (CD->chemtype[i].itype == MINERAL) ?
                        CD->Vcele[RT_RIVER(j)].t_conc[i] :
                        fabs(CD->Vcele[RT_RIVER(j)].t_conc[i] * 0.1);
                }
            }
        }

        if (!CD->RecFlg)
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (i = 0; i < nriver; i++)
            {
                Speciation(CD, RT_RIVER(i));
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (i = 0; i < CD->NumVol; i++)
            {
                if (CD->Vcele[i].type != VIRTUAL_VOL)
                {
                    Speciation(CD, i);
                }
            }
        }
    }

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < CD->NumVol; i++)
    {
        int             j;

        for (j = 0; j < CD->NumStc; j++)
        {
            CD->Vcele[i].log10_pconc[j] = log10(CD->Vcele[i].p_conc[j]);
        }
        for (j = 0; j < CD->NumSsc; j++)
        {
            CD->Vcele[i].log10_sconc[j] = log10(CD->Vcele[i].s_conc[j]);
        }
    }

    for (k = 0; k < CD->NumBTC; k++)
    {
        int             j;

        for (j = 0; j < CD->NumStc; j++)
        {
            if ((CD->BTC_loc[k] >= CD->pumps[0].Pump_Location - 1) &&
                (j == CD->pumps[0].Position_Species))
            {
                CD->Vcele[CD->BTC_loc[k]].btcv_pconc[j] =
                    log10((CD->Vcele[CD->BTC_loc[k]].p_conc[j] * CD->rivd +
                    CD->pumps[0].Injection_conc * CD->pumps[0].flow_rate) /
                    (CD->rivd + CD->pumps[0].flow_rate));
            }
            else
            {
                CD->Vcele[CD->BTC_loc[k]].btcv_pconc[j] =
                    CD->Vcele[CD->BTC_loc[k]].log10_pconc[j];
            }
        }
    }
}

void AdptTime(Chem_Data CD, double stepsize,
    double *t_duration_transp, double *t_duration_react)
{
    int             i, k;
    time_t          t_start_transp, t_end_transp;

    /* stepsize is in the unit of second */
    t_start_transp = time(NULL);

    time_t          t_start_react, t_end_react;

#if TEMP_DISABLED
    if (int_flg)
    {
        /* Do interpolation. Note that height_int always store the end time
         * height. */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (i = 0; i < CD->NumVol; i++)
        {
            if (CD->Vcele[i].type != VIRTUAL_VOL)
            {
                CD->Vcele[i].height_t =
                    CD->Vcele[i].height_o + CD->Vcele[i].height_sp * stepsize;
                CD->Vcele[i].vol = CD->Vcele[i].area * CD->Vcele[i].height_t;
            }
        }
    }
#endif

#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        for (j = 0; j < CD->NumSpc; j++)
        {
            if (CD->chemtype[j].mtype == 2)
            {
                for (k = 0; k < CD->NumSsc; k++)
                {
                    if ((CD->Totalconc[j][k + CD->NumStc] != 0) &&
                        (CD->chemtype[k + CD->NumStc].itype != AQUEOUS))
                    {
                        CD->Vcele[RT_GW(i)].t_conc[j] = CD->Vcele[RT_GW(i)].t_conc[j] -
                            CD->Totalconc[j][k + CD->NumStc] *
                            CD->Vcele[RT_GW(i)].s_conc[k] * CD->TimRiv;
                        CD->Vcele[RT_UNSAT(i)].t_conc[j] = CD->Vcele[RT_UNSAT(i)].t_conc[j] -
                            CD->Totalconc[j][k + CD->NumStc] *
                            CD->Vcele[RT_UNSAT(i)].s_conc[k] * CD->TimRiv;
                    }
                }
            }
        }
    }

    OS3D(stepsize, CD);

    /* Total concentration except for adsorptions have been transported and
     * adjusted by the volume. For example, if no transport but volume
     * increased by rain, the concentration need be decreased. However, the
     * adsorption part has not been treated yet, so they need be adjusted by
     * the volume change.
     * The porosity is not changed during the period, so the ratio between
     * pore space before and after OS3D is the same ratio between volume of
     * porous media before and after OS3D. */
#if defined(_OPENMP)
# pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int             j;

        for (j = 0; j < CD->NumSpc; j++)
        {
            if (CD->chemtype[j].mtype == 2)
            {
                for (k = 0; k < CD->NumSsc; k++)
                {
                    if ((CD->Totalconc[j][k + CD->NumStc] != 0) &&
                        (CD->chemtype[k + CD->NumStc].itype != AQUEOUS))
                    {
                        CD->Vcele[RT_GW(i)].t_conc[j] =
                            CD->Vcele[RT_GW(i)].t_conc[j] + CD->Totalconc[j][k +
                            CD->NumStc] * CD->Vcele[RT_GW(i)].s_conc[k] *
                            CD->TimRiv;
                        CD->Vcele[RT_UNSAT(i)].t_conc[j] =
                            CD->Vcele[RT_UNSAT(i)].t_conc[j] + CD->Totalconc[j][k +
                            CD->NumStc] * CD->Vcele[RT_UNSAT(i)].s_conc[k] *
                            CD->TimRiv;
                    }
                }
            }
        }
    }

    t_end_transp = time(NULL);
    *t_duration_transp += (t_end_transp - t_start_transp);

    t_start_react = time(NULL);

#if TEMP_DISABLED
    if (int_flg)
    {
#if defined(_OPENMP)
# pragma omp parallel for
#endif
        for (i = 0; i < CD->NumVol; i++)
        {
            CD->Vcele[i].height_o = CD->Vcele[i].height_t;
            CD->Vcele[i].vol_o = CD->Vcele[i].area * CD->Vcele[i].height_o;
        }
    }
#endif

    t_end_react = time(NULL);
    *t_duration_react += (t_end_react - t_start_react);
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
